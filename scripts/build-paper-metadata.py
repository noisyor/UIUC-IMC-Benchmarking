#!/usr/bin/env python3
"""Build full paper bylines and affiliations from DOI-indexed OpenAlex data."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import unicodedata
from datetime import date
from difflib import SequenceMatcher
from pathlib import Path
from urllib.parse import unquote, urlencode, urlparse
from urllib.request import Request, urlopen


OPENALEX_API = "https://api.openalex.org/works"


def load_javascript_object(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    return json.loads(text[text.index("{") : text.rindex("}") + 1])


def normalize(text: str) -> str:
    value = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    return " ".join(re.findall(r"[a-z0-9]+", value.lower()))


def title_score(expected: str, candidate: str) -> float:
    left = normalize(expected)
    right = normalize(candidate)
    sequence = SequenceMatcher(None, left, right).ratio()
    left_tokens = set(left.split())
    right_tokens = set(right.split())
    overlap = len(left_tokens & right_tokens) / max(1, len(left_tokens | right_tokens))
    return 0.62 * sequence + 0.38 * overlap


def doi_from_url(url: str) -> str:
    parsed = urlparse(url)
    if parsed.netloc.lower() != "doi.org":
        return ""
    return unquote(parsed.path.lstrip("/")).lower()


def fetch_openalex(dois: list[str], batch_size: int = 50) -> dict[str, dict]:
    works: dict[str, dict] = {}
    for offset in range(0, len(dois), batch_size):
        batch = dois[offset : offset + batch_size]
        query = urlencode(
            {
                "filter": "doi:" + "|".join(batch),
                "select": "id,doi,title,publication_year,authorships",
                "per-page": 100,
            }
        )
        request = Request(
            f"{OPENALEX_API}?{query}",
            headers={"User-Agent": "UIUC-IMC-Benchmarking-metadata-builder/1.0"},
        )
        with urlopen(request, timeout=90) as response:
            payload = json.load(response)
        for work in payload.get("results", []):
            doi = doi_from_url(work.get("doi", ""))
            if doi:
                works[doi] = work
        print(f"Fetched OpenAlex batch {offset // batch_size + 1}: {len(batch)} DOI records", flush=True)
        time.sleep(0.25)
    return works


def unique_strings(values: list[str]) -> list[str]:
    seen = set()
    result = []
    for value in values:
        cleaned = " ".join(value.split())
        key = cleaned.casefold()
        if cleaned and key not in seen:
            seen.add(key)
            result.append(cleaned)
    return result


def author_record(authorship: dict) -> dict:
    author = authorship.get("author") or {}
    name = authorship.get("raw_author_name") or author.get("display_name") or "Author unavailable"
    affiliations = unique_strings(authorship.get("raw_affiliation_strings") or [])
    if not affiliations:
        affiliations = unique_strings(
            institution.get("display_name", "") for institution in authorship.get("institutions") or []
        )
    return {"name": name, "affiliations": affiliations}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("links_path", type=Path)
    parser.add_argument("output_path", type=Path)
    args = parser.parse_args()

    papers = {}
    with args.csv_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            index = row["Index"].strip()
            if index and index not in papers:
                papers[index] = {
                    "title": row["Paper Title"].strip(),
                    "year": int(row["Year"].strip()),
                }

    links = load_javascript_object(args.links_path)
    index_to_doi = {index: doi_from_url(item.get("url", "")) for index, item in links.items()}
    missing_dois = [index for index in papers if not index_to_doi.get(index)]
    if missing_dois:
        raise SystemExit(f"Missing DOI links for benchmark indices: {', '.join(missing_dois)}")

    dois = sorted(set(index_to_doi.values()))
    works = fetch_openalex(dois)
    metadata = {}
    low_matches = []
    for index in sorted(papers, key=int):
        paper = papers[index]
        doi = index_to_doi[index]
        work = works.get(doi)
        if not work:
            metadata[index] = {
                "doi": doi,
                "authors": [],
                "status": "metadata-unavailable",
            }
            continue

        score = title_score(paper["title"], work.get("title", ""))
        if score < 0.79:
            low_matches.append((index, score, paper["title"], work.get("title", "")))
        authors = [author_record(item) for item in work.get("authorships") or []]
        metadata[index] = {
            "doi": doi,
            "openAlexUrl": work.get("id", ""),
            "matchedTitle": work.get("title", ""),
            "titleScore": round(score, 3),
            "authors": authors,
            "status": "complete" if authors else "metadata-unavailable",
        }

    output = {
        "source": "OpenAlex DOI metadata",
        "retrieved": date.today().isoformat(),
        "papers": metadata,
    }
    payload = "// Generated by scripts/build-paper-metadata.py from DOI-indexed OpenAlex metadata.\n"
    payload += "window.BENCHMARK_PAPER_METADATA = " + json.dumps(output, indent=2, ensure_ascii=False) + ";\n"
    args.output_path.write_text(payload, encoding="utf-8")

    complete = sum(item["status"] == "complete" for item in metadata.values())
    authors = sum(len(item["authors"]) for item in metadata.values())
    affiliated = sum(
        bool(author["affiliations"])
        for item in metadata.values()
        for author in item["authors"]
    )
    print(f"Wrote {len(metadata)} papers: {complete} with bylines, {authors} authors, {affiliated} with affiliation data.")
    if low_matches:
        print("Low title matches requiring review:")
        for index, score, expected, matched in low_matches:
            print(f"  #{index} ({score:.3f})\n    CSV: {expected}\n    API: {matched}")


if __name__ == "__main__":
    main()
