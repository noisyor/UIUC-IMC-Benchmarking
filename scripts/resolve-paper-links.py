#!/usr/bin/env python3
"""Resolve benchmark titles to DOI links, with exact-title search fallbacks."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import unicodedata
from difflib import SequenceMatcher
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


def normalized(text: str) -> str:
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    text = text.lower().replace("–", "-").replace("—", "-")
    return " ".join(re.findall(r"[a-z0-9]+", text))


def title_score(expected: str, candidate: str) -> float:
    left = normalized(expected)
    right = normalized(candidate)
    sequence = SequenceMatcher(None, left, right).ratio()
    left_tokens = set(left.split())
    right_tokens = set(right.split())
    overlap = len(left_tokens & right_tokens) / max(1, len(left_tokens | right_tokens))
    return 0.62 * sequence + 0.38 * overlap


def crossref_candidates(title: str, year: int) -> list[dict]:
    params = urlencode({
        "query.title": title,
        "filter": f"from-pub-date:{year - 1}-01-01,until-pub-date:{year + 1}-12-31",
        "rows": 5,
        "select": "DOI,title,published",
    })
    request = Request(
        f"https://api.crossref.org/works?{params}",
        headers={"User-Agent": "UIUC-IMC-Benchmarking-link-resolver/1.0"},
    )
    for attempt in range(4):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)["message"]["items"]
        except HTTPError as error:
            if error.code != 429 or attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    return []


def resolve(title: str, year: int) -> dict:
    best = None
    for item in crossref_candidates(title, year):
        candidate_title = (item.get("title") or [""])[0]
        score = title_score(title, candidate_title)
        if best is None or score > best[0]:
            best = (score, item, candidate_title)
    if best and best[0] >= 0.79 and best[1].get("DOI"):
        return {
            "url": f"https://doi.org/{best[1]['DOI']}",
            "kind": "doi",
            "score": round(best[0], 3),
            "matchedTitle": best[2],
        }
    return {
        "url": f"https://scholar.google.com/scholar?q={quote(chr(34) + title + chr(34))}",
        "kind": "search",
        "score": round(best[0], 3) if best else 0,
        "matchedTitle": best[2] if best else "",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
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

    resolved = {}
    if args.output_path.exists():
        existing_text = args.output_path.read_text(encoding="utf-8")
        resolved = json.loads(existing_text[existing_text.index("{"):existing_text.rindex("}") + 1])
    for position, (index, paper) in enumerate(sorted(papers.items(), key=lambda item: int(item[0])), 1):
        if resolved.get(index, {}).get("kind") == "doi":
            print(f"[{position:03}/{len(papers)}] #{index}: cached DOI", flush=True)
            continue
        try:
            result = resolve(paper["title"], paper["year"])
        except Exception as error:
            result = {
                "url": f"https://scholar.google.com/scholar?q={quote(chr(34) + paper['title'] + chr(34))}",
                "kind": "search",
                "score": 0,
                "matchedTitle": "",
                "error": str(error),
            }
        resolved[index] = result
        print(f"[{position:03}/{len(papers)}] #{index}: {result['kind']} {result['score']}", flush=True)
        time.sleep(0.55)

    payload = "// Generated from Benchmarking_Data.csv via Crossref title matching.\n"
    payload += "window.BENCHMARK_PAPER_LINKS = " + json.dumps(resolved, indent=2, ensure_ascii=False) + ";\n"
    args.output_path.write_text(payload, encoding="utf-8")
    doi_count = sum(item["kind"] == "doi" for item in resolved.values())
    print(f"Wrote {len(resolved)} links: {doi_count} DOI, {len(resolved) - doi_count} title-search fallbacks.")


if __name__ == "__main__":
    main()
