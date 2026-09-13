#!/usr/bin/env python3
"""Fit the builder's data-derived column-energy models to the benchmark CSV.

Analog IMC (one ADC column per invocation, as in the CSV's E_col):
    E_col = k1 * B_ADC * (node / 28 nm)^2 + a_family * N * VDD^2 * (node / 28 nm)^q
Digital IMC:
    E_col = a * (N * B_x * B_w)^b * VDD^2 * (node / 28 nm)^c

The ADC node exponent is fixed at 2: freeing it fits the data marginally better but
predicts implausible ADC energies outside the 28-65 nm bulk of the dataset. The 4^B ADC
term is not identifiable from these rows (the fit drives it to zero), so it is omitted.
Every family reports its leave-one-paper-out error so the builder can show how far a
new design is likely to sit from the fit.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ANALOG_FAMILIES = ("QR", "QS", "SRAM IS", "resistive IS")
NODE_REFERENCE_NM = 28.0
ADC_NODE_EXPONENT = 2.0
EXCLUDED = {"136": "E_col is about 29,000x larger than the row's own E_OP1 implies"}


def number(value: str | None) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def family(architecture: str, model: str) -> str | None:
    if model == "DIMC":
        return "DIMC"
    if architecture in ("eNVM", "eFlash") and model == "IS":
        return "resistive IS"
    return {"IS": "SRAM IS", "QR": "QR", "QS": "QS", "QS-QR": "QS"}.get(model)


def load_rows(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            index = (raw.get("Index") or "").strip()
            if not index or index in EXCLUDED:
                continue
            row = {
                "index": index,
                "family": family(raw["Architecture"].strip(), raw["Compute Model"].strip()),
                "adc": number(raw["B_ADC"]),
                "energy": number(raw["E_col (fJ)"]),
                "N": number(raw["N"]),
                "bx": number(raw["B_x"]),
                "bw": number(raw["B_w"]),
                "node": number(raw["Tech (nm)"]),
                "vdd": number(raw["Supply V(V)"]),
            }
            if not row["family"] or not all(row[key] for key in ("energy", "N", "node", "vdd")):
                continue
            if row["family"] == "DIMC" and row["bx"] and row["bw"]:
                rows.append(row)
            elif row["family"] != "DIMC" and row["adc"] and row["adc"] > 0:
                rows.append(row)
    return rows


def analog_energy(params: np.ndarray, row: dict) -> float:
    scale = row["node"] / NODE_REFERENCE_NM
    adc = 10 ** params[0] * row["adc"] * scale ** ADC_NODE_EXPONENT
    array = 10 ** params[1 + ANALOG_FAMILIES.index(row["family"])] * row["N"] * row["vdd"] ** 2 * scale ** params[5]
    return adc + array


def digital_energy(params: np.ndarray, row: dict) -> float:
    scale = row["node"] / NODE_REFERENCE_NM
    return 10 ** params[0] * (row["N"] * row["bx"] * row["bw"]) ** params[1] * row["vdd"] ** 2 * scale ** params[2]


def fit(model, rows: list[dict], start: list[float]) -> np.ndarray:
    measured = np.log10([row["energy"] for row in rows])
    return least_squares(lambda p: np.log10([max(model(p, row), 1e-12) for row in rows]) - measured, start).x


def leave_one_paper_out(model, rows: list[dict], start: list[float]) -> list[float]:
    errors = [0.0] * len(rows)
    for paper in sorted({row["index"] for row in rows}):
        params = fit(model, [row for row in rows if row["index"] != paper], start)
        for position, row in enumerate(rows):
            if row["index"] == paper:
                errors[position] = math.log10(max(model(params, row), 1e-12)) - math.log10(row["energy"])
    return errors


def summarize(rows: list[dict], errors: list[float]) -> dict:
    absolute = np.abs(errors)
    return {
        "rows": len(rows),
        "papers": len({row["index"] for row in rows}),
        "typicalFactor": round(float(10 ** np.median(absolute)), 3),
        "p80Factor": round(float(10 ** np.percentile(absolute, 80)), 3),
        "within2x": round(float(np.mean(absolute <= math.log10(2))), 3),
        "medianBias": round(float(10 ** np.median(errors)), 3),
    }


def rounded(value: float) -> float:
    return float(f"{value:.6g}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path, help="Benchmarking_Data.csv")
    parser.add_argument("output", type=Path, help="Generated JavaScript output")
    args = parser.parse_args()

    rows = load_rows(args.csv)
    analog_rows = [row for row in rows if row["family"] in ANALOG_FAMILIES]
    digital_rows = [row for row in rows if row["family"] == "DIMC"]
    analog_start = [0.8, 1.0, 1.0, 1.0, 1.0, 0.0]
    digital_start = [0.0, 1.0, 1.0]

    analog_params = fit(analog_energy, analog_rows, analog_start)
    digital_params = fit(digital_energy, digital_rows, digital_start)
    analog_errors = leave_one_paper_out(analog_energy, analog_rows, analog_start)
    digital_errors = leave_one_paper_out(digital_energy, digital_rows, digital_start)

    validation = {
        name: summarize(
            [row for row in analog_rows if row["family"] == name],
            [error for error, row in zip(analog_errors, analog_rows) if row["family"] == name],
        )
        for name in ANALOG_FAMILIES
    }
    validation["analog"] = summarize(analog_rows, analog_errors)
    validation["DIMC"] = summarize(digital_rows, digital_errors)

    payload = {
        "source": "Benchmarking_Data.csv",
        "generator": "scripts/build-energy-fit.py",
        "nodeReferenceNm": NODE_REFERENCE_NM,
        "analog": {
            "form": "E_col = k1·B_ADC·(node/28 nm)² + a_family·N·VDD²·(node/28 nm)^q",
            "k1FjPerBit": rounded(10 ** analog_params[0]),
            "adcNodeExponent": ADC_NODE_EXPONENT,
            "arrayFjPerRowVolt2": {name: rounded(10 ** analog_params[1 + i]) for i, name in enumerate(ANALOG_FAMILIES)},
            "arrayNodeExponent": rounded(analog_params[5]),
        },
        "dimc": {
            "form": "E_col = a·(N·B_x·B_w)^b·VDD²·(node/28 nm)^c",
            "aFj": rounded(10 ** digital_params[0]),
            "opsExponent": rounded(digital_params[1]),
            "nodeExponent": rounded(digital_params[2]),
        },
        "validation": {"method": "leave-one-paper-out on E_col", "families": validation},
        "excluded": EXCLUDED,
    }

    notice = """// Generated by scripts/build-energy-fit.py from Benchmarking_Data.csv.
// Data-derived column-energy fits with leave-one-paper-out validation.
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(notice + "window.IMC_ENERGY_FIT = " + json.dumps(payload, separators=(",", ":")) + ";\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
