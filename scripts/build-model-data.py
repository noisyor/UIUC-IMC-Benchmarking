#!/usr/bin/env python3
"""Build the browser-side eNVM model table from the published MIT-licensed data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def rounded(values: np.ndarray) -> list:
    return np.vectorize(lambda value: float(f"{value:.8g}"))(values).tolist()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Clone of calmyor/eNVM-IMC-Modeling")
    parser.add_argument("output", type=Path, help="Generated JavaScript output")
    args = parser.parse_args()

    model_dir = args.source / "SNDR-sim" / "SNDRd-vs-Energy"
    dimensions = np.logspace(1.5, 4, 30) * 0.5
    devices = {}
    for device in ("ReRAM", "MRAM", "FeFET"):
        vbl = np.logspace(-3, 0.7 if device == "FeFET" else -0.09, 40)
        devices[device] = {
            "vdd": 5.0 if device == "FeFET" else 0.8,
            "vbl": rounded(vbl),
            "sndrDb": rounded(np.load(model_dir / f"SNDRd_dB_{device}_vs_N.npy")),
            "averageGeq": rounded(np.load(model_dir / f"AGeq_{device}_vs_N.npy")),
        }

    payload = {
        "source": "https://github.com/calmyor/eNVM-IMC-Modeling",
        "sourceModel": "Energy-Accuracy Trade-offs for Resistive In-Memory Computing Architectures",
        "license": "MIT",
        "modelScope": "Published six-bit-ADC SNDR and energy surface",
        "dimensions": rounded(dimensions),
        "devices": devices,
        "energyConstants": {
            "adcLinearJPerBit": 100e-15,
            "adcExponentialJ": 1e-18,
            "sourceCorePeriodSeconds": 20e-9,
        },
    }

    notice = """// Generated from calmyor/eNVM-IMC-Modeling (MIT License).
// Copyright (c) 2023 Saion Kumar Roy. See the source repository for the full license.
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        notice + "window.IMC_MODEL_DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
