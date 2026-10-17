#!/usr/bin/env python3
"""Create a long-conditioning cross-model gradient-reversal campaign."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


MODELS = ("wave_pinning", "otsuji", "goryachev", "legi")
LENGTHS = tuple(float(value) for value in range(1, 40, 2))


def build_manifest() -> dict:
    tasks = []
    for model in MODELS:
        for length in LENGTHS:
            conditioning = max(2000.0, 400.0 * (length / 10.0) ** 2)
            post_reversal = max(4000.0, 800.0 * (length / 10.0) ** 2)
            tasks.append(
                {
                    "tier": "revision_audit", "model": model,
                    "assay": "gradient_reversal_long_conditioning",
                    "protocol": {
                        "name": "reversing_gradient", "slope": 0.1 / length,
                        "reversal_time": conditioning, "transition": 20.0,
                    },
                    "simulation": {
                        "length": length, "cells": max(101, round(length / 0.1)),
                        "end_time": conditioning + post_reversal,
                        "output_interval": max(2.0, (conditioning + post_reversal) / 2400.0),
                        "backend": "bdf", "max_step": 10.0,
                    },
                    "metrics": {
                        "analysis_start": 0.0, "polarization_threshold": 0.1,
                        "orientation_threshold": 0.05,
                        "sustained_duration": max(20.0, conditioning / 100.0),
                        "reversal_time": conditioning,
                    },
                    "replicate": 0,
                    "save_trajectory": model == "otsuji" or length in {5.0, 21.0, 39.0},
                }
            )
    return make_manifest(
        "jilkine_reversal_long_conditioning_v1", tasks,
        description=(
            "Cross-model reversal after a substantially longer conditioning phase, "
            "with explicit pre-reversal morphology and censoring."
        ),
        metadata={
            "conditioning_policy": "max(2000, 400*(L/10)^2) s",
            "post_reversal_policy": "max(4000, 800*(L/10)^2) s",
            "gradient_peak": 0.1,
            "purpose": "distinguish incomplete conditioning from reversal failure",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "jilkine_reversal_long_conditioning_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
