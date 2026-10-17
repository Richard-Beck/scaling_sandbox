#!/usr/bin/env python3
"""Create the focused Otsuji single-cue length/amplitude qualification."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


LENGTHS = tuple(float(value) for value in range(1, 40, 2))
AMPLITUDES = (0.04, 0.2)


def build_manifest() -> dict:
    tasks = []
    for amplitude in AMPLITUDES:
        for length in LENGTHS:
            end_time = max(2000.0, 500.0 * (length / 10.0) ** 2)
            tasks.append(
                {
                    "tier": "qualification",
                    "model": "otsuji",
                    "assay": "single_cue_polarization",
                    "protocol": {
                        "name": "localized_pulse",
                        "amplitude": amplitude,
                        "start": 0.0,
                        "stop": 40.0,
                        "taper": 20.0,
                        "width_fraction": 0.1,
                        "sides": ["left"],
                        "asymmetry": 0.0,
                    },
                    "simulation": {
                        "length": length,
                        "cells": max(101, round(length / 0.1)),
                        "end_time": end_time,
                        "output_interval": max(2.0, end_time / 800.0),
                        "backend": "bdf",
                        "max_step": 10.0,
                    },
                    "metrics": {
                        "analysis_start": 40.0,
                        "polarization_threshold": 0.1,
                        "orientation_threshold": 0.05,
                        "sustained_duration": 50.0,
                    },
                    "replicate": 0,
                    "seed": None,
                    "input_peak": amplitude,
                    "save_trajectory": length in {5.0, 13.0, 21.0, 31.0, 39.0},
                }
            )
    return make_manifest(
        "otsuji_single_cue_v1",
        tasks,
        description=(
            "Focused test of whether Otsuji's apparent large-cell failure under two "
            "competing cues persists for a single transient cue."
        ),
        metadata={
            "cue_amplitudes": {
                "figure_guided": 0.04,
                "production_common_cue": 0.2,
            },
            "cue_geometry": "leftmost 10% of the domain",
            "timing": "S/2 through 20 s, cosine taper to zero at 40 s",
            "initial_condition": "homogeneous Otsuji steady state; no stochastic perturbation",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=ROOT / "campaigns" / "otsuji_single_cue_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
