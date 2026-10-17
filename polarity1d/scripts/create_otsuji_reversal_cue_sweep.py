#!/usr/bin/env python3
"""Create an Otsuji-only reversal cue amplitude/conditioning-time sweep."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


LENGTH = 10.0
INPUT_PEAKS = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
CONDITIONING_TIMES = (100.0, 200.0, 500.0, 1000.0, 2000.0, 4000.0)
POST_REVERSAL_TIME = 4000.0


def build_manifest() -> dict:
    tasks = []
    for input_peak in INPUT_PEAKS:
        for conditioning in CONDITIONING_TIMES:
            end_time = conditioning + POST_REVERSAL_TIME
            tasks.append(
                {
                    "tier": "targeted_audit",
                    "model": "otsuji",
                    "assay": "gradient_reversal_cue_sweep",
                    "input_peak": input_peak,
                    "conditioning_duration_s": conditioning,
                    "protocol": {
                        "name": "reversing_gradient",
                        "slope": input_peak / LENGTH,
                        "reversal_time": conditioning,
                        "transition": 20.0,
                    },
                    "simulation": {
                        "length": LENGTH,
                        "cells": 101,
                        "end_time": end_time,
                        "output_interval": max(2.0, end_time / 2400.0),
                        "backend": "bdf",
                        "max_step": 10.0,
                    },
                    "metrics": {
                        "analysis_start": 0.0,
                        "polarization_threshold": 0.1,
                        "orientation_threshold": 0.05,
                        "sustained_duration": 20.0,
                        "reversal_time": conditioning,
                    },
                    "replicate": 0,
                    "save_trajectory": True,
                }
            )
    return make_manifest(
        "otsuji_reversal_cue_sweep_l10_v1",
        tasks,
        description=(
            "Otsuji reversal at L=10 um across gradient peak strength and "
            "pre-reversal conditioning duration, with a 4000-s follow-up."
        ),
        metadata={
            "length_um": LENGTH,
            "input_peaks_per_s": list(INPUT_PEAKS),
            "conditioning_times_s": list(CONDITIONING_TIMES),
            "post_reversal_time_s": POST_REVERSAL_TIME,
            "transition_time_s": 20.0,
            "primary_endpoint": "sustained opposite-side orientation",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=ROOT / "campaigns" / "otsuji_reversal_cue_sweep_l10_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
