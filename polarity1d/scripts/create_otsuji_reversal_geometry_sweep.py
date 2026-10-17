#!/usr/bin/env python3
"""Create an Otsuji reversal sweep over stimulus amplitude and spatial width."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


LENGTH = 15.0
AMPLITUDES = (0.01, 0.02, 0.04, 0.08, 0.12, 0.2, 0.35, 0.5)
WIDTH_FRACTIONS = (0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.75, 1.0)
REVERSAL_TIME = 2_000.0
FOLLOWUP = 20_000.0


def build_manifest() -> dict:
    tasks = []
    for amplitude in AMPLITUDES:
        for width in WIDTH_FRACTIONS:
            tasks.append(
                {
                    "tier": "targeted_audit", "model": "otsuji",
                    "assay": "otsuji_reversal_geometry_sweep",
                    "input_peak": amplitude, "stimulus_width_fraction": width,
                    "protocol": {
                        "name": "reversing_finite_gradient", "amplitude": amplitude,
                        "width_fraction": width, "reversal_time": REVERSAL_TIME,
                        "transition": 20.0,
                    },
                    "simulation": {
                        "length": LENGTH, "cells": 150,
                        "end_time": REVERSAL_TIME + FOLLOWUP,
                        "output_interval": 20.0, "backend": "bdf", "max_step": 10.0,
                    },
                    "metrics": {
                        "analysis_start": 0.0, "polarization_threshold": 0.1,
                        "orientation_threshold": 0.05, "sustained_duration": 20.0,
                        "reversal_time": REVERSAL_TIME,
                    },
                    "replicate": 0, "save_trajectory": True,
                }
            )
    return make_manifest(
        "otsuji_reversal_geometry_sweep_l15_v1", tasks,
        description=(
            "Otsuji reversal at L=15 um across independently varied peak amplitude "
            "and spatial gradient width, with 20,000 s post-reversal follow-up."
        ),
        metadata={
            "length_um": LENGTH, "amplitudes_per_s": list(AMPLITUDES),
            "width_fractions": list(WIDTH_FRACTIONS),
            "conditioning_time_s": REVERSAL_TIME, "followup_time_s": FOLLOWUP,
            "transition_time_s": 20.0,
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "otsuji_reversal_geometry_sweep_l15_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
