#!/usr/bin/env python3
"""Create the intact-first De Belly length-effects campaign manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_intact_length_v2"
LENGTHS = (
    12.0,
    16.0,
    20.0,
    2 ** (1 / 3) * 20.0,
    32.0,
    36.0,
    40.0,
    44.0,
    48.0,
    50.0,
    52.0,
    56.0,
    60.0,
    64.0,
)
REVERSAL_INPUTS = (0.2, 0.5, 1.0)
FIXED_GAUSSIAN_SIGMA_UM = 1.6


def _cells(length: float) -> int:
    return max(60, round(length / 0.4))


def _base_task(assay: str, length: float) -> dict:
    return {
        "tier": "intact_length",
        "model": "debelly",
        "assay": assay,
        "coupling_variant": "intact",
        "model_parameters": {},
        "protocol": {"name": "none"},
        "simulation": {
            "length": length,
            "cells": _cells(length),
            "end_time": 750.0,
            "output_interval": 2.0,
            "step": 0.5,
        },
        "metrics": {
            "analysis_start": 0.0,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 25.0,
        },
        "mechanochemical_metrics": {"stimulus_start": 0.0},
        "replicate": 0,
        "save_trajectory": False,
    }


def _competition_end_time(length: float) -> float:
    """Give L=50 at least 2,000 s and L=64 more than 3,000 s."""
    return max(1200.0, 2000.0 * (length / 50.0) ** 2)


def _conditioning_time(length: float) -> float:
    """Length-aware conditioning anchored to the original 250 s at L=20."""
    return max(250.0, 250.0 * (length / 20.0) ** 2)


def _post_reversal_time(length: float) -> float:
    """Length-aware follow-up anchored to the original 500 s at L=20."""
    return max(500.0, 500.0 * (length / 20.0) ** 2)


def campaign_tasks() -> list[dict]:
    tasks = []
    saved_lengths = {20.0, 40.0, 50.0, 64.0}

    for width_mode in ("fractional", "fixed_physical"):
        assay = f"dual_site_{width_mode}_width"
        for length in LENGTHS:
            width_fraction = (
                0.08
                if width_mode == "fractional"
                else FIXED_GAUSSIAN_SIGMA_UM / length
            )
            task = _base_task(assay, length)
            task["protocol"] = {
                "name": "gaussian_pulse",
                "amplitude": 1.0,
                "start": 50.0,
                "stop": 150.0,
                "width_fraction": width_fraction,
                "centers": [0.0, 1.0],
                "asymmetry": 0.02,
            }
            task["simulation"].update(
                {
                    "end_time": _competition_end_time(length),
                    "output_interval": 2.0,
                }
            )
            task["metrics"].update(
                {"analysis_start": 50.0, "sustained_duration": 50.0}
            )
            task["mechanochemical_metrics"] = {"stimulus_start": 50.0}
            task["save_trajectory"] = length in saved_lengths
            tasks.append(task)

    for input_peak in REVERSAL_INPUTS:
        for length in LENGTHS:
            conditioning = _conditioning_time(length)
            follow_up = _post_reversal_time(length)
            task = _base_task("gradient_reversal_length_scaled", length)
            task["input_peak"] = input_peak
            task["protocol"] = {
                "name": "reversing_gradient",
                "slope": input_peak / length,
                "reversal_time": conditioning,
                "transition": 25.0,
            }
            task["simulation"].update(
                {
                    "end_time": conditioning + follow_up,
                    "output_interval": 2.0,
                }
            )
            task["metrics"].update(
                {
                    "analysis_start": 0.0,
                    "sustained_duration": 25.0,
                    "reversal_time": conditioning,
                }
            )
            task["save_trajectory"] = length in saved_lengths
            tasks.append(task)
    return tasks


def build_manifest() -> dict:
    return make_manifest(
        CAMPAIGN_NAME,
        campaign_tasks(),
        description=(
            "Intact De Belly length-effects campaign with extended two-front "
            "follow-up, a dense large-cell grid, paired fractional/fixed-physical "
            "two-site cues, and length-aware reversal conditioning and follow-up."
        ),
        metadata={
            "focus": "cell-length effects in the intact De Belly model",
            "lengths_um": list(LENGTHS),
            "dual_site_width_modes": {
                "fractional": "Gaussian sigma = 0.08 L",
                "fixed_physical": f"Gaussian sigma = {FIXED_GAUSSIAN_SIGMA_UM} um",
            },
            "competition_end_time_policy": "max(1200, 2000*(L/50)^2) s",
            "reversal_inputs": list(REVERSAL_INPUTS),
            "reversal_conditioning_policy": "max(250, 250*(L/20)^2) s",
            "reversal_follow_up_policy": "max(500, 500*(L/20)^2) s",
            "reversal_endpoint": "sustained Rac side-orientation switch",
            "replication": "one deterministic run per condition",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
