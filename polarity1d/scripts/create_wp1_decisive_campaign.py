#!/usr/bin/env python3
"""Create the decisive WP1 size-by-mechanochemical-coupling campaign."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "wp1_decisive_v1"
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
REVERSAL_INPUTS = (0.2, 0.5)
FIXED_SIGMA_UM = 1.6
VARIANTS = {
    "intact": {},
    "rac_to_rho_half": {"rac_to_rho_coupling": 0.5},
    "rho_to_rac_half": {"rho_to_rac_coupling": 0.5},
    "rac_to_rho_off": {"rac_to_rho_coupling": 0.0},
    "rho_to_rac_off": {"rho_to_rac_coupling": 0.0},
    "both_off": {"rac_to_rho_coupling": 0.0, "rho_to_rac_coupling": 0.0},
}
SAVED_LENGTHS = {20.0, 2 ** (1 / 3) * 20.0, 40.0, 64.0}


def _cells(length: float) -> int:
    return max(60, round(length / 0.4))


def _base_task(assay: str, length: float, variant: str) -> dict:
    return {
        "tier": "decisive_wp1",
        "model": "debelly",
        "assay": assay,
        "coupling_variant": variant,
        "model_parameters": VARIANTS[variant],
        "protocol": {"name": "none"},
        "simulation": {
            "length": length,
            "cells": _cells(length),
            "end_time": 1200.0,
            "output_interval": 2.0,
            "step": 0.5,
        },
        "metrics": {
            "analysis_start": 0.0,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 50.0,
        },
        "mechanochemical_metrics": {"stimulus_start": 0.0},
        "replicate": 0,
        "save_trajectory": length in SAVED_LENGTHS,
    }


def _competition_end_time(length: float) -> float:
    """At least double the previous intact L=64 observation window."""
    return max(2400.0, 8000.0 * (length / 64.0) ** 2)


def _conditioning_time(length: float) -> float:
    return max(250.0, 250.0 * (length / 20.0) ** 2)


def _post_reversal_time(length: float) -> float:
    return max(750.0, 750.0 * (length / 20.0) ** 2)


def campaign_tasks() -> list[dict]:
    tasks = []
    for variant in VARIANTS:
        for length in LENGTHS:
            task = _base_task("dual_site_fixed_physical_width_extended", length, variant)
            task["protocol"] = {
                "name": "gaussian_pulse",
                "amplitude": 1.0,
                "start": 50.0,
                "stop": 150.0,
                "width_fraction": FIXED_SIGMA_UM / length,
                "centers": [0.0, 1.0],
                "asymmetry": 0.02,
            }
            task["simulation"]["end_time"] = _competition_end_time(length)
            task["metrics"]["analysis_start"] = 50.0
            task["mechanochemical_metrics"] = {"stimulus_start": 50.0}
            tasks.append(task)

    for input_peak in REVERSAL_INPUTS:
        for variant in VARIANTS:
            for length in LENGTHS:
                conditioning = _conditioning_time(length)
                task = _base_task("gradient_reversal_length_scaled", length, variant)
                task["input_peak"] = input_peak
                task["protocol"] = {
                    "name": "reversing_gradient",
                    "slope": input_peak / length,
                    "reversal_time": conditioning,
                    "transition": 25.0,
                }
                task["simulation"]["end_time"] = conditioning + _post_reversal_time(length)
                task["metrics"].update(
                    {"sustained_duration": 25.0, "reversal_time": conditioning}
                )
                tasks.append(task)

    for variant in VARIANTS:
        for length in LENGTHS:
            task = _base_task("single_site_fixed_physical_width", length, variant)
            task["protocol"] = {
                "name": "gaussian_pulse",
                "amplitude": 1.0,
                "start": 50.0,
                "stop": 150.0,
                "width_fraction": FIXED_SIGMA_UM / length,
                "centers": [0.0],
            }
            task["simulation"]["end_time"] = _competition_end_time(length)
            task["metrics"]["analysis_start"] = 50.0
            task["mechanochemical_metrics"] = {"stimulus_start": 50.0}
            tasks.append(task)
    return tasks


def build_manifest() -> dict:
    return make_manifest(
        CAMPAIGN_NAME,
        campaign_tasks(),
        description=(
            "Decisive WP1 campaign: long-time fixed-width two-front competition, "
            "length-scaled intermediate-input reversal, and fixed-width distal "
            "propagation across reciprocal mechanochemical arm perturbations."
        ),
        metadata={
            "steps": [
                "extended fixed-width competition across all coupling states",
                "persistent-versus-transient resolution endpoint",
                "diffusive rescaling analysis",
                "length-scaled reversal at inputs 0.2 and 0.5",
                "fixed-width single-site propagation control",
            ],
            "lengths_um": list(LENGTHS),
            "coupling_variants": list(VARIANTS),
            "fixed_gaussian_sigma_um": FIXED_SIGMA_UM,
            "competition_end_time_policy": "max(2400, 8000*(L/64)^2) s",
            "propagation_end_time_policy": "max(2400, 8000*(L/64)^2) s",
            "reversal_inputs": list(REVERSAL_INPUTS),
            "reversal_conditioning_policy": "max(250, 250*(L/20)^2) s",
            "reversal_follow_up_policy": "max(750, 750*(L/20)^2) s",
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
