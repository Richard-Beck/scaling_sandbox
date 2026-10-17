#!/usr/bin/env python3
"""Create a De Belly reversal/tracking input-strength titration manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_input_titration_v1"
INPUT_PEAKS = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
LENGTHS = (12.0, 16.0, 20.0, 2 ** (1 / 3) * 20.0, 32.0, 40.0, 50.0, 64.0)
TRACKING_LENGTHS = (20.0, 2 ** (1 / 3) * 20.0, 40.0, 64.0)
TRACKING_PERIODS = (100.0, 200.0, 400.0)
VARIANTS = {
    "intact": {},
    "rac_to_rho_half": {"rac_to_rho_coupling": 0.5},
    "rho_to_rac_half": {"rho_to_rac_coupling": 0.5},
    "rac_to_rho_off": {"rac_to_rho_coupling": 0.0},
    "rho_to_rac_off": {"rho_to_rac_coupling": 0.0},
    "both_off": {"rac_to_rho_coupling": 0.0, "rho_to_rac_coupling": 0.0},
}


def _cells(length: float) -> int:
    return max(60, round(length / 0.4))


def _base_task(assay: str, length: float, variant: str, input_peak: float) -> dict:
    return {
        "tier": "titration",
        "model": "debelly",
        "assay": assay,
        "coupling_variant": variant,
        "model_parameters": VARIANTS[variant],
        "input_peak": input_peak,
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


def titration_tasks() -> list[dict]:
    tasks = []
    tracking_saved_peaks = {0.05, 0.1, 0.2}
    for input_peak in INPUT_PEAKS:
        for variant in VARIANTS:
            for length in LENGTHS:
                reversal_time = 250.0
                task = _base_task(
                    "gradient_reversal_titration", length, variant, input_peak
                )
                task["protocol"] = {
                    "name": "reversing_gradient",
                    "slope": input_peak / length,
                    "reversal_time": reversal_time,
                    "transition": 25.0,
                }
                task["metrics"]["reversal_time"] = reversal_time
                task["save_trajectory"] = input_peak in {0.2, 0.5, 1.0} and (
                    (
                        variant == "intact"
                        and length
                        in {20.0, 2 ** (1 / 3) * 20.0, 32.0, 40.0, 64.0}
                    )
                    or (
                        variant == "rac_to_rho_half"
                        and length in {2 ** (1 / 3) * 20.0, 32.0, 64.0}
                    )
                    or (
                        variant in {"rho_to_rac_half", "both_off"}
                        and length in {32.0, 64.0}
                    )
                    or (
                        input_peak == 0.2
                        and variant in {"rac_to_rho_off", "rho_to_rac_off"}
                        and length == 32.0
                    )
                )
                tasks.append(task)

    for input_peak in INPUT_PEAKS:
        for variant in VARIANTS:
            for length in TRACKING_LENGTHS:
                for period in TRACKING_PERIODS:
                    task = _base_task(
                        "dynamic_tracking_titration", length, variant, input_peak
                    )
                    task["protocol"] = {
                        "name": "alternating_gradient",
                        "slope": input_peak / length,
                        "period": period,
                        "start": 0.0,
                    }
                    task["simulation"].update(
                        {
                            "end_time": 6.0 * period,
                            "output_interval": period / 20.0,
                        }
                    )
                    task["metrics"]["sustained_duration"] = period / 4.0
                    task["tracking"] = {
                        "period": period,
                        "start": period,
                        "minimum_amplitude": 0.05,
                        "maximum_lag_fraction": 0.25,
                    }
                    task["save_trajectory"] = (
                        length == 20.0
                        and period in {100.0, 400.0}
                        and input_peak in tracking_saved_peaks
                        and variant in {"intact", "both_off"}
                    )
                    tasks.append(task)
    return tasks


def build_manifest() -> dict:
    return make_manifest(
        CAMPAIGN_NAME,
        titration_tasks(),
        description=(
            "Deterministic De Belly reversal and alternating-gradient tracking "
            "titration across mechanochemical arm perturbations."
        ),
        metadata={
            "input_peak_definition": (
                "maximum Rac activation added at the high-input boundary; basal "
                "Rac activation is 1 and maximum mechanical increment is 0.2"
            ),
            "input_peaks": list(INPUT_PEAKS),
            "reversal_policy": (
                "same 250-s conditioning, 25-s cross-fade, and 500-s follow-up as "
                "the Phase 1 production assay"
            ),
            "tracking_policy": (
                "same lengths, periods, six cycles, and first-cycle exclusion as "
                "the Phase 1 production assay"
            ),
            "replication": "one deterministic run per condition",
            "purpose": (
                "locate input strengths at which reciprocal mechanochemical-arm "
                "perturbations become distinguishable"
            ),
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
