#!/usr/bin/env python3
"""Create the Phase 1 multi-assay production campaign manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


JILKINE_MODELS = ("wave_pinning", "otsuji", "goryachev", "legi")
JILKINE_LENGTHS = tuple(float(value) for value in range(1, 40, 2))
TRACKING_LENGTHS = (5.0, 15.0, 25.0, 39.0)
TRACKING_PERIODS = (50.0, 100.0, 200.0, 400.0)
DEBELLY_LENGTHS = (12.0, 16.0, 20.0, 2 ** (1 / 3) * 20.0, 32.0, 40.0, 50.0, 64.0)
DEBELLY_TRACKING_LENGTHS = (20.0, 2 ** (1 / 3) * 20.0, 40.0, 64.0)
DEBELLY_VARIANTS = {
    "intact": {},
    "rac_to_rho_half": {"rac_to_rho_coupling": 0.5},
    "rho_to_rac_half": {"rho_to_rac_coupling": 0.5},
    "rac_to_rho_off": {"rac_to_rho_coupling": 0.0},
    "rho_to_rac_off": {"rho_to_rac_coupling": 0.0},
    "both_off": {"rac_to_rho_coupling": 0.0, "rho_to_rac_coupling": 0.0},
}


def _jilkine_cells(length: float, multiplier: int = 1) -> int:
    return multiplier * max(101, round(length / 0.1))


def _debelly_cells(length: float, multiplier: int = 1) -> int:
    return multiplier * max(60, round(length / 0.4))


def _base_task(model: str, assay: str, length: float, cells: int) -> dict:
    return {
        "tier": "production",
        "model": model,
        "assay": assay,
        "protocol": {"name": "none"},
        "simulation": {
            "length": length,
            "cells": cells,
            "end_time": 300.0,
            "output_interval": 2.0,
        },
        "metrics": {
            "analysis_start": 0.0,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 20.0,
        },
        "replicate": 0,
        "save_trajectory": False,
    }


def jilkine_tasks() -> list[dict]:
    tasks = []
    selected_trajectories = {5.0, 21.0, 39.0}
    for model in JILKINE_MODELS:
        replicates = range(5) if model != "legi" else range(1)
        for length in JILKINE_LENGTHS:
            end_time = max(2000.0, 500.0 * (length / 10.0) ** 2)
            for replicate in replicates:
                task = _base_task(
                    model,
                    "competing_fronts",
                    length,
                    _jilkine_cells(length),
                )
                task["protocol"] = {
                    "name": "localized_pulse",
                    "amplitude": 0.2,
                    "start": 0.0,
                    "stop": 40.0,
                    "taper": 20.0,
                    "width_fraction": 0.1,
                    "sides": ["left", "right"],
                    "asymmetry": 0.02,
                }
                task["simulation"].update(
                    {
                        "end_time": end_time,
                        "output_interval": max(2.0, end_time / 800.0),
                        "backend": "bdf",
                        "max_step": 10.0,
                    }
                )
                task["metrics"].update(
                    {"analysis_start": 40.0, "sustained_duration": 50.0}
                )
                task.update(
                    {
                        "replicate": replicate,
                        "seed": replicate + 1,
                        "save_trajectory": (
                            replicate == 0 and length in selected_trajectories
                        ),
                    }
                )
                if model != "legi":
                    task["initial_noise"] = {
                        "amplitude": 0.01,
                        "seed": replicate + 1,
                    }
                tasks.append(task)

            conditioning = max(200.0, 100.0 * (length / 10.0) ** 2)
            post_reversal = max(400.0, 200.0 * (length / 10.0) ** 2)
            task = _base_task(
                model,
                "gradient_reversal",
                length,
                _jilkine_cells(length),
            )
            task["protocol"] = {
                "name": "reversing_gradient",
                # Hold the peak activation rate fixed at 0.1/s across sizes.
                "slope": 0.1 / length,
                "reversal_time": conditioning,
                "transition": min(20.0, 0.1 * conditioning),
            }
            task["simulation"].update(
                {
                    "end_time": conditioning + post_reversal,
                    "output_interval": max(2.0, conditioning / 100.0),
                    "backend": "bdf",
                    "max_step": 10.0,
                }
            )
            task["metrics"].update(
                {
                    "analysis_start": 0.0,
                    "sustained_duration": max(10.0, conditioning / 50.0),
                    "reversal_time": conditioning,
                }
            )
            task["save_trajectory"] = length in selected_trajectories
            tasks.append(task)

    for model in JILKINE_MODELS:
        for length in TRACKING_LENGTHS:
            for period in TRACKING_PERIODS:
                task = _base_task(
                    model,
                    "dynamic_tracking",
                    length,
                    _jilkine_cells(length),
                )
                task["protocol"] = {
                    "name": "alternating_gradient",
                    # Hold peak input fixed rather than strengthening the cue in larger cells.
                    "slope": 0.1 / length,
                    "period": period,
                    "start": 0.0,
                }
                task["simulation"].update(
                    {
                        "end_time": 6.0 * period,
                        "output_interval": period / 20.0,
                        "backend": "bdf",
                        "max_step": min(5.0, period / 20.0),
                    }
                )
                task["metrics"]["sustained_duration"] = period / 4.0
                task["tracking"] = {
                    "period": period,
                    "start": period,
                    "minimum_amplitude": 0.05,
                    "maximum_lag_fraction": 0.25,
                }
                task["save_trajectory"] = length in {5.0, 39.0} and period in {
                    50.0,
                    400.0,
                }
                tasks.append(task)
    return tasks


def debelly_tasks() -> list[dict]:
    tasks = []
    selected_trajectories = {20.0, 2 ** (1 / 3) * 20.0, 64.0}
    for variant, parameters in DEBELLY_VARIANTS.items():
        for length in DEBELLY_LENGTHS:
            common = {
                "coupling_variant": variant,
                "model_parameters": parameters,
                "mechanochemical_metrics": {"stimulus_start": 50.0},
            }
            task = _base_task(
                "debelly",
                "single_site_propagation",
                length,
                _debelly_cells(length),
            )
            task.update(common)
            task["protocol"] = {
                "name": "gaussian_pulse",
                "amplitude": 1.0,
                "start": 50.0,
                "stop": 150.0,
                "width_fraction": 0.2,
                "centers": [0.0],
            }
            task["simulation"].update(
                {"end_time": 700.0, "output_interval": 2.0, "step": 0.5}
            )
            task["metrics"].update(
                {"analysis_start": 50.0, "sustained_duration": 50.0}
            )
            task["save_trajectory"] = (
                length in selected_trajectories and variant in {"intact", "both_off"}
            )
            tasks.append(task)

            task = _base_task(
                "debelly",
                "dual_site_challenge",
                length,
                _debelly_cells(length),
            )
            task.update(common)
            task["protocol"] = {
                "name": "gaussian_pulse",
                "amplitude": 1.0,
                "start": 50.0,
                "stop": 150.0,
                "width_fraction": 0.08,
                "centers": [0.0, 1.0],
                "asymmetry": 0.02,
            }
            task["simulation"].update(
                {"end_time": 700.0, "output_interval": 2.0, "step": 0.5}
            )
            task["metrics"].update(
                {"analysis_start": 50.0, "sustained_duration": 50.0}
            )
            task["save_trajectory"] = (
                length in selected_trajectories and variant in {"intact", "both_off"}
            )
            tasks.append(task)

            reversal_time = 250.0
            task = _base_task(
                "debelly",
                "gradient_reversal",
                length,
                _debelly_cells(length),
            )
            task.update(
                {
                    "coupling_variant": variant,
                    "model_parameters": parameters,
                    "mechanochemical_metrics": {"stimulus_start": 0.0},
                }
            )
            task["protocol"] = {
                "name": "reversing_gradient",
                "slope": 1.0 / length,
                "reversal_time": reversal_time,
                "transition": 25.0,
            }
            task["simulation"].update(
                {"end_time": 750.0, "output_interval": 2.0, "step": 0.5}
            )
            task["metrics"].update(
                {
                    "analysis_start": 0.0,
                    "sustained_duration": 25.0,
                    "reversal_time": reversal_time,
                }
            )
            task["save_trajectory"] = (
                length in selected_trajectories and variant in {"intact", "both_off"}
            )
            tasks.append(task)

    for variant, parameters in DEBELLY_VARIANTS.items():
        for length in DEBELLY_TRACKING_LENGTHS:
            for period in (100.0, 200.0, 400.0):
                task = _base_task(
                    "debelly",
                    "dynamic_tracking",
                    length,
                    _debelly_cells(length),
                )
                task.update(
                    {
                        "coupling_variant": variant,
                        "model_parameters": parameters,
                        "mechanochemical_metrics": {"stimulus_start": 0.0},
                    }
                )
                task["protocol"] = {
                    "name": "alternating_gradient",
                    "slope": 1.0 / length,
                    "period": period,
                    "start": 0.0,
                }
                task["simulation"].update(
                    {
                        "end_time": 6.0 * period,
                        "output_interval": period / 20.0,
                        "step": 0.5,
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
                    length in {20.0, 64.0}
                    and period in {100.0, 400.0}
                    and variant in {"intact", "both_off"}
                )
                tasks.append(task)
    return tasks


def refinement_tasks() -> list[dict]:
    tasks = []
    for model in JILKINE_MODELS:
        length = 39.0
        end_time = 500.0 * (length / 10.0) ** 2
        task = _base_task(
            model,
            "competing_fronts_refinement",
            length,
            _jilkine_cells(length, multiplier=2),
        )
        task["tier"] = "refinement"
        task["protocol"] = {
            "name": "localized_pulse",
            "amplitude": 0.2,
            "start": 0.0,
            "stop": 40.0,
            "taper": 20.0,
            "width_fraction": 0.1,
            "sides": ["left", "right"],
            "asymmetry": 0.02,
        }
        task["simulation"].update(
            {
                "end_time": end_time,
                "output_interval": end_time / 800.0,
                "backend": "bdf",
                "rtol": 2e-7,
                "atol": 1e-10,
                "max_step": 5.0,
            }
        )
        task["metrics"].update(
            {"analysis_start": 40.0, "sustained_duration": 50.0}
        )
        if model != "legi":
            task["initial_noise"] = {"amplitude": 0.01, "seed": 1}
        task.update({"replicate": 0, "seed": 1, "save_trajectory": True})
        tasks.append(task)
    for variant in ("intact", "both_off"):
        for length in (20.0, 64.0):
            task = _base_task(
                "debelly",
                "single_site_propagation_refinement",
                length,
                _debelly_cells(length, multiplier=2),
            )
            task.update(
                {
                    "tier": "refinement",
                    "coupling_variant": variant,
                    "model_parameters": DEBELLY_VARIANTS[variant],
                    "protocol": {
                        "name": "gaussian_pulse",
                        "amplitude": 1.0,
                        "start": 50.0,
                        "stop": 150.0,
                        "width_fraction": 0.2,
                        "centers": [0.0],
                    },
                    "mechanochemical_metrics": {"stimulus_start": 50.0},
                    "save_trajectory": True,
                }
            )
            task["simulation"].update(
                {"end_time": 700.0, "output_interval": 2.0, "step": 0.25}
            )
            task["metrics"].update(
                {"analysis_start": 50.0, "sustained_duration": 50.0}
            )
            tasks.append(task)
    return tasks


def build_manifest() -> dict:
    tasks = jilkine_tasks() + debelly_tasks() + refinement_tasks()
    return make_manifest(
        "phase1_production_v1",
        tasks,
        description=(
            "Common-phenotype size, reversal, tracking, and De Belly coupling campaign "
            "on the frozen Jilkine/accepted De Belly implementation baseline."
        ),
        metadata={
            "jilkine_reproduction": "frozen qualitative equation-level reconstruction",
            "debelly_crosscheck": "not required; optimized port accepted for WP1 screening",
            "symmetric_cue_policy": (
                "production competition uses declared 2% cue asymmetry and seeded noise; "
                "solver-seeded winners are not interpreted biologically"
            ),
            "gradient_policy": "peak scalar input is fixed across lengths",
            "mesh_policy": "approximately fixed physical dx, with selected refinements",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=ROOT / "campaigns" / "phase1_production_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
