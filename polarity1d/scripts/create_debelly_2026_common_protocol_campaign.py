#!/usr/bin/env python3
"""Create the de Belly 2026 campaign under the exact Jilkine comparison protocol."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_jilkine_common_v1"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
TRACKING_PERIODS = (50.0, 100.0, 200.0, 400.0)
VARIANTS = {
    "full": {},
    "beta_zero": {"rho_mechanical_feedback": 0.0},
    "alpha_zero": {"rac_mechanical_feedback": 0.0},
    "local_inhibition": {
        "protrusion_speed": 0.0,
        "rac_mechanical_feedback": 0.0,
        "rho_mechanical_feedback": 0.0,
        "active_stress": 0.0,
    },
}
VALIDATION_CASES = {
    ("competing_fronts_common", 5.0, None),
    ("gradient_reversal_common", 15.0, None),
    ("dynamic_tracking_common", 39.0, 200.0),
}
PRODUCTION_INITIALIZATION = "algebraic_fixed_point"
IMPLEMENTATION_REVISION = "reflection_symmetric_jilkine_extension_v10_selected_regime"
SELECTED_ALPHA0 = 3.5
SELECTED_BETA0 = 4.0
SELECTED_ALPHA = 0.3 * SELECTED_ALPHA0
SELECTED_BETA = 0.3 * SELECTED_BETA0


def _cells(length: float) -> int:
    """Keep the finite-volume cell width at approximately 0.2 micrometres."""
    return max(5, round(length / 0.2))


def _base(assay: str, length: float, variant: str) -> dict:
    return {
        "tier": "debelly_2026_common_protocol",
        "implementation_revision": IMPLEMENTATION_REVISION,
        "model": "debelly_2026",
        "model_label": (
            "Reflection-symmetric extension of the 2026 de Belly "
            "mechanochemical model"
        ),
        "assay": assay,
        "coupling_variant": variant,
        "model_parameters": {
            "initialization_mode": PRODUCTION_INITIALIZATION,
            "assay_mode": "reflection_symmetric_extension",
            "input_sampling_mode": "cell_centred",
            "rac_basal_activation": SELECTED_ALPHA0,
            "rho_basal_activation": SELECTED_BETA0,
            "rac_mechanical_feedback": SELECTED_ALPHA,
            "rho_mechanical_feedback": SELECTED_BETA,
            **VARIANTS[variant],
        },
        "protocol": {"name": "none"},
        "simulation": {
            "length": length,
            "cells": _cells(length),
            "end_time": 20_000.0,
            "output_interval": 5.0,
            "step": 0.5,
        },
        "metrics": {
            "analysis_start": 0.0,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 50.0,
        },
        "mechanochemical_metrics": {"stimulus_start": 0.0},
        "field_audit": True,
        "replicate": 0,
        "save_trajectory": False,
    }


def campaign_tasks() -> list[dict]:
    tasks: list[dict] = []

    # Explicit L=15 one-sided response gate before freezing the amplitude at 1.
    calibration = _base("cue_response_calibration", 15.0, "full")
    calibration["protocol"] = {
        "name": "localized_pulse",
        "amplitude": SELECTED_ALPHA0,
        "start": 0.0,
        "stop": 40.0,
        "taper": 20.0,
        "width_fraction": 0.1,
        "sides": ["left"],
    }
    calibration["simulation"].update({"end_time": 200.0, "output_interval": 0.5})
    calibration["metrics"].update({"analysis_start": 0.0, "sustained_duration": 10.0})
    calibration["save_trajectory"] = True
    tasks.append(calibration)

    for variant in VARIANTS:
        for length in LENGTHS:
            competition = _base("competing_fronts_common", length, variant)
            competition["protocol"] = {
                "name": "localized_pulse",
                "amplitude": SELECTED_ALPHA0,
                "start": 0.0,
                "stop": 40.0,
                "taper": 20.0,
                "width_fraction": 0.1,
                "sides": ["left", "right"],
                "asymmetry": 0.02,
            }
            competition["metrics"]["analysis_start"] = 40.0
            competition["save_trajectory"] = length == 15.0 or (
                variant == "full"
                and (competition["assay"], length, None) in VALIDATION_CASES
            )
            tasks.append(competition)

            conditioning = max(2000.0, 400.0 * (length / 10.0) ** 2)
            reversal = _base("gradient_reversal_common", length, variant)
            reversal["protocol"] = {
                "name": "reversing_gradient",
                "slope": SELECTED_ALPHA0 / length,
                "reversal_time": conditioning,
                "transition": 20.0,
            }
            reversal["simulation"].update(
                {"end_time": conditioning + 20_000.0, "output_interval": 5.0}
            )
            reversal["metrics"].update(
                {
                    "reversal_time": conditioning,
                    "sustained_duration": max(20.0, conditioning / 100.0),
                }
            )
            reversal["save_trajectory"] = length == 15.0 or (
                variant == "full"
                and (reversal["assay"], length, None) in VALIDATION_CASES
            )
            tasks.append(reversal)

    for variant in VARIANTS:
        for length in LENGTHS:
            for period in TRACKING_PERIODS:
                tracking = _base("dynamic_tracking_common", length, variant)
                tracking["protocol"] = {
                    "name": "alternating_gradient",
                    "slope": SELECTED_ALPHA0 / length,
                    "period": period,
                    "start": 0.0,
                }
                tracking["simulation"].update(
                    {
                        "end_time": 6.0 * period,
                        "output_interval": period / (80.0 if length == 15.0 else 20.0),
                    }
                )
                tracking["metrics"]["sustained_duration"] = period / 4.0
                tracking["tracking"] = {
                    "period": period,
                    "start": period,
                    "minimum_amplitude": 0.05,
                    "maximum_lag_fraction": 0.25,
                }
                tracking["save_trajectory"] = length == 15.0 or (
                    variant == "full"
                    and (tracking["assay"], length, period) in VALIDATION_CASES
                )
                tasks.append(tracking)
    return tasks


def build_manifest() -> dict:
    return make_manifest(
        CAMPAIGN_NAME,
        campaign_tasks(),
        description=(
            "2026 de Belly mechanochemical model under the Jilkine comparison "
            "lengths, cue geometry/timing, horizons, and phenotype definitions."
        ),
        metadata={
            "model_provenance": (
                "VirtualEmbryo/mechanochemical_polarization, 2026 de Belly "
                "mechanochemical Rac-Rho-MCA model; reflection-symmetric "
                "boundary extension for Jilkine assays"
            ),
            "assay_mode": (
                "reflection-symmetric extension with matched Rac-dependent "
                "protrusion at both ends and removal of rigid translation; "
                "not the literally released one-sided assay"
            ),
            "lengths_um": list(LENGTHS),
            "mesh_policy": "cells=max(5, round(L/0.2)); approximately dx=0.2 um",
            "time_step_s": 0.5,
            "implementation_revision": IMPLEMENTATION_REVISION,
            "parameter_selection": {
                "length_um": 15.0,
                "rac_basal_activation_alpha0": SELECTED_ALPHA0,
                "rho_basal_activation_beta0": SELECTED_BETA0,
                "rac_mechanical_feedback_alpha": SELECTED_ALPHA,
                "rho_mechanical_feedback_beta": SELECTED_BETA,
                "selection_basis": (
                    "interior of the persistent spontaneous-polarization regime; "
                    "40-s asymmetric cue creates long-lived Rac-front/Rho-back "
                    "polarity; mechanics removal abolishes persistence; sustained "
                    "opposing gradient reverses in 75 s at L=15"
                ),
                "selected_before_length_sweep": True,
            },
            "initialization_policy": (
                "production uses the timestep-independent homogeneous algebraic "
                "Rac-Rho fixed point; Julia-matched finite warm-up is validation-only"
            ),
            "mca_policy": (
                "resting density fixed at 0.1; implied total abundance proportional to L"
            ),
            "cue_policy": (
                "outer 10%; amplitude alpha_cue/alpha0=1; full to 20 s; "
                "cosine taper to zero at 40 s"
            ),
            "competition_horizon_s": 20_000.0,
            "reversal_follow_up_s": 20_000.0,
            "tracking_periods_s": list(TRACKING_PERIODS),
            "variants": {
                "full": "published parameters",
                "beta_zero": "beta=0 only",
                "alpha_zero": "alpha=0 only",
                "local_inhibition": "v_CTE=alpha=beta=sigma_a0=0",
            },
            "cross_validation_cases": [
                "competition L=5", "reversal L=15", "tracking T=200 L=39"
            ],
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
