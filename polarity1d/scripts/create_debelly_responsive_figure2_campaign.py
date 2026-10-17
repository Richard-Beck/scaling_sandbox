#!/usr/bin/env python3
"""Create the approved responsive-regime production campaign for Figure 2."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_responsive_figure2_v2"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
PERIODS = (100.0, 200.0, 400.0, 800.0, 1_600.0)
ALPHA0 = 1.8
BETA0 = 1.8
ALPHA = 0.3 * ALPHA0
BETA = 0.3 * BETA0
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


def _cells(length: float) -> int:
    return max(5, round(length / 0.2))


def _pulse(side: str, start: float, duration: float, ratio: float, width: float, *, taper=True):
    return {
        "amplitude": ratio * ALPHA0,
        "start": start,
        "stop": start + duration,
        "taper": duration / 2.0 if taper else 0.0,
        "width_fraction": width,
        "sides": [side],
    }


def _base(assay: str, length: float, variant: str) -> dict:
    return {
        "tier": "debelly_2026_responsive_figure2",
        "implementation_revision": "reflection_symmetric_responsive_regime_v2",
        "model": "debelly_2026",
        "model_label": (
            "Reflection-symmetric extension of the 2026 de Belly mechanochemical model"
        ),
        "assay": assay,
        "coupling_variant": variant,
        "model_parameters": {
            "initialization_mode": "algebraic_fixed_point",
            "assay_mode": "reflection_symmetric_extension",
            "input_sampling_mode": "cell_centred",
            "rac_basal_activation": ALPHA0,
            "rho_basal_activation": BETA0,
            "rac_mechanical_feedback": ALPHA,
            "rho_mechanical_feedback": BETA,
            **VARIANTS[variant],
        },
        "protocol": {"name": "none"},
        "simulation": {
            "length": length,
            "cells": _cells(length),
            "end_time": 1_000.0,
            "output_interval": 5.0,
            "step": 0.5,
        },
        "metrics": {
            "analysis_start": 0.0,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 50.0,
        },
        "field_audit": True,
        "replicate": 0,
        "save_trajectory": False,
    }


def campaign_tasks() -> list[dict]:
    tasks = []
    for variant in VARIANTS:
        for length in LENGTHS:
            competition = _base("responsive_competing_fronts", length, variant)
            competition.update(
                input_peak=2.0,
                protocol={
                    "name": "localized_pulse",
                    "amplitude": 2.0 * ALPHA0,
                    "start": 0.0,
                    "stop": 80.0,
                    "taper": 40.0,
                    "width_fraction": 0.2,
                    "sides": ["left", "right"],
                    "asymmetry": 0.15,
                },
                strong_competition={
                    "cue_stop": 80.0,
                    "two_pole_deadline": 230.0,
                    "two_pole_strength_threshold": 0.05,
                    "two_pole_sustain_s": 20.0,
                    "contrast_threshold": 0.3,
                    "orientation_threshold": 0.25,
                    "resolution_sustain_s": 100.0,
                    "endpoint_sustain_s": 500.0,
                },
            )
            competition["simulation"].update({"end_time": 20_000.0, "output_interval": 5.0})
            competition["metrics"]["analysis_start"] = 80.0
            competition["save_trajectory"] = variant == "full" and length == 15.0
            tasks.append(competition)

            reversal = _base("responsive_finite_reversal", length, variant)
            reversal.update(
                input_peak=4.0,
                protocol={
                    "name": "finite_pulse_sequence",
                    "pulses": [
                        _pulse("left", 0.0, 160.0, 2.0, 0.2),
                        _pulse("right", 600.0, 320.0, 4.0, 0.4),
                    ],
                },
                strong_reversal={
                    "opposite_start": 600.0,
                    "opposite_stop": 920.0,
                    "response_deadline": 1_420.0,
                    "contrast_threshold": 0.3,
                    "orientation_threshold": 0.25,
                    "sustain_s": 100.0,
                    "minimum_persistence_s": 1_000.0,
                },
            )
            reversal["simulation"].update({"end_time": 2_500.0, "output_interval": 2.0})
            reversal["save_trajectory"] = variant == "full" and length == 15.0
            tasks.append(reversal)

            for period in PERIODS:
                half_period = period / 2.0
                train_start = 600.0
                half_cycles = 10
                pulses = [_pulse("left", 0.0, 160.0, 2.0, 0.2)]
                for event in range(half_cycles):
                    pulses.append(_pulse(
                        "right" if event % 2 == 0 else "left",
                        train_start + event * half_period,
                        half_period,
                        4.0,
                        0.4,
                        taper=False,
                    ))
                tracking = _base("responsive_pulse_train", length, variant)
                tracking.update(
                    input_peak=4.0,
                    tracking_period_s=period,
                    protocol={"name": "finite_pulse_sequence", "pulses": pulses},
                    pulse_train_frequency={
                        "start": train_start,
                        "period": period,
                        "half_cycles": half_cycles,
                        "contrast_threshold": 0.3,
                        "orientation_threshold": 0.25,
                    },
                )
                tracking["simulation"].update({
                    "end_time": train_start + half_cycles * half_period,
                    "output_interval": min(5.0, period / 50.0),
                })
                tracking["save_trajectory"] = (
                    variant == "full" and length == 15.0 and period == 400.0
                )
                tasks.append(tracking)
    return tasks


def build_manifest() -> dict:
    return make_manifest(
        CAMPAIGN_NAME,
        campaign_tasks(),
        description=(
            "Approved responsive-regime Figure 2 campaign: strict two-pole resolution, "
            "bounded finite-pulse reversal, and deterministic pulse-train frequency response."
        ),
        metadata={
            "model_provenance": (
                "2026 de Belly mechanochemical Rac-Rho-MCA model with the documented "
                "reflection-symmetric boundary extension"
            ),
            "lengths_um": list(LENGTHS),
            "mesh_policy": "cells=max(5, round(L/0.2)); approximately dx=0.2 um",
            "time_step_s": 0.5,
            "parameter_selection": {
                "alpha0": ALPHA0,
                "beta0": BETA0,
                "alpha": ALPHA,
                "beta": BETA,
                "selected_at_length_um": 15.0,
                "selected_before_production_length_sweep": True,
                "basis": (
                    "reviewed L=15 kymographs showing strong two-pole establishment and "
                    "resolution, autonomous finite-pulse reversal, and a frequency cutoff"
                ),
            },
            "competition_protocol": (
                "two outer-20% cues at 2*alpha0, 80 s with 40-s cosine taper, "
                "15% left/right asymmetry; 20,000-s common censoring horizon"
            ),
            "reversal_protocol": (
                "left induction 2*alpha0 outer 20% for 160 s; opposite right cue "
                "4*alpha0 outer 40% from 600-920 s; response must begin by 1420 s "
                "and persist at least 1000 s"
            ),
            "frequency_protocol": (
                "after the same induction, ten alternating square half-cycles at "
                "4*alpha0 over the outer 40%; full-cycle periods 100-1600 s"
            ),
            "variants": {
                "full": "intact model",
                "beta_zero": "beta=0 only",
                "alpha_zero": "alpha=0 only",
                "local_inhibition": "v_CTE=alpha=beta=sigma_a0=0",
            },
        },
    )


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
