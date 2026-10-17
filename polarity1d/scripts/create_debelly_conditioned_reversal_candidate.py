#!/usr/bin/env python3
"""Create the selected long-conditioning finite-reversal candidate sweep."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_conditioned_reversal_candidate_v1"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
ALPHA0 = BETA0 = 1.8
ALPHA = BETA = 0.54
CONDITIONING_RATIO = 2.0
OPPOSITE_RATIO = 8.0
OPPOSITE_START = 10_000.0
OPPOSITE_STOP = 10_320.0
RESPONSE_DEADLINE = 20_320.0
END_TIME = 21_320.0
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


def campaign_tasks() -> list[dict]:
    tasks = []
    for variant, intervention in VARIANTS.items():
        for length in LENGTHS:
            tasks.append({
                "tier": "debelly_2026_conditioned_reversal_candidate",
                "implementation_revision": "reflection_symmetric_responsive_regime_v2",
                "model": "debelly_2026",
                "model_label": "Reflection-symmetric extension of the 2026 de Belly mechanochemical model",
                "assay": "conditioned_finite_reversal_candidate",
                "coupling_variant": variant,
                "input_peak": OPPOSITE_RATIO,
                "model_parameters": {
                    "initialization_mode": "algebraic_fixed_point",
                    "assay_mode": "reflection_symmetric_extension",
                    "input_sampling_mode": "cell_centred",
                    "rac_basal_activation": ALPHA0,
                    "rho_basal_activation": BETA0,
                    "rac_mechanical_feedback": ALPHA,
                    "rho_mechanical_feedback": BETA,
                    **intervention,
                },
                "protocol": {
                    "name": "finite_pulse_sequence",
                    "pulses": [
                        {
                            "amplitude": CONDITIONING_RATIO * ALPHA0,
                            "start": 0.0,
                            "stop": OPPOSITE_START,
                            "taper": 0.0,
                            "width_fraction": 0.4,
                            "sides": ["left"],
                        },
                        {
                            "amplitude": OPPOSITE_RATIO * ALPHA0,
                            "start": OPPOSITE_START,
                            "stop": OPPOSITE_STOP,
                            "taper": 160.0,
                            "width_fraction": 0.4,
                            "sides": ["right"],
                        },
                    ],
                },
                "simulation": {
                    "length": length,
                    "cells": max(5, round(length / 0.2)),
                    "end_time": END_TIME,
                    "output_interval": 5.0,
                    "step": 0.5,
                },
                "metrics": {
                    "analysis_start": OPPOSITE_START,
                    "polarization_threshold": 0.1,
                    "orientation_threshold": 0.05,
                    "sustained_duration": 50.0,
                },
                "strong_reversal": {
                    "opposite_start": OPPOSITE_START,
                    "opposite_stop": OPPOSITE_STOP,
                    "response_deadline": RESPONSE_DEADLINE,
                    "contrast_threshold": 0.3,
                    "orientation_threshold": 0.25,
                    "sustain_s": 100.0,
                    "minimum_persistence_s": 1_000.0,
                },
                "field_audit": True,
                "replicate": 0,
                "save_trajectory": variant == "full" and length == 15.0,
            })
    return tasks


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json",
    )
    args = parser.parse_args(argv)
    manifest = make_manifest(
        CAMPAIGN_NAME,
        campaign_tasks(),
        description=(
            "Preselected candidate: 2*alpha0 left conditioning for 10 ks, finite "
            "8*alpha0 right pulse for 320 s, and 10 ks autonomous follow-up."
        ),
        metadata={
            "conditioning": "left, 2*alpha0, outer 40%, sustained from 0-10000 s",
            "opposite_pulse": "right, 8*alpha0, outer 40%, 10000-10320 s, 160-s taper",
            "autonomous_reversal_window_s": RESPONSE_DEADLINE - OPPOSITE_STOP,
            "minimum_persistence_s": 1_000.0,
            "selection": "selected at L=15 before the length sweep",
        },
    )
    destination = write_manifest(args.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
