#!/usr/bin/env python3
"""Create a focused trial of stronger finite pulses and long autonomous windows."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_long_window_reversal_trial_v1"
LENGTHS = (5.0, 15.0, 25.0, 39.0)
INDUCTION_RATIOS = (4.0, 8.0)
ALPHA0 = BETA0 = 1.8
ALPHA = BETA = 0.54
INDUCTION_STOP = 640.0
OPPOSITE_START = 10_000.0
OPPOSITE_STOP = 10_640.0
RESPONSE_DEADLINE = 20_640.0
END_TIME = 21_640.0
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


def _pulse(side: str, start: float, stop: float, ratio: float) -> dict:
    return {
        "amplitude": ratio * ALPHA0,
        "start": start,
        "stop": stop,
        "taper": (stop - start) / 2.0,
        "width_fraction": 0.4,
        "sides": [side],
    }


def campaign_tasks() -> list[dict]:
    tasks = []
    for ratio in INDUCTION_RATIOS:
        for variant, intervention in VARIANTS.items():
            for length in LENGTHS:
                tasks.append({
                    "tier": "debelly_2026_long_window_reversal_trial",
                    "implementation_revision": "reflection_symmetric_responsive_regime_v2",
                    "model": "debelly_2026",
                    "model_label": "Reflection-symmetric extension of the 2026 de Belly mechanochemical model",
                    "assay": "long_window_finite_reversal_trial",
                    "coupling_variant": variant,
                    "induction_ratio": ratio,
                    "input_peak": ratio,
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
                            _pulse("left", 0.0, INDUCTION_STOP, ratio),
                            _pulse("right", OPPOSITE_START, OPPOSITE_STOP, ratio),
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
                        "analysis_start": INDUCTION_STOP,
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
                    "save_trajectory": length == 15.0,
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
            "Focused sensitivity trial: stronger symmetric finite pulses, 10 ks for "
            "initial autonomous state formation and 10 ks for autonomous reversal."
        ),
        metadata={
            "induction_ratios": list(INDUCTION_RATIOS),
            "lengths_um": list(LENGTHS),
            "initial_pulse": "left outer 40%, 0-640 s, cosine taper over final 320 s",
            "initial_qualification_time_s": OPPOSITE_START,
            "opposite_pulse": "right outer 40%, 10000-10640 s, cosine taper over final 320 s",
            "reversal_onset_window_s": RESPONSE_DEADLINE - OPPOSITE_STOP,
            "minimum_persistence_s": 1_000.0,
        },
    )
    destination = write_manifest(args.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
