#!/usr/bin/env python3
"""Create a strong 10-ks conditioning / 10-ks reversal sensitivity campaign."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_long_conditioning_reversal_trial_v1"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
ALPHA0 = BETA0 = 1.8
ALPHA = BETA = 0.54
CUE_RATIO = 4.0
SWITCH_TIME = 10_000.0
REVERSE_CUE_STOP = 20_000.0
END_TIME = 21_000.0
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


def _pulse(side: str, start: float, stop: float) -> dict:
    return {
        "amplitude": CUE_RATIO * ALPHA0,
        "start": start,
        "stop": stop,
        "taper": 0.0,
        "width_fraction": 0.4,
        "sides": [side],
    }


def campaign_tasks() -> list[dict]:
    tasks = []
    for variant, intervention in VARIANTS.items():
        for length in LENGTHS:
            tasks.append({
                "tier": "debelly_2026_long_conditioning_reversal_trial",
                "implementation_revision": "reflection_symmetric_responsive_regime_v2",
                "model": "debelly_2026",
                "model_label": "Reflection-symmetric extension of the 2026 de Belly mechanochemical model",
                "assay": "long_conditioning_reversal_trial",
                "coupling_variant": variant,
                "input_peak": CUE_RATIO,
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
                        _pulse("left", 0.0, SWITCH_TIME),
                        _pulse("right", SWITCH_TIME, REVERSE_CUE_STOP),
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
                    "analysis_start": SWITCH_TIME,
                    "polarization_threshold": 0.1,
                    "orientation_threshold": 0.05,
                    "sustained_duration": 50.0,
                },
                # Setting opposite_stop to the cue-switch time makes the strict
                # search begin at reversal onset. The reverse cue itself lasts
                # until 20 ks; the final 1 ks tests persistence after withdrawal.
                "strong_reversal": {
                    "opposite_start": SWITCH_TIME,
                    "opposite_stop": SWITCH_TIME,
                    "response_deadline": REVERSE_CUE_STOP,
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
            "Sensitivity trial with 10 ks of strong left conditioning, 10 ks of "
            "strong right reversal cue and a 1-ks cue-free persistence check."
        ),
        metadata={
            "cue_ratio_alpha0": CUE_RATIO,
            "cue_width_fraction": 0.4,
            "conditioning_interval_s": SWITCH_TIME,
            "reversal_interval_s": REVERSE_CUE_STOP - SWITCH_TIME,
            "post_withdrawal_persistence_check_s": END_TIME - REVERSE_CUE_STOP,
            "scoring_note": (
                "The strong-reversal latency field is latency after cue switch; "
                "the reverse cue remains present until 20000 s."
            ),
        },
    )
    destination = write_manifest(args.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
