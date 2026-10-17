#!/usr/bin/env python3
"""Create the de Belly single-cue width/amplitude/duration trial campaign."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_single_cue_kinetics_trial_v1"
LENGTH = 15.0
CELLS = 75
WIDTHS = (0.025, 0.05, 0.1, 0.2, 0.3)
AMPLITUDES = (0.1, 0.25, 0.5, 1.0, 2.0)
DURATIONS = (10.0, 20.0, 40.0, 80.0, 160.0)
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
    for variant, overrides in VARIANTS.items():
        for width in WIDTHS:
            for amplitude in AMPLITUDES:
                for duration in DURATIONS:
                    tasks.append(
                        {
                            "tier": "tmp_trial_single_cue_kinetics",
                            "implementation_revision": "single_cue_kinetics_v1",
                            "model": "debelly_2026",
                            "model_label": (
                                "Reflection-symmetric extension of the 2026 de Belly "
                                "mechanochemical model"
                            ),
                            "assay": "single_rac_cue_kinetics",
                            "coupling_variant": variant,
                            "input_peak": amplitude,
                            "model_parameters": {
                                "assay_mode": "reflection_symmetric_extension",
                                "initialization_mode": "algebraic_fixed_point",
                                "input_sampling_mode": "cell_centred",
                                **overrides,
                            },
                            "protocol": {
                                "name": "localized_pulse",
                                "amplitude": amplitude,
                                "start": 0.0,
                                "stop": duration,
                                "taper": duration / 2.0,
                                "width_fraction": width,
                                "sides": ["left"],
                                "asymmetry": 0.0,
                            },
                            "simulation": {
                                "length": LENGTH,
                                "cells": CELLS,
                                "end_time": duration + 3_000.0,
                                "output_interval": 1.0,
                                "step": 0.5,
                            },
                            "metrics": {
                                "analysis_start": 0.0,
                                "polarization_threshold": 0.1,
                                "orientation_threshold": 0.05,
                                "sustained_duration": 10.0,
                            },
                            "single_cue_kinetics": {
                                "cue_start": 0.0,
                                "cue_stop": duration,
                                "cue_width_fraction": width,
                                "cue_amplitude": amplitude,
                                "orientation_threshold": 0.05,
                                "polarization_sustain_s": 10.0,
                                "depolarization_sustain_s": 50.0,
                                "cue_side": "left",
                            },
                            "field_audit": True,
                            "save_trajectory": False,
                            "replicate": 0,
                        }
                    )
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
            "Temporary L=15 single-left-Rac-cue sweep over width, amplitude and "
            "duration for four independent mechanochemical arm variants."
        ),
        metadata={
            "status": "temporary trial; does not feed manuscript Figure 2",
            "length_um": LENGTH,
            "dx_um": LENGTH / CELLS,
            "dt_s": 0.5,
            "width_fractions": list(WIDTHS),
            "amplitudes_alpha_cue_over_alpha0": list(AMPLITUDES),
            "total_cue_durations_s": list(DURATIONS),
            "cue_shape": "full for first half, cosine taper during second half",
            "post_withdrawal_observation_s": 3_000.0,
            "orientation_threshold": 0.05,
            "polarization_sustain_s": 10.0,
            "depolarization_sustain_s": 50.0,
            "variants": list(VARIANTS),
        },
    )
    write_manifest(args.output, manifest)
    print(f"{args.output}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
