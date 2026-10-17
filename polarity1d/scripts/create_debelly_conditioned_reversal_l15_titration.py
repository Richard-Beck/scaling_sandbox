#!/usr/bin/env python3
"""Titrate conditioning and finite opposite pulses at L=15 in the full model."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_conditioned_reversal_l15_titration_v1"
ALPHA0 = BETA0 = 1.8
ALPHA = BETA = 0.54
OPPOSITE_START = 10_000.0


def campaign_tasks() -> list[dict]:
    tasks = []
    for conditioning_ratio in (2.0, 4.0):
        for opposite_ratio in (4.0, 8.0):
            for opposite_duration in (320.0, 640.0):
                opposite_stop = OPPOSITE_START + opposite_duration
                deadline = opposite_stop + 10_000.0
                assay = (
                    f"conditioned_reversal_c{conditioning_ratio:g}_"
                    f"r{opposite_ratio:g}_d{opposite_duration:g}"
                )
                tasks.append({
                    "tier": "debelly_2026_conditioned_reversal_l15_titration",
                    "implementation_revision": "reflection_symmetric_responsive_regime_v2",
                    "model": "debelly_2026",
                    "model_label": "Reflection-symmetric extension of the 2026 de Belly mechanochemical model",
                    "assay": assay,
                    "coupling_variant": "full",
                    "input_peak": opposite_ratio,
                    "model_parameters": {
                        "initialization_mode": "algebraic_fixed_point",
                        "assay_mode": "reflection_symmetric_extension",
                        "input_sampling_mode": "cell_centred",
                        "rac_basal_activation": ALPHA0,
                        "rho_basal_activation": BETA0,
                        "rac_mechanical_feedback": ALPHA,
                        "rho_mechanical_feedback": BETA,
                    },
                    "protocol": {
                        "name": "finite_pulse_sequence",
                        "pulses": [
                            {
                                "amplitude": conditioning_ratio * ALPHA0,
                                "start": 0.0,
                                "stop": OPPOSITE_START,
                                "taper": 0.0,
                                "width_fraction": 0.4,
                                "sides": ["left"],
                            },
                            {
                                "amplitude": opposite_ratio * ALPHA0,
                                "start": OPPOSITE_START,
                                "stop": opposite_stop,
                                "taper": opposite_duration / 2.0,
                                "width_fraction": 0.4,
                                "sides": ["right"],
                            },
                        ],
                    },
                    "simulation": {
                        "length": 15.0,
                        "cells": 75,
                        "end_time": deadline + 1_000.0,
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
                        "opposite_stop": opposite_stop,
                        "response_deadline": deadline,
                        "contrast_threshold": 0.3,
                        "orientation_threshold": 0.25,
                        "sustain_s": 100.0,
                        "minimum_persistence_s": 1_000.0,
                    },
                    "field_audit": True,
                    "replicate": 0,
                    "save_trajectory": True,
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
        description="L=15 full-model titration for the long-conditioned finite reversal assay.",
    )
    destination = write_manifest(args.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
