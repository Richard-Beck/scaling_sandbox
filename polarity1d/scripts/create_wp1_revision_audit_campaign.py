#!/usr/bin/env python3
"""Create the focused campaign needed to revise WP1 Figures 1 and 2."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


JILKINE_MODELS = ("wave_pinning", "otsuji", "goryachev", "legi")
JILKINE_LENGTHS = tuple(float(value) for value in range(1, 40, 2))
DEBELLY_LENGTHS = (
    12.0, 16.0, 20.0, 2 ** (1 / 3) * 20.0, 32.0, 36.0, 40.0,
    44.0, 48.0, 50.0, 52.0, 56.0, 60.0, 64.0,
)


def _jilkine_cells(length: float) -> int:
    return max(101, round(length / 0.1))


def _debelly_cells(length: float) -> int:
    return max(60, round(length / 0.4))


def build_manifest() -> dict:
    tasks = []

    # Four times the previous physical and L^2-scaled observation policy.
    for model in JILKINE_MODELS:
        replicates = range(5) if model != "legi" else range(1)
        for length in JILKINE_LENGTHS:
            end_time = max(8000.0, 2000.0 * (length / 10.0) ** 2)
            for replicate in replicates:
                task = {
                    "tier": "revision_audit",
                    "model": model,
                    "assay": "competing_fronts_extended_4x",
                    "protocol": {
                        "name": "localized_pulse", "amplitude": 0.2,
                        "start": 0.0, "stop": 40.0, "taper": 20.0,
                        "width_fraction": 0.1, "sides": ["left", "right"],
                        "asymmetry": 0.02,
                    },
                    "simulation": {
                        "length": length, "cells": _jilkine_cells(length),
                        "end_time": end_time,
                        "output_interval": max(2.0, end_time / 1600.0),
                        "backend": "bdf", "max_step": 10.0,
                    },
                    "metrics": {
                        "analysis_start": 40.0, "polarization_threshold": 0.1,
                        "orientation_threshold": 0.05, "sustained_duration": 50.0,
                    },
                    "replicate": replicate, "seed": replicate + 1,
                    "save_trajectory": replicate == 0 and length in {5.0, 21.0, 39.0},
                }
                if model != "legi":
                    task["initial_noise"] = {"amplitude": 0.01, "seed": replicate + 1}
                tasks.append(task)

    # Re-run Otsuji reversal with the current pre-reversal diagnostics and all trajectories.
    for length in JILKINE_LENGTHS:
        conditioning = max(200.0, 100.0 * (length / 10.0) ** 2)
        post_reversal = max(400.0, 200.0 * (length / 10.0) ** 2)
        tasks.append(
            {
                "tier": "revision_audit", "model": "otsuji",
                "assay": "otsuji_gradient_reversal_diagnostic",
                "protocol": {
                    "name": "reversing_gradient", "slope": 0.1 / length,
                    "reversal_time": conditioning,
                    "transition": min(20.0, 0.1 * conditioning),
                },
                "simulation": {
                    "length": length, "cells": _jilkine_cells(length),
                    "end_time": conditioning + post_reversal,
                    "output_interval": max(1.0, conditioning / 200.0),
                    "backend": "bdf", "max_step": 5.0,
                },
                "metrics": {
                    "analysis_start": 0.0, "polarization_threshold": 0.1,
                    "orientation_threshold": 0.05,
                    "sustained_duration": max(10.0, conditioning / 50.0),
                    "reversal_time": conditioning,
                },
                "replicate": 0, "save_trajectory": True,
            }
        )

    # Intact De Belly fields across both cue strengths, including the 40–50 µm transition.
    for input_peak in (0.2, 0.5):
        for length in DEBELLY_LENGTHS:
            conditioning = max(250.0, 250.0 * (length / 20.0) ** 2)
            post_reversal = max(750.0, 750.0 * (length / 20.0) ** 2)
            tasks.append(
                {
                    "tier": "revision_audit", "model": "debelly",
                    "assay": "debelly_intact_reversal_field_audit",
                    "coupling_variant": "intact", "model_parameters": {},
                    "input_peak": input_peak,
                    "protocol": {
                        "name": "reversing_gradient", "slope": input_peak / length,
                        "reversal_time": conditioning, "transition": 25.0,
                    },
                    "simulation": {
                        "length": length, "cells": _debelly_cells(length),
                        "end_time": conditioning + post_reversal,
                        "output_interval": 2.0, "step": 0.5,
                    },
                    "metrics": {
                        "analysis_start": 0.0, "polarization_threshold": 0.1,
                        "orientation_threshold": 0.05, "sustained_duration": 25.0,
                        "reversal_time": conditioning,
                    },
                    "mechanochemical_metrics": {"stimulus_start": 0.0},
                    "replicate": 0, "save_trajectory": True,
                }
            )

    return make_manifest(
        "wp1_revision_audit_v1",
        tasks,
        description=(
            "Fourfold-longer Jilkine competition, explicit Otsuji conditioning "
            "diagnostics, and intact De Belly multi-field reversal trajectories."
        ),
        metadata={
            "competition_policy": "max(8000, 2000*(L/10)^2) s; 4x phase1_production_v1",
            "competition_censoring": "unresolved sustained multifront states are right-censored",
            "otsuji_purpose": "separate unpolarized, multifront, weakly oriented, and conditioned prestates",
            "debelly_purpose": "audit Rac-only reversal calls against Rho, MCA, and mechanics",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "wp1_revision_audit_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
