#!/usr/bin/env python3
"""Validate reflection symmetry and both mechanochemical paths of the extension."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly,
    DeBellyParameters,
    LocalizedPulse,
    simulate,
)


def _run(side: str, **overrides):
    return simulate(
        DeBelly(
            DeBellyParameters(
                assay_mode="reflection_symmetric_extension",
                **overrides,
            )
        ),
        protocol=LocalizedPulse(
            amplitude=1.0,
            start=0.0,
            stop=40.0,
            taper=20.0,
            width_fraction=0.1,
            sides=(side,),
        ),
        length=15.0,
        cells=75,
        end_time=500.0,
        output_interval=1.0,
        step=0.5,
    )


def build_report() -> dict:
    left = _run("left")
    right = _run("right")
    reflection = {}
    for field in ("rac", "rho", "tension", "mca_bound", "mca_unbound"):
        reflection[field] = float(
            np.max(np.abs(left.fields[field] - right.fields[field][:, ::-1]))
        )
    for field in ("velocity", "displacement"):
        reflection[field] = float(
            np.max(np.abs(left.fields[field] + right.fields[field][:, ::-1]))
        )

    pathway = {}
    length = 15.0
    for side, full in (("left", left), ("right", right)):
        beta0 = _run(side, rho_mechanical_feedback=0.0)
        alpha0 = _run(side, rac_mechanical_feedback=0.0)
        proximal = full.x <= 0.2 * length if side == "left" else full.x >= 0.8 * length
        distal = full.x >= 0.8 * length if side == "left" else full.x <= 0.2 * length

        def regional(trajectory, field, mask):
            return np.mean(trajectory.fields[field][:, mask], axis=1)

        distal_tension = regional(full, "tension", distal)
        distal_rho = regional(full, "rho", distal)
        beta0_distal_rho = regional(beta0, "rho", distal)
        proximal_mca = regional(full, "mca_bound", proximal)
        proximal_rac = regional(full, "rac", proximal)
        alpha0_proximal_rac = regional(alpha0, "rac", proximal)
        pathway[side] = {
            "distal_tension_peak_change": float(
                np.max(distal_tension - distal_tension[0])
            ),
            "distal_rho_peak_change_full": float(
                np.max(distal_rho - distal_rho[0])
            ),
            "distal_rho_peak_change_beta_zero": float(
                np.max(beta0_distal_rho - beta0_distal_rho[0])
            ),
            "proximal_bound_mca_minimum_change": float(
                np.min(proximal_mca - proximal_mca[0])
            ),
            "late_proximal_rac_change_full": float(
                proximal_rac[-1] - proximal_rac[0]
            ),
            "late_proximal_rac_change_alpha_zero": float(
                alpha0_proximal_rac[-1] - alpha0_proximal_rac[0]
            ),
        }
    maximum_reflection_error = max(reflection.values())
    pathway_passed = all(
        values["distal_tension_peak_change"] > 5.0
        and values["distal_rho_peak_change_full"] > 0.03
        and values["distal_rho_peak_change_beta_zero"] < 0.005
        and values["proximal_bound_mca_minimum_change"] < -5e-4
        and values["late_proximal_rac_change_full"] > 0.002
        and abs(values["late_proximal_rac_change_alpha_zero"]) < 5e-4
        for values in pathway.values()
    )
    return {
        "model": "reflection-symmetric extension of the 2026 de Belly model",
        "published_model_modified": False,
        "reflection": {
            "maximum_absolute_error": maximum_reflection_error,
            "field_errors": reflection,
            "passed": maximum_reflection_error <= 3e-13,
        },
        "pathway": {"sides": pathway, "passed": pathway_passed},
        "passed": maximum_reflection_error <= 3e-13 and pathway_passed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=(
            ROOT
            / "results"
            / "validation"
            / "debelly_2026_reflection_symmetric"
            / "validation.json"
        ),
    )
    arguments = parser.parse_args()
    report = build_report()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(arguments.output)
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
