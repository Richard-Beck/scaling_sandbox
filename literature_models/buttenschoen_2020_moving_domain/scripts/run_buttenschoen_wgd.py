#!/usr/bin/env python3
"""Run matched 1x/2x-volume pairs for the four Buttenschön Fig. 9 regimes."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from buttenschoen import (  # noqa: E402
    SimulationParameters,
    figure9_regimes,
    scale_wgd_condition,
    simulate,
)
from buttenschoen.render import plot_kymograph, plot_wgd_companion  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--volume-ratio", type=float, default=2.0)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/buttenschoen_spatial/wgd_figure9")
    parser.add_argument("--n-cells", type=int, default=161)
    parser.add_argument("--end-time", type=float, default=720.0)
    parser.add_argument("--output-interval", type=float, default=1.0)
    parser.add_argument("--rtol", type=float, default=2e-6)
    parser.add_argument("--atol", type=float, default=1e-8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.volume_ratio <= 0:
        raise SystemExit("--volume-ratio must be positive")
    simulation = SimulationParameters(
        n_cells=args.n_cells,
        end_time=args.end_time,
        output_interval=args.output_interval,
        rtol=args.rtol,
        atol=args.atol,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    scale = args.volume_ratio**(1.0/3.0)
    design = {
        "experiment": "Buttenschoen Figure 9 matched volume companion",
        "volume_ratio": args.volume_ratio,
        "linear_scale": scale,
        "scaling_rules": {
            "rest_length": "s L0",
            "elastic_modulus": "s E",
            "cell_viscosity": "s eta",
            "endpoint_friction": "s^2 mu",
            "maximum_force": "s^2 Fmax",
            "total_gtpase_amount": "s GT",
        },
        "chemistry_assumption": (
            "Preserve local active/inactive concentrations on the effective "
            "1D support; scale native conserved amount with line length."
        ),
        "unchanged": [
            "biochemical kinetic parameters",
            "physical active and inactive diffusion coefficients",
            "initial local concentrations and relative patch geometry",
        ],
        "initial_centroid_um": 0.0,
        "figure_axis_rule": "Matched physical-position limits within each 1x/2x pair",
    }
    (args.output/"design.json").write_text(json.dumps(design, indent=2)+"\n")
    paired_results = {}
    summary_rows = []
    for name, (baseline_mechanics, baseline_initial) in figure9_regimes().items():
        scaled_mechanics, scaled_initial = scale_wgd_condition(
            baseline_mechanics, baseline_initial, args.volume_ratio
        )
        conditions = {
            "1x": (baseline_mechanics, baseline_initial),
            "2x": (scaled_mechanics, scaled_initial),
        }
        paired_results[name] = {}
        for scale, (mechanics, initial) in conditions.items():
            started = time.perf_counter()
            result = simulate(mechanics=mechanics, initial=initial,
                              simulation=simulation)
            elapsed = time.perf_counter()-started
            directory = args.output/name/scale
            result.save(directory)
            plot_kymograph(
                result, directory/"active_kymograph.png",
                title=f"{name.replace('_', ' ').title()} · {scale}",
            )
            paired_results[name][scale] = result
            contrast = ((np.max(result.active, axis=1)-np.min(result.active, axis=1))
                        / np.maximum(np.max(result.active, axis=1)
                                     + np.min(result.active, axis=1), 1e-15))
            maximum_index = int(np.argmax(result.length))
            depolarized = np.flatnonzero(contrast < 0.1)
            depolarization_index = int(depolarized[0]) if depolarized.size else None
            summary_rows.append({
                "regime": name,
                "scale": scale,
                "initial_centroid_um": float(result.centroid[0]),
                "maximum_length_um": float(result.length[maximum_index]),
                "time_of_maximum_length_s": float(result.times[maximum_index]),
                "depolarization_time_s": (
                    float(result.times[depolarization_index])
                    if depolarization_index is not None else ""
                ),
                "centroid_at_depolarization_um": (
                    float(result.centroid[depolarization_index])
                    if depolarization_index is not None else ""
                ),
                "final_centroid_um": float(result.centroid[-1]),
                "final_left_edge_um": float(result.left_edge[-1]),
                "final_right_edge_um": float(result.right_edge[-1]),
                "final_length_um": float(result.length[-1]),
                "final_active_contrast": float(contrast[-1]),
                "maximum_relative_mass_error": float(
                    np.max(result.relative_mass_error)
                ),
            })
            print(
                f"{name} {scale}: {elapsed:.2f}s, "
                f"final length={result.length[-1]:.3f} um, "
                f"max mass error={result.relative_mass_error.max():.3e}"
            )
    plot_wgd_companion(paired_results, args.output/"wgd_companion_kymographs.png")
    with (args.output/"comparison_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
