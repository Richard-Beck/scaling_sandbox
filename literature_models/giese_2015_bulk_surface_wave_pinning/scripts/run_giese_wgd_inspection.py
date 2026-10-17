#!/usr/bin/env python3
"""Run the long two-site Giese pair used for WGD-focused inspection figures."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from giese_wp import GieseParameters, GieseSolver, Stimulus, make_circular_system  # noqa: E402
from giese_wp.metrics import summarize_trajectory  # noqa: E402

BASELINE_DIAMETER_UM = 10.0
LINEAR_WGD_SCALE = 2.0 ** (1.0 / 3.0)
DIAMETERS = (
    BASELINE_DIAMETER_UM,
    BASELINE_DIAMETER_UM * LINEAR_WGD_SCALE,
    4.5 * LINEAR_WGD_SCALE,
    6.0 * LINEAR_WGD_SCALE,
)
OUTPUT_LABELS = ("1x", "2x", "4p5_to_2x", "6_to_2x")


def stimuli() -> tuple[Stimulus, ...]:
    """Paper's two opposite, unequal, transient membrane stimuli."""
    return (
        Stimulus(0.44, 0.0, membrane_fraction=0.05, end_time=10.0),
        Stimulus(0.40, math.pi, membrane_fraction=0.05, end_time=10.0),
    )


def output_path(output: Path, task_index: int) -> Path:
    label = OUTPUT_LABELS[task_index]
    return output / "checkpoints" / f"two_site_wgd_{label}.npz"


def run(output: Path, diameter: float, target_h: float, dt: float,
        end_time: float, task_index: int) -> Path:
    system = make_circular_system(diameter, target_h=target_h)
    parameters = GieseParameters()
    selected_times = (0.0, 10.0, 50.0, 100.0, 300.0, 500.0, 1000.0,
                      1500.0, 2000.0, end_time)
    result = GieseSolver(system, parameters).run(
        end_time=end_time,
        dt=dt,
        stimuli=stimuli(),
        output_interval=1.0,
        selected_times=selected_times,
    )
    trajectory = summarize_trajectory(result.times, result.metrics)
    metrics = {
        key: np.asarray([row[key] for row in result.metrics])
        for key in result.metrics[0]
    }
    selected_keys = np.asarray(sorted(result.selected_bulk))
    initial_mass = float(result.total_mass[0])
    pair_baseline = (BASELINE_DIAMETER_UM, BASELINE_DIAMETER_UM, 4.5, 6.0)[task_index]
    metadata = {
        "diameter_um": diameter,
        "pair_baseline_diameter_um": pair_baseline,
        "volume_scale": (diameter / pair_baseline) ** 3,
        "linear_scale": diameter / pair_baseline,
        "target_h_um": target_h,
        "actual_max_h_um": system.max_edge_length,
        "dt_s": dt,
        "end_time_s": end_time,
        "triangles": int(system.mesh.t.shape[1]),
        "bulk_nodes": int(system.mesh.p.shape[1]),
        "boundary_nodes": int(system.boundary_nodes.size),
        "initial_native_mass": initial_mass,
        "maximum_relative_mass_error": float(np.max(np.abs(
            result.total_mass / initial_mass - 1.0
        ))),
        "parameters": {
            "k0": parameters.k0,
            "gamma": parameters.gamma,
            "K": parameters.K,
            "delta": parameters.delta,
            "Dm_um2_s": parameters.Dm,
            "Dc_um2_s": parameters.Dc,
            "u0_um_uM": parameters.u0,
            "v0_uM": parameters.v0,
        },
        "stimuli": [
            {
                "strength_um_s": stimulus.strength,
                "center_angle_rad": stimulus.center_angle,
                "membrane_fraction": stimulus.membrane_fraction,
                "start_time_s": stimulus.start_time,
                "end_time_s": stimulus.end_time,
            }
            for stimulus in stimuli()
        ],
        **trajectory,
    }
    path = output_path(output, task_index)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        metadata_json=json.dumps(metadata),
        times=result.times,
        membrane=result.membrane,
        total_mass=result.total_mass,
        mesh_points=system.mesh.p,
        mesh_triangles=system.mesh.t,
        boundary_nodes=system.boundary_nodes,
        boundary_angles=system.boundary_angles,
        edge_lengths=system.edge_lengths,
        selected_times=selected_keys,
        selected_bulk=np.asarray([result.selected_bulk[t] for t in selected_keys]),
        selected_membrane=np.asarray([
            result.selected_membrane[t] for t in selected_keys
        ]),
        **{f"metric_{key}": value for key, value in metrics.items()},
    )
    print(json.dumps(metadata, indent=2))
    print(path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-index", type=int, required=True,
                        choices=tuple(range(len(DIAMETERS))))
    parser.add_argument("--target-h", type=float, default=0.25)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--end-time", type=float, default=2000.0)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "reports/output/giese_wp/wgd_inspection",
    )
    args = parser.parse_args()
    run(args.output, DIAMETERS[args.task_index], args.target_h, args.dt,
        args.end_time, args.task_index)


if __name__ == "__main__":
    main()
