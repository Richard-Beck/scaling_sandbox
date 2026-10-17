#!/usr/bin/env python3
"""Run and collect the Giese two-site Dm x Dc x diameter screen."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from giese_wp import GieseParameters, GieseSolver, Stimulus, make_circular_system  # noqa: E402
from giese_wp.metrics import half_mass_pf  # noqa: E402

DM_VALUES = (0.015, 0.03, 0.05, 0.075, 0.10)
DC_VALUES = (3.0, 5.5, 10.0)
DIAMETERS = (3.0, 4.5, 6.0, 8.0, 10.0, 10.0 * 2.0 ** (1.0 / 3.0), 15.0)


def design() -> list[tuple[float, float, float]]:
    return [
        (dm, dc, diameter)
        for dc in DC_VALUES
        for dm in DM_VALUES
        for diameter in DIAMETERS
    ]


def checkpoint_path(output: Path, task_index: int) -> Path:
    return output / "checkpoints" / f"task_{task_index:03d}.npz"


def stimuli() -> tuple[Stimulus, ...]:
    return (
        Stimulus(0.44, 0.0, membrane_fraction=0.05, end_time=10.0),
        Stimulus(0.40, math.pi, membrane_fraction=0.05, end_time=10.0),
    )


def vertex_measure(edge_lengths: np.ndarray) -> np.ndarray:
    return 0.5 * (edge_lengths + np.roll(edge_lengths, 1))


def angular_distance(angle: np.ndarray, center: float) -> np.ndarray:
    return np.angle(np.exp(1j * (angle - center)))


def competition_metrics(membrane: np.ndarray, angles: np.ndarray,
                        edge_lengths: np.ndarray) -> dict[str, np.ndarray]:
    weights = vertex_measure(edge_lengths)
    strong = np.abs(angular_distance(angles, 0.0)) <= math.pi / 2.0
    weak = ~strong
    strong_height = membrane[:, strong].max(axis=1)
    weak_height = membrane[:, weak].max(axis=1)
    strong_mass = membrane[:, strong] @ weights[strong]
    weak_mass = membrane[:, weak] @ weights[weak]
    return {
        "strong_height": strong_height,
        "weak_height": weak_height,
        "height_bias": ((strong_height - weak_height)
                        / np.maximum(strong_height + weak_height, 1e-15)),
        "strong_mass": strong_mass,
        "weak_mass": weak_mass,
        "mass_bias": ((strong_mass - weak_mass)
                      / np.maximum(strong_mass + weak_mass, 1e-15)),
    }


def stable_winner_time(times: np.ndarray, domains: np.ndarray,
                       amplitude: np.ndarray) -> float:
    unique = (domains == 1) & (amplitude >= 0.05) & (times >= 10.0)
    suffix_all = np.logical_and.accumulate(unique[::-1])[::-1]
    reached = np.flatnonzero(suffix_all)
    return float(times[reached[0]]) if reached.size else float("nan")


def run_task(output: Path, task_index: int, target_h: float, dt: float,
             end_time: float) -> Path:
    dm, dc, diameter = design()[task_index]
    system = make_circular_system(diameter, target_h=target_h)
    parameters = replace(GieseParameters(), Dm=dm, Dc=dc)
    result = GieseSolver(system, parameters).run(
        end_time=end_time,
        dt=dt,
        stimuli=stimuli(),
        output_interval=1.0,
    )
    metrics = {
        key: np.asarray([row[key] for row in result.metrics])
        for key in result.metrics[0]
    }
    competition = competition_metrics(
        result.membrane, system.boundary_angles, system.edge_lengths
    )
    initial_mass = float(result.total_mass[0])
    winner_time = stable_winner_time(
        result.times,
        metrics["n_above_mean_clusters"],
        metrics["relative_amplitude"],
    )
    final_half_pf = half_mass_pf(result.membrane[-1], system.edge_lengths)
    summary = {
        "task_index": task_index,
        "Dm_um2_s": dm,
        "Dc_um2_s": dc,
        "diameter_um": diameter,
        "target_h_um": target_h,
        "actual_max_h_um": system.max_edge_length,
        "dt_s": dt,
        "end_time_s": end_time,
        "triangles": int(system.mesh.t.shape[1]),
        "boundary_nodes": int(system.boundary_nodes.size),
        "winner_selection_time_s": winner_time,
        "winner_selection_right_censored": bool(
            not np.isfinite(winner_time)
            and metrics["relative_amplitude"][-1] >= 0.05
        ),
        "final_domains": int(metrics["n_above_mean_clusters"][-1]),
        "final_relative_amplitude": float(metrics["relative_amplitude"][-1]),
        "final_POL": float(metrics["POL"][-1]),
        "final_PF_above_mean": float(metrics["PF"][-1]),
        "final_PF_half_mass": float(final_half_pf),
        "mass_bias_300_s": float(competition["mass_bias"][min(300, result.times.size - 1)]),
        "mass_bias_1000_s": float(competition["mass_bias"][min(1000, result.times.size - 1)]),
        "mass_bias_final": float(competition["mass_bias"][-1]),
        "height_bias_final": float(competition["height_bias"][-1]),
        "initial_native_mass": initial_mass,
        "maximum_relative_mass_error": float(np.max(np.abs(
            result.total_mass / initial_mass - 1.0
        ))),
    }
    path = checkpoint_path(output, task_index)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        summary_json=json.dumps(summary),
        times=result.times,
        membrane=result.membrane,
        total_mass=result.total_mass,
        boundary_angles=system.boundary_angles,
        edge_lengths=system.edge_lengths,
        **{f"metric_{key}": value for key, value in metrics.items()},
        **{f"competition_{key}": value for key, value in competition.items()},
    )
    print(json.dumps(summary, indent=2))
    return path


def collect(output: Path) -> Path:
    rows = []
    missing = []
    for task_index in range(len(design())):
        path = checkpoint_path(output, task_index)
        if not path.exists():
            missing.append(str(path))
            continue
        with np.load(path, allow_pickle=False) as raw:
            rows.append(json.loads(str(raw["summary_json"])))
    if missing:
        raise SystemExit("Missing diffusion-screen checkpoints:\n" + "\n".join(missing))
    rows.sort(key=lambda row: row["task_index"])
    path = output / "diffusion_screen_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "n_conditions": len(rows),
        "Dm_um2_s": DM_VALUES,
        "Dc_um2_s": DC_VALUES,
        "diameters_um": DIAMETERS,
        "summary": str(path),
    }
    (output / "design_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--target-h", type=float, default=0.25)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--end-time", type=float, default=2000.0)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "reports/output/giese_wp/diffusion_screen",
    )
    args = parser.parse_args()
    if args.collect:
        collect(args.output)
        return
    if args.task_index is None:
        raise SystemExit("Provide --task-index=0..104 or --collect")
    if not 0 <= args.task_index < len(design()):
        raise SystemExit(f"task index must be between 0 and {len(design()) - 1}")
    run_task(args.output, args.task_index, args.target_h, args.dt, args.end_time)


if __name__ == "__main__":
    main()
