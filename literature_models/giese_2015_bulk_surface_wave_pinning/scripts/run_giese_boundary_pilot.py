#!/usr/bin/env python3
"""Focused early-abort pilot near the Giese abundance/polarity boundary."""

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
sys.path.insert(0, str(ROOT / "scripts"))

from giese_wp import GieseParameters, GieseSolver, make_circular_system  # noqa: E402
from run_giese_competition_screen import (  # noqa: E402
    competition_metrics,
    reaction_balanced_state,
    stable_winner_time,
    stimuli,
)

MASS_SCALES = (0.92, 0.94, 0.96, 0.98, 1.00)
GAMMAS = (0.75, 1.00, 1.25)
DIAMETERS = (10.0, 10.0 * 2.0 ** (1.0 / 3.0))
DM = 0.10
DC = 20.0
WEAK_STIMULUS = 0.30


def design() -> list[dict]:
    return [
        {"mass_scale": mass, "gamma": gamma, "diameter_um": diameter}
        for gamma in GAMMAS
        for mass in MASS_SCALES
        for diameter in DIAMETERS
    ]


def checkpoint_path(output: Path, task_index: int) -> Path:
    return output / "checkpoints" / f"task_{task_index:03d}.npz"


def run_task(output: Path, task_index: int, end_time: float, dt: float,
             target_h: float) -> None:
    condition = design()[task_index]
    base = replace(GieseParameters(), gamma=condition["gamma"])
    u0, v0, target_mass = reaction_balanced_state(
        condition["diameter_um"], condition["mass_scale"], base
    )
    parameters = replace(base, Dm=DM, Dc=DC, u0=u0, v0=v0)
    system = make_circular_system(condition["diameter_um"], target_h=target_h)
    initial = (
        np.full(system.boundary_nodes.size, u0),
        np.full(system.mesh.p.shape[1], v0),
    )
    result = GieseSolver(system, parameters).run(
        end_time=end_time,
        dt=dt,
        stimuli=stimuli("two_site", WEAK_STIMULUS),
        output_interval=1.0,
        initial=initial,
    )
    metrics = {
        key: np.asarray([row[key] for row in result.metrics])
        for key in result.metrics[0]
    }
    competition = competition_metrics(
        result.membrane, system.boundary_angles, system.edge_lengths
    )
    winner_time = stable_winner_time(
        result.times, metrics["n_above_mean_clusters"],
        metrics["relative_amplitude"],
    )
    start_index = int(np.searchsorted(result.times, min(100.0, end_time / 2.0)))
    elapsed = max(result.times[-1] - result.times[start_index], 1.0)
    bias_growth_rate = (
        abs(competition["mass_bias"][-1])
        - abs(competition["mass_bias"][start_index])
    ) / elapsed
    initial_mass = float(result.total_mass[0])
    summary = {
        "task_index": task_index,
        **condition,
        "Dm_um2_s": DM,
        "Dc_um2_s": DC,
        "weak_stimulus_um_s": WEAK_STIMULUS,
        "u0": u0,
        "v0": v0,
        "target_native_mass": target_mass,
        "end_time_s": end_time,
        "winner_selection_time_s": winner_time,
        "winner_selection_right_censored": bool(
            not np.isfinite(winner_time)
            and metrics["relative_amplitude"][-1] >= 0.05
        ),
        "final_domains": int(metrics["n_above_mean_clusters"][-1]),
        "final_relative_amplitude": float(metrics["relative_amplitude"][-1]),
        "mass_bias_100_s": float(competition["mass_bias"][start_index]),
        "mass_bias_final": float(competition["mass_bias"][-1]),
        "absolute_bias_growth_rate_per_s": float(bias_growth_rate),
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
        boundary_angles=system.boundary_angles,
        metric_domains=metrics["n_above_mean_clusters"],
        metric_relative_amplitude=metrics["relative_amplitude"],
        competition_mass_bias=competition["mass_bias"],
    )
    print(json.dumps(summary, indent=2))


def collect(output: Path) -> None:
    rows = []
    for task_index in range(len(design())):
        path = checkpoint_path(output, task_index)
        if not path.exists():
            raise SystemExit(f"Missing {path}")
        with np.load(path, allow_pickle=False) as raw:
            rows.append(json.loads(str(raw["summary_json"])))
    path = output / "boundary_pilot_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"conditions": len(rows), "summary": str(path)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--end-time", type=float, default=500.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--target-h", type=float, default=0.25)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "reports/output/giese_wp/boundary_pilot",
    )
    args = parser.parse_args()
    if args.collect:
        collect(args.output)
        return
    if args.task_index is None or not 0 <= args.task_index < len(design()):
        raise SystemExit(f"Provide --task-index=0..{len(design()) - 1}")
    run_task(args.output, args.task_index, args.end_time, args.dt, args.target_h)


if __name__ == "__main__":
    main()
