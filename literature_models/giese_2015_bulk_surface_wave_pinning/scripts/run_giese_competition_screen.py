#!/usr/bin/env python3
"""Run and collect the coupled-diffusion, mass, and stimulus Giese screen."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from giese_wp import GieseParameters, GieseSolver, Stimulus, make_circular_system  # noqa: E402
from giese_wp.metrics import half_mass_pf  # noqa: E402

DIFFUSION_PAIRS = ((0.05, 10.0), (0.075, 15.0), (0.10, 20.0))
MASS_SCALES = (0.7, 0.8, 0.9, 1.0)
WEAK_STIMULI = (0.40, 0.35, 0.30)
DIAMETERS = (10.0, 10.0 * 2.0 ** (1.0 / 3.0))


def design() -> list[dict]:
    rows = []
    for dm, dc in DIFFUSION_PAIRS:
        for mass_scale in MASS_SCALES:
            for weak_stimulus in WEAK_STIMULI:
                for diameter in DIAMETERS:
                    rows.append({
                        "mode": "two_site", "Dm_um2_s": dm, "Dc_um2_s": dc,
                        "mass_scale": mass_scale, "weak_stimulus_um_s": weak_stimulus,
                        "diameter_um": diameter,
                    })
    # The single-site assay is independent of the absent weak stimulus, so run
    # each diffusion x mass x diameter condition only once.
    for dm, dc in DIFFUSION_PAIRS:
        for mass_scale in MASS_SCALES:
            for diameter in DIAMETERS:
                rows.append({
                    "mode": "single", "Dm_um2_s": dm, "Dc_um2_s": dc,
                    "mass_scale": mass_scale, "weak_stimulus_um_s": None,
                    "diameter_um": diameter,
                })
    return rows


def checkpoint_path(output: Path, task_index: int) -> Path:
    return output / "checkpoints" / f"task_{task_index:03d}.npz"


def reaction_balanced_state(diameter: float, mass_scale: float,
                            parameters: GieseParameters) -> tuple[float, float, float]:
    perimeter = math.pi * diameter
    area = math.pi * (diameter / 2.0) ** 2
    baseline_mass = perimeter * parameters.u0 + area * parameters.v0
    target_mass = mass_scale * baseline_mass

    def equilibrium_v(u: float) -> float:
        activation = (parameters.k0 + parameters.gamma * u**2
                      / (parameters.K**2 + u**2))
        return parameters.delta * u / activation

    def residual(u: float) -> float:
        return perimeter * u + area * equilibrium_v(u) - target_mass

    grid = np.geomspace(1e-10, 10.0, 20000)
    values = np.asarray([residual(value) for value in grid])
    roots = []
    for left, right, f_left, f_right in zip(
        grid[:-1], grid[1:], values[:-1], values[1:]
    ):
        if f_left * f_right < 0.0:
            root = brentq(residual, left, right)
            if not roots or not math.isclose(root, roots[-1], rel_tol=1e-8):
                roots.append(root)
    if not roots:
        raise RuntimeError(
            f"No positive homogeneous equilibrium for d={diameter}, mass={mass_scale}"
        )
    # Follow the low-u branch connected to the published homogeneous state.
    u0 = min(roots)
    v0 = equilibrium_v(u0)
    return float(u0), float(v0), float(target_mass)


def stimuli(mode: str, weak_stimulus: float | None) -> tuple[Stimulus, ...]:
    strong = Stimulus(0.44, 0.0, membrane_fraction=0.05, end_time=10.0)
    if mode == "single":
        return (strong,)
    return (
        strong,
        Stimulus(float(weak_stimulus), math.pi,
                 membrane_fraction=0.05, end_time=10.0),
    )


def competition_metrics(membrane: np.ndarray, angles: np.ndarray,
                        edge_lengths: np.ndarray) -> dict[str, np.ndarray]:
    weights = 0.5 * (edge_lengths + np.roll(edge_lengths, 1))
    distance = np.angle(np.exp(1j * angles))
    strong = np.abs(distance) <= math.pi / 2.0
    weak = ~strong
    strong_height = membrane[:, strong].max(axis=1)
    weak_height = membrane[:, weak].max(axis=1)
    strong_mass = membrane[:, strong] @ weights[strong]
    weak_mass = membrane[:, weak] @ weights[weak]
    return {
        "height_bias": ((strong_height - weak_height)
                        / np.maximum(strong_height + weak_height, 1e-15)),
        "mass_bias": ((strong_mass - weak_mass)
                      / np.maximum(strong_mass + weak_mass, 1e-15)),
    }


def stable_winner_time(times: np.ndarray, domains: np.ndarray,
                       amplitude: np.ndarray) -> float:
    unique = (times >= 10.0) & (domains == 1) & (amplitude >= 0.05)
    stable = np.logical_and.accumulate(unique[::-1])[::-1]
    reached = np.flatnonzero(stable)
    return float(times[reached[0]]) if reached.size else float("nan")


def threshold_time(times: np.ndarray, values: np.ndarray, target: float) -> float:
    reached = np.flatnonzero(values >= target - 1e-12)
    return float(times[reached[0]]) if reached.size else float("nan")


def run_task(output: Path, task_index: int, target_h: float, dt: float,
             end_time: float) -> Path:
    condition = design()[task_index]
    base = GieseParameters()
    u0, v0, target_mass = reaction_balanced_state(
        condition["diameter_um"], condition["mass_scale"], base
    )
    parameters = replace(
        base, Dm=condition["Dm_um2_s"], Dc=condition["Dc_um2_s"], u0=u0, v0=v0
    )
    system = make_circular_system(condition["diameter_um"], target_h=target_h)
    initial = (
        np.full(system.boundary_nodes.size, u0),
        np.full(system.mesh.p.shape[1], v0),
    )
    result = GieseSolver(system, parameters).run(
        end_time=end_time,
        dt=dt,
        stimuli=stimuli(condition["mode"], condition["weak_stimulus_um_s"]),
        output_interval=1.0,
        initial=initial,
    )
    metrics = {
        key: np.asarray([row[key] for row in result.metrics])
        for key in result.metrics[0]
    }
    half_pf = np.asarray([
        half_mass_pf(profile, system.edge_lengths) for profile in result.membrane
    ])
    competition = competition_metrics(
        result.membrane, system.boundary_angles, system.edge_lengths
    )
    winner_time = stable_winner_time(
        result.times, metrics["n_above_mean_clusters"], metrics["relative_amplitude"]
    )
    max_half_pf = float(np.max(half_pf))
    final_half_pf = float(half_pf[-1])
    time_90_max = threshold_time(result.times, half_pf, 0.9 * max_half_pf)
    time_90_final = threshold_time(result.times, half_pf, 0.9 * final_half_pf)
    initial_mass = float(result.total_mass[0])
    summary = {
        "task_index": task_index,
        **condition,
        "weak_stimulus_um_s": (condition["weak_stimulus_um_s"]
                                 if condition["weak_stimulus_um_s"] is not None
                                 else float("nan")),
        "u0": u0,
        "v0": v0,
        "target_native_mass": target_mass,
        "actual_initial_mass": initial_mass,
        "target_h_um": target_h,
        "actual_max_h_um": system.max_edge_length,
        "dt_s": dt,
        "end_time_s": end_time,
        "triangles": int(system.mesh.t.shape[1]),
        "boundary_nodes": int(system.boundary_nodes.size),
        "winner_selection_time_s": winner_time,
        "winner_selection_right_censored": bool(
            not np.isfinite(winner_time) and metrics["relative_amplitude"][-1] >= 0.05
        ),
        "final_domains": int(metrics["n_above_mean_clusters"][-1]),
        "final_relative_amplitude": float(metrics["relative_amplitude"][-1]),
        "final_POL": float(metrics["POL"][-1]),
        "final_half_mass_PF": final_half_pf,
        "maximal_half_mass_PF": max_half_pf,
        "polarization_time_90_max_half_mass_PF_s": time_90_max,
        "polarization_time_90_final_half_mass_PF_s": time_90_final,
        "mass_bias_300_s": float(competition["mass_bias"][min(300, result.times.size - 1)]),
        "mass_bias_1000_s": float(competition["mass_bias"][min(1000, result.times.size - 1)]),
        "mass_bias_final": float(competition["mass_bias"][-1]),
        "height_bias_final": float(competition["height_bias"][-1]),
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
        metric_half_mass_PF=half_pf,
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
        raise SystemExit("Missing competition-screen checkpoints:\n" + "\n".join(missing))
    rows.sort(key=lambda row: row["task_index"])
    path = output / "competition_screen_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "n_conditions": len(rows),
        "n_two_site": sum(row["mode"] == "two_site" for row in rows),
        "n_single": sum(row["mode"] == "single" for row in rows),
        "diffusion_pairs": DIFFUSION_PAIRS,
        "mass_scales": MASS_SCALES,
        "weak_stimuli_um_s": WEAK_STIMULI,
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
        default=ROOT / "reports/output/giese_wp/competition_screen",
    )
    args = parser.parse_args()
    if args.collect:
        collect(args.output)
        return
    if args.task_index is None:
        raise SystemExit("Provide --task-index=0..95 or --collect")
    if not 0 <= args.task_index < len(design()):
        raise SystemExit(f"task index must be between 0 and {len(design()) - 1}")
    run_task(args.output, args.task_index, args.target_h, args.dt, args.end_time)


if __name__ == "__main__":
    main()
