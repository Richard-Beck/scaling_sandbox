#!/usr/bin/env python3
"""Run or collect the native Giese cell-diameter experiments."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from giese_wp import GieseSolver, Stimulus, make_circular_system  # noqa: E402
from giese_wp.metrics import half_mass_pf, summarize_trajectory  # noqa: E402

DIAMETERS = (1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 15.0)
MODES = ("single", "two_site")


def task_design() -> list[tuple[str, float]]:
    return [(mode, diameter) for mode in MODES for diameter in DIAMETERS]


def checkpoint_path(output: Path, mode: str, diameter: float) -> Path:
    return output / "checkpoints" / f"{mode}_diameter_{diameter:04.1f}.npz"


def stimuli_for(mode: str) -> tuple[Stimulus, ...]:
    strong = Stimulus(0.44, 0.0, membrane_fraction=0.05, end_time=10.0)
    if mode == "single":
        return (strong,)
    return (strong, Stimulus(0.40, np.pi, membrane_fraction=0.05, end_time=10.0))


def run_task(output: Path, mode: str, diameter: float,
             target_h: float, dt: float) -> Path:
    end_time = 300.0 if mode == "single" else 2000.0
    selected_times = (
        (0.0, 10.0, 20.0, 50.0, 100.0, 200.0, 300.0)
        if mode == "single"
        else (0.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 1500.0, 2000.0)
    )
    system = make_circular_system(diameter, target_h=target_h)
    solver = GieseSolver(system)
    result = solver.run(
        end_time=end_time, dt=dt, stimuli=stimuli_for(mode),
        output_interval=1.0, selected_times=selected_times,
    )
    trajectory = summarize_trajectory(result.times, result.metrics)
    initial_mass = result.total_mass[0]
    summary = {
        "mode": mode,
        "diameter_um": diameter,
        "target_h_um": target_h,
        "actual_max_h_um": system.max_edge_length,
        "dt_s": dt,
        "end_time_s": end_time,
        "triangles": int(system.mesh.t.shape[1]),
        "bulk_nodes": int(system.mesh.p.shape[1]),
        "boundary_nodes": int(system.boundary_nodes.size),
        "maximum_relative_mass_error": float(np.max(np.abs(
            result.total_mass / initial_mass - 1.0
        ))),
        **trajectory,
    }
    metrics = {key: np.asarray([row[key] for row in result.metrics])
               for key in result.metrics[0]}
    selected_keys = np.asarray(sorted(result.selected_bulk))
    selected_bulk = np.asarray([result.selected_bulk[key] for key in selected_keys])
    selected_membrane = np.asarray([
        result.selected_membrane[key] for key in selected_keys
    ])
    path = checkpoint_path(output, mode, diameter)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        summary_json=json.dumps(summary),
        times=result.times,
        membrane=result.membrane,
        total_mass=result.total_mass,
        mesh_points=system.mesh.p,
        mesh_triangles=system.mesh.t,
        boundary_nodes=system.boundary_nodes,
        boundary_angles=system.boundary_angles,
        edge_lengths=system.edge_lengths,
        selected_times=selected_keys,
        selected_bulk=selected_bulk,
        selected_membrane=selected_membrane,
        **{f"metric_{key}": value for key, value in metrics.items()},
    )
    print(json.dumps(summary, indent=2))
    return path


def collect(output: Path) -> None:
    paths = [checkpoint_path(output, mode, diameter)
             for mode, diameter in task_design()]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise SystemExit("Missing native checkpoints:\n" + "\n".join(missing))
    summaries = []
    profile_rows = []
    field_manifest = []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            summary = json.loads(str(data["summary_json"]))
            half_pf = np.asarray([
                half_mass_pf(profile, data["edge_lengths"])
                for profile in data["membrane"]
            ])
            half_target = 0.9 * half_pf[-1]
            half_reached = np.flatnonzero(half_pf >= half_target)
            summary["final_half_mass_PF"] = float(half_pf[-1])
            summary["maximal_half_mass_PF"] = float(half_pf.max())
            summary["polarization_time_90_final_half_mass_PF_s"] = (
                float(data["times"][half_reached[0]])
                if half_reached.size else float("nan")
            )
            if summary["mode"] == "two_site":
                times = data["times"]
                clusters = data["metric_n_above_mean_clusters"]
                amplitude = data["metric_relative_amplitude"]
                unique = (clusters == 1) & (amplitude >= 0.05)
                stable = [index for index in range(times.size)
                          if times[index] >= 10.0 and unique[index]
                          and np.all(unique[index:])]
                summary["winner_selection_time_s"] = (
                    float(times[stable[0]]) if stable else float("nan")
                )
                summary["winner_selection_right_censored"] = (
                    not bool(stable) and bool(summary["polarized_final"])
                )
            else:
                summary["winner_selection_time_s"] = float("nan")
                summary["winner_selection_right_censored"] = False
            summaries.append(summary)
            angles = data["boundary_angles"]
            for time_index, time in enumerate(data["times"]):
                for angle, value in zip(angles, data["membrane"][time_index]):
                    profile_rows.append({
                        "mode": summary["mode"],
                        "diameter_um": summary["diameter_um"],
                        "time_s": float(time), "angle_rad": float(angle),
                        "u": float(value),
                    })
            for time in data["selected_times"]:
                field_manifest.append({
                    "mode": summary["mode"],
                    "diameter_um": summary["diameter_um"],
                    "time_s": float(time), "checkpoint": str(path),
                })

    def write_csv(path: Path, rows: list[dict]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(output / "native_size_summary.csv", summaries)
    write_csv(output / "membrane_profiles.csv", profile_rows)
    write_csv(output / "bulk_field_manifest.csv", field_manifest)
    print(f"Collected {len(paths)} checkpoints under {output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-index", type=int,
                        help="zero-based index in the 14-condition design")
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--target-h", type=float, default=0.25)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "reports/output/giese_wp/native")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.collect:
        collect(args.output)
        return
    if args.task_index is None:
        raise SystemExit("Provide --task-index=0..13 or --collect")
    design = task_design()
    if not 0 <= args.task_index < len(design):
        raise SystemExit(f"task index must be between 0 and {len(design) - 1}")
    mode, diameter = design[args.task_index]
    run_task(args.output, mode, diameter, args.target_h, args.dt)


if __name__ == "__main__":
    main()
