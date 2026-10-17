#!/usr/bin/env python3
"""Run and optionally render a moving Eroume/Maree phase-field cell."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from maree_migration import SimulationConfig, StimulusProtocol, run_simulation  # noqa: E402
from maree_migration.metrics import METRIC_COLUMNS, material_totals, trajectory_metrics  # noqa: E402
from maree_migration.render import render_outputs  # noqa: E402


def save_result(result, output_dir: Path, runtime: float) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / "trajectory.npz", times=result.times,
                       x=result.x, y=result.y, phi=result.phi,
                       concentrations=result.concentrations,
                       material=result.material,
                       exterior_material=result.exterior_material,
                       species=np.asarray(("cdc42_inact", "cdc42_act", "Rac_inact",
                                           "Rac_act", "Rho_inact", "Rho_act",
                                           "P1", "P2", "P3")))
    metrics = trajectory_metrics(result)
    with (output_dir / "trajectory.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(METRIC_COLUMNS)
        writer.writerows(metrics)
    totals = material_totals(result)
    pair_drift = {}
    for label, indices in (("cdc42", (0, 1)), ("rac", (2, 3)), ("rho", (4, 5))):
        total = totals[:, indices].sum(axis=1)
        pair_drift[label] = float(np.max(np.abs(total / total[0] - 1.0)))
    summary = {
        "runtime_seconds": runtime,
        "steps": int(np.ceil(result.config.duration / result.config.dt)),
        "frames": len(result.times),
        "area_relative_range": float(np.ptp(metrics[:, 4]) / metrics[0, 4]),
        "net_displacement": float(np.hypot(metrics[-1, 1] - metrics[0, 1],
                                            metrics[-1, 2] - metrics[0, 2])),
        "final_polarity_magnitude": float(metrics[-1, 7]),
        "gtpase_pair_max_relative_drift": pair_drift,
        "max_exterior_material_fraction_increase": float(np.max(
            result.exterior_material / result.material
            - result.exterior_material[0] / result.material[0])),
        "configuration": asdict(result.config),
        "model_status": {
            "chemistry": "exact transcription of R/vcell_model.cpp",
            "motility": "phenomenological Rac/Cdc42 protrusion minus Rho retraction",
            "scaling": "not implemented",
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("free", "gradient", "rotate", "transient", "legacy"),
                        default="free")
    parser.add_argument("--duration", type=float, default=600.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--gradient-strength", type=float, default=0.05)
    parser.add_argument("--gradient-angle", type=float, default=0.0)
    parser.add_argument("--rotate-at", type=float, default=300.0)
    parser.add_argument("--rotate-angle", type=float, default=90.0)
    parser.add_argument("--transient-duration", type=float, default=20.0)
    parser.add_argument("--transient-amplitude", type=float, default=1.0)
    parser.add_argument("--transient-width", type=float, default=5.0)
    parser.add_argument("--noise-amplitude", type=float, default=0.01,
                        help="set to zero for the deterministic symmetric control")
    parser.add_argument("--dt", type=float, default=0.004)
    parser.add_argument("--dx", type=float, default=1.0)
    parser.add_argument("--box-size", type=float, default=96.0)
    parser.add_argument("--cell-radius", type=float, default=12.5)
    parser.add_argument("--output-interval", type=float, default=2.0)
    parser.add_argument("--freeze-boundary", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--trajectory-overlay", action="store_true")
    parser.add_argument("--fps", type=float, default=15.0)
    parser.add_argument("--pixels", type=int, default=640)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    center = args.box_size / 2.0
    protocol = StimulusProtocol(
        mode=args.mode, gradient_strength=args.gradient_strength,
        gradient_angle=args.gradient_angle, rotate_at=args.rotate_at,
        rotate_angle=args.rotate_angle, pivot_x=center, pivot_y=center,
        cell_center_x=center, cell_center_y=center,
        transient_duration=args.transient_duration,
        transient_amplitude=args.transient_amplitude,
        transient_width=args.transient_width,
        transient_angle=args.gradient_angle,
    )
    config = SimulationConfig(
        duration=args.duration, seed=args.seed, dt=args.dt, dx=args.dx,
        box_size=args.box_size, cell_radius=args.cell_radius,
        output_interval=args.output_interval, noise_amplitude=args.noise_amplitude,
        freeze_boundary=args.freeze_boundary, stimulus=protocol,
    )
    output = args.output or ROOT / "reports/output/maree_migration" / f"{args.mode}_seed{args.seed}"
    last_percent = [-1]

    def progress(step: int, steps: int) -> None:
        percent = int(100 * step / steps)
        bucket = percent // 10
        if bucket > last_percent[0]:
            print(f"{percent:3d}%  t={min(step * config.dt, config.duration):.1f} s",
                  flush=True)
            last_percent[0] = bucket

    started = time.perf_counter()
    result = run_simulation(config, progress=progress)
    runtime = time.perf_counter() - started
    save_result(result, output, runtime)
    if args.render:
        for path in render_outputs(result, output / "rendered", fps=args.fps,
                                   pixels=args.pixels,
                                   trajectory=args.trajectory_overlay):
            print(path)


if __name__ == "__main__":
    main()
