#!/usr/bin/env python3
"""Validate homogeneity, conservation, and discretization convergence."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from giese_wp import GieseSolver, Stimulus, make_circular_system  # noqa: E402
from giese_wp.metrics import summarize_trajectory  # noqa: E402


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_case(target_h: float, dt: float, end_time: float = 60.0) -> dict:
    system = make_circular_system(3.0, target_h=target_h)
    result = GieseSolver(system).run(
        end_time=end_time, dt=dt, stimuli=(Stimulus(0.44, 0.0),),
        output_interval=1.0, selected_times=(end_time,),
    )
    summary = summarize_trajectory(result.times, result.metrics)
    return {
        "target_h_um": target_h,
        "actual_max_h_um": system.max_edge_length,
        "dt_s": dt,
        "triangles": system.mesh.t.shape[1],
        "boundary_nodes": system.boundary_nodes.size,
        "final_PF": summary["final_PF"],
        "final_POL": summary["final_POL"],
        "polarization_time_90_max_PF_s": summary["polarization_time_90_max_PF_s"],
        "relative_mass_error": float(np.max(np.abs(result.total_mass / result.total_mass[0] - 1.0))),
    }


def main() -> None:
    output = ROOT / "reports/output/giese_wp/validation"
    output.mkdir(parents=True, exist_ok=True)

    system = make_circular_system(3.0, target_h=0.25)
    solver = GieseSolver(system)
    homogeneous = solver.run(100.0, dt=0.1, output_interval=1.0)
    relative_mass_error = float(np.max(np.abs(
        homogeneous.total_mass / homogeneous.total_mass[0] - 1.0
    )))
    relative_spatial_range = float(max(
        row["relative_amplitude"] for row in homogeneous.metrics
    ))

    u0, v0 = solver.initial_state()
    reaction = solver.reaction_load(u0, v0, 0.0, (Stimulus(0.44, 0.0),))
    dt = 0.01
    u1, v1, integrated_flux = solver.step(
        u0, v0, 0.0, dt, (Stimulus(0.44, 0.0),)
    )
    membrane_change = float(
        np.ones(u1.size) @ system.membrane_mass @ (u1 - u0)
    )
    bulk_change = float(
        np.ones(v1.size) @ system.bulk_mass @ (v1 - v0)
    )
    flux_check = {
        "integrated_flux": integrated_flux,
        "expected_membrane_change": dt * integrated_flux,
        "membrane_change": membrane_change,
        "bulk_change": bulk_change,
        "opposite_sign_residual": membrane_change + bulk_change,
    }
    unit_bulk = np.ones_like(v0)
    coverage_load = solver.stimulus_reaction_load(
        unit_bulk, Stimulus(1.0, 0.0, membrane_fraction=0.05), 0.0
    )
    stimulus_fraction = float(coverage_load.sum() / system.perimeter)

    mesh_rows = [run_case(h, 0.05) for h in (0.5, 0.3, 0.15)]
    time_rows = [run_case(0.25, step) for step in (0.2, 0.1, 0.05)]
    write_rows(output / "mesh_convergence.csv", mesh_rows)
    write_rows(output / "time_convergence.csv", time_rows)

    validation = {
        "homogeneous_no_stimulus": {
            "duration_s": 100.0,
            "maximum_relative_spatial_range": relative_spatial_range,
            "maximum_relative_mass_error": relative_mass_error,
        },
        "shared_flux_one_step": flux_check,
        "stimulus_coverage": {
            "requested_membrane_fraction": 0.05,
            "integrated_membrane_fraction": stimulus_fraction,
        },
        "thresholds": {
            "maximum_relative_mass_error": 1e-10,
            "maximum_relative_spatial_range": 1e-4,
            "maximum_shared_flux_residual": 1e-12,
            "maximum_stimulus_fraction_error": 1e-12,
        },
    }
    validation["passed"] = bool(
        relative_mass_error <= 1e-10
        and relative_spatial_range <= 1e-4
        and abs(flux_check["opposite_sign_residual"]) <= 1e-12
        and abs(stimulus_fraction - 0.05) <= 1e-12
    )
    (output / "validation_summary.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )
    print(json.dumps(validation, indent=2))
    if not validation["passed"]:
        raise SystemExit("Giese validation thresholds were not met")


if __name__ == "__main__":
    main()
