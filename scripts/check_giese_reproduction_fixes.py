#!/usr/bin/env python3
"""Compare prioritized explanations for the Giese timing discrepancy at 3 um."""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from giese_wp import (GieseParameters, GieseSolver, Stimulus, half_mass_pf,  # noqa: E402
                      make_circular_system, make_gmsh_circular_system)


def summarize(label: str, system, parameters: GieseParameters,
              scheme: str = "imex") -> dict:
    started = time.time()
    result = GieseSolver(system, parameters).run(
        300.0, dt=0.05,
        stimuli=(Stimulus(0.44, 0.0), Stimulus(0.40, np.pi)),
        output_interval=1.0, scheme=scheme,
    )
    half_pf = np.asarray([
        half_mass_pf(profile, system.edge_lengths) for profile in result.membrane
    ])
    reached = np.flatnonzero(half_pf >= 0.9 * half_pf[-1])
    clusters = np.asarray([
        row["n_above_mean_clusters"] for row in result.metrics
    ])
    amplitude = np.asarray([row["relative_amplitude"] for row in result.metrics])
    unique = (clusters == 1) & (amplitude >= 0.05)
    stable = [index for index, run_time in enumerate(result.times)
              if run_time >= 10.0 and unique[index] and np.all(unique[index:])]
    return {
        "case": label,
        "scheme": scheme,
        "v0": parameters.v0,
        "triangles": system.mesh.t.shape[1],
        "bulk_nodes": system.mesh.p.shape[1],
        "boundary_nodes": system.boundary_nodes.size,
        "final_half_mass_PF": half_pf[-1],
        "half_mass_PF_time_90_final_s": (
            result.times[reached[0]] if reached.size else np.nan
        ),
        "stable_unique_site_time_s": (
            result.times[stable[0]] if stable else np.nan
        ),
        "maximum_relative_mass_error": np.max(np.abs(
            result.total_mass / result.total_mass[0] - 1.0
        )),
        "elapsed_s": time.time() - started,
    }


def main() -> None:
    printed = GieseParameters()
    exact_v = printed.delta * printed.u0 / (
        printed.k0 + printed.gamma * printed.u0**2 /
        (printed.K**2 + printed.u0**2)
    )
    exact = GieseParameters(v0=exact_v)
    rows = [
        summarize("regular_conservative", make_circular_system(3.0, 0.25), printed),
        summarize("regular_giese_split", make_circular_system(3.0, 0.25), printed,
                  scheme="giese_semi_implicit"),
        summarize("gmsh_boundary_refined", make_gmsh_circular_system(
            3.0, boundary_h=0.05, bulk_h=0.18), printed),
        summarize("gmsh_fine_boundary_refined", make_gmsh_circular_system(
            3.0, boundary_h=0.025, bulk_h=0.10), printed),
        summarize("gmsh_exact_equilibrium", make_gmsh_circular_system(
            3.0, boundary_h=0.05, bulk_h=0.18), exact),
    ]
    output = ROOT / "reports/output/giese_wp/fix_diagnostics.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
