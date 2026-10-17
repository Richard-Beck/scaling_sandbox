#!/usr/bin/env python3
"""Print compact numerical diagnostics for a completed Marée reconstruction run."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    with (args.run / "metrics.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit("empty metrics.csv")
    first, last = rows[0], rows[-1]
    pools = {}
    for name, active, inactive in (("C", "C_mass", "Ci_mass"), ("Rac", "R_mass", "Ri_mass"), ("Rho", "rho_mass", "rhoi_mass")):
        p0 = float(first[active]) + float(first[inactive])
        pools[name] = max(abs(float(r[active]) + float(r[inactive]) - p0) / p0 for r in rows)
    areas = [float(r["area_um2"]) for r in rows]
    summary = {
        "time_s": float(last["time_s"]),
        "steps": int(last["step"]),
        "final_sites": int(last["sites"]),
        "area_range_um2": [min(areas), max(areas)],
        "final_displacement_um": float(last["displacement_um"]),
        "mean_speed_nm_s": 1000 * float(last["displacement_um"]) / max(float(last["time_s"]), 1e-300),
        "max_pool_relative_drift": pools,
        "final_C_minus_Rho_x_um": float(last["C_minus_rho_x_um"]),
    }
    snap_path = args.run / "snapshots.csv"
    if snap_path.exists():
        final_step = int(last["step"])
        extrema = {key: [float("inf"), -float("inf")] for key in ("C", "R", "rho", "P_total")}
        count = 0
        with snap_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if int(row["step"]) != final_step:
                    continue
                count += 1
                for key, pair in extrema.items():
                    value = float(row[key]); pair[0] = min(pair[0], value); pair[1] = max(pair[1], value)
        summary["final_snapshot_sites"] = count
        summary["final_field_ranges"] = extrema
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
