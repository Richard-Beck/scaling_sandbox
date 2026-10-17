#!/usr/bin/env python3
"""Derive reversal and response-centroid metrics for the Otsuji geometry sweep."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.contract import Trajectory  # noqa: E402


def _centroid(profile: np.ndarray, relative_x: np.ndarray, *, excess: bool) -> float:
    weights = np.asarray(profile, dtype=float)
    if excess:
        weights = weights - np.min(weights)
    total = float(np.sum(weights))
    return 0.5 if total <= 1e-12 else float(np.sum(relative_x * weights) / total)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_root", type=Path)
    arguments = parser.parse_args()
    manifest = json.loads(arguments.manifest.read_text())
    rows = []
    for task in manifest["tasks"]:
        directory = arguments.output_root / "tasks" / task["task_id"]
        summary = json.loads((directory / "summary.json").read_text())
        trajectory = Trajectory.load(directory / "trajectory.npz")
        relative_x = trajectory.x / trajectory.x[-1]
        reversal_time = float(task["protocol"]["reversal_time"])
        before_index = int(np.flatnonzero(trajectory.time < reversal_time)[-1])
        final_window = trajectory.time >= trajectory.time[-1] - 500.0
        final_profile = np.mean(trajectory.observable[final_window], axis=0)
        pre_profile = trajectory.observable[before_index]
        rows.append(
            {
                "input_peak_per_s": task["input_peak"],
                "stimulus_width_fraction": task["stimulus_width_fraction"],
                "orientation_reversal_status": summary["orientation_reversal_status"],
                "strict_reversal_status": summary["reversal_status"],
                "reversal_latency_s": summary["reversal_latency_s"],
                "pre_orientation": summary["pre_reversal_orientation"],
                "final_orientation": summary["final_orientation"],
                "pre_mass_centroid": _centroid(pre_profile, relative_x, excess=False),
                "final_mass_centroid": _centroid(final_profile, relative_x, excess=False),
                "pre_excess_centroid": _centroid(pre_profile, relative_x, excess=True),
                "final_excess_centroid": _centroid(final_profile, relative_x, excess=True),
                "excess_centroid_shift": (
                    _centroid(final_profile, relative_x, excess=True)
                    - _centroid(pre_profile, relative_x, excess=True)
                ),
                "task_id": task["task_id"],
            }
        )
    destination = arguments.output_root / "reversal_geometry_sweep.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "task_count": len(rows),
        "orientation_reversal_success_count": sum(
            row["orientation_reversal_status"] == "success" for row in rows
        ),
        "strict_reversal_success_count": sum(
            row["strict_reversal_status"] == "success" for row in rows
        ),
    }
    (arguments.output_root / "reversal_geometry_sweep_summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
