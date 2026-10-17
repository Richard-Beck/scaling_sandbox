#!/usr/bin/env python3
"""Tabulate the Otsuji L=10 reversal cue sweep."""

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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_root", type=Path)
    arguments = parser.parse_args()

    manifest = json.loads(arguments.manifest.read_text())
    rows = []
    for task in manifest["tasks"]:
        summary_path = arguments.output_root / "tasks" / task["task_id"] / "summary.json"
        summary = json.loads(summary_path.read_text())
        trajectory = Trajectory.load(summary_path.with_name("trajectory.npz"))
        relative = (trajectory.x - trajectory.x[0]) / (trajectory.x[-1] - trajectory.x[0])
        left = np.mean(trajectory.observable[:, relative <= 0.2], axis=1)
        right = np.mean(trajectory.observable[:, relative >= 0.8], axis=1)
        orientation = (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)
        post_reversal = trajectory.time >= task["conditioning_duration_s"]
        rows.append(
            {
                "input_peak_per_s": task["input_peak"],
                "conditioning_duration_s": task["conditioning_duration_s"],
                "pre_reversal_state": summary["pre_reversal_state"],
                "pre_reversal_domain_count": summary["pre_reversal_domain_count"],
                "pre_reversal_orientation": summary["pre_reversal_orientation"],
                "orientation_reversal_status": summary["orientation_reversal_status"],
                "orientation_reversal_latency_s": summary["orientation_reversal_latency_s"],
                "strict_reversal_status": summary["reversal_status"],
                "strict_reversal_latency_s": summary["reversal_latency_s"],
                "old_front_decay_latency_s": summary["old_front_decay_latency_s"],
                "final_state": summary["endpoint_state"],
                "final_orientation": summary["final_orientation"],
                "minimum_post_reversal_orientation": float(np.min(orientation[post_reversal])),
                "maximum_post_reversal_orientation": float(np.max(orientation[post_reversal])),
            }
        )

    destination = arguments.output_root / "reversal_cue_sweep.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    conditioned = [row for row in rows if row["orientation_reversal_status"] != "not_conditioned"]
    successes = [row for row in rows if row["orientation_reversal_status"] == "success"]
    strict_successes = [row for row in rows if row["strict_reversal_status"] == "success"]
    report = {
        "task_count": len(rows),
        "orientation_conditioned_count": len(conditioned),
        "orientation_reversal_success_count": len(successes),
        "strict_single_front_reversal_success_count": len(strict_successes),
        "successful_conditions": successes,
    }
    report_path = arguments.output_root / "reversal_cue_sweep_summary.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
