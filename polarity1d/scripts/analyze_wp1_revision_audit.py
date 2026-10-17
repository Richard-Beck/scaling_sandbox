#!/usr/bin/env python3
"""Analyze the extended WP1 competition, conditioning, and field audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import Trajectory  # noqa: E402


AUDIT = ROOT / "results" / "campaigns" / "wp1_revision_audit_v1"
OUTPUT = AUDIT / "analysis"


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _number(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    if value in {"", "None", "null", "not_applicable"}:
        return None
    return float(value)


def _boolean(row: dict[str, str], key: str) -> bool:
    return row.get(key, "").lower() == "true"


def _first_sustained(condition: np.ndarray, times: np.ndarray, after: float, duration: float):
    start = int(np.searchsorted(times, after, side="left"))
    for index in range(start, len(times)):
        stop = int(np.searchsorted(times, times[index] + duration, side="left"))
        if stop < len(times) and np.all(condition[index : stop + 1]):
            return float(times[index])
    return None


def _orientation(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    relative = (x - x[0]) / (x[-1] - x[0])
    left = np.mean(values[:, relative <= 0.2], axis=1)
    right = np.mean(values[:, relative >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def competition_summary(rows: list[dict[str, str]]) -> list[dict]:
    selected = [row for row in rows if row["assay"] == "competing_fronts_extended_4x"]
    groups: dict[tuple[str, float], list[dict[str, str]]] = {}
    for row in selected:
        groups.setdefault((row["model"], float(row["length_um"])), []).append(row)
    result = []
    for (model, length), members in sorted(groups.items()):
        latencies = [
            _number(row, "competition_resolution_latency_s") for row in members
            if _number(row, "competition_resolution_latency_s") is not None
        ]
        result.append(
            {
                "model": model, "length_um": length, "n": len(members),
                "observation_end_time_s": _number(members[0], "observation_end_time_s"),
                "resolved_count": sum(row["competition_status"] == "resolved" for row in members),
                "right_censored_count": sum(row["competition_status"] == "right_censored" for row in members),
                "competition_not_established_count": sum(row["competition_status"] == "not_established" for row in members),
                "persistent_multifront_count": sum(
                    _boolean(row, "competition_unresolved_persistent_multifront") for row in members
                ),
                "resolution_latency_median_s": "" if not latencies else float(np.median(latencies)),
            }
        )
    return result


def otsuji_conditioning(rows: list[dict[str, str]]) -> list[dict]:
    result = []
    for row in rows:
        if row["assay"] != "otsuji_gradient_reversal_diagnostic":
            continue
        result.append(
            {
                "length_um": _number(row, "length_um"),
                "reversal_time_s": _number(row, "pre_reversal_time_s"),
                "pre_reversal_state": row.get("pre_reversal_state"),
                "pre_reversal_polarized": row.get("pre_reversal_polarized"),
                "pre_reversal_domain_count": _number(row, "pre_reversal_domain_count"),
                "pre_reversal_spatial_contrast": _number(row, "pre_reversal_spatial_contrast"),
                "pre_reversal_orientation": _number(row, "pre_reversal_orientation"),
                "morphology_conditioned": row.get("reversal_conditioned"),
                "orientation_conditioned": row.get("orientation_reversal_conditioned"),
                "orientation_reversal_status": row.get("orientation_reversal_status"),
                "morphology_reversal_status": row.get("reversal_status"),
                "task_id": row["task_id"],
            }
        )
    return sorted(result, key=lambda row: row["length_um"])


def debelly_fields(rows: list[dict[str, str]]) -> tuple[list[dict], list[dict]]:
    manifest = json.loads((ROOT / "campaigns" / "wp1_revision_audit_v1.json").read_text())
    tasks = {task["task_id"]: task for task in manifest["tasks"]}
    field_rows = []
    joint_rows = []
    thresholds = {"rac": 0.05, "rho": 0.05, "mca_bound": 0.005, "mca_unbound": 0.005}
    for row in rows:
        if row["assay"] != "debelly_intact_reversal_field_audit":
            continue
        task = tasks[row["task_id"]]
        trajectory = Trajectory.load(AUDIT / "tasks" / row["task_id"] / "trajectory.npz")
        reversal_time = float(task["protocol"]["reversal_time"])
        before_index = int(np.flatnonzero(trajectory.time < reversal_time)[-1])
        orientations = {}
        reversal_times = {}
        for field, threshold in thresholds.items():
            values = trajectory.fields[field]
            orientation = _orientation(values, trajectory.x)
            orientations[field] = orientation
            expected_after = orientation > threshold if field in {"rac", "mca_bound", "mca_unbound"} else orientation < -threshold
            event = _first_sustained(expected_after, trajectory.time, reversal_time, 25.0)
            reversal_times[field] = event
            field_rows.append(
                {
                    "input_peak": _number(row, "input_peak"),
                    "length_um": _number(row, "length_um"),
                    "field": field,
                    "pre_orientation": float(orientation[before_index]),
                    "final_orientation": float(orientation[-1]),
                    "post_switch_event_time_s": "" if event is None else event - reversal_time,
                    "event_observed": event is not None,
                    "task_id": row["task_id"],
                }
            )
        joint = (
            (orientations["rac"] > thresholds["rac"])
            & (orientations["rho"] < -thresholds["rho"])
            & (orientations["mca_bound"] > thresholds["mca_bound"])
            & (orientations["mca_unbound"] > thresholds["mca_unbound"])
        )
        joint_time = _first_sustained(joint, trajectory.time, reversal_time, 25.0)
        joint_rows.append(
            {
                "input_peak": _number(row, "input_peak"),
                "length_um": _number(row, "length_um"),
                "rac_only_status": row.get("orientation_reversal_status"),
                "rac_only_latency_s": _number(row, "orientation_reversal_latency_s"),
                "rac_rho_mca_joint_reversal": joint_time is not None,
                "joint_reversal_latency_s": "" if joint_time is None else joint_time - reversal_time,
                "final_joint_state": bool(joint[-1]),
                "task_id": row["task_id"],
            }
        )
    return field_rows, sorted(joint_rows, key=lambda row: (row["input_peak"], row["length_um"]))


def main() -> None:
    rows = _read(AUDIT / "task_summary.csv")
    competition = competition_summary(rows)
    otsuji = otsuji_conditioning(rows)
    fields, joint = debelly_fields(rows)
    _write(OUTPUT / "competition_extended_summary.csv", competition)
    _write(OUTPUT / "otsuji_conditioning_diagnostic.csv", otsuji)
    _write(OUTPUT / "debelly_field_orientation.csv", fields)
    _write(OUTPUT / "debelly_joint_reversal.csv", joint)
    print(OUTPUT)


if __name__ == "__main__":
    main()
