#!/usr/bin/env python3
"""Create conditioning and opposing-cue manifests for the Figure 2 trial panel."""

from __future__ import annotations

import argparse
import copy
import csv
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


SOURCE = ROOT / "campaigns/debelly_2026_continuous_polarity_alternatives_v1.json"
CONDITIONING_NAME = "debelly_2026_opposing_cue_duration_conditioning_trial"
RESPONSE_NAME = "debelly_2026_opposing_cue_duration_response_trial"
DOUBLE_RESPONSE_NAME = "debelly_2026_opposing_cue_duration_double_strength_trial"
CONDITIONING_RESULTS = ROOT / "results/campaigns" / CONDITIONING_NAME
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
BASE_DURATIONS = (60.0, 300.0, 1_200.0)
EXTENDED_DURATIONS = (7_200.0, 86_400.0)
DURATIONS = (*BASE_DURATIONS, *EXTENDED_DURATIONS)
CONDITIONING_STOP = 600.0
SETTLING_DEADLINE = CONDITIONING_STOP + 10_000.0
SETTLING_CONFIRMATION_S = 1_000.0
SETTLING_TOLERANCE = 0.002


def _templates() -> dict[float, dict]:
    manifest = load_manifest(SOURCE)
    selected = {
        float(task["simulation"]["length"]): task
        for task in manifest["tasks"]
        if task["assay"] == "figure2_B_beta_sweep"
        and task["coupling_variant"] == "beta_1x"
    }
    if selected.keys() != set(LENGTHS):
        raise RuntimeError("source campaign lacks the full beta=1x length grid")
    return selected


def _base_task(source: dict, assay: str) -> dict:
    task = copy.deepcopy(source)
    for key in ("task_id", "task_index", "task_digest"):
        task.pop(key, None)
    task.update(
        tier="debelly_2026_opposing_cue_duration_trial",
        implementation_revision="opposing_cue_duration_trial_v1",
        assay=assay,
        coupling_variant="beta_1x",
        source_campaign=load_manifest(SOURCE)["campaign_name"],
        source_assay="figure2_B_beta_sweep",
        save_trajectory=True,
    )
    task["simulation"].update(output_interval=5.0, step=0.5)
    return task


def _conditioning_manifest() -> dict:
    tasks = []
    for length, source in _templates().items():
        task = _base_task(source, "opposing_cue_duration_conditioning")
        task.update(
            input_peak=2.0,
            protocol={
                "name": "localized_pulse",
                "amplitude": 3.6,
                "start": 0.0,
                "stop": CONDITIONING_STOP,
                "taper": 0.0,
                "width_fraction": 0.2,
                "sides": ["left"],
                "asymmetry": 0.0,
            },
        )
        task["simulation"]["end_time"] = SETTLING_DEADLINE + SETTLING_CONFIRMATION_S
        task["metrics"]["analysis_start"] = CONDITIONING_STOP
        tasks.append(task)
    return make_manifest(
        CONDITIONING_NAME,
        tasks,
        description=(
            "Conditioning-only stage for the opposing-cue duration trial: a 600-s "
            "left cue in the intact beta=1x model, followed through the settling audit."
        ),
        metadata={
            "lengths_um": list(LENGTHS),
            "conditioning_stop_s": CONDITIONING_STOP,
            "settling_deadline_s": SETTLING_DEADLINE,
            "settling_confirmation_s": SETTLING_CONFIRMATION_S,
            "settling_tolerance_fraction_L": SETTLING_TOLERANCE,
        },
    )


def _half_mass(payload: dict[str, np.ndarray], field: str) -> np.ndarray:
    values = np.asarray(payload[f"field__{field}"], dtype=float)
    cumulative = np.cumsum(values, axis=1)
    half = 0.5 * cumulative[:, -1]
    crossing = np.argmax(cumulative >= half[:, None], axis=1)
    indices = np.arange(len(values))
    before = np.where(
        crossing > 0,
        cumulative[indices, np.maximum(crossing - 1, 0)],
        0.0,
    )
    fraction = np.divide(
        half - before,
        values[indices, crossing],
        out=np.full_like(half, 0.5),
        where=values[indices, crossing] > 0.0,
    )
    return np.clip((crossing + np.clip(fraction, 0.0, 1.0)) / values.shape[1], 0.0, 1.0)


def _settling_time(time: np.ndarray, score: np.ndarray) -> float | None:
    tail = time >= time[-1] - SETTLING_CONFIRMATION_S
    reference = float(np.median(score[tail]))
    if float(np.max(np.abs(score[tail] - reference))) > SETTLING_TOLERANCE:
        return None
    remains = np.logical_and.accumulate((np.abs(score - reference) <= SETTLING_TOLERANCE)[::-1])[::-1]
    eligible = remains & (time >= CONDITIONING_STOP) & (time <= SETTLING_DEADLINE)
    indices = np.flatnonzero(eligible)
    return None if len(indices) == 0 else float(time[indices[0]])


def _conditioning_audit() -> list[dict]:
    rows = list(csv.DictReader((CONDITIONING_RESULTS / "task_summary.csv").open(newline="")))
    audited = []
    for row in sorted(rows, key=lambda item: float(item["length_um"])):
        path = CONDITIONING_RESULTS / "tasks" / row["task_id"] / "trajectory.npz"
        with np.load(path, allow_pickle=False) as raw:
            payload = {key: raw[key] for key in raw.files}
        time = np.asarray(payload["time"], dtype=float)
        rac = _settling_time(time, _half_mass(payload, "rac"))
        mca = _settling_time(time, _half_mass(payload, "mca_bound"))
        accepted = rac is not None and mca is not None
        audited.append({
            "task_id": row["task_id"],
            "length_um": float(row["length_um"]),
            "rac_steady_state_time_s": "" if rac is None else rac,
            "mca_steady_state_time_s": "" if mca is None else mca,
            "both_reached_before_10000_s": accepted,
            "opposing_cue_start_s": "" if not accepted else max(rac, mca),
        })
    return audited


def _response_manifest(
    *,
    name: str = RESPONSE_NAME,
    durations: tuple[float, ...] = DURATIONS,
    opposing_amplitude: float = 3.6,
) -> tuple[dict, list[dict]]:
    audit = _conditioning_audit()
    templates = _templates()
    tasks = []
    accepted = [row for row in audit if row["both_reached_before_10000_s"]]

    def add_task(row: dict, duration: float) -> None:
        length = float(row["length_um"])
        start = float(row["opposing_cue_start_s"])
        task = _base_task(templates[length], "opposing_cue_duration_response")
        task.update(
            input_peak=opposing_amplitude / 1.8,
            opposing_cue_duration_s=duration,
            opposing_cue_start_s=start,
            conditioning_rac_steady_state_time_s=float(row["rac_steady_state_time_s"]),
            conditioning_mca_steady_state_time_s=float(row["mca_steady_state_time_s"]),
            protocol={
                "name": "finite_pulse_sequence",
                "pulses": [
                    {
                        "amplitude": 3.6,
                        "start": 0.0,
                        "stop": CONDITIONING_STOP,
                        "taper": 0.0,
                        "width_fraction": 0.2,
                        "sides": ["left"],
                    },
                    {
                        "amplitude": opposing_amplitude,
                        "start": start,
                        "stop": start + duration,
                        "taper": 0.0,
                        "width_fraction": 0.2,
                        "sides": ["right"],
                    },
                ],
            },
        )
        task["simulation"]["end_time"] = start + duration + 10_000.0
        task["metrics"]["analysis_start"] = start
        tasks.append(task)

    if durations == DURATIONS and np.isclose(opposing_amplitude, 3.6):
        # Retain the original task order and identifiers, then append new durations.
        for row in accepted:
            for duration in BASE_DURATIONS:
                add_task(row, duration)
        for duration in EXTENDED_DURATIONS:
            for row in accepted:
                add_task(row, duration)
    else:
        for row in accepted:
            for duration in durations:
                add_task(row, duration)
    manifest = make_manifest(
        name,
        tasks,
        description=(
            "Opposing-cue duration trial after length-specific Rac and bound-MCA "
            f"steady state: right-hand cues at amplitude {opposing_amplitude:g} "
            "in the beta=1x full model."
        ),
        metadata={
            "candidate_lengths_um": list(LENGTHS),
            "accepted_lengths_um": [row["length_um"] for row in audit if row["both_reached_before_10000_s"]],
            "discarded_lengths_um": [row["length_um"] for row in audit if not row["both_reached_before_10000_s"]],
            "opposing_cue_durations_s": list(durations),
            "opposing_cue_amplitude": opposing_amplitude,
            "opposing_cue_strength_over_alpha0": opposing_amplitude / 1.8,
            "settling_rule": (
                "both Rac and bound-MCA x50/L enter and remain within 0.002 of their "
                "final 1000-s median, no later than 10000 s after conditioning withdrawal"
            ),
        },
    )
    return manifest, audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("conditioning", "response", "response_double"))
    arguments = parser.parse_args()
    if arguments.stage == "conditioning":
        manifest = _conditioning_manifest()
        destination = ROOT / "campaigns" / f"{CONDITIONING_NAME}.json"
    elif arguments.stage == "response":
        manifest, audit = _response_manifest()
        destination = ROOT / "campaigns" / f"{RESPONSE_NAME}.json"
        audit_path = CONDITIONING_RESULTS / "conditioning_steady_state_audit.csv"
        with audit_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
            writer.writeheader()
            writer.writerows(audit)
    else:
        manifest, audit = _response_manifest(
            name=DOUBLE_RESPONSE_NAME,
            durations=BASE_DURATIONS,
            opposing_amplitude=7.2,
        )
        destination = ROOT / "campaigns" / f"{DOUBLE_RESPONSE_NAME}.json"
        audit_path = CONDITIONING_RESULTS / "conditioning_steady_state_audit.csv"
        with audit_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
            writer.writeheader()
            writer.writerows(audit)
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
