#!/usr/bin/env python3
"""Create conditioning and duration-sweep campaigns for Figure 2 panels E-F."""

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
SOURCE_NAME = "debelly_2026_continuous_polarity_alternatives_v1"
CONDITIONING_NAME = "debelly_2026_countercue_threshold_conditioning"
SWEEP_NAME = "debelly_2026_countercue_threshold_sweep"
CONDITIONING_RESULTS = ROOT / "results/campaigns" / CONDITIONING_NAME
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
BETA_VARIANTS = ("beta_0x", "beta_1x", "beta_2x")
DURATIONS = (
    300.0, 360.0, 430.0, 520.0, 630.0,
    750.0, 910.0, 1_090.0, 1_310.0, 1_580.0,
    1_900.0, 2_280.0, 2_750.0, 3_310.0, 3_980.0,
    4_790.0, 5_760.0, 6_930.0, 8_330.0, 10_000.0,
)
CONDITIONING_STOP = 600.0
SETTLING_DEADLINE = 10_600.0
CONFIRMATION = 1_000.0
SETTLING_TOLERANCE = 0.002
FLIP_TOLERANCE = 0.05


def _templates() -> dict[tuple[str, float], dict]:
    source = load_manifest(SOURCE)
    selected = {
        (task["coupling_variant"], float(task["simulation"]["length"])): task
        for task in source["tasks"]
        if task["assay"] == "figure2_B_beta_sweep"
        and task["coupling_variant"] in BETA_VARIANTS
    }
    expected = {(variant, length) for variant in BETA_VARIANTS for length in LENGTHS}
    if selected.keys() != expected:
        raise RuntimeError("source campaign lacks the requested beta/length grid")
    return selected


def _base(source: dict, assay: str) -> dict:
    task = copy.deepcopy(source)
    for key in ("task_id", "task_index", "task_digest"):
        task.pop(key, None)
    task.update(
        tier="debelly_2026_countercue_threshold",
        implementation_revision="countercue_threshold_v1",
        assay=assay,
        source_campaign=SOURCE_NAME,
        source_assay="figure2_B_beta_sweep",
    )
    task["simulation"].update(output_interval=5.0, step=0.5)
    return task


def _conditioning_manifest() -> dict:
    tasks = []
    templates = _templates()
    for variant in BETA_VARIANTS:
        for length in LENGTHS:
            task = _base(templates[(variant, length)], "countercue_threshold_conditioning")
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
                save_trajectory=True,
            )
            task["simulation"]["end_time"] = SETTLING_DEADLINE + CONFIRMATION
            task["metrics"]["analysis_start"] = CONDITIONING_STOP
            tasks.append(task)
    return make_manifest(
        CONDITIONING_NAME,
        tasks,
        description="Conditioning and joint Rac/MCA settling audit for the counter-cue threshold sweep.",
        metadata={
            "beta_variants": list(BETA_VARIANTS),
            "lengths_um": list(LENGTHS),
            "settling_deadline_s": SETTLING_DEADLINE,
            "settling_tolerance_fraction_L": SETTLING_TOLERANCE,
            "settling_confirmation_s": CONFIRMATION,
        },
    )


def _half_mass(payload: dict[str, np.ndarray], field: str) -> np.ndarray:
    values = np.asarray(payload[f"field__{field}"], dtype=float)
    cumulative = np.cumsum(values, axis=1)
    half = 0.5 * cumulative[:, -1]
    crossing = np.argmax(cumulative >= half[:, None], axis=1)
    rows = np.arange(len(values))
    before = np.where(crossing > 0, cumulative[rows, np.maximum(crossing - 1, 0)], 0.0)
    fraction = np.divide(
        half - before, values[rows, crossing], out=np.full_like(half, 0.5),
        where=values[rows, crossing] > 0.0,
    )
    return np.clip((crossing + np.clip(fraction, 0.0, 1.0)) / values.shape[1], 0.0, 1.0)


def _settling(time: np.ndarray, score: np.ndarray) -> tuple[float | None, float]:
    tail = time >= time[-1] - CONFIRMATION
    reference = float(np.median(score[tail]))
    if float(np.max(np.abs(score[tail] - reference))) > SETTLING_TOLERANCE:
        return None, reference
    remains = np.logical_and.accumulate((np.abs(score - reference) <= SETTLING_TOLERANCE)[::-1])[::-1]
    indices = np.flatnonzero(remains & (time >= CONDITIONING_STOP) & (time <= SETTLING_DEADLINE))
    return (None if len(indices) == 0 else float(time[indices[0]])), reference


def _audit() -> list[dict]:
    rows = list(csv.DictReader((CONDITIONING_RESULTS / "task_summary.csv").open(newline="")))
    audited = []
    for row in sorted(rows, key=lambda item: (BETA_VARIANTS.index(item["coupling_variant"]), float(item["length_um"]))):
        path = CONDITIONING_RESULTS / "tasks" / row["task_id"] / "trajectory.npz"
        with np.load(path, allow_pickle=False) as raw:
            payload = {key: raw[key] for key in raw.files}
        time = np.asarray(payload["time"], dtype=float)
        rac_time, rac_ref = _settling(time, _half_mass(payload, "rac"))
        mca_time, mca_ref = _settling(time, _half_mass(payload, "mca_bound"))
        accepted = rac_time is not None and mca_time is not None
        audited.append({
            "task_id": row["task_id"],
            "variant": row["coupling_variant"],
            "length_um": float(row["length_um"]),
            "rac_steady_state_time_s": "" if rac_time is None else rac_time,
            "mca_steady_state_time_s": "" if mca_time is None else mca_time,
            "rac_reference_x50_over_L": rac_ref,
            "mca_reference_x50_over_L": mca_ref,
            "both_reached_before_10000_s": accepted,
            "opposing_cue_start_s": "" if not accepted else max(rac_time, mca_time),
        })
    return audited


def _sweep_manifest() -> tuple[dict, list[dict]]:
    audit = _audit()
    templates = _templates()
    tasks = []
    for row in audit:
        if not row["both_reached_before_10000_s"]:
            continue
        variant = str(row["variant"])
        length = float(row["length_um"])
        start = float(row["opposing_cue_start_s"])
        for duration in DURATIONS:
            task = _base(templates[(variant, length)], "countercue_threshold_sweep")
            task.update(
                input_peak=2.0,
                countercue_duration_s=duration,
                opposing_cue_start_s=start,
                protocol={
                    "name": "finite_pulse_sequence",
                    "pulses": [
                        {"amplitude": 3.6, "start": 0.0, "stop": CONDITIONING_STOP,
                         "taper": 0.0, "width_fraction": 0.2, "sides": ["left"]},
                        {"amplitude": 3.6, "start": start, "stop": start + duration,
                         "taper": 0.0, "width_fraction": 0.2, "sides": ["right"]},
                    ],
                },
                polarity_reversal_audit={
                    "opposite_start": start,
                    "opposite_stop": start + duration,
                    "pre_rac_reference": float(row["rac_reference_x50_over_L"]),
                    "pre_mca_reference": float(row["mca_reference_x50_over_L"]),
                    "settling_tolerance": SETTLING_TOLERANCE,
                    "settling_confirmation_s": CONFIRMATION,
                    "flip_tolerance": FLIP_TOLERANCE,
                },
                save_trajectory=False,
            )
            task["simulation"]["end_time"] = start + duration + 10_000.0
            task["metrics"]["analysis_start"] = start
            tasks.append(task)
    return make_manifest(
        SWEEP_NAME,
        tasks,
        description="Minimum counter-cue duration sweep for steady sign-flipped Rac/MCA repolarization.",
        metadata={
            "beta_variants": list(BETA_VARIANTS),
            "durations_s": list(DURATIONS),
            "accepted_conditioning_cases": sum(row["both_reached_before_10000_s"] for row in audit),
            "discarded_conditioning_cases": sum(not row["both_reached_before_10000_s"] for row in audit),
            "flip_tolerance_fraction_L": FLIP_TOLERANCE,
            "final_followup_s": 10_000.0,
        },
    ), audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("conditioning", "sweep"))
    arguments = parser.parse_args()
    if arguments.stage == "conditioning":
        manifest = _conditioning_manifest()
        destination = ROOT / "campaigns" / f"{CONDITIONING_NAME}.json"
    else:
        manifest, audit = _sweep_manifest()
        destination = ROOT / "campaigns" / f"{SWEEP_NAME}.json"
        with (CONDITIONING_RESULTS / "conditioning_steady_state_audit.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
            writer.writeheader()
            writer.writerows(audit)
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
