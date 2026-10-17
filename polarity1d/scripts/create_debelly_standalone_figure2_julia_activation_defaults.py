#!/usr/bin/env python3
"""Create Figure 2 campaigns with only alpha0, beta0, alpha, and beta reset."""

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
BASE_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_base"
RESPONSE_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response"
BASE_RESULTS = ROOT / "results/campaigns" / BASE_NAME
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
MODEL_ARMS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
BETA_ARMS = ("beta_0x", "beta_1x", "beta_2x")
BETA_VALUES = {"beta_0x": 0.0, "beta_1x": 0.2, "beta_2x": 0.4}
TRAJECTORY_DURATIONS = (150.0, 300.0, 600.0)
SWEEP_DURATIONS = tuple(float(value) for value in np.linspace(10.0, 1_000.0, 50))
CONDITIONING_STOP = 600.0
SETTLING_DEADLINE = 10_600.0
CONFIRMATION = 1_000.0
SETTLING_TOLERANCE = 0.002
FLIP_TOLERANCE = 0.05


def _copy(source: dict) -> dict:
    task = copy.deepcopy(source)
    for key in ("task_id", "task_index", "task_digest"):
        task.pop(key, None)
    task.update(
        tier="debelly_2026_standalone_figure2_julia_activation_defaults",
        implementation_revision="julia_activation_defaults_v1",
        source_campaign="debelly_2026_continuous_polarity_alternatives_v1",
    )
    parameters = task["model_parameters"]
    parameters.update(rac_basal_activation=1.0, rho_basal_activation=1.0)
    return task


def _set_arm(task: dict, variant: str) -> None:
    parameters = task["model_parameters"]
    parameters["rac_mechanical_feedback"] = 0.0 if variant in {"alpha_zero", "local_inhibition"} else 0.2
    if variant in BETA_VALUES:
        parameters["rho_mechanical_feedback"] = BETA_VALUES[variant]
    else:
        parameters["rho_mechanical_feedback"] = 0.0 if variant in {"beta_zero", "local_inhibition"} else 0.2


def _base_manifest() -> dict:
    source = load_manifest(SOURCE)
    tasks = []
    singles = [
        task for task in source["tasks"]
        if task["assay"] == "single_stimulus_length_test"
        and task["coupling_variant"] in MODEL_ARMS
    ]
    competitions = [
        task for task in source["tasks"]
        if task["assay"] == "figure2_B_beta_sweep"
        and task["coupling_variant"] in BETA_ARMS
    ]
    if len(singles) != 80 or len(competitions) != 60:
        raise RuntimeError("source campaign lacks the required Figure 2 grids")

    for source_task in singles:
        task = _copy(source_task)
        _set_arm(task, task["coupling_variant"])
        task["source_assay"] = "single_stimulus_length_test"
        task["save_trajectory"] = True
        tasks.append(task)

    competition_templates = {}
    for source_task in competitions:
        task = _copy(source_task)
        variant = task["coupling_variant"]
        _set_arm(task, variant)
        task.update(
            beta_multiplier={"beta_0x": 0.0, "beta_1x": 1.0, "beta_2x": 2.0}[variant],
            beta_tension_feedback=BETA_VALUES[variant],
            source_assay="figure2_B_beta_sweep",
            save_trajectory=True,
        )
        tasks.append(task)
        competition_templates[(variant, float(task["simulation"]["length"]))] = task

    for variant in BETA_ARMS:
        for length in LENGTHS:
            task = _copy(competition_templates[(variant, length)])
            task.update(
                assay="countercue_threshold_conditioning",
                input_peak=2.0,
                source_assay="figure2_B_beta_sweep",
                protocol={
                    "name": "localized_pulse", "amplitude": 3.6,
                    "start": 0.0, "stop": CONDITIONING_STOP, "taper": 0.0,
                    "width_fraction": 0.2, "sides": ["left"], "asymmetry": 0.0,
                },
                save_trajectory=True,
            )
            task["simulation"].update(end_time=SETTLING_DEADLINE + CONFIRMATION, output_interval=5.0, step=0.5)
            task["metrics"]["analysis_start"] = CONDITIONING_STOP
            tasks.append(task)
    return make_manifest(
        BASE_NAME,
        tasks,
        description=(
            "Standalone Figure 2 base assays with only alpha0=beta0=1 and "
            "alpha=beta(default)=0.2 reset to the general 1D Julia values."
        ),
        metadata={
            "changed_parameters_only": {
                "rac_basal_activation": 1.0,
                "rho_basal_activation": 1.0,
                "rac_mechanical_feedback": 0.2,
                "rho_mechanical_feedback_0x_1x_2x": [0.0, 0.2, 0.4],
            },
            "all_other_parameters_and_protocols": "unchanged from Figure2_standalone_v2",
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
    rows = list(csv.DictReader((BASE_RESULTS / "task_summary.csv").open(newline="")))
    selected = [row for row in rows if row["assay"] == "countercue_threshold_conditioning"]
    audited = []
    for row in sorted(selected, key=lambda item: (BETA_ARMS.index(item["coupling_variant"]), float(item["length_um"]))):
        path = BASE_RESULTS / "tasks" / row["task_id"] / "trajectory.npz"
        with np.load(path, allow_pickle=False) as raw:
            payload = {key: raw[key] for key in raw.files}
        time = np.asarray(payload["time"], dtype=float)
        rac_time, rac_ref = _settling(time, _half_mass(payload, "rac"))
        mca_time, mca_ref = _settling(time, _half_mass(payload, "mca_bound"))
        accepted = rac_time is not None and mca_time is not None
        audited.append({
            "task_id": row["task_id"], "variant": row["coupling_variant"],
            "length_um": float(row["length_um"]),
            "rac_steady_state_time_s": "" if rac_time is None else rac_time,
            "mca_steady_state_time_s": "" if mca_time is None else mca_time,
            "rac_reference_x50_over_L": rac_ref, "mca_reference_x50_over_L": mca_ref,
            "both_reached_before_10000_s": accepted,
            "opposing_cue_start_s": "" if not accepted else max(rac_time, mca_time),
        })
    return audited


def _response_manifest() -> tuple[dict, list[dict]]:
    audit = _audit()
    base = load_manifest(ROOT / "campaigns" / f"{BASE_NAME}.json")
    templates = {
        (task["coupling_variant"], float(task["simulation"]["length"])): task
        for task in base["tasks"] if task["assay"] == "figure2_B_beta_sweep"
    }
    tasks = []

    def add(row: dict, duration: float, assay: str, save: bool) -> None:
        variant, length = str(row["variant"]), float(row["length_um"])
        start = float(row["opposing_cue_start_s"])
        task = _copy(templates[(variant, length)])
        _set_arm(task, variant)
        task.update(
            assay=assay, input_peak=2.0, countercue_duration_s=duration,
            opposing_cue_duration_s=duration, opposing_cue_start_s=start,
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
                "opposite_start": start, "opposite_stop": start + duration,
                "pre_rac_reference": float(row["rac_reference_x50_over_L"]),
                "pre_mca_reference": float(row["mca_reference_x50_over_L"]),
                "settling_tolerance": SETTLING_TOLERANCE,
                "settling_confirmation_s": CONFIRMATION,
                "flip_tolerance": FLIP_TOLERANCE,
            },
            source_campaign=BASE_NAME, source_assay="countercue_threshold_conditioning",
            save_trajectory=save,
        )
        task["simulation"].update(end_time=start + duration + 10_000.0, output_interval=5.0, step=0.5)
        task["metrics"]["analysis_start"] = start
        tasks.append(task)

    accepted = [row for row in audit if row["both_reached_before_10000_s"]]
    for row in accepted:
        if row["variant"] == "beta_1x":
            for duration in TRAJECTORY_DURATIONS:
                add(row, duration, "countercue_duration_trajectories", True)
    for row in accepted:
        for duration in SWEEP_DURATIONS:
            add(row, duration, "countercue_threshold_sweep", False)

    return make_manifest(
        RESPONSE_NAME,
        tasks,
        description="Figure 2 E-F response assays with Julia-default alpha0, beta0, alpha, and beta.",
        metadata={
            "trajectory_durations_s": list(TRAJECTORY_DURATIONS),
            "sweep_durations_s": list(SWEEP_DURATIONS),
            "accepted_conditioning_cases": len(accepted),
            "changed_parameters_only": {
                "rac_basal_activation": 1.0, "rho_basal_activation": 1.0,
                "rac_mechanical_feedback": 0.2,
                "rho_mechanical_feedback_0x_1x_2x": [0.0, 0.2, 0.4],
            },
        },
    ), audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("base", "response"))
    arguments = parser.parse_args()
    if arguments.stage == "base":
        manifest = _base_manifest()
        destination = ROOT / "campaigns" / f"{BASE_NAME}.json"
    else:
        manifest, audit = _response_manifest()
        destination = ROOT / "campaigns" / f"{RESPONSE_NAME}.json"
        with (BASE_RESULTS / "conditioning_steady_state_audit.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
            writer.writeheader(); writer.writerows(audit)
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
