"""Manifest-driven execution and collection for reproducible WP1 campaigns."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import traceback
from typing import Any, Iterable, Mapping

import numpy as np

from .contract import Grid1D
from .initial import noisy_initial_state
from .metrics import (
    mechanochemical_summary,
    phenotype_summary,
    profile_domain_count,
    tracking_summary,
)
from .protocols import (
    AlternatingGradient,
    GaussianPulse,
    FinitePulseSequence,
    LinearGradient,
    LocalizedPulse,
    NoStimulus,
    ReversingGradient,
    ReversingFiniteGradient,
)
from .registry import make_model
from .simulation import simulate


PROTOCOLS = {
    "none": NoStimulus,
    "localized_pulse": LocalizedPulse,
    "finite_pulse_sequence": FinitePulseSequence,
    "gaussian_pulse": GaussianPulse,
    "linear_gradient": LinearGradient,
    "reversing_gradient": ReversingGradient,
    "reversing_finite_gradient": ReversingFiniteGradient,
    "alternating_gradient": AlternatingGradient,
}

GROUP_KEYS = (
    "tier",
    "model",
    "assay",
    "length_um",
    "coupling_variant",
    "input_peak",
    "tracking_period_s",
)

BOOLEAN_METRICS = (
    "endpoint_state_sustained",
    "final_polarized",
    "ever_polarized",
    "front_establishment_censored",
    "competition_established",
    "competition_resolution_censored",
    "competition_resolution_persistent_at_end",
    "competition_resolution_transient_only",
    "competition_unresolved_persistent_multifront",
    "reversal_conditioned",
    "pre_reversal_polarized",
    "orientation_reversal_conditioned",
    "orientation_reversal_success",
    "orientation_reversal_censored",
    "reversal_success",
    "reversal_censored",
    "tracking_reliable",
)

NUMERIC_METRICS = (
    "final_spatial_contrast",
    "maximum_spatial_contrast",
    "final_domain_count",
    "maximum_domain_count",
    "front_establishment_time_s",
    "front_establishment_time_over_L2",
    "competition_establishment_time_s",
    "competition_resolution_time_s",
    "competition_resolution_time_over_L2",
    "competition_resolution_latency_s",
    "coexistence_duration_s",
    "final_orientation",
    "pre_reversal_time_s",
    "pre_reversal_domain_count",
    "pre_reversal_spatial_contrast",
    "pre_reversal_orientation",
    "orientation_reversal_time_s",
    "orientation_reversal_latency_s",
    "reversal_latency_s",
    "old_front_decay_latency_s",
    "new_front_establishment_latency_s",
    "response_amplitude",
    "phase_lag_seconds",
    "tracking_correlation",
    "rac_proximal_peak_change",
    "rac_distal_peak_change",
    "rac_distal_peak_time_s",
    "rac_distal_to_proximal_abs_ratio",
    "rho_proximal_peak_change",
    "rho_distal_peak_change",
    "rho_distal_peak_time_s",
    "rho_distal_to_proximal_abs_ratio",
    "maximum_absolute_tension",
    "maximum_absolute_velocity",
    "maximum_absolute_displacement",
)


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def task_digest(task: Mapping[str, Any]) -> str:
    """Return a stable digest excluding display/index fields."""
    payload = {
        key: value
        for key, value in task.items()
        if key not in {"task_id", "task_index", "task_digest"}
    }
    encoded = json.dumps(_jsonable(payload), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode()).hexdigest()


def make_manifest(
    name: str,
    tasks: Iterable[Mapping[str, Any]],
    *,
    description: str,
    metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assign deterministic identifiers and return a validated manifest."""
    finalized = []
    for index, raw in enumerate(tasks):
        task = _jsonable(dict(raw))
        digest = task_digest(task)
        task.update(
            {
                "task_index": index,
                "task_digest": digest,
                "task_id": f"{index:05d}-{digest[:12]}",
            }
        )
        finalized.append(task)
    manifest = {
        "schema_version": 1,
        "campaign_name": name,
        "description": description,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "metadata": {} if metadata is None else _jsonable(metadata),
        "task_count": len(finalized),
        "tasks": finalized,
    }
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    if manifest.get("schema_version") != 1:
        raise ValueError("campaign manifest schema_version must be 1")
    tasks = manifest.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ValueError("campaign manifest must contain at least one task")
    if manifest.get("task_count") != len(tasks):
        raise ValueError("task_count does not match tasks")
    identifiers = set()
    for index, task in enumerate(tasks):
        if task.get("task_index") != index:
            raise ValueError(f"task {index} has an invalid task_index")
        expected_digest = task_digest(task)
        if task.get("task_digest") != expected_digest:
            raise ValueError(f"task {index} digest does not match its configuration")
        if task.get("task_id") != f"{index:05d}-{expected_digest[:12]}":
            raise ValueError(f"task {index} has an invalid task_id")
        if task["task_id"] in identifiers:
            raise ValueError(f"duplicate task_id {task['task_id']}")
        identifiers.add(task["task_id"])
        for required in ("model", "assay", "protocol", "simulation", "metrics"):
            if required not in task:
                raise ValueError(f"task {index} is missing {required}")


def load_manifest(path: str | Path) -> dict[str, Any]:
    manifest = json.loads(Path(path).read_text())
    validate_manifest(manifest)
    return manifest


def write_manifest(path: str | Path, manifest: Mapping[str, Any]) -> Path:
    destination = Path(path)
    validate_manifest(manifest)
    _atomic_json(destination, manifest)
    return destination


def _make_protocol(configuration: Mapping[str, Any]):
    parameters = dict(configuration)
    name = parameters.pop("name")
    try:
        constructor = PROTOCOLS[name]
    except KeyError as error:
        raise ValueError(f"unknown campaign protocol {name!r}") from error
    return constructor(**parameters)


def _task_identity(task: Mapping[str, Any], trajectory_metadata=None) -> dict[str, Any]:
    simulation = task["simulation"]
    return {
        "task_id": task["task_id"],
        "task_index": task["task_index"],
        "task_digest": task["task_digest"],
        "tier": task.get("tier", "production"),
        "model": task["model"],
        "assay": task["assay"],
        "length_um": simulation["length"],
        "cells": simulation["cells"],
        "replicate": task.get("replicate", 0),
        "seed": task.get("seed"),
        "coupling_variant": task.get("coupling_variant", "not_applicable"),
        "input_peak": task.get("input_peak"),
        "tracking_period_s": task.get(
            "tracking_period_s", task.get("tracking", {}).get("period")
        ),
        "trajectory_metadata": trajectory_metadata,
    }


def run_manifest_task(
    manifest: Mapping[str, Any],
    task_index: int,
    output_root: str | Path,
    *,
    force: bool = False,
) -> Path:
    """Execute one indexed task with resumable, digest-checked output."""
    validate_manifest(manifest)
    tasks = manifest["tasks"]
    if task_index < 0 or task_index >= len(tasks):
        raise IndexError(f"task index {task_index} outside [0, {len(tasks) - 1}]")
    task = tasks[task_index]
    directory = Path(output_root) / "tasks" / task["task_id"]
    summary_path = directory / "summary.json"
    if summary_path.is_file() and not force:
        existing = json.loads(summary_path.read_text())
        if existing.get("status") == "success" and existing.get("task_digest") == task["task_digest"]:
            return summary_path
    directory.mkdir(parents=True, exist_ok=True)
    _atomic_json(directory / "task.json", task)
    try:
        model = make_model(task["model"], task.get("model_parameters"))
        protocol = _make_protocol(task["protocol"])
        simulation = dict(task["simulation"])
        noise = task.get("initial_noise")
        if noise is not None:
            grid = Grid1D(length=simulation["length"], cells=simulation["cells"])
            simulation["initial_state"] = noisy_initial_state(model, grid, **noise)
        trajectory = simulate(model, protocol=protocol, **simulation)
        summary = {
            "status": "success",
            "campaign_name": manifest["campaign_name"],
            **_task_identity(task, trajectory.metadata),
            **phenotype_summary(trajectory, **task["metrics"]),
        }
        if "tracking" in task:
            summary.update(tracking_summary(trajectory, **task["tracking"]))
        if task.get("mechanochemical_metrics") is not None:
            summary.update(
                mechanochemical_summary(
                    trajectory,
                    **task["mechanochemical_metrics"],
                )
            )
        if task.get("field_audit", False):
            summary.update(_field_audit_summary(trajectory))
        if task.get("single_cue_kinetics") is not None:
            summary.update(
                _single_cue_kinetics_summary(
                    trajectory, **task["single_cue_kinetics"]
                )
            )
        if task.get("pulse_switching") is not None:
            summary.update(
                _pulse_switching_summary(trajectory, **task["pulse_switching"])
            )
        if task.get("strong_competition") is not None:
            summary.update(
                _strong_competition_summary(
                    trajectory, **task["strong_competition"]
                )
            )
        if task.get("strong_reversal") is not None:
            summary.update(
                _strong_reversal_summary(trajectory, **task["strong_reversal"])
            )
        if task.get("polarity_reversal_audit") is not None:
            summary.update(
                _polarity_reversal_audit_summary(
                    trajectory, **task["polarity_reversal_audit"]
                )
            )
        if task.get("pulse_train_frequency") is not None:
            summary.update(
                _pulse_train_frequency_summary(
                    trajectory, **task["pulse_train_frequency"]
                )
            )
        if task.get("save_trajectory", False):
            trajectory.save(directory / "trajectory.npz")
        _atomic_json(summary_path, summary)
        return summary_path
    except Exception as error:
        failure = {
            "status": "failed",
            "campaign_name": manifest["campaign_name"],
            **_task_identity(task),
            "error_type": type(error).__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
        }
        _atomic_json(directory / "failure.json", failure)
        raise


def _signed_orientation(values: np.ndarray, relative: np.ndarray) -> np.ndarray:
    """Return the normalized right-minus-left orientation of every frame."""
    left = np.mean(values[:, relative <= 0.2], axis=1)
    right = np.mean(values[:, relative >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _first_sustained_time(
    condition: np.ndarray,
    times: np.ndarray,
    *,
    after: float,
    duration: float,
) -> float | None:
    start = int(np.searchsorted(times, after, side="left"))
    for index in range(start, len(times)):
        stop = int(np.searchsorted(times, times[index] + duration, side="left"))
        if stop < len(times) and np.all(condition[index : stop + 1]):
            return float(times[index])
    return None


def _single_cue_kinetics_summary(
    trajectory,
    *,
    cue_start: float,
    cue_stop: float,
    cue_width_fraction: float,
    cue_amplitude: float,
    orientation_threshold: float = 0.05,
    polarization_sustain_s: float = 10.0,
    depolarization_sustain_s: float = 50.0,
    cue_side: str = "left",
) -> dict[str, Any]:
    """Measure sustained orientation gain and loss around one withdrawn cue."""
    if cue_side not in {"left", "right"}:
        raise ValueError("cue_side must be left or right")
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    orientation = _signed_orientation(trajectory.observable, relative)
    directed = -orientation if cue_side == "left" else orientation
    polarized = directed >= orientation_threshold
    polarization_time = _first_sustained_time(
        polarized, trajectory.time,
        after=cue_start, duration=polarization_sustain_s,
    )
    depolarization_time = None
    if polarization_time is not None:
        depolarization_time = _first_sustained_time(
            np.abs(orientation) < orientation_threshold,
            trajectory.time,
            after=max(cue_stop, polarization_time),
            duration=depolarization_sustain_s,
        )
    if polarization_time is None:
        status = "not_polarized"
    elif depolarization_time is None:
        status = "depolarization_right_censored"
    else:
        status = "depolarized"
    withdrawal_index = int(np.argmin(np.abs(trajectory.time - cue_stop)))
    peak_index = int(np.argmax(directed))
    return {
        "cue_start_s": cue_start,
        "cue_stop_s": cue_stop,
        "cue_duration_s": cue_stop - cue_start,
        "cue_width_fraction": cue_width_fraction,
        "cue_amplitude": cue_amplitude,
        "cue_side": cue_side,
        "cue_kinetics_status": status,
        "cue_polarized": polarization_time is not None,
        "polarization_time_s": polarization_time,
        "polarization_latency_s": (
            None if polarization_time is None else polarization_time - cue_start
        ),
        "depolarization_time_s": depolarization_time,
        "depolarization_latency_s": (
            None if depolarization_time is None else depolarization_time - cue_stop
        ),
        "polarized_lifetime_s": (
            None
            if polarization_time is None or depolarization_time is None
            else depolarization_time - polarization_time
        ),
        "orientation_at_withdrawal": float(orientation[withdrawal_index]),
        "peak_directed_orientation": float(directed[peak_index]),
        "peak_orientation_time_s": float(trajectory.time[peak_index]),
        "final_cue_orientation": float(orientation[-1]),
        "observation_after_withdrawal_s": float(trajectory.time[-1] - cue_stop),
        "orientation_threshold": orientation_threshold,
        "polarization_sustain_s": polarization_sustain_s,
        "depolarization_sustain_s": depolarization_sustain_s,
    }


def _pulse_switching_summary(
    trajectory,
    *,
    events: list[Mapping[str, Any]],
    orientation_threshold: float = 0.05,
    sustain_s: float = 100.0,
    minimum_persistence_s: float = 500.0,
) -> dict[str, Any]:
    """Score autonomous switching only after each finite pulse is withdrawn."""
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    rac = _signed_orientation(trajectory.fields["rac"], relative)
    rho = _signed_orientation(trajectory.fields["rho"], relative)
    output: dict[str, Any] = {}
    successes = 0
    latencies: list[float] = []
    persistences: list[float] = []
    # Event zero establishes the initial orientation; later events are switch
    # opportunities and determine the reported switching probability.
    for event_index, event in enumerate(events):
        side = str(event["side"])
        start = float(event["start"])
        stop = float(event["stop"])
        if side not in {"left", "right"}:
            raise ValueError("pulse switching event side must be left or right")
        target_rac = -rac if side == "left" else rac
        target_rho = rho if side == "left" else -rho
        target = (target_rac >= orientation_threshold) & (
            target_rho >= orientation_threshold
        )
        next_start = (
            float(events[event_index + 1]["start"])
            if event_index + 1 < len(events)
            else float(trajectory.time[-1])
        )
        establishment = _first_sustained_time(
            target, trajectory.time, after=stop, duration=sustain_s
        )
        if establishment is not None and establishment >= next_start:
            establishment = None
        persistence = 0.0
        if establishment is not None:
            begin = int(np.searchsorted(trajectory.time, establishment, side="left"))
            end = int(np.searchsorted(trajectory.time, next_start, side="left"))
            false_after = np.flatnonzero(~target[begin:end])
            loss_index = begin + int(false_after[0]) if false_after.size else end
            loss_time = (
                float(trajectory.time[loss_index])
                if loss_index < len(trajectory.time)
                else float(trajectory.time[-1])
            )
            persistence = max(0.0, loss_time - establishment)
        success = bool(
            establishment is not None and persistence >= minimum_persistence_s
        )
        prefix = "induction" if event_index == 0 else f"switch_{event_index}"
        output.update(
            {
                f"{prefix}_side": side,
                f"{prefix}_success": success,
                f"{prefix}_establishment_time_s": establishment,
                f"{prefix}_latency_from_pulse_start_s": (
                    None if establishment is None else establishment - start
                ),
                f"{prefix}_latency_after_withdrawal_s": (
                    None if establishment is None else establishment - stop
                ),
                f"{prefix}_persistence_s": persistence,
            }
        )
        if event_index > 0:
            successes += int(success)
            if establishment is not None:
                latencies.append(establishment - start)
                persistences.append(persistence)
    opportunities = max(0, len(events) - 1)
    output.update(
        {
            "switch_opportunities": opportunities,
            "successful_switches": successes,
            "switch_probability": (
                successes / opportunities if opportunities else None
            ),
            "mean_switch_latency_s": (
                float(np.mean(latencies)) if latencies else None
            ),
            "minimum_post_pulse_persistence_s": (
                float(np.min(persistences)) if persistences else None
            ),
            "pulse_switch_orientation_threshold": orientation_threshold,
            "pulse_switch_sustain_s": sustain_s,
            "pulse_switch_minimum_persistence_s": minimum_persistence_s,
        }
    )
    return output


def _first_sustained_before(
    condition: np.ndarray,
    times: np.ndarray,
    *,
    after: float,
    before: float,
    duration: float,
) -> float | None:
    """Find a sustained interval wholly contained in ``[after, before)``."""
    start = int(np.searchsorted(times, after, side="left"))
    limit = int(np.searchsorted(times, before, side="left"))
    for index in range(start, limit):
        stop = int(np.searchsorted(times, times[index] + duration, side="left"))
        if stop < limit and np.all(condition[index : stop + 1]):
            return float(times[index])
    return None


def _contrast_series(values: np.ndarray) -> np.ndarray:
    high = np.max(values, axis=1)
    low = np.min(values, axis=1)
    return (high - low) / np.maximum(np.abs(high) + np.abs(low), 1e-12)


def _strong_competition_summary(
    trajectory,
    *,
    cue_stop: float,
    two_pole_deadline: float,
    two_pole_strength_threshold: float = 0.05,
    two_pole_sustain_s: float = 20.0,
    contrast_threshold: float = 0.3,
    orientation_threshold: float = 0.25,
    resolution_sustain_s: float = 100.0,
    endpoint_sustain_s: float = 500.0,
) -> dict[str, Any]:
    """Require real competing Rac peaks followed by a strong persistent winner."""
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    rac = np.asarray(trajectory.fields["rac"], dtype=float)
    left = np.mean(rac[:, relative <= 0.15], axis=1)
    right = np.mean(rac[:, relative >= 0.85], axis=1)
    centre = np.mean(rac[:, (relative >= 0.4) & (relative <= 0.6)], axis=1)
    two_strength = np.minimum(
        (left - centre) / np.maximum(np.abs(left) + np.abs(centre), 1e-12),
        (right - centre) / np.maximum(np.abs(right) + np.abs(centre), 1e-12),
    )
    orientation = _signed_orientation(rac, relative)
    contrast = _contrast_series(rac)
    domains = np.asarray([profile_domain_count(profile) for profile in rac], dtype=int)
    two_state = (domains >= 2) & (two_strength >= two_pole_strength_threshold)
    two_time = _first_sustained_before(
        two_state, trajectory.time,
        after=0.0, before=two_pole_deadline, duration=two_pole_sustain_s,
    )
    strong_single = (
        (domains == 1)
        & (contrast >= contrast_threshold)
        & (np.abs(orientation) >= orientation_threshold)
    )
    resolution = None
    if two_time is not None:
        resolution = _first_sustained_time(
            strong_single, trajectory.time,
            after=max(cue_stop, two_time), duration=resolution_sustain_s,
        )
    tail_start = max(float(trajectory.time[0]), float(trajectory.time[-1]) - endpoint_sustain_s)
    tail = trajectory.time >= tail_start
    persistent = bool(resolution is not None and np.all(strong_single[tail]))
    success = bool(two_time is not None and persistent)
    if two_time is None:
        status = "two_poles_not_established"
    elif resolution is None:
        status = "two_poles_unresolved"
    elif not persistent:
        status = "resolution_not_persistent"
    else:
        status = "resolved_persistent"
    final = rac[-1]
    low = float(np.min(final))
    amplitude = float(np.max(final) - low)
    return {
        "strong_competition_success": success,
        "strong_competition_status": status,
        "strong_two_pole_establishment_time_s": two_time,
        "strong_two_pole_peak_strength": float(np.max(two_strength[trajectory.time < two_pole_deadline])),
        "strong_resolution_time_s": resolution,
        "strong_resolution_latency_after_cue_s": (
            None if resolution is None else resolution - cue_stop
        ),
        "strong_resolution_persistent_at_end": persistent,
        "strong_final_rac_contrast": float(contrast[-1]),
        "strong_final_rac_orientation": float(orientation[-1]),
        "strong_final_rac_domain_count": int(domains[-1]),
        "strong_final_front_half_height_width_L": float(np.trapezoid(
            (final >= low + 0.5 * amplitude).astype(float), relative
        )),
        "strong_final_front_80pct_width_L": float(np.trapezoid(
            (final >= low + 0.8 * amplitude).astype(float), relative
        )),
    }


def _strong_reversal_summary(
    trajectory,
    *,
    opposite_start: float,
    opposite_stop: float,
    response_deadline: float,
    contrast_threshold: float = 0.3,
    orientation_threshold: float = 0.25,
    sustain_s: float = 100.0,
    minimum_persistence_s: float = 1_000.0,
) -> dict[str, Any]:
    """Score a strong autonomous reversal within a bounded causal window."""
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    rac = _signed_orientation(trajectory.fields["rac"], relative)
    rho = _signed_orientation(trajectory.fields["rho"], relative)
    contrast = _contrast_series(trajectory.fields["rac"])
    pre_index = max(0, int(np.searchsorted(trajectory.time, opposite_start, side="left")) - 1)
    preconditioned = bool(
        -rac[pre_index] >= orientation_threshold
        and rho[pre_index] >= orientation_threshold
        and contrast[pre_index] >= contrast_threshold
    )
    target = (
        (rac >= orientation_threshold)
        & (rho <= -orientation_threshold)
        & (contrast >= contrast_threshold)
    )
    establishment = _first_sustained_before(
        target, trajectory.time,
        after=opposite_stop, before=response_deadline, duration=sustain_s,
    )
    persistence = 0.0
    if establishment is not None:
        begin = int(np.searchsorted(trajectory.time, establishment, side="left"))
        false_after = np.flatnonzero(~target[begin:])
        loss_index = begin + int(false_after[0]) if false_after.size else len(trajectory.time) - 1
        persistence = max(0.0, float(trajectory.time[loss_index] - establishment))
    success = bool(
        preconditioned
        and establishment is not None
        and persistence >= minimum_persistence_s
    )
    if not preconditioned:
        status = "not_preconditioned"
    elif establishment is None:
        status = "no_causal_reversal"
    elif persistence < minimum_persistence_s:
        status = "reversal_not_persistent"
    else:
        status = "reversed_persistent"
    return {
        "strong_reversal_success": success,
        "strong_reversal_status": status,
        "strong_reversal_preconditioned": preconditioned,
        "strong_pre_opposite_rac_directed": float(-rac[pre_index]),
        "strong_pre_opposite_rho_directed": float(rho[pre_index]),
        "strong_pre_opposite_rac_contrast": float(contrast[pre_index]),
        "strong_reversal_establishment_time_s": establishment,
        "strong_reversal_latency_from_pulse_start_s": (
            None if establishment is None else establishment - opposite_start
        ),
        "strong_reversal_latency_after_withdrawal_s": (
            None if establishment is None else establishment - opposite_stop
        ),
        "strong_reversal_persistence_s": persistence,
        "strong_reversal_final_rac_directed": float(rac[-1]),
        "strong_reversal_final_rho_directed": float(-rho[-1]),
        "strong_reversal_final_rac_contrast": float(contrast[-1]),
    }


def _half_mass_position_series(trajectory, field: str) -> np.ndarray:
    """Return the interpolated relative half-mass position for a trajectory field."""
    values = np.asarray(trajectory.fields[field], dtype=float)
    cumulative = np.cumsum(values, axis=1)
    half = 0.5 * cumulative[:, -1]
    crossing = np.argmax(cumulative >= half[:, None], axis=1)
    rows = np.arange(len(values))
    before = np.where(
        crossing > 0,
        cumulative[rows, np.maximum(crossing - 1, 0)],
        0.0,
    )
    fraction = np.divide(
        half - before,
        values[rows, crossing],
        out=np.full_like(half, 0.5),
        where=values[rows, crossing] > 0.0,
    )
    return np.clip((crossing + np.clip(fraction, 0.0, 1.0)) / values.shape[1], 0.0, 1.0)


def _polarity_reversal_audit_summary(
    trajectory,
    *,
    opposite_start: float,
    opposite_stop: float,
    pre_rac_reference: float,
    pre_mca_reference: float,
    settling_tolerance: float = 0.002,
    settling_confirmation_s: float = 1_000.0,
    flip_tolerance: float = 0.05,
) -> dict[str, Any]:
    """Audit steady repolarization toward the sign-flipped pre-cue half-mass state."""
    time = np.asarray(trajectory.time, dtype=float)
    pre_index = max(0, int(np.searchsorted(time, opposite_start, side="left")) - 1)
    result: dict[str, Any] = {}
    successes = []
    for field, label, reference in (
        ("rac", "rac", float(pre_rac_reference)),
        ("mca_bound", "mca", float(pre_mca_reference)),
    ):
        score = _half_mass_position_series(trajectory, field)
        tail = time >= time[-1] - settling_confirmation_s
        final_reference = float(np.median(score[tail]))
        tail_stable = bool(
            np.max(np.abs(score[tail] - final_reference)) <= settling_tolerance
        )
        within = np.abs(score - final_reference) <= settling_tolerance
        remains = np.logical_and.accumulate(within[::-1])[::-1]
        eligible = (
            remains
            & (time >= opposite_stop)
            & (time <= time[-1] - settling_confirmation_s)
        )
        indices = np.flatnonzero(eligible)
        settling_time = None if not tail_stable or len(indices) == 0 else float(time[indices[0]])
        target = 1.0 - reference
        sign_flipped = bool((reference - 0.5) * (final_reference - 0.5) < 0.0)
        target_matched = bool(abs(final_reference - target) <= flip_tolerance)
        success = bool(settling_time is not None and sign_flipped and target_matched)
        successes.append(success)
        result.update({
            f"countercue_{label}_pre_x50_over_L": float(score[pre_index]),
            f"countercue_{label}_conditioning_reference_x50_over_L": reference,
            f"countercue_{label}_final_x50_over_L": final_reference,
            f"countercue_{label}_target_x50_over_L": target,
            f"countercue_{label}_final_steady": settling_time is not None,
            f"countercue_{label}_settling_time_s": settling_time,
            f"countercue_{label}_sign_flipped": sign_flipped,
            f"countercue_{label}_target_matched": target_matched,
            f"countercue_{label}_reversal_success": success,
        })
    result.update({
        "countercue_either_reversal_success": bool(any(successes)),
        "countercue_both_reversal_success": bool(all(successes)),
        "countercue_opposite_start_s": opposite_start,
        "countercue_opposite_stop_s": opposite_stop,
        "countercue_duration_s": opposite_stop - opposite_start,
        "countercue_settling_tolerance_fraction_L": settling_tolerance,
        "countercue_settling_confirmation_s": settling_confirmation_s,
        "countercue_flip_tolerance_fraction_L": flip_tolerance,
    })
    return result


def _pulse_train_frequency_summary(
    trajectory,
    *,
    start: float,
    period: float,
    half_cycles: int,
    contrast_threshold: float = 0.3,
    orientation_threshold: float = 0.25,
) -> dict[str, Any]:
    """Measure deterministic gain, phase, and event tracking of a square train."""
    if period <= 0 or half_cycles < 2:
        raise ValueError("pulse-train period and half-cycle count are invalid")
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    rac = _signed_orientation(trajectory.fields["rac"], relative)
    rho = _signed_orientation(trajectory.fields["rho"], relative)
    contrast = _contrast_series(trajectory.fields["rac"])
    half_period = period / 2.0
    tracked = 0
    peak_directed = []
    for event in range(half_cycles):
        event_start = start + event * half_period
        event_stop = event_start + half_period
        sign = 1.0 if event % 2 == 0 else -1.0
        target = (
            (sign * rac >= orientation_threshold)
            & (-sign * rho >= orientation_threshold)
            & (contrast >= contrast_threshold)
        )
        establishment = _first_sustained_before(
            target, trajectory.time,
            after=event_start, before=event_stop,
            duration=min(20.0, half_period / 4.0),
        )
        tracked += int(establishment is not None)
        window = (trajectory.time >= event_start) & (trajectory.time < event_stop)
        peak_directed.append(float(np.max(sign * rac[window])))

    fit_start = start + period
    fit_stop = start + half_cycles * half_period
    fit = (trajectory.time >= fit_start) & (trajectory.time < fit_stop)
    phase = 2.0 * np.pi * (trajectory.time[fit] - start) / period
    response = rac[fit]
    design = np.column_stack((np.ones(np.count_nonzero(fit)), np.sin(phase), np.cos(phase)))
    coefficients, *_ = np.linalg.lstsq(design, response, rcond=None)
    fitted = design @ coefficients
    residual = float(np.sum((response - fitted) ** 2))
    total = float(np.sum((response - np.mean(response)) ** 2))
    amplitude = float(np.hypot(coefficients[1], coefficients[2]))
    lag = float(np.arctan2(-coefficients[2], coefficients[1]))
    lag = ((lag + np.pi) % (2.0 * np.pi)) - np.pi
    lag_fraction = lag / (2.0 * np.pi)
    fraction = tracked / half_cycles
    reliable = bool(fraction >= 0.8 and amplitude >= 0.1)
    return {
        "pulse_train_tracked_events": tracked,
        "pulse_train_event_count": half_cycles,
        "pulse_train_tracking_fraction": fraction,
        "pulse_train_mean_peak_directed_rac": float(np.mean(peak_directed)),
        "pulse_train_minimum_peak_directed_rac": float(np.min(peak_directed)),
        "pulse_train_response_amplitude": amplitude,
        "pulse_train_phase_lag_fraction": lag_fraction,
        "pulse_train_harmonic_r2": 1.0 if total <= 1e-12 else 1.0 - residual / total,
        "pulse_train_reliable": reliable,
    }


def _field_audit_summary(trajectory) -> dict[str, Any]:
    """Numerical and multi-field diagnostics for mechanochemical campaigns."""
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    result: dict[str, Any] = {}
    for field in ("rac", "rho", "tension", "mca_bound", "mca_unbound"):
        if field not in trajectory.fields:
            continue
        values = np.asarray(trajectory.fields[field], dtype=float)
        orientation = _signed_orientation(values, relative)
        left = np.mean(values[:, relative <= 0.2], axis=1)
        right = np.mean(values[:, relative >= 0.8], axis=1)
        result.update(
            {
                f"{field}_minimum": float(np.min(values)),
                f"{field}_final_orientation": float(orientation[-1]),
                f"{field}_maximum_absolute_orientation": float(
                    np.max(np.abs(orientation))
                ),
                f"{field}_front_initial": float(left[0]),
                f"{field}_front_final": float(left[-1]),
                f"{field}_back_initial": float(right[0]),
                f"{field}_back_final": float(right[-1]),
            }
        )
    if {"mca_bound", "mca_unbound"} <= trajectory.fields.keys():
        total = np.sum(
            trajectory.fields["mca_bound"] + trajectory.fields["mca_unbound"],
            axis=1,
        )
        result["total_mca_initial"] = float(total[0])
        result["total_mca_final"] = float(total[-1])
        result["total_mca_max_relative_error"] = float(
            np.max(np.abs(total / total[0] - 1.0))
        )
    concentration_fields = ("rac", "rho", "mca_bound", "mca_unbound")
    finite = all(
        np.all(np.isfinite(trajectory.fields[field])) for field in concentration_fields
    )
    nonnegative = all(
        np.min(trajectory.fields[field]) >= -1e-10 for field in concentration_fields
    )
    conserved = result.get("total_mca_max_relative_error", np.inf) <= 1e-10
    result.update(
        {
            "field_audit_finite": bool(finite),
            "field_audit_nonnegative": bool(nonnegative),
            "field_audit_mca_conserved": bool(conserved),
            "field_audit_passed": bool(finite and nonnegative and conserved),
        }
    )
    return result


def _csv_value(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"))
    return value


def _write_csv(path: Path, rows: list[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key)) for key in fields})


def _group_value(row: Mapping[str, Any], key: str) -> Any:
    value = row.get(key)
    return "not_applicable" if value is None else value


def _aggregate_rows(rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[Mapping[str, Any]]] = {}
    for row in rows:
        key = tuple(_group_value(row, field) for field in GROUP_KEYS)
        groups.setdefault(key, []).append(row)
    aggregates = []
    for key, members in sorted(groups.items(), key=lambda item: tuple(map(str, item[0]))):
        successful = [row for row in members if row.get("status") == "success"]
        aggregate = dict(zip(GROUP_KEYS, key))
        aggregate.update(
            {
                "n_expected": len(members),
                "n_success": len(successful),
                "n_failed": len(members) - len(successful),
            }
        )
        for state in ("unpolarized", "single_front", "multifront", "polarized_no_resolved_peak"):
            count = sum(row.get("endpoint_state") == state for row in successful)
            aggregate[f"endpoint_{state}_count"] = count
            aggregate[f"endpoint_{state}_fraction"] = (
                None if not successful else count / len(successful)
            )
        for metric in BOOLEAN_METRICS:
            values = [row[metric] for row in successful if isinstance(row.get(metric), bool)]
            aggregate[f"{metric}_fraction"] = (
                None if not values else float(np.mean(values))
            )
        for metric in NUMERIC_METRICS:
            values = [
                float(row[metric])
                for row in successful
                if isinstance(row.get(metric), (int, float))
                and not isinstance(row.get(metric), bool)
                and np.isfinite(row[metric])
            ]
            aggregate[f"{metric}_median"] = (
                None if not values else float(np.median(values))
            )
        aggregates.append(aggregate)
    return aggregates


def collect_campaign(
    manifest: Mapping[str, Any],
    output_root: str | Path,
) -> dict[str, Any]:
    """Collect successful, failed, and missing tasks plus group summaries."""
    validate_manifest(manifest)
    root = Path(output_root)
    rows: list[dict[str, Any]] = []
    missing = []
    invalid = []
    for task in manifest["tasks"]:
        directory = root / "tasks" / task["task_id"]
        summary_path = directory / "summary.json"
        failure_path = directory / "failure.json"
        if summary_path.is_file():
            row = json.loads(summary_path.read_text())
        elif failure_path.is_file():
            row = json.loads(failure_path.read_text())
        else:
            missing.append(task["task_id"])
            continue
        if row.get("task_digest") != task["task_digest"]:
            invalid.append(task["task_id"])
            continue
        rows.append(row)
    aggregates = _aggregate_rows(rows)
    _write_csv(root / "task_summary.csv", rows)
    _write_csv(root / "group_summary.csv", aggregates)
    report = {
        "campaign_name": manifest["campaign_name"],
        "task_count": manifest["task_count"],
        "collected_count": len(rows),
        "success_count": sum(row.get("status") == "success" for row in rows),
        "failure_count": sum(row.get("status") == "failed" for row in rows),
        "missing_count": len(missing),
        "invalid_count": len(invalid),
        "missing_task_ids": missing,
        "invalid_task_ids": invalid,
        "complete": not missing and not invalid and all(
            row.get("status") == "success" for row in rows
        ),
        "group_count": len(aggregates),
    }
    _atomic_json(root / "collection_summary.json", report)
    return report
