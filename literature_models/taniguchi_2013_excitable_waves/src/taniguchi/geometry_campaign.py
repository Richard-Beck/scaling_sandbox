"""Streaming simulations and measurements for the Taniguchi geometry campaign.

The campaign uses the Fig. 4F phase-field chemistry in every arm.  A
``fixed_spontaneous`` or ``fixed_controlled`` run holds the initialized diffuse
circle fixed and advances only the conservative chemical fields; a
``deformable_spontaneous`` run advances the identical chemistry together with
the published boundary equation.

Only compact scalar/profile trajectories, wave-track summaries, stochastic
event coordinates, and eight activity-aware snapshots are retained.  Full
movies can therefore be recreated for selected deterministic seeds without
making every production replicate a multi-gigabyte file.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage, signal

from .phase_field import (
    PHASE_PARAMETER_SETS,
    _chemical_step_fast,
    _initial_phi,
    _phase_step_fast,
)


TRAJECTORY_KEYS = (
    "area", "active_area_fraction", "active_patch_count", "maximum_v",
    "centroid_x", "centroid_y", "aspect_ratio", "shape_anisotropy",
    "perimeter", "circularity", "boundary_radial_roughness",
    "protrusion_count", "dominant_patch_fraction", "activity_moment",
    "activity_angle", "boundary_active_fraction", "phase_singularity_count",
    "visited_area_fraction", "reexcited_area_fraction",
)


@dataclass
class ComponentState:
    track_id: int
    start_time: float
    last_time: float
    origin_x: float
    origin_y: float
    times: list[float]
    areas: list[float]
    centroid_x: list[float]
    centroid_y: list[float]
    lead_distance: list[float]
    boundary_contact: list[bool]
    reflected: bool = False
    ended_by_merge: bool = False
    right_censored: bool = False


def geometry_scales(volume_ratio: float) -> tuple[float, float]:
    linear = float(volume_ratio) ** (1.0 / 3.0)
    return linear, linear * linear


def _filled_radial_profile(
    boundary: np.ndarray,
    values: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    cx: float,
    cy: float,
    bins: int,
    *,
    maximum: bool,
) -> np.ndarray:
    yy, xx = np.nonzero(boundary)
    if yy.size == 0:
        return np.zeros(bins, dtype=float)
    angles = np.mod(np.arctan2(Y[yy, xx] - cy, X[yy, xx] - cx), 2.0 * np.pi)
    indices = np.minimum((angles * bins / (2.0 * np.pi)).astype(int), bins - 1)
    output = np.full(bins, np.nan)
    for index in range(bins):
        selected = indices == index
        if np.any(selected):
            selected_values = values[yy[selected], xx[selected]]
            output[index] = (
                float(np.max(selected_values)) if maximum
                else float(np.mean(selected_values))
            )
    good = np.flatnonzero(np.isfinite(output))
    if good.size == 0:
        return np.zeros(bins, dtype=float)
    if good.size == 1:
        output[:] = output[good[0]]
        return output
    extended_x = np.concatenate((good - bins, good, good + bins))
    extended_y = np.tile(output[good], 3)
    missing = np.flatnonzero(~np.isfinite(output))
    output[missing] = np.interp(missing, extended_x, extended_y)
    return output


def _phase_singularities(
    phi: np.ndarray,
    U: np.ndarray,
    V: np.ndarray,
    baseline_u: float,
    scale: float,
) -> int:
    """Count interior state-space phase-winding candidates.

    The phase is the polar angle about the homogeneous state in the local
    ``(U,V)`` plane.  A candidate is a plaquette with winding magnitude above
    pi whose corners include a chemically displaced state.  Nearby candidate
    plaquettes are counted once.  This is an operational spiral-core measure,
    not a claim that every one-frame winding is a stable spiral.
    """
    phase = np.arctan2(V / scale, (U - baseline_u) / scale)
    amplitude = np.hypot(V / scale, (U - baseline_u) / scale)

    def wrapped(delta: np.ndarray) -> np.ndarray:
        return (delta + np.pi) % (2.0 * np.pi) - np.pi

    p00 = phase[:-1, :-1]
    p01 = phase[:-1, 1:]
    p11 = phase[1:, 1:]
    p10 = phase[1:, :-1]
    winding = (wrapped(p01 - p00) + wrapped(p11 - p01)
               + wrapped(p10 - p11) + wrapped(p00 - p10))
    interior = ndimage.binary_erosion(phi >= 0.5, iterations=3)
    valid = (interior[:-1, :-1] & interior[:-1, 1:]
             & interior[1:, :-1] & interior[1:, 1:])
    displaced = np.maximum.reduce((
        amplitude[:-1, :-1], amplitude[:-1, 1:],
        amplitude[1:, :-1], amplitude[1:, 1:],
    )) >= 0.08
    candidates = valid & displaced & (np.abs(winding) > np.pi)
    return int(ndimage.label(candidates, structure=np.ones((3, 3), dtype=int))[1])


def _frame_geometry(
    phi: np.ndarray,
    V: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    dx: float,
    radial_bins: int,
) -> tuple[dict[str, float], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pixel_area = dx * dx
    weight = np.maximum(phi, 0.0)
    total = float(np.sum(weight))
    cx = float(np.sum(weight * X) / max(total, 1e-12))
    cy = float(np.sum(weight * Y) / max(total, 1e-12))
    xx = float(np.sum(weight * (X - cx) ** 2) / max(total, 1e-12))
    yy = float(np.sum(weight * (Y - cy) ** 2) / max(total, 1e-12))
    xy = float(np.sum(weight * (X - cx) * (Y - cy)) / max(total, 1e-12))
    eigenvalues = np.maximum(
        np.linalg.eigvalsh(np.asarray([[xx, xy], [xy, yy]])), 1e-12
    )
    anisotropy = float(
        (eigenvalues[-1] - eigenvalues[0]) /
        (eigenvalues[-1] + eigenvalues[0])
    )
    aspect_ratio = float(np.sqrt(eigenvalues[-1] / eigenvalues[0]))
    gradient_y, gradient_x = np.gradient(phi, dx, dx)
    perimeter = float(np.sum(np.hypot(gradient_x, gradient_y)) * pixel_area)
    area = total * pixel_area
    support = phi >= 0.5
    sharp_area = float(np.count_nonzero(support) * pixel_area)
    circularity = float(
        min(1.0, 4.0 * np.pi * sharp_area / max(perimeter * perimeter, 1e-12))
    )
    boundary = support & ~ndimage.binary_erosion(
        support, structure=np.ones((3, 3), dtype=bool), border_value=0
    )
    radius_values = np.hypot(X - cx, Y - cy)
    radial_profile = _filled_radial_profile(
        boundary, radius_values, X, Y, cx, cy, radial_bins, maximum=True
    )
    pip3_profile = _filled_radial_profile(
        boundary, V, X, Y, cx, cy, radial_bins, maximum=False
    )
    mean_radius = float(np.mean(radial_profile))
    roughness = float(np.std(radial_profile) / max(mean_radius, 1e-12))
    peaks, _ = signal.find_peaks(
        np.r_[radial_profile, radial_profile, radial_profile],
        prominence=max(0.02 * mean_radius, 2.0 * dx),
        distance=max(1, radial_bins // 12),
    )
    protrusion_count = int(np.count_nonzero((peaks >= radial_bins)
                                            & (peaks < 2 * radial_bins)))
    return ({
        "area": area,
        "centroid_x": cx,
        "centroid_y": cy,
        "aspect_ratio": aspect_ratio,
        "shape_anisotropy": anisotropy,
        "perimeter": perimeter,
        "circularity": circularity,
        "boundary_radial_roughness": roughness,
        "protrusion_count": float(protrusion_count),
    }, support, boundary, radial_profile, pip3_profile)


class WaveTracker:
    def __init__(self, shape: tuple[int, int], dx: float, minimum_area: float):
        self.dx = dx
        self.minimum_pixels = max(1, int(math.ceil(minimum_area / (dx * dx))))
        self.previous_tracks = np.zeros(shape, dtype=np.int32)
        self.next_track_id = 1
        self.tracks: dict[int, ComponentState] = {}
        self.finished: list[ComponentState] = []
        self.merges = 0
        self.splits = 0
        self.births = 0

    def update(
        self,
        time: float,
        active: np.ndarray,
        support: np.ndarray,
        X: np.ndarray,
        Y: np.ndarray,
    ) -> tuple[int, float, float, np.ndarray]:
        labels, count = ndimage.label(active, structure=np.ones((3, 3), dtype=int))
        if count:
            sizes = np.bincount(labels.ravel())
            keep = np.flatnonzero(sizes >= self.minimum_pixels)
            keep = keep[keep != 0]
        else:
            keep = np.empty(0, dtype=int)
        current_masks = [labels == label for label in keep]
        previous_to_current: dict[int, list[int]] = {}
        overlaps: list[list[tuple[int, int]]] = []
        for current_index, mask in enumerate(current_masks):
            ids, counts = np.unique(self.previous_tracks[mask], return_counts=True)
            candidates = [(int(track_id), int(n)) for track_id, n in zip(ids, counts)
                          if track_id > 0]
            candidates.sort(key=lambda item: item[1], reverse=True)
            overlaps.append(candidates)
            for track_id, _ in candidates:
                previous_to_current.setdefault(track_id, []).append(current_index)
        self.splits += sum(max(0, len(indices) - 1)
                           for indices in previous_to_current.values())

        assignments: dict[int, int] = {}
        used: set[int] = set()
        for current_index in sorted(range(len(current_masks)),
                                    key=lambda i: sum(current_masks[i].ravel()),
                                    reverse=True):
            available = [item for item in overlaps[current_index] if item[0] not in used]
            if available:
                track_id = available[0][0]
                assignments[current_index] = track_id
                used.add(track_id)
                if len(overlaps[current_index]) > 1:
                    merged = [item[0] for item in overlaps[current_index][1:]]
                    self.merges += len(merged)
                    for old_id in merged:
                        old = self.tracks.pop(old_id, None)
                        if old is not None:
                            old.ended_by_merge = True
                            self.finished.append(old)
            else:
                track_id = self.next_track_id
                self.next_track_id += 1
                assignments[current_index] = track_id
                self.births += 1
                yy, xx = np.nonzero(current_masks[current_index])
                origin_x = float(np.mean(X[yy, xx]))
                origin_y = float(np.mean(Y[yy, xx]))
                self.tracks[track_id] = ComponentState(
                    track_id, time, time, origin_x, origin_y,
                    [], [], [], [], [], [],
                )

        current_map = np.zeros_like(self.previous_tracks)
        boundary_distance = ndimage.distance_transform_edt(support) * self.dx
        boundary_zone = max(0.5, 1.5 * self.dx)
        areas: list[float] = []
        for current_index, mask in enumerate(current_masks):
            track_id = assignments[current_index]
            current_map[mask] = track_id
            yy, xx = np.nonzero(mask)
            area = float(yy.size * self.dx * self.dx)
            cx = float(np.mean(X[yy, xx]))
            cy = float(np.mean(Y[yy, xx]))
            track = self.tracks[track_id]
            lead = float(np.max(np.hypot(X[yy, xx] - track.origin_x,
                                        Y[yy, xx] - track.origin_y)))
            contact = bool(np.min(boundary_distance[yy, xx]) <= boundary_zone)
            if any(track.boundary_contact) and not contact and time - track.start_time >= 0.5:
                track.reflected = True
            track.last_time = time
            track.times.append(time)
            track.areas.append(area)
            track.centroid_x.append(cx)
            track.centroid_y.append(cy)
            track.lead_distance.append(lead)
            track.boundary_contact.append(contact)
            areas.append(area)

        present = set(assignments.values())
        for track_id in list(self.tracks):
            if track_id not in present:
                self.finished.append(self.tracks.pop(track_id))
        self.previous_tracks = current_map
        total_active = float(sum(areas))
        dominant_fraction = max(areas, default=0.0) / max(total_active, 1e-12)
        filtered_active = current_map > 0
        boundary_active = float(np.count_nonzero(
            filtered_active & (boundary_distance <= boundary_zone)
        ))
        boundary_fraction = boundary_active / max(
            float(np.count_nonzero(filtered_active)), 1.0
        )
        return len(current_masks), dominant_fraction, boundary_fraction, current_map > 0

    def finish(self) -> list[ComponentState]:
        for track in self.tracks.values():
            track.right_censored = True
            self.finished.append(track)
        self.tracks.clear()
        return self.finished


def _track_rows(
    tracks: list[ComponentState], diameter: float, output_interval: float
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for track in tracks:
        times = np.asarray(track.times, dtype=float)
        leads = np.asarray(track.lead_distance, dtype=float)
        if times.size >= 2:
            speed_steps = np.diff(leads) / np.diff(times)
            advancing = speed_steps[speed_steps > 0]
            speed = float(np.median(advancing)) if advancing.size else 0.0
        else:
            speed = 0.0
        lifetime = max(output_interval, track.last_time - track.start_time + output_interval)
        contacted = bool(any(track.boundary_contact))
        rows.append({
            "track_id": track.track_id,
            "start_time": track.start_time,
            "end_time": track.last_time,
            "lifetime": lifetime,
            "maximum_area": max(track.areas, default=0.0),
            "maximum_lead_distance": max(track.lead_distance, default=0.0),
            "travel_per_diameter": max(track.lead_distance, default=0.0) / diameter,
            "front_speed": speed,
            "boundary_contact": contacted,
            "reflection_candidate": bool(track.reflected),
            "ended_by_merge": bool(track.ended_by_merge),
            "right_censored": bool(track.right_censored),
            "extinction_candidate": bool(not track.ended_by_merge and not track.right_censored),
        })
    return rows


def _migration_summary(
    times: np.ndarray, x: np.ndarray, y: np.ndarray, radius: float
) -> dict[str, float | None]:
    coordinates = np.column_stack((x - x[0], y - y[0]))
    steps = np.diff(coordinates, axis=0)
    step_lengths = np.linalg.norm(steps, axis=1)
    duration = float(times[-1] - times[0])
    net = float(np.linalg.norm(coordinates[-1]))
    path = float(np.sum(step_lengths))
    lag = max(1, int(round(0.5 / np.median(np.diff(times)))))
    velocity = coordinates[lag:] - coordinates[:-lag]
    norms = np.linalg.norm(velocity, axis=1)
    good = norms > 1e-5 * radius
    turning: list[float] = []
    for first, second, good_first, good_second in zip(
        velocity[:-1], velocity[1:], good[:-1], good[1:]
    ):
        if good_first and good_second:
            cosine = float(np.dot(first, second) /
                           (np.linalg.norm(first) * np.linalg.norm(second)))
            turning.append(float(np.arccos(np.clip(cosine, -1.0, 1.0))))
    persistence_time: float | None = None
    if np.count_nonzero(good) >= 4:
        unit = velocity[good] / norms[good, None]
        maximum_lag = min(200, unit.shape[0] // 2)
        for index in range(1, maximum_lag + 1):
            correlation = float(np.mean(np.sum(unit[:-index] * unit[index:], axis=1)))
            if correlation <= math.exp(-1.0):
                persistence_time = index * float(np.median(np.diff(times)))
                break
    mean_speed = path / duration if duration > 0 else 0.0
    return {
        "net_displacement": net,
        "path_length": path,
        "mean_speed": mean_speed,
        "straightness": net / path if path > 0 else 0.0,
        "net_displacement_per_radius": net / radius,
        "path_length_per_radius": path / radius,
        "mean_speed_per_radius": mean_speed / radius,
        "mean_absolute_turning_angle": float(np.mean(turning)) if turning else None,
        "persistence_time": persistence_time,
        "persistence_distance": (
            mean_speed * persistence_time
            if persistence_time is not None else None
        ),
        "persistence_distance_per_radius": (
            mean_speed * persistence_time / radius
            if persistence_time is not None else None
        ),
    }


def _lagged_boundary_coupling(
    times: np.ndarray, radial: np.ndarray, activity: np.ndarray, lag_time: float = 0.5
) -> float | None:
    lag = max(1, int(round(lag_time / np.median(np.diff(times)))))
    if radial.shape[0] <= lag:
        return None
    future_extension = radial[lag:] - radial[:-lag]
    current_activity = activity[:-lag]
    current_activity = current_activity - np.mean(current_activity, axis=1, keepdims=True)
    future_extension = future_extension - np.mean(future_extension, axis=1, keepdims=True)
    first = current_activity.ravel()
    second = future_extension.ravel()
    if np.std(first) <= 1e-12 or np.std(second) <= 1e-12:
        return None
    return float(np.corrcoef(first, second)[0, 1])


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, np.ndarray):
        return [_json_ready(item) for item in value.tolist()]
    return value


def run_geometry_campaign(config: dict[str, Any], output: Path) -> dict[str, Any]:
    arm = str(config["arm"])
    if arm not in {"fixed_spontaneous", "fixed_controlled", "deformable_spontaneous"}:
        raise ValueError(f"unknown Taniguchi campaign arm {arm!r}")
    volume_ratio = float(config["volume_ratio"])
    linear_scale, area_scale = geometry_scales(volume_ratio)
    model = dict(config.get("model_config", {}))
    end_time = float(model.get("end_time", 400.0 if "spontaneous" in arm else 80.0))
    dt = float(model.get("dt", 5e-4))
    dx = float(model.get("dx", 0.1))
    output_interval = float(model.get("output_interval", 0.1))
    activity_threshold = float(model.get("activity_threshold", 0.5))
    minimum_domain_area = float(model.get("minimum_domain_area", 0.2))
    radial_bins = int(model.get("radial_bins", 72))
    snapshot_bins = int(model.get("snapshot_bins", 8))
    requested_snapshot_times = np.asarray(
        model.get("snapshot_times", []), dtype=float
    )
    seed = int(config.get("seed") if config.get("seed") is not None else 1)
    stochastic = arm != "fixed_controlled"
    deformable = arm == "deformable_spontaneous"
    baseline = PHASE_PARAMETER_SETS["4F"]
    parameters = replace(baseline, area0=baseline.area0 * area_scale)
    radius = math.sqrt(parameters.area0 / math.pi)
    box_half_width = float(model.get("baseline_box_half_width", 10.0)) * linear_scale
    x = np.arange(-box_half_width, box_half_width + 0.5 * dx, dx)
    y = x.copy()
    X, Y = np.meshgrid(x, y)
    phi = _initial_phi(X, Y, parameters)
    U = np.full_like(phi, 1.0 / parameters.gamma)
    V = np.zeros_like(phi)
    trigger: dict[str, float] | None = None
    if arm == "fixed_controlled":
        relative_x = float(model.get("trigger_relative_x", -0.5))
        relative_y = float(model.get("trigger_relative_y", 0.0))
        target_x, target_y = relative_x * radius, relative_y * radius
        inside = np.flatnonzero(phi.ravel() >= 0.5)
        distances = (X.ravel()[inside] - target_x) ** 2 + (Y.ravel()[inside] - target_y) ** 2
        center = int(inside[int(np.argmin(distances))])
        j, i = np.unravel_index(center, phi.shape)
        amplitude = float(model.get("trigger_amplitude", baseline.noise_sigma))
        width = float(model.get("trigger_width", baseline.noise_d))
        transfer = amplitude * np.exp(-((X - X[j, i]) ** 2 + (Y - Y[j, i]) ** 2)
                                      / (2.0 * width * width))
        U -= transfer
        V += transfer
        trigger = {
            "target_x": target_x, "target_y": target_y,
            "actual_x": float(X[j, i]), "actual_y": float(Y[j, i]),
            "relative_x": relative_x, "relative_y": relative_y,
            "amplitude": amplitude, "width": width,
        }
    q_u, q_v = phi * U, phi * V
    workspace = tuple(np.zeros_like(phi) for _ in range(15))
    values = (
        parameters.D_u, parameters.D_v, parameters.alpha, parameters.beta,
        parameters.K_k, parameters.K_p, parameters.S, parameters.gamma,
        parameters.mu, parameters.chi_u, parameters.chi_v, parameters.tau,
        parameters.epsilon, parameters.eta, parameters.area0,
        parameters.area_penalty, parameters.protrusion_a, parameters.retraction_b,
    )
    steps = int(np.ceil(end_time / dt))
    every = max(1, int(round(output_interval / dt)))
    save_steps = np.arange(0, steps + 1, every, dtype=int)
    if save_steps[-1] != steps:
        save_steps = np.r_[save_steps, steps]
    times = np.minimum(save_steps * dt, end_time)
    n_frames = times.size
    trajectory = {key: np.empty(n_frames, dtype=np.float64) for key in TRAJECTORY_KEYS}
    radial_profiles = np.empty((n_frames, radial_bins), dtype=np.float32)
    boundary_activity_profiles = np.empty_like(radial_profiles)
    tracker = WaveTracker(phi.shape, dx, minimum_domain_area)
    initial_support = phi >= 0.5
    visited = np.zeros_like(initial_support)
    activation_passes = np.zeros_like(phi, dtype=np.uint16)
    previous_active = np.zeros_like(initial_support)
    firing_rule = str(model.get("firing_rule", "per_area"))
    if firing_rule not in {"per_area", "per_cell_reference"}:
        raise ValueError(
            "firing_rule must be 'per_area' or 'per_cell_reference'"
        )
    reference_event_area = float(model.get("reference_event_area", baseline.area0))
    if reference_event_area <= 0:
        raise ValueError("reference_event_area must be positive")
    # Preserve the legacy single-stream process byte-for-byte for the primary
    # per-area campaign.  The per-cell sensitivity deliberately separates the
    # arrival/amplitude and spatial-placement streams.  Thus paired sizes have
    # identical attempted-event times, counts, and amplitudes even though the
    # eligible spatial supports differ.
    if firing_rule == "per_cell_reference":
        event_sequence, position_sequence = np.random.SeedSequence(seed).spawn(2)
        event_rng = np.random.default_rng(event_sequence)
        position_rng = np.random.default_rng(position_sequence)
    else:
        event_rng = np.random.default_rng(seed)
        position_rng = event_rng
    events: list[tuple[float, float, float, float]] = []
    inside_fixed = np.flatnonzero(initial_support.ravel())
    bin_best: list[tuple[float, float, np.ndarray, np.ndarray] | None] = [None] * snapshot_bins
    requested_snapshots: list[tuple[float, float, np.ndarray, np.ndarray]] = []
    phase_candidate_frames = 0
    phase_candidate_max = 0

    def record(index: int, time: float) -> None:
        nonlocal previous_active, phase_candidate_frames, phase_candidate_max
        current_U = q_u / np.maximum(phi, 1e-3)
        current_V = q_v / np.maximum(phi, 1e-3)
        geometry, support, boundary, radial, boundary_pip3 = _frame_geometry(
            phi, current_V, X, Y, dx, radial_bins
        )
        active_raw = support & (current_V >= activity_threshold)
        count, dominant, boundary_fraction, active = tracker.update(
            time, active_raw, support, X, Y
        )
        activation_passes[:] += (active & ~previous_active).astype(np.uint16)
        visited[:] |= active
        previous_active = active.copy()
        active_pixels = int(np.count_nonzero(active))
        support_pixels = max(1, int(np.count_nonzero(support)))
        active_weight = np.maximum(current_V, 0.0) * np.maximum(phi, 0.0)
        active_total = float(np.sum(active_weight))
        mx = float(np.sum(active_weight * (X - geometry["centroid_x"]))
                   / max(active_total, 1e-12))
        my = float(np.sum(active_weight * (Y - geometry["centroid_y"]))
                   / max(active_total, 1e-12))
        effective_radius = math.sqrt(geometry["area"] / np.pi)
        moment = min(1.0, math.hypot(mx, my) / max(effective_radius, 1e-12))
        singularities = _phase_singularities(
            phi, current_U, current_V, 1.0 / parameters.gamma,
            max(parameters.noise_sigma, 1e-12),
        )
        if singularities:
            phase_candidate_frames += 1
        phase_candidate_max = max(phase_candidate_max, singularities)
        values_now = {
            **geometry,
            "active_area_fraction": active_pixels / support_pixels,
            "active_patch_count": float(count),
            "maximum_v": float(np.max(current_V[support])) if np.any(support) else 0.0,
            "dominant_patch_fraction": dominant,
            "activity_moment": moment,
            "activity_angle": math.atan2(my, mx) if moment > 1e-8 else np.nan,
            "boundary_active_fraction": boundary_fraction,
            "phase_singularity_count": float(singularities),
            # Lab-frame pixels do not identify material points once the cell
            # translates or deforms.  These cumulative spatial-coverage
            # quantities are therefore defined only for fixed domains.
            "visited_area_fraction": (
                np.nan if deformable else np.count_nonzero(visited) / support_pixels
            ),
            "reexcited_area_fraction": (
                np.nan if deformable else
                np.count_nonzero(activation_passes >= 2) / support_pixels
            ),
        }
        for key in TRAJECTORY_KEYS:
            trajectory[key][index] = values_now[key]
        radial_profiles[index] = radial
        boundary_activity_profiles[index] = boundary_pip3
        score = values_now["active_area_fraction"] + 0.15 * count + 0.1 * singularities
        if requested_snapshot_times.size:
            if np.any(np.isclose(time, requested_snapshot_times,
                                 atol=0.51 * output_interval, rtol=0.0)):
                requested_snapshots.append((
                    score, time, phi.astype(np.float32).copy(),
                    current_V.astype(np.float32).copy(),
                ))
        else:
            bin_index = min(
                snapshot_bins - 1,
                int(time / max(end_time, dt) * snapshot_bins),
            )
            current = bin_best[bin_index]
            if current is None or score > current[0]:
                bin_best[bin_index] = (
                    score, time, phi.astype(np.float32).copy(),
                    current_V.astype(np.float32).copy(),
                )

    record(0, 0.0)
    next_save = 1
    for step in range(1, steps + 1):
        try:
            if deformable:
                _phase_step_fast(phi, q_u, q_v, dx, dt, values, workspace)
            else:
                _chemical_step_fast(phi, q_u, q_v, dx, dt, values, workspace)
        except ZeroDivisionError as exc:
            raise RuntimeError(
                "Taniguchi phase field collapsed to zero area at "
                f"t={step * dt:g} (arm={arm}, volume={volume_ratio:g}, "
                f"seed={seed})"
            ) from exc
        if stochastic:
            area = float(np.sum(phi) * dx * dx)
            event_area = (
                area if firing_rule == "per_area" else reference_event_area
            )
            n_events = int(event_rng.poisson(parameters.theta * event_area * dt))
            for _ in range(n_events):
                inside = (np.flatnonzero(phi.ravel() >= 0.5)
                          if deformable else inside_fixed)
                flat_index = int(inside[position_rng.integers(inside.size)])
                j, i = np.unravel_index(flat_index, phi.shape)
                amplitude = float(event_rng.exponential(parameters.noise_sigma))
                transfer = amplitude * np.exp(-((X - X[j, i]) ** 2
                                                + (Y - Y[j, i]) ** 2)
                                               / (2.0 * parameters.noise_d ** 2))
                q_u -= phi * transfer
                q_v += phi * transfer
                events.append((min(step * dt, end_time), float(X[j, i]),
                               float(Y[j, i]), amplitude))
        if next_save < n_frames and step == save_steps[next_save]:
            record(next_save, float(times[next_save]))
            next_save += 1

    tracks = tracker.finish()
    track_rows = _track_rows(tracks, 2.0 * radius, output_interval)
    area_time = float(np.trapezoid(trajectory["area"], times))
    lifetimes = [float(row["lifetime"]) for row in track_rows]
    speeds = [float(row["front_speed"]) for row in track_rows]
    travels = [float(row["travel_per_diameter"]) for row in track_rows]
    migration = _migration_summary(
        times, trajectory["centroid_x"], trajectory["centroid_y"], radius
    )
    summary: dict[str, Any] = {
        "campaign": str(config.get("campaign", "taniguchi_geometry_v1")),
        "arm": arm,
        "volume_ratio": volume_ratio,
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "seed": seed if stochastic else None,
        "stochastic": stochastic,
        "deformable": deformable,
        "end_time": end_time,
        "dt": dt,
        "dx": dx,
        "output_interval": output_interval,
        "parameter_set": "4F",
        "boundary_mode": "deformable" if deformable else "frozen_diffuse_circle",
        "target_area": parameters.area0,
        "initial_radius": radius,
        "firing_rule": firing_rule,
        "reference_event_area": (
            reference_event_area if firing_rule == "per_cell_reference" else None
        ),
        "prescribed_event_rate_per_cell_time": (
            parameters.theta * reference_event_area
            if firing_rule == "per_cell_reference" else None
        ),
        "event_count": len(events),
        "event_rate_per_cell_time": len(events) / end_time,
        "event_rate_per_area_time": len(events) / area_time,
        "successful_wave_initiation_count": tracker.births,
        "successful_wave_initiation_rate_per_cell_time": tracker.births / end_time,
        "successful_wave_initiation_rate_per_area_time": tracker.births / area_time,
        "mean_simultaneous_active_patches": float(np.mean(trajectory["active_patch_count"])),
        "maximum_simultaneous_active_patches": int(np.max(trajectory["active_patch_count"])),
        "fraction_time_active": float(np.mean(trajectory["active_patch_count"] > 0)),
        "mean_active_area_fraction": float(np.mean(trajectory["active_area_fraction"])),
        "maximum_active_area_fraction": float(np.max(trajectory["active_area_fraction"])),
        "visited_area_fraction": float(trajectory["visited_area_fraction"][-1]),
        "reexcited_area_fraction": float(trajectory["reexcited_area_fraction"][-1]),
        "wave_track_count": len(track_rows),
        "mean_wave_lifetime": float(np.mean(lifetimes)) if lifetimes else None,
        "median_wave_lifetime": float(np.median(lifetimes)) if lifetimes else None,
        "mean_wave_speed": float(np.mean(speeds)) if speeds else None,
        "median_wave_speed": float(np.median(speeds)) if speeds else None,
        "mean_wave_travel_per_diameter": float(np.mean(travels)) if travels else None,
        "boundary_contact_count": sum(bool(row["boundary_contact"]) for row in track_rows),
        "reflection_candidate_count": sum(bool(row["reflection_candidate"]) for row in track_rows),
        "extinction_candidate_count": sum(bool(row["extinction_candidate"]) for row in track_rows),
        "collision_merge_count": tracker.merges,
        "wave_split_count": tracker.splits,
        "phase_singularity_candidate_frames": phase_candidate_frames,
        "maximum_phase_singularity_candidates": phase_candidate_max,
        "fraction_frames_with_phase_singularity_candidate": phase_candidate_frames / n_frames,
        "mean_aspect_ratio": float(np.mean(trajectory["aspect_ratio"])),
        "sd_aspect_ratio": float(np.std(trajectory["aspect_ratio"], ddof=1)),
        "mean_circularity": float(np.mean(trajectory["circularity"])),
        "sd_circularity": float(np.std(trajectory["circularity"], ddof=1)),
        "mean_perimeter": float(np.mean(trajectory["perimeter"])),
        "sd_perimeter": float(np.std(trajectory["perimeter"], ddof=1)),
        "mean_protrusion_count": float(np.mean(trajectory["protrusion_count"])),
        "maximum_protrusion_count": int(np.max(trajectory["protrusion_count"])),
        "mean_dominant_patch_fraction_when_active": (
            float(np.mean(trajectory["dominant_patch_fraction"]
                          [trajectory["active_patch_count"] > 0]))
            if np.any(trajectory["active_patch_count"] > 0) else None
        ),
        "pip3_to_future_boundary_extension_correlation": (
            _lagged_boundary_coupling(times, radial_profiles,
                                      boundary_activity_profiles)
            if deformable else None
        ),
        "trigger": trigger,
        **migration,
    }
    if trigger is not None:
        active_indices = np.flatnonzero(trajectory["active_patch_count"] > 0)
        boundary_indices = np.flatnonzero(trajectory["boundary_active_fraction"] > 0)
        far_boundary = (X >= 0.85 * radius) & initial_support
        # The compact tracker does not retain the active mask at every frame;
        # reaching most of the disk and contacting its boundary are therefore
        # the predeclared crossing proxies.
        summary.update({
            "controlled_first_active_time": (
                float(times[active_indices[0]]) if active_indices.size else None
            ),
            "controlled_time_to_boundary": (
                float(times[boundary_indices[0]]) if boundary_indices.size else None
            ),
            "controlled_crossed_cell": bool(np.any(visited & far_boundary)),
        })

    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "trajectory.npz",
        format_version=np.asarray("taniguchi_geometry_compact_v1"),
        times=times,
        radial_profiles=radial_profiles,
        boundary_pip3_profiles=boundary_activity_profiles,
        **trajectory,
    )
    event_array = np.asarray(events, dtype=float).reshape((-1, 4))
    np.savez_compressed(
        output / "events_and_tracks.npz",
        event_time=event_array[:, 0] if event_array.size else np.empty(0),
        event_x=event_array[:, 1] if event_array.size else np.empty(0),
        event_y=event_array[:, 2] if event_array.size else np.empty(0),
        event_amplitude=event_array[:, 3] if event_array.size else np.empty(0),
        track_json=np.asarray(json.dumps(_json_ready(track_rows))),
    )
    selected = (
        sorted(requested_snapshots, key=lambda item: item[1])
        if requested_snapshot_times.size else
        [item for item in bin_best if item is not None]
    )
    np.savez_compressed(
        output / "snapshots.npz",
        times=np.asarray([item[1] for item in selected]),
        phi=np.stack([item[2] for item in selected]),
        V=np.stack([item[3] for item in selected]),
        x=x,
        y=y,
    )
    (output / "summary.json").write_text(
        json.dumps(_json_ready(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output / "resolved_config.json").write_text(
        json.dumps(_json_ready({
            **config,
            "resolved": {
                "linear_scale": linear_scale, "area_scale": area_scale,
                "target_area": parameters.area0, "initial_radius": radius,
                "box_half_width": box_half_width, "stochastic": stochastic,
                "deformable": deformable, "firing_rule": firing_rule,
                "reference_event_area": (
                    reference_event_area
                    if firing_rule == "per_cell_reference" else None
                ),
            },
        }), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
