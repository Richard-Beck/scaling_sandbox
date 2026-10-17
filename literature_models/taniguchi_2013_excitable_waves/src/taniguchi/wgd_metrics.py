"""Predeclared WGD statistics for deformable Taniguchi trajectories."""

from __future__ import annotations

import math

import numpy as np
from scipy import ndimage

from .phase_field import PhaseFieldResult


def _component_count(mask: np.ndarray, minimum_pixels: int) -> int:
    labels, count = ndimage.label(mask, structure=np.ones((3, 3), dtype=int))
    if count == 0:
        return 0
    sizes = np.bincount(labels.ravel())
    return int(np.count_nonzero(sizes[1:] >= minimum_pixels))


def _shape_and_direction(
    phi: np.ndarray, activity: np.ndarray, X: np.ndarray, Y: np.ndarray,
    pixel_area: float,
) -> tuple[float, float, float, float, float]:
    """Return centroid, shape anisotropy, activity moment, and its angle."""

    weight = np.maximum(phi, 0.0)
    total = float(np.sum(weight))
    cx = float(np.sum(weight * X) / max(total, 1e-12))
    cy = float(np.sum(weight * Y) / max(total, 1e-12))
    xx = float(np.sum(weight * (X - cx) ** 2) / max(total, 1e-12))
    yy = float(np.sum(weight * (Y - cy) ** 2) / max(total, 1e-12))
    xy = float(np.sum(weight * (X - cx) * (Y - cy)) / max(total, 1e-12))
    eigenvalues = np.linalg.eigvalsh(np.asarray([[xx, xy], [xy, yy]]))
    anisotropy = float(
        (eigenvalues[-1] - eigenvalues[0]) /
        max(eigenvalues[-1] + eigenvalues[0], 1e-12)
    )

    active_weight = weight * np.maximum(activity, 0.0)
    active_total = float(np.sum(active_weight))
    radius = math.sqrt(max(float(np.sum(weight)) * pixel_area, 0.0) / math.pi)
    mx = float(np.sum(active_weight * (X - cx)) / max(active_total, 1e-12))
    my = float(np.sum(active_weight * (Y - cy)) / max(active_total, 1e-12))
    moment = min(1.0, math.hypot(mx, my) / max(radius, 1e-12))
    angle = math.atan2(my, mx) if moment > 1e-10 else float("nan")
    return cx, cy, anisotropy, moment, angle


def extract_phase_wgd_observables(
    result: PhaseFieldResult,
    *,
    activity_threshold: float = 0.5,
    minimum_domain_area: float = 0.2,
) -> dict[str, np.ndarray]:
    """Extract the fixed per-frame observables needed by the WGD analyzer.

    A PIP3 domain is an eight-connected component within ``phi >= 0.5`` whose
    PIP3 concentration is at least ``activity_threshold`` and whose area is at
    least ``minimum_domain_area``. These operational definitions are fixed by
    the WGD manifest rather than inferred from a trajectory.
    """

    if activity_threshold <= 0 or minimum_domain_area <= 0:
        raise ValueError("activity threshold and minimum domain area must be positive")
    phi = np.asarray(result.phi, dtype=float)
    pip2 = np.asarray(result.U, dtype=float)
    pip3 = np.asarray(result.V, dtype=float)
    times = np.asarray(result.times, dtype=float)
    if phi.shape != pip2.shape or phi.shape != pip3.shape or phi.ndim != 3:
        raise ValueError("phi, U, and V must be equally shaped frame stacks")
    if times.ndim != 1 or times.size != phi.shape[0]:
        raise ValueError("times must contain one value per frame")
    if not (np.isfinite(phi).all() and np.isfinite(pip2).all()
            and np.isfinite(pip3).all() and np.isfinite(times).all()):
        raise ValueError("trajectory contains non-finite values")

    dx = float(result.x[1] - result.x[0])
    pixel_area = dx * dx
    X, Y = np.meshgrid(result.x, result.y)
    minimum_pixels = max(1, int(math.ceil(minimum_domain_area / pixel_area)))
    n_frames = times.size
    areas = np.empty(n_frames)
    active_fractions = np.empty(n_frames)
    contrasts = np.empty(n_frames)
    spatial_cvs = np.empty(n_frames)
    domain_counts = np.empty(n_frames, dtype=int)
    centroids = np.empty((n_frames, 2))
    anisotropy = np.empty(n_frames)
    moments = np.empty(n_frames)
    angles = np.empty(n_frames)
    pip2_integrals = np.empty(n_frames)
    pip3_integrals = np.empty(n_frames)
    box_edge_max = np.empty(n_frames)

    for index in range(n_frames):
        support = phi[index] >= 0.5
        count = max(int(np.count_nonzero(support)), 1)
        values = pip3[index][support]
        if values.size:
            maximum = float(np.max(values))
            minimum = float(np.min(values))
            mean = float(np.mean(values))
            contrasts[index] = (maximum - minimum) / max(maximum + minimum, 1e-12)
            spatial_cvs[index] = float(np.std(values)) / max(abs(mean), 1e-12)
        else:
            contrasts[index] = 0.0
            spatial_cvs[index] = 0.0
        active = support & (pip3[index] >= activity_threshold)
        domain_counts[index] = _component_count(active, minimum_pixels)
        active_fractions[index] = np.count_nonzero(active) / count
        areas[index] = float(np.sum(phi[index]) * pixel_area)
        cx, cy, shape, moment, angle = _shape_and_direction(
            phi[index], pip3[index], X, Y, pixel_area
        )
        centroids[index] = (cx, cy)
        anisotropy[index] = shape
        moments[index] = moment
        angles[index] = angle
        pip2_integrals[index] = float(np.sum(phi[index] * pip2[index]) * pixel_area)
        pip3_integrals[index] = float(np.sum(phi[index] * pip3[index]) * pixel_area)
        box_edge_max[index] = float(np.max(np.concatenate((
            phi[index, 0, :], phi[index, -1, :],
            phi[index, :, 0], phi[index, :, -1],
        ))))

    return {
        "times": times,
        "areas": areas,
        "active_fractions": active_fractions,
        "spatial_contrasts": contrasts,
        "spatial_cvs": spatial_cvs,
        "domain_counts": domain_counts,
        "centroid_x": centroids[:, 0],
        "centroid_y": centroids[:, 1],
        "shape_anisotropy": anisotropy,
        "directional_moments": moments,
        "directional_angles": angles,
        "pip2_native_integrals": pip2_integrals,
        "pip3_native_integrals": pip3_integrals,
        "box_edge_max_phi": box_edge_max,
    }


def summarize_phase_wgd_observables(
    observables: dict[str, np.ndarray],
    *,
    area0: float,
    event_count: int,
    activity_threshold: float = 0.5,
    minimum_domain_area: float = 0.2,
) -> dict[str, float | int | str | None]:
    """Derive terminal WGD statistics from retained per-frame observables."""

    required = {
        "times", "areas", "active_fractions", "spatial_contrasts",
        "spatial_cvs", "domain_counts", "centroid_x", "centroid_y",
        "shape_anisotropy", "directional_moments", "directional_angles",
        "pip2_native_integrals", "pip3_native_integrals", "box_edge_max_phi",
    }
    missing = required - set(observables)
    if missing:
        raise ValueError(f"missing Taniguchi observables: {sorted(missing)}")
    arrays = {key: np.asarray(observables[key]) for key in required}
    times = np.asarray(arrays["times"], dtype=float)
    if times.ndim != 1 or times.size == 0:
        raise ValueError("observable times must be a nonempty vector")
    n_frames = times.size
    for key, values in arrays.items():
        if values.ndim != 1 or values.size != times.size:
            raise ValueError(f"observable {key} must have one value per time")
        if key != "directional_angles" and not np.isfinite(values).all():
            raise ValueError(f"observable {key} contains non-finite values")
    if area0 <= 0 or event_count < 0:
        raise ValueError("area0 must be positive and event_count nonnegative")

    areas = np.asarray(arrays["areas"], dtype=float)
    active_fractions = np.asarray(arrays["active_fractions"], dtype=float)
    contrasts = np.asarray(arrays["spatial_contrasts"], dtype=float)
    spatial_cvs = np.asarray(arrays["spatial_cvs"], dtype=float)
    domain_counts = np.asarray(arrays["domain_counts"], dtype=int)
    centroids = np.column_stack((arrays["centroid_x"], arrays["centroid_y"])).astype(float)
    anisotropy = np.asarray(arrays["shape_anisotropy"], dtype=float)
    moments = np.asarray(arrays["directional_moments"], dtype=float)
    angles = np.asarray(arrays["directional_angles"], dtype=float)
    pip2_integrals = np.asarray(arrays["pip2_native_integrals"], dtype=float)
    pip3_integrals = np.asarray(arrays["pip3_native_integrals"], dtype=float)

    duration = float(times[-1] - times[0]) if n_frames > 1 else 0.0
    births = int(np.sum(np.maximum(np.diff(domain_counts), 0)))
    transitions = int(np.count_nonzero(np.diff(domain_counts)))
    target = 0.9 * float(np.max(moments))
    reached = np.flatnonzero(moments >= target)
    time_to_directional_peak = (
        float(times[reached[0]]) if target > 1e-10 and reached.size else None
    )
    valid_angles = np.isfinite(angles) & (moments > 1e-3)
    if np.any(valid_angles):
        weights = moments[valid_angles]
        resultant = np.sum(weights * np.exp(1j * angles[valid_angles])) / np.sum(weights)
        directional_persistence = float(abs(resultant))
    else:
        directional_persistence = 0.0
    displacement = float(np.linalg.norm(centroids[-1] - centroids[0]))
    integrated_area_time = (
        float(np.trapezoid(areas, times)) if n_frames > 1 else 0.0
    )
    event_rate = event_count / duration if duration > 0 else 0.0
    event_rate_per_area = (
        event_count / integrated_area_time if integrated_area_time > 0 else 0.0
    )
    box_edge_max = float(np.max(arrays["box_edge_max_phi"]))
    active_frames = domain_counts > 0
    state = "excitable" if np.any(active_frames) else "homogeneous"
    if domain_counts[-1] > 1:
        state = "multi_domain"

    return {
        "state_class": state,
        "terminal_spatial_contrast": float(contrasts[-1]),
        "maximum_spatial_contrast": float(np.max(contrasts)),
        "mean_spatial_contrast": float(np.mean(contrasts)),
        "terminal_spatial_cv": float(spatial_cvs[-1]),
        "terminal_n_domains": int(domain_counts[-1]),
        "maximum_n_domains": int(np.max(domain_counts)),
        "mean_n_domains": float(np.mean(domain_counts)),
        "domain_birth_count": births,
        "domain_birth_rate": births / duration if duration > 0 else 0.0,
        "domain_count_transition_rate": transitions / duration if duration > 0 else 0.0,
        "fraction_frames_with_active_domain": float(np.mean(active_frames)),
        "terminal_active_area_fraction": float(active_fractions[-1]),
        "maximum_active_area_fraction": float(np.max(active_fractions)),
        "mean_active_area_fraction": float(np.mean(active_fractions)),
        "terminal_directional_moment": float(moments[-1]),
        "maximum_directional_moment": float(np.max(moments)),
        "mean_directional_moment": float(np.mean(moments)),
        "directional_persistence": directional_persistence,
        "time_to_directional_peak": time_to_directional_peak,
        "maximum_shape_anisotropy": float(np.max(anisotropy)),
        "terminal_shape_anisotropy": float(anisotropy[-1]),
        "centroid_displacement": displacement,
        "initial_area": float(areas[0]),
        "final_area": float(areas[-1]),
        "minimum_area": float(np.min(areas)),
        "maximum_area": float(np.max(areas)),
        "maximum_relative_target_area_error": float(
            np.max(np.abs(areas / area0 - 1.0))
        ),
        "initial_pip2_native_integral": float(pip2_integrals[0]),
        "final_pip2_native_integral": float(pip2_integrals[-1]),
        "initial_pip3_native_integral": float(pip3_integrals[0]),
        "final_pip3_native_integral": float(pip3_integrals[-1]),
        "stochastic_event_count": int(event_count),
        "stochastic_event_rate_per_time": float(event_rate),
        "stochastic_event_rate_per_area_time": float(event_rate_per_area),
        "activity_threshold": float(activity_threshold),
        "minimum_domain_area": float(minimum_domain_area),
        "box_edge_max_phi": box_edge_max,
    }


def summarize_phase_wgd(
    result: PhaseFieldResult,
    *,
    activity_threshold: float = 0.5,
    minimum_domain_area: float = 0.2,
) -> dict[str, float | int | str | None]:
    """Extract and summarize a full in-memory trajectory."""

    observables = extract_phase_wgd_observables(
        result,
        activity_threshold=activity_threshold,
        minimum_domain_area=minimum_domain_area,
    )
    return summarize_phase_wgd_observables(
        observables,
        area0=result.parameters.area0,
        event_count=result.event_count,
        activity_threshold=activity_threshold,
        minimum_domain_area=minimum_domain_area,
    )
