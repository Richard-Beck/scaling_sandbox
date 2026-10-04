"""Model-independent polarity phenotypes for WP1 comparisons."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.signal import find_peaks

from .contract import Trajectory


ENDPOINT_STATES = (
    "unpolarized",
    "single_front",
    "multifront",
    "polarized_no_resolved_peak",
)


def diffusive_time(time: np.ndarray | float, length: float, *, reference_length: float = 10.0):
    """Map physical time to the reference-length diffusion clock ``t*(Lref/L)^2``."""
    if length <= 0 or reference_length <= 0:
        raise ValueError("lengths must be positive")
    scaled = np.asarray(time) * (reference_length / length) ** 2
    return float(scaled) if scaled.ndim == 0 else scaled


def spatial_contrast(profile: np.ndarray) -> float:
    profile = np.asarray(profile, dtype=float)
    maximum = float(np.max(profile))
    minimum = float(np.min(profile))
    return (maximum - minimum) / max(abs(maximum) + abs(minimum), 1e-12)


def profile_domain_count(
    profile: np.ndarray,
    *,
    relative_prominence: float = 0.08,
    min_distance_fraction: float = 0.12,
) -> int:
    """Count interior or boundary domains on a non-periodic interval."""
    values = np.asarray(profile, dtype=float)
    if values.ndim != 1 or len(values) < 5:
        raise ValueError("profile must be one-dimensional with at least five values")
    amplitude = float(np.max(values) - np.min(values))
    scale = max(float(np.max(np.abs(values))), 1.0)
    if amplitude <= 1e-8 * scale:
        return 0
    prominence = max(relative_prominence * amplitude, 1e-12)
    distance = max(1, round(min_distance_fraction * len(values)))
    peaks, _ = find_peaks(values, prominence=prominence, distance=distance)
    selected = list(peaks)
    height_threshold = float(np.min(values) + relative_prominence * amplitude)
    if values[0] > values[1] and values[0] >= height_threshold:
        selected.append(0)
    if values[-1] > values[-2] and values[-1] >= height_threshold:
        selected.append(len(values) - 1)
    return len(selected)


def _first_sustained(
    condition: np.ndarray,
    times: np.ndarray,
    *,
    after: float,
    duration: float,
) -> float | None:
    start = int(np.searchsorted(times, after, side="left"))
    for index in range(start, len(times)):
        stop_time = times[index] + duration
        stop = int(np.searchsorted(times, stop_time, side="left"))
        if stop < len(times) and np.all(condition[index : stop + 1]):
            return float(times[index])
    return None


def _endpoint_states(
    polarized: np.ndarray,
    domains: np.ndarray,
) -> np.ndarray:
    """Classify every frame without assigning biological meaning to its winner."""
    states = np.full(len(polarized), "unpolarized", dtype=object)
    states[polarized & (domains == 1)] = "single_front"
    states[polarized & (domains >= 2)] = "multifront"
    states[polarized & (domains == 0)] = "polarized_no_resolved_peak"
    return states


def _state_is_sustained_at_end(
    states: np.ndarray,
    times: np.ndarray,
    state: str,
    duration: float,
) -> bool:
    start = max(float(times[0]), float(times[-1]) - duration)
    selected = times >= start
    return bool(np.any(selected) and np.all(states[selected] == state))


def phenotype_summary(
    trajectory: Trajectory,
    *,
    analysis_start: float = 0.0,
    polarization_threshold: float = 0.1,
    orientation_threshold: float = 0.05,
    sustained_duration: float = 10.0,
    reversal_time: float | None = None,
) -> dict[str, Any]:
    """Score polarization, competition, and optional reversal uniformly.

    Resolution is only defined after a sustained multifront state has occurred.
    Unresolved event times are accompanied by explicit censoring/status fields;
    ``None`` therefore never silently means both "did not occur" and "was not
    observable".
    """
    if sustained_duration < 0:
        raise ValueError("sustained_duration must be nonnegative")
    if not 0 <= polarization_threshold <= 1 or not 0 <= orientation_threshold <= 1:
        raise ValueError("polarization/orientation thresholds must be in [0, 1]")
    values = trajectory.observable
    times = trajectory.time
    relative = (trajectory.x - trajectory.x[0]) / (trajectory.x[-1] - trajectory.x[0])
    left_region = relative <= 0.2
    right_region = relative >= 0.8
    maxima = np.max(values, axis=1)
    minima = np.min(values, axis=1)
    contrasts = (maxima - minima) / np.maximum(np.abs(maxima) + np.abs(minima), 1e-12)
    domains = np.asarray([profile_domain_count(profile) for profile in values], dtype=int)
    left = np.mean(values[:, left_region], axis=1)
    right = np.mean(values[:, right_region], axis=1)
    orientation = (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)
    polarized = contrasts >= polarization_threshold
    states = _endpoint_states(polarized, domains)
    establishment = _first_sustained(
        polarized, times, after=analysis_start, duration=sustained_duration
    )
    multifront = polarized & (domains >= 2)
    multifront_establishment = _first_sustained(
        multifront, times, after=analysis_start, duration=sustained_duration
    )
    resolution = None
    if multifront_establishment is not None:
        resolution = _first_sustained(
            states == "single_front",
            times,
            after=multifront_establishment,
            duration=sustained_duration,
        )
    active_window = times >= analysis_start
    if np.count_nonzero(active_window) > 1:
        coexistence = float(
            np.trapezoid(multifront[active_window].astype(float), times[active_window])
        )
    else:
        coexistence = 0.0
    winner = "none"
    if polarized[-1]:
        winner = (
            "right"
            if orientation[-1] > orientation_threshold
            else ("left" if orientation[-1] < -orientation_threshold else "center")
        )
    competition_status = "not_established"
    if multifront_establishment is not None:
        competition_status = "resolved" if resolution is not None else "right_censored"
    endpoint_state = str(states[-1])
    endpoint_state_sustained = _state_is_sustained_at_end(
        states, times, endpoint_state, sustained_duration
    )
    persistent_single_front = bool(
        multifront_establishment is not None
        and resolution is not None
        and endpoint_state == "single_front"
        and endpoint_state_sustained
    )
    transient_resolution_only = bool(
        multifront_establishment is not None
        and resolution is not None
        and not persistent_single_front
    )
    persistent_unresolved_multifront = bool(
        multifront_establishment is not None
        and resolution is None
        and endpoint_state == "multifront"
        and endpoint_state_sustained
    )
    if multifront_establishment is None:
        effective_competition_status = "not_established"
    elif persistent_single_front:
        effective_competition_status = "resolved_persistent"
    elif transient_resolution_only:
        effective_competition_status = f"resolved_transient_to_{endpoint_state}"
    elif persistent_unresolved_multifront:
        effective_competition_status = "unresolved_persistent_multifront"
    else:
        effective_competition_status = f"unresolved_to_{endpoint_state}"
    result: dict[str, Any] = {
        "model": trajectory.model,
        "observable": trajectory.observable_name,
        "observation_end_time_s": float(times[-1]),
        "endpoint_state": endpoint_state,
        "endpoint_state_sustained": endpoint_state_sustained,
        "final_polarized": bool(polarized[-1]),
        "ever_polarized": bool(np.any(polarized[active_window])),
        "final_spatial_contrast": float(contrasts[-1]),
        "maximum_spatial_contrast": float(np.max(contrasts)),
        "final_domain_count": int(domains[-1]),
        "maximum_domain_count": int(np.max(domains)),
        "front_establishment_time_s": establishment,
        "front_establishment_censored": establishment is None,
        "competition_established": multifront_establishment is not None,
        "competition_establishment_time_s": multifront_establishment,
        "competition_resolution_time_s": resolution,
        "competition_resolution_latency_s": (
            None
            if resolution is None or multifront_establishment is None
            else resolution - multifront_establishment
        ),
        "competition_resolution_censored": (
            multifront_establishment is not None and resolution is None
        ),
        "competition_status": competition_status,
        "competition_effective_status": effective_competition_status,
        "competition_resolution_persistent_at_end": persistent_single_front,
        "competition_resolution_transient_only": transient_resolution_only,
        "competition_unresolved_persistent_multifront": persistent_unresolved_multifront,
        "coexistence_duration_s": coexistence,
        "winning_side": winner,
        "final_orientation": float(orientation[-1]),
        "backend": trajectory.metadata.get("backend"),
    }
    length = trajectory.metadata.get("length_um")
    if length is not None:
        result["front_establishment_time_over_L2"] = (
            None if establishment is None else establishment / float(length) ** 2
        )
        result["competition_resolution_time_over_L2"] = (
            None if resolution is None else resolution / float(length) ** 2
        )
    if reversal_time is not None:
        result.update(
            {
                "pre_reversal_time_s": None,
                "pre_reversal_state": None,
                "pre_reversal_polarized": None,
                "pre_reversal_domain_count": None,
                "pre_reversal_spatial_contrast": None,
                "pre_reversal_orientation": None,
                "orientation_reversal_conditioned": False,
                "orientation_reversal_success": False,
                "orientation_reversal_status": "not_conditioned",
                "orientation_reversal_time_s": None,
                "orientation_reversal_latency_s": None,
                "orientation_reversal_censored": False,
                "reversal_conditioned": False,
                "reversal_success": False,
                "reversal_status": "not_conditioned",
                "reversal_time_s": None,
                "reversal_latency_s": None,
                "old_front_decay_time_s": None,
                "old_front_decay_latency_s": None,
                "new_front_establishment_time_s": None,
                "new_front_establishment_latency_s": None,
                "reversal_censored": False,
            }
        )
        before = np.flatnonzero(times < reversal_time)
        if before.size:
            before_index = int(before[-1])
            old_sign = np.sign(orientation[before_index])
            result.update(
                {
                    "pre_reversal_time_s": float(times[before_index]),
                    "pre_reversal_state": str(states[before_index]),
                    "pre_reversal_polarized": bool(polarized[before_index]),
                    "pre_reversal_domain_count": int(domains[before_index]),
                    "pre_reversal_spatial_contrast": float(contrasts[before_index]),
                    "pre_reversal_orientation": float(orientation[before_index]),
                }
            )
            conditioned = bool(
                polarized[before_index]
                and domains[before_index] <= 1
                and abs(orientation[before_index]) >= orientation_threshold
            )
            orientation_conditioned = bool(
                polarized[before_index]
                and abs(orientation[before_index]) >= orientation_threshold
            )
            result["reversal_conditioned"] = conditioned
            result["orientation_reversal_conditioned"] = orientation_conditioned
            if orientation_conditioned and old_sign != 0:
                orientation_reversed = polarized & (
                    old_sign * orientation < -orientation_threshold
                )
                orientation_reversal = _first_sustained(
                    orientation_reversed,
                    times,
                    after=reversal_time,
                    duration=sustained_duration,
                )
                result["orientation_reversal_time_s"] = orientation_reversal
                result["orientation_reversal_latency_s"] = (
                    None
                    if orientation_reversal is None
                    else orientation_reversal - reversal_time
                )
                result["orientation_reversal_success"] = (
                    orientation_reversal is not None
                )
                result["orientation_reversal_status"] = (
                    "success"
                    if orientation_reversal is not None
                    else "right_censored"
                )
                result["orientation_reversal_censored"] = (
                    orientation_reversal is None
                )
            if conditioned and old_sign != 0:
                reversed_condition = (
                    polarized
                    & (domains <= 1)
                    & (old_sign * orientation < -orientation_threshold)
                )
                reversal = _first_sustained(
                    reversed_condition,
                    times,
                    after=reversal_time,
                    duration=sustained_duration,
                )
                decay_condition = old_sign * orientation < orientation_threshold
                decay = _first_sustained(
                    decay_condition,
                    times,
                    after=reversal_time,
                    duration=sustained_duration,
                )
                result["reversal_time_s"] = reversal
                result["reversal_latency_s"] = None if reversal is None else reversal - reversal_time
                result["old_front_decay_time_s"] = decay
                result["old_front_decay_latency_s"] = (
                    None if decay is None else decay - reversal_time
                )
                result["new_front_establishment_time_s"] = reversal
                result["new_front_establishment_latency_s"] = (
                    None if reversal is None else reversal - reversal_time
                )
                result["reversal_success"] = reversal is not None
                result["reversal_status"] = (
                    "success" if reversal is not None else "right_censored"
                )
                result["reversal_censored"] = reversal is None
    return result


def tracking_summary(
    trajectory: Trajectory,
    period: float,
    *,
    start: float = 0.0,
    minimum_amplitude: float = 0.05,
    maximum_lag_fraction: float = 0.25,
) -> dict[str, float | bool]:
    """Estimate dimensionless orientation gain, lag, and cue tracking.

    ``AlternatingGradient`` points left at phase zero, so the expected signed
    orientation is ``-cos(phase)``. Positive lag means the response follows the
    cue after a delay. The reported amplitude is based on normalized left/right
    orientation and is therefore comparable across model concentration scales.
    """
    if period <= 0:
        raise ValueError("period must be positive")
    if minimum_amplitude < 0 or not 0 <= maximum_lag_fraction <= 0.5:
        raise ValueError("invalid tracking reliability thresholds")
    times = trajectory.time
    mask = times >= start
    if np.count_nonzero(mask) < 8:
        raise ValueError("too few post-start samples")
    relative = (trajectory.x - trajectory.x[0]) / (trajectory.x[-1] - trajectory.x[0])
    left = np.mean(trajectory.observable[:, relative <= 0.2], axis=1)
    right = np.mean(trajectory.observable[:, relative >= 0.8], axis=1)
    response = ((right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12))[mask]
    phase = 2.0 * np.pi * (times[mask] - start) / period
    response = response - np.mean(response)
    cosine = 2.0 * np.mean(response * np.cos(phase))
    sine = 2.0 * np.mean(response * np.sin(phase))
    amplitude = float(np.hypot(cosine, sine))
    lag = float(np.arctan2(-sine, -cosine))
    target = -np.cos(phase)
    response_scale = float(np.std(response))
    correlation = (
        0.0
        if response_scale <= 1e-12
        else float(np.corrcoef(response, target)[0, 1])
    )
    lag_seconds = lag * period / (2.0 * np.pi)
    reliable = bool(
        amplitude >= minimum_amplitude
        and abs(lag_seconds) <= maximum_lag_fraction * period
        and correlation > 0
    )
    return {
        "response_amplitude": amplitude,
        "phase_lag_radians": lag,
        "phase_lag_seconds": lag_seconds,
        "tracking_correlation": correlation,
        "tracking_reliable": reliable,
        "tracking_cycles_observed": float((times[mask][-1] - times[mask][0]) / period),
    }


def mechanochemical_summary(
    trajectory: Trajectory,
    *,
    stimulus_start: float,
    proximal_fraction: float = 0.2,
    distal_fraction: float = 0.2,
) -> dict[str, float]:
    """Score regional Rac/Rho propagation and mechanical response.

    The last frame strictly before ``stimulus_start`` is the baseline. Peak
    changes retain their sign but are selected by absolute magnitude. This
    avoids silently discarding reciprocal inhibition as a failed response.
    """
    if not 0 < proximal_fraction <= 0.5 or not 0 < distal_fraction <= 0.5:
        raise ValueError("regional fractions must be in (0, 0.5]")
    baseline_candidates = np.flatnonzero(trajectory.time < stimulus_start)
    baseline_index = int(baseline_candidates[-1]) if baseline_candidates.size else 0
    after = trajectory.time >= stimulus_start
    if not np.any(after):
        raise ValueError("stimulus_start falls after the trajectory")
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    proximal = relative <= proximal_fraction
    distal = relative >= 1.0 - distal_fraction
    result: dict[str, float] = {}
    for field in ("rac", "rho"):
        if field not in trajectory.fields:
            continue
        values = trajectory.fields[field]
        proximal_series = np.mean(values[:, proximal], axis=1)
        distal_series = np.mean(values[:, distal], axis=1)
        proximal_change = proximal_series - proximal_series[baseline_index]
        distal_change = distal_series - distal_series[baseline_index]
        post_indices = np.flatnonzero(after)
        proximal_peak = int(post_indices[np.argmax(np.abs(proximal_change[after]))])
        distal_peak = int(post_indices[np.argmax(np.abs(distal_change[after]))])
        proximal_value = float(proximal_change[proximal_peak])
        distal_value = float(distal_change[distal_peak])
        result.update(
            {
                f"{field}_proximal_peak_change": proximal_value,
                f"{field}_proximal_peak_time_s": float(trajectory.time[proximal_peak]),
                f"{field}_distal_peak_change": distal_value,
                f"{field}_distal_peak_time_s": float(trajectory.time[distal_peak]),
                f"{field}_distal_to_proximal_abs_ratio": abs(distal_value)
                / max(abs(proximal_value), 1e-12),
                f"{field}_final_proximal_change": float(proximal_change[-1]),
                f"{field}_final_distal_change": float(distal_change[-1]),
            }
        )
    for field in ("tension", "velocity", "displacement"):
        if field in trajectory.fields:
            result[f"maximum_absolute_{field}"] = float(
                np.max(np.abs(trajectory.fields[field][after]))
            )
    return result
