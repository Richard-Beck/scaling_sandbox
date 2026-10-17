"""Late-time spatial-state classification and directly observed domain counts."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ClassificationThresholds:
    amplitude: float = 0.5
    sustained_seconds: float = 10.0
    high_fraction: float = 0.5
    stationary_change: float = 0.01
    late_window_seconds: float = 50.0


def polarization_amplitude(profile: np.ndarray) -> float:
    mean = max(float(np.mean(profile)), 1.0e-12)
    return float((np.max(profile) - np.min(profile)) / mean)


def high_domain_count(profile: np.ndarray, high_fraction: float = 0.5) -> int:
    low, high = float(np.min(profile)), float(np.max(profile))
    if high - low <= 1.0e-10 * max(1.0, abs(high)):
        return 0
    mask = profile > low + high_fraction * (high - low)
    return int(mask[0]) + int(np.sum(mask[1:] & ~mask[:-1]))


def high_width_fraction(profile: np.ndarray, high_fraction: float = 0.5) -> float:
    low, high = float(np.min(profile)), float(np.max(profile))
    if high <= low:
        return 0.0
    return float(np.mean(profile > low + high_fraction * (high - low)))


def polarized_domain_count(
    profile: np.ndarray,
    *,
    amplitude_threshold: float = 0.5,
    high_fraction: float = 0.5,
) -> int:
    """Count high domains only when the profile is meaningfully polarized.

    A small numerical ripple can contain an above-midrange interval despite being
    effectively homogeneous. Requiring the same normalized-amplitude threshold
    used by the trajectory classifier makes zero mean "not polarized," one mean
    one spatially connected pole, and so on.
    """
    if polarization_amplitude(profile) < amplitude_threshold:
        return 0
    return high_domain_count(profile, high_fraction)


def _first_sustained_time(
    times: np.ndarray, values: np.ndarray, threshold: float, duration: float
) -> float | None:
    above = values >= threshold
    for index in np.flatnonzero(above):
        end = int(np.searchsorted(times, times[index] + duration, side="left"))
        if end < len(times) and np.all(above[index : end + 1]):
            return float(times[index])
    return None


def classify_trajectory(
    times: np.ndarray,
    profiles: np.ndarray,
    *,
    half_saturation: float = 1.0,
    thresholds: ClassificationThresholds | None = None,
) -> dict[str, object]:
    thresholds = thresholds or ClassificationThresholds()
    amplitudes = np.asarray([polarization_amplitude(row) for row in profiles])
    counts = np.asarray([
        polarized_domain_count(
            row,
            amplitude_threshold=thresholds.amplitude,
            high_fraction=thresholds.high_fraction,
        )
        for row in profiles
    ])
    widths = np.asarray([
        high_width_fraction(row, thresholds.high_fraction) for row in profiles
    ])
    late = times >= max(times[0], times[-1] - thresholds.late_window_seconds)
    late_profiles = profiles[late]
    scale = max(float(np.mean(np.abs(late_profiles))), 1.0e-12)
    if len(late_profiles) > 1:
        changes = np.sqrt(np.mean(np.diff(late_profiles, axis=0) ** 2, axis=1)) / scale
        maximum_change = float(np.max(changes))
    else:
        maximum_change = 0.0
    final_amplitude = float(amplitudes[-1])
    final_count = int(counts[-1])
    coarsened = bool(np.max(counts[late]) > final_count)
    if maximum_change > thresholds.stationary_change:
        state = "transient_coarsening" if coarsened else "oscillatory_dynamic"
    elif final_amplitude < thresholds.amplitude:
        state = "homogeneous_low" if float(np.mean(profiles[-1])) < half_saturation else "homogeneous_high"
    elif final_count == 1:
        state = "single_polarized"
    elif final_count == 2:
        state = "dual_polarized"
    else:
        state = "multi_peak"
    return {
        "state": state,
        "time_to_polarization": _first_sustained_time(
            times, amplitudes, thresholds.amplitude, thresholds.sustained_seconds
        ),
        "terminal_polarization_amplitude": final_amplitude,
        "terminal_pole_count": final_count,
        "terminal_high_width_fraction": float(widths[-1]),
        "late_maximum_normalized_change": maximum_change,
        "late_peak_count_maximum": int(np.max(counts[late])),
        "coarsened_in_late_window": coarsened,
        "amplitude_series": amplitudes,
        "peak_count_series": counts,
    }
