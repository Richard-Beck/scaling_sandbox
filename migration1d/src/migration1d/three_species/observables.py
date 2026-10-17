"""Composite Cdc42--Rac--Rho polarity metrics and regime labels."""

from __future__ import annotations

import numpy as np


def _range_contrast(profile: np.ndarray) -> float:
    high, low = float(np.max(profile)), float(np.min(profile))
    return (high - low) / max(high + low, 1.0e-12)


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    if np.std(left) <= 1.0e-12 or np.std(right) <= 1.0e-12:
        return 0.0
    return float(np.corrcoef(left, right)[0, 1])


def _moment(profile: np.ndarray) -> float:
    coordinate = np.linspace(-1.0, 1.0, len(profile))
    return float(np.trapezoid(coordinate * profile, coordinate) /
                 max(np.trapezoid(profile, coordinate), 1.0e-12))


def _domain_count(profile: np.ndarray, contrast: float, threshold: float) -> int:
    if contrast < threshold:
        return 0
    low, high = float(np.min(profile)), float(np.max(profile))
    mask = profile > 0.5 * (low + high)
    return int(mask[0]) + int(np.sum(mask[1:] & ~mask[:-1]))


def profile_metrics(
    profiles: np.ndarray, *, polarity_threshold: float = 0.2
) -> dict[str, float | int | bool]:
    """Metrics for one ``(3, nodes)`` active-GTPase profile."""
    cdc42, rac, rho = profiles
    c_contrast, r_contrast, q_contrast = (
        _range_contrast(profile) for profile in profiles
    )
    c_moment, r_moment, q_moment = (_moment(profile) for profile in profiles)
    cr_correlation = _correlation(cdc42, rac)
    cq_correlation = _correlation(cdc42, rho)
    directionally_coherent = (
        abs(c_moment) >= 0.01
        and np.sign(c_moment) == np.sign(r_moment)
        and np.sign(c_moment) == -np.sign(q_moment)
    )
    coherent = bool(
        c_contrast >= polarity_threshold
        and r_contrast >= polarity_threshold
        and cr_correlation >= 0.5
        and cq_correlation < 0.0
        and directionally_coherent
    )
    return {
        "cdc42_range_contrast": c_contrast,
        "rac_range_contrast": r_contrast,
        "rho_range_contrast": q_contrast,
        "cdc42_rac_correlation": cr_correlation,
        "cdc42_rho_correlation": cq_correlation,
        "cdc42_moment": c_moment,
        "rac_moment": r_moment,
        "rho_moment": q_moment,
        "rac_pole_count": _domain_count(rac, r_contrast, polarity_threshold),
        "coherent_polarity": coherent,
    }


def _first_sustained_time(
    times: np.ndarray, values: np.ndarray, duration: float
) -> float | None:
    for index in np.flatnonzero(values):
        end = int(np.searchsorted(times, times[index] + duration, side="left"))
        if end < len(times) and np.all(values[index : end + 1]):
            return float(times[index])
    return None


def _longest_true_duration(times: np.ndarray, values: np.ndarray) -> float:
    longest = 0.0
    start: int | None = None
    for index, value in enumerate(values):
        if value and start is None:
            start = index
        if start is not None and (not value or index == len(values) - 1):
            end = index if value else index - 1
            longest = max(longest, float(times[end] - times[start]))
            start = None
    return longest


def classify_three_species_trajectory(
    times: np.ndarray,
    profiles: np.ndarray,
    totals: np.ndarray,
    *,
    homogeneous_growth_rate: float,
    maximum_spatial_growth_rate: float,
    polarity_threshold: float = 0.2,
    sustained_seconds: float = 10.0,
    late_window_seconds: float = 50.0,
    stationary_change: float = 0.01,
) -> dict[str, object]:
    """Classify a stimulated trajectory using the parameter-mapping criteria."""
    metrics = [profile_metrics(row, polarity_threshold=polarity_threshold) for row in profiles]
    coherent = np.asarray([row["coherent_polarity"] for row in metrics], dtype=bool)
    late = times >= max(times[0], times[-1] - late_window_seconds)
    late_profiles = profiles[late]
    scale = np.maximum(np.mean(np.abs(late_profiles), axis=(0, 2)), 1.0e-12)
    if len(late_profiles) > 1:
        changes = np.sqrt(np.mean(
            (np.diff(late_profiles, axis=0) / scale[None, :, None]) ** 2,
            axis=(1, 2),
        ))
        maximum_change = float(np.max(changes))
    else:
        maximum_change = 0.0
    final = metrics[-1]
    unstable = max(homogeneous_growth_rate, maximum_spatial_growth_rate) > 1.0e-5
    terminal_pattern = max(
        float(final["cdc42_range_contrast"]),
        float(final["rac_range_contrast"]),
        float(final["rho_range_contrast"]),
    ) >= polarity_threshold
    rac_active_fraction = float(np.mean(profiles[-1, 1]) / totals[1])
    ever_responded = bool(np.any([
        max(float(row["cdc42_range_contrast"]), float(row["rac_range_contrast"]))
        >= polarity_threshold for row in metrics
    ]))
    if maximum_change > stationary_change:
        state = "oscillatory_dynamic"
    elif bool(final["coherent_polarity"]):
        state = "spontaneously_unstable" if unstable else "finitely_polarizable"
    elif unstable:
        state = "spontaneously_unstable"
    elif rac_active_fraction >= 0.65 and not terminal_pattern:
        state = "overstimulated_homogeneous"
    elif terminal_pattern:
        state = "incoherent_pattern"
    elif ever_responded:
        state = "transient_response"
    else:
        state = "insensitive_unpolarized"
    return {
        "terminal_polarity_status": state,
        "time_to_coherent_polarity": _first_sustained_time(
            times, coherent, sustained_seconds
        ),
        "maximum_coherent_polarity_duration": _longest_true_duration(times, coherent),
        "late_maximum_normalized_change": maximum_change,
        "terminal_rac_active_fraction": rac_active_fraction,
        **final,
    }
