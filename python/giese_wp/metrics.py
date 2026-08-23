"""Published polarization metrics for the Giese membrane solution."""

from __future__ import annotations

import numpy as np


def _above_mean_components(mask: np.ndarray) -> list[np.ndarray]:
    """Return connected true components on a periodic one-dimensional mesh."""
    mask = np.asarray(mask, dtype=bool)
    n = mask.size
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(n)]

    # Start immediately after a false entry so no component crosses the cut.
    start = (np.flatnonzero(~mask)[0] + 1) % n
    order = (start + np.arange(n)) % n
    components: list[np.ndarray] = []
    current: list[int] = []
    for index in order:
        if mask[index]:
            current.append(int(index))
        elif current:
            components.append(np.asarray(current, dtype=int))
            current = []
    if current:
        components.append(np.asarray(current, dtype=int))
    return components


def membrane_metrics(u: np.ndarray, edge_lengths: np.ndarray,
                     homogeneous_tolerance: float = 1e-5) -> dict[str, float]:
    """Compute POL and published PF on a periodic P1 membrane.

    PF's neighbourhood radius is operationalized as the connected above-mean
    cluster containing the global maximum.  Edge lengths are assigned half to
    each endpoint, matching the lumped P1 membrane measure.
    """
    u = np.asarray(u, dtype=float)
    edge_lengths = np.asarray(edge_lengths, dtype=float)
    if u.ndim != 1 or edge_lengths.shape != u.shape:
        raise ValueError("u and edge_lengths must be one-dimensional and equal-sized")
    vertex_measure = 0.5 * (edge_lengths + np.roll(edge_lengths, 1))
    perimeter = float(edge_lengths.sum())
    mean = float(np.dot(vertex_measure, u) / perimeter)
    peak_index = int(np.argmax(u))
    peak = float(u[peak_index])
    amplitude = peak - float(np.min(u))
    scale = max(abs(mean), np.finfo(float).eps)

    pol = (peak - mean) / (mean * perimeter) if mean != 0.0 else np.nan
    if amplitude <= homogeneous_tolerance * scale:
        pf = 0.0
        cluster_fraction = 1.0
        n_clusters = 0
    else:
        components = _above_mean_components(u > mean)
        containing = next((component for component in components
                           if peak_index in component), np.array([], dtype=int))
        cluster_length = float(vertex_measure[containing].sum())
        cluster_fraction = cluster_length / perimeter
        pf = 1.0 - cluster_fraction
        n_clusters = len(components)

    return {
        "membrane_mean": mean,
        "membrane_max": peak,
        "membrane_min": float(np.min(u)),
        "relative_amplitude": amplitude / mean if mean != 0.0 else np.nan,
        "POL": float(pol),
        "PF": float(pf),
        "peak_angle_index": peak_index,
        "winning_cluster_fraction": float(cluster_fraction),
        "n_above_mean_clusters": int(n_clusters),
    }


def half_mass_pf(u: np.ndarray, edge_lengths: np.ndarray,
                 minimum_samples: int = 2048,
                 samples_per_edge: int = 16) -> float:
    """Preprint/prose PF based on the smallest connected half-mass patch.

    The periodic P1 profile is sampled uniformly in physical arclength.  Among
    all connected patches containing the global maximum, this finds the
    shortest one containing at least half of the total membrane mass and
    returns ``1 - 2 * patch_length / perimeter``.  The final occupied sample is
    included fractionally, making the result substantially finer than the
    membrane-node spacing.
    """
    u = np.asarray(u, dtype=float)
    edge_lengths = np.asarray(edge_lengths, dtype=float)
    if u.ndim != 1 or edge_lengths.shape != u.shape:
        raise ValueError("u and edge_lengths must be one-dimensional and equal-sized")
    if np.any(u < 0.0):
        raise ValueError("half-mass PF requires nonnegative membrane concentration")
    perimeter = float(edge_lengths.sum())
    n_samples = max(int(minimum_samples), int(samples_per_edge * u.size))
    sample_ds = perimeter / n_samples
    sample_s = (np.arange(n_samples) + 0.5) * sample_ds
    edge_ends = np.cumsum(edge_lengths)
    edge = np.searchsorted(edge_ends, sample_s, side="right")
    edge_starts = np.concatenate(([0.0], edge_ends[:-1]))
    xi = (sample_s - edge_starts[edge]) / edge_lengths[edge]
    sampled_u = (1.0 - xi) * u[edge] + xi * u[(edge + 1) % u.size]
    sample_mass = sampled_u * sample_ds
    total_mass = float(sample_mass.sum())
    if total_mass <= np.finfo(float).eps:
        return 0.0

    peak = int(np.argmax(sampled_u))
    tiled_mass = np.tile(sample_mass, 3)
    prefix = np.concatenate(([0.0], np.cumsum(tiled_mass)))
    central_peak = peak + n_samples
    starts = np.arange(central_peak - n_samples + 1, central_peak + 1)
    targets = prefix[starts] + 0.5 * total_mass
    ends = np.searchsorted(prefix, targets, side="left")
    ends = np.maximum(ends, central_peak + 1)
    valid = ends <= starts + n_samples
    starts = starts[valid]
    ends = ends[valid]
    mass_before_last = prefix[ends - 1] - prefix[starts]
    remaining = 0.5 * total_mass - mass_before_last
    last_mass = tiled_mass[ends - 1]
    fraction_last = np.divide(
        remaining, last_mass, out=np.ones_like(remaining), where=last_mass > 0.0
    )
    lengths_in_samples = (ends - 1 - starts) + np.clip(fraction_last, 0.0, 1.0)
    patch_fraction = float(np.min(lengths_in_samples) / n_samples)
    return float(np.clip(1.0 - 2.0 * patch_fraction, 0.0, 1.0))


def summarize_trajectory(times: np.ndarray, metrics: list[dict[str, float]],
                         final_contrast_threshold: float = 0.05) -> dict[str, float]:
    """Summarize a run, including the paper's 90%-of-maximal-PF time."""
    times = np.asarray(times, dtype=float)
    pf = np.asarray([row["PF"] for row in metrics], dtype=float)
    final_relative_amplitude = float(metrics[-1]["relative_amplitude"])
    maximal_pf = float(np.nanmax(pf))
    target = 0.9 * maximal_pf
    reached = np.flatnonzero(pf >= target - 1e-14)
    time_90 = float(times[reached[0]]) if reached.size else np.nan
    return {
        "maximal_PF": maximal_pf,
        "final_PF": float(pf[-1]),
        "final_POL": float(metrics[-1]["POL"]),
        "final_relative_amplitude": final_relative_amplitude,
        "polarization_time_90_max_PF_s": time_90,
        "final_cluster_count": int(metrics[-1]["n_above_mean_clusters"]),
        "polarized_final": bool(final_relative_amplitude >= final_contrast_threshold),
        "unique_site_final": bool(
            final_relative_amplitude >= final_contrast_threshold
            and metrics[-1]["n_above_mean_clusters"] == 1
        ),
    }
