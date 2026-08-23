"""Diagnostics tied directly to the Taniguchi-port acceptance tests."""

from __future__ import annotations

import numpy as np

from .fixed import FixedResult


def fixed_diagnostics(result: FixedResult, threshold: float = 0.5) -> dict[str, float]:
    mask = result.grid.mask
    U = result.U[:, mask].astype(float)
    V = result.V[:, mask].astype(float)
    active_fraction = np.mean(V > threshold, axis=1)
    correlations = []
    for u, v in zip(U, V):
        if np.std(u) > 1e-8 and np.std(v) > 1e-8:
            correlations.append(float(np.corrcoef(u, v)[0, 1]))
    return {
        "event_count": float(result.event_count),
        "peak_V": float(np.max(V)),
        "peak_mean_V": float(np.max(np.mean(V, axis=1))),
        "peak_active_fraction": float(np.max(active_fraction)),
        "median_spatial_UV_correlation": (
            float(np.median(correlations)) if correlations else float("nan")
        ),
        "final_mean_V": float(np.mean(V[-1])),
    }


def directional_radius(phi: np.ndarray, x: np.ndarray, y: np.ndarray,
                       direction: str, level: float = 0.5) -> float:
    """Axis intercept of the phi level set, linearly interpolated."""
    if direction not in {"right", "left", "top", "bottom"}:
        raise ValueError("unknown direction")
    mid_x = int(np.argmin(np.abs(x)))
    mid_y = int(np.argmin(np.abs(y)))
    if direction in {"right", "left"}:
        coordinates, values = x, phi[mid_y]
        indices = np.arange(mid_x, len(x)) if direction == "right" else np.arange(mid_x, -1, -1)
    else:
        coordinates, values = y, phi[:, mid_x]
        indices = np.arange(mid_y, len(y)) if direction == "top" else np.arange(mid_y, -1, -1)
    previous = indices[0]
    for index in indices[1:]:
        if values[index] < level <= values[previous]:
            fraction = (level - values[index]) / (values[previous] - values[index])
            return float(coordinates[index] + fraction * (coordinates[previous] - coordinates[index]))
        previous = index
    return float("nan")
