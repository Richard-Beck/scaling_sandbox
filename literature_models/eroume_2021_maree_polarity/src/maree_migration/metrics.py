"""Trajectory observables for moving Eroume/Maree cells."""

from __future__ import annotations

import math

import numpy as np

from .model import SimulationResult


METRIC_COLUMNS = ("time", "centroid_x", "centroid_y", "speed", "area",
                  "perimeter", "aspect_ratio", "polarity_magnitude",
                  "polarity_direction_deg", "gradient_angle_deg",
                  "angle_to_gradient_deg")


def trajectory_metrics(result: SimulationResult) -> np.ndarray:
    rows = np.full((len(result.times), len(METRIC_COLUMNS)), np.nan, dtype=float)
    X, Y = np.meshgrid(result.x, result.y)
    dx = result.config.dx
    for index, (time, phi, fields) in enumerate(zip(
            result.times, result.phi.astype(float), result.concentrations.astype(float))):
        weight = phi
        mass = weight.sum()
        cx = float((weight * X).sum() / mass)
        cy = float((weight * Y).sum() / mass)
        xx = float((weight * (X - cx) ** 2).sum() / mass)
        yy = float((weight * (Y - cy) ** 2).sum() / mass)
        xy = float((weight * (X - cx) * (Y - cy)).sum() / mass)
        eigenvalues = np.linalg.eigvalsh(((xx, xy), (xy, yy)))
        aspect = math.sqrt(max(eigenvalues[-1], 0.0) / max(eigenvalues[0], 1e-15))
        gy, gx = np.gradient(phi, dx, edge_order=1)
        perimeter = float(np.hypot(gx, gy).sum() * dx * dx)
        # The environmental input acts directly on Cdc42, making active
        # Cdc42 the unambiguous cue-facing polarity reporter. Rac is retained
        # as the main rendered field.
        reporter = fields[1]
        reporter_mass = float((weight * reporter).sum())
        px = float((weight * reporter * (X - cx)).sum() / max(reporter_mass, 1e-15))
        py = float((weight * reporter * (Y - cy)).sum() / max(reporter_mass, 1e-15))
        polarity_magnitude = math.hypot(px, py)
        polarity_direction = math.degrees(math.atan2(py, px))
        angle = result.config.stimulus.angle_at(float(time)) if result.config.stimulus else None
        relative = np.nan
        if angle is not None and polarity_magnitude > 1e-12:
            relative = (polarity_direction - angle + 180.0) % 360.0 - 180.0
        rows[index] = (time, cx, cy, np.nan, mass * dx * dx, perimeter, aspect,
                       polarity_magnitude, polarity_direction,
                       np.nan if angle is None else angle, relative)
    if len(rows) > 1:
        displacement = np.hypot(np.diff(rows[:, 1]), np.diff(rows[:, 2]))
        rows[1:, 3] = displacement / np.diff(rows[:, 0])
        rows[0, 3] = 0.0
    else:
        rows[0, 3] = 0.0
    return rows


def material_totals(result: SimulationResult) -> np.ndarray:
    """Integral of each represented field, useful for leakage diagnostics."""
    return result.material.copy()
