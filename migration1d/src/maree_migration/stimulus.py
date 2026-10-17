"""Laboratory-coordinate Cdc42 stimulus protocols."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class StimulusProtocol:
    """Define the Eroume Cdc42 activation input ``Ic = S(x, y, t)``.

    Coordinates and the gradient pivot are fixed in the laboratory frame.
    ``free`` is uniform; ``transient`` adds a finite Gaussian pulse.
    """

    mode: str = "free"
    baseline: float = 2.95
    gradient_strength: float = 0.05
    gradient_angle: float = 0.0
    rotate_at: float = 300.0
    rotate_angle: float = 90.0
    pivot_x: float = 48.0
    pivot_y: float = 48.0
    transient_duration: float = 20.0
    transient_amplitude: float = 1.0
    transient_width: float = 5.0
    transient_angle: float = 0.0
    cell_center_x: float = 48.0
    cell_center_y: float = 48.0

    def __post_init__(self) -> None:
        if self.mode not in {"free", "gradient", "rotate", "transient", "legacy"}:
            raise ValueError("mode must be free, gradient, rotate, transient, or legacy")
        if self.baseline < 0 or self.gradient_strength < 0:
            raise ValueError("baseline and gradient_strength must be non-negative")
        if self.rotate_at < 0 or self.transient_duration < 0 or self.transient_width <= 0:
            raise ValueError("stimulus times must be non-negative and width positive")

    def angle_at(self, time: float) -> float | None:
        if self.mode == "gradient":
            return self.gradient_angle
        if self.mode == "rotate":
            return self.gradient_angle if time < self.rotate_at else self.rotate_angle
        if self.mode == "transient" and time < self.transient_duration:
            return self.transient_angle
        if self.mode == "legacy" and time <= 10.0:
            return 0.0
        return None

    def evaluate(self, X: np.ndarray, Y: np.ndarray, time: float) -> np.ndarray:
        if self.mode == "legacy":
            # Exact source stimulus from vcell_model.cpp. X is not recentered.
            return 2.6 + 0.05 * X if time <= 10.0 else np.full_like(X, 2.95)
        if self.mode in {"gradient", "rotate"}:
            angle = math.radians(float(self.angle_at(time)))
            projection = ((X - self.pivot_x) * math.cos(angle)
                          + (Y - self.pivot_y) * math.sin(angle))
            # A linear activation rate is physical only while non-negative.
            return np.maximum(0.0, self.baseline + self.gradient_strength * projection)
        result = np.full_like(X, self.baseline)
        if self.mode == "transient" and time < self.transient_duration:
            angle = math.radians(self.transient_angle)
            xc = self.cell_center_x + 0.75 * self.transient_width * math.cos(angle)
            yc = self.cell_center_y + 0.75 * self.transient_width * math.sin(angle)
            result += self.transient_amplitude * np.exp(
                -((X - xc) ** 2 + (Y - yc) ** 2) / (2.0 * self.transient_width ** 2)
            )
        return result
