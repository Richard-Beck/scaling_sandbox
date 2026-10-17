"""Stimulus and periodic membrane-profile utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def wrapped_gaussian_profile(
    membrane_fraction: np.ndarray,
    center_fraction: float,
    width_fraction: float,
) -> np.ndarray:
    """Return a unit-height Gaussian on a closed unit-period membrane."""
    if width_fraction <= 0:
        raise ValueError("width_fraction must be positive")
    center = center_fraction % 1.0
    separation = np.abs(np.asarray(membrane_fraction, dtype=float) - center)
    separation = np.minimum(separation, 1.0 - separation)
    return np.exp(-0.5 * (separation / width_fraction) ** 2)


@dataclass(frozen=True)
class NoStimulus:
    """A zero stimulus with the common in-place interface."""

    def values(self, time: float, membrane_fraction: np.ndarray,
               out: np.ndarray) -> None:
        _ = time, membrane_fraction
        out.fill(0.0)


@dataclass(frozen=True)
class WrappedGaussianStimulus:
    """A time-windowed Gaussian stimulus around the periodic membrane."""

    amplitude: float
    center_fraction: float
    width_fraction: float
    start: float = 0.0
    stop: float = np.inf

    def values(self, time: float, membrane_fraction: np.ndarray,
               out: np.ndarray) -> None:
        if self.start <= time < self.stop:
            out[:] = self.amplitude * wrapped_gaussian_profile(
                membrane_fraction, self.center_fraction, self.width_fraction
            )
        else:
            out.fill(0.0)


@dataclass(frozen=True)
class WrappedTopHatStimulus:
    """A time-windowed constant stimulus on a periodic membrane interval."""

    amplitude: float
    center_fraction: float
    width_fraction: float
    start: float = 0.0
    stop: float = np.inf

    def __post_init__(self) -> None:
        if not 0.0 < self.width_fraction <= 1.0:
            raise ValueError("width_fraction must be in (0, 1]")

    def values(self, time: float, membrane_fraction: np.ndarray,
               out: np.ndarray) -> None:
        if not self.start <= time < self.stop:
            out.fill(0.0)
            return
        center = self.center_fraction % 1.0
        separation = np.abs(np.asarray(membrane_fraction, dtype=float) - center)
        separation = np.minimum(separation, 1.0 - separation)
        tolerance = 8.0 * np.finfo(float).eps
        out[:] = self.amplitude * (
            separation <= 0.5 * self.width_fraction + tolerance
        )
