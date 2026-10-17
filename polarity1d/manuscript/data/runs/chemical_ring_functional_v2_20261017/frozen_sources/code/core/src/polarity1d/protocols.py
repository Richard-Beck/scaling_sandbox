"""Literature stimulation protocols expressed in physical coordinates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class NoStimulus:
    @property
    def event_times(self) -> tuple[float, ...]:
        return ()

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        del time, length
        return np.zeros_like(x)


def _rectangular_envelope(
    time: float, start: float, stop: float, taper: float
) -> float:
    if time < start or time >= stop:
        return 0.0
    if taper <= 0 or time <= stop - taper:
        return 1.0
    phase = np.pi * (time - (stop - taper)) / taper
    return 0.5 * (1.0 + np.cos(phase))


@dataclass(frozen=True)
class LocalizedPulse:
    """Jilkine cosine pulse at one or both ends of the domain.

    ``amplitude`` is the maximum conversion rate at the boundary.  A positive
    ``taper`` applies the half-cosine withdrawal used by Jilkine; a zero taper
    gives a rectangular optogenetic pulse.
    """

    amplitude: float = 0.1
    start: float = 0.0
    stop: float = 40.0
    taper: float = 20.0
    width_fraction: float = 0.1
    sides: tuple[str, ...] = ("left",)
    asymmetry: float = 0.0

    def __post_init__(self) -> None:
        if self.amplitude < 0 or self.start < 0 or self.stop <= self.start:
            raise ValueError("pulse amplitude/times are invalid")
        if self.taper < 0 or self.taper > self.stop - self.start:
            raise ValueError("taper must fit inside the pulse")
        if not 0 < self.width_fraction <= 0.5:
            raise ValueError("width_fraction must be in (0, 0.5]")
        if not self.sides or any(side not in {"left", "right"} for side in self.sides):
            raise ValueError("sides must contain left and/or right")
        if abs(self.asymmetry) >= 2:
            raise ValueError("absolute asymmetry must be below two")

    @property
    def event_times(self) -> tuple[float, ...]:
        events = [self.start, self.stop]
        if self.taper:
            events.append(self.stop - self.taper)
        return tuple(sorted(set(events)))

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        envelope = _rectangular_envelope(time, self.start, self.stop, self.taper)
        if envelope == 0:
            return np.zeros_like(x)
        width = self.width_fraction * length
        result = np.zeros_like(x)
        for side in self.sides:
            distance = x if side == "left" else length - x
            local = distance <= width
            shape = np.zeros_like(x)
            shape[local] = 0.5 * (1.0 + np.cos(np.pi * distance[local] / width))
            side_scale = 1.0 + (self.asymmetry / 2 if side == "left" else -self.asymmetry / 2)
            result += side_scale * shape
        return self.amplitude * envelope * result


@dataclass(frozen=True)
class FinitePulseSequence:
    """A sum of finite localized pulses with independently specified sides."""

    pulses: tuple[dict, ...]

    def __post_init__(self) -> None:
        normalized = tuple(dict(pulse) for pulse in self.pulses)
        if not normalized:
            raise ValueError("a finite pulse sequence requires at least one pulse")
        # Validate every component at construction time.
        for pulse in normalized:
            LocalizedPulse(**pulse)
        object.__setattr__(self, "pulses", normalized)

    @property
    def event_times(self) -> tuple[float, ...]:
        events = set()
        for configuration in self.pulses:
            events.update(LocalizedPulse(**configuration).event_times)
        return tuple(sorted(events))

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        result = np.zeros_like(x)
        for configuration in self.pulses:
            result += LocalizedPulse(**configuration)(time, x, length)
        return result


@dataclass(frozen=True)
class GaussianPulse:
    """Localized de Belly-style optogenetic input."""

    amplitude: float = 1.0
    start: float = 50.0
    stop: float = 150.0
    width_fraction: float = 0.2
    centers: tuple[float, ...] = (0.0,)
    asymmetry: float = 0.0

    def __post_init__(self) -> None:
        if self.amplitude < 0 or self.start < 0 or self.stop <= self.start:
            raise ValueError("pulse amplitude/times are invalid")
        if self.width_fraction <= 0:
            raise ValueError("width_fraction must be positive")
        if not self.centers or any(not 0 <= center <= 1 for center in self.centers):
            raise ValueError("centers are relative positions in [0, 1]")

    @property
    def event_times(self) -> tuple[float, ...]:
        return self.start, self.stop

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        if time < self.start or time >= self.stop:
            return np.zeros_like(x)
        width = self.width_fraction * length
        result = np.zeros_like(x)
        for index, center in enumerate(self.centers):
            scale = 1.0
            if len(self.centers) == 2:
                scale += self.asymmetry / 2 if index == 0 else -self.asymmetry / 2
            result += scale * np.exp(-0.5 * ((x - center * length) / width) ** 2)
        return self.amplitude * result


@dataclass(frozen=True)
class LinearGradient:
    """Persistent Jilkine gradient with slope in activation-rate/µm."""

    slope: float = 0.02
    direction: str = "left"
    start: float = 0.0
    stop: float = np.inf

    def __post_init__(self) -> None:
        if self.slope < 0 or self.direction not in {"left", "right"}:
            raise ValueError("slope must be nonnegative and direction left/right")
        if self.start < 0 or self.stop <= self.start:
            raise ValueError("gradient times are invalid")

    @property
    def event_times(self) -> tuple[float, ...]:
        values = [self.start]
        if np.isfinite(self.stop):
            values.append(self.stop)
        return tuple(values)

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        if time < self.start or time >= self.stop:
            return np.zeros_like(x)
        distance = length - x if self.direction == "left" else x
        return self.slope * distance


@dataclass(frozen=True)
class ReversingGradient:
    """A gradient that changes direction, optionally through a linear cross-fade."""

    slope: float = 0.02
    reversal_time: float = 100.0
    transition: float = 20.0

    def __post_init__(self) -> None:
        if self.slope < 0 or self.reversal_time <= 0 or self.transition < 0:
            raise ValueError("invalid reversal protocol")

    @property
    def event_times(self) -> tuple[float, ...]:
        return self.reversal_time, self.reversal_time + self.transition

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        left = self.slope * (length - x)
        right = self.slope * x
        if time < self.reversal_time:
            return left
        if self.transition == 0 or time >= self.reversal_time + self.transition:
            return right
        weight = (time - self.reversal_time) / self.transition
        return (1.0 - weight) * left + weight * right


@dataclass(frozen=True)
class ReversingFiniteGradient:
    """A reversing edge gradient with independently controlled peak and width.

    ``width_fraction=1`` recovers the full-domain linear stimulus used by
    :class:`ReversingGradient`; smaller values truncate the linear ramp at a
    specified fraction of the cell length.
    """

    amplitude: float = 0.1
    width_fraction: float = 1.0
    reversal_time: float = 100.0
    transition: float = 20.0

    def __post_init__(self) -> None:
        if self.amplitude < 0 or not 0 < self.width_fraction <= 1:
            raise ValueError("amplitude/width_fraction are invalid")
        if self.reversal_time <= 0 or self.transition < 0:
            raise ValueError("invalid reversal timing")

    @property
    def event_times(self) -> tuple[float, ...]:
        return self.reversal_time, self.reversal_time + self.transition

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        width = self.width_fraction * length
        left = self.amplitude * np.clip(1.0 - x / width, 0.0, 1.0)
        right = self.amplitude * np.clip(1.0 - (length - x) / width, 0.0, 1.0)
        if time < self.reversal_time:
            return left
        if self.transition == 0 or time >= self.reversal_time + self.transition:
            return right
        weight = (time - self.reversal_time) / self.transition
        return (1.0 - weight) * left + weight * right


@dataclass(frozen=True)
class AlternatingGradient:
    """Smoothly alternating input used to measure tracking bandwidth."""

    slope: float = 0.02
    period: float = 100.0
    start: float = 0.0

    def __post_init__(self) -> None:
        if self.slope < 0 or self.period <= 0 or self.start < 0:
            raise ValueError("invalid alternating gradient")

    @property
    def event_times(self) -> tuple[float, ...]:
        return (self.start,)

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        if time < self.start:
            return np.zeros_like(x)
        phase = np.cos(2.0 * np.pi * (time - self.start) / self.period)
        centered = 2.0 * x / length - 1.0
        # A nonnegative gradient with a fixed mean and alternating orientation.
        return 0.5 * self.slope * length * (1.0 - phase * centered)
