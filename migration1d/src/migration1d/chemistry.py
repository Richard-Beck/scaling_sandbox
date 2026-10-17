"""Reaction models and spatial cues for the 1D migration framework."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


class Cue(Protocol):
    """Dimensionless external signal sampled at material-node positions."""

    def __call__(self, time: float, position: Array) -> Array: ...


@dataclass(frozen=True)
class NoCue:
    def __call__(self, time: float, position: Array) -> Array:
        return np.zeros_like(position)


@dataclass(frozen=True)
class LinearCue:
    """A lab-frame linear cue; ``direction`` may be +1 or -1."""

    slope: float = 0.1
    direction: float = 1.0
    origin: float = 0.0

    def __call__(self, time: float, position: Array) -> Array:
        return self.direction * self.slope * (position - self.origin)


@dataclass(frozen=True)
class SwitchingCue:
    """Reverse a linear cue at ``switch_time`` to test polarity reversal."""

    switch_time: float
    slope: float = 0.1
    origin: float = 0.0

    def __call__(self, time: float, position: Array) -> Array:
        direction = 1.0 if time < self.switch_time else -1.0
        return direction * self.slope * (position - self.origin)


class Chemistry(Protocol):
    """Interface for any membrane/cytoplasm reaction network."""

    species: tuple[str, ...]
    membrane_diffusion: Array
    cytoplasm_diffusion: Array
    initial_membrane: Array
    initial_cytoplasm: Array

    def reaction(
        self,
        time: float,
        position: Array,
        membrane: Array,
        cytoplasm: Array,
        strain: Array,
        stress: Array,
    ) -> tuple[Array, Array]: ...


class TensionResponsiveChemistry(Chemistry, Protocol):
    """Optional interface required only when tension feedback is enabled."""

    def protrusive_activity(self, membrane: Array, cytoplasm: Array) -> Array:
        """Return the local non-negative Rac-equivalent source activity."""

        ...

    def tension_deactivation(
        self,
        membrane: Array,
        cytoplasm: Array,
        added_rate: Array,
    ) -> tuple[Array, Array]:
        """Return reaction corrections caused by an added local off-rate."""

        ...
