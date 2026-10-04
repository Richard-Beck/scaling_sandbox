"""The four canonical models compared by Jilkine & Edelstein-Keshet (2011)."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np

from .contract import Grid1D


class _ModelMetadata:
    def metadata(self) -> Mapping[str, Any]:
        return asdict(self)  # type: ignore[arg-type]


@dataclass(frozen=True)
class WavePinning(_ModelMetadata):
    """Mori/Jilkine mass-conserved wave-pinning model in physical units."""

    du: float = 0.1
    dv: float = 10.0
    k0: float = 0.067
    feedback: float = 1.0
    half_saturation: float = 1.0
    decay: float = 1.0
    u0: float = 0.2683
    v0: float = 2.0

    name = "jilkine_wave_pinning"
    state_names = ("active", "inactive")
    observable_name = "active"
    preferred_step = 0.05
    nonnegative = True
    conserved_groups = ((0, 1),)

    @property
    def diffusivities(self) -> tuple[float, float]:
        return self.du, self.dv

    def initial_state(self, grid: Grid1D) -> np.ndarray:
        return np.vstack((np.full(grid.cells, self.u0), np.full(grid.cells, self.v0)))

    def reaction(self, time, x, state, stimulus) -> np.ndarray:
        del time, x
        u, v = state
        transfer = v * (
            self.k0
            + self.feedback * u * u / (self.half_saturation**2 + u * u)
            + stimulus
        ) - self.decay * u
        return np.vstack((transfer, -transfer))


@dataclass(frozen=True)
class Otsuji(_ModelMetadata):
    """Jilkine's time-rescaled Otsuji mass-conserved polarity model."""

    du: float = 0.1
    dv: float = 10.0
    a1: float = 25.0
    a2: float = 0.7
    s: float = 1.0
    total: float = 2.0
    u0: float | None = None

    name = "jilkine_otsuji"
    state_names = ("active", "inactive")
    observable_name = "active"
    preferred_step = 0.025
    nonnegative = True
    conserved_groups = ((0, 1),)

    @property
    def diffusivities(self) -> tuple[float, float]:
        return self.du, self.dv

    def homogeneous_state(self) -> tuple[float, float]:
        v = self.total / (self.a2 * self.s * self.total + 1.0) ** 2
        u = self.total - v if self.u0 is None else self.u0
        return u, self.total - u

    def initial_state(self, grid: Grid1D) -> np.ndarray:
        u, v = self.homogeneous_state()
        return np.vstack((np.full(grid.cells, u), np.full(grid.cells, v)))

    def reaction(self, time, x, state, stimulus) -> np.ndarray:
        del time, x
        u, v = state
        total = u + v
        transfer = self.a1 * (
            v - total / (self.a2 * self.s * total + 1.0) ** 2
        ) + stimulus * v
        return np.vstack((transfer, -transfer))


@dataclass(frozen=True)
class Goryachev(_ModelMetadata):
    """Reduced Goryachev Cdc42 model as parameterized in Jilkine (2011).

    The review leaves the constant cytoplasmic Cdc24--Bem1 level and initial
    Cdc42 abundance unspecified.  ``ec=1`` and the positive homogeneous fixed
    point ``u0=1`` are explicit, documented reproduction assumptions.
    """

    du: float = 0.0025
    dv: float = 10.0
    phi: float = 100.0
    c: float = 0.01733
    ec: float = 1.0
    u0: float = 1.0
    v0: float | None = None

    name = "jilkine_goryachev"
    state_names = ("active", "inactive")
    observable_name = "active"
    preferred_step = 0.1
    nonnegative = True
    conserved_groups = ((0, 1),)

    @property
    def a(self) -> float:
        return 0.33 / (1.0 + self.phi)

    @property
    def b(self) -> float:
        return 0.67 / (1.0 + self.phi)

    @property
    def diffusivities(self) -> tuple[float, float]:
        return self.du, self.dv

    def homogeneous_state(self) -> tuple[float, float]:
        v = self.c / (self.ec * (self.a * self.u0 + self.b)) if self.v0 is None else self.v0
        return self.u0, v

    def initial_state(self, grid: Grid1D) -> np.ndarray:
        u, v = self.homogeneous_state()
        return np.vstack((np.full(grid.cells, u), np.full(grid.cells, v)))

    def reaction(self, time, x, state, stimulus) -> np.ndarray:
        del time, x
        u, v = state
        transfer = self.a * self.ec * u * u * v + self.b * self.ec * u * v - self.c * u
        transfer += stimulus * v
        return np.vstack((transfer, -transfer))


@dataclass(frozen=True)
class LEGI(_ModelMetadata):
    """Local-excitation/global-inhibition sensing model."""

    da: float = 0.0
    di: float = 10.0
    dr: float = 0.0
    k_a: float = 2.0
    k_minus_a: float = 2.0
    k_i: float = 1.0
    k_minus_i: float = 1.0
    k_r: float = 1.0
    k_minus_r: float = 1.0
    r_total: float = 1.0
    a0: float = 0.0
    i0: float = 0.0
    r0: float = 0.0

    name = "jilkine_legi"
    state_names = ("activator", "inhibitor", "response")
    observable_name = "response"
    preferred_step = 0.05
    nonnegative = True
    conserved_groups = ()

    @property
    def diffusivities(self) -> tuple[float, float, float]:
        return self.da, self.di, self.dr

    def initial_state(self, grid: Grid1D) -> np.ndarray:
        return np.vstack(
            tuple(np.full(grid.cells, value) for value in (self.a0, self.i0, self.r0))
        )

    def reaction(self, time, x, state, stimulus) -> np.ndarray:
        del time, x
        activator, inhibitor, response = state
        return np.vstack(
            (
                self.k_a * stimulus - self.k_minus_a * activator,
                self.k_i * stimulus - self.k_minus_i * inhibitor,
                self.k_r * activator * (self.r_total - response)
                - self.k_minus_r * inhibitor * response,
            )
        )


JILKINE_MODELS = {
    "wave_pinning": WavePinning,
    "otsuji": Otsuji,
    "goryachev": Goryachev,
    "legi": LEGI,
}
