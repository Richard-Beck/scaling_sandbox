from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from polarity1d.contract import Grid1D
from polarity1d.jilkine import _ModelMetadata

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


def make_model():
    return Otsuji()
