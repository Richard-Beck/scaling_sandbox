from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from polarity1d.contract import Grid1D
from polarity1d.jilkine import _ModelMetadata

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


def make_model():
    return Goryachev()
