from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from polarity1d.contract import Grid1D
from polarity1d.jilkine import _ModelMetadata

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


def make_model():
    return LEGI()
