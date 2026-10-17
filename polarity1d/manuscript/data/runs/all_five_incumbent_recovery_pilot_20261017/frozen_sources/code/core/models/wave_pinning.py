from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from polarity1d.contract import Grid1D
from polarity1d.jilkine import _ModelMetadata

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


def make_model():
    return WavePinning()
