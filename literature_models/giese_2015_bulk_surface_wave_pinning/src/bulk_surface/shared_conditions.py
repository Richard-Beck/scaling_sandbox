"""Mori-like stimulation and Chiou-like competition with shared chemistry."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import pi

import numpy as np
from numba import njit

from .core import ConservationPool, InitialState, ReactionSystem, Simulator
from .examples import ExampleCase
from .geometry import FixedGeometry
from .protocols import NoStimulus, WrappedTopHatStimulus, wrapped_gaussian_profile


@dataclass(frozen=True)
class SharedChemistry:
    """Common local response and transport parameters for both conditions."""

    basal: float = 0.067
    feedback: float = 1.0
    half_level: float = 1.0
    saturation: float = 1.0
    inactivation: float = 1.0
    bulk_diffusion: float = 10.0
    surface_diffusion: float = 0.1
    initial_bulk: float = 2.0
    initial_surface: float = 0.2683

    def kernel_parameters(self) -> np.ndarray:
        return np.asarray([
            self.basal,
            self.feedback,
            self.half_level,
            self.saturation,
            self.inactivation,
        ])

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@njit(cache=True)
def shared_local_response_kernel(
    time, surface, bulk_trace, stimulus, parameters,
    surface_local_out, exchange_out,
):
    """Evaluate C[k0 + gamma*A^2/(K^2+kappa*A^2) + S] - delta*A."""
    _ = time
    surface_local_out.fill(0.0)
    basal, feedback, half_level, saturation, inactivation = parameters
    half_squared = half_level * half_level
    for node in range(surface.shape[1]):
        active = surface[0, node]
        inactive = bulk_trace[0, node]
        active_squared = active * active
        activation = basal + feedback * active_squared / (
            half_squared + saturation * active_squared
        ) + stimulus[node]
        exchange_out[0, node] = activation * inactive - inactivation * active


def baseline_geometry(*, perimeter: float = 10.0, dx: float = 0.1,
                      n_surface: int = 256) -> FixedGeometry:
    """Return the shared 1x circle, parameterized by membrane perimeter."""
    return FixedGeometry.circle(
        radius=perimeter / (2.0 * pi), dx=dx, n_surface=n_surface
    )


def _simulator(geometry: FixedGeometry, chemistry: SharedChemistry,
               stimulus) -> Simulator:
    reactions = ReactionSystem(
        bulk_species=("inactive_bulk",),
        surface_species=("active_surface",),
        parameters=chemistry.kernel_parameters(),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=shared_local_response_kernel,
        stimulus=stimulus,
    )
    pool = ConservationPool(
        "total_protein", np.asarray([1.0]), np.asarray([1.0])
    )
    return Simulator(
        geometry,
        reactions,
        bulk_diffusion=np.asarray([chemistry.bulk_diffusion]),
        surface_diffusion=np.asarray([chemistry.surface_diffusion]),
        conservation_pools=(pool,),
    )


def mori_condition(
    geometry: FixedGeometry | None = None,
    chemistry: SharedChemistry | None = None,
) -> ExampleCase:
    """Homogeneous state followed by the Mori transient four-percent cue."""
    g = geometry or baseline_geometry()
    p = chemistry or SharedChemistry()
    cue = WrappedTopHatStimulus(
        amplitude=10.0,
        center_fraction=0.0,
        width_fraction=0.04,
        start=50.0,
        stop=70.0,
    )
    initial = InitialState(
        bulk=np.full((1, g.n_bulk), p.initial_bulk),
        surface=np.full((1, g.n_surface), p.initial_surface),
    )
    return ExampleCase(_simulator(g, p, cue), initial)


def chiou_condition(
    geometry: FixedGeometry | None = None,
    chemistry: SharedChemistry | None = None,
    *,
    peak_amplitude: float = 0.8,
    peak_width_fraction: float = 0.04,
) -> ExampleCase:
    """Two unforced 60:40 surface seeds competing for the shared bulk pool."""
    if peak_amplitude <= 0.0:
        raise ValueError("peak_amplitude must be positive")
    g = geometry or baseline_geometry()
    p = chemistry or SharedChemistry()
    strong = wrapped_gaussian_profile(
        g.membrane_fraction, 0.25, peak_width_fraction
    )
    weak = wrapped_gaussian_profile(
        g.membrane_fraction, 0.75, peak_width_fraction
    )
    surface = p.initial_surface + peak_amplitude * (0.6 * strong + 0.4 * weak)
    added_surface_mass = float(
        (surface - p.initial_surface) @ g.surface_weights
    )
    bulk_level = p.initial_bulk - added_surface_mass / g.area
    if bulk_level <= 0.0:
        raise ValueError("surface seeds exhaust the shared bulk pool")
    initial = InitialState(
        bulk=np.full((1, g.n_bulk), bulk_level),
        surface=surface[None, :],
    )
    return ExampleCase(_simulator(g, p, NoStimulus()), initial)


def local_response(
    active: np.ndarray,
    inactive_trace: np.ndarray,
    stimulus: np.ndarray,
    chemistry: SharedChemistry | None = None,
) -> np.ndarray:
    """Vectorized form of the shared response for diagnostics and plotting."""
    p = chemistry or SharedChemistry()
    active_squared = np.asarray(active) ** 2
    activation = (
        p.basal
        + p.feedback * active_squared
        / (p.half_level**2 + p.saturation * active_squared)
        + np.asarray(stimulus)
    )
    return np.asarray(inactive_trace) * activation - p.inactivation * active
