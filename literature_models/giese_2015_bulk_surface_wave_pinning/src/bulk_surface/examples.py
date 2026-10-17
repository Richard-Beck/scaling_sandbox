"""Small model configurations exercising the common framework."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from .core import ConservationPool, InitialState, ReactionSystem, Simulator
from .geometry import FixedGeometry
from .protocols import NoStimulus, WrappedGaussianStimulus, wrapped_gaussian_profile


@njit(cache=True)
def _exchange_kernel(time, surface, bulk_trace, stimulus, parameters,
                     surface_local_out, exchange_out):
    _ = time, stimulus
    surface_local_out.fill(0.0)
    kon, koff = parameters
    for node in range(surface.shape[1]):
        exchange_out[0, node] = (
            kon * bulk_trace[0, node] - koff * surface[0, node]
        )


@njit(cache=True)
def _wave_pinning_kernel(time, surface, bulk_trace, stimulus, parameters,
                         surface_local_out, exchange_out):
    _ = time
    surface_local_out.fill(0.0)
    basal, feedback, half_level, inactivation = parameters
    half_squared = half_level * half_level
    for node in range(surface.shape[1]):
        active = surface[0, node]
        inactive = bulk_trace[0, node]
        activation = basal + feedback * active * active / (
            half_squared + active * active
        ) + stimulus[node]
        exchange_out[0, node] = activation * inactive - inactivation * active


@dataclass(frozen=True)
class ExampleCase:
    simulator: Simulator
    initial: InitialState


def exchange_case(
    geometry: FixedGeometry | None = None,
    *,
    kon: float = 0.8,
    koff: float = 0.3,
) -> ExampleCase:
    """One bulk and one surface species undergoing reversible exchange."""
    g = geometry or FixedGeometry.circle(radius=2.0, dx=0.25, n_surface=64)
    reactions = ReactionSystem(
        bulk_species=("bulk",),
        surface_species=("surface",),
        parameters=np.asarray([kon, koff]),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=_exchange_kernel,
        stimulus=NoStimulus(),
    )
    pool = ConservationPool("combined", np.asarray([1.0]), np.asarray([1.0]))
    simulator = Simulator(
        g,
        reactions,
        bulk_diffusion=np.asarray([1.0]),
        surface_diffusion=np.asarray([0.05]),
        conservation_pools=(pool,),
    )
    initial = InitialState(
        bulk=np.full((1, g.n_bulk), 1.2),
        surface=np.full((1, g.n_surface), 0.1),
    )
    return ExampleCase(simulator, initial)


def wave_pinning_case(
    geometry: FixedGeometry | None = None,
    *,
    stimulus: WrappedGaussianStimulus | None = None,
) -> ExampleCase:
    """Membrane-active/bulk-inactive positive-feedback test configuration."""
    g = geometry or FixedGeometry.circle(radius=2.0, dx=0.2, n_surface=96)
    cue = stimulus or WrappedGaussianStimulus(
        amplitude=1.5,
        center_fraction=0.98,
        width_fraction=0.035,
        start=1.0,
        stop=3.0,
    )
    reactions = ReactionSystem(
        bulk_species=("inactive",),
        surface_species=("active",),
        parameters=np.asarray([0.067, 1.0, 1.0, 1.0]),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=_wave_pinning_kernel,
        stimulus=cue,
    )
    pool = ConservationPool("gtpase", np.asarray([1.0]), np.asarray([1.0]))
    simulator = Simulator(
        g,
        reactions,
        bulk_diffusion=np.asarray([2.0]),
        surface_diffusion=np.asarray([0.02]),
        conservation_pools=(pool,),
    )
    initial_surface = 0.18 + 0.002 * wrapped_gaussian_profile(
        g.membrane_fraction, 0.98, 0.04
    )
    initial = InitialState(
        bulk=np.full((1, g.n_bulk), 1.5),
        surface=initial_surface[None, :],
    )
    return ExampleCase(simulator, initial)
