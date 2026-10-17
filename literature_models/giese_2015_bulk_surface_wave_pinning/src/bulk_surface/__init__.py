"""Optimized fixed-geometry bulk--surface reaction--diffusion framework."""

from .core import (
    ConservationPool,
    InitialState,
    ReactionSystem,
    SimulationConfig,
    SimulationResult,
    Simulator,
)
from .geometry import FixedGeometry
from .diffuse import (
    DiffuseDomainSimulator,
    DiffuseGeometry,
    DiffuseGeometryFields,
    DiffuseInitialState,
    DiffuseSimulationResult,
    PhaseFieldMechanics,
    geometry_fields,
)
from .moving_curve import (
    MovingCurveInitialState,
    MovingCurveResult,
    MovingCurveSimulator,
    MovingMembrane,
)
from .protocols import (
    NoStimulus,
    WrappedGaussianStimulus,
    WrappedTopHatStimulus,
    wrapped_gaussian_profile,
)

__all__ = [
    "ConservationPool",
    "DiffuseDomainSimulator",
    "DiffuseGeometry",
    "DiffuseGeometryFields",
    "DiffuseInitialState",
    "DiffuseSimulationResult",
    "FixedGeometry",
    "InitialState",
    "MovingCurveInitialState",
    "MovingCurveResult",
    "MovingCurveSimulator",
    "MovingMembrane",
    "NoStimulus",
    "PhaseFieldMechanics",
    "ReactionSystem",
    "SimulationConfig",
    "SimulationResult",
    "Simulator",
    "WrappedGaussianStimulus",
    "WrappedTopHatStimulus",
    "geometry_fields",
    "wrapped_gaussian_profile",
]
