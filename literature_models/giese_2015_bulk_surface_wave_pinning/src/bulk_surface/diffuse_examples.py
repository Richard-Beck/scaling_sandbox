"""Small mechanochemical configurations for the diffuse-domain backend."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi
from pathlib import Path

import numpy as np
from numba import njit
from scipy.interpolate import RegularGridInterpolator

from .core import ConservationPool, ReactionSystem
from .diffuse import (
    DiffuseDomainSimulator,
    DiffuseGeometry,
    DiffuseInitialState,
    PhaseFieldMechanics,
)
from .protocols import NoStimulus, WrappedGaussianStimulus, wrapped_gaussian_profile
from .shared_conditions import SharedChemistry, shared_local_response_kernel
from .moving_curve import (
    MovingCurveInitialState,
    MovingCurveSimulator,
    MovingMembrane,
)


@njit(cache=True)
def _diffuse_wave_pinning_kernel(
    time, surface, bulk_local, stimulus, parameters,
    surface_local_out, exchange_out,
):
    _ = time
    surface_local_out.fill(0.0)
    basal, feedback, half_level, inactivation = parameters
    half_squared = half_level * half_level
    for node in range(surface.shape[1]):
        active = max(surface[0, node], 0.0)
        inactive = max(bulk_local[0, node], 0.0)
        activation = basal + feedback * active * active / (
            half_squared + active * active
        ) + stimulus[node]
        exchange_out[0, node] = activation * inactive - inactivation * active


@dataclass(frozen=True)
class DiffuseExampleCase:
    simulator: DiffuseDomainSimulator
    initial: DiffuseInitialState


@dataclass(frozen=True)
class MovingCurveExampleCase:
    simulator: MovingCurveSimulator
    initial: MovingCurveInitialState


def diffuse_wave_pinning_case(
    geometry: DiffuseGeometry | None = None,
    *,
    protrusion: float = 0.35,
    stimulus: WrappedGaussianStimulus | None = None,
    freeze_boundary: bool = False,
) -> DiffuseExampleCase:
    """Membrane-active/bulk-inactive wave pinning driving protrusion."""
    g = geometry or DiffuseGeometry.circle(
        radius=2.0, dx=0.2, interface_width=0.3, box_half_width=4.5
    )
    cue = stimulus or WrappedGaussianStimulus(
        amplitude=1.5,
        center_fraction=0.0,
        width_fraction=0.04,
        start=0.25,
        stop=1.0,
    )
    reactions = ReactionSystem(
        bulk_species=("inactive",),
        surface_species=("active",),
        parameters=np.asarray([0.067, 1.0, 1.0, 1.0]),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=_diffuse_wave_pinning_kernel,
        stimulus=cue,
    )
    pool = ConservationPool("gtpase", np.asarray([1.0]), np.asarray([1.0]))
    mechanics = PhaseFieldMechanics(
        surface_force_coefficients=np.asarray([protrusion]),
        relaxation_time=1.0,
        surface_tension=1.0,
        area_penalty=0.2,
        remove_mean_force=True,
        enforce_area=True,
        frozen=freeze_boundary,
    )
    simulator = DiffuseDomainSimulator(
        g,
        reactions,
        bulk_diffusion=np.asarray([1.0]),
        surface_diffusion=np.asarray([0.02]),
        mechanics=mechanics,
        conservation_pools=(pool,),
    )
    fields = g.fields()
    localized = wrapped_gaussian_profile(
        fields.membrane_fraction.reshape(g.shape), 0.0, 0.04
    )
    initial = DiffuseInitialState(
        bulk=np.full((1, *g.shape), 1.5),
        surface=(0.18 + 0.002 * localized)[None, :, :],
    )
    return DiffuseExampleCase(simulator=simulator, initial=initial)


def passive_membrane_case(
    geometry: DiffuseGeometry,
    initial_surface: np.ndarray,
    *,
    force_coefficient: float = 0.0,
    surface_diffusion: float = 0.0,
    freeze_boundary: bool = False,
) -> DiffuseExampleCase:
    """A reaction-free membrane tracer useful for mechanics validation."""
    reactions = ReactionSystem(
        bulk_species=(),
        surface_species=("tracer",),
        bulk_stoichiometry=np.empty((0, 0)),
        surface_stoichiometry=np.empty((1, 0)),
        stimulus=NoStimulus(),
    )
    mechanics = PhaseFieldMechanics(
        surface_force_coefficients=np.asarray([force_coefficient]),
        remove_mean_force=True,
        enforce_area=True,
        frozen=freeze_boundary,
    )
    simulator = DiffuseDomainSimulator(
        geometry,
        reactions,
        bulk_diffusion=np.empty(0),
        surface_diffusion=np.asarray([surface_diffusion]),
        mechanics=mechanics,
    )
    surface = np.asarray(initial_surface, dtype=float)
    if surface.shape != geometry.shape:
        raise ValueError("initial_surface must match the diffuse grid")
    initial = DiffuseInitialState(
        bulk=np.empty((0, *geometry.shape)),
        surface=surface[None, :, :],
    )
    return DiffuseExampleCase(simulator=simulator, initial=initial)


def released_mori_case(
    source: str | Path,
    geometry: DiffuseGeometry,
    *,
    source_time: float | None = None,
    protrusion: float = 0.35,
    surface_tension: float = 1.0,
    area_penalty: float = 0.2,
    perimeter_elasticity: float = 0.0,
    perimeter_viscosity: float = 0.0,
) -> tuple[DiffuseExampleCase, dict[str, float | str]]:
    """Release a polarized fixed-circle Mori state into the moving backend.

    The bulk field is bilinearly interpolated in laboratory coordinates and
    the ordered membrane field is periodically interpolated by polar angle.
    Bulk and membrane amounts are then matched separately to the source run,
    so switching numerical representations does not create or destroy protein.
    The original transient cue is absent after release.
    """
    source = Path(source)
    trajectory = source / "trajectory.npz" if source.is_dir() else source
    with np.load(trajectory, allow_pickle=False) as raw:
        times = raw["times"]
        index = times.size - 1 if source_time is None else int(
            np.argmin(np.abs(times - source_time))
        )
        selected_time = float(times[index])
        source_x = raw["x"].copy()
        source_y = raw["y"].copy()
        source_phi = raw["phi"].copy()
        source_bulk = raw["bulk"][index, 0].reshape(source_phi.shape).copy()
        source_surface = raw["surface"][index, 0].copy()
        source_fraction = raw["membrane_fraction"].copy()
        source_points = raw["membrane_points"].copy()

    source_dx = float(source_x[1] - source_x[0])
    # The fixed geometry applies a constant diffuse-area quadrature correction.
    radius = float(np.mean(np.linalg.norm(source_points, axis=1)))
    source_area = pi * radius * radius
    bulk_scale = source_area / (float(np.sum(source_phi)) * source_dx**2)
    source_bulk_amount = float(
        np.sum(source_bulk * source_phi) * source_dx**2 * bulk_scale
    )
    edges = np.linalg.norm(np.roll(source_points, -1, axis=0) - source_points,
                           axis=1)
    source_surface_weights = 0.5 * (edges + np.roll(edges, 1))
    source_surface_amount = float(source_surface @ source_surface_weights)

    xx, yy = np.meshgrid(geometry.x, geometry.y)
    interpolate_bulk = RegularGridInterpolator(
        (source_y, source_x), source_bulk,
        bounds_error=False, fill_value=float(np.mean(source_bulk)),
    )
    mapped_bulk = interpolate_bulk(np.column_stack((yy.ravel(), xx.ravel()))).reshape(
        geometry.shape
    )
    fields = geometry.fields()
    target_fraction = fields.membrane_fraction.reshape(geometry.shape)
    periodic_fraction = np.concatenate((source_fraction, [source_fraction[0] + 1.0]))
    periodic_surface = np.concatenate((source_surface, [source_surface[0]]))
    mapped_surface = np.interp(target_fraction, periodic_fraction, periodic_surface)

    current_bulk_amount = float(
        np.sum(geometry.phi * mapped_bulk) * geometry.dx**2
    )
    current_surface_amount = float(
        np.sum(fields.delta * mapped_surface) * geometry.dx**2
    )
    mapped_bulk *= source_bulk_amount / current_bulk_amount
    mapped_surface *= source_surface_amount / current_surface_amount

    chemistry = SharedChemistry()
    reactions = ReactionSystem(
        bulk_species=("inactive_bulk",),
        surface_species=("active_surface",),
        parameters=chemistry.kernel_parameters(),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=shared_local_response_kernel,
        stimulus=NoStimulus(),
    )
    pool = ConservationPool("total_protein", np.ones(1), np.ones(1))
    mechanics = PhaseFieldMechanics(
        surface_force_coefficients=np.asarray([protrusion]),
        relaxation_time=1.0,
        surface_tension=surface_tension,
        area_penalty=area_penalty,
        remove_mean_force=True,
        enforce_area=True,
        perimeter_elasticity=perimeter_elasticity,
        perimeter_viscosity=perimeter_viscosity,
    )
    simulator = DiffuseDomainSimulator(
        geometry,
        reactions,
        bulk_diffusion=np.asarray([chemistry.bulk_diffusion]),
        surface_diffusion=np.asarray([chemistry.surface_diffusion]),
        mechanics=mechanics,
        conservation_pools=(pool,),
    )
    case = DiffuseExampleCase(
        simulator=simulator,
        initial=DiffuseInitialState(
            bulk=mapped_bulk[None, :, :],
            surface=mapped_surface[None, :, :],
        ),
    )
    provenance: dict[str, float | str] = {
        "source_trajectory": str(trajectory.resolve()),
        "source_time": selected_time,
        "source_bulk_amount": source_bulk_amount,
        "source_surface_amount": source_surface_amount,
        "source_total_amount": source_bulk_amount + source_surface_amount,
        "mapped_bulk_amount": float(
            np.sum(geometry.phi * mapped_bulk) * geometry.dx**2
        ),
        "mapped_surface_amount": float(
            np.sum(fields.delta * mapped_surface) * geometry.dx**2
        ),
    }
    return case, provenance


def released_mori_curve_case(
    source: str | Path,
    geometry: DiffuseGeometry,
    *,
    source_time: float | None = None,
    protrusion: float = 0.35,
    surface_tension: float = 1.0,
    area_penalty: float = 0.2,
    perimeter_elasticity: float = 0.0,
    perimeter_viscosity: float = 0.0,
) -> tuple[MovingCurveExampleCase, dict[str, float | str]]:
    """Release fixed Mori fields using an explicitly moving 1D membrane."""
    source = Path(source)
    trajectory = source / "trajectory.npz" if source.is_dir() else source
    with np.load(trajectory, allow_pickle=False) as raw:
        times = raw["times"]
        index = times.size - 1 if source_time is None else int(
            np.argmin(np.abs(times - source_time))
        )
        selected_time = float(times[index])
        source_x = raw["x"].copy()
        source_y = raw["y"].copy()
        source_phi = raw["phi"].copy()
        source_bulk = raw["bulk"][index, 0].reshape(source_phi.shape).copy()
        source_surface = raw["surface"][index, 0].copy()
        source_points = raw["membrane_points"].copy()

    source_dx = float(source_x[1] - source_x[0])
    radius = float(np.mean(np.linalg.norm(source_points, axis=1)))
    source_area = pi * radius * radius
    bulk_scale = source_area / (float(np.sum(source_phi)) * source_dx**2)
    source_bulk_amount = float(
        np.sum(source_bulk * source_phi) * source_dx**2 * bulk_scale
    )
    membrane = MovingMembrane(source_points)
    source_surface_amount = float(source_surface @ membrane.weights)

    xx, yy = np.meshgrid(geometry.x, geometry.y)
    interpolate_bulk = RegularGridInterpolator(
        (source_y, source_x), source_bulk,
        bounds_error=False, fill_value=float(np.mean(source_bulk)),
    )
    mapped_bulk = interpolate_bulk(
        np.column_stack((yy.ravel(), xx.ravel()))
    ).reshape(geometry.shape)
    mapped_bulk *= source_bulk_amount / float(
        np.sum(geometry.phi * mapped_bulk) * geometry.dx**2
    )

    chemistry = SharedChemistry()
    reactions = ReactionSystem(
        bulk_species=("inactive_bulk",),
        surface_species=("active_surface",),
        parameters=chemistry.kernel_parameters(),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=shared_local_response_kernel,
        stimulus=NoStimulus(),
    )
    pool = ConservationPool("total_protein", np.ones(1), np.ones(1))
    mechanics = PhaseFieldMechanics(
        surface_force_coefficients=np.asarray([protrusion]),
        relaxation_time=1.0,
        surface_tension=surface_tension,
        area_penalty=area_penalty,
        remove_mean_force=True,
        enforce_area=True,
        perimeter_elasticity=perimeter_elasticity,
        perimeter_viscosity=perimeter_viscosity,
    )
    simulator = MovingCurveSimulator(
        geometry, membrane, reactions,
        bulk_diffusion=np.asarray([chemistry.bulk_diffusion]),
        surface_diffusion=np.asarray([chemistry.surface_diffusion]),
        mechanics=mechanics,
        conservation_pools=(pool,),
    )
    case = MovingCurveExampleCase(
        simulator=simulator,
        initial=MovingCurveInitialState(
            bulk=mapped_bulk[None, :, :],
            surface=source_surface[None, :],
        ),
    )
    provenance: dict[str, float | str] = {
        "source_trajectory": str(trajectory.resolve()),
        "source_time": selected_time,
        "source_bulk_amount": source_bulk_amount,
        "source_surface_amount": source_surface_amount,
        "source_total_amount": source_bulk_amount + source_surface_amount,
        "mapped_bulk_amount": float(
            np.sum(geometry.phi * mapped_bulk) * geometry.dx**2
        ),
        "mapped_surface_amount": float(source_surface @ membrane.weights),
        "membrane_representation": "ordered_material_1d_curve",
    }
    return case, provenance
