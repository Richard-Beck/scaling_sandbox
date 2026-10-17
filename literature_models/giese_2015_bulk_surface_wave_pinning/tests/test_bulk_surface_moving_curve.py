from __future__ import annotations

import numpy as np

from bulk_surface import (
    ConservationPool,
    DiffuseGeometry,
    MovingCurveInitialState,
    MovingCurveSimulator,
    MovingMembrane,
    PhaseFieldMechanics,
    ReactionSystem,
    SimulationConfig,
)
from bulk_surface.geometry import _bilinear_trace


def _geometry_and_membrane():
    geometry = DiffuseGeometry.circle(
        radius=1.0, dx=0.2, interface_width=0.2, box_half_width=2.2
    )
    return geometry, MovingMembrane.circle(1.0, 64)


def _exchange_kernel(time, surface, bulk_trace, stimulus, parameters,
                     surface_local_out, exchange_out):
    _ = time, stimulus
    surface_local_out.fill(0.0)
    exchange_out[0] = parameters[0] * bulk_trace[0] - parameters[1] * surface[0]


def test_explicit_curve_exchange_conserves_bulk_surface_pool():
    geometry, membrane = _geometry_and_membrane()
    reactions = ReactionSystem(
        bulk_species=("bulk",), surface_species=("surface",),
        parameters=np.asarray([0.8, 0.3]),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=_exchange_kernel,
    )
    simulator = MovingCurveSimulator(
        geometry, membrane, reactions,
        bulk_diffusion=np.zeros(1), surface_diffusion=np.zeros(1),
        mechanics=PhaseFieldMechanics(
            surface_force_coefficients=np.zeros(1), frozen=True
        ),
        conservation_pools=(
            ConservationPool("combined", np.ones(1), np.ones(1)),
        ),
    )
    initial = MovingCurveInitialState(
        bulk=np.full((1, *geometry.shape), 1.2),
        surface=np.full((1, membrane.n_points), 0.1),
    )
    result = simulator.run(
        initial, SimulationConfig(end_time=0.04, dt=0.002, output_interval=0.01)
    )
    total = result.pool_totals["combined"]
    assert np.max(np.abs(total / total[0] - 1.0)) < 2e-12
    assert result.surface_material[-1, 0] > result.surface_material[0, 0]


def test_one_dimensional_membrane_diffusion_smooths_and_conserves():
    geometry, membrane = _geometry_and_membrane()
    reactions = ReactionSystem(
        bulk_species=(), surface_species=("tracer",),
        bulk_stoichiometry=np.empty((0, 0)),
        surface_stoichiometry=np.empty((1, 0)),
    )
    simulator = MovingCurveSimulator(
        geometry, membrane, reactions,
        bulk_diffusion=np.empty(0), surface_diffusion=np.asarray([0.05]),
        mechanics=PhaseFieldMechanics(
            surface_force_coefficients=np.zeros(1), frozen=True
        ),
    )
    theta = 2.0 * np.pi * np.arange(membrane.n_points) / membrane.n_points
    initial = MovingCurveInitialState(
        bulk=np.empty((0, *geometry.shape)),
        surface=(1.0 + 0.3 * np.cos(theta))[None, :],
    )
    result = simulator.run(
        initial, SimulationConfig(end_time=0.02, dt=0.001, output_interval=0.01)
    )
    mass = result.surface_material[:, 0]
    assert np.max(np.abs(mass / mass[0] - 1.0)) < 2e-12
    assert np.ptp(result.surface[-1, 0]) < np.ptp(result.surface[0, 0])


def test_material_membrane_tracks_phi_half_contour_during_migration():
    geometry, membrane = _geometry_and_membrane()
    reactions = ReactionSystem(
        bulk_species=(), surface_species=("active",),
        bulk_stoichiometry=np.empty((0, 0)),
        surface_stoichiometry=np.empty((1, 0)),
    )
    simulator = MovingCurveSimulator(
        geometry, membrane, reactions,
        bulk_diffusion=np.empty(0), surface_diffusion=np.zeros(1),
        mechanics=PhaseFieldMechanics(
            surface_force_coefficients=np.asarray([2.0]),
            enforce_area=True,
        ),
    )
    theta = 2.0 * np.pi * np.arange(membrane.n_points) / membrane.n_points
    active = 0.05 + 2.0 * np.exp(-0.5 * (np.angle(np.exp(1j * theta)) / 0.35) ** 2)
    initial = MovingCurveInitialState(
        bulk=np.empty((0, *geometry.shape)), surface=active[None, :]
    )
    result = simulator.run(
        initial, SimulationConfig(end_time=0.02, dt=0.001, output_interval=0.01)
    )
    final_trace = _bilinear_trace(
        result.x, result.y, result.membrane_points[-1]
    ) @ result.phi[-1].ravel()
    assert np.max(np.abs(final_trace - 0.5)) < 2e-5
    assert result.centroid[-1, 0] > result.centroid[0, 0]
    assert np.max(np.abs(result.area / result.area[0] - 1.0)) < 2e-10
    mass = result.surface_material[:, 0]
    assert np.max(np.abs(mass / mass[0] - 1.0)) < 2e-12
    assert np.max(result.maximum_edge_ratio) < 1.1

