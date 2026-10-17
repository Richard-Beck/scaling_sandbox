from __future__ import annotations

import numpy as np

from bulk_surface import (
    ConservationPool,
    DiffuseDomainSimulator,
    DiffuseGeometry,
    DiffuseInitialState,
    PhaseFieldMechanics,
    ReactionSystem,
    SimulationConfig,
)
from bulk_surface.diffuse_examples import passive_membrane_case


def _exchange_kernel(time, surface, bulk_local, stimulus, parameters,
                     surface_local_out, exchange_out):
    _ = time, stimulus
    surface_local_out.fill(0.0)
    exchange_out[0] = parameters[0] * bulk_local[0] - parameters[1] * surface[0]


def small_geometry() -> DiffuseGeometry:
    return DiffuseGeometry.circle(
        radius=1.0, dx=0.2, interface_width=0.3, box_half_width=2.2
    )


def test_diffuse_circle_measures_and_outward_normal_are_consistent():
    geometry = small_geometry()
    fields = geometry.fields()
    assert abs(fields.area / np.pi - 1.0) < 0.15
    assert abs(fields.perimeter / (2.0 * np.pi) - 1.0) < 0.12
    row = np.argmin(np.abs(geometry.y))
    right = np.argmin(np.abs(geometry.x - 1.0))
    assert fields.normal_x[row, right] > 0.9
    interface = fields.delta > 0.5 * fields.delta.max()
    assert np.median(fields.curvature[interface]) > 0


def test_fixed_diffuse_exchange_conserves_paired_pool():
    geometry = small_geometry()
    reactions = ReactionSystem(
        bulk_species=("bulk",),
        surface_species=("surface",),
        parameters=np.asarray([0.8, 0.3]),
        bulk_stoichiometry=np.asarray([[-1.0]]),
        surface_stoichiometry=np.asarray([[1.0]]),
        surface_kernel=_exchange_kernel,
    )
    pool = ConservationPool("combined", np.ones(1), np.ones(1))
    simulator = DiffuseDomainSimulator(
        geometry,
        reactions,
        bulk_diffusion=np.zeros(1),
        surface_diffusion=np.zeros(1),
        mechanics=PhaseFieldMechanics(
            surface_force_coefficients=np.zeros(1), frozen=True
        ),
        conservation_pools=(pool,),
    )
    initial = DiffuseInitialState(
        bulk=np.full((1, *geometry.shape), 1.2),
        surface=np.full((1, *geometry.shape), 0.1),
    )
    result = simulator.run(
        initial, SimulationConfig(end_time=0.04, dt=0.002, output_interval=0.01)
    )
    total = result.pool_totals["combined"]
    assert np.max(np.abs(total / total[0] - 1.0)) < 2e-12
    assert result.surface_material[-1, 0] > result.surface_material[0, 0]


def test_projected_surface_diffusion_conserves_membrane_material():
    geometry = small_geometry()
    fields = geometry.fields()
    angle = 2.0 * np.pi * fields.membrane_fraction.reshape(geometry.shape)
    initial_surface = 1.0 + 0.2 * np.cos(angle)
    case = passive_membrane_case(
        geometry, initial_surface, surface_diffusion=0.01, freeze_boundary=True
    )
    result = case.simulator.run(
        case.initial,
        SimulationConfig(end_time=0.002, dt=0.0001, output_interval=0.001),
    )
    material = result.surface_material[:, 0]
    assert np.max(np.abs(material / material[0] - 1.0)) < 2e-12
    interface = result.phi[-1] > 0.2
    assert np.isfinite(result.surface[-1, 0, interface]).all()
    delta = fields.delta
    initial_mode = np.sum(delta * result.surface[0, 0] * np.cos(angle))
    final_mode = np.sum(delta * result.surface[-1, 0] * np.cos(angle))
    assert 0 < final_mode < initial_mode


def test_localized_protrusion_moves_cell_and_hard_area_control_holds():
    geometry = small_geometry()
    fields = geometry.fields()
    angle = 2.0 * np.pi * fields.membrane_fraction.reshape(geometry.shape)
    initial_surface = 0.05 + 2.0 * np.exp(-0.5 * (np.angle(np.exp(1j * angle)) / 0.35) ** 2)
    case = passive_membrane_case(
        geometry, initial_surface, force_coefficient=2.0, freeze_boundary=False
    )
    result = case.simulator.run(
        case.initial,
        SimulationConfig(end_time=0.02, dt=0.001, output_interval=0.01),
    )
    assert result.centroid[-1, 0] > result.centroid[0, 0]
    assert abs(result.centroid[-1, 1] - result.centroid[0, 1]) < 2e-3
    assert np.max(np.abs(result.area / result.area[0] - 1.0)) < 2e-10
    material = result.surface_material[:, 0]
    assert np.max(np.abs(material / material[0] - 1.0)) < 2e-12


def test_kelvin_voigt_perimeter_tension_relaxes_an_elongated_boundary():
    geometry = small_geometry()
    xx, yy = np.meshgrid(geometry.x, geometry.y)
    semi_major, semi_minor = 1.25, 0.8
    radial = np.sqrt((xx / semi_major) ** 2 + (yy / semi_minor) ** 2)
    geometry.phi = 0.5 * (
        1.0 - np.tanh((radial - 1.0) / (np.sqrt(2.0) * geometry.epsilon))
    )
    geometry.target_area = float(np.sum(geometry.phi) * geometry.dx**2)
    initial_perimeter = geometry.fields().perimeter
    rest_perimeter = 2.0 * np.pi * np.sqrt(semi_major * semi_minor)
    reactions = ReactionSystem(
        bulk_species=(), surface_species=(),
        bulk_stoichiometry=np.empty((0, 0)),
        surface_stoichiometry=np.empty((0, 0)),
    )
    simulator = DiffuseDomainSimulator(
        geometry,
        reactions,
        bulk_diffusion=np.empty(0),
        surface_diffusion=np.empty(0),
        mechanics=PhaseFieldMechanics(
            surface_force_coefficients=np.empty(0),
            surface_tension=0.0,
            perimeter_elasticity=3.0,
            perimeter_viscosity=0.2,
            rest_perimeter=rest_perimeter,
            enforce_area=True,
        ),
    )
    empty = DiffuseInitialState(
        bulk=np.empty((0, *geometry.shape)),
        surface=np.empty((0, *geometry.shape)),
    )
    result = simulator.run(
        empty, SimulationConfig(end_time=0.02, dt=0.001, output_interval=0.01)
    )
    assert result.perimeter[-1] < initial_perimeter
    assert np.max(np.abs(result.area / result.area[0] - 1.0)) < 2e-10
