from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from migration1d import (  # noqa: E402
    FieldForce,
    MechanicsParameters,
    MigrationModel,
    SimulationConfig,
    SumForces,
    run_simulation,
    solve_endpoint_mechanics,
)
from migration1d.model import (  # noqa: E402
    control_volumes,
    diffusion_amount_rate,
    reaction_diffusion_jacobian_sparsity,
)


class ThreeFieldChemistry:
    species = ("a", "b", "inhibitor")
    membrane_diffusion = np.asarray([0.0, 0.0, 0.0])
    cytoplasm_diffusion = np.asarray([0.0, 0.0, 0.0])
    initial_membrane = np.asarray([0.2, 0.4, 0.6])
    initial_cytoplasm = np.asarray([1.0, 1.2, 1.4])

    def reaction(self, time, position, membrane, cytoplasm, strain, stress):
        return np.zeros_like(membrane), np.zeros_like(cytoplasm)


class ProfileChemistry(ThreeFieldChemistry):
    def initial_profiles(self, position):
        fraction = (position-position[0])/(position[-1]-position[0])
        membrane = np.repeat(self.initial_membrane[:, None], len(position), axis=1)
        cytoplasm = np.repeat(self.initial_cytoplasm[:, None], len(position), axis=1)
        membrane[0] += 0.1*(2.0*fraction-1.0)
        cytoplasm[0] -= 0.1*(2.0*fraction-1.0)
        return membrane, cytoplasm


def arrays(nodes=11, species=1):
    x = np.linspace(0.0, 10.0, nodes)
    membrane = np.ones((species, nodes))
    cytoplasm = np.ones_like(membrane)
    return x, membrane, cytoplasm


def test_no_flux_diffusion_conserves_each_field_amount():
    x = np.asarray([0.0, 0.3, 1.1, 2.0])
    concentration = np.asarray([[1.0, 4.0, 0.5, 2.0], [2.0, 1.0, 3.0, 0.0]])
    rate = diffusion_amount_rate(concentration, x, np.asarray([0.2, 1.3]))
    assert np.allclose(rate.sum(axis=1), 0.0)


def test_field_force_samples_named_field_at_each_endpoint():
    chemistry = ThreeFieldChemistry()
    mapping = FieldForce(
        field="b", compartment="membrane", gain=3.0,
        response=lambda value: value**2, direction="outward",
    )
    model = MigrationModel(chemistry=chemistry, left_force=mapping, right_force=mapping)
    x, membrane, cytoplasm = arrays(species=3)
    membrane[1] = np.linspace(2.0, 4.0, len(x))
    _, applied, _, _ = solve_endpoint_mechanics(model, 0.0, x, membrane, cytoplasm)
    assert np.allclose(applied, [-12.0, 48.0])


def test_force_mappings_can_be_composed_and_can_inspect_full_state():
    chemistry = ThreeFieldChemistry()
    mapped = FieldForce("a", gain=2.0, direction="lab")

    def difference(state):
        return state.sample("inhibitor", "cytoplasm") - state.sample("b")

    model = MigrationModel(
        chemistry=chemistry,
        left_force=SumForces((mapped, difference)),
        right_force=SumForces((mapped, difference)),
    )
    x, membrane, cytoplasm = arrays(species=3)
    membrane[0] = 3.0
    membrane[1] = 0.5
    cytoplasm[2] = 2.5
    _, applied, _, _ = solve_endpoint_mechanics(model, 0.0, x, membrane, cytoplasm)
    assert np.allclose(applied, 8.0)


def test_force_free_relaxed_element_is_stationary():
    model = MigrationModel()
    x, membrane, cytoplasm = arrays()
    velocity, applied, spring, dashpot = solve_endpoint_mechanics(
        model, 0.0, x, membrane, cytoplasm
    )
    assert np.allclose(velocity, 0.0)
    assert np.allclose(applied, 0.0)
    assert spring == 0.0
    assert dashpot == 0.0


def test_equal_lab_frame_forces_translate_without_deforming():
    model = MigrationModel(
        mechanics=MechanicsParameters(left_friction=2.0, right_friction=2.0),
        left_force=lambda state: 3.0,
        right_force=lambda state: 3.0,
    )
    x, membrane, cytoplasm = arrays()
    velocity, _, spring, dashpot = solve_endpoint_mechanics(
        model, 0.0, x, membrane, cytoplasm
    )
    assert np.allclose(velocity, 1.5)
    assert spring == 0.0
    assert abs(dashpot) < 1e-14


def test_kelvin_voigt_dashpot_resists_length_rate():
    mechanics = MechanicsParameters(
        spring_stiffness=2.0, dashpot_viscosity=3.0,
        left_friction=1.0, right_friction=1.0,
    )
    model = MigrationModel(
        mechanics=mechanics,
        left_force=lambda state: -2.0,
        right_force=lambda state: 2.0,
    )
    x, membrane, cytoplasm = arrays()
    velocity, _, spring, dashpot = solve_endpoint_mechanics(
        model, 0.0, x, membrane, cytoplasm
    )
    assert velocity[0] < 0.0 < velocity[-1]
    assert spring == 0.0
    assert np.isclose(dashpot, 12.0 / 7.0)
    # The dashpot force plus endpoint drag balances the imposed force.
    assert np.isclose(-velocity[0] + dashpot, 2.0)


def test_spring_restoring_force_opposes_extension():
    model = MigrationModel(mechanics=MechanicsParameters(spring_stiffness=4.0))
    x, membrane, cytoplasm = arrays()
    x = np.linspace(0.0, 11.0, len(x))
    velocity, _, spring, _ = solve_endpoint_mechanics(
        model, 0.0, x, membrane, cytoplasm
    )
    assert spring == 4.0
    assert velocity[0] > 0.0 > velocity[-1]


def test_arbitrary_field_count_and_moving_domain_conserve_material():
    chemistry = ThreeFieldChemistry()
    model = MigrationModel(
        chemistry=chemistry,
        left_force=FieldForce("b", gain=0.2, direction="lab"),
        right_force=FieldForce("b", gain=0.2, direction="lab"),
    )
    result = run_simulation(
        model,
        SimulationConfig(
            nodes=15, duration=0.5, output_interval=0.1,
            max_step=0.02, noise_amplitude=0.0,
        ),
    )
    membrane_total, cytoplasm_total = result.compartment_totals()
    assert result.membrane.shape[1] == 3
    assert np.max(np.abs(membrane_total - membrane_total[0])) < 1e-10
    assert np.max(np.abs(cytoplasm_total - cytoplasm_total[0])) < 1e-10
    assert result.center[-1] > result.center[0]
    assert np.allclose(result.length, result.length[0], atol=1e-9)


def test_moving_control_volume_representation_dilutes_uniform_field():
    x0 = np.linspace(0.0, 10.0, 11)
    x1 = np.linspace(0.0, 20.0, 11)
    amount = 3.0 * control_volumes(x0)
    concentration = amount / control_volumes(x1)
    assert np.allclose(concentration, 1.5)
    assert np.isclose(np.sum(amount), 30.0)


def test_chemistry_can_supply_mass_balanced_spatial_initial_profiles():
    model = MigrationModel(chemistry=ProfileChemistry())
    result = run_simulation(
        model,
        SimulationConfig(nodes=11, duration=0.0, noise_amplitude=0.0),
    )
    assert result.membrane[0, 0, -1] > result.membrane[0, 0, 0]
    np.testing.assert_allclose(
        result.membrane[0, 0] + result.cytoplasm[0, 0], 1.2
    )


def test_reaction_diffusion_jacobian_pattern_has_expected_shape_and_diagonal():
    pattern = reaction_diffusion_jacobian_sparsity(species=3, nodes=11)
    assert pattern.shape == (68, 68)
    assert np.all(pattern.diagonal())
    assert pattern.nnz < 68 * 68
