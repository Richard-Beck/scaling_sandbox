from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from maree_migration.model import reaction_rates  # noqa: E402
from migration1d import (  # noqa: E402
    EroumeChemistry,
    EroumeTwoEdgeCue,
    MechanicsParameters,
    MigrationModel,
    RacRhoForce,
    SimulationConfig,
    run_simulation,
)


def contrast(values):
    return ((values.max(axis=1)-values.min(axis=1))
            / (values.max(axis=1)+values.min(axis=1)))


def test_eroume_reaction_and_native_initial_values_are_exact():
    chemistry = EroumeChemistry()
    position = np.linspace(0.0, 10.0, 5)
    fields = np.repeat(chemistry.initial_membrane[:, None], 5, axis=1)
    rate, companion = chemistry.reaction(
        0.0, position, fields, np.zeros_like(fields), np.zeros(5), np.zeros(5)
    )
    np.testing.assert_allclose(
        rate, reaction_rates(fields, chemistry.stimulus(0.0, position))
    )
    np.testing.assert_array_equal(companion, 0.0)
    np.testing.assert_array_equal(
        chemistry.initial_membrane, np.asarray(chemistry.parameters.initial)
    )


def test_eroume_coordinate_adapter_preserves_diffusive_timescale():
    native = EroumeChemistry(reference_length=70.54013)
    normalized = EroumeChemistry(reference_length=10.0)
    np.testing.assert_allclose(
        normalized.membrane_diffusion,
        native.membrane_diffusion*(10.0/70.54013)**2,
    )
    x_native = np.linspace(0.0, 70.54013, 11)
    x_normalized = np.linspace(0.0, 10.0, 11)
    np.testing.assert_allclose(
        native.stimulus(0.0, x_native), normalized.stimulus(0.0, x_normalized)
    )
    np.testing.assert_allclose(normalized.stimulus(10.1, x_normalized), 2.95)


def test_eroume_two_edge_cue_has_finite_endpoint_support_and_turns_off():
    cue = EroumeTwoEdgeCue(
        left_level=5.1, right_level=4.65, basal_level=2.95,
        width_fraction=0.06, start=0.0, stop=10.0,
    )
    position = np.linspace(0.0, 10.0, 101)
    values = cue.values(5.0, position)
    assert values[0] == 5.1
    assert values[-1] == 4.65
    assert np.all(values[7:-7] == 2.95)
    assert np.all(cue.values(10.0, position) == 2.95)


def test_native_protocol_retains_polarity_after_cue_removal_while_moving():
    chemistry = EroumeChemistry(reference_length=10.0)
    force = RacRhoForce(
        rac_field="Rac_act", rho_field="Rho_act",
        rac_gain=0.08, rho_gain=0.08,
        rac_half_level=3.0, rho_half_level=1.25,
    )
    model = MigrationModel(
        mechanics=MechanicsParameters(
            reference_length=10.0, spring_stiffness=1.0,
            dashpot_viscosity=2.0, left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry, left_force=force, right_force=force,
    )
    result = run_simulation(
        model,
        SimulationConfig(
            nodes=31, duration=30.0, output_interval=2.0, max_step=0.2,
            noise_amplitude=0.0, mechanics_activation_time=10.0,
        ),
    )
    assert contrast(result.membrane[:, 3])[-1] > 0.8
    assert result.center[-1] > result.center[0]
    membrane_total, _ = result.compartment_totals()
    for pair in ((0, 1), (2, 3), (4, 5)):
        total = membrane_total[:, list(pair)].sum(axis=1)
        assert np.max(np.abs(total/total[0]-1.0)) < 1e-8
