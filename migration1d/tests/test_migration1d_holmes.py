from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from migration1d import (  # noqa: E402
    HolmesChemistry,
    LinearCue,
    MechanicsParameters,
    MigrationModel,
    RacRhoForce,
    SimulationConfig,
    run_simulation,
    solve_endpoint_mechanics,
)


def test_holmes_reaction_matches_reduced_reference_formula():
    chemistry = HolmesChemistry()
    position = np.linspace(0.0, 10.0, 5)
    active = np.asarray((0.6, 1.1, 0.8))[:, None] * np.ones((3, 5))
    inactive = np.asarray((1.8, 6.4, 2.3))[:, None] * np.ones((3, 5))
    rate, opposite = chemistry.reaction(
        0.0, position, active, inactive, np.zeros(5), np.zeros(5)
    )
    p = chemistry.parameters
    expected = np.asarray((
        p.cdc42_activation / (1.0 + (0.8 / p.rho_inhibits_cdc42_half_level) ** 3)
        * 1.8 / p.cdc42_total - 0.6,
        (p.rac_basal_activation + p.cdc42_to_rac_feedback * 0.6
         + p.signal_baseline) * 6.4 / p.rac_total - 1.1,
        p.rho_activation / (1.0 + (1.1 / p.rac_inhibits_rho_half_level) ** 3)
        * 2.3 / p.rho_total - 0.8,
    ))
    np.testing.assert_allclose(rate, np.repeat(expected[:, None], 5, axis=1))
    np.testing.assert_array_equal(opposite, -rate)


def test_initial_state_is_the_prestimulus_reference_equilibrium():
    parameters = replace(HolmesChemistry().parameters, signal_baseline=0.0)
    chemistry = HolmesChemistry(parameters=parameters)
    position = np.linspace(0.0, 10.0, 7)
    active = np.repeat(chemistry.initial_membrane[:, None], 7, axis=1)
    inactive = np.repeat(chemistry.initial_cytoplasm[:, None], 7, axis=1)
    rate, _ = chemistry.reaction(
        0.0, position, active, inactive, np.zeros(7), np.zeros(7)
    )
    assert np.max(np.abs(rate)) < 2e-8


def test_front_rac_and_rear_rho_cooperate_to_translate_without_deforming():
    chemistry = HolmesChemistry()
    force = RacRhoForce(rac_gain=1.0, rho_gain=1.0)
    model = MigrationModel(
        mechanics=MechanicsParameters(
            spring_stiffness=2.0, dashpot_viscosity=3.0,
            left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry,
        left_force=force,
        right_force=force,
    )
    position = np.linspace(0.0, 10.0, 11)
    active = np.zeros((3, 11))
    inactive = np.ones((3, 11))
    active[1, -1] = 2.0  # Rac at the right/front: force to the right.
    active[2, 0] = 2.0   # Rho at the left/rear: force to the right.
    velocity, applied, _, dashpot = solve_endpoint_mechanics(
        model, 0.0, position, active, inactive
    )
    assert np.all(applied > 0.0)
    np.testing.assert_allclose(velocity, velocity[0])
    assert velocity[0] > 0.0
    assert abs(dashpot) < 1e-14


def test_rac_rho_force_reverses_when_polarity_reverses():
    chemistry = HolmesChemistry()
    force = RacRhoForce(rac_gain=1.0, rho_gain=1.0)
    model = MigrationModel(chemistry=chemistry, left_force=force, right_force=force)
    position = np.linspace(0.0, 10.0, 11)
    inactive = np.ones((3, 11))
    rightward = np.zeros((3, 11))
    rightward[1, -1] = rightward[2, 0] = 2.0
    leftward = rightward[:, ::-1].copy()
    v_right, _, _, _ = solve_endpoint_mechanics(
        model, 0.0, position, rightward, inactive
    )
    v_left, _, _, _ = solve_endpoint_mechanics(
        model, 0.0, position, leftward, inactive
    )
    assert np.mean(v_right) > 0.0 > np.mean(v_left)
    np.testing.assert_allclose(v_left, -v_right)


def test_six_field_moving_run_conserves_each_gtpase_pool():
    chemistry = HolmesChemistry(cue=LinearCue(slope=0.04, origin=5.0))
    force = RacRhoForce(rac_gain=0.08, rho_gain=0.08)
    model = MigrationModel(
        mechanics=MechanicsParameters(
            spring_stiffness=0.1, dashpot_viscosity=2.0,
            left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry,
        left_force=force,
        right_force=force,
    )
    result = run_simulation(
        model,
        SimulationConfig(
            nodes=21, duration=2.0, output_interval=0.2,
            max_step=0.05, noise_amplitude=0.0,
        ),
    )
    active_total, inactive_total = result.compartment_totals()
    total = active_total + inactive_total
    relative_drift = np.abs(total / total[0] - 1.0)
    assert result.membrane.shape[1:] == (3, 21)
    assert result.cytoplasm.shape[1:] == (3, 21)
    assert np.max(relative_drift) < 2e-8
    assert result.center[-1] > result.center[0]
