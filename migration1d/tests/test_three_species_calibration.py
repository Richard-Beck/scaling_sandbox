from dataclasses import replace

import numpy as np

from migration1d.holmes import HolmesChemistry
from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from migration1d.three_species.model import (
    GaussianRacPulse,
    active_reaction,
    canonical_holmes_parameters,
    homogeneous_equilibrium,
    homogeneous_stability,
)
from migration1d.three_species.observables import (
    classify_three_species_trajectory,
    profile_metrics,
)
from migration1d.three_species.wgd import all_tasks, task_key, transition_group
from migration1d.three_species.wgd_extended import all_tasks as all_extended_tasks


def test_canonical_mapping_uses_documented_holmes_controls():
    parameters = canonical_holmes_parameters()
    assert np.isclose(parameters.rac_basal_activation, 0.2)
    assert np.isclose(parameters.cdc42_to_rac_feedback, 0.65)
    assert parameters.signal_baseline == 0.0
    assert parameters.active_diffusivity == 0.1
    assert parameters.inactive_diffusivity == 50.0


def test_relaxed_homogeneous_state_is_an_equilibrium():
    parameters = canonical_holmes_parameters()
    active = homogeneous_equilibrium(parameters)
    residual = active_reaction(active, parameters.totals - active, parameters)
    assert np.max(np.abs(residual)) < 1.0e-7
    chemistry = HolmesChemistry(
        replace(parameters, initial_active=tuple(float(value) for value in active))
    )
    np.testing.assert_allclose(chemistry.initial_membrane + chemistry.initial_cytoplasm,
                               parameters.totals)


def test_homogeneous_stability_reports_finite_growth_rates_and_mode():
    parameters = canonical_holmes_parameters()
    active = homogeneous_equilibrium(parameters)
    result = homogeneous_stability(parameters, active, domain_length=20.0)
    assert np.isfinite(result["homogeneous_growth_rate"])
    assert np.isfinite(result["maximum_spatial_growth_rate"])
    assert 1 <= result["fastest_spatial_mode"] <= 64


def test_gaussian_rac_pulse_is_local_and_transient():
    pulse = GaussianRacPulse(start=2.0, stop=4.0, center_fraction=0.1)
    position = np.linspace(0.0, 20.0, 101)
    assert np.all(pulse(1.0, position) == 0.0)
    active = pulse(3.0, position)
    assert np.argmax(active) == 10
    assert np.all(pulse(4.0, position) == 0.0)


def test_composite_metrics_recognize_cdc42_rac_front_and_rho_rear():
    coordinate = np.linspace(0.0, 1.0, 101)
    front = 0.1 + 1.5 * np.exp(-0.5 * ((coordinate - 0.15) / 0.08) ** 2)
    rear = front[::-1]
    profiles = np.stack((front, 1.5 * front, rear))
    metrics = profile_metrics(profiles)
    assert metrics["coherent_polarity"]
    assert metrics["cdc42_rac_correlation"] > 0.99
    assert metrics["cdc42_rho_correlation"] < 0.0
    assert metrics["rac_pole_count"] == 1


def test_stable_coherent_endpoint_is_finitely_polarizable():
    coordinate = np.linspace(0.0, 1.0, 101)
    front = 0.1 + 1.5 * np.exp(-0.5 * ((coordinate - 0.15) / 0.08) ** 2)
    endpoint = np.stack((front, 1.5 * front, front[::-1]))
    times = np.linspace(0.0, 100.0, 101)
    trajectory = np.repeat(endpoint[None, :, :], len(times), axis=0)
    result = classify_three_species_trajectory(
        times, trajectory, np.asarray((2.4, 7.5, 3.1)),
        homogeneous_growth_rate=-0.1, maximum_spatial_growth_rate=-0.01,
    )
    assert result["terminal_polarity_status"] == "finitely_polarizable"


def test_wgd_task_grid_is_paired_and_shardable():
    class Args:
        lambda_min = 10.0
        lambda_max = 20.0
        lambda_count = 2
        kappa_min = 0.0
        kappa_max = 0.1
        kappa_count = 2
        hill_values = "2,3"
        volume_multipliers = "0.5,1,2"

    tasks = all_tasks(Args())
    assert len(tasks) == 24
    keys = {task_key(task) for task in tasks}
    assert len(keys) == len(tasks)
    assert {key[0] for key in keys} == {0.5, 1.0, 2.0}


def test_wgd_transition_groups_match_single_species_primary_outcomes():
    assert transition_group("single_polarized", "dual_polarized") == "single_to_dual"
    assert transition_group("homogeneous_low", "single_polarized") == "gained_polarity"
    assert transition_group("single_polarized", "homogeneous_low") == "lost_polarity"


def test_extended_wgd_suite_matches_latest_single_species_task_shape():
    class Args:
        volume_multipliers = "0.5,1,2"
        hill_values = "2,3,4"
        asymmetries = "0,0.001,0.01,0.05"
        dual_kappa_min = 0.0
        dual_kappa_max = 0.08
        dual_kappa_count = 13
        dual_lambda_min = 31.6227766017
        dual_lambda_max = 94.8683298051
        dual_lambda_count = 13
        reversal_kappa_min = 0.0
        reversal_kappa_max = 0.08
        reversal_kappa_count = 17
        reversal_lambda_min = 31.6227766017
        reversal_lambda_max = 94.8683298051
        reversal_lambda_count = 17

    tasks = all_extended_tasks(Args())
    assert len(tasks) == 8685
    assert sum(task[0] == "dual" for task in tasks) == 6084
    assert sum(task[0] == "reversal" for task in tasks) == 2601


def test_sparse_and_dense_bdf_paths_agree_for_six_field_model():
    parameters = canonical_holmes_parameters()
    active = homogeneous_equilibrium(parameters)
    chemistry = HolmesChemistry(
        replace(parameters, initial_active=tuple(float(value) for value in active)),
        cue=GaussianRacPulse(start=0.5, stop=1.0),
    )
    model = MigrationModel(
        mechanics=MechanicsParameters(reference_length=20.0), chemistry=chemistry
    )
    common = dict(
        nodes=21, duration=2.0, output_interval=0.5, max_step=0.1,
        noise_amplitude=0.0, mechanics_activation_time=None,
    )
    dense = run_simulation(model, SimulationConfig(**common))
    sparse = run_simulation(
        model, SimulationConfig(**common, jacobian_sparsity=True)
    )
    np.testing.assert_allclose(sparse.membrane, dense.membrane, rtol=2e-6, atol=2e-8)
    np.testing.assert_allclose(sparse.cytoplasm, dense.cytoplasm, rtol=2e-6, atol=2e-8)
