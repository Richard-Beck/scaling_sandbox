import numpy as np

from migration1d import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from migration1d.single_species import (
    ActivationPulse,
    DualPulseInitialCondition,
    SingleSpeciesParameters,
    StripWavePinningChemistry,
    classify_trajectory,
    homogeneous_equilibria,
    local_reaction_roots,
    polarized_domain_count,
)


def test_dimensionless_length_sets_active_diffusion_and_ratio():
    parameters = SingleSpeciesParameters(
        domain_length=10.0, dimensionless_length=20.0,
        inactivation_rate=2.0, diffusion_ratio=250.0,
    )
    assert np.isclose(parameters.membrane_diffusivity, 0.5)
    assert np.isclose(parameters.cytoplasm_diffusivity, 125.0)


def test_baseline_local_kinetics_have_three_roots_at_v_two():
    roots = local_reaction_roots(2.0, SingleSpeciesParameters())
    assert len(roots) == 3
    assert [root["stable"] for root in roots] == [True, False, True]
    np.testing.assert_allclose(
        [root["active"] for root in roots],
        [0.2683312, 0.32390323, 1.54176557], rtol=1e-6,
    )


def test_baseline_homogeneous_state_matches_mapped_resting_state():
    equilibria = homogeneous_equilibria(SingleSpeciesParameters())
    resting = min((item for item in equilibria if item["stable"]), key=lambda x: x["active"])
    assert np.isclose(resting["active"], 0.2683312, rtol=1e-6)
    assert np.isclose(resting["inactive"], 2.0, rtol=1e-6)


def test_strip_reaction_conserves_u_plus_hv():
    parameters = SingleSpeciesParameters(cytoplasm_depth=4.0)
    chemistry = StripWavePinningChemistry(
        parameters,
        activation=ActivationPulse(start=0.0, stop=10.0),
    )
    position = np.linspace(0.0, 10.0, 21)
    membrane = np.full((1, 21), chemistry.initial_membrane[0])
    cytoplasm = np.full((1, 21), chemistry.initial_cytoplasm[0])
    membrane_rate, cytoplasm_rate = chemistry.reaction(
        5.0, position, membrane, cytoplasm,
        np.zeros(21), np.zeros(21),
    )
    np.testing.assert_allclose(
        membrane_rate + parameters.cytoplasm_depth * cytoplasm_rate, 0.0
    )


def test_dual_initial_condition_is_mass_conserving_and_asymmetric():
    parameters = SingleSpeciesParameters(cytoplasm_depth=2.0)
    chemistry = StripWavePinningChemistry(
        parameters,
        initial_condition=DualPulseInitialCondition(asymmetry=0.05),
    )
    position = np.linspace(0.0, 10.0, 101)
    active, inactive = chemistry.initial_profiles(position)
    total = active + parameters.cytoplasm_depth * inactive
    np.testing.assert_allclose(total, parameters.total_density)
    assert active[0, 20] > active[0, 80]


def test_late_trajectory_classifier_distinguishes_one_and_two_domains():
    times = np.linspace(0.0, 100.0, 101)
    x = np.linspace(0.0, 1.0, 101)
    one = 0.2 + 1.5 * np.exp(-0.5 * ((x - 0.15) / 0.06) ** 2)
    two = one + 1.5 * np.exp(-0.5 * ((x - 0.85) / 0.06) ** 2)
    one_result = classify_trajectory(times, np.repeat(one[None, :], len(times), axis=0))
    two_result = classify_trajectory(times, np.repeat(two[None, :], len(times), axis=0))
    assert one_result["state"] == "single_polarized"
    assert two_result["state"] == "dual_polarized"


def test_polarized_domain_count_rejects_weak_ripples():
    x = np.linspace(0.0, 1.0, 101)
    weak_two_peak_profile = 1.0 + 0.01 * (
        np.exp(-0.5 * ((x - 0.2) / 0.04) ** 2)
        + np.exp(-0.5 * ((x - 0.8) / 0.04) ** 2)
    )
    assert polarized_domain_count(weak_two_peak_profile) == 0


def test_short_strip_simulation_conserves_geometry_weighted_amount():
    parameters = SingleSpeciesParameters(cytoplasm_depth=3.0)
    chemistry = StripWavePinningChemistry(
        parameters,
        activation=ActivationPulse(start=1.0, stop=3.0),
    )
    result = run_simulation(
        MigrationModel(
            mechanics=MechanicsParameters(reference_length=10.0),
            chemistry=chemistry,
        ),
        SimulationConfig(
            nodes=31, duration=8.0, output_interval=1.0, max_step=0.2,
            noise_amplitude=0.0, mechanics_activation_time=None,
        ),
    )
    membrane, cytoplasm = result.compartment_totals()
    total = membrane[:, 0] + parameters.cytoplasm_depth * cytoplasm[:, 0]
    assert np.max(np.abs(total - total[0])) / total[0] < 1e-8
