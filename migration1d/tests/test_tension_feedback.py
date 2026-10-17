import numpy as np

from migration1d import (
    EroumeChemistry,
    GlobalTensionFeedback,
    HolmesChemistry,
    MechanicsParameters,
    MigrationModel,
    SimulationConfig,
    SpatialTensionFeedback,
    StripWavePinningChemistry,
    run_simulation,
)
from migration1d.model import reaction_diffusion_jacobian_sparsity


class ConservedPairChemistry:
    """Minimal active/inactive pair used to isolate feedback from base kinetics."""

    species = ("nominal_polarity_species",)
    membrane_diffusion = np.asarray([0.0])
    cytoplasm_diffusion = np.asarray([0.0])
    initial_membrane = np.asarray([1.0])
    initial_cytoplasm = np.asarray([0.0])

    def reaction(self, time, position, membrane, cytoplasm, strain, stress):
        return np.zeros_like(membrane), np.zeros_like(cytoplasm)

    def protrusive_activity(self, membrane, cytoplasm):
        return np.maximum(membrane[0], 0.0)

    def tension_deactivation(self, membrane, cytoplasm, added_rate):
        loss = np.asarray(added_rate) * np.maximum(membrane[0], 0.0)
        return -loss[None, :], loss[None, :]


class LocalizedPairChemistry(ConservedPairChemistry):
    def initial_profiles(self, position):
        fraction = (position - position[0]) / (position[-1] - position[0])
        active = np.exp(-0.5 * ((fraction - 0.5) / 0.06) ** 2)[None, :]
        return active, 1.0 - active


def run_pair(*, left_force, right_force, spring, dashpot, feedback):
    return run_simulation(
        MigrationModel(
            mechanics=MechanicsParameters(
                reference_length=10.0,
                spring_stiffness=spring,
                dashpot_viscosity=dashpot,
            ),
            chemistry=ConservedPairChemistry(),
            left_force=left_force,
            right_force=right_force,
            tension_feedback=feedback,
        ),
        SimulationConfig(
            nodes=11,
            duration=0.5,
            output_interval=0.1,
            max_step=0.01,
            noise_amplitude=0.0,
        ),
    )


def test_global_feedback_responds_to_extension_but_not_compression():
    feedback = GlobalTensionFeedback(alpha_t=2.0)
    extended = run_pair(
        left_force=lambda state: -2.0,
        right_force=lambda state: 2.0,
        spring=3.0,
        dashpot=1.0,
        feedback=feedback,
    )
    compressed = run_pair(
        left_force=lambda state: 2.0,
        right_force=lambda state: -2.0,
        spring=3.0,
        dashpot=1.0,
        feedback=feedback,
    )
    extended_active, extended_inactive = extended.compartment_totals()
    compressed_active, _ = compressed.compartment_totals()

    assert extended.elastic_tension[-1] > 0.0
    assert extended_active[-1, 0] < extended_active[0, 0]
    assert compressed.spring_force[-1] < 0.0
    assert compressed.elastic_tension[-1] == 0.0
    np.testing.assert_allclose(compressed_active[:, 0], compressed_active[0, 0])
    np.testing.assert_allclose(
        extended_active + extended_inactive,
        np.repeat(
            (extended_active[0] + extended_inactive[0])[None, :],
            len(extended_active),
            axis=0,
        ),
        rtol=2e-8,
    )


def test_dashpot_tension_alone_does_not_deactivate_active_species():
    result = run_pair(
        left_force=lambda state: -2.0,
        right_force=lambda state: 2.0,
        spring=0.0,
        dashpot=3.0,
        feedback=GlobalTensionFeedback(alpha_t=10.0),
    )
    active, _ = result.compartment_totals()
    assert np.max(result.dashpot_force) > 0.0
    np.testing.assert_array_equal(result.elastic_tension, 0.0)
    np.testing.assert_allclose(active[:, 0], active[0, 0])


def test_uniform_spatial_tension_matches_generation_relaxation_solution():
    feedback = SpatialTensionFeedback(
        diffusivity=4.0,
        relaxation_time=2.0,
        generation_rate=0.5,
        alpha_t=0.0,
        baseline=3.0,
    )
    result = run_simulation(
        MigrationModel(
            chemistry=ConservedPairChemistry(), tension_feedback=feedback
        ),
        SimulationConfig(
            nodes=21,
            duration=4.0,
            output_interval=0.5,
            max_step=0.05,
            noise_amplitude=0.0,
            mechanics_activation_time=None,
            jacobian_sparsity=True,
        ),
    )
    expected = (
        feedback.generation_rate
        * feedback.relaxation_time
        * (1.0 - np.exp(-result.times / feedback.relaxation_time))
    )
    expected_field = np.repeat(expected[:, None], result.config.nodes, axis=1)
    np.testing.assert_allclose(result.tension_excess, expected_field, rtol=2e-5)
    np.testing.assert_allclose(
        result.tension, feedback.baseline + expected_field, rtol=2e-5
    )
    assert np.isclose(feedback.characteristic_length, np.sqrt(8.0))


def test_spatial_diffusion_propagates_tension_without_changing_its_integral():
    common = dict(
        relaxation_time=10.0,
        generation_rate=1.0,
        alpha_t=0.0,
    )
    config = SimulationConfig(
        nodes=41,
        duration=0.5,
        output_interval=0.5,
        max_step=0.01,
        noise_amplitude=0.0,
        mechanics_activation_time=None,
    )
    local = run_simulation(
        MigrationModel(
            chemistry=LocalizedPairChemistry(),
            tension_feedback=SpatialTensionFeedback(diffusivity=0.0, **common),
        ),
        config,
    )
    propagated = run_simulation(
        MigrationModel(
            chemistry=LocalizedPairChemistry(),
            tension_feedback=SpatialTensionFeedback(diffusivity=2.0, **common),
        ),
        config,
    )
    center = config.nodes // 2
    assert propagated.tension_excess[-1, center] < local.tension_excess[-1, center]
    assert propagated.tension_excess[-1, center + 6] > local.tension_excess[-1, center + 6]
    weights = np.diff(propagated.position[-1]).mean()
    # Trapezoidal endpoint weights are half-width under this uniform grid.
    quadrature = np.ones(config.nodes) * weights
    quadrature[[0, -1]] *= 0.5
    np.testing.assert_allclose(
        np.sum(propagated.tension_excess[-1] * quadrature),
        np.sum(local.tension_excess[-1] * quadrature),
        rtol=2e-6,
    )


def test_spatial_tension_suppresses_activity_and_preserves_total_pool():
    config = SimulationConfig(
        nodes=11,
        duration=1.0,
        output_interval=0.2,
        max_step=0.02,
        noise_amplitude=0.0,
        mechanics_activation_time=None,
    )
    control = run_simulation(
        MigrationModel(
            chemistry=ConservedPairChemistry(),
            tension_feedback=SpatialTensionFeedback(
                diffusivity=1.0,
                relaxation_time=2.0,
                generation_rate=1.0,
                alpha_t=0.0,
            ),
        ),
        config,
    )
    inhibited = run_simulation(
        MigrationModel(
            chemistry=ConservedPairChemistry(),
            tension_feedback=SpatialTensionFeedback(
                diffusivity=1.0,
                relaxation_time=2.0,
                generation_rate=1.0,
                alpha_t=3.0,
            ),
        ),
        config,
    )
    active_control, _ = control.compartment_totals()
    active, inactive = inhibited.compartment_totals()
    assert active[-1, 0] < active_control[-1, 0]
    np.testing.assert_allclose(
        active[:, 0] + inactive[:, 0],
        np.full(len(active), active[0, 0] + inactive[0, 0]),
        rtol=2e-8,
    )


def test_builtin_chemistries_apply_mass_balanced_tension_deactivation():
    nodes = 5
    rate = np.linspace(0.0, 0.4, nodes)

    single = StripWavePinningChemistry()
    single_m = np.repeat(single.initial_membrane[:, None], nodes, axis=1)
    single_c = np.repeat(single.initial_cytoplasm[:, None], nodes, axis=1)
    dm, dc = single.tension_deactivation(single_m, single_c, rate)
    np.testing.assert_allclose(
        dm + single.parameters.cytoplasm_depth * dc, 0.0
    )

    holmes = HolmesChemistry()
    holmes_m = np.repeat(holmes.initial_membrane[:, None], nodes, axis=1)
    holmes_c = np.repeat(holmes.initial_cytoplasm[:, None], nodes, axis=1)
    dm, dc = holmes.tension_deactivation(holmes_m, holmes_c, rate)
    np.testing.assert_allclose(dm + dc, 0.0)
    np.testing.assert_array_equal(dm[[0, 2]], 0.0)

    eroume = EroumeChemistry()
    eroume_m = np.repeat(eroume.initial_membrane[:, None], nodes, axis=1)
    eroume_c = np.zeros_like(eroume_m)
    dm, dc = eroume.tension_deactivation(eroume_m, eroume_c, rate)
    np.testing.assert_allclose(dm[2] + dm[3], 0.0)
    np.testing.assert_array_equal(dc, 0.0)


def test_spatial_tension_extends_jacobian_sparsity_by_one_field():
    baseline = reaction_diffusion_jacobian_sparsity(species=3, nodes=11)
    spatial = reaction_diffusion_jacobian_sparsity(
        species=3, nodes=11, spatial_tension=True
    )
    assert baseline.shape == (68, 68)
    assert spatial.shape == (79, 79)
    assert np.all(spatial.diagonal())
