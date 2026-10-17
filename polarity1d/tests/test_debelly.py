import numpy as np

from polarity1d import (
    DeBelly,
    DeBellyParameters,
    GaussianPulse,
    LocalizedPulse,
    ReversingGradient,
    simulate,
)


def test_debelly_baseline_is_symmetric_and_nonnegative():
    trajectory = simulate(
        DeBelly(),
        length=20.0,
        cells=60,
        end_time=20.0,
        output_interval=2.0,
        step=0.25,
    )
    assert np.min(trajectory.fields["rac"]) >= 0
    assert np.min(trajectory.fields["rho"]) >= 0
    assert np.max(np.abs(trajectory.fields["rac"] - trajectory.fields["rac"][:, ::-1])) < 1e-8
    assert np.max(np.abs(trajectory.fields["rho"] - trajectory.fields["rho"][:, ::-1])) < 1e-8


def test_debelly_mca_is_conserved_to_roundoff_under_mechanical_transport():
    trajectory = simulate(
        DeBelly(),
        protocol=GaussianPulse(start=5.0, stop=15.0, width_fraction=0.1),
        length=20.0,
        cells=60,
        end_time=50.0,
        output_interval=1.0,
        step=0.25,
    )
    total = trajectory.fields["mca_bound"] + trajectory.fields["mca_unbound"]
    mass = np.sum(total, axis=1)
    assert np.max(np.abs(mass / mass[0] - 1.0)) < 2e-14
    assert trajectory.metadata["mca_conservation_error"] < 2e-14


def test_debelly_arm_ablation_removes_tension_and_changes_distal_rac():
    protocol = GaussianPulse(start=10.0, stop=30.0, width_fraction=0.1)
    intact = simulate(
        DeBelly(DeBellyParameters(protrusion_speed=0.0)),
        protocol=protocol,
        length=20.0,
        cells=60,
        end_time=100.0,
        output_interval=2.0,
        step=0.25,
    )
    ablated = simulate(
        DeBelly(
            DeBellyParameters(
                protrusion_speed=0.0,
                rac_to_rho_coupling=0.0,
                rho_to_rac_coupling=0.0,
            )
        ),
        protocol=protocol,
        length=20.0,
        cells=60,
        end_time=100.0,
        output_interval=2.0,
        step=0.25,
    )
    distal = intact.x >= 0.8 * 20.0
    assert np.max(np.abs(intact.fields["tension"])) > 1e-3
    assert np.max(np.abs(ablated.fields["tension"])) < 1e-12
    assert abs(np.mean(intact.observable[-1, distal]) - np.mean(ablated.observable[-1, distal])) > 1e-3


def test_debelly_spatial_refinement_is_close():
    protocol = GaussianPulse(start=5.0, stop=15.0, width_fraction=0.15)
    coarse = simulate(
        DeBelly(), protocol=protocol, cells=50, end_time=40.0, output_interval=5.0, step=0.25
    )
    fine = simulate(
        DeBelly(), protocol=protocol, cells=100, end_time=40.0, output_interval=5.0, step=0.25
    )
    interpolated = np.interp(fine.x, coarse.x, coarse.observable[-1])
    nrmse = np.sqrt(np.mean((interpolated - fine.observable[-1]) ** 2)) / np.ptp(fine.observable[-1])
    assert nrmse < 0.025


def test_debelly_common_stimulus_can_target_rho():
    protocol = GaussianPulse(start=2.0, stop=6.0, width_fraction=0.1)
    trajectory = simulate(
        DeBelly(DeBellyParameters(opto_target="rho", protrusion_speed=0.0)),
        protocol=protocol,
        cells=40,
        end_time=10.0,
        output_interval=1.0,
        step=0.5,
    )
    assert np.max(trajectory.fields["rho"][-1] - trajectory.fields["rho"][0]) > 1e-3


def test_debelly_production_initialization_is_algebraic_and_step_independent():
    coarse = simulate(
        DeBelly(), cells=20, end_time=1.0, output_interval=1.0, step=0.5
    )
    fine = simulate(
        DeBelly(), cells=20, end_time=1.0, output_interval=1.0, step=0.1
    )
    assert coarse.metadata["parameters"]["initialization_mode"] == "algebraic_fixed_point"
    assert np.allclose(coarse.fields["rac"][0], fine.fields["rac"][0], atol=1e-14)
    assert np.allclose(coarse.fields["rho"][0], fine.fields["rho"][0], atol=1e-14)


def test_debelly_accepts_explicit_rac_rho_initial_fields():
    cells = 20
    rac = np.linspace(2.0, 0.1, cells)
    rho = np.linspace(0.1, 2.0, cells)
    trajectory = simulate(
        DeBelly(
            DeBellyParameters(assay_mode="reflection_symmetric_extension")
        ),
        length=4.0,
        cells=cells,
        end_time=1.0,
        output_interval=1.0,
        initial_fields={"rac": rac, "rho": rho},
    )
    assert np.allclose(trajectory.fields["rac"][0], rac)
    assert np.allclose(trajectory.fields["rho"][0], rho)
    assert trajectory.metadata["custom_initial_fields"] == ["rac", "rho"]


def test_debelly_algebraic_initialization_resolves_picard_two_cycle():
    trajectory = simulate(
        DeBelly(
            DeBellyParameters(
                rac_basal_activation=2.0,
                rho_basal_activation=2.0,
                rac_mechanical_feedback=0.0,
                rho_mechanical_feedback=0.0,
            )
        ),
        cells=20,
        end_time=0.5,
        output_interval=0.5,
        step=0.5,
    )
    rac = float(trajectory.fields["rac"][0, 0])
    rho = float(trajectory.fields["rho"][0, 0])
    assert abs(rac - 2.0 / (1.0 + rho**2)) < 1e-10
    assert abs(rho - 2.0 / (1.0 + rac**2)) < 1e-10
    assert np.ptp(trajectory.fields["rac"][0]) < 1e-14
    assert np.ptp(trajectory.fields["rho"][0]) < 1e-14


def test_debelly_julia_validation_initialization_is_explicit_and_step_dependent():
    model = DeBelly(DeBellyParameters(initialization_mode="julia_warmup"))
    coarse = simulate(model, cells=20, end_time=1.0, output_interval=1.0, step=0.5)
    fine = simulate(model, cells=20, end_time=1.0, output_interval=1.0, step=0.1)
    assert coarse.metadata["parameters"]["initialization_mode"] == "julia_warmup"
    assert not np.allclose(coarse.fields["rac"][0], fine.fields["rac"][0])


def test_debelly_julia_validation_input_sampling_matches_element_right_edges():
    length, cells = 5.0, 25
    protocol = LocalizedPulse(
        amplitude=1.0, start=0.0, stop=40.0, taper=20.0,
        width_fraction=0.1, sides=("left", "right"), asymmetry=0.02,
    )
    centred = simulate(
        DeBelly(), protocol=protocol, length=length, cells=cells,
        end_time=0.5, output_interval=0.5, step=0.5,
    )
    julia_sampled = simulate(
        DeBelly(DeBellyParameters(input_sampling_mode="julia_element_right_edge")),
        protocol=protocol, length=length, cells=cells,
        end_time=0.5, output_interval=0.5, step=0.5,
    )
    dx = length / cells
    expected_difference = (
        protocol(0.0, np.minimum(centred.x + 0.5 * dx, length), length)
        - protocol(0.0, centred.x, length)
    )
    observed_difference = (
        julia_sampled.fields["rac_activation"][0]
        - centred.fields["rac_activation"][0]
    )
    assert np.allclose(observed_difference, expected_difference, atol=1e-14)
    assert julia_sampled.metadata["parameters"]["input_sampling_mode"] == (
        "julia_element_right_edge"
    )


def test_debelly_boundary_relaxation_is_the_released_per_step_factor():
    parameters = DeBellyParameters(boundary_relaxation=0.98)
    assert parameters.boundary_relaxation == 0.98


def test_debelly_moving_boundary_depletes_front_mca_without_losing_total_mca():
    length = 15.0
    trajectory = simulate(
        DeBelly(DeBellyParameters(initialization_mode="julia_warmup")),
        protocol=ReversingGradient(
            slope=1.0 / length,
            reversal_time=2_000.0,
            transition=20.0,
        ),
        length=length,
        cells=75,
        end_time=200.0,
        output_interval=20.0,
        step=0.5,
    )
    bound = trajectory.fields["mca_bound"][-1]
    left = np.mean(bound[trajectory.x <= 0.2 * length])
    right = np.mean(bound[trajectory.x >= 0.8 * length])
    orientation = (right - left) / (right + left)
    total = trajectory.fields["mca_bound"] + trajectory.fields["mca_unbound"]
    mass = np.sum(total, axis=1)
    assert orientation > 0.1
    assert np.max(np.abs(mass / mass[0] - 1.0)) < 2e-14


def _symmetric_pulse(side: str, **parameter_overrides):
    length = 15.0
    return simulate(
        DeBelly(
            DeBellyParameters(
                assay_mode="reflection_symmetric_extension",
                **parameter_overrides,
            )
        ),
        protocol=LocalizedPulse(
            amplitude=1.0,
            start=0.0,
            stop=40.0,
            taper=20.0,
            width_fraction=0.1,
            sides=(side,),
        ),
        length=length,
        cells=75,
        end_time=500.0,
        output_interval=1.0,
        step=0.5,
    )


def test_debelly_symmetric_extension_is_reflection_equivariant():
    left = _symmetric_pulse("left")
    right = _symmetric_pulse("right")
    for field in ("rac", "rho", "tension", "mca_bound", "mca_unbound"):
        assert np.allclose(
            left.fields[field], right.fields[field][:, ::-1], atol=3e-13, rtol=0
        )
    for field in ("velocity", "displacement"):
        assert np.allclose(
            left.fields[field], -right.fields[field][:, ::-1], atol=3e-13, rtol=0
        )
    for field in ("boundary_separation", "mean_strain", "tension_range"):
        assert np.allclose(left.fields[field], right.fields[field], atol=3e-13, rtol=0)


def test_debelly_exports_actual_mechanical_feedback_increments():
    trajectory = _symmetric_pulse("left")
    parameters = DeBellyParameters(assay_mode="reflection_symmetric_extension")
    expected_alpha = parameters.rac_mechanical_feedback * 0.5 * (
        1.0
        - np.tanh(
            (trajectory.fields["mca_bound"] - parameters.mca_rac_threshold)
            / parameters.mca_switch_scale
        )
    )
    expected_beta = parameters.rho_mechanical_feedback * 0.5 * (
        1.0
        + np.tanh(
            (trajectory.fields["tension"] - parameters.tension_rho_threshold)
            / parameters.tension_switch_scale
        )
    )
    assert np.allclose(trajectory.fields["delta_alpha_mca"], expected_alpha)
    assert np.allclose(trajectory.fields["delta_beta_tension"], expected_beta)


def test_debelly_symmetric_extension_has_both_mirrored_mechanical_paths():
    length = 15.0
    for side in ("left", "right"):
        full = _symmetric_pulse(side)
        no_tension_to_rho = _symmetric_pulse(side, rho_mechanical_feedback=0.0)
        no_mca_to_rac = _symmetric_pulse(side, rac_mechanical_feedback=0.0)
        proximal = full.x <= 0.2 * length if side == "left" else full.x >= 0.8 * length
        distal = full.x >= 0.8 * length if side == "left" else full.x <= 0.2 * length

        full_distal_rho = np.mean(full.fields["rho"][:, distal], axis=1)
        beta0_distal_rho = np.mean(no_tension_to_rho.fields["rho"][:, distal], axis=1)
        distal_tension = np.mean(full.fields["tension"][:, distal], axis=1)
        proximal_mca = np.mean(full.fields["mca_bound"][:, proximal], axis=1)
        full_proximal_rac = np.mean(full.fields["rac"][:, proximal], axis=1)
        alpha0_proximal_rac = np.mean(
            no_mca_to_rac.fields["rac"][:, proximal], axis=1
        )

        assert np.max(distal_tension - distal_tension[0]) > 5.0
        assert np.max(full_distal_rho - full_distal_rho[0]) > 0.03
        assert np.max(beta0_distal_rho - beta0_distal_rho[0]) < 0.005
        assert np.min(proximal_mca - proximal_mca[0]) < -5e-4
        assert full_proximal_rac[-1] - full_proximal_rac[0] > 0.002
        assert abs(alpha0_proximal_rac[-1] - alpha0_proximal_rac[0]) < 5e-4
