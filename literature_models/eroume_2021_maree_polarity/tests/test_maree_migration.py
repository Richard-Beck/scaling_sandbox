from __future__ import annotations

from dataclasses import replace

import numpy as np

from maree_migration import SimulationConfig, StimulusProtocol, reaction_rates, run_simulation
from maree_migration.metrics import material_totals, trajectory_metrics


def small_config(**kwargs):
    base = SimulationConfig(duration=0.04, dt=0.004, dx=1.0, box_size=32.0,
                            cell_radius=6.0, output_interval=0.02,
                            noise_amplitude=0.0)
    return replace(base, **kwargs)


def test_exact_reaction_stoichiometry_and_pi_source_are_preserved():
    c = np.asarray((1.4, 1.0, 4.5, 3.0, 1.85, 1.25, 55.0, 30.0, 0.05))[:, None, None]
    rates = reaction_rates(c, 2.95)
    np.testing.assert_allclose(rates[0] + rates[1], 0.0)
    np.testing.assert_allclose(rates[2] + rates[3], 0.0)
    np.testing.assert_allclose(rates[4] + rates[5], 0.0)
    # At P1=55, the intentionally retained source is 10.5 - .21*55 = -1.05.
    np.testing.assert_allclose(rates[6] + rates[7] + rates[8], -1.05)
    # Values independently evaluated from the C++/R equations at the VCML
    # initial state and uniform Ic=2.95.
    expected = (0.139583333333333, -0.139583333333333, 0.0, 0.0,
                -0.003225806451613, 0.003225806451613, -5.25, 4.2, 0.0)
    np.testing.assert_allclose(rates[:, 0, 0], expected, rtol=1e-13, atol=1e-14)


def test_stimulus_is_in_fixed_laboratory_coordinates_and_rotates():
    X, Y = np.meshgrid(np.asarray((40.0, 48.0, 56.0)), np.asarray((40.0, 48.0, 56.0)))
    protocol = StimulusProtocol(mode="rotate", pivot_x=48, pivot_y=48,
                                gradient_strength=0.05, gradient_angle=0,
                                rotate_at=5, rotate_angle=90)
    before = protocol.evaluate(X, Y, 4.9)
    after = protocol.evaluate(X, Y, 5.0)
    assert before[1, 2] > before[1, 0]
    assert after[2, 1] > after[0, 1]
    assert before[1, 1] == after[1, 1] == 2.95


def test_short_moving_run_is_finite_bounded_and_area_controlled():
    result = run_simulation(small_config(noise_amplitude=0.005))
    assert np.isfinite(result.phi).all()
    assert np.isfinite(result.concentrations).all()
    assert result.phi.min() >= 0 and result.phi.max() <= 1
    metrics = trajectory_metrics(result)
    assert np.ptp(metrics[:, 4]) / metrics[0, 4] < 0.01
    fraction = result.exterior_material / result.material
    assert np.max(fraction - fraction[0]) < 1e-3


def test_frozen_symmetric_domain_preserves_uniformity_and_gtpase_material():
    result = run_simulation(small_config(freeze_boundary=True))
    inside = result.phi[-1] >= 0.5
    for species in range(9):
        # Saved movie fields are intentionally float32.
        assert np.ptp(result.concentrations[-1, species][inside]) < 5e-5
    totals = material_totals(result)
    for pair in ((0, 1), (2, 3), (4, 5)):
        pair_total = totals[:, pair].sum(axis=1)
        # The q update is conservative; reconstructing totals from float32
        # movie concentrations incurs a small diffuse-tail cutoff error.
        np.testing.assert_allclose(pair_total, pair_total[0], rtol=2e-4)


def test_no_scaling_control_is_exposed():
    fields = SimulationConfig.__dataclass_fields__
    assert "linear_scale" not in fields and "size_scale" not in fields


def test_frozen_legacy_gradient_reproduces_source_direction():
    protocol = StimulusProtocol(mode="legacy")
    result = run_simulation(small_config(duration=1.0, output_interval=1.0,
                                         freeze_boundary=True, stimulus=protocol))
    active = result.concentrations[-1, 1]
    inside = result.phi[-1] >= 0.5
    X = np.broadcast_to(result.x, active.shape)
    assert active[inside & (X > 17)].mean() > active[inside & (X < 15)].mean()
