from __future__ import annotations

import numpy as np

from bulk_surface import (
    FixedGeometry,
    SimulationConfig,
    WrappedTopHatStimulus,
    wrapped_gaussian_profile,
)
from bulk_surface.examples import exchange_case, wave_pinning_case
from bulk_surface.render import extract_surface_profiles, render_run
from bulk_surface.shared_conditions import (
    SharedChemistry,
    chiou_condition,
    mori_condition,
)


def test_circle_geometry_and_operators_preserve_constants():
    geometry = FixedGeometry.circle(radius=2.0, dx=0.2, n_surface=96)
    assert abs(geometry.area / (np.pi * 2.0**2) - 1.0) < 0.01
    assert abs(geometry.perimeter / (2.0 * np.pi * 2.0) - 1.0) < 0.001
    assert np.allclose(
        geometry.trace_bulk(np.full(geometry.n_bulk, 3.25)), 3.25, atol=1e-14
    )
    assert np.max(np.abs(geometry.bulk_stiffness @ np.ones(geometry.n_bulk))) < 1e-12
    assert (
        np.max(np.abs(geometry.surface_stiffness @ np.ones(geometry.n_surface)))
        < 1e-12
    )


def test_wrapped_stimulus_crosses_membrane_storage_boundary():
    fraction = np.arange(100) / 100
    profile = wrapped_gaussian_profile(fraction, 0.99, 0.03)
    assert profile[99] == profile.max()
    assert profile[0] > 0.9
    assert profile[50] < 1e-20


def test_wrapped_top_hat_has_requested_periodic_support():
    fraction = np.arange(100) / 100
    values = np.empty(100)
    stimulus = WrappedTopHatStimulus(3.0, 0.0, 0.04, start=1.0, stop=2.0)
    stimulus.values(0.5, fraction, values)
    assert np.all(values == 0.0)
    stimulus.values(1.5, fraction, values)
    assert np.array_equal(np.flatnonzero(values), np.asarray([0, 1, 2, 98, 99]))


def test_reversible_exchange_is_conservative():
    case = exchange_case(FixedGeometry.circle(radius=1.5, dx=0.3, n_surface=48))
    result = case.simulator.run(
        case.initial,
        SimulationConfig(end_time=1.0, dt=0.02, output_interval=0.1),
    )
    total = result.pool_totals["combined"]
    assert np.max(np.abs(total / total[0] - 1.0)) < 2e-12
    assert np.mean(result.surface[-1, 0]) > np.mean(result.surface[0, 0])
    assert np.mean(result.bulk[-1, 0]) < np.mean(result.bulk[0, 0])


def test_final_state_is_saved_when_output_interval_does_not_divide_duration():
    case = exchange_case(FixedGeometry.circle(radius=1.0, dx=0.4, n_surface=32))
    result = case.simulator.run(
        case.initial,
        SimulationConfig(end_time=0.6, dt=0.1, output_interval=0.4),
    )
    assert np.allclose(result.times, np.asarray([0.0, 0.4, 0.6]), atol=1e-15)


def test_wave_pinning_configuration_builds_localized_membrane_response():
    case = wave_pinning_case(
        FixedGeometry.circle(radius=1.5, dx=0.25, n_surface=64)
    )
    result = case.simulator.run(
        case.initial,
        SimulationConfig(end_time=4.0, dt=0.02, output_interval=0.2),
    )
    total = result.pool_totals["gtpase"]
    initial_contrast = np.ptp(result.surface[0, 0])
    final_contrast = np.ptp(result.surface[-1, 0])
    assert np.max(np.abs(total / total[0] - 1.0)) < 3e-12
    assert final_contrast > 10 * initial_contrast


def test_shared_conditions_start_with_equal_total_amount_and_are_conservative(tmp_path):
    geometry = FixedGeometry.circle(radius=1.0, dx=0.3, n_surface=48)
    chemistry = SharedChemistry()
    mori = mori_condition(geometry, chemistry)
    chiou = chiou_condition(geometry, chemistry)
    config = SimulationConfig(end_time=0.2, dt=0.02, output_interval=0.1)
    mori_result = mori.simulator.run(mori.initial, config)
    chiou_result = chiou.simulator.run(chiou.initial, config)
    mori_total = mori_result.pool_totals["total_protein"]
    chiou_total = chiou_result.pool_totals["total_protein"]
    assert np.isclose(mori_total[0], chiou_total[0], atol=1e-13)
    assert np.max(np.abs(mori_total / mori_total[0] - 1.0)) < 1e-11
    assert np.max(np.abs(chiou_total / chiou_total[0] - 1.0)) < 1e-11

    profiles = extract_surface_profiles(chiou_result)
    assert profiles["bulk_trace"].shape == chiou_result.surface.shape
    diagnostics = render_run(
        chiou_result, tmp_path / "rendered", condition="chiou",
        snapshot_times=(0.0, 0.2),
    )
    assert diagnostics["maximum_relative_pool_error"] < 1e-11
    assert (tmp_path / "rendered" / "profiles.npz").is_file()
    assert len(list((tmp_path / "rendered" / "snapshots").glob("*.png"))) == 2
