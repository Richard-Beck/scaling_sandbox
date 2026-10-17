from __future__ import annotations

import numpy as np

from taniguchi.geometry_campaign import geometry_scales, run_geometry_campaign
from taniguchi.fixed import (
    FixedParameters,
    TriangularGrid,
    gaussian_transfer,
    simulate_fixed,
)
from taniguchi.phase_field import (
    PHASE_PARAMETER_SETS,
    PhaseFieldResult,
    simulate_phase_field,
)
from taniguchi.render import (
    MicroscopyRenderer,
    RenderScale,
    determine_scale,
    representative_indices,
    still_indices,
    video_indices,
)
from taniguchi.wgd_metrics import summarize_phase_wgd


def test_table_s1_parameter_sets_are_complete():
    assert set(PHASE_PARAMETER_SETS) == {"4F", "4G", "4H", "4I"}
    assert PHASE_PARAMETER_SETS["4G"].chi_u == 80.0
    assert PHASE_PARAMETER_SETS["4I"].theta == 0.85e-3
    assert np.isclose(PHASE_PARAMETER_SETS["4F"].area0, 25 * np.pi)


def test_gaussian_impulse_is_exactly_balanced():
    grid = TriangularGrid.circle(radius=1.0, dx=0.2)
    U = np.full(grid.mask.shape, 4.0)
    V = np.zeros_like(U)
    before = (U + V).copy()
    gaussian_transfer(U, V, grid, 0.1, -0.1, 0.7, 0.2)
    np.testing.assert_allclose((U + V)[grid.mask], before[grid.mask], atol=1e-14)


def test_homogeneous_rest_is_stable_without_firing():
    p = FixedParameters()
    result = simulate_fixed(0.05, parameters=p,
                            grid=TriangularGrid.circle(radius=1.0, dx=0.2),
                            stochastic=False, output_interval=0.05)
    np.testing.assert_allclose(result.U[-1][result.grid.mask], 1 / p.gamma,
                               atol=1e-12)
    np.testing.assert_allclose(result.V[-1][result.grid.mask], 0.0, atol=1e-12)


def test_subthreshold_and_suprathreshold_responses_differ():
    grid = TriangularGrid.circle(radius=2.0, dx=0.2)
    low = simulate_fixed(1.0, grid=grid, stochastic=False, output_interval=0.05,
                         initial_transfers=((0, 0, 0.2, 0.3),))
    high = simulate_fixed(1.0, grid=grid, stochastic=False, output_interval=0.05,
                          initial_transfers=((0, 0, 1.0, 0.3),))
    assert low.V.max() < 0.25
    assert high.V.max() > 2.0


def test_phase_field_short_run_is_finite_and_bounded():
    result = simulate_phase_field(8e-4, dx=0.2, dt=8e-5,
                                  box_half_width=6.5, stochastic=False,
                                  output_interval=8e-4)
    assert np.isfinite(result.phi).all()
    assert np.isfinite(result.U).all()
    assert np.isfinite(result.V).all()
    assert result.phi.min() >= 0.0
    assert result.phi.max() <= 1.0


def test_fast_phase_kernel_matches_reference_kernel():
    kwargs = dict(
        end_time=8e-4,
        dx=0.2,
        dt=8e-5,
        box_half_width=6.5,
        stochastic=False,
        output_interval=8e-4,
        initial_transfers=((1.0, 0.0, 1.0, 0.3),),
    )
    reference = simulate_phase_field(**kwargs, backend="reference")
    fast = simulate_phase_field(**kwargs, backend="fast")
    np.testing.assert_allclose(fast.phi, reference.phi, atol=1e-12)
    np.testing.assert_allclose(fast.U, reference.U, atol=1e-12)
    np.testing.assert_allclose(fast.V, reference.V, atol=1e-12)


def test_phase_backend_is_validated():
    try:
        simulate_phase_field(0.0, backend="not-a-backend")
    except ValueError as exc:
        assert "backend" in str(exc)
    else:
        raise AssertionError("invalid backend was accepted")


def test_phase_time_and_spacing_are_validated():
    for kwargs in ({"end_time": -1.0}, {"end_time": 0.0, "dt": 0.0},
                   {"end_time": 0.0, "dx": 0.0},
                   {"end_time": 0.0, "output_interval": 0.0}):
        try:
            simulate_phase_field(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid arguments were accepted: {kwargs}")


def test_render_scale_ignores_exterior_and_is_movie_wide():
    phi = np.zeros((2, 4, 4))
    phi[:, 1:3, 1:3] = 1
    U = np.full_like(phi, 5.0)
    V = np.zeros_like(phi)
    V[1, 2, 2] = 4.0
    U[:, 0, 0] = 1e6
    scale = determine_scale(U, V, phi, percentile=100)
    assert scale.pip2_max == 5.0
    assert scale.pip3_max == 4.0


def test_frame_selection_includes_trajectory_endpoints():
    np.testing.assert_array_equal(video_indices(101, fps=2, duration=3),
                                  [0, 20, 40, 60, 80, 100])
    np.testing.assert_array_equal(still_indices(101, 5), [0, 25, 50, 75, 100])


def test_representative_frames_select_separated_activity_peaks():
    phi = np.ones((11, 2, 2))
    V = np.zeros_like(phi)
    V[2] = 2.0
    V[7] = 3.0
    selected = representative_indices(V, phi, 4)
    np.testing.assert_array_equal(selected, [0, 2, 7, 10])


def test_microscopy_renderer_produces_even_rgb_frame():
    x = np.linspace(-2, 2, 21)
    X, Y = np.meshgrid(x, x)
    phi = (np.hypot(X, Y) < 1.5).astype(float)
    U = 1.0 + np.zeros_like(phi)
    V = np.exp(-((X - 0.7) ** 2 + Y ** 2) / 0.2)
    renderer = MicroscopyRenderer(x, x, condition="4F", mode="pip3",
                                  scale=RenderScale(1.0, 1.0), pixels=128)
    frame = renderer.render(phi, U, V, 0.5)
    assert frame.shape == (128, 128, 3)
    assert frame.dtype == np.uint8
    assert frame.max() > 0


def test_wgd_metrics_count_domains_and_directional_activity():
    x = np.linspace(-2, 2, 9)
    X, Y = np.meshgrid(x, x)
    phi = np.repeat((np.hypot(X, Y) <= 1.8)[None, :, :], 3, axis=0).astype(float)
    U = np.full_like(phi, 1.0)
    V = np.zeros_like(phi)
    V[1, 3:6, 5:7] = 2.0
    V[2, 3:6, 5:7] = 2.0
    result = PhaseFieldResult(
        x=x, y=x, parameters=PHASE_PARAMETER_SETS["4F"],
        times=np.asarray([0.0, 1.0, 2.0]), phi=phi, U=U, V=V,
        event_count=2,
    )
    metrics = summarize_phase_wgd(
        result, activity_threshold=0.5, minimum_domain_area=0.25
    )
    assert metrics["state_class"] == "excitable"
    assert metrics["terminal_n_domains"] == 1
    assert metrics["domain_birth_count"] == 1
    assert metrics["stochastic_event_rate_per_time"] == 1.0
    assert metrics["terminal_directional_moment"] > 0


def test_geometry_campaign_uses_volume_to_linear_and_area_scaling():
    linear, area = geometry_scales(2.0)
    assert np.isclose(linear, 2 ** (1 / 3))
    assert np.isclose(area, 2 ** (2 / 3))


def test_fixed_campaign_holds_phi_and_writes_compact_outputs(tmp_path):
    config = {
        "arm": "fixed_controlled",
        "volume_ratio": 1,
        "seed": None,
        "model_config": {
            "end_time": 0.001,
            "dt": 0.0005,
            "dx": 0.5,
            "output_interval": 0.0005,
            "baseline_box_half_width": 6.5,
            "radial_bins": 24,
            "snapshot_bins": 2,
            "trigger_width": 0.5,
        },
    }
    summary = run_geometry_campaign(config, tmp_path)
    assert summary["boundary_mode"] == "frozen_diffuse_circle"
    assert summary["event_count"] == 0
    with np.load(tmp_path / "snapshots.npz", allow_pickle=False) as snapshots:
        np.testing.assert_allclose(
            snapshots["phi"],
            np.repeat(snapshots["phi"][[0]], snapshots["phi"].shape[0], axis=0),
        )
    with np.load(tmp_path / "trajectory.npz", allow_pickle=False) as trajectory:
        assert str(trajectory["format_version"]) == "taniguchi_geometry_compact_v1"
        assert trajectory["times"].size == 3


def test_per_cell_firing_schedule_is_identical_across_paired_sizes(tmp_path):
    common = {
        "end_time": 0.001,
        "dt": 0.0005,
        "dx": 0.5,
        "output_interval": 0.0005,
        "baseline_box_half_width": 6.5,
        "radial_bins": 24,
        "snapshot_bins": 2,
        "firing_rule": "per_cell_reference",
        # Force several events in this very short unit test.
        "reference_event_area": 1.0e7,
    }
    event_data = []
    for volume in (1, 2):
        output = tmp_path / f"{volume}x"
        summary = run_geometry_campaign({
            "campaign": "test_per_cell",
            "arm": "fixed_spontaneous",
            "volume_ratio": volume,
            "seed": 19,
            "model_config": common,
        }, output)
        assert summary["firing_rule"] == "per_cell_reference"
        assert summary["prescribed_event_rate_per_cell_time"] == 1.0e4
        with np.load(output / "events_and_tracks.npz", allow_pickle=False) as events:
            event_data.append((events["event_time"].copy(),
                               events["event_amplitude"].copy()))
    np.testing.assert_array_equal(event_data[0][0], event_data[1][0])
    np.testing.assert_array_equal(event_data[0][1], event_data[1][1])
    assert event_data[0][0].size > 0
