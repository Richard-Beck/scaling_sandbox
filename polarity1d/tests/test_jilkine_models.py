import numpy as np
import pytest

from polarity1d import (
    Goryachev,
    Grid1D,
    LEGI,
    LocalizedPulse,
    Otsuji,
    WavePinning,
    noisy_initial_state,
    simulate,
)


@pytest.mark.parametrize("model", [WavePinning(), Otsuji(), Goryachev(), LEGI()])
def test_all_jilkine_models_run_through_one_contract(model):
    trajectory = simulate(
        model,
        protocol=LocalizedPulse(amplitude=0.1, stop=10.0, taper=5.0),
        length=10.0,
        cells=41,
        end_time=20.0,
        output_interval=2.0,
        backend="bdf",
    )
    assert trajectory.model == model.name
    assert trajectory.observable.shape == (11, 41)
    assert np.min(trajectory.observable) >= -1e-8
    assert np.ptp(trajectory.observable[-1]) > 1e-4


@pytest.mark.parametrize("model", [WavePinning(), Otsuji(), Goryachev()])
def test_mass_conserved_jilkine_models_conserve_total_material(model):
    trajectory = simulate(
        model,
        protocol=LocalizedPulse(amplitude=0.15, stop=15.0, taper=5.0),
        cells=61,
        end_time=50.0,
        output_interval=1.0,
        backend="bdf",
        rtol=2e-7,
        atol=1e-10,
    )
    total = trajectory.fields["active"] + trajectory.fields["inactive"]
    mass = np.sum(total, axis=1)
    assert np.max(np.abs(mass / mass[0] - 1.0)) < 2e-8


def test_split_solver_independently_matches_sparse_bdf_wave_pinning():
    protocol = LocalizedPulse(amplitude=0.15, stop=20.0, taper=10.0)
    common = dict(
        protocol=protocol,
        cells=61,
        end_time=60.0,
        output_interval=5.0,
    )
    adaptive = simulate(WavePinning(), backend="bdf", **common)
    split = simulate(WavePinning(), backend="split", **common)
    amplitude = np.ptp(adaptive.observable[-1])
    nrmse = np.sqrt(np.mean((adaptive.observable[-1] - split.observable[-1]) ** 2)) / amplitude
    assert nrmse < 0.004


def test_goryachev_assumed_initial_state_is_an_exact_homogeneous_fixed_point():
    model = Goryachev()
    u, v = model.homogeneous_state()
    reaction = model.reaction(
        0.0,
        np.array([0.5]),
        np.array([[u], [v]]),
        np.zeros(1),
    )
    assert np.max(np.abs(reaction)) < 1e-14


def test_noise_initialization_is_reproducible_and_preserves_local_mass():
    model = WavePinning()
    grid = Grid1D(cells=51)
    first = noisy_initial_state(model, grid, amplitude=0.05, seed=7)
    second = noisy_initial_state(model, grid, amplitude=0.05, seed=7)
    baseline = model.initial_state(grid)
    assert np.array_equal(first, second)
    assert np.isclose(np.mean(first[0]), np.mean(baseline[0]))
    assert np.allclose(np.sum(first, axis=0), np.sum(baseline, axis=0))
    trajectory = simulate(
        model,
        initial_state=first,
        cells=51,
        end_time=5.0,
        output_interval=1.0,
        backend="bdf",
    )
    assert trajectory.observable.shape == (6, 51)
