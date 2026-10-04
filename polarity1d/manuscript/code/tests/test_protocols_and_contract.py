import numpy as np

from polarity1d import (
    FinitePulseSequence,
    Grid1D, LocalizedPulse, ReversingFiniteGradient, ReversingGradient, Trajectory,
)


def test_finite_pulse_sequence_has_cue_free_interval_and_both_sides():
    sequence = FinitePulseSequence(
        pulses=(
            {"amplitude": 1.0, "start": 0.0, "stop": 20.0, "taper": 0.0,
             "width_fraction": 0.1, "sides": ("left",)},
            {"amplitude": 1.0, "start": 40.0, "stop": 60.0, "taper": 0.0,
             "width_fraction": 0.1, "sides": ("right",)},
        )
    )
    grid = Grid1D(length=10.0, cells=50)
    assert np.all(sequence(30.0, grid.x, grid.length) == 0.0)
    assert np.argmax(sequence(10.0, grid.x, grid.length)) == 0
    assert np.argmax(sequence(50.0, grid.x, grid.length)) == grid.cells - 1
    assert sequence.event_times == (0.0, 20.0, 40.0, 60.0)


def test_dual_localized_pulse_has_jilkine_support_and_conserves_symmetry():
    grid = Grid1D(length=10.0, cells=100)
    pulse = LocalizedPulse(
        amplitude=0.2,
        stop=40.0,
        taper=20.0,
        sides=("left", "right"),
    )
    initial = pulse(10.0, grid.x, grid.length)
    tapered = pulse(30.0, grid.x, grid.length)
    assert np.allclose(initial, initial[::-1])
    assert np.max(initial) <= 0.2
    assert np.all(initial[(grid.x > 1.0) & (grid.x < 9.0)] == 0)
    assert np.max(tapered) < np.max(initial)
    assert np.all(pulse(40.0, grid.x, grid.length) == 0)


def test_localized_pulse_matches_jilkine_spatial_and_temporal_normalization():
    x = np.array([0.0, 0.5, 1.0, 5.0, 9.0, 9.5, 10.0])
    pulse = LocalizedPulse(
        amplitude=0.2,
        stop=40.0,
        taper=20.0,
        width_fraction=0.1,
        sides=("left", "right"),
        asymmetry=0.0,
    )
    full = pulse(20.0, x, 10.0)
    half_withdrawn = pulse(30.0, x, 10.0)
    assert np.allclose(full, [0.2, 0.1, 0.0, 0.0, 0.0, 0.1, 0.2])
    assert np.allclose(half_withdrawn, 0.5 * full)
    assert np.array_equal(full, full[::-1])


def test_dual_pulse_asymmetry_is_the_total_left_right_fractional_difference():
    pulse = LocalizedPulse(
        amplitude=0.2,
        stop=40.0,
        taper=20.0,
        sides=("left", "right"),
        asymmetry=0.01,
    )
    signal = pulse(0.0, np.array([0.0, 10.0]), 10.0)
    assert np.allclose(signal, [0.201, 0.199])


def test_reversing_gradient_changes_orientation():
    grid = Grid1D(length=10.0, cells=21)
    protocol = ReversingGradient(slope=0.02, reversal_time=10.0, transition=0.0)
    before = protocol(9.0, grid.x, grid.length)
    after = protocol(10.0, grid.x, grid.length)
    assert before[0] > before[-1]
    assert after[0] < after[-1]
    assert np.allclose(before, after[::-1])


def test_finite_reversing_gradient_controls_peak_and_spatial_width():
    grid = Grid1D(length=10.0, cells=21)
    protocol = ReversingFiniteGradient(
        amplitude=0.3, width_fraction=0.25, reversal_time=10.0, transition=0.0,
    )
    before = protocol(9.0, grid.x, grid.length)
    after = protocol(10.0, grid.x, grid.length)
    assert np.isclose(before[0], 0.3 * (1.0 - grid.x[0] / 2.5))
    assert np.all(before[grid.x >= 2.5] == 0)
    assert np.allclose(before, after[::-1])


def test_trajectory_round_trip(tmp_path):
    trajectory = Trajectory(
        model="test",
        time=np.array([0.0, 1.0]),
        x=np.array([0.25, 0.75]),
        fields={"active": np.array([[1.0, 1.0], [2.0, 1.0]])},
        observable_name="active",
        metadata={"backend": "unit"},
    )
    restored = Trajectory.load(trajectory.save(tmp_path / "trajectory.npz"))
    assert restored.model == trajectory.model
    assert restored.metadata == trajectory.metadata
    assert np.array_equal(restored.observable, trajectory.observable)
