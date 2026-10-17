import numpy as np

from polarity1d import (
    Trajectory,
    mechanochemical_summary,
    phenotype_summary,
    profile_domain_count,
    tracking_summary,
)


def test_domain_counter_includes_both_boundaries():
    profile = np.array([2.0, 1.2, 0.2, 0.1, 0.2, 1.2, 2.0])
    assert profile_domain_count(profile) == 2


def test_common_summary_detects_resolution_and_winner():
    time = np.arange(0.0, 31.0, 5.0)
    x = (np.arange(21) + 0.5) / 21
    left = np.exp(-x / 0.1)
    right = np.exp(-(1 - x) / 0.1)
    fields = np.vstack(
        [
            np.zeros_like(x),
            left + right,
            left + right,
            1.2 * left + 0.5 * right,
            1.4 * left + 0.2 * right,
            1.5 * left,
            1.5 * left,
        ]
    )
    trajectory = Trajectory(
        model="synthetic",
        time=time,
        x=x,
        fields={"active": fields},
        observable_name="active",
    )
    summary = phenotype_summary(
        trajectory,
        analysis_start=5.0,
        polarization_threshold=0.1,
        sustained_duration=5.0,
    )
    assert summary["maximum_domain_count"] == 2
    assert summary["final_domain_count"] == 1
    assert summary["winning_side"] == "left"
    assert summary["competition_established"]
    assert summary["competition_status"] == "resolved"
    assert summary["competition_resolution_time_s"] == 25.0
    assert summary["competition_resolution_latency_s"] == 20.0
    assert summary["endpoint_state"] == "single_front"
    assert summary["competition_effective_status"] == "resolved_persistent"
    assert summary["competition_resolution_persistent_at_end"]
    assert not summary["competition_resolution_transient_only"]


def test_competition_distinguishes_transient_from_persistent_resolution():
    time = np.arange(0.0, 41.0, 5.0)
    x = (np.arange(21) + 0.5) / 21
    left = np.exp(-x / 0.1)
    right = np.exp(-(1 - x) / 0.1)
    fields = np.vstack(
        [
            left + right,
            left + right,
            1.5 * left,
            1.5 * left,
            np.zeros_like(x),
            np.zeros_like(x),
            np.zeros_like(x),
            np.zeros_like(x),
            np.zeros_like(x),
        ]
    )
    trajectory = Trajectory(
        model="synthetic_transient_resolution",
        time=time,
        x=x,
        fields={"active": fields},
        observable_name="active",
    )
    summary = phenotype_summary(
        trajectory,
        polarization_threshold=0.1,
        sustained_duration=5.0,
    )
    assert summary["competition_status"] == "resolved"
    assert summary["endpoint_state"] == "unpolarized"
    assert summary["competition_effective_status"] == "resolved_transient_to_unpolarized"
    assert not summary["competition_resolution_persistent_at_end"]
    assert summary["competition_resolution_transient_only"]


def test_single_front_is_not_misclassified_as_resolved_competition():
    time = np.arange(0.0, 31.0, 5.0)
    x = (np.arange(21) + 0.5) / 21
    front = np.exp(-x / 0.1)
    trajectory = Trajectory(
        model="synthetic",
        time=time,
        x=x,
        fields={"active": np.vstack([front for _ in time])},
        observable_name="active",
    )
    summary = phenotype_summary(
        trajectory,
        polarization_threshold=0.1,
        sustained_duration=5.0,
    )
    assert not summary["competition_established"]
    assert summary["competition_status"] == "not_established"
    assert summary["competition_resolution_time_s"] is None
    assert not summary["competition_resolution_censored"]


def test_reversal_reports_conditioning_decay_and_new_front_establishment():
    time = np.arange(0.0, 101.0, 5.0)
    x = (np.arange(41) + 0.5) / 41
    left = np.exp(-x / 0.1)
    right = np.exp(-(1 - x) / 0.1)
    profiles = []
    for at in time:
        if at < 50:
            profiles.append(left)
        elif at < 60:
            profiles.append(left + right)
        else:
            profiles.append(right)
    trajectory = Trajectory(
        model="synthetic",
        time=time,
        x=x,
        fields={"active": np.asarray(profiles)},
        observable_name="active",
    )
    summary = phenotype_summary(
        trajectory,
        polarization_threshold=0.1,
        sustained_duration=10.0,
        reversal_time=50.0,
    )
    assert summary["reversal_conditioned"]
    assert summary["pre_reversal_time_s"] == 45.0
    assert summary["pre_reversal_state"] == "single_front"
    assert summary["pre_reversal_polarized"]
    assert summary["pre_reversal_domain_count"] == 1
    assert summary["pre_reversal_orientation"] < -0.05
    assert summary["orientation_reversal_conditioned"]
    assert summary["orientation_reversal_success"]
    assert summary["orientation_reversal_latency_s"] == 10.0
    assert summary["reversal_success"]
    assert summary["reversal_status"] == "success"
    assert summary["new_front_establishment_time_s"] == 60.0
    assert summary["new_front_establishment_latency_s"] == 10.0


def test_orientation_reversal_accepts_two_peaks_on_the_same_side():
    time = np.arange(0.0, 101.0, 5.0)
    x = (np.arange(81) + 0.5) / 81
    left = np.exp(-x / 0.035) + 0.9 * np.exp(
        -0.5 * ((x - 0.18) / 0.035) ** 2
    )
    right = left[::-1]
    profiles = np.asarray([left if at < 60 else right for at in time])
    trajectory = Trajectory(
        model="synthetic_clustered_front",
        time=time,
        x=x,
        fields={"active": profiles},
        observable_name="active",
    )
    summary = phenotype_summary(
        trajectory,
        polarization_threshold=0.1,
        sustained_duration=10.0,
        reversal_time=50.0,
    )
    assert summary["pre_reversal_state"] == "multifront"
    assert summary["pre_reversal_domain_count"] == 2
    assert not summary["reversal_conditioned"]
    assert summary["reversal_status"] == "not_conditioned"
    assert summary["orientation_reversal_conditioned"]
    assert summary["orientation_reversal_success"]
    assert summary["orientation_reversal_latency_s"] == 10.0


def test_tracking_summary_uses_dimensionless_orientation_and_positive_delay():
    time = np.arange(0.0, 401.0, 5.0)
    x = (np.arange(41) + 0.5) / 41
    centered = 2.0 * x - 1.0
    period = 100.0
    delay = 12.5
    profiles = np.asarray(
        [
            1.0
            + 0.3
            * (-np.cos(2.0 * np.pi * (at - delay) / period))
            * centered
            for at in time
        ]
    )
    trajectory = Trajectory(
        model="synthetic",
        time=time,
        x=x,
        fields={"active": profiles},
        observable_name="active",
    )
    summary = tracking_summary(trajectory, period, start=0.0)
    assert summary["tracking_reliable"]
    assert summary["tracking_correlation"] > 0.6
    assert abs(summary["phase_lag_seconds"] - delay) < 2.0


def test_mechanochemical_summary_retains_signed_distal_response():
    time = np.array([0.0, 5.0, 10.0, 15.0])
    x = (np.arange(10) + 0.5) / 10
    rac = np.ones((4, 10))
    rho = np.ones((4, 10))
    rac[2:, :2] += 2.0
    rac[2:, -2:] += 1.0
    rho[2:, :2] -= 0.5
    rho[2:, -2:] -= 0.25
    trajectory = Trajectory(
        model="synthetic_mechanical",
        time=time,
        x=x,
        fields={
            "rac": rac,
            "rho": rho,
            "tension": np.tile(np.arange(4)[:, None], (1, 10)),
        },
        observable_name="rac",
    )
    summary = mechanochemical_summary(trajectory, stimulus_start=10.0)
    assert summary["rac_distal_peak_change"] == 1.0
    assert summary["rho_distal_peak_change"] == -0.25
    assert summary["rac_distal_to_proximal_abs_ratio"] == 0.5
    assert summary["maximum_absolute_tension"] == 3.0
