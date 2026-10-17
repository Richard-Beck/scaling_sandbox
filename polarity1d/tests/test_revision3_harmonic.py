import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from reports.four_test_suite.revision3.harmonic import edge_orientation, score


def test_known_delayed_harmonic_has_positive_lag_and_is_reliable():
    period = 100.0
    time = np.arange(0.0, 6.0 * period + 0.625, 0.625)
    delay = 10.0
    orientation = -0.35 * np.cos(2.0 * np.pi * (time - delay) / period)
    result = score(time, orientation, period)

    assert result["cycle_count"] == 5
    assert abs(result["amplitude"] - 0.35) < 0.003
    # The archived inclusive endpoint convention adds one sample at 6P.
    assert abs(result["phase_fraction"] - delay / period) < 3e-4
    assert result["correlation"] > 0
    assert result["r_squared"] > 0.999
    assert result["reliable"] is True
    assert result["failure_reasons"] == []
    assert result["scoring_window"] == {"start": period, "end": 6 * period, "inclusive": True}


def test_edge_orientation_uses_outer_twenty_percent_means():
    x = np.linspace(0.0, 10.0, 11)
    active = np.tile(np.arange(11.0), (2, 1))
    result = edge_orientation(x, active)
    left = np.mean(np.arange(11.0)[x / x[-1] <= 0.2])
    right = np.mean(np.arange(11.0)[x / x[-1] >= 0.8])
    assert np.allclose(result, (right - left) / (abs(right) + abs(left)))


def test_values_outside_scoring_window_do_not_affect_fit():
    period = 100.0
    time = np.arange(0.0, 7.0 * period + 0.625, 0.625)
    orientation = -0.3 * np.cos(2.0 * np.pi * (time - 8.0) / period)
    baseline = score(time, orientation, period)
    altered = orientation.copy()
    altered[time < period] = 1000.0
    altered[time > 6.0 * period] = -1000.0
    result = score(time, altered, period)
    for key in ("amplitude", "phase_fraction", "correlation", "r_squared",
                "cycle_amplitude_cv", "cycle_phase_dispersion", "reliable"):
        assert result[key] == baseline[key]


def test_archived_audit_condition_reproduces_saved_harmonic_metrics():
    archive = ROOT / "polarity1d/results/campaigns/wp1_tracking_reliability_audit_v1"
    with (archive / "tracking_reliability_audit.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    reference = next(row for row in rows if row["model"] == "jilkine_wave_pinning"
                     and float(row["length_um"]) == 5.0 and float(row["period_s"]) == 50.0)
    task_path = next((archive / "tasks").glob("*/task.json"))
    selected = None
    for candidate in (archive / "tasks").glob("*/task.json"):
        task = json.loads(candidate.read_text())
        if (task["model"] == "wave_pinning" and task["simulation"]["length"] == 5.0
                and task["tracking"]["period"] == 50.0):
            selected = (candidate.parent, task)
            break
    assert selected is not None, f"archived trajectory missing; searched from {task_path}"
    directory, task = selected
    with np.load(directory / "trajectory.npz") as data:
        time = data["time"]
        orientation = edge_orientation(data["x"], data["field__active"])
    result = score(time, orientation, task["tracking"]["period"])

    assert abs(result["amplitude"] - float(reference["audit_amplitude"])) < 1e-12
    assert abs(result["phase_fraction"] - float(reference["audit_phase_fraction"])) < 1e-12
    assert abs(result["correlation"] - float(reference["audit_correlation"])) < 1e-12
    assert abs(result["r_squared"] - float(reference["harmonic_r_squared"])) < 1e-12
    assert result["reliable"] == (reference["audited_reliable"] == "True")
