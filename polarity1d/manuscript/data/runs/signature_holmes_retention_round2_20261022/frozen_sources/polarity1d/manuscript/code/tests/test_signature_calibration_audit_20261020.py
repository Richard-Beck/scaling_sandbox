"""Check independent calibration gates against contradictory saved outcomes."""
from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_calibration_audit_20261020 import audit_saved_tail, validate_duration_search


def example():
    case = dict(model="wave_pinning")
    config = dict(calibration_horizon=dict(wave_pinning=4.), calibration=dict(
        scans=dict(wave_pinning=[1., 4.]), absolute_tolerance=dict(wave_pinning=1.),
        relative_tolerance=0.))
    trials, evaluations = [], []
    for t, d in [(0., -1.), (1., .4), (2., .6), (4., .7)]:
        folder = f"workers/0000/evaluations/{t:g}"
        trials.append(dict(trial_directory=folder, final_d=d, readout_eligible=True))
        evaluations.append(dict(duration=t, trial_directory=folder, final_d=d,
            readout_eligible=True, replacement=d >= .5, field_tail_pass=False))
    calibration = dict(evaluations=evaluations, finite_duration_horizon=4.,
        scan_durations=[0., 1., 4.], observed_nonmonotone_commitment=False,
        lower_failing_duration=1., upper_committing_duration=2., bracket_resolved=True,
        train_ready=True, field_tail_used_to_define_fate=False,
        status="sampled_commitment_bracket_resolved")
    return case, dict(calibration=calibration), dict(trials=trials), config


def test_failed_field_tails_do_not_change_postwithdrawal_calibration_fate():
    case, result, native, config = example()
    receipt = validate_duration_search(case, result, native, config)
    assert receipt["passed"] and receipt["train_ready"]
    assert receipt["global_minimality_proven"] is False


def test_saved_fate_must_match_independently_reconstructed_eligible_activity():
    case, result, native, config = example()
    native["trials"][2]["readout_eligible"] = False
    result["calibration"]["evaluations"][2]["readout_eligible"] = False
    with pytest.raises(ValueError, match="sampled postwithdrawal fate"):
        validate_duration_search(case, result, native, config)


def test_omitted_scan_and_invented_bracket_are_rejected():
    case, result, native, config = example()
    bad = deepcopy(result)
    bad["calibration"]["lower_failing_duration"] = .5
    with pytest.raises(ValueError, match="lower failing bracket"):
        validate_duration_search(case, bad, native, config)
    config["calibration"]["scans"]["wave_pinning"].append(3.)
    with pytest.raises(ValueError, match="omitted a declared scan"):
        validate_duration_search(case, result, native, config)


def test_later_loss_of_commitment_prevents_train_ready_even_with_early_bracket():
    case, result, native, config = example()
    native["trials"][-1]["final_d"] = .2
    last = result["calibration"]["evaluations"][-1]
    last.update(final_d=.2, replacement=False)
    result["calibration"].update(observed_nonmonotone_commitment=True, train_ready=False)
    assert validate_duration_search(case, result, native, config)["train_ready"] is False


def test_full_field_tail_failure_is_recomputed_even_when_scalar_fate_is_constant(tmp_path):
    case = dict(model="wave_pinning", native_length=2., fixed_scales=dict(active=1.))
    folder = tmp_path/"workers/0000/evaluations/1"
    folder.mkdir(parents=True)
    np.savez(folder/"trace.npz", time=np.array([0., 1., 2.]),
             field_active=np.array([[1., 0.], [1., 0.], [1.001, 0.]]))
    evaluation = dict(duration=1., trial_directory="workers/0000/evaluations/1",
        field_tail_pass=False, field_tail=dict(available=True, window=2., tail_samples=3,
            tolerance=.001, field_tail_pass=False, field_tail_max_scaled=.002,
            field_tail_metrics=dict(active=dict(range_scaled=.001, projected_drift_500=.002))))
    receipt = audit_saved_tail(case, evaluation, tmp_path, None)
    assert receipt["passed"] and receipt["native_field_tail_pass"] is False
    evaluation["field_tail"]["field_tail_pass"] = True
    with pytest.raises(ValueError, match="full-field tail flag"):
        audit_saved_tail(case, evaluation, tmp_path, None)
