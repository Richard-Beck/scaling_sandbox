"""Independent saved-field and duration-bracket audit of commitment calibration.

No integration runs here. Frozen native phases and independently reconstructed
readouts establish each sampled fate; full-field tail flags remain separate.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path

import numpy as np
import signature_analysis_20261020 as analysis


def _equal(actual, expected, label):
    if actual is None or expected is None or isinstance(expected, (str, bool)):
        if actual != expected:
            raise ValueError(f"Calibration audit mismatch at {label}: {actual!r} != {expected!r}")
    elif not math.isclose(float(actual), float(expected), rel_tol=2e-10, abs_tol=2e-11):
        raise ValueError(f"Calibration audit mismatch at {label}: {actual!r} != {expected!r}")


def validate_duration_search(case, result, reconstructed, config):
    """Check native sampled fates and bracket arithmetic without rerunning search."""
    calibration = result["calibration"]
    evaluations = calibration["evaluations"]
    by_folder = {trial["trial_directory"]: trial for trial in reconstructed["trials"]}
    if len(by_folder) != len(evaluations):
        raise ValueError("Calibration native evaluation count differs from saved search")
    records = {}
    for evaluation in evaluations:
        duration = float(evaluation["duration"])
        if duration in records:
            raise ValueError("Repeated duration in saved calibration evaluations")
        native = by_folder[evaluation["trial_directory"]]
        _equal(evaluation["final_d"], native["final_d"], "native final D")
        _equal(evaluation["readout_eligible"], native["readout_eligible"], "native eligibility")
        replacement = bool(native["readout_eligible"] and native["final_d"] is not None
                           and native["final_d"] >= .5)
        _equal(evaluation["replacement"], replacement, "sampled postwithdrawal fate")
        records[duration] = replacement
    horizon = float(config["calibration_horizon"][case["model"]])
    scan = sorted({0., horizon, *map(float, config["calibration"]["scans"][case["model"]])})
    if any(t < 0 or t > horizon for t in records) or not set(scan).issubset(records):
        raise ValueError("Calibration omitted a declared scan or exceeded its horizon")
    _equal(calibration["finite_duration_horizon"], horizon, "search horizon")
    if list(calibration["scan_durations"]) != scan:
        raise ValueError("Recorded calibration scan differs from frozen configuration")
    if records[0.]:
        raise ValueError("Zero-duration control already replaces; incumbent is invalid")
    times = sorted(records)
    nonmonotone = any(records[a] and not records[b] for a, b in zip(times, times[1:]))
    _equal(calibration["observed_nonmonotone_commitment"], nonmonotone, "sampled monotonicity")
    successes = [t for t in times if records[t]]
    lower = upper = None
    resolved = False
    if successes:
        upper = min(successes)
        lower = max(t for t in times if t < upper and not records[t])
        settings = config["calibration"]
        tolerance = max(settings["absolute_tolerance"][case["model"]],
                        settings["relative_tolerance"] * upper)
        resolved = upper-lower <= tolerance
        if case["model"] == "debelly":
            if not all(math.isclose(t/.5, round(t/.5), abs_tol=1e-9) for t in times):
                raise ValueError("deBelly calibration left its native half-unit lattice")
            resolved = resolved or round((upper-lower)/.5) <= 1
    _equal(calibration["lower_failing_duration"], lower, "lower failing bracket")
    _equal(calibration["upper_committing_duration"], upper, "upper committing bracket")
    _equal(calibration["bracket_resolved"], resolved, "duration resolution")
    _equal(calibration["train_ready"], bool(resolved and not nonmonotone), "train readiness")
    _equal(calibration["field_tail_used_to_define_fate"], False, "independent field-tail flag")
    expected_status = ("no_sampled_commitment_through_finite_horizon" if not successes else
                       "sampled_commitment_bracket_resolved" if resolved else
                       "sampled_commitment_bracket_iteration_limit")
    _equal(calibration["status"], expected_status, "search outcome")
    return dict(passed=True, evaluations=len(records), finite_duration_horizon=horizon,
                lower_failing_duration=lower, upper_committing_duration=upper,
                train_ready=calibration["train_ready"],
                observed_nonmonotone_commitment=nonmonotone,
                global_minimality_proven=False, field_tail_defines_fate=False)


def _field(archive, name, model, length):
    if name == "physical_length" and model == "debelly":
        values = length + archive["field_boundary_separation"]
    elif "field_" + name in archive:
        values = archive["field_" + name]
    elif name == "length":
        values = archive["physical_length"]
    else:
        raise ValueError("Missing monitored native field: " + name)
    if name == "displacement" and model == "debelly":
        values = values-values.mean(axis=1, keepdims=True)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite monitored native field: " + name)
    return values


def audit_saved_tail(case, evaluation, root, basal):
    """Recompute finite-tail range/drift and zero-cue excursions from full fields."""
    path = root / evaluation["trial_directory"] / "trace.npz"
    tail = evaluation["field_tail"]
    if tail.get("available") is not True:
        raise ValueError("Calibration recovery lacks its declared full-field tail")
    with np.load(path) as trace:
        time = trace["time"]
        keep = time >= time[-1]-tail["window"]-1e-10
        tt = time[keep]
        if len(tt) < 2 or len(tt) != tail["tail_samples"]:
            raise ValueError("Calibration tail sample count mismatch")
        maxima, excursions = [], {}
        metrics = tail["field_tail_metrics"]
        if set(metrics) != set(case["fixed_scales"]):
            raise ValueError("Tail does not cover every declared monitored native field")
        for name, values in metrics.items():
            field = _field(trace, name, case["model"], case["native_length"])
            scale = float(case["fixed_scales"][name])
            yy = field[keep]
            variation = float(np.max(np.max(yy, axis=0)-np.min(yy, axis=0))/scale)
            drift = float(np.max(np.abs(yy[-1]-yy[-2]))/(tt[-1]-tt[-2])*tail["window"]/scale)
            _equal(values["range_scaled"], variation, name+" tail range")
            _equal(values["projected_drift_500"], drift, name+" projected tail drift")
            maxima.extend([variation, drift])
            if evaluation["duration"] == 0:
                reference = _field(basal, name, case["model"], case["native_length"])[-1]
                excursions[name] = float(np.max(np.abs(field-reference))/scale)
        maximum = max(maxima)
        _equal(tail["field_tail_max_scaled"], maximum, "full-field tail maximum")
        _equal(tail["field_tail_pass"], maximum <= tail["tolerance"], "full-field tail flag")
        _equal(evaluation["field_tail_pass"], tail["field_tail_pass"], "evaluation tail flag")
        if excursions:
            control = evaluation["no_cue_control"]
            for name, value in excursions.items():
                _equal(control["no_cue_scaled_excursions"][name], value, name+" zero-cue excursion")
            _equal(control["no_cue_control_pass"], max(excursions.values()) <= control["tolerance"],
                   "zero-cue field stability")
            if control["passed"] is not True or control["spontaneous_D_replacement"] is not False:
                raise ValueError("Calibration zero-cue control failed")
    return dict(duration=evaluation["duration"], native_field_tail_pass=tail["field_tail_pass"],
                maximum_scaled_tail=maximum, monitored_fields=sorted(metrics),
                zero_cue_excursions=excursions, passed=True)


def audit_worker(root, task, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if not output.is_relative_to(root/"results"):
        raise ValueError("Audit outputs must be separate from original native workers")
    manifest = analysis._read(root/"manifest.json")
    case = next(c for c in manifest["cases"] if c["task"] == task)
    if case["stage"] != "calibration":
        raise ValueError("Expected a frozen calibration worker")
    result_path = root/"workers"/f"{task:04d}"/"result.json"
    result = analysis._read(result_path)
    if not all(e.get("events") and all("archive" in event for event in e["events"])
               for e in result["calibration"]["evaluations"]):
        raise ValueError("Every calibration phase must declare its full native archive")
    config = analysis._read(root/"config.json")
    if analysis._sha(root/"config.json") != manifest["config_sha256"]:
        raise ValueError("Frozen calibration configuration hash mismatch")
    reconstructed = analysis.analyze_worker_output(case, result, run_root=root, result_path=result_path)
    if not all(trial["native_phase_join_audit"].get("passed") is True
               for trial in reconstructed["trials"]):
        raise ValueError("Calibration native state phase continuity is not independently verified")
    search = validate_duration_search(case, result, reconstructed, config)
    with np.load(root/result["basal"]["checkpoint"]) as basal:
        tails = [audit_saved_tail(case, evaluation, root, basal)
                 for evaluation in result["calibration"]["evaluations"]]
    receipt = dict(generated_utc=datetime.now(timezone.utc).isoformat(), task=task,
        model=case["model"], reference_length=case["requested_reference_length"], passed=True,
        run_root=str(root), manifest_sha256=analysis._sha(root/"manifest.json"),
        result_sha256=analysis._sha(result_path), audit_source_sha256=analysis._sha(Path(__file__)),
        analysis_source_sha256=analysis._sha(Path(analysis.__file__)),
        search=search, tails=tails, native_trial_count=len(reconstructed["trials"]),
        numeric_readouts_reconstructed=True, complete_state_phase_continuity=True,
        field_tail_defines_fate=False)
    analysis._write_json(output/f"worker_{task:04d}.json", receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--task", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Full native-field audit requires SLURM")
    receipt = audit_worker(args.run_root, args.task, args.output)
    print(json.dumps({k:receipt[k] for k in ("task", "passed", "native_trial_count", "search")}))


if __name__ == "__main__":
    main()
