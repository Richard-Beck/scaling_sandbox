"""Audit frozen signature workers and extract finite, sampled assay outcomes.

No solver or figure renderer is called. Continuous cue-ON events, recovery
events, final fate and field-tail qualification remain separate. Missing work
never becomes finite censoring; observed low-activity states remain known
nonreplacement under the eligibility gate. Ratios propagate sampled brackets.
"""
from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from datetime import datetime, timezone

import numpy as np


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path):
    return json.loads(Path(path).read_text())


def _inside(root, name):
    path = Path(name)
    path = path.resolve() if path.is_absolute() else (root / path).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"Scientific artifact lies outside frozen run root: {name}")
    return path


def _finite(value):
    return value is not None and not isinstance(value, bool) and math.isfinite(float(value))


def threshold_event(series, start, stop, threshold, *, value_key="D_half",
                    eligibility_key="D_half_eligible", time_tolerance=1e-9):
    """First eligible sampled threshold hit within one exact phase interval.

    Returns absolute and phase-elapsed brackets. Explicit failed eligibility is
    a known-false replacement predicate, even when D is undefined. Unexplained
    missing D or unavailable native samples remain unscorable. Interpolation is
    omitted when the lower bracket has known failed eligibility and no D value.
    """
    start, stop, threshold = float(start), float(stop), float(threshold)
    if not all(math.isfinite(v) for v in (start, stop, threshold)) or stop < start:
        raise ValueError("Finite ordered interval and threshold required")
    time = np.asarray(series["time"], dtype=float)
    values = series[value_key]
    eligibility = series.get(eligibility_key, [True] * len(time))
    available = series.get("native_sample_available", [True] * len(time))
    if (time.ndim != 1 or not len(time) or not np.isfinite(time).all()
            or np.any(np.diff(time) <= 0)
            or len(values) != len(time) or len(eligibility) != len(time)
            or len(available) != len(time)):
        raise ValueError("Aligned scalar series and strictly increasing finite time required")
    keep = np.flatnonzero((time >= start - time_tolerance) & (time <= stop + time_tolerance))
    boundary_complete = bool(len(keep) and abs(time[keep[0]] - start) <= time_tolerance
                             and abs(time[keep[-1]] - stop) <= time_tolerance)
    valid = [bool(available[i]) and bool(eligibility[i]) and _finite(values[i]) for i in keep]
    known_false = [bool(available[i]) and isinstance(eligibility[i], (bool, np.bool_))
                   and not bool(eligibility[i]) for i in keep]
    known = [a or b for a, b in zip(valid, known_false)]
    base = dict(threshold=threshold, value_key=value_key, interval=[start, stop],
                samples=len(keep), eligible_samples=sum(valid),
                ineligible_samples=sum(known_false), missing_readout_samples=len(known) - sum(known),
                exact_boundary_coverage=boundary_complete,
                lower_time=None, upper_time=None, lower_elapsed=None, upper_elapsed=None,
                estimated_elapsed=None, right_censored_at=None,
                lower_D_unavailable=False,
                definition="First eligible sampled threshold hit; interpolation is descriptive.")
    if not boundary_complete:
        return dict(base, status="incomplete_window", observed=False)
    hits = [j for j, i in enumerate(keep) if valid[j] and float(values[i]) >= threshold]
    if not hits:
        if not all(known):
            return dict(base, status="unscorable_samples", observed=False)
        return dict(base, status="right_censored", observed=False,
                    right_censored_at=stop - start)
    j = hits[0]
    i = int(keep[j])
    upper = float(time[i])
    if j == 0:
        return dict(base, status="already_above_at_interval_start", observed=True,
                    upper_time=upper, upper_elapsed=upper - start)
    if not all(known[:j]):
        return dict(base, status="observed_after_unscorable_samples", observed=True,
                    upper_time=upper, upper_elapsed=upper - start)
    previous = int(keep[j - 1])
    lower = float(time[previous])
    estimate = None
    if valid[j - 1]:
        fraction = (threshold - float(values[previous])) / (float(values[i]) - float(values[previous]))
        estimate = lower + fraction * (upper - lower) - start
    return dict(base, status="observed", observed=True, lower_time=lower,
                upper_time=upper, lower_elapsed=lower - start,
                upper_elapsed=upper - start, estimated_elapsed=estimate,
                lower_D_unavailable=not valid[j - 1])


def replacement_time_ratio(narrow, broad):
    """Bound narrow/broad event time using sampled intervals and finite censoring."""
    allowed = {"observed", "right_censored"}
    if narrow.get("status") not in allowed or broad.get("status") not in allowed:
        return dict(status="unavailable", lower=None, upper=None,
                    reason="Both inputs require an observed sampled bracket or complete finite censoring.")

    def limits(event):
        if event["status"] == "right_censored":
            value = event.get("right_censored_at")
            if not _finite(value) or float(value) <= 0:
                raise ValueError("Finite positive censoring horizon required")
            return float(value), math.inf, True
        lo, hi = event.get("lower_elapsed"), event.get("upper_elapsed")
        if not _finite(lo) or not _finite(hi) or not 0 <= float(lo) < float(hi):
            raise ValueError("Observed time requires an ordered nonnegative sampled bracket")
        return float(lo), float(hi), False

    nlo, nhi, nc = limits(narrow)
    blo, bhi, bc = limits(broad)
    lower = 0. if math.isinf(bhi) else nlo / bhi
    upper = math.inf if blo <= 0 or math.isinf(nhi) else nhi / blo
    estimate = None
    if (not nc and not bc and _finite(narrow.get("estimated_elapsed"))
            and _finite(broad.get("estimated_elapsed")) and broad["estimated_elapsed"] > 0):
        estimate = narrow["estimated_elapsed"] / broad["estimated_elapsed"]
    return dict(status="bounded" if math.isfinite(upper) else "lower_bound_only",
                lower=float(lower), upper=None if math.isinf(upper) else float(upper),
                lower_strict=bool(nc and not bc), upper_strict=bool(bc and not nc),
                estimate_from_interpolation=estimate,
                narrow_censored=nc, broad_censored=bc,
                informative=bool(lower > 0 or math.isfinite(upper)),
                interpretation="Sampled-event/censoring bounds; no finite event time assigned to censoring.")


def _endpoint(series, when, tolerance=1e-9):
    time = np.asarray(series["time"], dtype=float)
    matches = np.flatnonzero(np.abs(time - when) <= tolerance)
    if len(matches) != 1:
        return dict(time=float(when), status="missing_exact_endpoint")
    index = int(matches[0])
    result = {key: values[index] for key, values in series.items()
              if isinstance(values, list) and len(values) == len(time)}
    result.update(time=float(when), status="observed_exact_endpoint")
    return result


def pulse_gap_endpoints(series, events, *, gain_tolerance=1e-12):
    """Extract all exact endpoints and continuous gained-share gap retention."""
    initial = _endpoint(series, float(series["time"][0]))
    basal_share = initial.get("share_half")
    rows, last_on = [], None
    expanded = []
    for event in events:
        is_on = event["name"] == "cue_on" or str(event["name"]).endswith("on")
        if expanded and is_on:
            previous = expanded[-1]
            previous_on = previous["name"] == "cue_on" or str(previous["name"]).endswith("on")
            if previous_on and math.isclose(previous["stop"], event["start"], rel_tol=0, abs_tol=1e-9):
                expanded.append(dict(name="implicit_zero_gap", start=event["start"],
                                     stop=event["start"], zero_duration=True))
        expanded.append(event)
    for event in expanded:
        name = str(event["name"])
        kind = ("recovery" if name == "recovery" else
                "gap" if name.endswith("gap") else
                "on" if name == "cue_on" or name.endswith("on") else "other")
        start = _endpoint(series, float(event["start"]))
        end = _endpoint(series, float(event["stop"]))
        row = dict(kind=kind, phase_name=name, phase_index=event.get("index"),
                   start=start, end=end, retention_ratio=None,
                   retention_eligible=False, retention_reason=None,
                   zero_duration=bool(event.get("zero_duration", False)),
                   retention_evidence_state=None,
                   retention_probability_denominator_eligible=False,
                   observed_retention_loss=False,
                   primary_retention_pulse_eligible=False,
                   primary_retention_denominator_eligible=False,
                   primary_retention_ratio=None,
                   primary_retention_success=None,
                   primary_retention_success_by_threshold={"0.7": None, "0.8": None, "0.9": None})
        if kind == "gap":
            on_share = None if last_on is None else last_on["end"].get("share_half")
            gap_share = end.get("share_half")
            pulse = {} if last_on is None else last_on["end"]
            pulse_D, gap_D = pulse.get("D_half"), end.get("D_half")
            pulse_eligible = bool(pulse.get("readout_eligible", pulse.get("D_half_eligible", True)))
            pulse_ok = (pulse_eligible and _finite(pulse_D) and pulse_D >= 0
                        and _finite(on_share) and on_share > 0)
            gap_eligibility = end.get("readout_eligible", end.get("D_half_eligible", True))
            gap_known_ineligible = gap_eligibility is False
            gap_known = (end["status"] == "observed_exact_endpoint"
                         and end.get("native_sample_available", True)
                         and (gap_known_ineligible or (_finite(gap_D) and _finite(gap_share))))
            row["primary_retention_pulse_eligible"] = bool(pulse_ok)
            row["primary_retention_denominator_eligible"] = bool(pulse_ok and gap_known)
            if pulse_ok and gap_known:
                ratio = gap_share / on_share if not gap_known_ineligible else None
                row["primary_retention_ratio"] = ratio
                for threshold in (.7, .8, .9):
                    row["primary_retention_success_by_threshold"][str(threshold)] = bool(
                        not gap_known_ineligible and gap_D >= 0 and ratio >= threshold)
                row["primary_retention_success"] = row["primary_retention_success_by_threshold"]["0.8"]
            if not all(_finite(v) for v in (basal_share, on_share)):
                row["retention_reason"] = "missing_or_ineligible_baseline_or_pulse_share"
                row["retention_evidence_state"] = "uninformative_pulse"
            elif float(on_share) - float(basal_share) <= gain_tolerance:
                row["retention_reason"] = "pulse_did_not_gain_resolvable_challenger_share"
                row["retention_evidence_state"] = "uninformative_pulse"
            elif end["status"] != "observed_exact_endpoint":
                row["retention_reason"] = "missing_exact_gap_endpoint"
                row["retention_evidence_state"] = "missing"
            elif end.get("readout_eligible", end.get("D_half_eligible")) is False:
                row["retention_reason"] = "observed_gap_activity_below_eligibility_gate"
                row["retention_evidence_state"] = "observed_activity_loss"
                row["retention_probability_denominator_eligible"] = True
                row["observed_retention_loss"] = True
            elif not _finite(gap_share):
                row["retention_reason"] = "missing_gap_share_readout"
                row["retention_evidence_state"] = "missing"
            else:
                row["retention_ratio"] = ((gap_share - basal_share) / (on_share - basal_share))
                row["retention_eligible"] = True
                row["retention_probability_denominator_eligible"] = True
                row["retention_evidence_state"] = "observed_numeric"
                row["challenger_amount_change"] = (
                    end.get("E_challenger_half", 0.) - last_on["end"].get("E_challenger_half", 0.))
                row["D_change"] = end["D_half"] - last_on["end"]["D_half"]
        rows.append(row)
        if kind == "on":
            last_on = row
    return rows


def _trace_mapping(path, model):
    with np.load(path, allow_pickle=False) as data:
        trace = dict(time=data["time"].copy(), x=data["x"].copy(), active=data["active"].copy(),
                     sample_kind="nodes" if model == "spring" else "cells")
        if model == "spring":
            length = data["field_length"] if "field_length" in data else data["physical_length"]
            trace["physical_x"] = np.asarray(length)[:, None] * np.linspace(0, 1, len(trace["x"]))
        else:
            x = trace["x"]
            dx = x[1] - x[0]
            trace["reference_edges"] = np.r_[x[0] - .5 * dx, x + .5 * dx]
    return trace


def _readout_parameters(case, basal):
    length = float(case["native_length"])
    physical_length = (float(basal["physical_x"][-1, -1] - basal["physical_x"][-1, 0])
                       if case["model"] == "spring" else length)
    fixed_width = 15.2 / 38. if case["model"] == "spring" else 15.2
    own_ref = float(case["cue_width"]) * length / physical_length
    fixed_ref = fixed_width * length / physical_length
    if case.get("assay") == "position":
        return dict(challenger_own_intervals=[[case["cue_center"] - case["cue_width"] / 2,
                                              case["cue_center"] + case["cue_width"] / 2]],
                    challenger_fixed_intervals=[[max(0., length - fixed_ref), length]],
                    challenger_side="right", denominator_atol=1e-12, contrast_atol=.05)
    return dict(challenger_own_intervals=[(max(0., length - own_ref), length)],
                challenger_fixed_intervals=[(max(0., length - fixed_ref), length)],
                challenger_own_width_physical=float(case["cue_width"]),
                challenger_side="right", denominator_atol=1e-12, contrast_atol=.05)


def _same_numeric_tree(first, second, path="root"):
    if isinstance(first, dict) and isinstance(second, dict):
        if set(first) != set(second):
            raise ValueError(f"Recomputed readout keys differ at {path}")
        for key in first:
            _same_numeric_tree(first[key], second[key], path + "." + key)
    elif isinstance(first, (list, tuple)) and isinstance(second, (list, tuple)):
        if len(first) != len(second):
            raise ValueError(f"Recomputed readout length differs at {path}")
        for i, (a, b) in enumerate(zip(first, second)):
            _same_numeric_tree(a, b, f"{path}[{i}]")
    elif isinstance(first, bool) or isinstance(second, bool) or first is None or second is None:
        if type(first) is not type(second) or first != second:
            raise ValueError(f"Recomputed readout value differs at {path}")
    elif isinstance(first, (int, float)) and isinstance(second, (int, float)):
        if not math.isclose(first, second, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError(f"Recomputed numeric readout differs at {path}")
    elif first != second:
        raise ValueError(f"Recomputed readout value differs at {path}")


def _regional_peak(trace, intervals_by_frame):
    values = []
    x = trace["x"]
    for i, active in enumerate(trace["active"]):
        excess = active - np.min(active)
        maxima = []
        for lo, hi in intervals_by_frame[i]:
            if hi <= lo:
                continue
            if trace["sample_kind"] == "cells":
                edges = trace["reference_edges"]
                mask = (edges[:-1] < hi) & (edges[1:] > lo)
                if np.any(mask):
                    maxima.append(float(np.max(excess[mask])))
            else:
                inside = excess[(x > lo) & (x < hi)]
                maxima.extend(np.r_[np.interp([lo, hi], x, excess), inside].tolist())
        values.append(max(maxima) if maxima else None)
    return values


def _trial_items(result):
    for name in ("zero_control", "challenge", "trial"):
        if name in result:
            yield name, result[name], name if name != "zero_control" else "zero"
    for evaluation in result.get("calibration", {}).get("evaluations", []):
        label = format(evaluation["duration"], ".12g").replace(".", "p")
        yield "calibration_" + label, evaluation, "evaluations/" + label


def declared_state_fields(model, root, manifest):
    """Read frozen native declarations without importing or running a solver."""
    if model == "spring":
        name = "polarity1d/manuscript/code/mechanisms/signature_native_20261020.py"
        # The frozen adapter's state_at constructs interleaved active/inactive
        # followed by both raw endpoint coordinates; length is derived.
        return ("field_active", "field_inactive", "raw_boundaries"), dict(
            source=name, sha256=manifest["source_sha256"][name], declaration="state_at spring branch")
    if model == "debelly":
        suffix, declaration = "code/core/native.py", "DB_NAMES"
    elif model in ("wave_pinning", "goryachev"):
        suffix, declaration = f"code/core/models/{model}.py", "state_names"
    elif model == "holmes_model3":
        suffix, declaration = "code/core/src/polarity1d/holmes.py", "state_names"
    else:
        raise ValueError("Unknown native-state declaration: " + model)
    matches = [name for name in manifest["source_sha256"] if name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError("Unique frozen native-state declaration required")
    name = matches[0]
    tree = ast.parse((root / "frozen_sources" / name).read_text())
    values = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == declaration for t in node.targets):
            values.append(tuple(ast.literal_eval(node.value)))
    if len(values) != 1 or not values[0] or not all(isinstance(v, str) for v in values[0]):
        raise ValueError("Native-state declaration is not one literal field list")
    return tuple("field_" + field for field in values[0]), dict(
        source=name, sha256=manifest["source_sha256"][name], declaration=declaration)


def audit_native_phase_join(trace_path, events, active_trace, *, state_fields=None):
    """Compare the retained full native fields with their exact phase archives.

    Hashes establish artifact identity; this independent check also establishes
    that the joined trace and every phase endpoint contain the same native
    samples. Legacy fixtures lacking phase archives explicitly lack this gate.
    """
    if not all("archive" in event for event in events):
        return dict(available=False, passed=None, reason="phase_archives_not_declared")
    with np.load(trace_path, allow_pickle=False) as joined:
        fields = [key for key in joined if key.startswith("field_")]
        if "raw_boundaries" in joined:
            fields.append("raw_boundaries")
        if not fields:
            raise ValueError("Declared native phases require full field histories")
        values = {key: joined[key].copy() for key in fields}
        values["active"] = active_trace["active"]
        state_fields = set(fields if state_fields is None else state_fields) | {"active"}
        if state_fields - set(values):
            raise ValueError("Declared native state fields absent from full trace")
        time = active_trace["time"]
        previous = None
        covered = np.zeros(len(time), dtype=bool)
        checks = []
        for event in events:
            path = Path(event["archive"]["path"])
            if not path.is_absolute():
                # Worker events currently carry absolute archive paths. Relative
                # paths are interpreted by the caller before this function.
                raise ValueError("Phase audit requires resolved archive paths")
            with np.load(path, allow_pickle=False) as phase:
                native_time = phase["time"] + float(event["start"])
                if (abs(native_time[0] - event["start"]) > 1e-9 or
                        abs(native_time[-1] - event["stop"]) > 1e-9):
                    raise ValueError("Native phase does not span its declared exact interval")
                index = np.searchsorted(time, native_time)
                if (np.any(index >= len(time)) or
                        not np.allclose(time[index], native_time, atol=1e-9, rtol=0)):
                    raise ValueError("Native phase samples are missing from joined trace")
                covered[index] = True
                endpoints, diagnostic_jumps = {}, {}
                for key, full in values.items():
                    if key not in phase:
                        raise ValueError("Full native field absent from phase: " + key)
                    samples = phase[key]
                    diagnostic = key not in state_fields
                    retained = slice(1, None) if previous is not None and diagnostic else slice(None)
                    tolerance = 32 * np.finfo(float).eps * np.maximum(1., np.abs(samples[retained]))
                    if not np.all(np.abs(full[index[retained]] - samples[retained]) <= tolerance):
                        raise ValueError("Joined trace differs from native phase: " + key)
                    if previous is not None and key in previous and not diagnostic:
                        tolerance = 32 * np.finfo(float).eps * np.maximum(1., np.abs(samples[0]))
                        if not np.all(np.abs(previous[key] - samples[0]) <= tolerance):
                            raise ValueError("Native field discontinuity at phase boundary: " + key)
                    elif previous is not None and diagnostic:
                        diagnostic_jumps[key] = dict(
                            preceding_phase_endpoint_min=float(np.min(previous[key])),
                            preceding_phase_endpoint_max=float(np.max(previous[key])),
                            current_phase_start_min=float(np.min(samples[0])),
                            current_phase_start_max=float(np.max(samples[0])),
                            maximum_absolute_one_sided_difference=float(np.max(np.abs(previous[key] - samples[0]))),
                            joined_trace_keeps="preceding_phase_endpoint",
                            exact_one_sided_arrays_retained_in_native_phase_archives=True)
                    endpoints[key] = samples[-1].copy()
                previous = endpoints
                checks.append(dict(index=event.get("index"), interval=[event["start"], event["stop"]],
                                   samples=len(index), native_fields=fields, passed=True,
                                   one_sided_diagnostic_boundaries=diagnostic_jumps))
        if not covered.all():
            raise ValueError("Joined trace contains samples outside declared native phases")
    return dict(available=True, passed=True, phases=checks,
                all_joined_samples_covered=True, roundoff_tolerance_binary64_eps=32,
                declared_state_fields=sorted(state_fields),
                diagnostic_boundary_convention="Joined trace retains preceding endpoint; later-phase duplicate diagnostic sample remains in its own phase archive.")


def analyze_worker_output(case, result, *, run_root, result_path=None,
                          verify_hashes=True, recompute_readouts=True, shared_hash_cache=None):
    """Audit and extract a completed worker, accepting legacy and schema-2 output.

    ``result`` may be its JSON path or an already loaded mapping. Missing files
    raise errors, rather than producing scientific censoring. Schema 2 requires
    declared complete scientific artifact hash coverage. Legacy readout JSONs
    without original digests are numerically recomputed from hashed native data
    using the frozen, hash-verified readout source; identity is never invented.
    """
    root = Path(run_root).resolve()
    if isinstance(result, (str, Path)):
        result_path = Path(result)
        result = _read(result_path)
    case = dict(case)
    worker_dir = root / "workers" / f"{int(case['task']):04d}"
    result_path = worker_dir / "result.json" if result_path is None else _inside(root, result_path)
    manifest_path = root / "manifest.json"
    manifest = _read(manifest_path)
    schema = int(result.get("result_schema", manifest.get("result_schema", 1)))
    checks, unhashed, measured_hashes = [], [], {}

    def check(path, expected, role):
        path = _inside(root, path)
        if not path.is_file():
            raise FileNotFoundError(f"Missing scientific {role}: {path}")
        if path not in measured_hashes:
            shared = role in {"frozen_source", "exact_basal", "frozen_case_basal", "frozen_basal_input",
                              "worker_manifest", "frozen_manifest"}
            stat = path.stat()
            identity = (str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            if shared and shared_hash_cache is not None and identity in shared_hash_cache:
                measured_hashes[path] = shared_hash_cache[identity]
            else:
                measured_hashes[path] = _sha(path)
                if shared and shared_hash_cache is not None:
                    shared_hash_cache[identity] = measured_hashes[path]
        measured = measured_hashes[path]
        if expected is not None and verify_hashes and measured != expected:
            raise ValueError(f"Scientific {role} hash mismatch: {path}")
        record = dict(path=str(path.relative_to(root)), role=role, actual_sha256=measured,
                      expected_sha256=expected, identity_verified=bool(expected and verify_hashes))
        checks.append(record)
        if expected is None:
            unhashed.append(record["path"])
        return path

    manifest_digest = _sha(manifest_path)
    if result_path.is_file():
        _same_numeric_tree(result, _read(result_path), "worker_result_on_disk")
        check(result_path, None, "worker_result_json")
    check(manifest_path, result.get("manifest_sha256"), "worker_manifest")
    if (root / "freeze.json").is_file():
        check(manifest_path, _read(root / "freeze.json")["manifest_sha256"], "frozen_manifest")
    declared_cases = manifest["cases"]
    matches = [row for row in declared_cases if row["task"] == case["task"]]
    if len(matches) != 1:
        raise ValueError("Worker task missing or duplicated in frozen case manifest")
    for key in ("task", "model", "stage", "requested_reference_length", "native_length", "cells", "cue_width", "cue_area"):
        if key in matches[0] and case.get(key) != matches[0][key]:
            raise ValueError("Worker case differs from frozen manifest: " + key)
        if key in case and result.get("case", {}).get(key) != case[key]:
            raise ValueError("Worker result differs from requested case: " + key)
    readout_source = None
    for name, digest in manifest.get("source_sha256", {}).items():
        source = check(root / "frozen_sources" / name, digest, "frozen_source")
        if source.name == "signature_readouts_20261020.py":
            readout_source = source
    declared_artifacts = result.get("artifact_sha256", {})
    for name, digest in declared_artifacts.items():
        check(name, digest, "worker_scientific_artifact")
    if schema >= 2:
        scientific = {str(p.relative_to(root)) for p in worker_dir.rglob("*")
                      if p.is_file() and (p.suffix == ".npz" or p.suffix == ".json")
                      and p.name not in ("result.json", "status.json", "failed.json")}
        missing = scientific - set(declared_artifacts)
        if missing:
            raise ValueError("Schema-2 scientific artifacts lack hash coverage: " + ", ".join(sorted(missing)))
    basal_record = result["basal"]
    basal_path = check(basal_record["checkpoint"], basal_record.get("checkpoint_sha256") or
                       basal_record.get("actual_sha256"), "exact_basal")
    frozen_basal = matches[0].get("checkpoint")
    if frozen_basal:
        expected_basal = (matches[0].get("checkpoint_sha256") or
                          manifest.get("input_sha256", {}).get(frozen_basal))
        check(basal_path, expected_basal, "frozen_case_basal")
        check(frozen_basal, manifest.get("input_sha256", {}).get(frozen_basal, expected_basal),
              "frozen_basal_input")
    frozen_module = None
    if recompute_readouts:
        if readout_source is None:
            raise ValueError("Numeric readout audit requires the declared frozen readout source")
        spec = importlib.util.spec_from_file_location("signature_frozen_readout_audit_" + _sha(readout_source)[:12], readout_source)
        frozen_module = importlib.util.module_from_spec(spec)
        # Compile verified source directly: importing for analysis must not add
        # bytecode cache files to an immutable frozen executable directory.
        exec(compile(readout_source.read_bytes(), str(readout_source), "exec"), frozen_module.__dict__)
    trials = []
    for label, trial, legacy_dir in _trial_items(result):
        folder = _inside(root, trial["trial_directory"]) if "trial_directory" in trial else worker_dir / legacy_dir
        artifacts = trial.get("artifacts", {})
        readout_record = artifacts.get("readouts", {})
        trace_record = artifacts.get("trace", {})
        readouts_name = _inside(root, readout_record.get("path", folder / trial.get("readouts", "readouts.json")))
        original_readout_digest = (readout_record.get("sha256") or
                                   declared_artifacts.get(str(readouts_name.relative_to(root))))
        readouts_path = check(readouts_name, original_readout_digest, "readout_json")
        trace_name = _inside(root, trace_record.get("path", folder / "trace.npz"))
        trace_path = check(trace_name, trace_record.get("sha256") or trial.get("full_trace_sha256") or
                           declared_artifacts.get(str(trace_name.relative_to(root))), "full_native_trace")
        tail_path = folder / "native_tail.json"
        if tail_path.is_file() and "native_tail" not in artifacts:
            check(tail_path, None, "legacy_native_tail_json")
        for event in trial["events"]:
            if "archive" in event:
                archive = event["archive"]
                check(archive["path"], archive.get("sha256"), "native_phase")
        for name, record in artifacts.items():
            if name not in ("readouts", "trace"):
                check(record["path"], record.get("sha256"), name)
        measurement = _read(readouts_path)
        native_trace = _trace_mapping(trace_path, case["model"])
        basal_trace = _trace_mapping(basal_path, case["model"])
        if recompute_readouts:
            regenerated = frozen_module.analyze_trace(native_trace, _readout_parameters(case, basal_trace), basal_trace)
            _same_numeric_tree(measurement["series"], regenerated["series"], "series")
            _same_numeric_tree(measurement["basal_descriptor"], regenerated["basal_descriptor"], "basal_descriptor")
            _same_numeric_tree(measurement["final"], regenerated["final"], "final")
            for key in ("incumbent_intervals", "challenger_own_intervals", "challenger_fixed_intervals",
                        "challenger_own_intervals_by_frame", "challenger_fixed_intervals_by_frame",
                        "final_d", "final_replacement", "readout_eligible", "ineligibility_reason"):
                _same_numeric_tree(measurement[key], regenerated[key], key)
        series = measurement["series"]
        if not np.array_equal(np.asarray(series["time"]), native_trace["time"]):
            raise ValueError("Readout and native trace times differ")
        for key in ("final_d", "readout_eligible", "ineligibility_reason"):
            if key in trial:
                _same_numeric_tree(trial[key], measurement[key], "trial." + key)
        events = trial["events"]
        resolved_events = [dict(event, archive=dict(event["archive"],
                           path=str(_inside(root, event["archive"]["path"]))))
                           if "archive" in event else event for event in events]
        state_fields, state_declaration = (declared_state_fields(case["model"], root, manifest)
            if all("archive" in event for event in events) else (None, None))
        phase_audit = audit_native_phase_join(trace_path, resolved_events, native_trace, state_fields=state_fields)
        phase_audit["state_declaration"] = state_declaration
        on_intervals = trial.get("phase_on_intervals", [[e["start"], e["stop"]] for e in events
                            if e["name"] == "cue_on" or str(e["name"]).endswith("on")])
        recovery = trial.get("recovery_interval", next(([e["start"], e["stop"]] for e in events
                                                        if e["name"] == "recovery"), None))
        on_events = [dict(interval=[float(a), float(b)],
                          D0=threshold_event(series, a, b, 0.),
                          D0_5=threshold_event(series, a, b, .5)) for a, b in on_intervals]
        recovery_events = (dict(D0=threshold_event(series, *recovery, 0.),
                               D0_5=threshold_event(series, *recovery, .5)) if recovery else None)
        own = measurement["challenger_own_intervals_by_frame"]
        peak = _regional_peak(native_trace, own)
        activation_threshold = .5 * float(measurement["basal_descriptor"]["excess_peak"])
        activation_series = dict(time=series["time"], own_excess_peak=peak,
                                 peak_eligible=[_finite(v) for v in peak])
        activation = [threshold_event(activation_series, a, b, activation_threshold,
                                     value_key="own_excess_peak", eligibility_key="peak_eligible")
                      for a, b in on_intervals]
        continuous = len(on_intervals) == 1 and not trial.get("schedule", {}).get("pulse_count")
        tail = trial.get("field_tail", _read(folder / "native_tail.json") if (folder / "native_tail.json").is_file() else {})
        if tail_path.is_file():
            _same_numeric_tree(tail, _read(tail_path), "native_tail_sidecar")
        if "field_tail_pass" in tail:
            _same_numeric_tree(trial.get("field_tail_pass"), tail["field_tail_pass"], "field_tail_flag")
        trials.append(dict(label=label, trial_directory=str(folder.relative_to(root)),
            assay_kind="continuous_cue" if continuous else "train" if on_intervals else "zero_control",
            continuous_cue_on_events=on_events[0] if continuous else None,
            phase_on_events=on_events, recovery_events=recovery_events,
            activation_events=activation, activation_threshold=activation_threshold,
            activation_definition="Own-cue max(current-min positive excess) >= half frozen basal excess peak; cue-ON only.",
            native_readout_numeric_audit=bool(recompute_readouts),
            full_trace_sha256=measured_hashes[trace_path], readouts_current_sha256=measured_hashes[readouts_path],
            native_phase_join_audit=phase_audit,
            basal_descriptor=measurement["basal_descriptor"], schedule=trial.get("schedule", {}),
            events=events, phase_on_intervals=on_intervals, recovery_interval=recovery,
            original_readout_hash_available=bool(original_readout_digest),
            field_tail=tail, field_tail_pass=trial.get("field_tail_pass"),
            final_d=measurement["final_d"], final_replacement=measurement["final_replacement"],
            readout_eligible=measurement["readout_eligible"],
            ineligibility_reason=measurement["ineligibility_reason"],
            endpoints=pulse_gap_endpoints(series, events),
            regional_amount_series={key: value for key, value in series.items() if key == "time" or
                key.startswith("E_") or key.startswith("roi_") or key in ("D_half", "D_half_raw", "share_half", "contrast", "current_min_floor", "readout_eligible")},
            own_excess_peak_series=peak))
    return dict(task=case["task"], model=case["model"], case=case,
                basal_sha256=measured_hashes[basal_path], manifest_sha256=manifest_digest,
                analysis_source_sha256=_sha(Path(__file__)),
                result_schema=schema, trials=trials,
                provenance=dict(declared_hash_verification_enabled=bool(verify_hashes),
                    numeric_readout_recomputation=bool(recompute_readouts), checks=checks,
                    files_without_original_declared_hash=sorted(set(unhashed)),
                    scientific_scope="Verified declared identities and independent numeric recomputation; missing historical digests are not reconstructed."))


def width_ratio_table(worker_analyses, *, narrow_width=.95, broad_width=7.6, event_key="D0_5"):
    """Pair completed width assays by model, length and area; retain missing pairs."""
    groups = {}
    if event_key not in ("D0", "D0_5"):
        raise ValueError("Width ratios support D0 or D0_5 sampled event times")
    for analysis in worker_analyses:
        case = analysis["case"]
        if case.get("stage") != "width":
            continue
        width = float(case["requested_width_label"])
        if not (math.isclose(width, narrow_width) or math.isclose(width, broad_width)):
            continue
        key = (case["model"], case["requested_reference_length"], case["integrated_area_factor"])
        slot = "narrow" if math.isclose(width, narrow_width) else "broad"
        if slot in groups.setdefault(key, {}):
            raise ValueError("Duplicate width condition in paired scientific table")
        trial = next((t for t in analysis["trials"] if t["label"] == "trial"), None)
        event = (trial.get("continuous_cue_on_events") or {}).get(event_key) if trial else None
        groups[key][slot] = dict(event=event, basal_sha256=analysis["basal_sha256"], task=analysis["task"])
    rows = []
    for (model, length, area), pair in sorted(groups.items()):
        narrow, broad = pair.get("narrow"), pair.get("broad")
        if not narrow or not broad or narrow["event"] is None or broad["event"] is None:
            ratio = dict(status="unavailable", lower=None, upper=None, reason="Missing completed paired continuous-cue output.")
        elif narrow["basal_sha256"] != broad["basal_sha256"]:
            raise ValueError("Narrow and broad width arms use different exact basals")
        else:
            ratio = replacement_time_ratio(narrow["event"], broad["event"])
        rows.append(dict(model=model, requested_reference_length=length, integrated_area_factor=area,
                         narrow_width=narrow_width, broad_width=broad_width, event_key=event_key,
                         pair=pair, ratio=ratio))
    return rows


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def analyze_completed_workers(run_root, output, *, task_indices=None):
    """Stream a read-only worker audit into a separate results metadata folder.

    The snapshot reports missing workers explicitly and makes no whole-brief
    completion claim. Existing worker scientific directories are never changed.
    """
    root, output = Path(run_root).resolve(), Path(output).resolve()
    if not output.is_relative_to(root / "results"):
        raise ValueError("Analysis output must be beneath the campaign results directory")
    manifest = _read(root / "manifest.json")
    compact, event_rows, missing, errors, worker_digests, phase_checks = [], [], [], [], {}, []
    all_cases = manifest["cases"]
    if task_indices is None:
        cases = all_cases
    else:
        selected = set(map(int, task_indices))
        known = {int(case["task"]) for case in all_cases}
        if selected - known:
            raise ValueError("Requested task indices absent from frozen manifest")
        cases = [case for case in all_cases if int(case["task"]) in selected]
    shared_hash_cache = {}
    source_hash = _sha(Path(__file__))
    summary = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                   producer_state="running", analysis_source_sha256=source_hash,
                   manifest_sha256=_sha(root / "manifest.json"), stage=manifest.get("stage"),
                   planned_workers=len(cases), manifest_workers=len(all_cases), analyzed_workers=0,
                   selected_tasks=[int(case["task"]) for case in cases],
                   missing_tasks=[], errors=[], declared_subset_complete=False,
                   full_brief_complete=False,
                   scope="Independent finite-readout/source audit of this frozen campaign only; no figure rendering or manuscript integration.")
    _write_json(output / "status.json", summary)
    for case in cases:
        task = int(case["task"])
        worker = root / "workers" / f"{task:04d}"
        result_path = worker / "result.json"
        if not result_path.is_file():
            missing.append(task)
            continue
        status_path = worker / "status.json"
        if status_path.is_file() and _read(status_path).get("state") != "complete":
            missing.append(task)
            continue
        try:
            analysis = analyze_worker_output(case, result_path, run_root=root,
                                             shared_hash_cache=shared_hash_cache)
            phase_checks.extend(trial.get("native_phase_join_audit", {}).get("passed") is True
                                for trial in analysis["trials"])
            worker_analysis_path = output / "workers" / f"{task:04d}.json"
            _write_json(worker_analysis_path, analysis)
            worker_digests[str(task)] = dict(path=str(worker_analysis_path.relative_to(root)),
                                            sha256=_sha(worker_analysis_path),
                                            result_sha256=_sha(result_path))
            compact.append({key: analysis[key] for key in ("task", "model", "case", "basal_sha256")} |
                           {"trials": [{key: trial[key] for key in ("label", "continuous_cue_on_events")}
                                       for trial in analysis["trials"]]})
            for trial in analysis["trials"]:
                event_rows.append(dict(task=task, model=case["model"],
                    requested_reference_length=case.get("requested_reference_length"),
                    requested_width_label=case.get("requested_width_label"),
                    integrated_area_factor=case.get("integrated_area_factor"),
                    pulse_count=case.get("pulse_count"),
                    gap_multiple_of_pulse=case.get("gap_multiple_of_pulse"),
                    total_on_factor_of_calibration=case.get("total_on_factor_of_calibration"),
                    label=trial["label"], assay_kind=trial["assay_kind"],
                    continuous_cue_on_events=trial["continuous_cue_on_events"],
                    recovery_events=trial["recovery_events"], activation_events=trial["activation_events"],
                    final_d=trial["final_d"], final_replacement=trial["final_replacement"],
                    readout_eligible=trial["readout_eligible"], field_tail_pass=trial["field_tail_pass"],
                    own_roi_overlap_any=any(trial["regional_amount_series"].get("roi_own_overlap", [])),
                    fixed_roi_overlap_any=any(trial["regional_amount_series"].get("roi_fixed_overlap", []))))
            summary["analyzed_workers"] += 1
        except Exception as error:
            errors.append(dict(task=task, error_type=type(error).__name__, error=str(error)))
        summary.update(generated_utc=datetime.now(timezone.utc).isoformat(),
                       missing_tasks=missing, errors=errors)
        _write_json(output / "status.json", summary)
    _write_json(output / "event_rows.json", event_rows)
    for event_key in ("D0", "D0_5"):
        _write_json(output / f"width_ratios_{event_key}.json", width_ratio_table(compact, event_key=event_key))
    required_width_grid = {(model, length, width, area)
        for model in ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
        for length in (19, 23, 27, 31, 35, 39) for width in (.95, 1.9, 3.8, 7.6, 15.2)
        for area in (.5, 1., 2.)}
    observed_width_grid = {(a["case"]["model"], a["case"].get("requested_reference_length"),
                            a["case"].get("requested_width_label"), a["case"].get("integrated_area_factor"))
                           for a in compact if a["case"].get("stage") == "width"}
    summary.update(generated_utc=datetime.now(timezone.utc).isoformat(), producer_state="stopped",
                   missing_tasks=missing, errors=errors,
                   declared_subset_complete=not missing and not errors,
                   full_width_factor_coverage=observed_width_grid == required_width_grid,
                   width_cases_audited=len(observed_width_grid), full_width_required_cases=450,
                   worker_analysis_sha256=worker_digests,
                   all_native_phase_join_audits_pass=bool(phase_checks) and all(phase_checks),
                   model_coverage=sorted({a["model"] for a in compact}),
                   outputs={str(path.relative_to(output)): _sha(path) for path in
                            (output / "event_rows.json", output / "width_ratios_D0.json",
                             output / "width_ratios_D0_5.json")})
    _write_json(output / "status.json", summary)
    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--task-start", type=int, default=0)
    parser.add_argument("--task-stop", type=int)
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--shard-count", type=int)
    arguments = parser.parse_args()
    cases = _read(arguments.run_root / "manifest.json")["cases"]
    selected = [int(case["task"]) for case in cases if int(case["task"]) >= arguments.task_start
                and (arguments.task_stop is None or int(case["task"]) < arguments.task_stop)]
    if arguments.shard_count is not None or arguments.shard_index is not None:
        if (arguments.shard_count is None or arguments.shard_index is None or
                not 0 <= arguments.shard_index < arguments.shard_count):
            parser.error("Require 0 <= shard-index < shard-count")
        selected = [task for task in selected if task % arguments.shard_count == arguments.shard_index]
    outcome = analyze_completed_workers(arguments.run_root, arguments.output, task_indices=selected)
    print(json.dumps(outcome, indent=2, allow_nan=False))
    raise SystemExit(1 if outcome["errors"] else 0)
