"""Coverage-gated Figure 2 assembly from independently audited native assays.

This program launches no integrations. The default export is explicitly
provisional. Final export additionally requires complete factor coverage and a
hash-linked independent completion review of stages 0–5 and panels A–F.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
from pathlib import Path
import textwrap

import numpy as np

from signature_analysis_20261020 import replacement_time_ratio, threshold_event

MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
NAMES = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "deBelly", "Holmes 3", "Spring")))
COLORS = dict(zip(MODELS, ("#3975A4", "#B85942", "#82599C", "#318875", "#B68A22")))
LENGTHS = (19, 23, 27, 31, 35, 39)
WIDTHS = (.95, 1.9, 3.8, 7.6, 15.2)
AREAS = (.5, 1., 2.)
GAPS = (0., .1, .5, 1., 4., 16.)
PULSES = (2, 4, 8)
ON_FACTORS = (1., 1.5, 2.)
EXAMPLE = dict(pulse_count=4, gap_multiple_of_pulse=4., total_on_factor_of_calibration=1.5)


def example_length(model, calibration_unavailable):
    """Choose by calibration availability alone, never by native response fate."""
    eligible = [length for length in LENGTHS
                if ("train", model, length, 4, 4., 1.5) not in calibration_unavailable]
    return max(eligible) if eligible else None


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def condition_key(case):
    model, length = case["model"], int(case["requested_reference_length"])
    if case["stage"] == "width":
        return ("width", model, length, float(case["requested_width_label"]),
                float(case["integrated_area_factor"]))
    if case["stage"] == "train":
        return ("train", model, length, int(case["pulse_count"]),
                float(case["gap_multiple_of_pulse"]), float(case["total_on_factor_of_calibration"]))
    raise ValueError("Figure 2 core supports only width and train scientific records")


def required_conditions():
    widths = {("width", m, l, w, a) for m in MODELS for l in LENGTHS for w in WIDTHS for a in AREAS}
    trains = {("train", m, l, n, g, o) for m in MODELS for l in LENGTHS
              for n in PULSES for g in GAPS for o in ON_FACTORS}
    return widths, trains


def load_audited_records(analysis_directories):
    """Read compact audit tables, verify their receipts and preserve omissions.

    Missing/unscorable work never creates a synthetic D or a censoring time.
    Duplicate shard inputs are accepted only when the scientific row agrees.
    """
    records, receipts = {}, []
    for directory in map(lambda p: Path(p).resolve(), analysis_directories):
        root = next((p for p in directory.parents if (p / "manifest.json").is_file()
                     and (p / "workers").is_dir()), None)
        if root is None:
            raise ValueError("Audit directory must belong to an immutable native run root")
        status_path = directory / "status.json"
        status = read(status_path)
        if status.get("producer_state") != "stopped":
            raise ValueError("Audit producer is not terminal: " + str(directory))
        if status.get("manifest_sha256") != sha(root / "manifest.json"):
            raise ValueError("Audit manifest identity no longer matches run root")
        table = directory / "event_rows.json"
        expected = status.get("outputs", {}).get("event_rows.json")
        if expected is None or sha(table) != expected:
            raise ValueError("Independent scientific table lacks its verified digest")
        manifest = read(root / "manifest.json")
        cases = {int(case["task"]): case for case in manifest["cases"]}
        selected = set(status.get("selected_tasks", cases))
        accepted = 0
        for row in read(table):
            if row["label"] != "trial":
                continue
            task = int(row["task"])
            if task not in selected or task not in cases:
                raise ValueError("Scientific table contains task outside declared audit selection")
            case = cases[task]
            if case["stage"] not in ("width", "train"):
                continue
            for key in ("model", "requested_reference_length", "requested_width_label", "integrated_area_factor"):
                if key in row and row[key] is not None and case.get(key) != row[key]:
                    raise ValueError("Scientific table condition differs from frozen manifest")
            key = condition_key(case)
            payload = dict(case=case, measurement=row,
                           source=dict(run_root=str(root), task=task, analysis_directory=str(directory),
                                       audit_status_sha256=sha(status_path), event_table_sha256=expected,
                                       manifest_sha256=status["manifest_sha256"],
                                       worker_analysis_receipt=status.get("worker_analysis_sha256", {}).get(str(task))))
            if key in records:
                if records[key]["measurement"] != row:
                    raise ValueError("Conflicting independent measurements for the same condition")
                continue
            records[key] = payload
            accepted += 1
        receipts.append(dict(path=str(directory), run_root=str(root), status_sha256=sha(status_path),
                             scientific_table_sha256=expected, accepted_conditions=accepted,
                             missing_tasks=status.get("missing_tasks", []), errors=status.get("errors", []),
                             declared_subset_complete=status.get("declared_subset_complete", False),
                             native_phase_join_gate=status.get("all_native_phase_join_audits_pass")))
    return records, receipts


def audited_worker_metadata(record):
    """Load a hash-linked scientific audit without rereading large native arrays."""
    source = record["source"]
    root = Path(source["run_root"])
    receipt = source.get("worker_analysis_receipt")
    path = (root / receipt["path"] if receipt else
            Path(source["analysis_directory"]) / "workers" / f"{source['task']:04d}.json")
    if receipt and sha(path) != receipt["sha256"]:
        raise ValueError("Selected detailed scientific audit changed")
    worker = read(path)
    if worker["manifest_sha256"] != source["manifest_sha256"]:
        raise ValueError("Selected audit manifest mismatch")
    trial = next(t for t in worker["trials"] if t["label"] == "trial")
    folder = root / trial["trial_directory"]
    if not trial["native_readout_numeric_audit"]:
        raise ValueError("Selected trace has no independent numeric reconstruction")
    return trial, folder, dict(path=str(path), current_sha256=sha(path),
                               originally_hashed=bool(receipt))


def scientific_master(config_path, output):
    """Freeze all main-train receipts and planned calibration-unavailable cells."""
    config_path, output = Path(config_path).resolve(), Path(output).resolve()
    config=read(config_path)
    records,receipts=load_audited_records(config["analysis_directories"])
    unavailable,unavailable_receipt=load_calibration_unavailable(config["unavailable_trains"])
    gate=coverage(records,receipts,calibration_unavailable=unavailable)
    if (gate["audited_width_conditions"]!=450 or gate["audited_train_conditions"]!=1458 or
            gate["calibration_unavailable_train_conditions"]!=162 or gate["audit_errors"] or
            gate["missing_width_conditions"] or gate["missing_train_conditions"]):
        raise ValueError("Scientific master requires all completed native assay conditions")
    selected=[]
    for key in sorted(k for k in records if k[0]=="train"):
        record=records[key]
        trial,_,reference=audited_worker_metadata(record)
        native_join=trial.get("native_phase_join_audit",{})
        if not native_join.get("passed"):
            raise ValueError("Scientific master selected a failed native phase join")
        worker=read(reference["path"])
        selected.append(dict(condition=list(key),case=record["case"],source=record["source"],
            independent_worker_audit=reference,analysis_source_sha256=worker["analysis_source_sha256"],
            all_native_state_joins_and_phase_samples_pass=True,
            native_readout_numeric_audit_pass=bool(trial["native_readout_numeric_audit"]),
            native_trace_sha256=trial["full_trace_sha256"],native_readout_sha256=trial["readouts_current_sha256"],
            final_d=trial["final_d"],final_replacement=trial["final_replacement"],
            readout_eligible=trial["readout_eligible"],field_tail_pass=trial["field_tail_pass"],
            exact_phase_endpoints=len(trial["endpoints"])))
    report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,
        all_native_readout_recomputations_pass=True,all_native_phase_join_audits_pass=True,
        scope="All 1458 executed main trains; width readout coverage is contextual and has separate audit receipts.",
        width_readout_conditions_in_config=450,native_train_conditions=1458,calibration_unavailable_train_conditions=162,
        planned_train_conditions=1620,all_planned_conditions_scientifically_accounted=True,
        scalar_fate_and_field_tail_separate=True,source_config=dict(path=str(config_path),sha256=sha(config_path)),
        selected_worker_audits=selected,source_audit_groups=receipts,
        calibration_unavailable_evidence=unavailable_receipt,
        calibration_unavailable_conditions=list(unavailable.values()),
        checker_history="Valid original chemical audits are retained; deBelly and Spring use the corrected declared-state checker. Superseded diagnostic-continuity errors remain in original immutable audit receipts and are not native worker failures.",
        renderer_source_sha256=sha(Path(__file__)))
    write(output,report)
    return dict(path=str(output),sha256=sha(output),passed=True,widths=450,trains=1458,unavailable=162)


def audit_worker(record):
    """Load one selected detailed audit; verify its native/readout identities."""
    trial, folder, _ = audited_worker_metadata(record)
    for name, expected in (("readouts.json", trial["readouts_current_sha256"]),
                           ("trace.npz", trial["full_trace_sha256"])):
        if sha(folder / name) != expected:
            raise ValueError("Selected independently audited native evidence changed: " + name)
    return trial, folder


def load_relative_activation(directories):
    """Load same-peak relative-width diagnostics from terminal S4 native audit."""
    if not directories:
        return {}, None
    if isinstance(directories, (str, Path)):
        directories = [directories]
    result, receipts = {}, []
    for directory in map(lambda p: Path(p).resolve(), directories):
        root = next((p for p in directory.parents if (p / "manifest.json").is_file()
                     and (p / "workers").is_dir()), None)
        if root is None:
            raise ValueError("Supporting audit must belong to its frozen native run root")
        status = read(directory / "status.json")
        if status.get("producer_state") != "stopped":
            raise ValueError("Relative-width overlay requires terminal supporting audit snapshots")
        if status["manifest_sha256"] != sha(root / "manifest.json"):
            raise ValueError("Supporting manifest identity mismatch")
        table = directory / "event_rows.json"
        if sha(table) != status["outputs"]["event_rows.json"]:
            raise ValueError("Supporting scientific table identity mismatch")
        cases = {int(c["task"]): c for c in read(root / "manifest.json")["cases"]}
        for row in read(table):
            case = cases[int(row["task"])]
            if case.get("assay") != "relative_width_initiation" or row["label"] != "trial":
                continue
            key = (case["model"], int(case["requested_reference_length"]), float(case["strength_factor"]))
            if key in result and result[key]["measurement"] != row:
                raise ValueError("Conflicting relative-width scientific activation condition")
            result[key] = dict(case=case, measurement=row)
        receipts.append(dict(path=str(directory), status_sha256=sha(directory / "status.json"),
                             event_rows_sha256=sha(table)))
    expected = {(m, l, a) for m in MODELS[:3] for l in LENGTHS for a in AREAS}
    if set(result) != expected:
        raise ValueError("Supporting relative-width activation grid is incomplete")
    return result, dict(audits=receipts, conditions=len(result),
                        interpretation="Same peak as each fixed-width strength; width=20%L changes commanded area. Focal three-model sensitivity only.")


def load_calibration_unavailable(path):
    """Account for planned trains whose native calibration has no sampled hit."""
    if path is None:
        return {}, None
    path = Path(path).resolve()
    table = read(path)
    rows = table["conditions"]
    result, verified = {}, {}
    _, required = required_conditions()
    for row in rows:
        case = dict(row, stage="train")
        key = condition_key(case)
        if key not in required or key in result:
            raise ValueError("Invalid or duplicate calibration-unavailable train condition")
        if (row.get("status") != "calibration_unavailable" or row.get("train_executed") is not False
                or row.get("calibration_duration") is not None
                or row.get("calibration_artifact_hashes_verified") is not True):
            raise ValueError("Unavailable train requires an audited absent calibration, never a substituted duration")
        source = Path(row["calibration_result_path"])
        if not source.is_absolute():
            source = path.parent / source
        if str(source) not in verified:
            if sha(source) != row["calibration_result_sha256"]:
                raise ValueError("Unavailable calibration scientific identity mismatch")
            verified[str(source)] = read(source)
        native = verified[str(source)]
        calibration = native["calibration"]
        if (native["case"]["model"] != row["model"] or
                native["case"]["requested_reference_length"] != row["requested_reference_length"] or
                native["basal"]["checkpoint_sha256"] != row["exact_basal_sha256"] or
                native["manifest_sha256"] != row["calibration_manifest_sha256"] or
                calibration["status"] != "no_sampled_commitment_through_finite_horizon" or
                calibration["train_ready"] is not False or
                calibration["finite_duration_horizon"] != row["searched_on_cap"] or
                not calibration["evaluations"] or
                max(e["duration"] for e in calibration["evaluations"]) != row["searched_on_cap"] or
                any(e["replacement"] for e in calibration["evaluations"])):
            raise ValueError("Unavailable train record contradicts its native calibration")
        result[key] = row
    return result, dict(path=str(path), sha256=sha(path), conditions=len(result),
                        no_finite_train_duration_assigned=True)


def coverage(records, receipts, completion_review=None, *, calibration_unavailable=None):
    widths, trains = required_conditions()
    present = set(records)
    unavailable = set(calibration_unavailable or {})
    if unavailable & present or unavailable - trains:
        raise ValueError("Calibration-unavailable conditions must be planned trains without executed output")
    resolved = {(r["source"]["run_root"], int(r["source"]["task"])) for r in records.values() if "source" in r}
    unresolved_audits = []
    for receipt in receipts:
        missing = [task for task in receipt["missing_tasks"]
                   if (receipt.get("run_root"), int(task)) not in resolved]
        errors = [error for error in receipt["errors"]
                  if (receipt.get("run_root"), int(error["task"])) not in resolved]
        if missing or errors:
            unresolved_audits.append(dict(path=receipt["path"], missing_tasks=missing, errors=errors))
    report = dict(required_width_conditions=len(widths), audited_width_conditions=len(widths & present),
                  required_train_conditions=len(trains), audited_train_conditions=len(trains & present),
                  calibration_unavailable_train_conditions=len(unavailable),
                  scientifically_accounted_train_conditions=len((trains & present) | unavailable),
                  missing_width_conditions=[list(k) for k in sorted(widths - present)],
                  missing_train_conditions=[list(k) for k in sorted(trains - present - unavailable)],
                  unexpected_conditions=[list(k) for k in sorted(present - widths - trains)],
                  audit_errors=unresolved_audits,
                  final_export_eligible=False, finite_field_tail_filters_fate=False)
    if completion_review:
        review_path = Path(completion_review).resolve()
        review = read(review_path)
        required_panels = set("ABCDEF")
        required_stages = set(map(str, range(6)))
        if (review.get("independent_review_passed") is not True or
                set(review.get("panel_completion", {})) != required_panels or
                not all(review["panel_completion"].values()) or
                set(review.get("stage_completion", {})) != required_stages or
                not all(review["stage_completion"].values())):
            raise ValueError("Final export requires the independent all-panel/stage scientific review")
        evidence = review.get("evidence", [])
        if not evidence or set(e["stage"] for e in evidence) != required_stages:
            raise ValueError("Final completion review must link actual evidence for all six stages")
        for item in evidence:
            path = Path(item["path"])
            if not path.is_absolute():
                path = review_path.parent / path
            if sha(path) != item["sha256"]:
                raise ValueError("Independent completion evidence hash mismatch")
        report["independent_completion_review"] = dict(path=str(review_path), sha256=sha(review_path))
        report["final_export_eligible"] = bool(not (widths - present) and not (trains - present - unavailable)
                                               and not report["unexpected_conditions"]
                                               and not report["audit_errors"])
    return report


def ratio_records(records, *, area=1.):
    rows = []
    for model in MODELS:
        for length in LENGTHS:
            pair = [records.get(("width", model, length, width, area)) for width in (.95, 7.6)]
            if all(pair):
                identities = [p["case"].get("checkpoint_sha256") for p in pair]
                if not all(identities) or identities[0] != identities[1]:
                    raise ValueError("Width ratio requires the same exact frozen basal identity")
                events = [(p["measurement"].get("continuous_cue_on_events") or {}).get("D0_5") for p in pair]
                ratio = replacement_time_ratio(*events) if all(events) else dict(status="unavailable")
            else:
                ratio = dict(status="unavailable", reason="Missing independently audited width arm")
            rows.append(dict(model=model, length=length, area=area, ratio=ratio))
    return rows


def ratio_text(ratio):
    if ratio.get("status") == "unavailable":
        return "unavailable"
    lo, hi = ratio.get("lower"), ratio.get("upper")
    if hi is None:
        return f">{lo:.3g}" if lo and lo > 0 else "uninformative"
    if lo == 0:
        return f"<{hi:.3g}"
    return f"{lo:.3g}–{hi:.3g}"


def regional_peak_vectorized(active, x, intervals, *, sample_kind):
    """Exact native-cell maximum or piecewise-linear nodal footprint maximum."""
    active, x = np.asarray(active), np.asarray(x)
    bounds = np.asarray(intervals, dtype=float)
    if bounds.shape != (len(active), 1, 2):
        raise ValueError("Single contiguous own-cue support required for reporter sensitivity")
    lo, hi = bounds[:, 0, 0], bounds[:, 0, 1]
    if sample_kind == "cells":
        dx = x[1] - x[0]
        edges = np.r_[x[0] - .5 * dx, x + .5 * dx]
        mask = (edges[:-1][None, :] < hi[:, None]) & (edges[1:][None, :] > lo[:, None])
        maximum = np.max(np.where(mask, active, -np.inf), axis=1)
    else:
        mask = (x[None, :] > lo[:, None]) & (x[None, :] < hi[:, None])
        maximum = np.max(np.where(mask, active, -np.inf), axis=1)
        rows = np.arange(len(active))
        for edge in (lo, hi):
            index = np.clip(np.searchsorted(x, edge, side="right") - 1, 0, len(x) - 2)
            fraction = (edge - x[index]) / (x[index + 1] - x[index])
            interpolated = active[rows, index] * (1 - fraction) + active[rows, index + 1] * fraction
            maximum = np.maximum(maximum, interpolated)
    if not np.isfinite(maximum).all():
        raise ValueError("Native grid has no finite own-support reporter maximum")
    return maximum - np.min(active, axis=1)


def reporter_threshold_sensitivity(records, output):
    """Audit all fixed reporter thresholds before a main initiation definition.

    The four thresholds are 1/5/10/20% of each model's exact L19 basal excess
    peak, held constant across lengths, widths and strengths. The original
    per-basal half-peak event remains a separately labelled height diagnostic.
    """
    widths, _ = required_conditions()
    if not widths <= set(records):
        raise ValueError("Reporter sensitivity requires all 450 independently audited width histories")
    fractions = (.01, .05, .1, .2)
    refs, rows, cache = {}, [], {}
    for model in MODELS:
        record = records[("width", model, 19, 3.8, 1.)]
        trial, folder = audit_worker(record)
        measurement = read(folder / "readouts.json")
        refs[model] = dict(excess_peak=float(measurement["basal_descriptor"]["excess_peak"]),
                           active_peak=float(measurement["basal_descriptor"]["active_peak"]),
                           basal_sha256=record["case"]["checkpoint_sha256"],
                           source=record["source"], readout_sha256=sha(folder / "readouts.json"),
                           units="Model-native active concentration above current minimum; not optical reporter intensity.")
        cache[(model, 19, 3.8, 1.)] = (trial, folder)
    for key in sorted(widths):
        record = records[key]
        _, model, length, width, area = key
        trial, folder = cache.get((model, length, width, area), (None, None))
        if trial is None:
            trial, folder = audit_worker(record)
        measurement = read(folder / "readouts.json")
        with np.load(folder / "trace.npz", allow_pickle=False) as native:
            time = native["time"].copy()
            peak = regional_peak_vectorized(native["active"], native["x"],
                measurement["challenger_own_intervals_by_frame"],
                sample_kind="nodes" if model == "spring" else "cells")
        if not np.allclose(peak, trial["own_excess_peak_series"], rtol=1e-10, atol=1e-12):
            raise ValueError("Independent vectorized native reporter maximum differs from original audit")
        continuous = record["measurement"]["continuous_cue_on_events"]
        interval = continuous["interval"]
        series = dict(time=time.tolist(), peak=peak.tolist(), peak_eligible=[True] * len(time))
        events = {str(fraction): threshold_event(series, *interval, refs[model]["excess_peak"] * fraction,
                                                value_key="peak", eligibility_key="peak_eligible")
                  for fraction in fractions}
        rows.append(dict(model=model, length=length, width=width, area=area,
                         fixed_reference_peak=refs[model]["excess_peak"],
                         fixed_thresholds={str(f): f * refs[model]["excess_peak"] for f in fractions},
                         fixed_reporter_events=events,
                         original_half_per_basal_peak_height_event=record["measurement"]["activation_events"][0],
                         source=record["source"], native_trace_sha256=trial["full_trace_sha256"],
                         native_peak_independent_recomputation_pass=True))
    summaries = []
    for model in MODELS:
        for area in AREAS:
            subset = sorted((r for r in rows if r["model"] == model and r["area"] == area and r["width"] == 3.8),
                            key=lambda r: r["length"])
            for fraction in fractions:
                events = [row["fixed_reporter_events"][str(fraction)] for row in subset]
                summaries.append(dict(model=model, area=area, fixed_width=3.8, reference_fraction=fraction,
                                      reference_threshold=refs[model]["excess_peak"] * fraction,
                                      lengths=[r["length"] for r in subset], events=events,
                                      L39_div_L19=replacement_time_ratio(events[-1], events[0]),
                                      all_six_events_observed=all(e["status"] == "observed" for e in events),
                                      strictly_increasing_nonoverlapping_brackets=bool(
                                          all(e["status"] == "observed" for e in events) and
                                          all(a["upper_elapsed"] < b["lower_elapsed"] for a, b in zip(events, events[1:])))))
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    write(output / "reporter_thresholds.json", refs)
    write(output / "all450_fixed_reporter_events.json", rows)
    write(output / "fixed_width_length_strength_sensitivity.json", summaries)
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), producer_state="stopped", passed=True,
                  native_histories=450, independent_native_peak_recomputations=450,
                  reference_fractions=list(fractions), reference_length=19,
                  threshold_policy="One fixed threshold per model, independent of length/width/strength; complete predeclared1/5/10/20% sensitivity before main selection.",
                  original_readout="Half each run's basal peak is retained as height maturation, not initiation.",
                  main_threshold_selected=False, renderer_source_sha256=sha(Path(__file__)),
                  output_sha256={name: sha(output / name) for name in ("reporter_thresholds.json", "all450_fixed_reporter_events.json", "fixed_width_length_strength_sensitivity.json")})
    write(output / "reporter_sensitivity_audit.json", report)
    return report


def load_reporter_definition(directory, fraction):
    if directory is None and fraction is None:
        return {}, None
    if directory is None or fraction not in (.01, .05, .1, .2):
        raise ValueError("Reporter definition requires the complete sensitivity audit and one declared fraction")
    directory = Path(directory).resolve()
    report = read(directory / "reporter_sensitivity_audit.json")
    if (report.get("passed") is not True or report.get("native_histories") != 450 or
            report.get("independent_native_peak_recomputations") != 450 or
            report.get("reference_fractions") != [.01, .05, .1, .2]):
        raise ValueError("Complete independently reconstructed reporter sensitivity is required")
    for name, expected in report["output_sha256"].items():
        if sha(directory / name) != expected:
            raise ValueError("Reporter sensitivity artifact identity mismatch")
    rows = read(directory / "all450_fixed_reporter_events.json")
    events = {("width", r["model"], int(r["length"]), float(r["width"]), float(r["area"])):
              r["fixed_reporter_events"][str(fraction)] for r in rows}
    required, _ = required_conditions()
    if set(events) != required:
        raise ValueError("Reporter scientific grid does not cover all450 width histories")
    reference_peaks = read(directory / "reporter_thresholds.json")
    for model, field in zip(MODELS, ("active", "active", "rac", "rac_active", "active")):
        reference_peaks[model]["native_field"] = field
        reference_peaks[model]["units"] = "Native active polarity-signal units; no molecular-count, optical or cross-model calibration."
        reference_peaks[model]["selected_fixed_threshold"] = fraction * reference_peaks[model]["excess_peak"]
    return events, dict(path=str(directory), audit_sha256=sha(directory / "reporter_sensitivity_audit.json"),
                        fraction=fraction, reference_length=19, thresholds=reference_peaks,
                        fixed_reference_thresholds_sha256=sha(directory / "reporter_thresholds.json"),
                        selection="Declared main reporter fraction after full1/5/10/20% sensitivity; all sensitivity values remain available.")


def relative_reporter_sensitivity(index_path, reference_directory, output):
    """Independently extend the same four reporter thresholds to all54 S4 cues."""
    index_path, reference_directory = Path(index_path).resolve(), Path(reference_directory).resolve()
    index = read(index_path)
    if not index.get("passed") or not index.get("declared_subset_complete") or index.get("audited_workers") != 270:
        raise ValueError("Relative reporter sensitivity requires the complete S4 independent audit")
    root = index_path.parent.parent
    if sha(root / "manifest.json") != index["manifest_sha256"]:
        raise ValueError("S4 manifest mismatch")
    refs = read(reference_directory / "reporter_thresholds.json")
    reference_audit = read(reference_directory / "reporter_sensitivity_audit.json")
    if not reference_audit.get("passed") or sha(reference_directory / "reporter_thresholds.json") != reference_audit["output_sha256"]["reporter_thresholds.json"]:
        raise ValueError("Fixed reporter reference peak identity mismatch")
    rows = []
    for receipt in index["records"]:
        if receipt["assay"] != "relative_width_initiation":
            continue
        path = root / receipt["audited_worker"]
        if sha(path) != receipt["audited_worker_sha256"]:
            raise ValueError("S4 detailed worker audit identity mismatch")
        worker = read(path)
        case = worker["case"]
        record = dict(case=case, source=dict(run_root=str(root), task=case["task"],
            manifest_sha256=index["manifest_sha256"],
            worker_analysis_receipt=dict(path=receipt["audited_worker"], sha256=receipt["audited_worker_sha256"])))
        trial, folder = audit_worker(record)
        measurement = read(folder / "readouts.json")
        with np.load(folder / "trace.npz", allow_pickle=False) as native:
            time = native["time"].copy()
            peak = regional_peak_vectorized(native["active"], native["x"],
                measurement["challenger_own_intervals_by_frame"], sample_kind="cells")
        if not np.allclose(peak, trial["own_excess_peak_series"], rtol=1e-10, atol=1e-12):
            raise ValueError("Relative cue native reporter reconstruction differs from original audit")
        fractions = (.01, .05, .1, .2)
        series = dict(time=time.tolist(), peak=peak.tolist(), peak_eligible=[True] * len(time))
        model = case["model"]
        events = {str(f): threshold_event(series, *trial["continuous_cue_on_events"]["interval"],
            refs[model]["excess_peak"] * f, value_key="peak", eligibility_key="peak_eligible") for f in fractions}
        rows.append(dict(model=model, length=int(case["requested_reference_length"]),
                         area=case["strength_factor"], cue_width=case["cue_width"],
                         fixed_reporter_events=events, case=case,
                         own_peak_native_recomputation_pass=True, source=record["source"],
                         dose_changes_with_width=True))
    expected = {(m, l, a) for m in MODELS[:3] for l in LENGTHS for a in AREAS}
    if {(r["model"], r["length"], r["area"]) for r in rows} != expected:
        raise ValueError("Relative reporter sensitivity does not cover all54 cases")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    write(output / "all54_relative_reporter_events.json", rows)
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True, native_histories=54,
                  reference_fractions=list(fractions), main_threshold_selected=False,
                  index_sha256=sha(index_path), source_index=str(index_path),
                  native_peak_independent_recomputations=54,
                  fixed_reference_thresholds_sha256=sha(reference_directory / "reporter_thresholds.json"),
                  output_sha256={"all54_relative_reporter_events.json": sha(output / "all54_relative_reporter_events.json")},
                  renderer_source_sha256=sha(Path(__file__)))
    write(output / "relative_reporter_sensitivity_audit.json", report)
    return report


def load_relative_reporter_definition(directory, fraction):
    if directory is None:
        return {}, None
    directory = Path(directory).resolve()
    report = read(directory / "relative_reporter_sensitivity_audit.json")
    path = directory / "all54_relative_reporter_events.json"
    if (fraction not in report["reference_fractions"] or not report.get("passed") or
            report.get("native_histories") != 54 or
            sha(path) != report["output_sha256"][path.name]):
        raise ValueError("Incomplete or changed relative reporter sensitivity")
    rows = {(r["model"], r["length"], r["area"]):
            dict(case=r["case"], measurement=dict(activation_events=[r["fixed_reporter_events"][str(fraction)]]))
            for r in read(path)}
    return rows, dict(path=str(directory), audit_sha256=sha(directory / "relative_reporter_sensitivity_audit.json"),
                      reference_fraction=fraction, fixed_reference_thresholds_sha256=report["fixed_reference_thresholds_sha256"])


def plot_reporter_sensitivity(directory, output):
    """Publication-style full sensitivity PDF; no threshold is chosen here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    directory, output = Path(directory).resolve(), Path(output).resolve()
    audit = read(directory / "reporter_sensitivity_audit.json")
    for name, expected in audit["output_sha256"].items():
        if sha(directory / name) != expected:
            raise ValueError("Reporter plot source identity mismatch")
    rows, refs = read(directory / "all450_fixed_reporter_events.json"), read(directory / "reporter_thresholds.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.get_cmap("viridis")(np.linspace(.05, .9, len(WIDTHS)))
    with PdfPages(output) as pdf:
        for fraction in [*audit["reference_fractions"], None]:
            fig, axes = plt.subplots(3, 5, figsize=(13.5, 9.), layout="constrained")
            description = (f"Fixed reporter threshold = {100 * fraction:g}% of L19 excess peak" if fraction is not None else
                           "Original half-per-basal peak: height maturation diagnostic")
            fig.suptitle("Figure 2 initiation-readout sensitivity\n" + description, fontsize=13)
            for i, area in enumerate(AREAS):
                for j, model in enumerate(MODELS):
                    ax = axes[i, j]
                    for width, color in zip(WIDTHS, colors):
                        selected = sorted((r for r in rows if r["model"] == model and r["area"] == area and r["width"] == width), key=lambda r: r["length"])
                        points = []
                        for row in selected:
                            event = (row["fixed_reporter_events"][str(fraction)] if fraction is not None else
                                     row["original_half_per_basal_peak_height_event"])
                            y, lo, hi, kind = _event_value(event)
                            if y is not None:
                                ax.errorbar(row["length"], y, yerr=[[y - lo], [hi - y]], fmt="^" if kind == "lower_bound" else "o",
                                            color=color, markersize=2.6, capsize=1, linewidth=.6)
                                if kind == "observed":
                                    points.append((row["length"], y))
                            elif kind == "already_above_at_interval_start":
                                ax.plot(row["length"], 0, "D", markerfacecolor="none", color=color, markersize=3)
                        if points:
                            ax.plot(*zip(*points), color=color, linewidth=.8, label=f"w={width:g}")
                    threshold = refs[model]["excess_peak"] * fraction if fraction is not None else None
                    ax.set_title(NAMES[model] + (f"; level={threshold:.4g}" if threshold is not None else ""), fontsize=8.6)
                    ax.set_yscale("symlog", linthresh=.5 if model == "debelly" else .001)
                    ax.set_xlabel("Length label L", fontsize=8)
                    ax.set_xticks([19, 29, 39])
                    ax.tick_params(labelsize=7)
                    if j == 0:
                        ax.set_ylabel(f"Area {area:g}×\nNative threshold time", fontsize=8)
                    if i == 0 and j == 4:
                        ax.legend(fontsize=6.5, frameon=False)
            fig.supxlabel(textwrap.fill("Bars are sampled crossing brackets; triangles are finite lower bounds; open diamonds are already above the reporter level at cue onset. Native units differ by model.", 180), fontsize=8)
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)
    report = dict(path=str(output), sha256=sha(output), pages=5, main_threshold_selected=False,
                  passed=True, source_directory=str(directory),
                  source_audit_sha256=sha(directory / "reporter_sensitivity_audit.json"),
                  renderer_source_sha256=sha(Path(__file__)))
    write(output.parents[1] / "metadata/S4_reporter_sensitivity_display_audit.json", report)
    return report


def verified_organization(path):
    path = Path(path).resolve()
    validation = path.parents[1] / "metadata/validation.json"
    checks = read(validation).get("source_hashes", {})
    declared = [digest for name, digest in checks.items() if name.endswith("/organization_all97.csv")]
    if len(declared) != 1 or sha(path) != declared[0]:
        raise ValueError("Retained organization table fails original display/source identity")
    rows = list(csv.DictReader(path.open()))
    if len(rows) != 97 or set(row["model"] for row in rows) != set(MODELS):
        raise ValueError("Retained all-five organization coverage differs from 97 states")
    return rows, dict(path=str(path), sha256=sha(path), validation=str(validation), validation_sha256=sha(validation))


def _event_value(event):
    if event.get("status") == "observed":
        lo, hi = event["lower_elapsed"], event["upper_elapsed"]
        return .5 * (lo + hi), lo, hi, "observed"
    if event.get("status") == "right_censored":
        value = event["right_censored_at"]
        return value, value, value, "lower_bound"
    return None, None, None, event.get("status", "missing")


def render(records, receipts, output, *, organization, final=False, completion_review=None,
           relative_activation=None, relative_receipt=None,
           calibration_unavailable=None, unavailable_receipt=None,
           reporter_events=None, reporter_definition=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    output = Path(output).resolve()
    calibration_unavailable = calibration_unavailable or {}
    report = coverage(records, receipts, completion_review, calibration_unavailable=calibration_unavailable)
    relative_activation = relative_activation or {}
    if final and not report["final_export_eligible"]:
        raise ValueError("Final export blocked by incomplete independent scientific coverage")
    if final and reporter_definition is None:
        raise ValueError("Final panelB requires the reviewed frozen reporter-initiation definition")
    matching_relative_reporter = bool(reporter_definition and relative_receipt and
        relative_receipt.get("reference_fraction") == reporter_definition["fraction"] and
        relative_receipt.get("fixed_reference_thresholds_sha256") == reporter_definition["fixed_reference_thresholds_sha256"])
    if final and not matching_relative_reporter:
        raise ValueError("Final panelB requires matching audited relative-width reporter sensitivity")
    output.mkdir(parents=True, exist_ok=True)
    data, vectors, metadata = (output / "supporting" / name for name in ("data", "vectors", "metadata"))
    for path in (data, vectors, metadata):
        path.mkdir(parents=True, exist_ok=True)
    org, org_receipt = verified_organization(organization)
    plt.rcParams.update({"font.size": 8.3, "axes.titlesize": 9.5, "axes.labelsize": 8.,
                         "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(12.1, 14.6), layout="constrained")
    grid = fig.add_gridspec(6, 1, height_ratios=(1.7, 1.65, 2.1, 1.65, 1.75, 1.6))
    state = "AUDITED FINAL" if final else "PROVISIONAL • independently audited completed arms only"
    fig.suptitle("Figure 2. Shape and timing distinguish replacement responses\n" + state,
                 fontsize=13, fontweight="bold")

    # A: the retained 97-state organization curves preserve their original units,
    # reference lengths, half-height region and all-five context.
    organization_grid = grid[0].subgridspec(1, 3)
    activation_grid = grid[1].subgridspec(1, 5)
    fate_grid = grid[3].subgridspec(1, 6, width_ratios=(1, 1, 1, 1, 1, .07))
    trace_grid = grid[4].subgridspec(1, 5)
    for j, (metric, title) in enumerate((("excess_amount", "Incumbent excess amount"),
                                        ("width", "Incumbent half-height extent"),
                                        ("peak", "Incumbent peak signal"))):
        ax = fig.add_subplot(organization_grid[0, j])
        for model in MODELS:
            rows = sorted((r for r in org if r["model"] == model), key=lambda r: float(r["size_native"]))
            ax.plot([float(r["size_over_reference"]) for r in rows],
                    [float(r[metric + "_over_reference"]) for r in rows], color=COLORS[model], label=NAMES[model])
        ax.set(title=("A  " if j == 0 else "") + title, xlabel="Native length / own reference", ylabel="Fold of own reference")
        ax.spines[["top", "right"]].set_visible(False)
        if j == 2:
            ax.legend(fontsize=6.8, loc="upper left", frameon=False)

    # B: fixed physical width, with actual commanded-area strengths. Native
    # clocks receive separate axes; activation is never called replacement.
    for j, model in enumerate(MODELS):
        ax = fig.add_subplot(activation_grid[0, j])
        for area, style in ((.5, ":"), (1., "-"), (2., "--")):
            points = []
            for length in LENGTHS:
                row = records.get(("width", model, length, 3.8, area))
                event = ([reporter_events[("width", model, length, 3.8, area)]]
                         if reporter_events and row else row["measurement"].get("activation_events", []) if row else [])
                if event:
                    y, lo, hi, kind = _event_value(event[0])
                    if y is not None:
                        if kind == "observed":
                            points.append((length, y))
                        ax.errorbar(length, y, yerr=[[y - lo], [hi - y]], fmt="^" if kind == "lower_bound" else "o",
                                    color=COLORS[model], markersize=3, capsize=1)
                    elif kind == "already_above_at_interval_start":
                        ax.plot(length, 0, "D", markerfacecolor="none", color=COLORS[model], markersize=4)
            if points:
                ax.plot(*zip(*points), style, color=COLORS[model], linewidth=1, label=f"{area:g}× area")
        ylabel = "Reporter initiation time (native)" if reporter_definition else "Half-basal height time (native)"
        ax.set(title=("B  " if j == 0 else "") + NAMES[model], xlabel="Length label L", ylabel=ylabel,
               yscale="linear" if reporter_definition else "log")
        ax.set_xticks([19, 29, 39])
        for area, style in (() if reporter_definition and not matching_relative_reporter else
                            ((.5, ":"), (1., "-"), (2., "--"))):
            points = []
            for length in LENGTHS:
                relative = relative_activation.get((model, length, area))
                events = relative["measurement"].get("activation_events", []) if relative else []
                if events:
                    y, lo, hi, kind = _event_value(events[0])
                    if y is not None:
                        ax.plot(length, y, "v" if kind == "lower_bound" else "s", markerfacecolor="none",
                                color=".35", markersize=3)
                        if kind == "observed":
                            points.append((length, y))
                    elif kind == "already_above_at_interval_start":
                        ax.plot(length, 0, "D", markerfacecolor="none", color=".35", markersize=4)
            if points:
                ax.plot(*zip(*points), style, color=".4", linewidth=.7,
                        label="20% L, same peak" if area == 1 else None)
        if j == 0:
            ax.legend(fontsize=6.2, frameon=False)
        if reporter_definition:
            # Apply limits after fixed- and relative-width overlays; fixing the
            # lower limit earlier would disable autoscaling and clip Goryachev.
            ax.relim()
            ax.autoscale_view()
            _, top = ax.get_ylim()
            ax.set_ylim(0, top * (1.55 if j == 0 else 1.02))

    ratios = ratio_records(records)
    middle = grid[2].subgridspec(1, 2, width_ratios=(1.1, 1))
    ax = fig.add_subplot(middle[0, 0])
    for model in MODELS:
        points = []
        for row in (r for r in ratios if r["model"] == model):
            r, length = row["ratio"], row["length"]
            if r.get("lower") is None or not r.get("informative", True):
                continue
            lo, hi = r["lower"], r.get("upper")
            if hi is None and lo > 0:
                ax.plot(length, lo, "^", color=COLORS[model], markersize=5)
                ax.annotate("", (length, lo * 1.32), (length, lo), arrowprops=dict(arrowstyle="->", color=COLORS[model]))
            elif hi is not None and lo > 0:
                y = math.sqrt(lo * hi)
                ax.errorbar(length, y, yerr=[[y - lo], [hi - y]], color=COLORS[model], fmt="o", capsize=2, markersize=4)
                points.append((length, y))
            elif hi is not None:
                ax.plot(length, hi, "v", color=COLORS[model])
        if points:
            ax.plot(*zip(*points), color=COLORS[model], label=NAMES[model])
    ax.axhline(1, color=".5", linestyle=":", linewidth=.8)
    ax.set(title="C  Equal commanded area: narrow / broad", xlabel="Length label L", ylabel="Replacement-time ratio", yscale="log")
    current_bottom, current_top = ax.get_ylim()
    ax.set_ylim(current_bottom, current_top * 1.35)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    profile = fig.add_subplot(middle[0, 1])
    profile_sources = []
    for width, style in ((.95, "--"), (7.6, "-")):
        row = records.get(("width", "wave_pinning", 39, width, 1.))
        if row:
            trial, folder = audit_worker(row)
            with np.load(folder / "trace.npz", allow_pickle=False) as native:
                # Same predeclared native time, 1500, in both shape arms. No
                # matching by outcome or independently rescaled concentration.
                time = native["time"]
                index = int(np.argmin(np.abs(time - 1500.)))
                profile.plot(native["x"], native["active"][index], style, color=COLORS["wave_pinning"],
                             label=f"w={width:g}, t={time[index]:g}")
                profile_sources.append(dict(source=row["source"], width=width, sample_index=index,
                                            native_time=float(time[index]), trace_sha256=trial["full_trace_sha256"]))
    profile.set(title="Wave-pinning L39: same native time", xlabel="Native position", ylabel="Active concentration")
    if profile_sources:
        profile.legend(fontsize=7, frameon=False)

    # D: predeclared N=4, ON=1.5 calibration slice; ALL six gaps and lengths.
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#dedede")
    norm = Normalize(-1, 1)
    for j, model in enumerate(MODELS):
        ax = fig.add_subplot(fate_grid[0, j])
        values = np.full((len(GAPS), len(LENGTHS)), np.nan)
        for yi, gap in enumerate(GAPS):
            for xi, length in enumerate(LENGTHS):
                row = records.get(("train", model, length, 4, gap, 1.5))
                if ("train", model, length, 4, gap, 1.5) in calibration_unavailable:
                    ax.text(xi, yi, "†", ha="center", va="center", fontsize=8)
                if row:
                    m = row["measurement"]
                    if m["readout_eligible"] and m["final_d"] is not None:
                        values[yi, xi] = m["final_d"]
                    else:
                        ax.text(xi, yi, "×", ha="center", va="center", fontsize=8)
                    if m["field_tail_pass"] is False:
                        ax.plot(xi, yi, "s", fillstyle="none", color="black", markersize=9, markeredgewidth=.7)
        image = ax.imshow(values, origin="lower", cmap=cmap, norm=norm, aspect="auto")
        ax.set(title=("D  " if j == 0 else "") + NAMES[model], xlabel="Length label L")
        ax.set_xticks([0, 2, 5], [19, 27, 39])
        ax.set_yticks(range(len(GAPS)), [f"{g:g}" for g in GAPS])
        if j == 0:
            ax.set_ylabel("Gap / pulse duration")
    color_axis = fig.add_subplot(fate_grid[0, 5])
    fig.colorbar(image, cax=color_axis, label="Final D after recovery")

    # E: one fixed example in every model; normalized pulse time makes phase
    # positions comparable while preserving each actual native pulse duration.
    examples = []
    for j, model in enumerate(MODELS):
        ax = fig.add_subplot(trace_grid[0, j])
        selected_length = example_length(model, calibration_unavailable)
        row = records.get(("train", model, selected_length, 4, 4., 1.5))
        if row:
            trial, folder = audit_worker(row)
            measurement = read(folder / "readouts.json")
            native_trial = read(Path(row["source"]["run_root"]) / "workers" / f"{row['source']['task']:04d}" / "result.json")["trial"]
            pulse = float(native_trial["schedule"]["pulse_duration"])
            stop = float(native_trial["recovery_interval"][0])
            series = measurement["series"]
            keep = np.asarray(series["time"]) <= stop + pulse
            d = np.array([np.nan if v is None else v for v in series["D_half"]])
            ax.plot(np.asarray(series["time"])[keep] / pulse, d[keep], color=COLORS[model], linewidth=1)
            for start, end in native_trial["phase_on_intervals"]:
                ax.axvspan(start / pulse, end / pulse, color=COLORS[model], alpha=.15)
            endpoints = [r for r in trial["endpoints"] if r["kind"] in ("on", "gap")]
            for endpoint in endpoints:
                value = endpoint["end"].get("D_half")
                if value is not None:
                    ax.plot(endpoint["end"]["time"] / pulse, value,
                            "o" if endpoint["kind"] == "on" else "s", color=COLORS[model], markersize=3)
            ax.text(.02, .02, f"Pulse = {pulse:.3g} native", transform=ax.transAxes, fontsize=6.6)
            examples.append(dict(model=model, requested_length=39, displayed_length=selected_length,
                                 calibration_fallback=selected_length != 39,
                                 source=row["source"], schedule=native_trial["schedule"],
                                 final_d=trial["final_d"], field_tail_pass=trial["field_tail_pass"], endpoints=endpoints))
        else:
            unavailable = calibration_unavailable.get(("train", model, 39, 4, 4., 1.5))
            message = (f"† No sampled commitment\nthrough {unavailable['searched_on_cap']:g} ON\nTrain calibration unavailable"
                       if unavailable else "Awaiting audited\ntrain condition")
            ax.text(.5, .5, message, transform=ax.transAxes, ha="center", color=".4", fontsize=7)
        ax.axhline(.5, color=".5", linestyle=":", linewidth=.7)
        ax.set(title=("E  " if j == 0 else "") + NAMES[model] +
               (f" L{selected_length}" if selected_length is not None else "") +
               ("*" if model == "debelly" else ""),
               xlabel="Time / pulse duration", ylim=(-1.05, 1.05))
        if selected_length is not None and selected_length != 39:
            ax.text(.02, .88, "L39 calibration unavailable", transform=ax.transAxes, fontsize=6.5,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.8, pad=1))
        if j == 0:
            ax.set_ylabel("Activity dominance D")

    # F: measured signatures, rather than hard-coded model labels or asserted
    # causal mechanisms. All five models remain visible when evidence is absent.
    summary_rows = []
    for model in MODELS:
        ratio = {r["length"]: r["ratio"] for r in ratios if r["model"] == model}
        example = next((e for e in examples if e["model"] == model), None)
        gap = next((r for r in example["endpoints"] if r["kind"] == "gap"), None) if example else None
        pulse_D = None if gap is None else gap["start"].get("D_half")
        gap_D = None if gap is None else gap["end"].get("D_half")
        gap_text = f"{pulse_D:.2f} → {gap_D:.2f}" if pulse_D is not None and gap_D is not None else "unavailable"
        final_text = f"{example['final_d']:.2f}" if example and example["final_d"] is not None else "unavailable"
        initiation = "under review"
        if reporter_events:
            short = reporter_events[("width", model, 19, 3.8, 1.)]
            long = reporter_events[("width", model, 39, 3.8, 1.)]
            initiation = ratio_text(replacement_time_ratio(long, short))
            if model == "debelly" and short.get("lower_elapsed") == long.get("lower_elapsed") == 0:
                initiation = "both in (0, 0.5]"
        summary_rows.append([NAMES[model] + ("*" if model == "debelly" else ""),
                             f"L{example['displayed_length']}" if example else "unavailable", initiation,
                             ratio_text(ratio[19]) + " → " + ratio_text(ratio[39]), gap_text, final_text])
    ax = fig.add_subplot(grid[5])
    ax.axis("off")
    ax.set_title("F  Signatures tied to measurements; causal interpretation requires native-field support", loc="left")
    table = ax.table(cellText=summary_rows,
                     colLabels=["Model", "Trace L", "Initiation L39/L19", "Narrow / broad, L19 → L39", "Pulse 1 → gap 1: D", "Final D"],
                     cellLoc="center", bbox=[0, .20, 1, .73], colWidths=[.14, .07, .17, .25, .25, .12])
    table.auto_set_font_size(False)
    table.set_fontsize(7.4)
    ax.text(0, .09, "Shape: beam width at fixed input area   •   Timing: equal total ON pulse trains   •   Reporter: activity share after each gap",
            transform=ax.transAxes, fontsize=7.8)
    ax.text(0, .005, "↑ triangles: finite-cap ratio lower bounds; bounds cannot be ranked. †: no calibration; ×: low activity; □: failed field tail.",
            transform=ax.transAxes, fontsize=7.1)
    ax.text(0, -.075, "* deBelly later-gap timing and endpoint fields have numerical qualifications (S5). Native clocks differ; optical dose is uncalibrated.",
            transform=ax.transAxes, fontsize=7.1)

    prefix = "figure2_signature_assays" if final else "figure2_signature_assays_provisional"
    png = output / (prefix + ".png")
    fig.savefig(png, dpi=180, facecolor="white")
    for suffix in ("pdf", "svg"):
        fig.savefig(vectors / (prefix + "." + suffix), facecolor="white")
    plt.close(fig)
    generated = datetime.now(timezone.utc).isoformat()
    report.update(generated_utc=generated, status="final" if final else "provisional",
                  producer_state="stopped", renderer_source_sha256=sha(Path(__file__)),
                  retained_organization=org_receipt, selected_example=EXAMPLE,
                  native_profile_selection="Wave-pinning L39 area1 widths0.95/7.6; nearest saved frame to native t=1500 in each arm.",
                  separate_native_clock_per_model=True, pilot_numbers_used=False,
                  audits=receipts, relative_activation_overlay=relative_receipt,
                  calibration_unavailable_evidence=unavailable_receipt,
                  reporter_initiation_definition=reporter_definition,
                  panelB_current_measurement="fixed reporter threshold" if reporter_definition else "height maturation diagnostic; initiation definition under review")
    write(data / "width_ratio_sampled_bounds.json", ratios)
    write(data / "selected_train_endpoints.json", examples)
    write(metadata / "selected_profile_sources.json", profile_sources)
    write(metadata / "coverage.json", report)
    legend = dict(title="Figure 2. Shape and timing distinguish replacement responses",
                  status=report["status"], generated_utc=generated,
                  paragraphs=[
                      "A, Retained organization of all 97 qualified basal states. Excess amount is integrated inside the interpolated basal half-height core. Amount, extent and native peak signal are divided by each model's own reference values; Spring retains its physical domain and native resting-length convention.",
                      ("B, Reporter initiation: first cue-site positive-excess peak reaching " +
                       f"{100 * reporter_definition['fraction']:g}% of the exact L19 basal excess peak, held fixed within each model across all lengths, widths and strengths. "
                       if reporter_definition else
                       "B, Height-maturation diagnostic pending review of the initiation readout: first cue-site positive-excess peak reaching half that run's basal excess peak. This can remain below threshold despite replacement and is not interpreted as initiation. ") +
                      "The fixed physical-width raised-cosine cue has support 3.8 chemical units or 0.1 Spring units. Dotted, solid and dashed coloured curves have 0.5, 1 and 2 times the reference commanded spatial input area. Grey open-square curves retain the same peak at relative support 20% of length for the three focal models; this diagnostic changes spatial area with length. Bars are saved-sample crossing brackets. deBelly events in the first native update remain in the interval (0, 0.5]; substep interpolation is descriptive only. Open diamonds at zero denote a footprint already above the reporter threshold before cue application. Native clocks are shown separately.",
                      "C, Time to eligible D ≥ 0.5 during continuous cue ON for width 0.95 divided by width 7.6 (Spring widths 0.025 and 0.2), at the central fixed commanded spatial area equal to native reference peak times 3.8 (Spring times 0.1). Bars propagate both sampled crossing brackets; upward triangles and arrows are strict finite-cap lower bounds when the narrow arm remains below replacement at its native horizon. Censored bounds are not finite replacement times and cannot be ranked as measured ratios. The profile inset compares exact native L39 fields at the predeclared saved time nearest 1500 with one common concentration scale.",
                      "D, Final activity dominance after four equal pulses with total ON 1.5 times each model–length's postwithdrawal commitment calibration, then 10,000 native cue-free units (4000 for deBelly). Gap durations are multiples of one pulse. E, Four pulses with gaps four times one pulse at L39 where calibrated; otherwise the largest calibrated length is used. This deterministic fallback shows Holmes 3 at L27, explicitly labelled, because L39 has no sampled continuous-cue commitment through 100,000 ON units. Selection depends on calibration availability, not train fate. Shading marks cue ON; circles and squares mark pulse and gap endpoints. Time is divided by each model's actual pulse duration, printed in native units. F, The trace-length column preserves this exception beside the numerical signature summary. Shared responses at calibrated lengths do not establish a matched-length comparison at L39.",
                      "Dominance D is the positive active-signal excess in the challenger half minus that in the incumbent half, divided by their sum. The instantaneous profile minimum is the common floor. Eligible D = 0.5 means 75% challenger share; low-activity states require explicit amount and contrast gates and receive no invented D. A square marks a failed finite native-field tail independently of final fate. Grey dagger cells have no sampled continuous-cue commitment through the explicitly searched finite cap, so no calibration duration or train trajectory is assigned. Plain grey cells lack a completed independent audit; neither is censoring of an executed train.",
                      "Commitment calibration is the shortest sampled continuous ON duration, bracketed by native-duration bisection, that yields eligible D ≥ 0.5 after withdrawal and the declared model-native recovery period. The finite native-field tail is evaluated independently and does not redefine scalar replacement. Native clocks differ between models; no absolute cross-model time ranking is used.",
                      "The deBelly trace and summary carry an asterisk because later-gap commitment timing is numerically sensitive. At L39, the third-gap endpoint is D = −0.586 with the native 0.5-unit update and D = +0.190 in the separately qualified 0.25-unit companion. The first two gaps still collapse and final replacement agrees, but later-gap timing and field-tail status are not resolution invariant. The L39 401-cell train also has native endpoint-field sensitivity despite phase-D and fate agreement. S5 retains exact-native-start and separately qualified quarter-step preparations, every endpoint, native-field errors and all gate outcomes. Noisy deBelly ensembles use the native 0.5-unit update and remain conditional on that numerical scheme.",
                      "All five models keep their native kinetics. Commanded spatial input area per ON time is not calibrated photon dose or realized activated molecular amount. A width penalty or gap retention is an observable signature, not a causal allocation to a particular feedback term or unique molecular identification. Holmes 3 and Spring determine whether supported signatures are shared. Polarity replacement is distinct from observed migration and conditional OU persistence estimates."
                  ])
    noise_display = metadata / "S3_display_audit.json"
    if noise_display.is_file():
        noise_receipt = read(noise_display)
        count_path = data / "figure2_S3_retention_counts.csv"
        if (noise_receipt.get("passed") is not True or
                noise_receipt.get("native_trajectories") != 486 or
                noise_receipt["artifact_sha256"].get(str(count_path.relative_to(output))) != sha(count_path)):
            raise ValueError("Main legend requires the complete verified S3 count table")
        with count_path.open(newline="") as stream:
            count_rows = list(csv.DictReader(stream))
        descriptions = []
        for model in ("wave_pinning", "debelly"):
            for length in (19, 29, 39):
                selected = sorted((r for r in count_rows if r["model"] == model and int(r["length"]) == length
                    and r["kind"] == "noise" and float(r["gap_multiple"]) == 16
                    and int(r["gap_index"]) == 3 and float(r["threshold"]) == .8),
                    key=lambda r: float(r["mean_fraction"]))
                if len(selected) != 2:
                    raise ValueError("S3 illustrative condition coverage changed")
                descriptions.append(f"{NAMES[model]} L{length}: " + " and ".join(
                    f"{r['successes']}/{r['eligible']}" for r in selected))
        legend["paragraphs"].append("S3 reports conditional retention through weakly driven gaps, separately from final replacement. At the third gap with gap duration 16 times one pulse, using the predeclared 80% raw-share criterion with D ≥ 0, counts at mean input 0.5% and 2% of native peak are " +
            "; ".join(descriptions) + ". These are independent-seed denominators, with exact 95% intervals shown in S3; the three gaps of one trajectory are not independent replicates. Mean and zero controls remain visible. Weak describes instantaneous bounded input (maximum 1% or 4% of native peak); at L39 with gap multiple 16 and mean input 2%, whole-domain gap exposure is about 2.46 times targeted ON exposure. This commanded native-coordinate coefficient is not photon dose or realized activation. deBelly retention remains conditional on the native 0.5-unit update.")
    write(metadata / "legend.json", legend)
    if reporter_definition:
        write(metadata / "reporter_initiation_definition.json", dict(
            **reporter_definition,
            observable="Maximum native active polarity signal inside actual cue support minus instantaneous whole-profile minimum.",
            threshold_rule=f"{100 * reporter_definition['fraction']:g}% of each model's exact L19 basal excess peak, fixed across length, width and strength.",
            baseline_active_policy="If already above threshold at cue onset, record baseline-active, never a fabricated activation time.",
            event_rule="First eligible saved crossing, retaining the preceding/current sampled bracket; interpolation is descriptive only.",
            debelly_native_clock="0.5-unit updates; first-step activation remains (0,0.5].",
            sensitivity="Full1/5/20% alternatives and original half-per-basal height-maturation diagnostic retained.",
            calibration_scope="Native signal units, not photon dose or biological reporter calibration."))
    repo = next((p for p in Path(__file__).resolve().parents
                 if (p / "polarity1d/manuscript").is_dir() and (p / "manuscript_outline").is_dir()),
                Path.cwd().resolve())
    rel = lambda p: str(p.relative_to(repo)) if p.is_relative_to(repo) else str(p)
    content = dict(id="figure2", title=legend["title"],
                   role="Test whether shape and timing of a measured challenge distinguish explanations for size-dependent replacement resistance.",
                   status=report["status"],
                   assets=[dict(path=rel(png), caption="Five-model shape and timing signature assays",
                                page_size=[11.8, 16.1], minimum_printed_width_inches=10.0,
                                page_title=legend["title"])],
                   legend=legend["paragraphs"],
                   narrative=["Narrow-to-broad replacement-time ratios at L19 and L39 are " +
                              "; ".join(r[0] + ": " + r[3] for r in summary_rows) +
                              ". These are sampled bounds, not fitted scaling laws.",
                              "Activity-share trajectories and their exact gap endpoints test whether an early challenged state persists through cue-free gaps. A later retained replacement does not establish that early activity itself persisted."],
                   measurements=["Vary physical cue support at fixed commanded spatial input area and measure time to eligible75% challenger activity share.",
                                 "Apply equal-total-ON pulse trains, measure activity share at every pulse and gap endpoint, then report postwithdrawal fate and field stability independently."],
                   limitations=["Native-to-optical calibration is absent. Native clocks, reporter species, amount scales and stimulus channels differ across models.",
                                "Observable signatures do not uniquely identify a molecular mechanism or allocate causation to one feedback term. All five models determine whether responses are shared.",
                                "deBelly early activity collapse and final replacement survive the qualified quarter-step check, but the third-gap dominance sign and finite-field-tail status change; the L39 refined-mesh train also has endpoint-field sensitivity. Later commitment timing and native-step noise probabilities are not claimed to be numerically invariant.",
                                "The current export remains provisional until all named native assays, supporting sensitivity checks, numerical companions and independent source/display audits are complete." if not final else
                                "Finite native-field tails qualify only the declared horizon; deterministic activity replacement does not measure migration persistence."],
                   sources=[dict(label="Independent coverage and scientific source receipts", path=rel(metadata / "coverage.json")),
                            dict(label="Sampled width-ratio bounds", path=rel(data / "width_ratio_sampled_bounds.json")),
                            dict(label="Selected pulse/gap endpoint measurements", path=rel(data / "selected_train_endpoints.json"))])
    write(metadata / "content.json", content)
    if (data / "figure2_S5_debelly_temporal_comparisons.json").is_file():
        content["sources"].append(dict(label="Material deBelly temporal qualification and all phase gates",
                                       path=rel(data / "figure2_S5_debelly_temporal_comparisons.json")))
        write(metadata / "content.json", content)
    for filename, label in (("S3_display_audit.json", "Complete bounded-noise retention and finite-seed uncertainty"),
                            ("S5_combined_display_audit.json", "All numerical companions and material qualifications")):
        if (metadata / filename).is_file():
            content["sources"].append(dict(label=label, path=rel(metadata / filename)))
    write(metadata / "content.json", content)
    artifacts = [png, *vectors.glob(prefix + ".*"), data / "width_ratio_sampled_bounds.json",
                 data / "selected_train_endpoints.json", metadata / "selected_profile_sources.json",
                 metadata / "coverage.json", metadata / "legend.json", metadata / "content.json"]
    if reporter_definition:
        artifacts.append(metadata / "reporter_initiation_definition.json")
    paras = "".join("<p>" + html.escape(p) + "</p>" for p in legend["paragraphs"])
    supporting_links = "".join("<li><a href='" + html.escape(str(p.relative_to(output)), quote=True) + "'>" +
        html.escape(p.stem.replace("figure2_", "").replace("_", " ")) + "</a></li>"
        for p in sorted(vectors.glob("figure2_S*.pdf")))
    (output / "index.html").write_text("<!doctype html><meta charset='utf-8'><title>Figure2 signature assays</title>"
        "<style>body{max-width:1100px;margin:2em auto;font:16px/1.5 system-ui}img{width:100%}pre{white-space:pre-wrap}</style>"
        "<h1>" + html.escape(legend["title"]) + "</h1><p><strong>" + state + "</strong></p>"
        "<p>Generated " + generated + "; renderer stopped. Opening or manually refreshing reads this saved report.</p>"
        "<p>Audited widths " + str(report["audited_width_conditions"]) + "/450; audited trains "
        + str(report["audited_train_conditions"]) + " native trains plus "
        + str(report["calibration_unavailable_train_conditions"]) + " explicitly calibration-unavailable conditions; "
        + str(report["scientifically_accounted_train_conditions"]) + "/1620 planned train conditions accounted for. "
        "<a href='supporting/metadata/coverage.json'>Full coverage record</a>.</p>"
        "<img src='" + png.name + "' alt='Figure2 A–F audited signature assays'>"
        "<h2>Supporting assays</h2><ul>" + supporting_links + "</ul>"
        "<div style='display:grid;grid-template-columns:1fr auto;align-items:center'>"
        "<math display='block'><mrow><mi>D</mi><mo>=</mo><mfrac><mrow><msub><mi>E</mi><mi>C</mi></msub><mo>−</mo><msub><mi>E</mi><mi>I</mi></msub></mrow>"
        "<mrow><msub><mi>E</mi><mi>C</mi></msub><mo>+</mo><msub><mi>E</mi><mi>I</mi></msub></mrow></mfrac></mrow></math><span>(1)</span></div>"
        "<div style='display:grid;grid-template-columns:1fr auto;align-items:center'><math display='block'><mrow>"
        "<msub><mi>E</mi><mi>j</mi></msub><mo>(</mo><mi>t</mi><mo>)</mo><mo>=</mo>"
        "<msub><mo>∫</mo><msub><mi>Ω</mi><mi>j</mi></msub></msub>"
        "<msub><mrow><mo>[</mo><mi>u</mi><mo>(</mo><mi>x</mi><mo>,</mo><mi>t</mi><mo>)</mo><mo>−</mo>"
        "<munder><mo>min</mo><mi>y</mi></munder><mi>u</mi><mo>(</mo><mi>y</mi><mo>,</mo><mi>t</mi><mo>)</mo><mo>]</mo></mrow><mo>+</mo></msub>"
        "<mi>d</mi><mi>x</mi></mrow></math><span>(2)</span></div>"
        "<p>Here Ω denotes the native challenger or incumbent half, and u is the model's active polarity signal. Amount is a native one-dimensional signal integral; molecular counts and optical intensities are uncalibrated.</p>" + paras)
    artifacts.append(output / "index.html")
    write(metadata / "artifact_sha256.json", {str(p.relative_to(output)): sha(p) for p in artifacts})
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, help="Saved audit/source selection, without launching native work")
    parser.add_argument("--analysis-dir", action="append", default=[], type=Path)
    parser.add_argument("--organization", type=Path)
    parser.add_argument("--relative-activation-audit", action="append", default=[], type=Path)
    parser.add_argument("--unavailable-trains", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--completion-review", type=Path)
    parser.add_argument("--final", action="store_true")
    parser.add_argument("--reporter-sensitivity", action="store_true")
    parser.add_argument("--reporter-sensitivity-dir", type=Path)
    parser.add_argument("--activation-fraction", type=float, choices=[.01, .05, .1, .2])
    parser.add_argument("--relative-reporter-index", type=Path)
    parser.add_argument("--relative-reporter-sensitivity-dir", type=Path)
    parser.add_argument("--plot-reporter-sensitivity", type=Path)
    parser.add_argument("--scientific-master", type=Path)
    args = parser.parse_args()
    if args.scientific_master:
        if not args.inputs: parser.error("--scientific-master requires --inputs")
        print(json.dumps(scientific_master(args.inputs,args.scientific_master),indent=2))
        raise SystemExit(0)
    if args.inputs:
        config = read(args.inputs)
        if not args.analysis_dir:
            args.analysis_dir = list(map(Path, config.get("analysis_directories", [])))
        for name in ("organization", "output", "reporter_sensitivity_dir", "relative_reporter_sensitivity_dir", "unavailable_trains"):
            if getattr(args, name) is None and config.get(name):
                setattr(args, name, Path(config[name]))
        if args.activation_fraction is None:
            args.activation_fraction = config.get("activation_fraction")
    if args.output is None:
        parser.error("--output or a saved inputs output is required")
    if args.relative_reporter_index:
        if args.reporter_sensitivity_dir is None:
            parser.error("Relative sensitivity requires the fixed reporter sensitivity reference")
        print(json.dumps(relative_reporter_sensitivity(args.relative_reporter_index,
            args.reporter_sensitivity_dir, args.output), indent=2))
        raise SystemExit(0)
    if args.plot_reporter_sensitivity:
        print(json.dumps(plot_reporter_sensitivity(args.plot_reporter_sensitivity, args.output), indent=2))
        raise SystemExit(0)
    records, receipts = load_audited_records(args.analysis_dir)
    if args.reporter_sensitivity:
        print(json.dumps(reporter_threshold_sensitivity(records, args.output), indent=2))
        raise SystemExit(0)
    if args.organization is None:
        parser.error("--organization is required for a figure export")
    relative, relative_receipt = load_relative_activation(args.relative_activation_audit)
    unavailable, unavailable_receipt = load_calibration_unavailable(args.unavailable_trains)
    reporter_events, reporter_definition = load_reporter_definition(args.reporter_sensitivity_dir, args.activation_fraction)
    if args.relative_reporter_sensitivity_dir:
        relative, relative_receipt = load_relative_reporter_definition(args.relative_reporter_sensitivity_dir, args.activation_fraction)
    report = render(records, receipts, args.output, organization=args.organization,
                    final=args.final, completion_review=args.completion_review,
                    relative_activation=relative, relative_receipt=relative_receipt,
                    calibration_unavailable=unavailable, unavailable_receipt=unavailable_receipt,
                    reporter_events=reporter_events, reporter_definition=reporter_definition)
    print(json.dumps({k: report[k] for k in ("status", "audited_width_conditions", "audited_train_conditions", "final_export_eligible")}, indent=2))
