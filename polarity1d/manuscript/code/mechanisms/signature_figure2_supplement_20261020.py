"""Complete S1--S5 display assembly from audited Figure 2 records; no solvers.

Full source time series remain in the immutable native/readout and audit roots.
Display reduction retains every phase endpoint and each time-bin extremum.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import math
import json
from pathlib import Path
import textwrap
import shutil

import numpy as np

from signature_figure2_render_20261020 import (
    AREAS, GAPS, LENGTHS, MODELS, NAMES, ON_FACTORS, PULSES, WIDTHS,
    _event_value, audited_worker_metadata, coverage, load_audited_records,
    load_calibration_unavailable, read, sha, write,
)


def display_indices(time, values, endpoints=(), bins=500):
    """Preserve all extrema and named phase endpoints while reducing artwork."""
    time, values = np.asarray(time), np.asarray(values, dtype=float)
    if len(time) != len(values) or not len(time):
        raise ValueError("Aligned nonempty display series required")
    if len(time) <= 2 * bins:
        return np.arange(len(time))
    selected = {0, len(time) - 1}
    for chunk in np.array_split(np.arange(len(time)), bins):
        if len(chunk):
            selected.add(int(chunk[np.nanargmin(values[chunk])]))
            selected.add(int(chunk[np.nanargmax(values[chunk])]))
    for endpoint in endpoints:
        index = int(np.argmin(np.abs(time - endpoint)))
        if abs(time[index] - endpoint) > 1e-9:
            raise ValueError("Named phase endpoint is missing from source display series")
        selected.add(index)
    return np.array(sorted(selected))


def exact_binomial_interval(successes, eligible, alpha=.05):
    """Two-sided finite-seed interval; three gaps never become three trials."""
    from scipy.stats import beta
    if eligible == 0:
        return None, None
    if not 0 <= successes <= eligible:
        raise ValueError("Invalid independent-trajectory retention count")
    return (0. if successes == 0 else float(beta.ppf(alpha / 2, successes, eligible - successes + 1)),
            1. if successes == eligible else float(beta.ppf(1 - alpha / 2, successes + 1, eligible - successes)))


def _csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build_noise(summary_path, output):
    """S3 from all 486 independently audited trajectories, including controls."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    summary_path, output = Path(summary_path).resolve(), Path(output).resolve()
    root = summary_path.parent.parent
    proof_path = summary_path.parent / "independent_native_readout_audit.json"
    summary, proof = read(summary_path), read(proof_path)
    if (proof.get("passed") is not True or proof.get("workers") != 486 or
            proof.get("noise_summary_sha256") != sha(summary_path) or
            proof.get("manifest_sha256") != sha(root / "manifest.json") or
            not proof.get("all_native_phase_join_audits_pass") or
            not proof.get("all_readout_recomputations_pass") or
            not summary.get("full_ensemble_complete") or
            not summary.get("all_density_conservation_guards_pass") or
            not summary.get("all_hard_input_bounds_pass")):
        raise ValueError("S3 requires the complete independently audited native/noise ensemble")
    by_task = {int(r["task"]): r for r in summary["trajectories"]}
    rows = {(int(r["task"]), int(r["gap_index"])): r for r in summary["rows"]}
    if len(by_task) != 486 or len(rows) != 1458 or set(by_task) != set(range(486)):
        raise ValueError("S3 trajectory/gap coverage is incomplete or duplicated")
    for task, trajectory in by_task.items():
        receipt = proof["worker_analysis_sha256"][str(task)]
        audit_path = root / receipt["path"]
        if sha(audit_path) != receipt["sha256"] or sha(trajectory["result_path"]) != trajectory["result_sha256"]:
            raise ValueError("S3 audited worker or native result changed")
        audit = read(audit_path)
        trial = next(t for t in audit["trials"] if t["label"] == "trial")
        gaps = [p for p in trial["endpoints"] if p["kind"] == "gap"]
        if len(gaps) != 3:
            raise ValueError("S3 must preserve all three exact gap endpoints")
        for index, gap in enumerate(gaps, 1):
            row = rows[task, index]
            for source, destination in (("primary_retention_pulse_eligible", "pulse_eligible"),
                    ("primary_retention_denominator_eligible", "denominator_eligible"),
                    ("primary_retention_ratio", "retention_ratio"), ("primary_retention_success", "success"),
                    ("primary_retention_success_by_threshold", "sensitivity")):
                if row[destination] != gap[source]:
                    raise ValueError("S3 summary differs from independently reconstructed native endpoint")
        if trajectory["final_d"] != trial["final_d"] or trajectory["final_field_tail_pass"] != trial["field_tail_pass"]:
            raise ValueError("S3 final fate/tail summary differs from independent native audit")
    folders = {name: output / "supporting" / name for name in ("vectors", "data", "metadata")}
    for folder in folders.values(): folder.mkdir(parents=True, exist_ok=True)
    summary_copy = folders["data"] / "figure2_S3_all_native_trajectories.json"
    shutil.copy2(summary_path, summary_copy)
    grouped = []
    for model in ("wave_pinning", "debelly"):
        for length in (19, 29, 39):
            for gap in (1., 4., 16.):
                for kind, mean, expected in (("noise", .005, 12), ("noise", .02, 12),
                                             ("mean", .005, 1), ("mean", .02, 1), ("zero", 0., 1)):
                    for index in (1, 2, 3):
                        selected = [r for r in rows.values() if (r["model"], r["reference_length"], r["gap_multiple"],
                            r["kind"], r["mean_fraction"], r["gap_index"]) == (model, length, gap, kind, mean, index)]
                        if len(selected) != expected or len({r["seed"] for r in selected}) != expected:
                            raise ValueError("S3 design coverage or independent seed identity failed")
                        for threshold in (.7, .8, .9):
                            eligible = [r for r in selected if r["denominator_eligible"]]
                            count = sum(r["sensitivity"][str(threshold)] is True for r in eligible)
                            low, high = exact_binomial_interval(count, len(eligible)) if kind == "noise" else (None, None)
                            grouped.append(dict(model=model, length=length, gap_multiple=gap, kind=kind,
                                mean_fraction=mean, gap_index=index, threshold=threshold, independent_seeds=expected,
                                eligible=len(eligible), successes=count,
                                probability=count / len(eligible) if eligible else None, ci_lower=low, ci_upper=high,
                                uncertainty="95% exact binomial, conditional on observed eligibility" if kind == "noise" else "deterministic control; no seed CI"))
    table_path = folders["data"] / "figure2_S3_retention_counts.csv"
    _csv(table_path, grouped)
    pdf_path = folders["vectors"] / "figure2_S3_bounded_noise_retention.pdf"
    colors = {0.: ".25", .005: "#2878b5", .02: "#cb6635"}
    pages = 0
    with PdfPages(pdf_path) as pdf:
        for threshold in (.7, .8, .9):
            for index in (1, 2, 3):
                fig, axes = plt.subplots(2, 3, figsize=(11.8, 7.8), layout="constrained")
                fig.suptitle(f"S3  Gap {index}: retain ≥{threshold:.0%} of pulse-end challenger share and D≥0", fontsize=14)
                for row_index, model in enumerate(("wave_pinning", "debelly")):
                    for column, length in enumerate((19, 29, 39)):
                        ax = axes[row_index, column]
                        selected = [r for r in grouped if (r["model"], r["length"], r["gap_index"], r["threshold"]) ==
                            (model, length, index, threshold)]
                        for kind, mean in (("noise", .005), ("noise", .02), ("mean", .005), ("mean", .02), ("zero", 0.)):
                            series = sorted([r for r in selected if r["kind"] == kind and r["mean_fraction"] == mean], key=lambda r: r["gap_multiple"])
                            offset = -.10 if mean == .005 else .10 if mean == .02 else 0.
                            x = np.arange(3) + offset
                            y = [np.nan if r["probability"] is None else r["probability"] for r in series]
                            ax.plot(x, y, marker="o" if kind == "noise" else "s" if kind == "mean" else "x",
                                color=colors[mean], markerfacecolor=colors[mean] if kind == "noise" else "none",
                                ls="-" if kind == "noise" else "--" if kind == "mean" else ":", lw=.8, ms=4,
                                label=f"{mean:g}× peak, {kind}")
                            for xx, yy, record in zip(x, y, series):
                                if kind == "noise":
                                    if record["eligible"]:
                                        ax.vlines(xx, record["ci_lower"], record["ci_upper"], color=colors[mean], lw=.8)
                                        ax.annotate(f"{record['successes']}/{record['eligible']}", (xx, yy),
                                            xytext=(0, 7 if mean == .005 else -14), textcoords="offset points", ha="center", fontsize=6)
                                    else:
                                        ax.text(xx, .53 if mean == .005 else .43, "0/12 eligible", color=colors[mean], ha="center", fontsize=6)
                        ax.set(title=f"{NAMES[model]} L{length}", xlabel="Gap / pulse duration", ylim=(-.08, 1.14), ylabel="Conditional retention fraction")
                        ax.set_xticks(range(3), ("1", "4", "16")); ax.axhline(0, color=".7", lw=.5)
                axes[0, 0].legend(fontsize=6, frameon=False, loc="center left")
                fig.supxlabel(textwrap.fill("Each noisy denominator contains at most 12 independent trajectories; three gaps from one trajectory are not independent replicates. Eligibility requires pulse-end challenger bias, not eventual replacement. Known gap-end bias/activity loss counts as failure. Error bars: exact 95% finite-seed intervals. deBelly probabilities are conditional on its native0.5 update; later-gap timing is sensitive in S5.", 165), fontsize=7)
                pdf.savefig(fig); plt.close(fig); pages += 1
        fig, axes = plt.subplots(2, 3, figsize=(11.8, 7.8), layout="constrained")
        fig.suptitle("S3  Bounded instantaneous input can accumulate substantial gap exposure", fontsize=14)
        for i, model in enumerate(("wave_pinning", "debelly")):
            for j, length in enumerate((19, 29, 39)):
                ax = axes[i, j]
                for mean in (.005, .02):
                    selected = [r for r in by_task.values() if (r["model"], r["reference_length"], r["kind"], r["mean_fraction"]) == (model, length, "noise", mean)]
                    for gap in (1., 4., 16.):
                        vals = [r["gap_exposure_over_on_exposure"] for r in selected if r["gap_multiple"] == gap]
                        ax.plot([gap] * len(vals), vals, "o", color=colors[mean], ms=3, alpha=.6)
                ax.set(title=f"{NAMES[model]} L{length}", xlabel="Gap / pulse duration", ylabel="Integrated gap / ON coefficient", xscale="log")
                ax.set_xticks((1, 4, 16), ("1", "4", "16")); ax.axhline(1, color=".5", ls=":", lw=.8)
        fig.supxlabel(textwrap.fill("Weak describes the instantaneous bounded coefficient (maximum1% or4% of native peak), not cumulative exposure. Spatially distributed gap input is integrated over the whole native domain; at L39/gap16/mean2% its expected commanded exposure is about2.46× targeted ON exposure. This is not calibrated photon dose or biochemical activation. Native guards, all exposure summaries, final fate and field-tail flags remain separately in the source JSON.", 165), fontsize=8)
        pdf.savefig(fig); plt.close(fig); pages += 1
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True, pages=pages,
        native_trajectories=486, independent_noise_trajectories=432, deterministic_mean_controls=36,
        deterministic_zero_controls=18, exact_gap_records=1458, summary_sha256=sha(summary_path),
        native_audit_path=str(proof_path), native_audit_sha256=sha(proof_path),
        retention_is_conditional_not_final_fate=True, field_tail_filters_retention=False,
        uncertainty="Independent-seed exact binomial intervals; no pooled gaps",
        native_integrations_run=0, source_sha256=sha(Path(__file__)),
        artifact_sha256={str(p.relative_to(output)): sha(p) for p in (pdf_path, summary_copy, table_path)})
    write(folders["metadata"] / "S3_display_audit.json", report)
    return report


def build_numerics(summary_path, output):
    """Display all mesh/time-step comparisons, including known qualifications."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    summary_path, output = Path(summary_path).resolve(), Path(output).resolve()
    summary = read(summary_path)
    if not summary.get("independent_native_identity_audits_passed") or summary.get("actual_protocol_comparisons") != 41:
        raise ValueError("S5 requires all 41 independent numerical protocol comparisons")
    records = summary["records"]
    if len(records) != 41 or len({(r["run"], r["task"], r["variant"]) for r in records}) != 41:
        raise ValueError("S5 numerical comparison coverage is incomplete or duplicated")
    for source in summary["sources"]:
        root = Path(source["run"])
        audit_path = root / "results/independent_protocol_audit.json"
        if sha(root / "manifest.json") != source["manifest_sha256"] or sha(audit_path) != source["native_audit_sha256"]:
            raise ValueError("S5 source manifest/native audit changed")
        if not read(audit_path).get("passed"):
            raise ValueError("S5 native identity audit failed")
    folders = {name: output / "supporting" / name for name in ("vectors", "data", "metadata")}
    for folder in folders.values(): folder.mkdir(parents=True, exist_ok=True)
    copies = []
    for filename, digest in summary["outputs"].items():
        original = summary_path.parent / filename
        if sha(original) != digest:
            raise ValueError("S5 comparison CSV identity changed")
        copied = folders["data"] / ("figure2_S5_combined_" + filename)
        shutil.copy2(original, copied); copies.append(copied)
    copied_summary = folders["data"] / "figure2_S5_combined_comparisons.json"
    shutil.copy2(summary_path, copied_summary); copies.append(copied_summary)
    flat = []
    for r in records:
        flat.append(dict(model=r["model"], length=r["requested_reference_length"], assay=r["assay"], variant=r["variant"],
            companion_no_cue_control_pass=r["companion_no_cue_control_pass"],
            maximum_endpoint_D_difference=r["phase_comparison"]["maximum_absolute_eligible_D_difference"],
            maximum_endpoint_native_field_scaled_difference=r["endpoint_native_fields"]["maximum_scaled_difference"],
            first_replacement_event_comparison=r["first_replacement_event_comparison"]["classification"],
            main_final_D=r["main201_final"]["final_d"], companion_final_D=r["companion_final_D"],
            main_tail=r["main201_final"]["field_tail_pass"], companion_tail=r["companion_field_tail_pass"],
            all_requested_readouts_agree=r["all_declared_readout_comparisons_agree"],
            issues="; ".join(r["issues"]), source_run=r["run"], source_task=r["task"],
            native_worker_result_sha256=r["native_worker_result_sha256"]))
    flat_path = folders["data"] / "figure2_S5_all_comparison_gates.csv"
    _csv(flat_path, flat); copies.append(flat_path)
    pdf_path = folders["vectors"] / "figure2_S5_all_numerical_qualifications.pdf"
    pages = 0
    with PdfPages(pdf_path) as pdf:
        for model in MODELS:
            selected = [r for r in flat if r["model"] == model]
            fig, ax = plt.subplots(figsize=(13.5, 9.5)); ax.axis("off")
            ax.set_title(f"S5  {NAMES[model]}: all mesh/time-step comparison gates", fontsize=14, loc="left")
            def value(x): return "unavailable" if x is None else f"{x:.4g}"
            issue_codes = {"companion_no_cue_control_failed": "NC", "finite_field_tail_status_changed": "tail",
                "phase_endpoint_eligibility_changed": "eligibility", "phase_endpoint_D_difference_above_tolerance": "D",
                "native_endpoint_field_difference_above_tolerance": "fields", "phase_endpoint_challenger_bias_changed": "bias"}
            cells = [[r["length"], r["assay"].replace("train_4_gap4_on1p5", "4 pulses, gap4, ON1.5τ"),
                r["variant"].replace("_", " "), str(r["companion_no_cue_control_pass"]),
                value(r["maximum_endpoint_D_difference"]), value(r["maximum_endpoint_native_field_scaled_difference"]),
                value(r["main_final_D"]) + " / " + value(r["companion_final_D"]),
                f"{r['main_tail']} / {r['companion_tail']}",
                ", ".join(issue_codes.get(code, code) for code in r["issues"].split("; ")) if r["issues"] else "none"] for r in selected]
            table = ax.table(cellText=cells, colLabels=["L", "Assay", "Companion", "No-cue\ngate", "Max phase\n|ΔD|",
                "Max native\nfield error", "Final D:\nmain / companion", "Tail:\nmain / companion", "Qualification codes"],
                colWidths=[.035, .14, .13, .075, .085, .09, .115, .105, .225],
                cellLoc="center", bbox=[0, .21, 1, .69])
            table.auto_set_font_size(False); table.set_fontsize(6.7 if len(selected) > 10 else 8)
            message = "Codes: NC, companion no-cue gate failed; tail, finite-tail flag changed; eligibility, phase eligibility changed; D, phase |ΔD| exceeds tolerance; fields, native-field error exceeds tolerance; bias, challenger-bias sign changed. Full warning text is retained in CSV. Independent native identity/readout audits pass for every actual comparison. Numerical agreement is assessed separately for every requested readout; final fate agreement does not establish pulse/gap timing or full-field agreement. The mesh comparison is201 versus401 cells, using the same frozen native kinetics and realized schedules."
            if model == "debelly":
                message += " Both exact-native-start quarter-step diagnostics and separately qualified quarter-step basals are retained. The L39 mesh401 train has native-field endpoint sensitivity; qualified quarter-step trains preserve early bias loss and final fate but later-gap timing, eligibility and tail status differ."
            if model == "holmes_model3":
                message += " The requested L39 train has no sampled commitment calibration through100,000 ON; its numerical companion is explicitly unavailable. It is not replaced with a forced-cap schedule."
            ax.text(0, .17, textwrap.fill(message, 195), transform=ax.transAxes, fontsize=8.5, va="top")
            pdf.savefig(fig); plt.close(fig); pages += 1
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        pages=pages, independent_native_identity_audits_passed=True, actual_protocol_comparisons=41,
        numerical_agreement_for_all_requested_readouts=summary["numerical_agreement_for_all_requested_readouts"],
        unavailable_reference_cases=summary["unavailable_reference_cases"],
        source_summary=str(summary_path), source_summary_sha256=sha(summary_path),
        native_integrations_run=0, source_sha256=sha(Path(__file__)),
        artifact_sha256={str(p.relative_to(output)): sha(p) for p in [pdf_path, *copies]})
    write(folders["metadata"] / "S5_combined_display_audit.json", report)
    return report


def combine_supplements(output):
    """Join the complete display sections while retaining their source hashes."""
    from pypdf import PdfReader, PdfWriter
    output = Path(output).resolve()
    vectors = output / "supporting/vectors"
    metadata = output / "supporting/metadata"
    receipts = {}
    for name in ("S1_S2_display_audit.json", "S3_display_audit.json", "S4_display_audit.json",
                 "S5_combined_display_audit.json", "S5_debelly_temporal_display_audit.json"):
        path = metadata / name
        value = read(path)
        if value.get("passed") is False:
            raise ValueError("A scientific display section failed its completion gate")
        for relative, digest in value["artifact_sha256"].items():
            if sha(output / relative) != digest:
                raise ValueError("A supplement section changed after its display audit")
        receipts[name] = dict(path=str(path), sha256=sha(path))
    s1s2 = read(metadata / "S1_S2_display_audit.json")
    if (s1s2["S1"]["audited_width_histories"] != 450 or
            not s1s2["S2"]["complete_native_or_calibration_coverage"] or
            s1s2["S2"]["native_train_histories"] != 1458):
        raise ValueError("The main assay supplement grid is incomplete")
    reporter_path=metadata/"S4_reporter_sensitivity_display_audit.json"
    reporter=read(reporter_path)
    if (not reporter.get("passed") or reporter["sha256"]!=sha(reporter["path"]) or
            reporter["source_audit_sha256"]!=sha(Path(reporter["source_directory"])/"reporter_sensitivity_audit.json")):
        raise ValueError("Reporter sensitivity source/display identities changed")
    receipts[reporter_path.name]=dict(path=str(reporter_path),sha256=sha(reporter_path))
    names = ["figure2_S1_equal_area_widths.pdf", "figure2_S2_all_train_fate_maps.pdf",
        "figure2_S3_bounded_noise_retention.pdf", "figure2_S4_position_share_dwell.pdf",
        "figure2_S4_reporter_threshold_sensitivity.pdf", "figure2_S5_all_numerical_qualifications.pdf",
        "figure2_S5_debelly_temporal_qualification.pdf"]
    writer, sections = PdfWriter(), []
    for name in names:
        path = vectors / name
        document = PdfReader(path)
        start = len(writer.pages) + 1
        writer.append(document)
        sections.append(dict(path=str(path), sha256=sha(path), first_page=start, pages=len(document.pages)))
    target = vectors / "figure2_signature_assays_complete_supplement.pdf"
    with target.open("wb") as stream: writer.write(stream)
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        all_named_supplement_sections_complete=True, final_main_export_status="requires separate independent panel/stage review",
        pages=len(writer.pages), sections=sections, source_display_receipts=receipts,
        numerical_agreement_for_all_requested_readouts=False,
        pdf=dict(path=str(target), sha256=sha(target)), source_sha256=sha(Path(__file__)))
    write(metadata / "combined_supplement_audit.json", report)
    return report


def _verified_history(record):
    trial, folder, receipt = audited_worker_metadata(record)
    path = folder / "readouts.json"
    if sha(path) != trial["readouts_current_sha256"]:
        raise ValueError("S1 source readout identity changed")
    measurement = read(path)
    for key, value in trial["regional_amount_series"].items():
        if key not in measurement["series"] or value != measurement["series"][key]:
            raise ValueError("S1 native-audited history differs from original readout series")
    return trial, measurement, receipt


def build_supplements(records, receipts, unavailable, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import Normalize

    output = Path(output).resolve()
    vectors = output / "supporting/vectors"
    data = output / "supporting/data"
    metadata = output / "supporting/metadata"
    for folder in (vectors, data, metadata):
        folder.mkdir(parents=True, exist_ok=True)
    gate = coverage(records, receipts, calibration_unavailable=unavailable)
    if gate["audited_width_conditions"] != 450:
        raise ValueError("S1 requires the complete independently audited450-condition grid")
    plt.rcParams.update({"font.size": 8., "axes.titlesize": 9., "axes.labelsize": 8.,
                         "xtick.labelsize": 7., "ytick.labelsize": 7., "pdf.fonttype": 42})
    colors = plt.get_cmap("viridis")(np.linspace(.05, .9, len(WIDTHS)))
    s1_path = vectors / "figure2_S1_equal_area_widths.pdf"
    histories, s1_pages = [], 0
    with PdfPages(s1_path) as pdf:
        for model in MODELS:
            for area in AREAS:
                for pair in ((19, 23), (27, 31), (35, 39)):
                    fig, axes = plt.subplots(3, 2, figsize=(11.8, 10.0), layout="constrained")
                    fig.suptitle(f"S1  {NAMES[model]} • commanded area {area:g}× • native kinetics", fontsize=13)
                    for column, length in enumerate(pair):
                        overlap = dict(own=False, fixed=False)
                        for width, color in zip(WIDTHS, colors):
                            record = records[("width", model, length, width, area)]
                            trial, measurement, audit_receipt = _verified_history(record)
                            event = record["measurement"]["continuous_cue_on_events"]["D0_5"]
                            y, lo, hi, kind = _event_value(event)
                            x = width / 38 if model == "spring" else width
                            if y is not None:
                                axes[0, column].errorbar(x, y, yerr=[[y - lo], [hi - y]], color=color,
                                    fmt="^" if kind == "lower_bound" else "o", capsize=2, markersize=5)
                            series = measurement["series"]
                            time = np.asarray(series["time"])
                            on_stop = event["interval"][1]
                            endpoints = [e["end"]["time"] for e in trial["endpoints"]]
                            for row, (amount_key, overlap_key, label) in enumerate((
                                ("E_challenger_own", "roi_own_overlap", "Own cue-support excess amount"),
                                ("E_challenger_fixed", "roi_fixed_overlap", "Fixed widest-zone excess amount")), start=1):
                                values = np.asarray(series[amount_key], dtype=float)
                                index = display_indices(time, values, endpoints)
                                axes[row, column].plot(time[index], values[index], color=color, linewidth=.9,
                                                       label=f"w={x:g}")
                                axes[row, column].set_ylabel(label)
                                axes[row, column].set_xscale("symlog", linthresh=.01 if model == "spring" else .5)
                                axes[row, column].set_xlabel("Native time; ON then cue-free recovery")
                                overlap["own" if row == 1 else "fixed"] |= any(series[overlap_key])
                            histories.append(dict(case=record["case"], source=record["source"],
                                audited_metadata=audit_receipt, native_trace_sha256=trial["full_trace_sha256"],
                                readout_path=str((Path(record["source"]["run_root"]) / trial["trial_directory"] / "readouts.json")),
                                readout_sha256=trial["readouts_current_sha256"], on_stop=on_stop,
                                full_native_recovery_interval=trial.get("recovery_interval"),
                                final_d=trial["final_d"], field_tail_pass=trial["field_tail_pass"],
                                own_roi_overlap_any=any(series["roi_own_overlap"]),
                                fixed_roi_overlap_any=any(series["roi_fixed_overlap"])))
                        axes[0, column].set(title=f"L{length}: continuous-cue replacement",
                            xlabel="Native physical cue support", ylabel="Time to eligible D≥0.5", xscale="log", yscale="log")
                        for row in (1, 2):
                            axes[row, column].axvline(on_stop, color=".3", linestyle=":", linewidth=.8)
                            axes[row, column].set_title(f"L{length}; " +
                                ("ROI overlaps initial incumbent core" if overlap["own" if row == 1 else "fixed"] else "No initial incumbent-core overlap"))
                        if column == 1:
                            axes[1, column].legend(fontsize=6.5, frameon=False)
                    fig.supxlabel(textwrap.fill("Triangles are finite lower bounds. Regional excess uses the instantaneous spatial minimum; ROI conventions can overlap the initial incumbent and do not define the primary disjoint-half D. Every field-tail flag remains in the source table.", 165), fontsize=8)
                    pdf.savefig(fig, facecolor="white")
                    plt.close(fig)
                    s1_pages += 1
    write(data / "figure2_S1_complete_source_index.json", histories)
    s2_path = vectors / "figure2_S2_all_train_fate_maps.pdf"
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#dedede")
    s2_pages = 0
    with PdfPages(s2_path) as pdf:
        for pulses in PULSES:
            for factor in ON_FACTORS:
                fig, axes = plt.subplots(1, 5, figsize=(13.5, 4.0), layout="constrained")
                fig.suptitle(f"S2  {pulses} pulses; total ON = {factor:g}× postwithdrawal calibration", fontsize=13)
                for ax, model in zip(axes, MODELS):
                    values = np.full((len(GAPS), len(LENGTHS)), np.nan)
                    for yi, gap in enumerate(GAPS):
                        for xi, length in enumerate(LENGTHS):
                            key = ("train", model, length, pulses, gap, factor)
                            if key in records:
                                row = records[key]["measurement"]
                                if row["readout_eligible"] and row["final_d"] is not None:
                                    values[yi, xi] = row["final_d"]
                                else:
                                    ax.text(xi, yi, "×", ha="center", va="center")
                                if row["field_tail_pass"] is False:
                                    ax.plot(xi, yi, "s", markerfacecolor="none", color="black", markersize=10, markeredgewidth=.7)
                            elif key in unavailable:
                                ax.text(xi, yi, "†", ha="center", va="center")
                    image = ax.imshow(values, origin="lower", cmap=cmap, norm=Normalize(-1, 1), aspect="auto")
                    ax.set(title=NAMES[model], xlabel="Length label L")
                    ax.set_xticks(range(6), LENGTHS)
                    ax.set_yticks(range(6), [f"{g:g}" for g in GAPS])
                axes[0].set_ylabel("Gap / pulse duration")
                fig.colorbar(image, ax=axes.tolist(), label="Final eligible D after native recovery", fraction=.015)
                fig.supxlabel(textwrap.fill("□: failed native-field tail, independent of fate; ×: known low activity; †: no sampled calibration at the finite cap; plain grey: missing independent native audit. Native numerical schemes are retained; deBelly uses the 0.5-unit update and later-gap timing is sensitive to its qualified quarter-step companion (S5).", 190), fontsize=8)
                pdf.savefig(fig, facecolor="white")
                plt.close(fig)
                s2_pages += 1
    endpoint_rows = []
    for key in sorted(k for k in records if k[0] == "train"):
        record = records[key]
        trial, _, receipt = audited_worker_metadata(record)
        endpoint_rows.append(dict(case=record["case"], source=record["source"], audited_metadata=receipt,
                                  schedule=trial["schedule"], endpoints=trial["endpoints"],
                                  final_d=trial["final_d"], final_replacement=trial["final_replacement"],
                                  readout_eligible=trial["readout_eligible"], field_tail_pass=trial["field_tail_pass"],
                                  native_phase_join_audit_pass=trial.get("native_phase_join_audit", {}).get("passed")))
    write(data / "figure2_S2_all_native_pulse_gap_endpoints.json", endpoint_rows)
    write(data / "figure2_S2_calibration_unavailable_conditions.json", list(unavailable.values()))
    artifacts = [s1_path, s2_path, data / "figure2_S1_complete_source_index.json",
                 data / "figure2_S2_all_native_pulse_gap_endpoints.json",
                 data / "figure2_S2_calibration_unavailable_conditions.json"]
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), producer_state="stopped",
                  source_sha256=sha(Path(__file__)), native_integrations_run=0,
                  S1=dict(pages=s1_pages, audited_width_histories=450, both_roi_conventions_retained=True),
                  S2=dict(pages=s2_pages, native_train_histories=len(endpoint_rows),
                          calibration_unavailable_conditions=len(unavailable),
                          planned_conditions=1620, complete_native_or_calibration_coverage=not gate["missing_train_conditions"]),
                  source_audits=receipts, missing_train_conditions=gate["missing_train_conditions"],
                  display_reduction="At most500 index bins per curve, each min/max plus every phase endpoint; full exact histories retained by hash-linked source paths.",
                  artifact_sha256={str(p.relative_to(output)): sha(p) for p in artifacts})
    write(metadata / "S1_S2_display_audit.json", report)
    return report


def build_supporting_assays(index_path, output):
    """S4 placement and share/dwell maps, retaining all finite nonreaches."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import Normalize
    index_path, output = Path(index_path).resolve(), Path(output).resolve()
    index = read(index_path)
    if (not index.get("passed") or not index.get("declared_subset_complete") or
            index.get("audited_workers") != 270):
        raise ValueError("S4 requires the complete combined independent native audit")
    root = index_path.parent.parent
    if sha(root / "manifest.json") != index["manifest_sha256"]:
        raise ValueError("S4 frozen manifest identity changed")
    summaries = []
    for receipt in index["records"]:
        path = root / receipt["audited_worker"]
        if sha(path) != receipt["audited_worker_sha256"]:
            raise ValueError("S4 selected independent worker metadata identity changed")
        audit = read(path)
        case = audit["case"]
        trial = next((t for t in audit["trials"] if t["label"] == "trial"), None)
        geometry = None
        if case["assay"] == "position" and trial:
            edge = max(interval[1] for interval in trial["basal_descriptor"]["incumbent_intervals"])
            geometry = dict(incumbent_initial_halfheight_edge=edge,
                            cue_near_edge=case["cue_center"] - case["cue_width"] / 2,
                            edge_gap=case["cue_center"] - case["cue_width"] / 2 - edge,
                            definition="Native near-support edge minus initial finite-volume basal half-height core edge; negative means overlap.")
        summaries.append(dict(case=case, executed=receipt["withdrawal_executed"],
                              status=receipt["support_status"], geometry=geometry,
                              continuous_cue_on_events=trial["continuous_cue_on_events"] if trial else None,
                              final_d=trial["final_d"] if trial else None,
                              readout_eligible=trial["readout_eligible"] if trial else None,
                              field_tail_pass=trial["field_tail_pass"] if trial else None,
                              source=receipt))
    folders = {name: output / "supporting" / name for name in ("vectors", "data", "metadata")}
    for path in folders.values():
        path.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 8, "axes.titlesize": 10, "pdf.fonttype": 42})
    path = folders["vectors"] / "figure2_S4_position_share_dwell.pdf"
    colors = plt.get_cmap("viridis")(np.linspace(.05, .9, len(LENGTHS)))
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#dedede")
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(12.5, 8.), layout="constrained")
        fig.suptitle("S4  Centered equal-area cue placement at fixed native length", fontsize=14)
        for j, model in enumerate(MODELS[:3]):
            for length, color in zip(LENGTHS, colors):
                rows = sorted((r for r in summaries if r["case"]["model"] == model and
                               r["case"]["requested_reference_length"] == length and r["case"]["assay"] == "position"),
                              key=lambda r: r["geometry"]["edge_gap"])
                points = []
                for row in rows:
                    x = row["geometry"]["edge_gap"]
                    y, lo, hi, kind = _event_value(row["continuous_cue_on_events"]["D0_5"])
                    if y is not None:
                        axes[0, j].errorbar(x, y, yerr=[[y - lo], [hi - y]], color=color,
                            fmt="^" if kind == "lower_bound" else "o", capsize=2, markersize=4)
                        if kind == "observed":
                            points.append((x, y))
                    if row["readout_eligible"] and row["final_d"] is not None:
                        axes[1, j].plot(x, row["final_d"], "o", color=color, markersize=4)
                        if row["field_tail_pass"] is False:
                            axes[1, j].plot(x, row["final_d"], "s", color="black", markerfacecolor="none", markersize=9)
                if points:
                    axes[0, j].plot(*zip(*points), color=color, label=f"L{length}")
                if rows:
                    axes[1, j].plot([r["geometry"]["edge_gap"] for r in rows],
                                   [r["final_d"] if r["final_d"] is not None else np.nan for r in rows], color=color)
            for row in (0, 1):
                axes[row, j].set_xlabel("Native gap from initial half-height core")
                axes[row, j].axvspan(min(-1., axes[row, j].get_xlim()[0]), 0, color=".92", zorder=-1)
            axes[0, j].set(title=NAMES[model], ylabel="Cue-ON time to D≥0.5", yscale="log")
            axes[1, j].set(ylabel="Final D after recovery", ylim=(-1.05, 1.05))
        axes[0, 2].legend(fontsize=7, frameon=False, ncol=2)
        fig.supxlabel(textwrap.fill("Each length uses three centered full-support cosine cues with the same commanded spatial area. This is a supporting native placement test, distinct from the Figure 3 ring experiment. Triangles are lower bounds; open squares mark failed finite field tails. deBelly distances are native nominal solver coordinates.", 180), fontsize=8)
        pdf.savefig(fig, facecolor="white")
        plt.close(fig)
        for model in MODELS[:3]:
            fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), layout="constrained")
            fig.suptitle("S4  " + NAMES[model] + ": share-triggered withdrawal and added dwell", fontsize=14)
            for ax, target in zip(axes, (.5, .75, .9)):
                values = np.full((3, 6), np.nan)
                for yi, dwell in enumerate((0., .25, 1.)):
                    for xi, length in enumerate(LENGTHS):
                        row = next(r for r in summaries if r["case"]["model"] == model and
                            r["case"]["assay"] == "share_withdrawal" and
                            r["case"]["requested_reference_length"] == length and
                            r["case"]["target_share"] == target and r["case"]["dwell_factor_of_calibration"] == dwell)
                        if not row["executed"]:
                            ax.text(xi, yi, "†", ha="center", va="center")
                        elif row["readout_eligible"] and row["final_d"] is not None:
                            values[yi, xi] = row["final_d"]
                        else:
                            ax.text(xi, yi, "×", ha="center", va="center")
                        if row["field_tail_pass"] is False:
                            ax.plot(xi, yi, "s", color="black", markerfacecolor="none", markersize=12)
                image = ax.imshow(values, origin="lower", cmap=cmap, norm=Normalize(-1, 1), aspect="auto")
                ax.set(title=f"First sampled challenger share ≥ {target:g}", xlabel="Length label L")
                ax.set_xticks(range(6), LENGTHS)
                ax.set_yticks(range(3), ["0", "0.25", "1"])
            axes[0].set_ylabel("Additional dwell / commitment calibration")
            fig.colorbar(image, ax=axes.tolist(), fraction=.02, label="Final eligible D after recovery")
            fig.supxlabel(textwrap.fill("†: share target not reached in the finite shared native prefix; no withdrawal was executed. ×: observed low activity. □: failed finite native-field tail, independent of fate. Exact first-hit native states are reused without interpolation; threshold and added dwell remain separate. deBelly retains its native 0.5-unit clock.", 180), fontsize=8)
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)
    table = folders["data"] / "figure2_S4_all270_results.json"
    write(table, summaries)
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True, native_integrations_run=0,
                  pages=4, source_index=str(index_path), source_index_sha256=sha(index_path),
                  audited_workers=270, native_trials=sum(r["executed"] for r in summaries),
                  finite_prefix_nonreaches=sum(not r["executed"] for r in summaries),
                  artifact_sha256={str(p.relative_to(output)): sha(p) for p in (path, table)},
                  source_sha256=sha(Path(__file__)))
    write(folders["metadata"] / "S4_display_audit.json", report)
    return report


def build_debelly_numerics(summary_path, output):
    """S5 material temporal qualifications; successful fate cannot hide timing."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    summary_path, output = Path(summary_path).resolve(), Path(output).resolve()
    summary = read(summary_path)
    if summary.get("independent_native_identity_audits_passed") is not True or summary.get("actual_protocol_comparisons") != 6:
        raise ValueError("S5 temporal display requires all six independently audited qualified companions")
    for name, expected in summary["outputs"].items():
        if sha(summary_path.parent / name) != expected:
            raise ValueError("S5 temporal scientific table identity changed")
    with (summary_path.parent / "debelly_native_endpoint_profiles.csv").open() as stream:
        profiles = list(csv.DictReader(stream))
    folders = {name: output / "supporting" / name for name in ("vectors", "data", "metadata")}
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)
    copies = []
    for name in summary["outputs"]:
        target = folders["data"] / ("figure2_S5_debelly_" + name)
        shutil.copyfile(summary_path.parent / name, target)
        copies.append(target)
    copy = folders["data"] / "figure2_S5_debelly_temporal_comparisons.json"
    shutil.copyfile(summary_path, copy)
    copies.append(copy)
    pdf_path = folders["vectors"] / "figure2_S5_debelly_temporal_qualification.pdf"
    plt.rcParams.update({"font.size": 8.5, "axes.labelsize": 8, "axes.titlesize": 10, "pdf.fonttype": 42})
    records = summary["records"]
    if not all(record["identical_realized_phase_schedule"] for record in records):
        raise ValueError("Temporal display requires identical realized phase schedules")
    trains = [r for r in records if r["assay"].startswith("train_")]
    colors = ("#82599C", "#C67535")
    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(12.5, 7.5), layout="constrained")
        ax.axis("off")
        ax.set_title("S5  deBelly: native identity audits pass; numerical agreement is incomplete", fontsize=14, loc="left")
        cells = []
        for record in records:
            maximum_D = record["phase_comparison"]["maximum_absolute_eligible_D_difference"]
            cells.append([record["requested_reference_length"], record["assay"],
                f"{100 * record['source_basal_field_shift']['maximum_scaled_difference']:.3g}%",
                "unavailable" if maximum_D is None else f"{maximum_D:.3g}",
                f"{100 * record['endpoint_native_fields']['maximum_scaled_difference']:.3g}%",
                f"{record['main201_final']['final_d']:.3g} / {record['companion_final_D']:.3g}",
                f"{record['main201_final']['field_tail_pass']} / {record['companion_field_tail_pass']}",
                "\n".join(record["issues"]) if record["issues"] else "none"])
        table = ax.table(cellText=cells, colLabels=["L", "Assay", "Basal field shift", "Max endpoint |ΔD|",
            "Max endpoint field error", "Final D:0.5 /0.25", "Tail:0.5 /0.25", "Nonexclusive qualifications"],
            bbox=[0, .23, 1, .65], cellLoc="center", colWidths=[.04, .13, .11, .11, .13, .12, .11, .25])
        table.auto_set_font_size(False)
        table.set_fontsize(7.0)
        ax.text(0, .13, textwrap.fill("The 0.25-unit companions use separately qualified basal states and a controller relaxation rate adjusted to preserve its physical rate. The native 0.5-unit campaign remains unchanged. Final fate agreement does not establish agreement in phase timing, eligibility, coupled native fields or finite-tail status.", 170), transform=ax.transAxes, fontsize=9)
        ax.text(0, .025, textwrap.fill("At L39 the first two gaps still lose challenger bias, but third-gap D changes from −0.58556 to +0.18947. Final D changes from 0.91646 to 0.92085. Tension and velocity endpoint differences reach 8.341% and 6.802% of native scales. The first L19 gap also loses bias in both schemes while its eligibility flag changes. These are explicit numerical limits, not omitted failed simulations.", 170), transform=ax.transAxes, fontsize=9)
        pdf.savefig(fig, facecolor="white")
        plt.close(fig)
        fig, axes = plt.subplots(2, 1, figsize=(12.5, 8.), layout="constrained")
        fig.suptitle("S5  Every native pulse/gap endpoint under the identical realized cue schedule", fontsize=14)
        for ax, record in zip(axes, trains):
            points = record["phase_comparison"]["endpoints"]
            for key, color, label in (("main201", colors[0], "Native0.5 update"),
                                      ("companion", colors[1], "Qualified0.25 update")):
                eligible = [np.nan if p[key]["D_half"] is None else p[key]["D_half"] for p in points]
                ax.plot(range(len(points)), eligible, "o-", color=color, linewidth=1, label=label)
                for i, p in enumerate(points):
                    if p[key]["D_half"] is None and p[key]["D_half_raw"] is not None:
                        ax.plot(i, p[key]["D_half_raw"], "^", color=color, markerfacecolor="none", markersize=8)
            labels = []
            on, gap = 0, 0
            for point in points:
                if point["name"].endswith("on"):
                    on += 1
                    labels.append(f"ON{on}")
                elif point["name"].endswith("gap"):
                    gap += 1
                    labels.append(f"Gap{gap}")
                else:
                    labels.append("Recovery")
            ax.set(title=f"L{record['requested_reference_length']}", ylabel="Eligible activity dominance D", ylim=(-1.05, 1.05))
            ax.set_xticks(range(len(points)), labels)
            ax.axhline(0, color=".4", linewidth=.6)
            ax.axhline(.5, color=".5", linestyle=":", linewidth=.8)
            ax.legend(frameon=False)
        fig.supxlabel(textwrap.fill("Points are exact phase-end samples; connecting segments do not represent continuous trajectories or equal time spacing. Open triangles show raw D when the amount/contrast eligibility gate fails; raw values are not substituted for eligible D. Full native times, eligibility, contrast, amounts and gate results are retained in the companion CSV.", 180), fontsize=8)
        pdf.savefig(fig, facecolor="white")
        plt.close(fig)
        fields = ("rac", "rho", "mca_bound", "mca_unbound", "displacement", "velocity", "tension",
                  "raw_left_displacement", "raw_right_displacement", "left_boundary_displacement", "right_boundary_displacement")
        for record in trains:
            length = record["requested_reference_length"]
            fig, axes = plt.subplots(3, 4, figsize=(12.5, 9.5), layout="constrained")
            fig.suptitle(f"S5  deBelly L{length}: all 11 native endpoint state components after recovery", fontsize=14)
            for ax, field in zip(axes.flat, fields):
                for profile, color, label in (("main201", colors[0], "Native0.5"),
                                               (record["variant"], colors[1], "Qualified0.25")):
                    selected = sorted((p for p in profiles if p["profile"] == profile and
                        p["assay"] == record["assay"] and float(p["requested_reference_length"]) == length), key=lambda p: int(p["cell"]))
                    ax.plot([float(p["nominal_x"]) for p in selected], [float(p[field]) for p in selected],
                            color=color, linewidth=1, label=label)
                ax.set(title=field.replace("_", " "), xlabel="Nominal native position", ylabel="Native field units")
            axes[0, 0].legend(fontsize=7, frameon=False)
            axes[-1, -1].axis("off")
            axes[-1, -1].text(0, .85, textwrap.fill("No fitted spatial shift or amplitude normalization. Final fate, native-field error and finite-tail gate are separate readouts. All six width/train comparisons and every state profile remain in the copied source tables.", 35), fontsize=9, va="top")
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
                  source_summary=str(summary_path), source_summary_sha256=sha(summary_path),
                  independent_native_identity_audits_passed=True,
                  numerical_agreement_for_all_requested_readouts=summary["numerical_agreement_for_all_requested_readouts"],
                  native_integrations_run=0, pages=4,
                  artifact_sha256={str(p.relative_to(output)): sha(p) for p in [pdf_path, *copies]},
                  source_sha256=sha(Path(__file__)))
    write(folders["metadata"] / "S5_debelly_temporal_display_audit.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, help="Renderer review_inputs.json with audited directories")
    parser.add_argument("--support-index", type=Path)
    parser.add_argument("--debelly-numerics-summary", type=Path)
    parser.add_argument("--noise-summary", type=Path)
    parser.add_argument("--numerics-summary", type=Path)
    parser.add_argument("--combine", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.combine:
        print(json.dumps(combine_supplements(args.output), indent=2))
        raise SystemExit(0)
    if args.numerics_summary:
        print(json.dumps(build_numerics(args.numerics_summary, args.output), indent=2))
        raise SystemExit(0)
    if args.noise_summary:
        print(json.dumps(build_noise(args.noise_summary, args.output), indent=2))
        raise SystemExit(0)
    if args.debelly_numerics_summary:
        print(json.dumps(build_debelly_numerics(args.debelly_numerics_summary, args.output), indent=2))
        raise SystemExit(0)
    if args.support_index:
        print(json.dumps(build_supporting_assays(args.support_index, args.output), indent=2))
        raise SystemExit(0)
    if args.inputs is None:
        parser.error("--inputs is required for S1/S2")
    config = read(args.inputs)
    records, receipts = load_audited_records(config["analysis_directories"])
    unavailable, _ = load_calibration_unavailable(config.get("unavailable_trains"))
    print(json.dumps(build_supplements(records, receipts, unavailable, args.output), indent=2))
