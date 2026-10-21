"""Render the independently audited Figure 4 E0 interim supplement.

No native integrations or scientific estimators live here.  Each PDF page has
an unabridged CSV of its plotted source rows.  E1--E4-dependent decisions stay
provisional, and finite-window or cadence limits never become fitted values.
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

# Final H3 records retain every preset pair in a nested CSV cell. Capacity is
# an engineering parser setting, not truncation or a scientific threshold.
csv.field_size_limit(64*1024*1024)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from matplotlib.lines import Line2D


MANUSCRIPT = Path(__file__).resolve().parents[2]
DEFAULT_RUN = MANUSCRIPT / "data/runs/signature_figure4_e0_corrected_basals_20261020"
DEFAULT_OUT = MANUSCRIPT / "figures/length_mechanisms_20261008/dev/signature_figure4_20261020"
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring_ring")
LABELS = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "Holmes 3", "Spring")))
COLORS = ("#2878b5", "#34974a", "#e58b29")
FILES = ("condition_metrics.csv", "heading_curves.csv", "site_curves.csv",
         "size_ratios.csv", "noise_tolerance.csv", "decisions.csv", "case_metrics.csv")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def read_csv(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Missing or duplicate CSV columns: {path}")
        return list(reader)


def display_rows(rows):
    """Add display aliases without discarding any original scientific column."""
    aliases = dict(tau="tau_p", tau_status="tau_p_status", tau_lower="tau_p_lower",
                   tau_upper="tau_p_upper", tau_coverage="tau_p_coverage",
                   tau_ci_lower="tau_p_ci_low", tau_ci_upper="tau_p_ci_high",
                   tau_ci_upper_unbounded="tau_p_ci_upper_unbounded",
                   Dtheta_ci_lower="Dtheta_ci_low", Dtheta_ci_upper="Dtheta_ci_high",
                   Dx_ci_lower="Dx_ci_low", Dx_ci_upper="Dx_ci_high",
                   C_ci_lower="C_ci_low", C_ci_upper="C_ci_high",
                   angular_MSD="angular_msd", angular_MSD_ci_lower="angular_msd_ci_low",
                   angular_MSD_ci_upper="angular_msd_ci_high", MSD="msd",
                   MSD_ci_lower="msd_ci_low", MSD_ci_upper="msd_ci_high",
                   single_front_fraction="single_fraction", multiple_fronts_fraction="multiple_fraction",
                   actual_geometry_ratio="geometric_ratio")
    result=[]
    for original in rows:
        row=dict(original)
        for source,target in aliases.items():
            if source in row and target not in row:row[target]=row[source]
        geometry=number(row.get("actual_geometry_ratio"))
        if geometry>0:
            for side in ("lower","upper"):
                value=number(row.get("tau_ratio_"+side))
                if not math.isnan(value):row.setdefault("excess_persistence_"+side,str(value/geometry))
        for prefix in ("tau_ratio","excess_persistence","fhalf"):
            for side in ("lower","upper"):
                source=prefix+"_ci_"+side
                if source in row:row.setdefault(prefix+"_ci_"+("low" if side=="lower" else "high"),row[source])
            if prefix in row and prefix+"_status" not in row:
                point=number(row[prefix]);low=number(row.get(prefix+"_lower"));high=number(row.get(prefix+"_upper"))
                row[prefix+"_status"]=("resolved" if math.isfinite(point) else
                    "bounded_interval" if math.isfinite(low) and low>0 else "unavailable")
        result.append(row)
    return result


def number(value, default=math.nan):
    if value is None or str(value).strip().lower() in ("", "none", "null", "nan"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def get(row, *keys, default=None):
    for key in keys:
        if key in row and row[key] not in (None, ""):
            return row[key]
    return default


def num(row, *keys, default=math.nan):
    return number(get(row, *keys), default)


def boolean(value):
    if value is True or str(value).lower() in ("true", "1"):
        return True
    if value is False or str(value).lower() in ("false", "0"):
        return False
    return None


def outcome(value):
    labels={"pass_provisional":"criterion provisionally met", "pass_model_criterion_provisional":"model criterion provisionally met",
            "inconclusive":"inconclusive", "fail":"criterion not met"}
    return labels.get(str(value),str(value).replace("_"," "))


def yesno(value):
    value=boolean(value)
    return "yes" if value is True else "no" if value is False else "unavailable"


def interval_text(value):
    try:
        values=json.loads(str(value))
    except (ValueError,TypeError):return str(value)
    if not isinstance(values,list) or len(values)!=2:return str(value)
    return "["+", ".join("∞" if x is None else f"{float(x):g}" for x in values)+"]"


def reversal_text(value):
    try:
        values=json.loads(str(value))
    except (ValueError,TypeError):return str(value)
    return ", ".join(f"{float(x):g}" for x in values) if values else "none observed"


def model(row):
    value = row.get("model", "")
    return "spring_ring" if value == "spring" else value


def main_rows(rows, view="main"):
    return [r for r in rows if get(r, "view", default="main") == view
            and get(r, "kind", default="noise") == "noise"]


def subset(rows, *, which_model=None, strength=None, size=None, view=None):
    return [r for r in rows if
            (which_model is None or model(r) == which_model) and
            (strength is None or math.isclose(num(r, "f"), strength, abs_tol=1e-14)) and
            (size is None or math.isclose(num(r, "C_label"), size, abs_tol=1e-10)) and
            (view is None or get(r, "view", default="main") == view)]


def condition_key(row):
    return model(row), num(row, "C_label"), num(row, "f"), get(row, "view", default="main")


def size_color(which_model, size, forced=False):
    if forced:
        return "#888888"
    sizes = (1., 2., 3.) if which_model == "spring_ring" else (20., 40., 60.)
    position = max(0., min(2., 2 * (size - sizes[0]) / (sizes[-1] - sizes[0])))
    left = min(1, int(position)); weight = position - left
    rgb = [[int(color[start:start + 2], 16) for start in (1, 3, 5)] for color in COLORS]
    return "#" + "".join(f"{round((1-weight)*a+weight*b):02x}" for a, b in zip(rgb[left], rgb[left+1]))


def forced(row):
    if "display_forced" in row:
        return boolean(row["display_forced"])
    return boolean(get(row, "forced", "forced_state", "matched_mean_forced"))


def load_inputs(run_root):
    """Reject incomplete or mutated analysis before creating any display file."""
    run_root = Path(run_root).resolve()
    analysis = run_root / "results/analysis"
    status_path = analysis / "status.json"
    status = json.loads(status_path.read_text())
    complete = status.get("complete", status.get("declared_condition_groups_complete", False))
    if status.get("passed") is not True or complete is not True:
        raise ValueError("E0 scientific analysis is not declared complete and passed")
    if status.get("errors") or status.get("missing") or status.get("missing_conditions"):
        raise ValueError("E0 status retains scientific errors or missing conditions")
    native_path = run_root / "results/native_identity_verification.json"
    native = json.loads(native_path.read_text())
    inventory_path = run_root / "inventory.json"
    if native.get("passed") is not True or native.get("inventory_sha256") != sha(inventory_path):
        raise ValueError("E0 native identity inventory gate failed")
    correction_path=run_root/"results/basal_layout_correction.json"
    correction=json.loads(correction_path.read_text())
    if (correction.get("passed") is not True or correction.get("corrected_inventory_sha256")!=sha(inventory_path)
            or not correction.get("all_raw_contrasts_identical") or len(correction.get("rows",[]))!=12
            or not all(row.get("native_initial_primary_exact") for row in correction["rows"])):
        raise ValueError("E0 requires the corrected native Spring basal adapter and all12 exact profiles")
    sources = dict(inventory=dict(path=str(inventory_path), sha256=sha(inventory_path)),
                   native_identity=dict(path=str(native_path), sha256=sha(native_path)),
                   status=dict(path=str(status_path), sha256=sha(status_path)),
                   corrected_basal_layout=dict(path=str(correction_path),sha256=sha(correction_path)))
    invariant_path=run_root/"results/scalar_invariance_after_basal_correction.json"
    invariant=json.loads(invariant_path.read_text())
    if (not invariant.get("passed") or not invariant.get("complete") or invariant.get("condition_groups")!=110
            or not invariant.get("all_metrics_bit_identical") or not invariant.get("all_npz_file_hashes_identical")):
        raise ValueError("Corrected-basal E0 reanalysis/invariance receipt is incomplete")
    sources["corrected_scalar_reanalysis"]=dict(path=str(invariant_path),sha256=sha(invariant_path))
    independent_path=analysis/"independent_estimand_audit.json"
    independent=json.loads(independent_path.read_text())
    if (independent.get("passed") is not True or independent.get("complete") is not True
            or independent.get("condition_groups")!=110 or not independent.get("all_seven_views")):
        raise ValueError("E0 independent estimand audit failed or remains incomplete")
    manifest_sha=sha(run_root/"analysis_manifest.json")
    if independent.get("analysis_manifest_sha256")!=manifest_sha or status.get("analysis_manifest_sha256")!=manifest_sha:
        raise ValueError("E0 independent estimand/collector manifests disagree")
    if {int(r["task"]) for r in independent.get("receipts",[])}!=set(range(110)):
        raise ValueError("E0 independent estimand receipts do not cover all conditions")
    for reference in independent["receipts"]:
        receipt_path=run_root/f"estimand_audit_workers/{int(reference['task']):04d}.json"
        if sha(receipt_path)!=reference["receipt_sha256"]:
            raise ValueError("E0 independent estimand worker receipt changed")
        receipt=json.loads(receipt_path.read_text())
        producer=run_root/f"analysis_workers/{int(reference['task']):04d}/result.json"
        if receipt.get("passed") is not True or receipt.get("producer_result_sha256")!=sha(producer):
            raise ValueError("E0 producer changed after independent estimand verification")
    sources["independent_estimands"]=dict(path=str(independent_path),sha256=sha(independent_path))
    tables = {}
    for name in FILES:
        path = analysis / name
        claimed = status.get("outputs", {}).get(name)
        if isinstance(claimed, dict):
            claimed = claimed.get("sha256")
        if not claimed or claimed != sha(path):
            raise ValueError(f"E0 CSV lacks its verified digest: {name}")
        tables[name] = display_rows(read_csv(path))
        sources[name] = dict(path=str(path), sha256=claimed, rows=len(tables[name]))
    estimands_path = run_root / "inputs/estimands.json"
    if status.get("estimands_sha256")!=sha(estimands_path) or status.get("inventory_sha256")!=sha(inventory_path):
        raise ValueError("E0 scientific estimands or inventory changed after collection")
    if status.get("native_identity_sha256")!=sha(native_path):
        raise ValueError("E0 collected native-identity receipt changed")
    sources["estimands"] = dict(path=str(estimands_path), sha256=sha(estimands_path))
    conditions = main_rows(tables["condition_metrics.csv"])
    keys = [condition_key(r) for r in conditions]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate main noisy E0 condition")
    inventory = read_csv(run_root / "inventory.csv")
    expected = {condition_key(dict(r, view="main")) for r in inventory if r["kind"] == "noise"}
    if set(keys) != expected:
        raise ValueError("E0 main condition table does not cover the frozen noisy inventory")
    if len(expected) != 110:
        raise ValueError("The declared E0 design requires all 110 noisy conditions")
    if set(model(r) for r in conditions) != set(MODELS):
        raise ValueError("E0 does not cover all four native ring models")
    views={"main","analysis_first_half","analysis_second_half","late_quarter","cadence_x2","cadence_x4","cadence_x8"}
    all_keys=[condition_key(r) for r in tables["condition_metrics.csv"]]
    expected_keys={(m,c,f,view) for m,c,f,_ in expected for view in views}
    if len(all_keys)!=len(expected_keys) or set(all_keys)!=expected_keys:
        raise ValueError("E0 sensitivity windows/cadences do not cover all 110 conditions")
    if len(tables["case_metrics.csv"])!=962*7:
        raise ValueError("E0 case table omits or duplicates native realizations/controls across views")
    memo_path=run_root/"results/interim_handoff/memo.json"
    memo=json.loads(memo_path.read_text())
    if (not memo.get("passed") or memo.get("E0_status_sha256")!=sha(status_path) or
            memo.get("independent_estimand_audit_sha256")!=sha(independent_path)):
        raise ValueError("E0 interim decision memo does not match the audited source")
    for name,digest in memo["outputs"].items():
        path=memo_path.parent/name
        if sha(path)!=digest:raise ValueError("E0 interim decision/trigger table changed")
        tables[name]=read_csv(path)
        sources[name]=dict(path=str(path),sha256=digest,rows=len(tables[name]))
    tables["memo"]=memo
    sources["memo"]=dict(path=str(memo_path),sha256=sha(memo_path))
    return tables, inventory, sources, status


def interval_mark(ax, x, row, prefix, *, color, marker="o", label=None):
    """Plot point/bound/unresolved intervals without converting a bound to a fit."""
    state = str(get(row, prefix + "_status", "status", default="unavailable"))
    point = num(row, prefix, prefix + "_estimate", prefix + "_point")
    low = num(row, prefix + "_lower", prefix + "_lower_bound")
    high = num(row, prefix + "_upper", prefix + "_upper_bound")
    ci_low = num(row, prefix + "_ci_low", prefix + "_ci_lower")
    ci_high = num(row, prefix + "_ci_high", prefix + "_ci_upper")
    unavailable = any(word in state for word in ("unavailable", "inconclusive", "unresolved"))
    resolved = "resolved" in state and not unavailable and "censored" not in state
    if resolved and math.isfinite(point):
        ax.plot(x, point, marker=marker, color=color, ms=4, label=label,zorder=4)
        if math.isfinite(ci_low) and math.isfinite(ci_high):
            ax.vlines(x, ci_low, ci_high, color=color, lw=1)
        elif boolean(row.get(prefix+"_ci_upper_unbounded")) is True or ci_high==math.inf:
            # The empirical point can be identified while its whole-seed
            # uncertainty has no finite upper limit. Never substitute the
            # observation cap or drop this uncertainty silently.
            if math.isfinite(ci_low) and ci_low>0:
                ax.vlines(x,ci_low,point,color=color,lw=1,ls=":")
            text="CI[0,∞)" if ci_low==0 else "CI∞"
            ax.annotate(text,xy=(x,point),xytext=(3,4),textcoords="offset points",
                        color=color,fontsize=6,annotation_clip=False)
            if isinstance(row,dict):
                row["display_"+prefix+"_CI_marker"]=text
                row["display_"+prefix+"_finite_point_upper_unbounded_CI"]=True
                row["display_"+prefix+"_CI_lower_endpoint"]=get(row,prefix+"_ci_low",prefix+"_ci_lower")
        return "resolved"
    if math.isfinite(low) and low > 0 and ("censored" in state or "bound" in state or
                                         "inconclusive" in state and high >= low):
        ax.plot(x, low, marker="^", color=color, ms=5, fillstyle="none", label=label)
        if math.isfinite(high) and high > low:
            ax.vlines(x, low, high, color=color, ls=":", lw=1)
        return "bounded"
    if unavailable and math.isfinite(high) and high > 0:
        ax.plot(x, high, marker="v", color=color, ms=5, fillstyle="none", label=label)
        return "unresolved"
    ax.plot(x,.025,"x",color=color,ms=3,transform=ax.get_xaxis_transform(),label=label)
    return "unavailable"


def curve(ax, rows, xkey, ykeys, *, color, label, lowkeys=(), highkeys=(), ls="-"):
    rows = sorted(rows, key=lambda r: num(r, xkey))
    x = np.array([num(r, xkey) for r in rows])
    y = np.array([num(r, *ykeys) for r in rows])
    valid = np.isfinite(x) & np.isfinite(y)
    ax.plot(x, np.where(valid, y, np.nan), color=color, lw=.9, label=label, ls=ls)
    if lowkeys and highkeys:
        lo = np.array([num(r, *lowkeys) for r in rows])
        hi = np.array([num(r, *highkeys) for r in rows])
        ax.fill_between(x, lo, hi, where=valid & np.isfinite(lo) & np.isfinite(hi),
                        color=color, alpha=.13, linewidth=0)


def polish(ax, *, xlabel=None, ylabel=None, logx=False, logy=False):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7, length=2)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")


class Pages:
    def __init__(self, output, sources):
        self.output = Path(output)
        self.vectors = self.output / "supporting/vectors"
        self.data = self.output / "supporting/data"
        self.metadata = self.output / "supporting/metadata"
        for path in (self.vectors, self.data, self.metadata):
            path.mkdir(parents=True, exist_ok=True)
        self.path = self.vectors / "figure4_E0_interim_supplement.pdf"
        self.pdf = PdfPages(self.path)
        self.sources = sources
        self.pages = []
        self.generated = datetime.now(timezone.utc).isoformat()

    def save(self, fig, page, title, source_rows, note):
        if not source_rows:
            raise ValueError(f"Page {page} has no source rows")
        sequence = len(self.pages) + 1
        csv_path = self.data / f"figure4_E0_{sequence:02d}_{page}.csv"
        fields = list(dict.fromkeys(key for row in source_rows for key in row))
        with csv_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(source_rows)
        display_page=page.split("_")[0]
        fig.suptitle(f"{display_page}  {title}", x=.055, ha="left", y=.985, fontsize=13)
        footer = "E0 INTERIM · A1–A4/A7 · provisional decisions · " + note
        lines = textwrap.wrap(footer, 150)
        fig.text(.055, .018, "\n".join(lines) + f"\nCSV: {csv_path.name}",
                 fontsize=6.3, va="bottom", ha="left")
        self.pdf.savefig(fig)
        if sequence == 1:
            fig.savefig(self.output / "figure4_E0_interim_preview.png", dpi=150)
        plt.close(fig)
        self.pages.append(dict(page=page, title=title, pdf_page=sequence,
                               csv=str(csv_path), csv_sha256=sha(csv_path), rows=len(source_rows)))

    def close(self):
        self.pdf.close()
        report = dict(generated_utc=self.generated, producer_state="stopped", interim=True,
                      final_main_framing_selected=False, passed=True,
                      pdf=dict(path=str(self.path), sha256=sha(self.path), pages=len(self.pages)),
                      source_files=self.sources, pages=self.pages,
                      renderer_source=dict(path=str(Path(__file__).resolve()),sha256=sha(Path(__file__))),
                      excluded_pending_sections=["S7", "S8", "S9", "S10", "S11"],
                      pending=["E1", "E2", "E3", "E4", "A5", "A6", "final H1/H2/H3 criteria"])
        write_json(self.metadata / "interim_display_audit.json", report)
        return report


def decision_pages(pages, tables):
    rows = tables["decisions.csv"]
    fig, axes = plt.subplots(4, 1, figsize=(11.7, 8.3))
    for ax, which_model in zip(axes, MODELS):
        ax.axis("off")
        selected = subset(rows, which_model=which_model, view="main")
        ax.text(0, .96, LABELS[which_model], va="top", fontsize=11, weight="bold")
        lines = []
        for row in selected:
            inclusion="Forced-inclusive" if boolean(row.get("forced_inclusive")) else "Primary non-forced"
            evidence=(f"Evaluable {row.get('evaluable_conditions','?')}; gain CIs exclude 1: {row.get('gain_CI_excludes_one_count','?')}; "
                f"reversal: {reversal_text(row.get('reversal_f','?'))}; f½ small CI {interval_text(row.get('fhalf_small_ci','?'))}, "
                f"large CI {interval_text(row.get('fhalf_large_ci','?'))}; f½ increase CIs separate: {yesno(row.get('fhalf_increase_nonoverlap'))}.")
            sentence=f"{inclusion}: H1 {outcome(row.get('H1','inconclusive'))}; H2 {outcome(row.get('H2','inconclusive'))}. {evidence}"
            lines.extend(textwrap.wrap(sentence.replace("null","∞").replace("_"," "),145))
        ax.text(0, .73, "\n".join(lines), va="top", fontsize=8, linespacing=1.4)
    fig.subplots_adjust(left=.06, right=.97, top=.92, bottom=.14, hspace=.35)
    pages.save(fig, "S0", "Provisional per-model hypothesis decisions", rows,
               "No final overall H1/H2/H3 claim: A6, fixed-mean sigma sweeps, fine profiles and extended horizons remain pending.")
    memo=tables["memo"]
    fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis("off")
    paragraphs=[("Provisional H1",memo["provisional_H1"]),("Provisional H2",memo["provisional_H2"]),
                ("Physical scope",memo["physical_scope"])]
    holmes_windows=[r for r in rows if model(r)=="holmes_model3" and boolean(r.get("forced_inclusive")) is False]
    unconfirmed=[r["view"].replace("analysis_","").replace("_"," ") for r in holmes_windows if not str(r.get("H2","")).startswith("pass")]
    if unconfirmed:
        paragraphs[1]=(paragraphs[1][0],paragraphs[1][1]+" Holmes' criterion is not robust across all declared windows: "+", ".join(unconfirmed)+" remain inconclusive.")
    source_rows=[dict(section=title,text=body) for title,body in paragraphs]
    y=.93
    for title,body in paragraphs:
        ax.text(0,y,title,fontsize=10,weight="bold",va="top")
        ax.text(0,y-.04,"\n".join(textwrap.wrap(body,145)),fontsize=9,va="top")
        y-=.14
    for alternative in memo["alternatives"]:
        ax.text(0,y,alternative["framing"],fontsize=10,weight="bold",va="top")
        text="Support: "+alternative["support"]+" Weakness: "+alternative["weakness"]+" Needed: "+alternative["next"]
        ax.text(0,y-.04,"\n".join(textwrap.wrap(text,145)),fontsize=9,va="top")
        source_rows.append(dict(section="Alternative framing",**alternative));y-=.16
    fig.subplots_adjust(left=.06,right=.97,top=.93,bottom=.13)
    pages.save(fig,"S0_memo","Decision memo and possible later framings",source_rows,
               "Alternatives describe evidence and dependencies; none is selected as the main-figure framing. Holmes' local peak-drop criterion does not establish overall H2.")


def inventory_page(pages, inventory, tables):
    fig = plt.figure(figsize=(11.7, 8.3))
    ax = fig.add_axes([.06, .44, .88, .46]); ax.axis("off")
    ax.text(.5, .98, r"$I=f I_{\rm ref}\exp(\sigma Z-\sigma^2/2),\quad \sigma=1$", fontsize=17, va="top",ha="center")
    ax.text(.99,.98,"(1)",fontsize=12,va="top",ha="right")
    description = [
        "Four native ring models; chemical C = 20, 40, 60; Spring resting C = 1, 2, 3.",
        "Original Gaussian-OU spatial forcing retained. Kernel length 1.5 (Spring 0.15); relaxation time 50 (Spring 5).",
        "Native horizons 32,000 (Spring 3,200); discard the first half. Heading/scalar cadence 5 (Spring 0.5).",
        "Saved profiles are coarser: 40 (Spring 4). E0 region lifetimes and peak-position sensitivities have that resolution.",
        "Cosine endpoints require R ≥ 0.1. MSD tracks never bridge undefined headings. CIs resample whole independent seeds (2,000 draws).",
        "Persistence is the empirical first 1/e crossing. Upward triangles denote finite-window lower bounds; downward triangles denote cadence limits.",
        "Grey indicates an abnormal matched-mean control. Missing forced-state classification stays unavailable.",
        "Spring physical increments use the measured circumference at each frame. Native times are compared only within a model.",
        "Spring's forcing stays attached to material nodes: its physical covariance width is 0.15 × C(t)/C(basal), rather than a retuned fixed optical width.",
        "Blue: smallest ring; green: intermediate ring; orange: largest ring. Grey: forced matched-mean condition, retained in every display.",
    ]
    ax.text(0, .78, "\n\n".join("\n".join(textwrap.wrap(line, 135)) for line in description),
            va="top", fontsize=9)
    ax2 = fig.add_axes([.06, .14, .88, .25]); ax2.axis("off")
    summaries = []
    for which_model in MODELS:
        selected = subset(inventory, which_model=which_model)
        strengths = sorted({num(r, "f") for r in selected if r["kind"] == "noise"})
        counts = [sum(r["kind"] == kind for r in selected) for kind in ("noise", "mean", "zero")]
        seeds = sorted({sum(r["kind"] == "noise" and num(r, "f") == f and num(r, "C_label") == c for r in selected)
                        for f in strengths for c in sorted({num(r, "C_label") for r in selected})
                        if any(r["kind"] == "noise" and num(r, "f") == f and num(r, "C_label") == c for r in selected)})
        summaries.append([LABELS[which_model], "/".join(map(str, counts)), ", ".join(map(str, seeds)),
                          ", ".join(f"{f:g}" for f in strengths)])
    table = ax2.table(cellText=summaries, colLabels=["Model", "Noisy / mean / zero", "Seeds / condition", "Observed input fractions f"],
                      colWidths=[.16, .16, .14, .54], loc="upper left", cellLoc="left")
    table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1, 2)
    pages.save(fig, "S1", "Design, forcing and complete existing-case inventory", inventory,
               "All native source identities rehashed independently; the two-seed outer strengths retain only two independent realizations.")


def regime_pages(pages, tables):
    rows = main_rows(tables["condition_metrics.csv"])
    metric_sets = [(("tau_p", r"Empirical $\tau_p$", True), ("tau_coverage", "Heading-pair coverage", False),
                    ("forced", "Declared forced-condition flag", False)),
                   (("single_fraction", "Single-region fraction", False), ("multiple_fraction", "Multiple-region fraction", False),
                    ("low_contrast_fraction", "Low-contrast fraction", False))]
    for index, metrics in enumerate(metric_sets, 1):
        fig, axes = plt.subplots(4, 3, figsize=(11.7, 8.3), squeeze=False)
        page_rows = []
        for i, which_model in enumerate(MODELS):
            selected = subset(rows, which_model=which_model)
            strengths = sorted({num(r, "f") for r in selected}); sizes = sorted({num(r, "C_label") for r in selected})
            for j, (key, title, logarithmic) in enumerate(metrics):
                ax = axes[i, j]; data = np.full((len(sizes), len(strengths)), np.nan)
                if not selected:
                    ax.text(.5,.5,"Outside this arm's declared grid",transform=ax.transAxes,ha="center",fontsize=8)
                    ax.axis("off");continue
                for row in selected:
                    y, x = sizes.index(num(row, "C_label")), strengths.index(num(row, "f"))
                    if key == "tau_p":
                        state = str(get(row, "tau_p_status", default=""))
                        value = num(row, "tau_p", "tau_p_estimate", "tau_p_point")
                        if not math.isfinite(value): value = num(row, "tau_p_lower")
                        data[y, x] = math.log10(value) if value > 0 else math.nan
                        if "censored" in state: ax.text(x, y, "↑", ha="center", va="center", fontsize=10)
                        if "unresolved" in state: ax.text(x, y, "↓", ha="center", va="center", fontsize=10)
                    elif key == "forced":
                        state = forced(row); data[y, x] = math.nan if state is None else float(state)
                    else:
                        aliases = {"tau_coverage": ("tau_p_coverage", "coverage"), "single_fraction": ("fraction_single",),
                                   "multiple_fraction": ("fraction_multiple",), "low_contrast_fraction": ("fraction_low_contrast",)}
                        data[y, x] = num(row, key, *aliases[key])
                    if forced(row) is True and key != "forced":
                        ax.add_patch(plt.Rectangle((x-.47, y-.47), .94, .94, facecolor="#cccccc", alpha=.8, ec="#777777", lw=.8))
                image = ax.imshow(data, aspect="auto", origin="lower", cmap="viridis", vmin=None if logarithmic else 0,
                                  vmax=None if logarithmic else 1)
                ax.set_xticks(range(len(strengths)), [f"{f:g}" for f in strengths], rotation=45, fontsize=6)
                ax.set_yticks(range(len(sizes)), [f"{c:g}" for c in sizes], fontsize=7)
                ax.set_title((LABELS[which_model] + " · " if j == 0 else "") + title, fontsize=8)
                fig.colorbar(image, ax=ax, fraction=.025, pad=.02).ax.tick_params(labelsize=6)
                if i == 3: ax.set_xlabel("Input fraction f", fontsize=8)
                if j == 0: ax.set_ylabel("Resting C", fontsize=8)
            page_rows.extend(selected)
        fig.subplots_adjust(left=.065, right=.96, top=.92, bottom=.15, hspace=.8, wspace=.45)
        pages.save(fig, f"S2_{index}", "Condition overview" + (" · log₁₀ native persistence" if index == 1 else " · state fractions"), page_rows,
                   "All declared input strengths in this observation arm are shown. Grey overlays retain forced conditions; the exact display flag and its components remain in CSV. Arrows preserve censoring and cadence limits.")


def strength_batches(tables):
    values = sorted({num(r, "f") for r in main_rows(tables["condition_metrics.csv"])})
    return [values[start:start + 3] for start in range(0, len(values), 3)]


def heading_pages(pages, tables):
    rows = main_rows(tables["heading_curves.csv"])
    conditions = {condition_key(r): r for r in main_rows(tables["condition_metrics.csv"])}
    for batch_index, strengths in enumerate(strength_batches(tables), 1):
        for metric in ("cosine", "angular_msd"):
            fig, axes = plt.subplots(4, 3, figsize=(11.7, 8.3), squeeze=False)
            page_rows = []
            for i, which_model in enumerate(MODELS):
                for j in range(3):
                    ax = axes[i, j]
                    if j >= len(strengths): ax.axis("off"); continue
                    strength = strengths[j]; selected = subset(rows, which_model=which_model, strength=strength)
                    if not selected: ax.text(.5, .5, "Outside declared grid", transform=ax.transAxes, ha="center", fontsize=8)
                    for size in sorted({num(r, "C_label") for r in selected}):
                        series = subset(selected, size=size); reference = conditions.get(condition_key(series[0]), {})
                        color = size_color(which_model, size, forced(reference) is True)
                        if metric == "cosine":
                            curve(ax, series, "lag", ("C", "cosine", "heading_cosine"), color=color, label=f"C {size:g}",
                                  lowkeys=("C_ci_low", "cosine_ci_low"), highkeys=("C_ci_high", "cosine_ci_high"))
                            coverage_rows = [dict(r, strip=-1.25 + .10 * num(r, "coverage", "pair_coverage")) for r in series]
                            curve(ax, coverage_rows, "lag", ("strip",), color=color, label=None)
                            fitted=num(reference,"exponential_tau")
                            if fitted>0:
                                xx=np.array(sorted({num(r,"lag") for r in series}))
                                ax.plot(xx,np.exp(-xx/fitted),color=color,lw=.6,ls=":",alpha=.7)
                            ax.axhline(math.exp(-1), color=".6", lw=.5, ls=":"); ax.set_ylim(-1.29, 1.05)
                        else:
                            short = [r for r in series if num(r, "lag") <= 10 * num(reference, "cadence", "scalar_cadence", default=.5 if which_model == "spring_ring" else 5)]
                            curve(ax, short, "lag", ("angular_msd", "heading_msd", "msd"), color=color, label=f"C {size:g}",
                                  lowkeys=("angular_msd_ci_low", "msd_ci_low"), highkeys=("angular_msd_ci_high", "msd_ci_high"))
                    ax.set_title(f"{LABELS[which_model]} · f={strength:g}", fontsize=8)
                    polish(ax, xlabel="Lag (native time)" if i == 3 else None,
                           ylabel=(r"$C(\Delta t)$ / coverage strip" if metric == "cosine" else r"Angular MSD (rad$^2$)") if j == 0 else None)
                    if i == 0 and ax.get_legend_handles_labels()[0]: ax.legend(fontsize=6, frameon=False, loc="best")
                    page_rows.extend(selected)
            fig.subplots_adjust(left=.08, right=.97, top=.92, bottom=.15, hspace=.65, wspace=.38)
            pages.save(fig, f"S3_{batch_index}_{metric}", "All heading correlations" if metric == "cosine" else "Short-lag angular MSD", page_rows,
                       "Lines/intervals aggregate whole seeds; dotted exponentials are comparisons, never empirical-bound replacements. All lag rows remain in CSV. The full cosine range [−1,1] remains visible; coverage strips span −1.25 to −1.15 for coverage 0 to 1. MSD tracks never bridge invalid headings.")


def size_pages(pages, tables):
    rows = main_rows(tables["condition_metrics.csv"])
    for batch_index, strengths in enumerate(strength_batches(tables), 1):
        for prefix, label in (("tau_p", r"Empirical $\tau_p$"), ("Dtheta", r"Short-lag $D_\theta$ (rad$^2$/time)")):
            fig, axes = plt.subplots(4, 3, figsize=(11.7, 8.3), squeeze=False); page_rows=[]
            for i, which_model in enumerate(MODELS):
                for j in range(3):
                    ax=axes[i,j]
                    if j>=len(strengths): ax.axis("off"); continue
                    selected=sorted(subset(rows,which_model=which_model,strength=strengths[j]),key=lambda r:num(r,"C_label"))
                    for row in selected:
                        x=num(row,"actual_C_mean","actualC_mean","C_mean","C_label")
                        color=size_color(which_model,num(row,"C_label"),forced(row) is True)
                        if prefix=="tau_p": interval_mark(ax,x,row,prefix,color=color)
                        else:
                            y=num(row,"Dtheta","D_theta");lo=num(row,"Dtheta_ci_low","D_theta_ci_low");hi=num(row,"Dtheta_ci_high","D_theta_ci_high")
                            if math.isfinite(y):
                                ax.plot(x,y,"o",color=color,ms=4)
                                if math.isfinite(lo) and math.isfinite(hi):ax.vlines(x,lo,hi,color=color,lw=1)
                    if prefix=="tau_p" and selected:
                        first=selected[0];anchor=num(first,"tau_p","tau_p_estimate","tau_p_point")
                        x0=num(first,"actual_C_mean","actualC_mean","C_mean","C_label")
                        if "resolved"==get(first,"tau_p_status") and anchor>0:
                            xs=np.array([num(r,"actual_C_mean","actualC_mean","C_mean","C_label") for r in selected])
                            ax.plot(xs,anchor*(xs/x0)**2,color=".65",ls=":",lw=.8,label=r"$C^2$ reference")
                    if not selected:ax.text(.5,.5,"Outside declared grid",transform=ax.transAxes,ha="center",fontsize=8)
                    ax.set_title(f"{LABELS[which_model]} · f={strengths[j]:g}",fontsize=8)
                    polish(ax,xlabel="Actual mean circumference" if i==3 else None,ylabel=label if j==0 else None,logy=prefix=="tau_p")
                    page_rows.extend(selected)
            fig.subplots_adjust(left=.08,right=.97,top=.92,bottom=.15,hspace=.65,wspace=.38)
            pages.save(fig,f"S4_{batch_index}_{prefix}","Persistence and angular mobility across size",page_rows,
                       "Circumference is measured in the analysis window. Dθ is a signed short-lag mobility estimate, not an assumption of free diffusion. Dotted C² reference is anchored only to a resolved smallest-size persistence.")


def ratio_page(pages,tables):
    rows=main_rows(tables["size_ratios.csv"])
    fig,axes=plt.subplots(4,2,figsize=(11.7,8.3),squeeze=False)
    for i,which_model in enumerate(MODELS):
        selected=subset(rows,which_model=which_model)
        for j,prefix in enumerate(("tau_ratio","excess_persistence")):
            ax=axes[i,j]
            for row in selected:
                interval_mark(ax,num(row,"f"),row,prefix,color="#888888" if forced(row) is True else "#333333")
                if j==0:
                    geometric=num(row,"geometric_ratio","actual_geometric_ratio")
                    if math.isfinite(geometric):ax.plot(num(row,"f"),geometric,"_",color="#34974a",ms=8)
                    proxy=num(row,"Dtheta_mobility_proxy")
                    if math.isfinite(proxy) and proxy>0 and get(row,"tau_ratio_status")!="resolved":
                        proxy_color="#888888" if forced(row) is True else "#2878b5"
                        ax.plot(num(row,"f"),proxy,"D",mec=proxy_color,mfc="none",ms=4)
                        lo=num(row,"Dtheta_mobility_proxy_ci_lower");hi=num(row,"Dtheta_mobility_proxy_ci_upper")
                        if math.isfinite(lo) and math.isfinite(hi):ax.vlines(num(row,"f"),lo,hi,color=proxy_color,lw=.7)
            ax.axhline(1,color=".6",ls=":",lw=.7)
            if j==0:ax.axhline(9,color=".75",ls="--",lw=.6)
            polish(ax,xlabel="Input fraction f" if i==3 else None,ylabel=LABELS[which_model],logx=True,logy=True)
            if i==0:ax.set_title("Large/small empirical persistence" if j==0 else "Ratio / actual geometric expectation",fontsize=9)
    fig.subplots_adjust(left=.08,right=.97,top=.92,bottom=.15,hspace=.55,wspace=.3)
    pages.save(fig,"S4_ratios","Size advantage and geometric reference",rows,
               "Triangles preserve finite-window or resolution bounds. Two censored times do not yield a finite ratio. Open diamonds show the short-lag global activity-heading small/large Dθ proxy: blue for unforced conditions and grey for forced conditions. This proxy does not resolve a direct persistence ratio or measure one tracked site's diffusion in multi-region states. Green marks: measured circumference-ratio squared; dashed 9: resting-size reference.")


def site_pages(pages,tables):
    rows=main_rows(tables["site_curves.csv"])
    references={condition_key(r):r for r in main_rows(tables["condition_metrics.csv"])}
    for batch_index,strengths in enumerate(strength_batches(tables),1):
        fig,axes=plt.subplots(4,3,figsize=(11.7,8.3),squeeze=False);page_rows=[]
        for i,which_model in enumerate(MODELS):
            for j in range(3):
                ax=axes[i,j]
                if j>=len(strengths):ax.axis("off");continue
                selected=subset(rows,which_model=which_model,strength=strengths[j])
                for size in sorted({num(r,"C_label") for r in selected}):
                    series=subset(selected,size=size)
                    color=size_color(which_model,size,forced(references.get(condition_key(series[0]),{})) is True)
                    for method,ls in (("heading_arc","-"),("heading_direct","--"),("profile_peak",":")):
                        method_rows=[r for r in series if get(r,"method",default="heading_arc")==method]
                        if not method_rows:continue
                        curve(ax,method_rows,"lag",("msd","site_msd","physical_msd"),color=color,
                              label=f"C {size:g}" if method=="heading_arc" else None,
                              lowkeys=("msd_ci_low","site_msd_ci_low"),highkeys=("msd_ci_high","site_msd_ci_high"),ls=ls)
                if not selected:ax.text(.5,.5,"Outside declared grid",transform=ax.transAxes,ha="center",fontsize=8)
                ax.set_title(f"{LABELS[which_model]} · f={strengths[j]:g}",fontsize=8)
                polish(ax,xlabel="Lag (native time)" if i==3 else None,ylabel="Physical site MSD" if j==0 else None)
                if i==0 and ax.get_legend_handles_labels()[0]:ax.legend(fontsize=6,frameon=False)
                page_rows.extend(selected)
        fig.subplots_adjust(left=.08,right=.97,top=.92,bottom=.15,hspace=.65,wspace=.38)
        pages.save(fig,f"S5_{batch_index}_siteMSD","Physical-coordinate estimator sensitivities",page_rows,
                   "These inherited estimators are sensitivities, not the literal fixed-material-origin A3 primary. Solid: midpoint-C increment integration. Dashed: track-reanchored (θ−θstart)C(t)/(2π). Dotted: profile-maximum midpoint-C increments at the retained profile cadence. Exact cadence is in CSV and can equal scalar cadence in dense companions. Tracks reset at invalid frames; a plateau can indicate confinement. Literal native-M1 θC(t)/(2π) and antipodal ambiguity checks are displayed separately.")
    physical_mobility_page(pages,tables)


def physical_mobility_page(pages,tables):
    """Show the declared angular-to-length conversion, separately from MSD."""
    conditions=main_rows(tables["condition_metrics.csv"])
    fig,axes=plt.subplots(2,2,figsize=(11.7,8.3));allrows=[]
    for ax,which_model in zip(axes.flat,MODELS):
        selected=subset(conditions,which_model=which_model)
        for size in sorted({num(r,"C_label") for r in selected}):
            series=sorted(subset(selected,size=size),key=lambda r:num(r,"f"));x=[num(r,"f") for r in series];y=[num(r,"Dx","D_x") for r in series]
            color=size_color(which_model,size);ax.plot(x,y,"-",color=color,lw=.8,label=f"C {size:g}")
            for xx,yy,row in zip(x,y,series):ax.plot(xx,yy,"o",color=size_color(which_model,size,forced(row)is True),ms=3)
            lo=[num(r,"Dx_ci_low","D_x_ci_low") for r in series];hi=[num(r,"Dx_ci_high","D_x_ci_high") for r in series]
            ax.fill_between(x,lo,hi,color=color,alpha=.13)
        ax.set_title(LABELS[which_model],fontsize=10);polish(ax,xlabel="Input fraction f",ylabel=r"$D_x$ global-heading proxy"+"\n(native length²/time)",logx=True)
        ax.legend(fontsize=7,frameon=False);allrows.extend(selected)
    fig.subplots_adjust(left=.09,right=.97,top=.84,bottom=.22,hspace=.4,wspace=.35)
    pages.save(fig,"S5_Dx","Global-heading mobility proxy",allrows,
               "Dx = Dθ × (actual mean C / 2π)² converts the signed short-lag global activity-heading estimate using each window's actual mean circumference. In multi-region states this proxy is not diffusion of one tracked site. It is distinct from qualified single-region literal heading/peak coordinate MSD and is not a slope fit to θC(t)/(2π). Increased Dx and angular persistence advantages are both retained when they coexist. Grey retains forced conditions. Native length/time units differ across models. No free-diffusion or migration velocity law is assumed.")


def noise_tolerance_pages(pages,tables):
    rows=main_rows(tables["condition_metrics.csv"])
    fig,axes=plt.subplots(2,2,figsize=(11.7,8.3))
    for ax,which_model in zip(axes.flat,MODELS):
        selected=subset(rows,which_model=which_model)
        bounds=[]
        for size in sorted({num(r,"C_label") for r in selected}):
            for row in sorted(subset(selected,size=size),key=lambda r:num(r,"f")):
                first=len(ax.lines)
                state=interval_mark(ax,num(row,"f"),row,"tau_p",color=size_color(which_model,size,forced(row)is True))
                if state in ('bounded','unresolved'):
                    bounds.append((row,state,ax.lines[first]))
        coincident_bound_markers(ax,bounds,'tau_p')
        ax.set_title(LABELS[which_model],fontsize=10);polish(ax,xlabel="Input fraction f",ylabel=r"Empirical $\tau_p$ (native time)",logx=True,logy=True)
        ax.margins(y=.25)
    legend_mapping=shared_size_legend(fig,rows,y=.135)
    for row in rows:row['display_shared_size_legend']=json.dumps(legend_mapping)
    fig.subplots_adjust(left=.08,right=.97,top=.86,bottom=.27,hspace=.65,wspace=.3)
    pages.save(fig,"S6_tau","Empirical persistence against input strength",rows,
               "Upward triangles are empirical lower bounds; downward triangles are unresolved below two cadences. Nested triangles at one unchanged coordinate retain coincident sizes; ×n gives their count. Their exact sizes, glyph sizes and coordinates remain in CSV. The shared legend maps chemical C to Spring resting C₀; grey retains forced observations. Bootstrap CI markers are separate from empirical bounds. Exponential fits never replace empirical measurements.")
    tolerance=[r for r in main_rows(tables["noise_tolerance.csv"]) if boolean(r.get("forced_inclusive")) is not True]
    fig,axes=plt.subplots(2,2,figsize=(11.7,8.3))
    for ax,which_model in zip(axes.flat,MODELS):
        selected=subset(tolerance,which_model=which_model)
        for row in selected:
            interval_mark(ax,num(row,"C_label"),row,"fhalf",color=size_color(which_model,num(row,"C_label")))
        ax.set_title(LABELS[which_model],fontsize=10);polish(ax,xlabel="Resting circumference",ylabel=r"Half-persistence input $f_{1/2}$",logy=True)
        if all(not math.isfinite(num(r,"fhalf","fhalf_estimate","fhalf_point")) for r in selected):
            ax.text(.05,.9,"No resolved point estimate; see bounds/status in CSV",transform=ax.transAxes,fontsize=8,va="top")
    fig.subplots_adjust(left=.08,right=.97,top=.91,bottom=.14,hspace=.4,wspace=.3)
    pages.save(fig,"S6_fhalf","Noise tolerance with finite-window uncertainty",tolerance,
               "First crossing of half the actual f=0.0001 persistence, using log-f interpolation only when identified. Censored references remain bounded or inconclusive. Intermediate-size grids lacking the declared f=0.0001 baseline retain unavailable markers/reasons; f=0.001 is never substituted and no tolerance curve bridges them.")


def coincident_bound_markers(ax,entries,prefix,*,x_field='f'):
    """Expose exactly coincident bounds with unchanged scientific coordinates."""
    groups={}
    for row,state,artist in entries:
        x=num(row,x_field);y=num(row,prefix+'_lower' if state=='bounded' else prefix+'_upper')
        groups.setdefault((x,y,state),[]).append((row,artist))
    for (x,y,state),items in groups.items():
        if len(items)<2:continue
        ordered=sorted(items,key=lambda item:num(item[0],'C_label'))
        sizes=[num(row,'C_label') for row,_ in ordered]
        if len(sizes)!=len(set(sizes)):raise ValueError('One condition/size has duplicate empirical bounds')
        for index,(row,artist) in enumerate(ordered):
            artist.remove()
            scale=4+2*index
            ax.plot(x,y,'^' if state=='bounded' else 'v',ms=scale,mfc='none',
                    mec=size_color(model(row),num(row,'C_label'),forced(row)is True),mew=.8,label=artist.get_label())
            row['display_'+prefix+'_coincident_bound_count']=len(ordered)
            row['display_'+prefix+'_coincident_bound_sizes']=json.dumps(sizes)
            row['display_'+prefix+'_bound_marker_points']=scale
            row['display_'+prefix+'_bound_x_field']=x_field
            row['display_'+prefix+'_scientific_coordinate_unchanged']=True
        ax.annotate('×'+str(len(ordered)),xy=(x,y),xytext=(0,7+len(ordered)),textcoords='offset points',
                    ha='center',fontsize=6,color='#333333',annotation_clip=False)


def shared_size_legend(fig,rows,y):
    """Use the actual size mapping; no size is inferred from a response."""
    grids={which:sorted({num(row,'C_label') for row in subset(rows,which_model=which)}) for which in MODELS}
    chemical=[grids[which] for which in MODELS if which!='spring_ring']
    if any(grid!=chemical[0] for grid in chemical) or len(chemical[0])!=len(grids['spring_ring']):
        raise ValueError('A shared chemical/Spring size legend needs the actual matching designed grids')
    handles=[Line2D([],[],marker='o',lw=0,color=size_color('wave_pinning',size),
                    label=f'{size:g} ({spring:g})') for size,spring in zip(chemical[0],grids['spring_ring'])]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.52,y),ncol=len(handles),
               frameon=False,fontsize=7.5,title='Circumference: chemical C (Spring resting C₀)',title_fontsize=7.5)
    return [dict(chemical_C=size,Spring_resting_C0=spring,color=size_color('wave_pinning',size))
            for size,spring in zip(chemical[0],grids['spring_ring'])]


def robustness_pages(pages,tables):
    rows=tables["size_ratios.csv"]
    view_groups=(("main","analysis_first_half","analysis_second_half","late_quarter"),
                 ("main","cadence_x2","cadence_x4","cadence_x8"))
    styles=("-","--",":","-.")
    for index,views in enumerate(view_groups,1):
        fig,axes=plt.subplots(2,2,figsize=(11.7,8.3));page_rows=[]
        for ax,which_model in zip(axes.flat,MODELS):
            for view,ls in zip(views,styles):
                selected=sorted(subset(rows,which_model=which_model,view=view),key=lambda r:num(r,"f"))
                for row in selected:
                    result=interval_mark(ax,num(row,"f"),row,"tau_ratio",color="#888888" if forced(row)is True else "#333333")
                resolved=[r for r in selected if get(r,"tau_ratio_status")=="resolved"]
                ax.plot([num(r,"f") for r in resolved],[num(r,"tau_ratio","tau_ratio_estimate") for r in resolved],ls=ls,lw=.8,color=".3",label=view.replace("analysis_",""))
                page_rows.extend(selected)
            ax.axhline(1,color=".6",ls=":",lw=.6);ax.set_title(LABELS[which_model],fontsize=10)
            polish(ax,xlabel="Input fraction f",ylabel="Large/small empirical persistence",logx=True,logy=True);ax.legend(fontsize=6,frameon=False)
        fig.subplots_adjust(left=.08,right=.97,top=.91,bottom=.14,hspace=.4,wspace=.3)
        pages.save(fig,f"S12_{index}","Window robustness" if index==1 else "Scalar-cadence robustness",page_rows,
                   "Each view is recomputed from its own valid pairs; the bootstrap unit remains the whole seed. Grey marks are retained forced conditions. E4 extended horizons are pending, so no extrapolated correction is shown.")
    decisions=tables["decisions.csv"]
    for inclusive in (False,True):
        selected=[r for r in decisions if boolean(r.get("forced_inclusive")) is inclusive]
        if not selected:continue
        fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis("off")
        cells=[[LABELS.get(model(r),model(r)),r.get("view","main").replace("analysis_","").replace("_"," "),outcome(r.get("H1","inconclusive")).replace("criterion provisionally met","provisional met"),outcome(r.get("H2","inconclusive")).replace("model criterion provisionally met","provisional met"),
                r.get("evaluable_conditions",""),reversal_text(r.get("reversal_f","")),yesno(r.get("fhalf_increase_nonoverlap"))]
               for r in selected]
        table=ax.table(cellText=cells,colLabels=["Model","Window / cadence","H1","H2","Evaluable n","Reversal f","f½ CIs separate"],
                       bbox=[0,.08,1,.82],cellLoc="left",colWidths=[.13,.23,.12,.12,.11,.16,.13])
        table.auto_set_font_size(False);table.set_fontsize(7)
        fig.subplots_adjust(left=.04,right=.97,top=.94,bottom=.10)
        pages.save(fig,"S12_decisions_"+("inclusive" if inclusive else "nonforced"),
                   "Forced-inclusive decision sensitivity" if inclusive else "Non-forced decision sensitivity",selected,
                   "Provisional criteria are recomputed for every window/cadence. E1 forced-state continuation and the A6 overall H2 dependency are unavailable in E0.")
    triggers=tables["e4_trigger_inventory.csv"]
    fig,axes=plt.subplots(2,2,figsize=(11.7,8.3))
    for ax,which_model in zip(axes.flat,MODELS):
        selected=subset(triggers,which_model=which_model)
        for row in selected:
            size=num(row,"C_label");is_trigger=boolean(row.get("trigger"))
            mark="o" if is_trigger else "x" if "no_usable" in row.get("status","") else "s"
            color=size_color(which_model,size,forced(row)is True)
            ax.plot(num(row,"f"),size,marker=mark,mec=color,mfc=color if is_trigger else "none",ms=5,ls="none")
        n=sum(boolean(r.get("trigger")) is True for r in selected)
        ax.set_title(f"{LABELS[which_model]} · {n} horizon-extension triggers",fontsize=10)
        polish(ax,xlabel="Input fraction f",ylabel="Resting circumference",logx=True)
        ax.set_yticks(sorted({num(r,"C_label") for r in selected}))
    fig.subplots_adjust(left=.08,right=.97,top=.91,bottom=.14,hspace=.4,wspace=.3)
    pages.save(fig,"S12_E4_pending","Finite-horizon limitations and pending extension inventory",triggers,
               "Filled circles: cosine >0.8 at the longest usable lag; open squares: no trigger; crosses: no usable lag. All 110 groups remain shown, including forced conditions in grey. The 52 trigger groups require 624 proposed integrations at 12 seeds, and fourfold horizons do not guarantee uncensoring. These integrations have not entered E0.")


def write_reader(output,report,tables):
    output=Path(output)
    links="\n".join(f'<li>Page {p["pdf_page"]}: {html.escape(p["page"])} — {html.escape(p["title"])} · <a href="supporting/data/{Path(p["csv"]).name}">CSV</a></li>' for p in report["pages"])
    content=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 4 E0 interim supplement</title>
<style>body{{font:17px/1.5 system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem}}a{{color:#16699d}}.note{{background:#fff2ce;padding:1rem}}li{{margin:.35rem 0}}</style>
<h1>Figure 4 · E0 interim supplement</h1><p class="note">A1–A4 and A7 from all existing native noisy-ring cases. Hypothesis decisions are provisional. E1–E4, A5 and A6 remain pending; a final main-figure framing has not been selected.</p>
<p>Generated {html.escape(report["generated_utc"])}. Renderer stopped. <a href="supporting/vectors/figure4_E0_interim_supplement.pdf">Open the {len(report["pages"])}-page interim PDF</a> · <a href="supporting/metadata/interim_display_audit.json">Input and display hashes</a></p>
<p>Every page links its complete plotted rows, including unavailable and censored measurements. Whole-seed uncertainty, valid-heading coverage, native clocks and actual Spring circumferences remain explicit.</p><ol>{links}</ol></html>'''
    (output/"index.html").write_text(content)
    (output/"README.md").write_text("# Figure 4 E0 interim supplement\n\n[Reader](index.html) · [Interim PDF](supporting/vectors/figure4_E0_interim_supplement.pdf)\n\nGenerated "+report["generated_utc"]+". Renderer stopped. All outputs are an interim E0 analysis; final hypotheses and Figure 4 framing await E1–E4, A5 and A6. Every page has an unabridged plotted-value CSV.\n\nThe original lognormal sigma=1 forcing is retained. Grey conditions have abnormal matched-mean controls; empirical finite-window and cadence bounds remain explicit. CIs use 2,000 whole-seed bootstrap draws. Spring physical motion uses its actual circumference. Scalar heading and coarse profile cadences are distinct.\n")


def render(run_root=DEFAULT_RUN,output=DEFAULT_OUT):
    tables,inventory,sources,status=load_inputs(run_root)
    plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"svg.fonttype":"none","axes.titlesize":9})
    pages=Pages(output,sources)
    try:
        decision_pages(pages,tables)
        inventory_page(pages,inventory,tables)
        regime_pages(pages,tables)
        heading_pages(pages,tables)
        size_pages(pages,tables)
        ratio_page(pages,tables)
        site_pages(pages,tables)
        noise_tolerance_pages(pages,tables)
        robustness_pages(pages,tables)
        report=pages.close()
    except Exception:
        pages.pdf.close()
        raise
    write_reader(output,report,tables)
    files=[Path(report["pdf"]["path"]),*(Path(p["csv"]) for p in report["pages"]),
           Path(output)/"index.html",Path(output)/"figure4_E0_interim_preview.png",
           Path(output)/"supporting/metadata/interim_display_audit.json"]
    write_json(Path(output)/"supporting/metadata/artifact_sha256.json",{str(path.relative_to(Path(output))):sha(path) for path in files})
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root",type=Path,default=DEFAULT_RUN)
    parser.add_argument("--output",type=Path,default=DEFAULT_OUT)
    args=parser.parse_args()
    report=render(args.run_root,args.output)
    print(json.dumps(dict(pdf=report["pdf"],interim=True),indent=2))


if __name__=="__main__":
    main()
