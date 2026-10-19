"""Complete-condition supplement from immutable native experiment readouts.

This renderer neither integrates native models nor changes the guarded terminal
pipeline. Pending conditions remain in every requested grid. The authoritative
tables/NPZ retain exact values; display JSON uses eight significant digits.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import html
import json
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[2]
GROUP = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
TOPIC = GROUP / "narrow_arc_experiment"
OUT = TOPIC / "supplement"
RUN = ROOT / "data/runs/revision_arc_narrow_grid_20261019"
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes", spring="Spring")
REGIMES = {"inside": "I", "straddling": "S", "outside": "O"}
COLORS = dict(inside="#ba5b36", straddling="#b79a42", outside="#237b83")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recorded_source_digest(recorded, path):
    """Resolve the same exported source after a repository checkout moves."""
    relative = Path(path).resolve().relative_to(ROOT).as_posix()
    matches = [digest for name, digest in recorded.items()
               if name == relative or name.endswith("/" + relative)]
    assert len(matches) == 1, ("ambiguous or missing source identity", relative)
    return matches[0]


def atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(text)
    temp.replace(path)


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(format(float(value), ".8g")) if np.isfinite(value) else None
    return value


def is_complete(row):
    return row["status"] == "complete"


def size_label(model, size):
    return ("C₀=" if model == "spring" else "C=") + f"{size:g}"


def footer(fig, status, extra=""):
    fig.text(.045, .047, f"All requested placements and durations retained · {status['complete']}/645 conditions available · data {status['generated_at'][:19]} UTC", fontsize=9)
    fig.text(.045, .027, "I / S / O = entire initial cue footprint inside / straddling / outside the 10%-contrast region. Centre membership is recorded separately.", fontsize=8.7)
    if extra:
        fig.text(.045, .010, extra, fontsize=8.5)


def title(fig, model, subtitle):
    fig.text(.045, .966, f"{NAMES[model]} · {subtitle}", fontsize=18, weight="bold")


def format_grid(ax, part, durations, distances):
    first = part.drop_duplicates("physical_distance_at_onset").set_index("physical_distance_at_onset")
    ax.set_xlim(-.5, len(durations)-.5)
    ax.set_ylim(len(distances)-.5, -.5)
    ax.set_xticks(range(len(durations)), [f"{t:g}" for t in durations])
    ax.set_yticks(range(len(distances)), [f"{d:.5g}  {REGIMES[first.loc[d, 'footprint_relative_to_core_10']]}" for d in distances])
    ax.tick_params(length=0, labelsize=8.5)
    for spine in ax.spines.values():
        spine.set_color("#cbd2d6")


def vector_page(model, table, status):
    selected = table[table.model == model]
    sizes = sorted(selected.size_parameter.unique())
    durations = sorted(selected.cue_duration.unique())
    distances = sorted(selected.physical_distance_at_onset.unique())
    fig, axes = plt.subplots(1, 3, figsize=(15, 10))
    fig.subplots_adjust(left=.09, right=.99, bottom=.15, top=.82, wspace=.30)
    title(fig, model, "Every pulse: withdrawal and retained activity vectors")
    fig.text(.045, .925, "Orange: cue withdrawal. Black: after 13,000 cue-free units. Each arrow is the actual P; one common scale in every cell.", fontsize=10)
    fig.text(.045, .900, "Background: fraction of initial target-angle error removed at the retained endpoint (colour clipped to 0–1 only). Grey: pending / weak final R.", fontsize=9.5)
    quivers, weak, pending = 0, 0, 0
    for ax, size in zip(axes, sizes):
        part = selected[selected.size_parameter == size]
        format_grid(ax, part, durations, distances)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(size_label(model, size), loc="left", fontsize=12)
        for r in part.to_dict("records"):
            x, y = durations.index(r["cue_duration"]), distances.index(r["physical_distance_at_onset"])
            done = is_complete(r)
            val = r["retained_error_closed"] if done else np.nan
            face = plt.get_cmap("Blues")(.10+.75*np.clip(val, 0, 1)) if np.isfinite(val) else "#f0f0f0"
            ax.add_patch(Rectangle((x-.48, y-.48), .96, .96, facecolor=face, edgecolor="white", linewidth=.5))
            if not done:
                ax.text(x, y, "·", ha="center", va="center", color="#777")
                pending += 1
                continue
            ax.plot([x, x+.30], [y, y], color="#a9b1b7", lw=.55, zorder=2)
            for endpoint, color, width in (("withdrawal", "#db852d", .010), ("final", "#17242b", .014)):
                re, im, R = r[f"P_{endpoint}_real"], r[f"P_{endpoint}_imag"], r[f"R_{endpoint}"]
                if R >= .1:
                    ax.quiver(x, y, .30*re, -.30*im, angles="xy", scale_units="xy", scale=1, color=color, width=width, headwidth=3.5, headlength=4, zorder=4)
                    quivers += 1
                else:
                    ax.add_patch(Circle((x, y), .30*R, fill=False, edgecolor=color, linewidth=1.0, zorder=5))
                    weak += 1
            tail = r["full_field_tail_status"]
            if tail != "passes_finite_tail_gate":
                ax.text(x+.34, y-.32, "?" if "unavailable" in str(tail) else "!", fontsize=8, color="#a92323", ha="center", va="center")
        ax.set_xlabel("Tested cue duration (categorical native units)", fontsize=9)
    axes[0].set_ylabel("Initial physical centre distance d₀ and footprint class", fontsize=10)
    footer(fig, status, "Rightward guide is P(0)=1; circle = measured R<0.1 without a heading. ? full-field tail unavailable; ! failed finite-tail gate. No minimum-pulse inference.")
    return fig, dict(quivers=quivers, weak_circles=weak, pending_cells=pending, vector_scale=.30, equal_xy=True)


MATRIX_GROUPS = [
    ("cue_response", "Every pulse: cue-on response and sampled coverage", [
        ("withdrawal_error_closed", "Target-angle error removed at withdrawal", "Blues", (0, 1)),
        ("R_cue_on_sampled_minimum", "Minimum saved cue-on R", "viridis", (0, 1)),
        ("rhs_withdrawal", "Native RHS angular rate at withdrawal (° / native unit)", "coolwarm", None),
        ("secant_coverage", "Fraction of saved cue-on intervals with usable phase", "viridis", (0, 1))]),
    ("amount_region", "Every pulse: vector magnitude, active amount and fixed-region response", [
        ("R_withdrawal", "R at cue withdrawal", "viridis", (0, 1.5)),
        ("M0_ratio_withdrawal", "Total active amount / initial amount at withdrawal", "coolwarm", None),
        ("fixed_region_excess_cue_on_sampled_minimum", "Minimum saved cue-on original-region excess", "viridis", (0, 1)),
        ("fixed_region_excess_final", "Original-region excess after 13,000 cue-free units", "viridis", (0, 1))]),
    ("retained_organization", "Every pulse: finite cue-free endpoint organization", [
        ("retained_error_closed", "Target-angle error removed after 13,000 cue-free units", "Blues", (0, 1)),
        ("R_final", "R at finite cue-free endpoint", "viridis", (0, 1.5)),
        ("final_resolved_peaks", "Resolved final active-profile peaks", "viridis", (0, 3)),
        ("final_halfcontrast_regions", "Final half-contrast connected regions", "viridis", (0, 3))])]


def matrix_page(model, table, status, group):
    slug, subtitle, metrics = group
    selected = table[table.model == model]
    sizes = sorted(selected.size_parameter.unique())
    durations = sorted(selected.cue_duration.unique())
    distances = sorted(selected.physical_distance_at_onset.unique())
    fig, axes = plt.subplots(len(metrics), 3, figsize=(15, 17))
    fig.subplots_adjust(left=.105, right=.925, top=.89, bottom=.11, wspace=.35, hspace=.46)
    title(fig, model, subtitle)
    fig.text(.045, .932, "Five tested cue durations × every prescribed distance × three exact qualified basals. Grey = pending or the specific quantity is unavailable.", fontsize=10)
    for row, (metric, label, cmap, limits) in enumerate(metrics):
        values = pd.to_numeric(selected.loc[selected.status == "complete", metric], errors="coerce").dropna().to_numpy()
        if limits is None:
            if metric == "M0_ratio_withdrawal":
                bound = max(.01, np.max(abs(values-1)) if len(values) else .01)
                norm = TwoSlopeNorm(vmin=1-bound, vcenter=1, vmax=1+bound)
            else:
                bound = max(.001, np.max(abs(values)) if len(values) else .001)
                norm = TwoSlopeNorm(vmin=-bound, vcenter=0, vmax=bound)
        else:
            if metric in ("R_withdrawal", "R_final"):
                limits = (0, max(1.1, float(np.max(values))*1.01 if len(values) else 1.1))
            if metric in ("final_resolved_peaks", "final_halfcontrast_regions"):
                limits = (0, max(3, float(np.max(values)) if len(values) else 3))
            norm = Normalize(*limits)
        for col, size in enumerate(sizes):
            ax = axes[row, col]
            part = selected[selected.size_parameter == size]
            arr = np.full((len(distances), len(durations)), np.nan)
            for r in part.to_dict("records"):
                if is_complete(r):
                    arr[distances.index(r["physical_distance_at_onset"]), durations.index(r["cue_duration"])] = r[metric]
            cm = plt.get_cmap(cmap).copy(); cm.set_bad("#e9ecee")
            image = ax.imshow(arr, cmap=cm, norm=norm, aspect="auto", interpolation="nearest")
            format_grid(ax, part, durations, distances)
            ax.set_title((size_label(model, size) + " · " if row == 0 else "") + label, fontsize=9.5, loc="left", wrap=True)
            ax.set_xlabel("Cue duration (categorical)", fontsize=8.7)
            if col == 0:
                ax.set_ylabel("d₀ and footprint I/S/O", fontsize=9)
        cbax = fig.add_axes([.945, axes[row, 2].get_position().y0, .010, axes[row, 2].get_position().height])
        fig.colorbar(image, cax=cbax).ax.tick_params(labelsize=8)
    footer(fig, status, "Error-closure colours clip to 0–1; exact raw values remain in CSV. Fixed original-region loss is redistribution-compatible, not site ancestry or destruction.")
    return fig


def trajectory_page(model, table, traces, status):
    chosen = table[(table.model == model) & (table.cue_duration == 1000.)]
    sizes = sorted(chosen.size_parameter.unique())
    distances = sorted(chosen.physical_distance_at_onset.unique())
    colors = {d: plt.get_cmap("turbo")(.04+.92*j/max(1, len(distances)-1)) for j, d in enumerate(distances)}
    fig, axes = plt.subplots(5, 3, figsize=(15, 17))
    fig.subplots_adjust(left=.09, right=.985, bottom=.12, top=.85, hspace=.32, wspace=.24)
    title(fig, model, "Every placement: full longest-cue time courses")
    fig.legend([Line2D([0], [0], color=colors[d], label=f"d₀={d:.5g}") for d in distances], [f"d₀={d:.5g}" for d in distances], loc="upper left", bbox_to_anchor=(.04, .935), ncol=6, frameon=False, fontsize=9)
    data_keys = ["angle_masked_deg", "R", "M0_ratio", "fixed_region_transported_background_excess"]
    common = {key: [] for key in data_keys[1:]+["cue_sampled_angular_secant"]}
    for r in chosen.to_dict("records"):
        if is_complete(r):
            prefix = f"c{int(r['condition_id'])}_"; mask = traces[prefix+"phase"] == 0
            for key in common:
                values = traces[prefix+key] if key == "cue_sampled_angular_secant" else traces[prefix+key][mask]
                common[key].extend(values[np.isfinite(values)])
    for col, size in enumerate(sizes):
        part = chosen[chosen.size_parameter == size]
        axes[0, col].set_title(size_label(model, size)+f" · {(part.status=='complete').sum()}/{len(part)} histories", loc="left", fontsize=11)
        for r in part.to_dict("records"):
            if not is_complete(r):
                continue
            cid = int(r["condition_id"]); prefix = f"c{cid}_"; mask = traces[prefix+"phase"] == 0
            time = traces[prefix+"time"][mask]
            style = {"inside": "-", "straddling": "--", "outside": ":"}[r["footprint_relative_to_core_10"]]
            for k, key in enumerate(data_keys):
                axes[k, col].plot(time, traces[prefix+key][mask], color=colors[r["physical_distance_at_onset"]], ls=style, lw=1.0)
            left, right, rate = (traces[prefix+k] for k in ("cue_rate_left", "cue_rate_right", "cue_sampled_angular_secant"))
            seg = [[(a, v), (b, v)] for a, b, v in zip(left, right, rate) if np.isfinite(v)]
            if seg:
                axes[4, col].add_collection(LineCollection(seg, colors=[colors[r["physical_distance_at_onset"]]], linewidths=.7))
                axes[4, col].autoscale_view()
        for row in range(5):
            ax = axes[row, col]
            ax.set_xscale("symlog", linthresh=.1); ax.set_xlim(0, 1000)
            ax.set_xticks([0, .1, 1, 10, 100, 1000], ["0", ".1", "1", "10", "100", "1000"])
            ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=8.5); ax.minorticks_off()
        axes[0, col].set_ylim(-5, 185)
        axes[1, col].axhline(.1, color="#888", lw=.6, ls="--")
        axes[1, col].set_ylim(0, max(1.2, max(common["R"], default=1)*1.04))
        amount_bound = max(.01, max((abs(v-1) for v in common["M0_ratio"]), default=0)*1.06)
        axes[2, col].set_ylim(1-amount_bound, 1+amount_bound)
        region_values = common["fixed_region_transported_background_excess"]
        axes[3, col].set_ylim(min(0, min(region_values, default=0))-.05, max(1, max(region_values, default=1))+.05)
        axes[4, col].set_yscale("symlog", linthresh=.001)
        rate_bound = max(.001, max((abs(v) for v in common["cue_sampled_angular_secant"]), default=0)*1.05)
        axes[4, col].set_ylim(-rate_bound, rate_bound)
        axes[4, col].set_xlabel("Cue-on time (native units)", fontsize=10)
    for row, label in enumerate(["Heading (degrees)", "R = |P|", "Total active amount / initial", "Fixed original-region excess", "Sampled phase secant\n(degrees / native unit)"]):
        axes[row, 0].set_ylabel(label, fontsize=10)
    footer(fig, status, "Solid/dashed/dotted = I/S/O at initial 10% edge. Time linear to 0.1 then log; rates linear near ±0.001 then log. Each rate segment spans its actual saved interval.")
    return fig


def profile_pages(model, table, profiles, basal, pdf, status):
    selected = table[table.model == model]
    all_values = [v["u"] for (m, _), v in basal.items() if m == model]
    complete_ids = set(selected.loc[selected.status == "complete", "condition_id"].astype(int))
    subset = profiles[profiles.condition_id.isin(complete_ids)]
    ymin = min(np.min(v) for v in all_values)
    ymax = max(np.max(v) for v in all_values)
    if len(subset):
        ymin, ymax = min(ymin, subset.final_active.min()), max(ymax, subset.final_active.max())
    pad = .06*(ymax-ymin)
    by_id = {int(cid): part for cid, part in subset.groupby("condition_id")}
    for size in sorted(selected.size_parameter.unique()):
        initial = basal[(model, size)]
        for duration in sorted(selected.cue_duration.unique()):
            part = selected[(selected.size_parameter == size) & (selected.cue_duration == duration)]
            fig, axes = plt.subplots(3, 4, figsize=(15, 11))
            fig.subplots_adjust(left=.065, right=.99, bottom=.12, top=.85, hspace=.55, wspace=.27)
            title(fig, model, f"Every final profile · {size_label(model, size)} · cue T={duration:g}")
            fig.text(.045, .925, "Dotted grey: initial native basal. Black: measured final active field. Blue: basal rotated by final first-moment phase; this is a shape comparison, not lineage.", fontsize=9.5)
            fig.text(.045, .902, "Coloured support bands show the prescribed initial footprint (I/S/O at the 10% edge). Full-ring coordinate is physical arc on the initial basal, also for Spring.", fontsize=9)
            for ax, r in zip(axes.flat, part.to_dict("records")):
                cid = int(r["condition_id"]); C = float(r["circumference"])
                q = initial["q"]; order = np.argsort((q-r["preparation_axis_fraction"]+.5) % 1-.5)
                x = C*((q-r["preparation_axis_fraction"]+.5) % 1-.5)
                ax.plot(x[order], initial["u"][order], color="#999", lw=1.1, ls=":")
                for shift in (-C, 0, C):
                    a = max(-C/2, r["physical_distance_at_onset"]-r["full_physical_support"]/2+shift)
                    b = min(C/2, r["physical_distance_at_onset"]+r["full_physical_support"]/2+shift)
                    if b > a:
                        ax.axvspan(a, b, color=COLORS[r["footprint_relative_to_core_10"]], alpha=.12)
                line = f"d₀={r['physical_distance_at_onset']:.5g}; θ*={r['target_angle_in_P_coordinates_deg']:.3g}°; {REGIMES[r['footprint_relative_to_core_10']]}"
                if is_complete(r):
                    p = by_id[cid]
                    px = C*((p.material_position.to_numpy()-r["preparation_axis_fraction"]+.5) % 1-.5); idx = np.argsort(px)
                    ax.plot(px[idx], p.final_active.to_numpy()[idx], color="#16252d", lw=1.3)
                    ax.plot(px[idx], p.rotated_basal_active.to_numpy()[idx], color="#517fa8", lw=.8, alpha=.8)
                    line += f"\nR={r['R_final']:.3g}; peaks={r['final_resolved_peaks']:.0f}; regions={r['final_halfcontrast_regions']:.0f}"
                    tail = "tail PASS" if r["full_field_tail_status"] == "passes_finite_tail_gate" else "tail unavailable" if "unavailable" in str(r["full_field_tail_status"]) else "tail FAIL"
                    line += f"; {tail}"
                else:
                    ax.text(.5, .5, "Native trajectory pending", transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#777")
                ax.set_title(line, fontsize=8.5, loc="left")
                ax.set_ylim(ymin-pad, ymax+pad); ax.set_xlim(-C/2, C/2)
                ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=8)
            for ax in axes.flat[len(part):]:
                ax.axis("off")
            for ax in axes[-1, :]:
                if ax.axison: ax.set_xlabel("Initial physical arc from preparation axis", fontsize=8)
            for ax in axes[:, 0]: ax.set_ylabel("Native active concentration", fontsize=9)
            footer(fig, status, "Endpoint is 13,000 units after withdrawal; final orientation alone does not establish stationarity. All profiles use common native-amplitude limits within this model.")
            pdf.savefig(fig); plt.close(fig)


def html_report(table, status, source, artifacts, stamp):
    metrics = [m for _, _, items in MATRIX_GROUPS for m, _, _, _ in items]
    metrics += ["angle_withdrawal_deg", "angle_final_deg", "M0_ratio_final", "R_sampled_minimum", "physical_C_withdrawal", "physical_C_final", "fixed_region_excess_withdrawal", "fixed_region_rebound_from_sampled_minimum"]
    cols = ["condition_id", "model", "size_parameter", "cue_duration", "physical_distance_at_onset", "target_angle_in_P_coordinates_deg", "status", "execution", "centre_relative_to_core_10", "footprint_relative_to_core_10", "centre_relative_to_core_50", "footprint_relative_to_core_50", "fraction_cue_support_inside_core_10", "full_field_tail_status", "actual_cue_path", "actual_recovery_path"] + list(dict.fromkeys(metrics))
    embedded = clean(table[cols].to_dict("records"))
    cards = "".join(f'<article><h3>{name}</h3><p><a href="supporting/vectors/{model}_core_supplement.pdf">Readable core supplement: 5 pages</a> · <a href="supporting/vectors/{model}_native_profile_atlas.pdf">Full native-profile appendix: 15 pages</a> · <a href="supporting/vectors/{model}_complete_supplement.pdf">Combined archive PDF</a></p><a href="{model}_all_vectors.png"><img loading="lazy" src="{model}_all_vectors.png" alt="{name}: every requested withdrawal and retained vector"></a><p><a href="{model}_cue_response.png">Cue-on response / coverage matrices</a> · <a href="{model}_amount_region.png">Amount and original-region matrices</a> · <a href="{model}_retained_organization.png">Finite endpoint organization matrices</a> · <a href="{model}_all_longest_cue.png">All longest-cue histories</a></p></article>' for model, name in NAMES.items())
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 3 — Complete narrow-grid stationary-cue supplement</title><style>
body{font:16px/1.55 system-ui,sans-serif;color:#20313b;background:#f5f6f7;margin:0}main{max-width:1250px;margin:auto;padding:28px}h1,h2,h3{line-height:1.2}a{color:#176380}section,article{background:white;border:1px solid #d7dfe3;border-radius:6px;padding:22px;margin:22px 0}img{width:100%;height:auto}p{max-width:1050px}.notice{background:#fff4db;border-left:5px solid #b2761e;padding:16px}label{display:inline-block;margin:8px 18px 8px 0}select,input{font:inherit;padding:5px}canvas{width:100%;height:auto;background:#fff}.scroll{overflow:auto;max-height:650px}table{border-collapse:collapse;font-size:13px;white-space:nowrap}th,td{padding:6px 9px;border-bottom:1px solid #dbe2e5;text-align:right}th{position:sticky;top:0;background:#eaf0f3}code{font-size:13px;overflow-wrap:anywhere}.equation{text-align:center;font-family:serif;font-size:22px;margin:24px}.small{font-size:13px;color:#536773}footer{font-size:13px}</style><main>
<h1>Complete stationary-cue experiment supplement</h1><p>Figure 3 · narrow native size range, denser common initial physical arc grid.</p>
<p><b>Generated STAMP.</b> Authoritative collector: DATASTAMP. <a href="../../arc_review/index.html">Scientific report and compact examples</a> · <a href="../narrow_initial_design.png">All initial profiles and cue footprints</a>.</p>
<div class="notice"><b>Coverage: COMPLETE / 645 requested conditions.</b> This is an explicitly partial snapshot until every native branch has completed and terminal checks pass. Pending cells are retained. The design contains 175 proved exact reuses and 470 fresh branches; all initial overlaps remain included. Cue-on reorientation is the main experimental readout. Cue-free endpoints provide separate finite-horizon context.</div>
<section><h2>What the full supplement contains</h2><p>Every requested distance, duration and basal appears in the vector grids and numerical matrix selector below. Four model PDFs add cue-on response/coverage matrices, total activity and original-region readouts, finite endpoint organization, all longest-cue time courses, and one native final-profile panel for every requested condition. Pending profiles show the actual initial basal and a pending label. These are comprehensive supporting views, not a proposed main-figure catalogue.</p>
<p><a href="../supporting/data/narrow_arc_summary.csv">Full 166-column authoritative condition table</a> · <a href="../supporting/data/narrow_arc_initial_context.csv">All 129 placement contexts</a> · <a href="../supporting/data/narrow_arc_traces.npz">All completed physical vector/amount/regional histories</a> · <a href="../supporting/data/narrow_arc_cue_on_sampled_secants.csv">Every saved cue-on rate interval</a> · <a href="../supporting/data/narrow_arc_cue_on_RHS_rates.csv">Selected exact-state native RHS rates</a> · <a href="../supporting/data/narrow_arc_crossings.csv">All sampled regional-entry brackets and rebound flags</a> · <a href="../supporting/data/narrow_arc_final_profiles.csv">All final native active profiles</a>.</p>
<p><a href="supporting/data/displayed_conditions.csv">Exact plotted condition values</a> · <a href="supporting/metadata/manifest.json">Renderer, source and artifact hashes</a> · <a href="README.md">Methods and scope</a> · <a href="LEGENDS.md">Supporting figure legends</a> · <a href="supporting/data/initial_basal_profiles.csv">Compact initial basal profiles</a>.</p></section>
<section><h2>Every condition: choose the quantity</h2><label>Model <select id="model"></select></label><label>Quantity <select id="metric"></select></label><label>Initial footprint <select id="regime"><option value="all">All</option><option>inside</option><option>straddling</option><option>outside</option></select></label><label>Operational edge <select id="edge"><option value="10">10% contrast</option><option value="50">50% contrast</option></select></label><canvas id="matrix" width="1200" height="680"></canvas><p id="matrixnote" class="small"></p><p>Grey: pending, unavailable quantity, or excluded by the selected footprint filter. Exact values below remain unclipped. Duration columns are tested categories, not uniformly spaced elapsed time. R is a moment-magnitude ratio; it can exceed one.</p></section>
<section><h2>Searchable full-grid values and actual source cases</h2><label>Filter model / size / status / condition <input id="search" type="search" placeholder="e.g. holmes_model3 or pending"></label><p id="tablecount"></p><div class="scroll"><table><thead id="head"></thead><tbody id="body"></tbody></table></div><p class="small">Display JSON is rounded to eight significant digits. The linked CSV/NPZ files retain authoritative floating-point precision. Source paths identify the actual duration-specific cue and recovery archives, including reused native histories.</p></section>
CARDS
<section><h2>Definitions and interpretation limits</h2><div class="equation">M<sub>n</sub>(t) = C(t)〈u(q,t)e<sup>2πinq</sup>〉, &nbsp; P(t) = M<sub>1</sub>(t)/M<sub>1</sub>(0), &nbsp; R(t) = |P(t)|.</div><p>The raw physical moment uses fixed ring reference axes. Complex division by the initial first moment rotates the normalized P axes once, so P(0)=1. Activity-vector phase is a global concentration-weighted normal direction. It is not cell velocity, peak travel, or molecular/site ancestry.</p>
<div class="equation">Error removed = 1 − |wrap(φ − θ*)| / |wrap(−θ*)|.</div><p>The actual target θ* includes the measured initial first-moment offset from the declared preparation axis. Phase-dependent values are unavailable when R&lt;0.1. Negative error-removal values are retained in the table; only figure background colours are clipped to 0–1. This endpoint measure is not trajectory speed or a binary success criterion.</p>
<p>Sampled angular secants use their actual left/right saved times; weak endpoint R creates gaps and no phase is unwrapped across them. Exact-antipodal angular rates are withheld. Near-antipodal vectors require their magnitude and profiles. Selected native RHS rates are available only where the actual full state is retained; inactive fields are never inferred for legacy primary-only histories. Initial rate differences can arise from the angular projection of the finite-width input and the initial moment magnitude, without demonstrating incumbent remodeling.</p>
<p>The initial 50%/10% edges are linearly interpolated min-plus-contrast thresholds. Cue-centre membership and whole-footprint intersection are separate. I/S/O are spatial readout classes, not demonstrated internal causal pathways. Spring uses the qualified native mean-force global-extension closure: the material cue centre stays fixed while its physical arc changes with C(t), and its full physical support stays 0.2. Other native supports are 3; native clocks are not calibrated across models.</p>
<p>Original-region excess integrates the fixed initial material region and subtracts the transported initial background amount. Its decrease can reflect redistribution; entry times are sampled brackets with right censoring and rebound context. They do not prove site death or provide refined minimum effective pulse durations. Final profiles are measured after 13,000 cue-free units; a passing finite-tail gate and an unavailable full-field tail are distinct from an asymptotic claim.</p>
<p><a href="../../standard_arclength_sizes/README.md">Earlier complete broad-size fixed-arc arm</a> · <a href="../../standard_angles_sizes/README.md">Earlier complete fixed-angle arm</a> · <a href="../../paired_stationary_review/README.md">Preserved paired-arm vector review</a>. These preceding designs remain distinct; none is silently relabelled as fresh narrow-grid data.</p></section>
<footer>No new native integrations or model parameters are introduced by this supplement. Opening or manually refreshing loads the latest generated snapshot; there is no timed self-refresh.</footer>
<script id="data" type="application/json">DATA</script><script>
const rows=JSON.parse(document.getElementById('data').textContent), names={wave_pinning:'Wave-pinning',goryachev:'Goryachev',holmes_model3:'Holmes',spring:'Spring'};
const metrics=METRICS, el=id=>document.getElementById(id), fmt=v=>v===null||v===undefined?'—':typeof v==='number'?Number(v.toPrecision(6)).toString():String(v);
Object.entries(names).forEach(([v,n])=>{const o=document.createElement('option');o.value=v;o.textContent=n;el('model').appendChild(o)});
metrics.forEach(v=>{const o=document.createElement('option');o.value=v;o.textContent=v.replaceAll('_',' ');el('metric').appendChild(o)});
function draw(){const all=rows.filter(r=>r.model===el('model').value), metric=el('metric').value, edge=el('edge').value, regime=el('regime').value, sizes=[...new Set(all.map(r=>r.size_parameter))].sort((a,b)=>a-b), ds=[...new Set(all.map(r=>r.physical_distance_at_onset))].sort((a,b)=>a-b), ts=[...new Set(all.map(r=>r.cue_duration))].sort((a,b)=>a-b), canvas=el('matrix'), c=canvas.getContext('2d'); c.clearRect(0,0,1200,680);const vals=all.filter(r=>r.status==='complete'&&r[metric]!==null&&typeof r[metric]==='number').map(r=>r[metric]);let lo=vals.length?Math.min(...vals):0,hi=vals.length?Math.max(...vals):1;if(metric.includes('error_closed')||metric==='secant_coverage'){lo=0;hi=1}if(hi===lo)hi=lo+1;const plotw=300,ploth=480,cellw=plotw/ts.length,cellh=ploth/ds.length;c.font='15px system-ui';c.fillStyle='#22343d';c.fillText(metric.replaceAll('_',' '),22,26);for(let s=0;s<sizes.length;s++){const x0=80+s*395,y0=92;c.font='bold 16px system-ui';c.fillText((all[0].model==='spring'?'C₀=':'C=')+sizes[s],x0,y0-27);for(let j=0;j<ds.length;j++){const context=all.find(r=>r.size_parameter===sizes[s]&&r.physical_distance_at_onset===ds[j]);c.font='12px system-ui';c.fillStyle='#22343d';c.fillText(fmt(ds[j])+' '+context['footprint_relative_to_core_'+edge][0].toUpperCase(),x0-68,y0+(j+.62)*cellh);for(let k=0;k<ts.length;k++){const r=all.find(r=>r.size_parameter===sizes[s]&&r.physical_distance_at_onset===ds[j]&&r.cue_duration===ts[k]),v=r[metric],shown=regime==='all'||r['footprint_relative_to_core_'+edge]===regime,valid=shown&&r.status==='complete'&&v!==null&&typeof v==='number';let p=valid?Math.min(1,Math.max(0,(v-lo)/(hi-lo))):0;c.fillStyle=valid?'hsl(200,55%,'+(96-57*p)+'%)':'#e9ecee';c.fillRect(x0+k*cellw,y0+j*cellh,cellw-1,cellh-1);c.font='10px system-ui';c.fillStyle=p>.65&&valid?'white':'#233741';c.textAlign='center';c.fillText(valid?fmt(v):r.status==='complete'?'—':'pending',x0+(k+.5)*cellw,y0+(j+.6)*cellh);c.textAlign='left'}}c.fillStyle='#22343d';c.font='12px system-ui';for(let k=0;k<ts.length;k++)c.fillText(fmt(ts[k]),x0+(k+.35)*cellw,y0+ploth+22);c.fillText('Cue duration (categorical native units)',x0,y0+ploth+47)}el('matrixnote').textContent='Colour bounds '+fmt(lo)+' to '+fmt(hi)+'. '+all.filter(r=>r.status==='complete').length+'/'+all.length+' complete in this model. I/S/O applies to the selected initial operational edge.';drawTable()}
const headers=['condition_id','model','size_parameter','physical_distance_at_onset','cue_duration','status','footprint_relative_to_core_10','angle_withdrawal_deg','angle_final_deg','R_withdrawal','R_final','withdrawal_error_closed','retained_error_closed','full_field_tail_status'];
function drawTable(){const q=el('search').value.toLowerCase(),selected=rows.filter(r=>r.model===el('model').value).filter(r=>!q||JSON.stringify(r).toLowerCase().includes(q));el('head').innerHTML='<tr>'+headers.map(h=>'<th>'+h+'</th>').join('')+'</tr>';el('body').innerHTML=selected.map(r=>'<tr>'+headers.map(h=>'<td>'+fmt(r[h])+'</td>').join('')+'</tr>').join('');el('tablecount').textContent=selected.length+' condition rows shown. Actual source paths and all derived quantities are in the full CSV.'}
['model','metric','regime','edge'].forEach(id=>el(id).addEventListener('change',draw));el('search').addEventListener('input',drawTable);draw();
</script></main></html>'''.replace("STAMP", stamp, 1).replace("DATASTAMP", status["generated_at"]).replace("COMPLETE", str(status["complete"]), 1).replace("CARDS", cards).replace("METRICS", json.dumps(metrics)).replace("DATA</script>", json.dumps(embedded, separators=(",", ":"), allow_nan=False)+"</script>")

    if status["complete"] == 645:
        page = page.replace("This is an explicitly partial snapshot until every native branch has completed and terminal checks pass. Pending cells are retained.", "All 645 native outcomes are present in this readout snapshot. Native terminal validation and independent checks remain separate evidence receipts; the scientific report records their current state.")
    return page


def main():
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42, "svg.fonttype": "none"})
    (OUT/"supporting/data").mkdir(parents=True, exist_ok=True)
    (OUT/"supporting/metadata").mkdir(parents=True, exist_ok=True)
    (OUT/"supporting/vectors").mkdir(parents=True, exist_ok=True)
    status = json.loads((TOPIC/"supporting/metadata/latest_status.json").read_text())
    data = TOPIC/"supporting/data"
    filenames = ["narrow_arc_summary.csv", "narrow_arc_final_profiles.csv", "narrow_arc_cue_on_RHS_rates.csv"]
    sources = {str(data/name): sha(data/name) for name in filenames}
    for name in filenames:
        assert sources[str(data/name)] == status["artifact_sha256"]["supporting/data/"+name]
    table = pd.read_csv(data/filenames[0]).copy(); profiles = pd.read_csv(data/filenames[1]); rhs = pd.read_csv(data/filenames[2])
    interactive = data/"narrow_cue_on129_interactive.json"
    presentation_path = TOPIC/"supporting/metadata/presentation.json"
    presentation = json.loads(presentation_path.read_text())
    assert presentation["interactive_sha256"] == sha(interactive)
    assert recorded_source_digest(presentation["source_sha256"], data/filenames[0]) == sources[str(data/filenames[0])]
    assert presentation["collector_generation"] == status["generated_at"]
    sources[str(interactive)] = sha(interactive)
    sources[str(presentation_path)] = sha(presentation_path)
    assert len(table) == 645 and table.condition_id.nunique() == 645
    table["secant_coverage"] = table.cue_secant_available_intervals / table.cue_secant_intervals
    selected_rhs = rhs.merge(table[["condition_id", "cue_duration"]], on="condition_id")
    rates = selected_rhs[np.isclose(selected_rhs.requested_time, selected_rhs.cue_duration)].set_index("condition_id")
    assert rates.index.is_unique
    table["rhs_withdrawal"] = table.condition_id.map(rates.native_RHS_angular_rate_deg_per_native_time)
    atomic(OUT/"supporting/data/displayed_conditions.csv", table.to_csv(index=False))
    basal = {}; basal_path = OUT/"supporting/data/initial_basal_profiles.csv"
    lineage_path = OUT/"supporting/metadata/initial_basal_lineage.json"
    if not basal_path.exists():
        rows = []; lineage = []
        for case in table.drop_duplicates(["model", "size_parameter"]).to_dict("records"):
            path = RUN/case["input"]
            assert sha(path) == case["input_sha256"]
            with np.load(path) as z:
                n = case["cells"]
                u = z["state"][:n] if case["model"] == "spring" else z["state"][0]
                q = (np.arange(n)+.5)/n if case["model"] == "spring" else z["x"]/case["circumference"]
                rows.extend(dict(model=case["model"], size_parameter=case["size_parameter"], material_position=float(a), initial_active=float(b)) for a,b in zip(q,u))
            lineage.append({key:case[key] for key in ("model","size_parameter","circumference","cells","input","input_sha256","preparation_axis_fraction")})
        atomic(basal_path, pd.DataFrame(rows).to_csv(index=False))
        atomic(lineage_path, json.dumps(dict(source_manifest_path=str(RUN/"manifest.json"), source_manifest_sha256=sha(RUN/"manifest.json"), basals=lineage), indent=2)+"\n")
    for key, part in pd.read_csv(basal_path).groupby(["model", "size_parameter"]):
        expected = table[(table.model==key[0])&(table.size_parameter==key[1])].iloc[0]
        assert len(part) == expected.cells
        basal[key] = dict(q=part.material_position.to_numpy(), u=part.initial_active.to_numpy())
    sources[str(basal_path)] = sha(basal_path)
    sources[str(lineage_path)] = sha(lineage_path)
    # All longest-cue values are already exported at full precision by the
    # guarded presenter. Using this compact JSON avoids a plotting dependency
    # on ignored large NPZ archives in a compact source checkout.
    traces = {}
    for curve in json.loads(interactive.read_text()):
        if curve["status"] != "complete":
            continue
        prefix = f"c{int(curve['condition_id'])}_"
        traces[prefix+"phase"] = np.zeros(len(curve["time"]), dtype=int)
        for old, new in (("time","time"),("heading","angle_masked_deg"),("R","R"),("M0_ratio","M0_ratio"),("regional_excess","fixed_region_transported_background_excess"),("rate_left","cue_rate_left"),("rate_right","cue_rate_right"),("rate","cue_sampled_angular_secant")):
            traces[prefix+new] = np.array([np.nan if v is None else v for v in curve[old]], dtype=float)
    stamp = datetime.now(timezone.utc).isoformat()
    artist_checks = {}
    artifacts = []
    if traces or status["complete"] == 0:
        for model in NAMES:
            path = OUT/"supporting/vectors"/(model+"_complete_supplement.pdf")
            with PdfPages(path) as pdf:
                fig, checks = vector_page(model, table, status)
                artist_checks[model] = checks
                png = OUT/(model+"_all_vectors.png"); fig.savefig(png, dpi=150, facecolor="white"); pdf.savefig(fig); plt.close(fig); artifacts.append(png)
                for group in MATRIX_GROUPS:
                    fig = matrix_page(model, table, status, group)
                    png = OUT/(model+"_"+group[0]+".png"); fig.savefig(png, dpi=135, facecolor="white"); pdf.savefig(fig); plt.close(fig); artifacts.append(png)
                fig = trajectory_page(model, table, traces, status)
                png = OUT/(model+"_all_longest_cue.png"); fig.savefig(png, dpi=135, facecolor="white"); pdf.savefig(fig); plt.close(fig); artifacts.append(png)
                profile_pages(model, table, profiles, basal, pdf, status)
            artifacts.append(path)
            reader = PdfReader(path)
            assert len(reader.pages) == 20
            for label, pages in (("core_supplement", range(5)), ("native_profile_atlas", range(5, 20))):
                split = OUT/"supporting/vectors"/(model+"_"+label+".pdf")
                writer = PdfWriter()
                for page_index in pages:
                    writer.add_page(reader.pages[page_index])
                with split.open("wb") as stream:
                    writer.write(stream)
                artifacts.append(split)
            print(str(OUT/(model+"_all_vectors.png")), flush=True)
    htmlpath = OUT/"index.html"
    atomic(htmlpath, html_report(table, status, sources, artifacts, stamp)); artifacts.append(htmlpath)
    text = f"""# Complete-condition narrow-grid Figure 3 supplement

Generated {stamp}; authoritative collector {status['generated_at']}. Coverage {status['complete']}/645; all pending cells retained. Native producer state: {status['producer']}.

[HTML supplement](index.html) · [current scientific report](../../arc_review/index.html) · [all placement contexts](../supporting/data/narrow_arc_initial_context.csv) · [exact plotted values](supporting/data/displayed_conditions.csv).

Each model has a readable five-page core supplement: one complete withdrawal/retained vector grid, three complete condition-matrix pages and all longest-cue time courses. The separate fifteen-page native-profile appendix covers every size × duration × placement. A combined twenty-page archival PDF is also retained. These full appendices are linked rather than forced into the main manuscript packet. Every requested placement is included. Pending profile panels show only the exact initial basal. Native amplitude limits and each trajectory quantity's y limits are common across sizes within a model. I/S/O refers to whole initial cue footprints at the 10% min-plus-contrast edge; centre membership remains separate in the table. Both 50%/10% classes are available in the HTML selector. Spatial classes do not assign a causal pathway.

Physical moments use instantaneous C times native-cell means. P=M1/M1(0) rotates normalized axes once; R=|P| can exceed1. Arrows use one .30-data-unit-per-P scale with equal x/y coordinates, so actual directions and lengths are preserved; R<.1 is a magnitude-only circle. Final background is branch-safe fraction of initial target-angle error removed, clipped only for colour; raw values are exported. Failed and unavailable full-field tails are independent flags. Native endpoint is13000 cue-free units after withdrawal.

Cue-on curves use every saved time for T1000. Rate segments span their actual saved intervals and never bridge weak endpoint R; no phase unwrapping. Selected native RHS withdrawal rates use exact full retained states/snapshots only. No inactive field is inferred for primary-only histories; rate availability and all endpoint R remain exported. Source times and cadence are in the authoritative rate tables. An observed angular rate is not local-profile travel, nor a uniquely identified incumbent-remodeling mechanism. Initial angular projection and initial moment magnitude matter.

Original-region excess retains the fixed initial material region and transported initial background amount. Depletion can be redistribution; first/retained sampled entries, censoring and rebound are in the full regional table. No minimum effective pulse duration or site ancestry is inferred. Rotated basal overlays are shape comparisons, and finite tail gates do not establish an asymptotic attractor.

Each PDF is complete with respect to the requested grid, not complete with respect to still-pending native outputs. There is no timed HTML refresh. The renderer adds no integrations and does not modify the native campaign or guarded terminal pipeline. Display JSON rounds to8 significant digits; CSV/NPZ sources retain full precision. The supplement redraws from compact CSVs, full-precision longest-cue JSON and the exported initial-basal CSV; ignored native NPZ files are not needed to rebuild the plots. Native source lineage and complete raw histories remain separately retained.
"""
    atomic(OUT/"README.md", text); artifacts.append(OUT/"README.md")
    # Ensure the atomic source generation did not change underneath rendering.
    for path, digest in sources.items():
        assert sha(path) == digest, "Source generation changed: "+path
    result = dict(generated_at=stamp, collector_generation=status["generated_at"], requested=645, complete=int(status["complete"]),
                  producer=status["producer"], pending=645-int(status["complete"]), total_profiles=len(profiles),
                  complete_for_requested_grid=True, all_native_outputs_complete=bool(status["complete"] == 645), new_integrations=0,
                  renderer=str(Path(__file__).resolve()), renderer_sha256=sha(__file__), source_sha256=sources,
                  compact_rebuild_inputs="Authoritative CSVs, full-precision longest-cue JSON and exported initial-basal CSV; ignored native archives are not required to redraw this supplement.",
                  vector_artist_checks=artist_checks, pdf_pages_per_model=20,
                  display_precision="8 significant digits in embedded JSON only; exact table and authoritative sources linked",
                  artifact_sha256={str(p.relative_to(OUT)): sha(p) for p in artifacts})
    result["artifact_sha256"]["supporting/data/displayed_conditions.csv"] = sha(OUT/"supporting/data/displayed_conditions.csv")
    atomic(OUT/"supporting/metadata/manifest.json", json.dumps(result, indent=2)+"\n")
    print(json.dumps({"generated_at": stamp, "complete": int(status["complete"]), "html": str(htmlpath)}, indent=2))


if __name__ == "__main__":
    main()
