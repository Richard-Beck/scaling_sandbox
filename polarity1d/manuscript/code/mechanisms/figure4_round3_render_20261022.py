"""Portable Figure 4 round-3 composition from compact records.

No simulation or native archive is opened by this renderer. Undefined estimates
remain absent, and the complete omission inventory is exported with the figure.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
OUT = ROOT / "figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4"
MODELS = ("wave_pinning", "holmes_model3", "spring_ring", "goryachev")
NAMES = dict(wave_pinning="Wave-pinning", holmes_model3="Holmes", spring_ring="Spring", goryachev="Goryachev")
COLORS = dict(wave_pinning="#17619c", holmes_model3="#28783b", spring_ring="#aa7a00", goryachev="#b93434")
NULLS = ("drift", "same_rate", "perimeter")
NULL_NAMES = dict(drift="Drifting peak\n(same physical $D_s$)", same_rate="Jumping peak\n(same rate)", perimeter="Jumping peak\n(rate $\\propto$ perimeter)")
NULL_B_NAMES = dict(drift="Drift only", same_rate="Drift + same-rate jumps", perimeter="Drift + perimeter-rate jumps")
SIZES = (20., 30., 40., 50., 60.)
STYLE = dict(fontsize=6.5, titlesize=7., labelsize=6.7, markersize=2., linewidth=.72)


def number(value, default=math.nan):
    try:
        return float(value) if value not in (None, "", "None") else default
    except (TypeError, ValueError):
        return default


def finite(value):
    return math.isfinite(number(value))


def yes(value):
    return value is True or str(value).lower() in ("true", "1")


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=True) + "\n")


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def first(row, *keys, default=None):
    return next((row[key] for key in keys if key in row and row[key] not in (None, "")), default)


def n(row, *keys, default=math.nan):
    return number(first(row, *keys), default)


def null_kind(row):
    raw = str(first(row, "null", "null_type", "kind", "model", default=""))
    if "perimeter" in raw:
        return "perimeter"
    if "same" in raw or raw in ("jump", "jumps", "drift_jump", "constant_jump", "jump_only"):
        return "same_rate"
    if "drift" in raw:
        return "drift"
    return raw


def primary(rows):
    """Keep the declared main arm/view without suppressing forced conditions."""
    return [row for row in rows if row.get("arm", "common") in ("common","",None) and row.get("view", "main") in ("main","",None)]


def equivalent_size(row):
    value = n(row, "C_label", "C", "circumference", "size")
    return 20 * value if row.get("model") == "spring_ring" and value <= 3.1 else value


def size_color(model, value):
    index = int(np.argmin(np.abs(np.asarray(SIZES) - value)))
    if model in NULLS or model == "null":
        return str(np.linspace(.64, .03, 5)[index])
    rgb = np.asarray(to_rgb(COLORS[model]))
    return tuple(1 - np.linspace(.38, 1., 5)[index] * (1 - rgb))


def style_axis(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(.48)
    ax.spines["bottom"].set_linewidth(.48)
    ax.tick_params(axis="both", which="major", width=.45, length=2., pad=1.6)
    ax.tick_params(axis="both", which="minor", width=.35, length=1.1)
    ax.set_axisbelow(True)


def scientific_tick(value, pos=None):
    if value == 0:
        return "0"
    exponent = int(round(math.log10(value))) if value > 0 else 0
    if value > 0 and math.isclose(value, 10 ** exponent, rel_tol=1e-6):
        return rf"$10^{{{exponent}}}$"
    return f"{value:g}"


def log_f(ax):
    ax.set_xscale("log")
    ax.set_xlim(7e-5, 1.4)
    ax.set_xticks([1e-4, 1e-2, 1], [r"$10^{-4}$", r"$10^{-2}$", "1"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"Input fraction $f$", labelpad=1.5)


def facet(ax, text):
    ax.set_title(text, color="black", fontsize=STYLE["titlesize"], pad=3.5, fontweight="normal")


def panel_label(fig, letter, y):
    fig.text(.017, y, letter, fontsize=10, fontweight="bold", ha="left", va="top")


def ramp_key(fig, models, y, *, x=.74, width=.225, labels=True):
    """One local five-step family key directly above the plots it serves."""
    step = width / len(models)
    for j, model in enumerate(models):
        ax = fig.add_axes([x + j * step, y, step * .83, .006])
        colors = [size_color(model, value) for value in SIZES]
        ax.imshow(np.arange(5)[None], aspect="auto", cmap=LinearSegmentedColormap.from_list(model, colors, N=5), extent=(0, 5, 0, 1), interpolation="nearest")
        ax.set_axis_off()
        if labels:
            for k,value in enumerate((1,1.5,2,2.5,3) if model=="spring_ring" else SIZES):
                fig.text(x + j*step + (k+.5)*step*.83/5,y+.008,f"{value:g}",ha="center",va="bottom",fontsize=5.1)


def bounded_point(ax, x, y, lo, hi, color, *, marker="o", open_marker=False, ms=2.1, lw=.46):
    """Intervals clipped to limits use outward edge triangles, never hidden tails."""
    ymin, ymax = ax.get_ylim()
    point_below = finite(y) and y < ymin
    point_above = finite(y) and y > ymax
    if finite(y):
        if ymin <= y <= ymax:
            ax.plot(x, y, marker, ms=ms, color=color, mec=color, mew=.55, mfc="white" if open_marker else color, zorder=4,clip_on=False)
        else:
            ax.plot(x, ymin if y < ymin else ymax, "v" if y < ymin else "^", color=color, mfc="white" if open_marker else color, ms=3.0, mew=.6, clip_on=False, zorder=5)
    if not math.isnan(lo) and not math.isnan(hi):
        bottom, top = max(lo, ymin), min(hi, ymax)
        if top >= bottom:
            ax.plot([x, x], [bottom, top], color=color, lw=lw, alpha=.8, zorder=2)
        # The point triangle already identifies its occupied edge. An open CI
        # triangle there would cover its fill and falsely imply a forced point.
        if lo < ymin and not point_below:
            ax.plot(x, ymin, "v", color=color, mfc="white", ms=2.6, mew=.55, clip_on=False, zorder=5)
        if hi > ymax and not point_above:
            ax.plot(x, ymax, "^", color=color, mfc="white", ms=2.6, mew=.55, clip_on=False, zorder=5)


def rows_for(rows, **kwargs):
    return [row for row in rows if all(str(row.get(key)) == str(value) for key, value in kwargs.items())]


def load_kymographs(data, prefix):
    """Index CSV contract: file, model/null, C_label, f; arrays time/field/C."""
    entries = read_csv(data / (prefix + "_kymographs.csv"))
    out = []
    for entry in entries:
        path = Path(first(entry, "file", "path", "npz", default=""))
        if not path.is_absolute():
            path = data / path
        if not path.is_file():
            continue
        with np.load(path, allow_pickle=False) as archive:
            t = archive["time"].astype(float)
            field = archive["field" if "field" in archive else "primary"].astype(float)
            circumference = archive["circumference"] if "circumference" in archive else n(entry, "C_mean", "C_label")
            space_edges=archive['space_edges'].astype(float) if 'space_edges' in archive else np.linspace(0,float(np.mean(circumference)),field.shape[1]+1)
            out.append(dict(entry, time=t, field=field, circ=float(np.mean(circumference)), space_edges=space_edges,npz_path=path))
    return out


def kymo(ax, record, vmax, *, ylabel=False, label="", switches=(),elapsed=False):
    field, t, circ = record["field"], record["time"], record["circ"]
    origin=float(t[0]) if elapsed else 0.
    t=t-origin
    # Pixel centres retain native sample times; time always increases upward.
    middle=(t[:-1]+t[1:])/2
    time_edges=np.r_[t[0],middle,t[-1]] if len(t)>1 else np.array([t[0]-.5,t[0]+.5])
    ax.pcolormesh(record['space_edges'],time_edges/1000,field,cmap="Greys",vmin=0,vmax=vmax,rasterized=True,shading='flat',linewidth=0)
    ax.set_ylim(time_edges[0]/1000,time_edges[-1]/1000)
    ax.set_xlim(0, circ)
    ax.set_xticks([])
    ax.set_xlabel(label, fontsize=6., labelpad=1.4)
    if ylabel:
        ax.set_ylabel(r"Time ($\times 10^3$)", fontsize=6.2, labelpad=2)
    else:
        ax.tick_params(labelleft=False)
    for when in switches:
        ax.plot([-.10 * circ, -.025 * circ], [(when-origin) / 1000] * 2, color="black", lw=.65, clip_on=False)
    style_axis(ax)


def correlation(ax, rows, model, *, annotation=None, null=False, legend=False):
    relevant = [row for row in rows if (null_kind(row) == model if null else row.get("model") == model)]
    sizes = sorted({equivalent_size(row) for row in relevant if finite(equivalent_size(row))})
    if len(sizes) > 2:
        sizes = [sizes[0], sizes[-1]]
    for j, size in enumerate(sizes):
        selected = sorted([row for row in relevant if math.isclose(equivalent_size(row), size)], key=lambda row: n(row, "lag"))
        x = np.array([n(row, "lag") for row in selected]) / 1000
        y = np.array([n(row, "C", "correlation") for row in selected])
        color = size_color("null" if null else model, size)
        ax.plot(x, y, "-" if j == 0 else ":", lw=.8, color=color, label=f"C = {size:g}")
    ax.set_ylim(-.08, 1.05)
    if null:ax.set_xlim(0,24)
    ax.set_xlabel(r"Lag ($\times 10^3$)", labelpad=1.4)
    if annotation:
        ax.text(.98, .94, annotation, transform=ax.transAxes, ha="right", va="top", fontsize=6.3)
    if legend and relevant:
        ax.legend(frameon=False, fontsize=5.9, loc="lower left", handlelength=1.7, borderpad=.2, labelspacing=.1)
    style_axis(ax)


def null_panels(fig, grids, data, exported):
    examples = load_kymographs(data, "null_A")
    acf = read_csv(data / "null_A_acf.csv")
    for j, kind in enumerate(NULLS):
        outer = grids[j].subgridspec(2, 1, height_ratios=[1.45, 1.], hspace=.32)
        pair = outer[0].subgridspec(1, 2, width_ratios=[1, 3], wspace=.28)
        selected = [row for row in examples if null_kind(row) == kind]
        vmax = max((row["field"].max() for row in selected), default=1.)
        for k, size in enumerate((20, 60)):
            ax = fig.add_subplot(pair[k])
            record = next((row for row in selected if math.isclose(equivalent_size(row), size)), None)
            if record is not None:
                kymo(ax, record, vmax, ylabel=k == 0, label=f"C{size}")
            else:
                ax.text(.5, .5, "Awaiting null record", transform=ax.transAxes, ha="center", fontsize=6)
                ax.set_axis_off()
        fig.text((grids[j].get_position(fig).x0 + grids[j].get_position(fig).x1) / 2, grids[j].get_position(fig).y1 + .007, NULL_NAMES[kind], ha="center", fontsize=7)
        ax = fig.add_subplot(outer[1])
        annotation = rf"$\alpha = {dict(drift='2',same_rate='0',perimeter='-1')[kind]}$"
        correlation(ax, acf, kind, annotation=annotation, null=True, legend=j == 0)
        if j == 0:
            ax.set_ylabel(r"Heading $C(\ell)$", labelpad=2)
        if kind == "same_rate":
            ax.text(.98, .80, "(curves coincide)", transform=ax.transAxes, ha="right", fontsize=6.)
    exported["A_acf"] = acf


def alpha_panel(ax, rows, model, *, null=False, first_axis=False, ylabel=False):
    log_f(ax)
    ax.set_ylim((-1.5, 2.5) if null else (-1.5, 4.))
    for guide in (2, 0, -1):
        ax.axhline(guide, color=".79", lw=.4, ls="--", zorder=0)
        if first_axis and null:
            ax.text(.98, guide, {2:"drift",0:"same-rate jumps",-1:"perimeter-rate jumps"}[guide], transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=5.3)
    relevant = [row for row in rows if (null_kind(row) == model if null else row.get("model") == model)]
    for index, size in enumerate(SIZES):
        selected = sorted([row for row in relevant if math.isclose(equivalent_size(row), size)], key=lambda row: n(row, "f"))
        color = size_color("null" if null else model, size)
        offset = np.exp(np.linspace(-.04, .04, 5)[index]) if not null else 1.
        # A missing designed point interrupts a curve rather than bridging it.
        xs = [n(row, "f") * offset for row in selected]
        ys = [n(row, "alpha") for row in selected]
        ax.plot(xs, ys, color=color, lw=.67, zorder=1)
        if null:continue
        for row, x, y in zip(selected, xs, ys):
            if not finite(y):
                continue
            bounded_point(ax, x, y, n(row, "alpha_ci_lower", "ci_lower"), n(row, "alpha_ci_upper", "ci_upper"), color, open_marker=yes(first(row, "forced_combined", "forced", default=False)), ms=1.7 if null else 1.9)
    if ylabel:
        ax.set_ylabel(r"Local exponent $\alpha$", labelpad=2.)
    facet(ax, NULL_B_NAMES[model] if null else NAMES[model])
    style_axis(ax)


def fold_panel(ax, rows, *, null=False):
    ax.set(xscale="log", yscale="log", xlim=(18, 67), ylim=(1., 4.))
    ax.set_xticks([20, 40, 60], ["20", "40", "60"])
    ax.set_yticks([1, 2, 4], ["1", "2", "4"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.plot([20, 60], [1, 3], "--", color=".7", lw=.55)
    ax.text(43,2.3,r"$\propto C$",color='.55',fontsize=5.6,ha='left',va='bottom')
    groups = NULLS if null else MODELS
    for index, group in enumerate(groups):
        selected = [row for row in rows if (null_kind(row) == group if null else row.get("model") == group)]
        # The main has the requested reference lag; other reference lags remain
        # available in compact CSV and the supplement.
        selected = [row for row in selected if n(row, "reference_lag", "tau_ref", default=200 if group == "spring_ring" else 2000) == (200 if group == "spring_ring" else 2000)]
        selected.sort(key=equivalent_size)
        x = [equivalent_size(row) for row in selected]
        y = [n(row, "fold", "fhalf_fold", "ratio") for row in selected]
        color = "black" if null else COLORS[group]
        style = ("-", "--", ":")[index] if null else "-"
        marker = None if null else ("o", "s", "^", "D")[index]
        legend_name=dict(drift="Drift",same_rate="+ same-rate jumps",perimeter="+ perimeter-rate jumps")[group] if null else NAMES[group]
        ax.plot(x, y, ls=style, color=color, lw=.8, marker=marker, ms=2., mfc="white" if group == "goryachev" else color, mew=.5, label=legend_name)
        if null:continue
        for row, xx, yy in zip(selected, x, y):
            if finite(yy):
                bounded_point(ax, xx, yy, n(row, "fold_ci_lower", "fhalf_fold_ci_lower", "ci_lower"), n(row, "fold_ci_upper", "fhalf_fold_ci_upper", "ci_upper"), color, marker=marker or "o", open_marker=group == "goryachev", ms=0 if null else 2.)
    ax.set_xlabel("Circumference C", labelpad=1.5)
    ax.set_ylabel(r"$f_{1/2}/f_{1/2}(C_{20})$", labelpad=2)
    ax.legend(frameon=False, fontsize=5.15 if null else 5.5, loc="upper left", handlelength=1.6, labelspacing=.15, borderpad=.1)
    style_axis(ax)


def model_examples(fig, grids, data, exported):
    examples = load_kymographs(data, "model_D")
    selection = read_csv(data / "model_D_selection.csv")
    correlations = primary(read_csv(data / "condition_lag_metrics.csv"))
    end_alpha = primary(read_csv(data / "end_to_end_alpha.csv"))
    plotted = []
    for j, model in enumerate(MODELS):
        selected = [row for row in examples if row.get("model") == model]
        selection_row = next((row for row in selection if row.get("model") == model), {})
        chosen_f = n(selection_row, "f", default=n(selected[0], "f") if selected else math.nan)
        outer = grids[j].subgridspec(2, 1, height_ratios=[1.45, 1.], hspace=.40)
        widths=[next((row["circ"] for row in selected if math.isclose(equivalent_size(row),size)),size) for size in (20,60)]
        pair = outer[0].subgridspec(1, 2, width_ratios=widths, wspace=.32)
        vmax = max((row["field"].max() for row in selected), default=1.)
        for k, size in enumerate((20, 60)):
            ax = fig.add_subplot(pair[k])
            record = next((row for row in selected if math.isclose(equivalent_size(row), size)), None)
            if record is not None:
                switch_string = str(first(record, "switch_times", default="[]"))
                try:
                    switches = json.loads(switch_string)
                except ValueError:
                    switches = []
                native_size = size / 20 if model == "spring_ring" else size
                kymo(ax, record, vmax, ylabel=k == 0, label=f"C{native_size:g}", switches=switches,elapsed=True)
                if j==0 and k==1:
                    ax.legend([Line2D([0],[0],color='black',lw=.65)],['switch'],frameon=False,fontsize=5.2,loc='upper right',handlelength=.7,borderpad=.1,labelspacing=.1)
            else:
                ax.text(.5, .5, "No retained\nexemplar", transform=ax.transAxes, ha="center", fontsize=6)
                ax.set_axis_off()
        bounds = grids[j].get_position(fig)
        label = NAMES[model] + (rf", $f={chosen_f:g}$" if finite(chosen_f) else "")
        fig.text((bounds.x0 + bounds.x1) / 2, bounds.y1 + .007, label, ha="center", fontsize=6.8)
        ax = fig.add_subplot(outer[1])
        rr = [row for row in correlations if row.get("model") == model and math.isclose(n(row, "f"), chosen_f, rel_tol=1e-8) and any(math.isclose(equivalent_size(row),size) for size in (20,60))]
        aa = next((row for row in end_alpha if row.get("model") == model and math.isclose(n(row, "f"), chosen_f, rel_tol=1e-8)), {})
        value = n(aa, "alpha", "end_to_end_alpha")
        annotation = rf"$\alpha_{{20,60}}={value:.2g}$" if finite(value) else r"$\alpha_{20,60}$ unavailable"
        lo, hi = n(aa, "alpha_ci_lower"), n(aa, "alpha_ci_upper")
        if not math.isnan(lo) and not math.isnan(hi):
            bound=lambda value:('−∞' if value<0 else '∞') if math.isinf(value) else f'{value:.2g}'
            annotation += f"\n95% CI [{bound(lo)}, {bound(hi)}]"
        correlation(ax, rr, model, annotation=annotation, legend=j == 0)
        if j == 0:
            ax.set_ylabel(r"Heading $C(\ell)$", labelpad=2)
        plotted.extend(rr)
    exported["D_acf"] = plotted


def decomposition(fig, grid, data, exported, omissions):
    rows = [r for r in primary(read_csv(data / "drift_decomposition.csv")) if n(r,"flicker_frames",default=5)==5 and math.isclose(n(r,"MSD_range_fraction",default=.25),.25)]
    # Descriptive multiple-front fractions may be available even when native
    # tracked-region archives are inaccessible. These remain separately named.
    source=primary(read_csv(data/'analysis_lag_metrics.csv'))
    front_rows=[dict(model=r['model'],C_label=r['C_label'],f=r['f'],forced_combined=yes(first(r,'forced_combined','forced',default=False)),multiple_front_fraction=r.get('multiple_front_fraction'),multiple_front_fraction_ci_lower=first(r,'multiple_front_fraction_ci_lower','multiple_fronts_ci_lower'),multiple_front_fraction_ci_upper=first(r,'multiple_front_fraction_ci_upper','multiple_fronts_ci_upper'),n_seeds=r.get('n_seeds'),source=r.get('source')) for r in source]
    if not front_rows:front_rows=primary(read_csv(data/'multiple_front_fraction.csv'))
    write_csv(data/'multiple_front_fraction.csv',front_rows)
    keys = ("Ds", "lambda_eff", "kpred_over_kobs", "multiple_front_fraction")
    labels = (r"$D_s$ (physical)", r"$\lambda_{\mathrm{eff}}$", r"$k_{\mathrm{pred}}/k_{\mathrm{obs}}$", "Multiple-front\nfraction")
    plotted = []
    for i, (key, label) in enumerate(zip(keys, labels)):
        axes = []
        values = [n(row, key) for row in rows if n(row, key) > 0]
        if key == "kpred_over_kobs":
            limits = (.1, 10.)
        elif key == "multiple_front_fraction":
            limits = (0., 1.)
        else:
            limits = (min(values) / 2, max(values) * 2) if values else ((1e-7, 1.) if key == "Ds" else (1e-7, 1e-1))
        for j, model in enumerate(MODELS):
            ax = fig.add_subplot(grid[i, j])
            axes.append(ax)
            ax.set_ylim(limits)
            if i < 3:
                ax.set_yscale("log")
            if key == "kpred_over_kobs":
                ax.axhline(1, color=".65", lw=.5)
            for size in SIZES:
                relevant = [row for row in rows if row.get("model") == model and math.isclose(equivalent_size(row), size)]
                if key == "multiple_front_fraction" and not any(finite(n(row, key, "multiple_fronts_fraction")) for row in relevant):
                    relevant = [row for row in front_rows if row.get("model") == model and math.isclose(equivalent_size(row), size)]
                relevant.sort(key=lambda row: n(row, "f"))
                color = size_color(model, size)
                x, y = [], []
                for row in relevant:
                    value = n(row, key, "multiple_fronts_fraction" if key == "multiple_front_fraction" else key)
                    if key == "kpred_over_kobs" and not finite(value) and n(row, "kobs") > 0:
                        value = n(row, "kpred") / n(row, "kobs")
                    x.append(n(row, "f"))
                    y.append(value if finite(value) and (value >= 0 if i == 3 else value > 0) else math.nan)
                    if not finite(y[-1]):
                        reason=('Empirical zero omitted from logarithmic axis; Poisson zero-count reference bound and whole-seed CI remain distinct in source' if key=='lambda_eff' else 'Empirical zero omitted from logarithmic axis; exact zero preserved in source') if value==0 and i<3 else first(row, key + "_status", "status", "reason", default="Undefined in compact record")
                        if key=='kpred_over_kobs' and (not finite(n(row,'kobs')) or n(row,'kobs')<=0):reason='Observed k unavailable at any eligible single-condition analysis lag'
                        elif key=='kpred_over_kobs' and not finite(n(row,'kpred')):reason='Prediction unavailable because a drift or switch ingredient is undefined'
                        elif str(reason)=='resolved':reason='Quantity is undefined despite other resolved tracked ingredients; inspect source fields'
                        omissions.append(dict(panel="E", model=model, C_label=row.get("C_label"), f=row.get("f"), quantity=key, reason=reason))
                    else:
                        lower=n(row,key+'_ci_lower','multiple_fronts_ci_lower' if i==3 else key+'_ci_lower')
                        upper=n(row,key+'_ci_upper','multiple_fronts_ci_upper' if i==3 else key+'_ci_upper')
                        forced=yes(first(row,'forced_combined','forced',default=False))
                        plotted.append(dict(panel="E", model=model, C_label=row.get("C_label"), f=row.get("f"), forced_combined=forced,quantity=key, value=value,ci_lower=lower,ci_upper=upper,analysis_lag=row.get('analysis_lag'),n_seeds=row.get('n_seeds'),source="multiple_front_fraction.csv" if relevant and row in front_rows else "drift_decomposition.csv"))
                        bounded_point(ax,n(row,'f'),value,lower,upper,color,open_marker=forced,ms=1.6,lw=.35)
                # Draw only the connecting curve here; the point/CI routine
                # above preserves forced-marker fill even without an interval.
                ax.plot(x, y, "-", color=color, lw=.65,zorder=1)
            if not ax.lines or not any(np.isfinite(line.get_ydata()).any() for line in ax.lines if len(line.get_ydata()) > 2):
                if i != 3 and not any(finite(n(row, key)) for row in rows if row.get("model") == model):
                    ax.text(.5, .62, "Not estimable", transform=ax.transAxes, ha="center", va="center", fontsize=5.8, color="black")
            log_f(ax)
            if i != 3:
                ax.tick_params(labelbottom=False)
                ax.set_xlabel("")
            if j == 0:
                ax.set_ylabel(label, labelpad=2)
            else:
                # Every model column uses the same limits within this row.
                # Repeated long log labels otherwise reach the previous facet.
                ax.tick_params(axis='y',which='both',labelleft=False)
                ax.yaxis.get_offset_text().set_visible(False)
            if i == 0:
                facet(ax, NAMES[model])
                if j == 0:
                    handles=[Line2D([],[],marker='o',ls='none',color='black',mfc=fill,ms=2.6,mew=.55,label=label) for fill,label in [('black','Unforced'),('white','Forced')]]
                    ax.legend(handles=handles,frameon=False,loc='upper left',ncol=2,fontsize=5.1,handlelength=.7,handletextpad=.3,columnspacing=.7,borderpad=.15)
            style_axis(ax)
    exported["E"] = plotted


def render(out=OUT):
    data, meta = out / "supporting/data", out / "supporting/metadata"
    data.mkdir(parents=True, exist_ok=True);meta.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":STYLE["fontsize"], "axes.labelsize":STYLE["labelsize"], "axes.titlesize":STYLE["titlesize"], "axes.linewidth":.5, "xtick.labelsize":5.9, "ytick.labelsize":5.9, "svg.fonttype":"none", "pdf.fonttype":42, "savefig.facecolor":"white"})
    fig = plt.figure(figsize=(180/25.4, 500/25.4))
    # Dedicated horizontal space before the fold facets prevents ylabel clashes.
    outer = fig.add_gridspec(5,1,left=.082,right=.977,bottom=.03,top=.965,height_ratios=[2.7,1.3,1.5,2.7,3.7],hspace=.28)
    exported, omissions = {}, []
    null_panels(fig, outer[0].subgridspec(1,3,wspace=.33), data, exported)
    bgrid = outer[1].subgridspec(1,5,width_ratios=[1,1,1,.15,1.08],wspace=.32)
    b_alpha = read_csv(data / "null_B_local_alpha.csv")
    for j, kind in enumerate(NULLS):
        alpha_panel(fig.add_subplot(bgrid[j]),b_alpha,kind,null=True,first_axis=j==0,ylabel=j==0)
    bfolds=read_csv(data/"null_B_noise_tolerance_folds.csv") or read_csv(data/"null_B_fhalf.csv")
    fold_panel(fig.add_subplot(bgrid[4]),bfolds,null=True)
    exported["B_alpha"] = [dict(r,f_display=n(r,'f'),C_equivalent=equivalent_size(r)) for r in b_alpha]
    exported['B_folds']=[dict(r,C_display=equivalent_size(r)) for r in bfolds if n(r,'reference_lag','tau_ref',default=2000)==2000]
    cgrid = outer[2].subgridspec(1,6,width_ratios=[1,1,1,1,.17,1.13],wspace=.35)
    c_alpha = primary(read_csv(data/"local_alpha.csv"))
    for j, model in enumerate(MODELS):
        alpha_panel(fig.add_subplot(cgrid[j]),c_alpha,model,ylabel=j==0)
    folds = primary(read_csv(data/"noise_tolerance_folds.csv"))
    fold_panel(fig.add_subplot(cgrid[5]),folds)
    exported["C_alpha"]=[dict(r,f_display=n(r,'f')*np.exp(np.linspace(-.04,.04,5)[int(np.argmin(abs(np.asarray(SIZES)-equivalent_size(r))))]),C_equivalent=equivalent_size(r)) for r in c_alpha]
    exported["C_folds"]=[dict(r,C_display=equivalent_size(r)) for r in folds if n(r,'reference_lag')==(200 if r.get('model')=='spring_ring' else 2000)]
    model_examples(fig,outer[3].subgridspec(1,4,wspace=.39),data,exported)
    decomposition(fig,outer[4].subgridspec(4,4,wspace=.31,hspace=.20),data,exported,omissions)
    for letter,index in zip("ABCDE",range(5)):
        bounds=outer[index].get_position(fig)
        panel_label(fig,letter,bounds.y1+.018)
    # Family keys sit in the dedicated white space immediately above their row.
    by=outer[1].get_position(fig).y1+.012
    ramp_key(fig,["null"],by,x=.52,width=.12)
    fig.text(.649,by+.002,"C",fontsize=5.8,ha="left",va="center")
    cy=outer[2].get_position(fig).y1+.012
    for j,model in enumerate(MODELS):
        bounds=cgrid[j].get_position(fig)
        ramp_key(fig,[model],cy,x=bounds.x0+.12*bounds.width,width=bounds.width*.87)
    fig.text(.81,cy+.002,"C; Spring ×20 in fold panel",fontsize=5.5,ha="left",va="center")
    for extension in ("png","pdf","svg"):
        fig.savefig(out/("figure4_main."+extension),dpi=450 if extension=="png" else 300)
    plt.close(fig)
    for panel,rows in exported.items():
        write_csv(data/("plotted_"+panel+".csv"),rows)
    alpha_failures=[dict(panel="C",**r) for r in c_alpha if not finite(n(r,"alpha"))]
    omissions.extend(alpha_failures)
    write_csv(data/"omitted_estimates.csv",omissions)
    outliers=[dict(panel='C',**r) for r in c_alpha if finite(n(r,'alpha')) and not -1.5<=n(r,'alpha')<=4.]
    outliers.extend(r for r in exported.get('E',[]) if r.get('quantity')=='kpred_over_kobs' and finite(n(r,'value')) and not .1<=n(r,'value')<=10.)
    write_csv(data/'display_outliers.csv',outliers)
    from figure4_round3_supplement_20261022 import supplement
    supplement(data,out/"figure4_supplement.pdf")
    inputs={str(p.relative_to(out)):sha(p) for p in data.iterdir() if p.is_file() and not p.name.startswith("plotted_") and p.name!="omitted_estimates.csv"}
    try: commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip()
    except Exception: commit="unavailable"
    source_files=sorted(Path(__file__).parent.glob('figure4_round3_*_20261022.py'))+[Path(__file__).parent/'rebuild_figure4_round3_20261022.py']
    write_json(meta/"render_receipt.json",dict(generated_utc=datetime.now(timezone.utc).isoformat(),code_commit=commit,code_commit_status='repository HEAD at build; exact source hashes below include working-tree code',code_sha256={str(p.relative_to(REPO)):sha(p) for p in source_files if p.exists()},renderer=str(Path(__file__).relative_to(REPO)),renderer_sha256=sha(Path(__file__)),width_mm=180,height_mm=500,bootstrap_resamples=2000,inputs=inputs,outputs={p.name:sha(p) for p in out.glob("figure4_*.*") if p.is_file()},omission_count=len(omissions),no_imputation=True,all_main_titles_removed=True,time_increases_upward=True))
    return exported,omissions


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=OUT)
    args=parser.parse_args()
    render(args.out)
