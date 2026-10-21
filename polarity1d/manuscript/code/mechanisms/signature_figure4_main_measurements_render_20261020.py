"""Purpose-built compact Figure 4 measurements from completed audited tables.

The caller declares each model, input and measurement. This source never chooses
a main framing, an input that supports a claim, or a native example. Each single
vector page is designed for a main-figure slot, with all original scientific rows
in its own CSV. Root review of the actual completed candidate remains required.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

import signature_figure4_decisions_render_20261020 as D
import signature_figure4_interim_render_20261020 as P
import signature_figure4_stage_render_20261020 as S


MEASUREMENTS = {
    "persistence_ratio": ("size_ratios.csv", "tau_ratio", "Persistence, large / small"),
    "geometry_normalized_persistence": ("size_ratios.csv", "excess_persistence", "Persistence / actual geometry"),
    "global_heading_mobility_ratio": ("size_ratios.csv", "Dx_large_over_small", "Global-heading mobility, large / small"),
    "fingerprint_region_count": ("fingerprint_effects.csv", "estimate", "Observed region count, large / small"),
    "fingerprint_normalized_contrast": ("fingerprint_effects.csv", "estimate", "Normalized contrast, large / small"),
    "native_multiregion_fraction": ("condition_metrics.csv", "multiple_fronts_fraction", "Native multiregion fraction"),
    "native_low_contrast_fraction": ("condition_metrics.csv", "low_contrast_fraction", "Native low-contrast fraction"),
    "persistence_ratio_overview": ("size_ratios.csv", "tau_ratio", "Full persistence-ratio overview"),
    "persistence_ratio_all_views": ("size_ratios.csv", "tau_ratio", "Persistence across all declared views"),
}
FINGERPRINTS = {"fingerprint_region_count": "region_count",
                "fingerprint_normalized_contrast": "normalized_contrast"}
MARKERS = ("o", "s")


def one(rows, predicate, description):
    selected = [dict(row) for row in rows if predicate(row)]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one declared {description}; found {len(selected)}")
    return selected[0]


def selected_rows(tables, panel):
    """Project declared exact rows; never drop forced, bounded or missing values."""
    measurement = panel["measurement"]
    if measurement not in MEASUREMENTS:
        raise ValueError("Unknown compact measurement")
    arm, view = panel["arm"], panel.get("view", "main")
    inclusive = panel.get("forced_inclusive", False)
    if type(inclusive) is not bool:
        raise ValueError("The forcing view must be an explicit boolean")
    if measurement == "persistence_ratio_all_views":
        if panel.get("views",list(D.VIEWS)) != list(D.VIEWS):
            raise ValueError("All seven predeclared views must remain in the main robustness comparison")
        rows=[]
        for which_view in D.VIEWS:
            projection=dict(panel,measurement="persistence_ratio",view=which_view)
            for row in selected_rows(tables,projection):
                row.update(display_measurement=measurement,display_view=which_view,
                           display_all_seven_views=True)
                rows.append(row)
        return rows
    if measurement == "persistence_ratio_overview":
        if panel.get("conditions"):
            raise ValueError("The full overview cannot select a favorable subset of model/input cells")
        ratios = [dict(row) for row in tables["size_ratios.csv"] if row.get("arm") == arm and row.get("view") == view]
        keys = [(P.model(row), P.num(row,"f")) for row in ratios]
        if set(model for model,_ in keys) != set(P.MODELS) or len(keys) != len(set(keys)):
            raise ValueError("Full overview requires every model and unique complete source model/input cells")
        for row in ratios:
            model, strength = P.model(row), P.num(row,"f")
            matches=[dict(item) for item in tables["condition_metrics.csv"] if item.get("arm") == arm and item.get("view") == view
                     and P.model(item) == model and math.isclose(P.num(item,"f"),strength,rel_tol=0,abs_tol=1e-14)
                     and math.isclose(P.num(item,"C_label"),P.num(row,"large_C_label"),rel_tol=0,abs_tol=1e-10)]
            if len(matches)==1:
                large=matches[0]
            elif not matches and row.get("comparison_status")=="unavailable_unmatched_designed_input" and math.isclose(
                    P.num(row,"missing_size"),P.num(row,"large_C_label"),rel_tol=0,abs_tol=1e-10):
                large=None
            else:
                raise ValueError("Overview large-cell source is unexpectedly missing or duplicated")
            forced=row.get("forced_combined")
            if P.boolean((large or {}).get("forced_combined")) is True:forced=True
            row.update(display_measurement=measurement,display_source_table="size_ratios.csv",display_arm=arm,display_view=view,
                       display_forced_inclusive=inclusive,display_large_condition=json.dumps(large,sort_keys=True),
                       display_large_native_regime=(large or {}).get("regime"),display_operationally_forced=forced)
        return ratios
    declarations = panel["conditions"]
    keys = [(row["model"], float(row["f"])) for row in declarations]
    if not keys or len(keys) > 8 or len(set(keys)) != len(keys):
        raise ValueError("Declare one to eight unique model/input comparisons")
    table_name = MEASUREMENTS[measurement][0]
    result = []
    for group, (model, strength) in enumerate(keys):
        if model not in P.MODELS or not math.isfinite(strength) or strength <= 0:
            raise ValueError("Every declared comparison needs a known model and positive native input")
        def matches(row):
            return (row.get("arm") == arm and row.get("view") == view and P.model(row) == model
                    and math.isclose(P.num(row, "f"), strength, rel_tol=0, abs_tol=1e-14))
        ratio = one(tables["size_ratios.csv"], matches, "small/large source comparison")
        if table_name == "condition_metrics.csv":
            sizes = [(side, P.num(ratio, side + "_C_label")) for side in ("small", "large")]
            chosen = [dict(one(tables[table_name], lambda row:matches(row) and
                          math.isclose(P.num(row, "C_label"), size, rel_tol=0, abs_tol=1e-10),
                          side + " native condition"), display_size_role=side) for side, size in sizes]
        elif table_name == "fingerprint_effects.csv":
            chosen = [one(tables[table_name], lambda row:matches(row) and
                          row.get("statistic") == FINGERPRINTS[measurement], "within-model optical quotient")]
        else:
            chosen = [ratio]
        for row in chosen:
            # Final ratio/fingerprint tables may contain both forcing views.
            # Where the scientific producer declares that axis, it is decisive.
            if "forced_inclusive" in row and P.boolean(row["forced_inclusive"]) is not inclusive:
                raise ValueError("The selected source row uses a different forcing view")
            row.update(display_group=group, display_measurement=measurement,
                       display_arm=arm, display_view=view, display_forced_inclusive=inclusive,
                       display_source_table=table_name, display_small_C_label=ratio["small_C_label"],
                       display_large_C_label=ratio["large_C_label"],
                       display_operationally_forced=row.get("forced_combined", ratio.get("forced_combined")),
                       display_comparison_horizons=json.dumps({side:ratio.get("native_horizon_"+side)
                                                              for side in ("small", "large")}, sort_keys=True))
            result.append(row)
    return result


def endpoints(row, measurement):
    """Keep empirical bounds and whole-seed confidence bounds distinct."""
    _, field, _ = MEASUREMENTS[measurement]
    if measurement in FINGERPRINTS:
        ci = (P.num(row, "ci_lower"), P.num(row, "ci_upper"))
        infinite = P.boolean(row.get("ci_upper_unbounded")) is True
        bounds = (math.nan, math.nan)
    elif measurement == "global_heading_mobility_ratio":
        ci = (P.num(row, "Dx_ratio_ci_lower"), P.num(row, "Dx_ratio_ci_upper"))
        infinite = math.isinf(ci[1])
        bounds = (math.nan, math.nan)
    elif measurement.startswith("native_"):
        stem = field.removesuffix("_fraction")
        ci = (P.num(row, stem+"_ci_lower"), P.num(row, stem+"_ci_upper"))
        infinite = False
        bounds = (math.nan, math.nan)
    else:
        ci = (P.num(row, field+"_ci_lower"), P.num(row, field+"_ci_upper"))
        infinite = P.boolean(row.get(field+"_ci_upper_unbounded")) is True or math.isinf(ci[1])
        bounds = (P.num(row, field+"_lower"), P.num(row, field+"_upper"))
        if P.boolean(row.get(field+"_upper_unbounded")) is True:
            bounds = (bounds[0], math.inf)
        if measurement == "geometry_normalized_persistence":
            geometry = P.num(row, "actual_geometry_ratio")
            if geometry > 0 and math.isfinite(geometry):
                bounds = tuple(P.num(row, "tau_ratio_"+side)/geometry for side in ("lower", "upper"))
                if P.boolean(row.get("tau_ratio_upper_unbounded")) is True:
                    bounds = (bounds[0], math.inf)
                row["display_empirical_bound_identity"] = "tau_ratio empirical bounds / actual_geometry_ratio"
    return P.num(row, field), bounds, (ci[0], math.inf if infinite else ci[1])


def measured_mark(ax, x, row, measurement, color, marker):
    point, (lo, hi), (ci_lo, ci_hi) = endpoints(row, measurement)
    fraction = measurement.startswith("native_")
    signed = measurement == "global_heading_mobility_ratio"
    if point < 0 and not signed:
        raise ValueError("A negative persistence/fingerprint value cannot silently become unavailable")
    valid = lambda value:math.isfinite(value) and (signed or (value >= 0 if fraction else value > 0))
    annotation = ""
    if valid(point):
        ax.plot(x, point, marker, ms=4.5, color=color, mfc="white" if marker == "o" else color, zorder=4)
        if valid(ci_lo) and valid(ci_hi):
            ax.vlines(x, ci_lo, ci_hi, color=color, lw=.8)
            ax.plot([x, x], [ci_lo, ci_hi], "_", color=color, ms=8, mew=1.2, zorder=6)
        elif math.isinf(ci_hi):
            if valid(ci_lo):ax.vlines(x, ci_lo, point, color=color, lw=.8, ls=":")
            annotation = "CI[0,∞)" if ci_lo == 0 else "CI∞"
            ax.annotate(annotation, (x, point), xytext=(3, 4), textcoords="offset points", fontsize=10.5, color=color)
        elif not fraction and ci_lo == 0 and valid(ci_hi):
            ax.vlines(x, point, ci_hi, color=color, lw=.8)
            annotation = "CI starts at 0"
            ax.annotate(annotation, (x, point), xytext=(3, 4), textcoords="offset points", fontsize=10.5, color=color)
        elif not math.isfinite(ci_lo) and not math.isfinite(ci_hi):
            annotation = "CI unavailable"
            ax.annotate(annotation, (x, point), xytext=(3, 4), textcoords="offset points", fontsize=10.5, color=color)
        status = "native point"
    elif point == 0 and not fraction:
        label = "0"
        if math.isfinite(ci_lo) and (math.isfinite(ci_hi) or math.isinf(ci_hi)):
            annotation = "CI["+f"{ci_lo:g}"+","+("∞" if math.isinf(ci_hi) else f"{ci_hi:g}")+"]"
            label += "\n"+annotation
        ax.text(x, .04, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=10.5, color=color)
        status = "native zero; text at axis foot, no substituted log coordinate"
    elif valid(lo) and math.isinf(hi):
        ax.plot(x, lo, "^", color=color, mfc="white", ms=5)
        status = "empirical lower bound"
    elif lo == 0 and valid(hi):
        ax.plot(x, hi, "v", color=color, mfc="white", ms=5)
        status = "empirical upper bound"
    elif valid(lo) and valid(hi):
        ax.vlines(x, lo, hi, color=color, lw=1.1)
        ax.plot([x, x], [lo, hi], "_", color=color, ms=5)
        status = "empirical interval, no point"
    else:
        label = "[0,∞)" if lo == 0 and math.isinf(hi) else "[0,0]" if lo == hi == 0 else "unavailable"
        ax.text(x, .04, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=10.5, color=color)
        status = label
    if not valid(point) and math.isinf(ci_hi):
        # An empirical lower bound and a bootstrap confidence lower endpoint
        # are distinct observations. Retain both at their actual coordinates.
        anchor=lo if valid(lo) else hi if valid(hi) else ci_lo
        if valid(ci_lo):
            ax.plot(x,ci_lo,'_',color=color,ms=5)
            if valid(anchor):ax.vlines(x,ci_lo,anchor,color=color,lw=.8,ls=':')
        annotation='CI[0,∞)' if ci_lo==0 else 'CI∞'
        if valid(anchor):
            ax.annotate(annotation,(x,anchor),xytext=(3,4),textcoords='offset points',fontsize=10.5,color=color)
        row['display_CI_lower_endpoint_retained']=ci_lo if math.isfinite(ci_lo) else ''
    row.update(display_point_status=status, display_CI_annotation=annotation,
               display_x=x, display_scientific_values_unchanged=True)


def legend_notes(measurement):
    notes=["Grey marks retain operationally forced conditions. Uncertainty uses whole-seed 95% bootstrap intervals where available.",
           "Open upward/downward triangles retain native empirical lower/upper bounds. CI∞ marks an upper-unbounded bootstrap confidence interval; CI[0,∞) spans the full nonnegative range. Native zero or unidentified [0,∞) values are text at the logarithmic axis foot; no finite positive endpoint is substituted."]
    if measurement.startswith("native_"):
        notes.append("Open circles denote small cells and squares large cells. Native time-window fractions are not event rates.")
    if measurement == "global_heading_mobility_ratio":
        notes.append("The physical-mobility quantity is a short-lag global activity-heading proxy using actual circumference. In multiple-region states it is not diffusion of a tracked single region; qualified literal A3 site MSDs are separate.")
    if measurement in FINGERPRINTS:
        notes.append("Optical quantities are within-model large/small quotients, not cross-model comparisons of raw counts. The primary optical sampling and smoothing use their declared physical scales.")
    if measurement == "persistence_ratio_overview":
        notes=[notes[0],"Ρ = τₚ(largest) / τₚ(smallest), from direct measured heading persistence. The full overview includes every designed model/input comparison. Cell colour shows only a finite measured direct persistence ratio, on a logarithmic colour scale centred at one. Bound or unidentified cells receive no fabricated ratio colour. Up/down arrows denote empirical one-sided bounds, a bar a finite empirical interval, and ? an unidentified or unavailable ratio or regime. A dagger denotes a finite point with upper-unbounded bootstrap confidence interval. Hatching denotes an input not designed for that model. Grey retains operationally forced cells. Letters show the largest-cell native raw regime: S single front, M multiple fronts, L low contrast, X mixed. All empirical/CI bounds, native regime fractions, forcing flags and unavailable reasons remain in the panel CSV."]
    if measurement == "persistence_ratio_all_views":
        notes.append("Every declared temporal window and sampling cadence is shown. Circles are native finite ratios; right-pointing open triangles are empirical lower bounds, without a substituted ratio point. Horizontal bars retain whole-seed 95% intervals. A dagger denotes an upper-unbounded confidence interval; its finite confidence lower endpoint remains a vertical cap, distinct from the empirical bound. Main, first half, second half and late quarter refer to the declared analysis windows. Sampling ×2/×4/×8 retains the same independent seed histories.")
    return notes


def all_views_page(rows,panel,output):
    """Three compact horizontal interval columns retain all seven views."""
    declarations=panel['conditions']
    if not 1<=len(declarations)<=3:raise ValueError('All-view main columns permit at most three reviewed model/input conditions')
    expected={(i,view) for i in range(len(declarations)) for view in D.VIEWS}
    if {(int(row['display_group']),row['display_view']) for row in rows}!=expected or len(rows)!=len(expected):
        raise ValueError('Main robustness display lost or duplicated a declared view')
    fig,axes=plt.subplots(1,len(declarations),figsize=(6.25,2.7),squeeze=False,sharey=True)
    finite=[]
    for row in rows:
        point,bounds,ci=endpoints(row,'persistence_ratio')
        finite.extend(value for value in (point,*bounds,*ci) if value>0 and math.isfinite(value))
    low=min([1.,*finite])*.86;high=max([1.,*finite])*1.28
    labels=('Main','First half','Second half','Late quarter','Sampling ×2','Sampling ×4','Sampling ×8')
    for i,decl in enumerate(declarations):
        ax=axes[0,i]
        for row in rows:
            if int(row['display_group'])!=i:continue
            y=D.VIEWS.index(row['display_view']);point,(lo,hi),(ci_lo,ci_hi)=endpoints(row,'persistence_ratio')
            color='#888888' if P.boolean(row.get('display_operationally_forced')) is True else '#24649a'
            valid=lambda value:value>0 and math.isfinite(value)
            if valid(point):
                ax.plot(point,y,'o',ms=4.5,color=color,mfc='white');anchor=point;status='native point'
            elif valid(lo) and math.isinf(hi):
                ax.plot(lo,y,'>',ms=5,color=color,mfc='white');anchor=lo;status='empirical lower bound'
            elif lo==0 and valid(hi):
                ax.plot(hi,y,'<',ms=5,color=color,mfc='white');anchor=hi;status='empirical upper bound'
            elif valid(lo) and valid(hi):
                ax.hlines(y,lo,hi,color=color);anchor=lo;status='empirical interval without point'
            else:
                ax.text(.5,y,'?',transform=ax.get_yaxis_transform(),ha='center',fontsize=10.5);anchor=math.nan;status='unidentified or unavailable'
            if valid(ci_lo) and valid(ci_hi):
                ax.hlines(y,ci_lo,ci_hi,color=color,lw=.8)
                ax.plot([ci_lo,ci_hi],[y,y],'|',color=color,ms=5)
            elif math.isinf(ci_hi):
                if valid(ci_lo):
                    ax.plot(ci_lo,y,'|',color=color,ms=5)
                    if valid(anchor):ax.hlines(y,ci_lo,anchor,color=color,lw=.8,ls=':')
                if valid(anchor):ax.annotate('†',(anchor,y),xytext=(4,1),textcoords='offset points',fontsize=10.5,color=color)
            row.update(display_point_status=status,display_CI_annotation='upper-unbounded †' if math.isinf(ci_hi) else 'finite confidence interval',
                       display_CI_lower_endpoint_retained=ci_lo if math.isfinite(ci_lo) else '',
                       display_y=y,display_scientific_values_unchanged=True)
        ax.set_xscale('log');ax.set_xlim(low,high);ax.axvline(1,color='.65',lw=.7,ls='--')
        ax.set_xticks([v for v in (1,2,4,8,16,32) if low<=v<=high]);ax.xaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
        ax.minorticks_off();ax.tick_params(labelsize=10.5);ax.set_title(P.LABELS[decl['model']],fontsize=11.5,pad=5)
        ax.set_ylim(len(D.VIEWS)-.5,-.5);ax.set_yticks(range(len(D.VIEWS)),labels,fontsize=10.5)
    fig.text(.64,.035,'Ρ = persistence, large / small',ha='center',fontsize=11.5)
    fig.text(.03,.95,panel.get('title','All seven views · f = 0.05'),fontsize=12)
    fig.subplots_adjust(left=.23,right=.97,top=.75,bottom=.23,wspace=.25)
    bounds=S.presentation_bounds(fig)
    if not bounds['passed']:plt.close(fig);raise ValueError('All-view main comparison exceeds actual compact page bounds')
    fig.savefig(output);plt.close(fig);return bounds


def compact_page(rows, panel, output):
    """One readable source page; detailed qualifications belong in its legend."""
    measurement = panel["measurement"]
    fraction = measurement.startswith("native_")
    groups = sorted({int(row["display_group"]) for row in rows})
    if groups != list(range(len(panel["conditions"]))):
        raise ValueError("Compact measurement group coverage changed")
    fig, ax = plt.subplots(figsize=(6.25, 2.7))
    for row in rows:
        role = row.get("display_size_role", "ratio")
        shift = (-.12 if role == "small" else .12) if fraction else 0.
        forced = P.boolean(row.get("display_operationally_forced")) is True
        color = "#888888" if forced else "#24649a"
        measured_mark(ax, int(row["display_group"])+shift, row, measurement, color,
                      MARKERS[1] if role == "large" else MARKERS[0])
    names=dict(wave_pinning="Wave-\npinning",goryachev="Gory-\nachev",holmes_model3="Holmes 3",spring_ring="Spring")
    labels = [names[row["model"]]+"\nf="+f"{float(row['f']):g}" for row in panel["conditions"]]
    ax.set_xticks(groups, labels, fontsize=10.5)
    ax.tick_params(axis="y", labelsize=10.5)
    ax.set_xlim(-.5, len(groups)-.5)
    ax.set_title(panel.get("title", MEASUREMENTS[measurement][2]), loc="left", fontsize=12, pad=5)
    ax.set_ylabel("Native fraction" if fraction else "Within-model ratio", fontsize=11.5)
    if fraction:
        ax.set_ylim(-.03, 1.08)
        ax.plot([], [], "o", ms=4, color="#24649a", mfc="white", label="Small")
        ax.plot([], [], "s", ms=4, color="#24649a", label="Large")
        ax.legend(fontsize=10.5, loc="upper right", frameon=False, ncol=2)
    else:
        ax.set_yscale("symlog",linthresh=.1) if measurement=="global_heading_mobility_ratio" else ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, position:f"{value:g}"))
        ax.axhline(1, color=".65", lw=.6, ls="--")
        ax.margins(y=.3)
    if panel.get('show_qualification_footer',True):
        footer=fig.text(.16, .02, "Grey: operationally forced. Native bounds and CI∞: see legend.", fontsize=10.5, va="bottom")
        footer.set_gid("document_footer")
    fig.subplots_adjust(left=.16, right=.95, top=.84, bottom=.40)
    bounds = S.presentation_bounds(fig)
    if not bounds["passed"]:
        plt.close(fig)
        raise ValueError("Compact measurement text exceeds its actual page bounds")
    fig.savefig(output)
    plt.close(fig)
    return bounds


def overview_page(rows, panel, output):
    """Full model-by-input direct-ratio and native-regime context, without selection."""
    strengths=sorted({P.num(row,"f") for row in rows})
    lookup={(P.model(row),P.num(row,"f")):row for row in rows}
    if len(lookup)!=len(rows) or set(P.model(row) for row in rows)!=set(P.MODELS):
        raise ValueError("Full ratio overview needs all four models and every unique declared input")
    missing=[dict(model=model,f=strength,source_status="input not designed for this model")
             for model in P.MODELS for strength in strengths if (model,strength) not in lookup]
    for row in rows:
        row.update(display_grid_model_order=json.dumps(P.MODELS),display_grid_inputs=json.dumps(strengths),
                   display_not_designed_cells=json.dumps(missing,sort_keys=True))
    values=[math.log2(P.num(row,"tau_ratio")) for row in rows if P.num(row,"tau_ratio")>0 and math.isfinite(P.num(row,"tau_ratio"))]
    low=min(-1.,min(values,default=-1.));high=max(1.,max(values,default=1.))
    norm=TwoSlopeNorm(vmin=low,vcenter=0.,vmax=high);cmap=plt.get_cmap("RdBu_r")
    page_width=float(panel.get('width_inches',6.25))
    if not 6.25<=page_width<=12.:raise ValueError('Overview width must preserve the declared readable main scale')
    fig,ax=plt.subplots(figsize=(page_width,2.7))
    regimes={"single_front":"S","single":"S","multiple_fronts":"M","multiple":"M","low_contrast":"L","mixed":"X"}
    for i,model in enumerate(P.MODELS):
        for j,strength in enumerate(strengths):
            row=lookup.get((model,strength))
            if row is None:
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor="white",edgecolor=".8",hatch="///",lw=.5))
                continue
            point,(lo,hi),ci=endpoints(row,"persistence_ratio")
            forced=P.boolean(row.get("display_operationally_forced")) is True
            face="#aaaaaa" if forced else cmap(norm(math.log2(point))) if point>0 and math.isfinite(point) else "#f3f3f3"
            ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=face,edgecolor="white",lw=.7))
            if point>0 and math.isfinite(point):
                status="finite direct ratio colour";mark="†" if math.isinf(ci[1]) else ""
            elif point==0:
                status="native zero";mark="0"
            elif lo>0 and math.isfinite(lo) and math.isinf(hi):
                status="empirical lower bound";mark="↑"
            elif lo==0 and hi>0 and math.isfinite(hi):
                status="empirical upper bound";mark="↓"
            elif lo>0 and hi>0 and math.isfinite(lo) and math.isfinite(hi):
                status="empirical interval without point";mark="↕"
            else:
                status="unidentified or unavailable";mark="?"
            regime=regimes.get(row.get("display_large_native_regime"),"?")
            label=mark+"\n"+regime if mark else regime
            text_color="white" if not forced and point>0 and math.isfinite(point) and abs(norm(math.log2(point))-.5)>.32 else "black"
            ax.text(j,i,label,ha="center",va="center",fontsize=10.5,color=text_color,linespacing=.9)
            row.update(display_point_status=status,display_cell_row=i,display_cell_column=j,
                       display_ratio_colour_finite_point_only=True,display_colour_log2_limits=json.dumps([low,high]),
                       display_regime_symbol=regime,display_bound_or_CI_symbol=mark,display_scientific_values_unchanged=True)
    ax.set_xlim(-.5,len(strengths)-.5);ax.set_ylim(len(P.MODELS)-.5,-.5)
    ax.set_yticks(range(len(P.MODELS)),["Wave-pinning","Goryachev","Holmes 3","Spring"],fontsize=10.5)
    ax.set_xticks(range(len(strengths)),[f"{value:g}" for value in strengths],fontsize=10.5,rotation=60,ha="right")
    ax.tick_params(length=0)
    ax.set_xlabel("Native input f",fontsize=11.5,labelpad=1)
    ax.set_title(panel.get("title","All inputs · persistence ratio and large-cell regime"),loc="left",fontsize=12,pad=32)
    cax=fig.add_axes([.25,.84,.62,.025])
    bar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal")
    ticks=[low,0.,high];bar.set_ticks(ticks,labels=[f"{2**value:.2g}" for value in ticks]);bar.ax.tick_params(labelsize=10.5,pad=1,length=2)
    fig.text(.03,.845,"Ρ",fontsize=12,va="center")
    fig.subplots_adjust(left=.13 if page_width>=10 else .20,right=.97,top=.75,bottom=.36)
    bounds=S.presentation_bounds(fig)
    if not bounds["passed"]:
        plt.close(fig);raise ValueError("Full overview exceeds its actual compact page bounds")
    fig.savefig(output);plt.close(fig)
    return bounds


def render(spec_path, output, name):
    spec_path, output = Path(spec_path).resolve(), Path(output).resolve()
    spec = json.loads(spec_path.read_text())
    if spec.get("schema") != "figure4_compact_measurements_v1":
        raise ValueError("Unknown compact-main source contract")
    tables, sources = D.load_final(spec["final_statistics_root"], spec["producer_sha256"], spec["independent_audit_sha256"])
    panels = spec["panels"]
    identities=[row['id'] for row in panels]
    if not panels or len(panels)>4 or identities!=sorted(set(identities)) or not set(identities)<=set('BCDE'):
        raise ValueError("Declare one to four unique ordered compact measurement components")
    records = []
    for panel in panels:
        stem = name+"_panel"+panel["id"]
        pdf = output/"supporting/vectors"/(stem+".pdf")
        csv_path = output/"supporting/data"/(stem+".csv")
        audit = output/"supporting/metadata"/(stem+"_display_audit.json")
        if any(path.exists() for path in (pdf, csv_path, audit)):
            raise ValueError("Use a fresh named compact source closure")
        for path in (pdf, csv_path, audit):path.parent.mkdir(parents=True, exist_ok=True)
        rows = selected_rows(tables, panel)
        bounds = (overview_page(rows,panel,pdf) if panel["measurement"]=="persistence_ratio_overview" else
                  all_views_page(rows,panel,pdf) if panel["measurement"]=="persistence_ratio_all_views" else compact_page(rows, panel, pdf))
        columns = list(dict.fromkeys(key for row in rows for key in row))
        with csv_path.open("w", newline="") as stream:
            writer=csv.DictWriter(stream, fieldnames=columns);writer.writeheader();writer.writerows(rows)
        record = dict(passed=True, complete=True, status="candidate_component", final_export_eligible=False,
            generated_utc=datetime.now(timezone.utc).isoformat(), id=panel["id"], measurement=panel["measurement"],
            small_large_comparison=True, panel_spec=panel, scientific_sources=sources,
            pdf=dict(path=str(pdf), sha256=P.sha(pdf), pages=1),
            csvs=[dict(path=str(csv_path), sha256=P.sha(csv_path), rows=len(rows))], presentation_bounds=bounds,
            required_legend_notes=legend_notes(panel["measurement"]),minimum_printed_main_width_inches=12.,
            recommended_printed_main_width_inches=14.,minimum_core_source_font_pt=10.5,
            source_spec=dict(path=str(spec_path), sha256=P.sha(spec_path)), renderer_sha256=P.sha(__file__),
            shared_source_sha256={str(Path(module.__file__).resolve()):P.sha(module.__file__) for module in (D, P, S)},
            no_framing_or_condition_selection_by_renderer=True, no_scientific_estimator_change=True)
        P.write_json(audit, record)
        records.append(dict(id=panel["id"], measurement=panel["measurement"], small_large_comparison=True,
            required_legend_notes=record["required_legend_notes"],minimum_printed_main_width_inches=12.,
            pdf=record["pdf"], csvs=record["csvs"], display_audit=dict(path=str(audit),sha256=P.sha(audit),assertions=dict(passed=True,complete=True))))
    return records


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", required=True)
    args=parser.parse_args()
    print(json.dumps(render(args.spec,args.output,args.name), ensure_ascii=False))
