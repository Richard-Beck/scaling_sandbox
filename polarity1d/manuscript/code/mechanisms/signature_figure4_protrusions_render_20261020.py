"""Render A5 protrusion measurements without selecting a Figure 4 framing.

The source is an audited tidy-table closure.  Every page retains all its source
rows in CSV, including unavailable estimands, independent-seed denominators,
native horizons and the complete threshold/cadence/optics factorial.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from signature_figure4_interim_render_20261020 import (
    MODELS, LABELS, boolean, curve, main_rows, model, num, polish, read_csv,
    sha, size_color, subset, write_json,
)


NAMES = {
    "region_count_mean": "Mean simultaneous qualified-region count",
    "single_front_fraction": "Time fraction with one qualified region",
    "multiple_front_fraction": "Time fraction with multiple qualified regions",
    "low_contrast_time_fraction": "Time fraction below the contrast gate",
    "births_per_time_per_perimeter": "Region births per native time per physical perimeter",
    "mean_region_width": "Mean physical region width",
    "mean_region_width_fraction": "Mean region width / actual circumference",
    "mean_birth_angle": "Mean newborn angle from the previous dominant region",
    "mean_birth_physical_distance": "Mean newborn physical arc distance",
    "opposite_birth_fraction": "Fraction of births opposite the previous dominant region",
    "observed_births": "Observed region births per trajectory",
    "observed_takeovers": "Observed globally dominant newborn takeovers per trajectory",
    "observed_takeovers_per_birth": "Globally dominant newborn takeovers per birth",
    "observed_pairwise_takeovers": "Literal pairwise takeover sensitivity per trajectory",
    "observed_pairwise_takeovers_per_birth": "Literal pairwise takeover sensitivity per birth",
    "unscorable_takeover_tracks": "Tracks whose takeover cannot be scored",
    "ever_unscorable_takeover_fraction": "Fraction of trajectories with an unscorable takeover track",
    "restricted_mean_lifetime_lower": "Lower-censoring restricted mean lifetime",
    "restricted_mean_lifetime_observed": "Observed-frame restricted mean lifetime",
    "restricted_mean_lifetime_upper": "Upper-censoring restricted mean lifetime",
    **{f"mean_absolute_turn_{i}_frames": f"Mean absolute dominant-heading turn over {i} observation frames" for i in (1,2,4,8)},
}


def identity(row):
    return (model(row), num(row,"C_label"), num(row,"f"), row["kind"],
            num(row,"sigma"), row["experiment"])


def variant(row):
    return num(row,"threshold"), int(num(row,"downsample_factor")), boolean(row["smoothed"])


def validate_metric_identities(rows):
    """A sink/table name must never replace a native condition's kind."""
    if not rows or {row["kind"] for row in rows}!={"noise","mean","zero"}:
        raise ValueError("A5 condition kinds must remain noise, mean and zero in every variant")
    keys=[(identity(row),variant(row),row["statistic"]) for row in rows]
    if len(keys)!=len(set(keys)):
        raise ValueError("A5 tidy metric identities are duplicated")
    primary=[r for r in rows if boolean(r["primary"]) is True]
    for row in primary:
        cadence=.8*10 if model(row)=="spring_ring" else 80.
        if (not math.isclose(num(row,"observation_cadence"),cadence) or
                num(row,"threshold")!=.5 or boolean(row["smoothed"]) is not True):
            raise ValueError("A5 primary must use the declared common physical optics/cadence")
    for which_variant in {variant(r) for r in rows}:
        selected=[r for r in rows if variant(r)==which_variant]
        if {r["kind"] for r in selected}!={"noise","mean","zero"}:
            raise ValueError("A5 variant lost its native noise/mean/zero identities")
    return primary


def load_audited_tables(run_root, heading_root, analysis_subdir="analysis"):
    run_root,heading_root=Path(run_root).resolve(),Path(heading_root).resolve()
    complete_path=run_root/"results"/analysis_subdir/"complete.json"
    complete=json.loads(complete_path.read_text())
    native_path=run_root/"results/independent_primary_native_frames_audit.json"
    native=json.loads(native_path.read_text())
    if (complete.get("passed") is not True or complete.get("cases")!=962 or
            native.get("passed") is not True or native.get("cases")!=962 or
            native.get("primary_native_frames")!=193362 or
            complete.get("condition_identity_checked_against_inventory") is not True or
            complete.get("group_kind_counts")!={"noise":110,"mean":110,"zero":12}):
        raise ValueError("E0 A5 native/statistic coverage is incomplete")
    inventory_path=heading_root/"inventory.json"
    if complete["inventory_sha256"]!=sha(inventory_path):
        raise ValueError("A5 must use the corrected heading inventory")
    correction=json.loads((heading_root/"results/basal_layout_correction.json").read_text())
    if not correction.get("passed") or correction["corrected_inventory_sha256"]!=sha(inventory_path):
        raise ValueError("A5 source uses a superseded basal adapter")
    sources=dict(complete=dict(path=str(complete_path),sha256=sha(complete_path)),
        primary_native_audit=dict(path=str(native_path),sha256=sha(native_path)),
        corrected_inventory=dict(path=str(inventory_path),sha256=sha(inventory_path)))
    for shard in native["source_shards"]:
        if sha(shard["path"])!=shard["sha256"]:
            raise ValueError("Independent A5 primary-native audit changed")
    allviews_path=run_root/"results/independent_all_views_audit.json"
    allviews=json.loads(allviews_path.read_text())
    if (not allviews.get("passed") or allviews.get("native_cases")!=962 or allviews.get("views_per_case")!=24
            or allviews.get("declared_view_audits")!=23088 or allviews.get("native_observation_frames_checked")!=4352088):
        raise ValueError("A5 requires the complete independent24-view native audit")
    for shard in allviews["shards"]:
        if sha(shard["path"])!=shard["sha256"]:raise ValueError("A5 all-view independent source audit changed")
    sources["all24_native_views"]=dict(path=str(allviews_path),sha256=sha(allviews_path),scope=allviews["scope"])
    seed=complete["seed_receipts"]
    if sha(seed["path"])!=seed["sha256"]:raise ValueError("A5 seed source receipts changed")
    tables={}
    for name in ("metrics","distributions","survival"):
        source=complete["outputs"][name]
        if sha(source["path"])!=source["sha256"]:raise ValueError("A5 tidy source table changed")
        sources[name]=source
        # All factorial metric rows are retained. Curves use the declared
        # imaging primary; full curve-source files remain hash-linked above.
        with Path(source["path"]).open(newline="") as stream:
            tables[name]=[r for r in csv.DictReader(stream) if name=="metrics" or boolean(r["primary"]) is True]
    primary=validate_metric_identities(tables["metrics"])
    expected_variants={(threshold,factor,smooth) for threshold in (.3,.5,.7) for factor in (1,2,4,8) for smooth in (False,True)}
    identities={identity(r) for r in tables["metrics"]}
    if (len(identities)!=232 or {r["statistic"] for r in tables["metrics"]}!=set(NAMES)
            or {variant(r) for r in tables["metrics"]}!=expected_variants):
        raise ValueError("E0 A5 requires all232 condition identities,24 statistics and24 variants")
    for current in expected_variants:
        selected={(identity(r),r["statistic"]) for r in tables["metrics"] if variant(r)==current}
        if selected!={(key,statistic) for key in identities for statistic in NAMES}:
            raise ValueError("An A5 factorial sensitivity omits a condition or statistic")
    if len(primary)!=232*len(NAMES):
        raise ValueError("Every A5 metric condition requires one declared imaging primary")
    for name in ("distributions","survival"):
        if any(identity(r) not in identities or r["kind"] not in ("noise","mean","zero") for r in tables[name]):
            raise ValueError("A5 curve identities differ from corrected metric identities")
    status=json.loads((heading_root/"results/analysis/status.json").read_text())
    heading_path=heading_root/"results/analysis/condition_metrics.csv"
    if not status.get("complete") or sha(heading_path)!=status["outputs"][heading_path.name]:
        raise ValueError("A5 forced-state annotation lacks verified heading analysis")
    heading={identity(dict(r,kind="noise",sigma="1",experiment="E0")):r for r in main_rows(read_csv(heading_path))}
    sources["heading_condition_metrics"]=dict(path=str(heading_path),sha256=sha(heading_path))
    return tables,heading,sources


def annotate_forced(row,heading):
    if "display_forced" in row:
        return boolean(row["display_forced"]) is True
    key=(model(row),num(row,"C_label"),num(row,"f"),"noise",1.,"E0")
    return boolean(heading.get(key,{}).get("forced")) is True


def size_line_style(index):
    # Preserve the three E0 styles while supporting both new intermediate sizes.
    return ("-", "--", ":", "-.", (0, (3, 1, 1, 1, 1, 1)))[index % 5]


def write_csv(path,rows):
    columns=list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)


class Supplement:
    def __init__(self,output,sources):
        self.output=Path(output).resolve();self.sources=sources;self.pages=[]
        for name in ("data","vectors","metadata"):(self.output/"supporting"/name).mkdir(parents=True,exist_ok=True)
        self.path=self.output/"supporting/vectors/figure4_S9_E0_protrusion_statistics.pdf"
        self.pdf=PdfPages(self.path)

    def save(self,fig,label,title,rows,note):
        if not rows:raise ValueError("A5 supplement page has no source rows")
        csv_path=self.output/f"supporting/data/figure4_S9_E0_{len(self.pages)+1:03d}_{label}.csv"
        write_csv(csv_path,rows)
        fig.suptitle("S9 · E0 provisional  "+title,fontsize=12,x=.06,ha="left",y=.98)
        footer="E0 only; later fine-profile experiments pending. "+note
        fig.text(.06,.016,"\n".join(textwrap.wrap(footer,160))+"\nCSV: "+csv_path.name,fontsize=6.2,va="bottom")
        self.pdf.savefig(fig);plt.close(fig)
        self.pages.append(dict(label=label,title=title,pdf_page=len(self.pages)+1,csv=str(csv_path),csv_sha256=sha(csv_path),rows=len(rows)))

    def close(self):
        self.pdf.close()
        report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,provisional=True,experiment="E0",
            producer_state="stopped",native_integrations_run=0,final_Figure4_framing_selected=False,
            pdf=dict(path=str(self.path),sha256=sha(self.path),pages=len(self.pages)),pages=self.pages,sources=self.sources,
            native_all24_views_pass=True,
            primary_native_audit_scope="193362 primary native frames including strict globally dominant event gates; independently checked24-view sensitivity closure covers4352088 native observation frames. Latent identity completeness is not claimed.",
            source_sha256=sha(Path(__file__)))
        write_json(self.output/"supporting/metadata/S9_E0_display_audit.json",report)
        return report


def metric_series(ax,rows,heading,*,ls="-",label=None,show_mean=False,point_marker="o",label_suffix=""):
    noisy=[r for r in rows if r["kind"]=="noise"]
    for size in sorted({num(r,"C_label") for r in noisy}):
        selected=sorted([r for r in noisy if num(r,"C_label")==size],key=lambda r:num(r,"f"))
        x=np.array([num(r,"f") for r in selected]);y=np.array([num(r,"estimate") for r in selected])
        color=size_color(model(selected[0]),size)
        ax.plot(x,y,ls=ls,color=color,lw=.8,label=(f"C {size:g}" if label is None else label)+label_suffix)
        for xx,yy,row in zip(x,y,selected):
            current=".6" if annotate_forced(row,heading) else color
            if math.isfinite(yy):
                ax.plot(xx,yy,point_marker,color=current,ms=2.4)
                lo,hi=num(row,"ci_lower"),num(row,"ci_upper")
                if math.isfinite(lo) and math.isfinite(hi):ax.vlines(xx,lo,hi,color=current,lw=.65,alpha=.6)
            else:
                lower,upper=num(row,"lower_bound"),num(row,"upper_bound")
                if math.isfinite(lower):
                    ax.plot(xx,lower,"^",color=current,mfc="none",ms=4)
                    if math.isfinite(upper):
                        ax.vlines(xx,lower,upper,color=current,lw=.9)
                        ax.plot(xx,upper,"v",color=current,mfc="none",ms=4)
                else:
                    ax.plot(xx,.025,"x",color=current,ms=3,transform=ax.get_xaxis_transform())
        if show_mean:
            controls=sorted([r for r in rows if r["kind"]=="mean" and num(r,"C_label")==size],key=lambda r:num(r,"f"))
            ax.plot([num(r,"f") for r in controls],[num(r,"estimate") for r in controls],ls="--",color=color,lw=.6,alpha=.5)
            zero=[r for r in rows if r["kind"]=="zero" and num(r,"C_label")==size]
            if len(zero)==1 and math.isfinite(num(zero[0],"estimate")):ax.axhline(num(zero[0],"estimate"),ls=":",color=color,lw=.5,alpha=.5)


def metric_pages(book,tables,heading):
    rows=tables["metrics"]
    for statistic in sorted({r["statistic"] for r in rows}):
        all_rows=[r for r in rows if r["statistic"]==statistic]
        primary=[r for r in all_rows if boolean(r["primary"]) is True]
        fig,axes=plt.subplots(4,1,figsize=(11.7,8.3))
        for ax,which_model in zip(axes,MODELS):
            selected=subset(primary,which_model=which_model)
            metric_series(ax,selected,heading,show_mean=True)
            polish(ax,xlabel="Input fraction f" if which_model==MODELS[-1] else None,ylabel=LABELS[which_model],logx=True)
            ax.legend(fontsize=6,frameon=False,ncol=3)
        fig.subplots_adjust(left=.09,right=.97,top=.91,bottom=.15,hspace=.5)
        takeover_note=("Literal pairwise events are a separate sensitivity and can occur while a third region remains dominant. "
            if statistic.startswith("observed_pairwise_") else
            "Primary takeover requires the newborn to be globally dominant among qualified regions and have pairwise D>0.5 against the previous dominant for two observed frames. ")
        book.save(fig,statistic,NAMES.get(statistic,statistic.replace("_"," ")),primary,
            "Primary optics: threshold 0.5; fixed physical Gaussian σ=1.5 (Spring 0.15); common cadence 80 (Spring 8). Grey points are forced conditions. Dashed mean controls and dotted zero-input references remain shown. CIs resample whole seeds; open triangles preserve identified lifetime bounds and × marks no identified estimate/bound. "+takeover_note+"Complete native peak basins define amounts; lifetimes remain cadence/censoring limited. Units are model-native; absolute cross-model comparison is not used.")
        fig,axes=plt.subplots(4,3,figsize=(12.5,9.2),squeeze=False)
        for i,which_model in enumerate(MODELS):
            selected=subset(all_rows,which_model=which_model)
            primary_factor=2
            choices=(("threshold",(.3,.5,.7)),("downsample_factor",(1,2,4,8)),("smoothed",(False,True)))
            for j,(factor,values) in enumerate(choices):
                ax=axes[i,j]
                for value,ls in zip(values,("-","--",":","-.")):
                    chosen=[]
                    for row in selected:
                        th=num(row,"threshold");ds=int(num(row,"downsample_factor"));sm=boolean(row["smoothed"])
                        if factor=="threshold":keep=th==value and ds==primary_factor and sm is True
                        elif factor=="downsample_factor":keep=th==.5 and ds==value and sm is True
                        else:keep=th==.5 and ds==primary_factor and sm is value
                        if keep:chosen.append(row)
                    metric_series(ax,chosen,heading,ls=ls,label=str(value))
                if i==0:ax.set_title(factor.replace("_"," "),fontsize=9)
                polish(ax,xlabel="Input fraction f" if i==3 else None,ylabel=LABELS[which_model] if j==0 else None,logx=True)
        fig.subplots_adjust(left=.07,right=.97,top=.92,bottom=.16,hspace=.6,wspace=.3)
        book.save(fig,statistic+"_sensitivity",NAMES.get(statistic,statistic)+" · sensitivity",all_rows,
            "Size colours are fixed. One-factor contrasts: threshold 0.3/0.5/0.7 uses solid/dashed/dotted; native-profile downsampling ×1/×2/×4/×8 uses solid/dashed/dotted/dash-dot; optics raw/smoothed uses solid/dashed. Other factors use the common primary. The full 24-way factorial, controls, CIs and unavailable estimands remain in this page's CSV; no combined-frame bootstrap is used.")


def distribution_pages(book,tables,heading):
    rows=[r for r in tables["distributions"] if r["kind"]=="noise"]
    strengths=sorted({num(r,"f") for r in rows})
    for statistic in sorted({r["statistic"] for r in rows}):
        intervals=sorted({num(r,"frame_interval") for r in rows if r["statistic"]==statistic and math.isfinite(num(r,"frame_interval"))}) or [None]
        for interval in intervals:
            selected=[r for r in rows if r["statistic"]==statistic and (interval is None or num(r,"frame_interval")==interval)]
            for start in range(0,len(strengths),3):
                batch=strengths[start:start+3];fig,axes=plt.subplots(4,3,figsize=(11.7,8.3),squeeze=False);source=[]
                for i,which_model in enumerate(MODELS):
                    for j in range(3):
                        ax=axes[i,j]
                        if j>=len(batch):ax.axis("off");continue
                        group=subset(selected,which_model=which_model,strength=batch[j])
                        if not group:ax.text(.5,.5,"Not observed",transform=ax.transAxes,ha="center",fontsize=8)
                        for size in sorted({num(r,"C_label") for r in group}):
                            series=subset(group,size=size)
                            for row in series:
                                row["display_x"]=(num(row,"bin_lower")+num(row,"bin_upper"))/2 if math.isfinite(num(row,"bin_lower")) else num(row,"value")
                            color=".6" if annotate_forced(series[0],heading) else size_color(which_model,size)
                            size_index=sorted({num(r,"C_label") for r in group}).index(size)
                            curve(ax,series,"display_x",("probability",),color=color,label=f"C {size:g}",lowkeys=("ci_lower",),highkeys=("ci_upper",),ls=size_line_style(size_index))
                            if not any(math.isfinite(num(r,"probability")) for r in series):
                                ax.text(.02,.08+.10*sorted({num(r,"C_label") for r in group}).index(size),f"C {size:g}: no eligible items",transform=ax.transAxes,fontsize=6,color=color)
                        ax.set_title(f"{LABELS[which_model]} · f={batch[j]:g}",fontsize=8)
                        polish(ax,xlabel="Count / angle (radians)" if i==3 else None,ylabel="Probability" if j==0 else None)
                        lower=[num(r,"bin_lower") for r in group if math.isfinite(num(r,"bin_lower"))]
                        upper=[num(r,"bin_upper") for r in group if math.isfinite(num(r,"bin_upper"))]
                        if lower and upper:ax.set_xlim(min(lower),max(upper))
                        elif group:
                            values=[num(r,"value") for r in group if math.isfinite(num(r,"value"))]
                            if values:ax.set_xlim(min(values)-.5,max(values)+.5)
                        ax.set_ylim(0,1.03 if not any(math.isfinite(num(r,"probability")) for r in group)
                                    else max(.05,ax.get_ylim()[1]))
                        if j==0 and ax.get_legend_handles_labels()[0]:ax.legend(fontsize=5,frameon=False)
                        source.extend(group)
                fig.subplots_adjust(left=.08,right=.97,top=.92,bottom=.15,hspace=.65,wspace=.33)
                book.save(fig,f"{statistic}_{interval}_{start}",statistic.replace("_"," ")+("" if interval is None else f" · {interval:g} observation frames"),source,
                    "Imaging primary only; size colours match S1 and forced conditions are grey. Solid/dashed/dotted also distinguish small/intermediate/large sizes. Histogram bin centres are displayed; exact bins, finite-seed uncertainty, observed totals and frame intervals remain in CSV. Native coarsely saved profiles limit birth/turn observability. Missing densities remain missing, not zero.")


def survival_pages(book,tables,heading):
    rows=[r for r in tables["survival"] if r["kind"]=="noise"]
    strengths=sorted({num(r,"f") for r in rows})
    for clock in ("lower","observed","upper"):
        for start in range(0,len(strengths),3):
            batch=strengths[start:start+3];fig,axes=plt.subplots(4,3,figsize=(11.7,8.3),squeeze=False);source=[]
            for i,which_model in enumerate(MODELS):
                for j in range(3):
                    ax=axes[i,j]
                    if j>=len(batch):ax.axis("off");continue
                    group=[r for r in subset(rows,which_model=which_model,strength=batch[j]) if r["lifetime_clock"]==clock]
                    if not group:ax.text(.5,.5,"Not observed",transform=ax.transAxes,ha="center",fontsize=8)
                    for size in sorted({num(r,"C_label") for r in group}):
                        series=subset(group,size=size)
                        color=".6" if annotate_forced(series[0],heading) else size_color(which_model,size)
                        size_index=sorted({num(r,"C_label") for r in group}).index(size)
                        curve(ax,series,"time",("survival",),color=color,label=f"C {size:g}",lowkeys=("ci_lower",),highkeys=("ci_upper",),ls=size_line_style(size_index))
                        if not any(math.isfinite(num(r,"survival")) for r in series):
                            ax.text(.02,.08+.10*sorted({num(r,"C_label") for r in group}).index(size),f"C {size:g}: no incident-cohort estimate",transform=ax.transAxes,fontsize=6,color=color)
                    ax.set_title(f"{LABELS[which_model]} · f={batch[j]:g}",fontsize=8)
                    polish(ax,xlabel="Lifetime (native time)" if i==3 else None,ylabel="Survival probability" if j==0 else None)
                    ax.set_ylim(-.02,1.03);source.extend(group)
                    if group:
                        cadence=8 if which_model=="spring_ring" else 80
                        maximum=max(num(row,"time") for row in group)
                        ax.set_xscale("symlog",linthresh=cadence)
                        ax.set_xlim(0,maximum)
                        ticks=[0.,float(cadence)]
                        while ticks[-1]*10<=maximum:ticks.append(ticks[-1]*10)
                        ax.set_xticks(ticks,labels=[f"{value:g}" for value in ticks])
                        ax.set_xticks([],minor=True)
                    if j==0 and ax.get_legend_handles_labels()[0]:ax.legend(fontsize=5,frameon=False)
            fig.subplots_adjust(left=.08,right=.97,top=.92,bottom=.15,hspace=.65,wspace=.33)
            book.save(fig,f"lifetime_{clock}_{start}",f"Incident-region lifetime survival · {clock} observation clock",source,
                "Kaplan–Meier curves use incident tracks; preexisting left-censored regions are excluded. End-of-window and unresolved states retain right censoring. Lower/observed/upper clocks preserve cadence-limited times; symlog lifetime axes are linear through one observation cadence. Solid/dashed/dotted distinguish small/intermediate/large; forced conditions are grey. Whole-seed CIs and at-risk/censoring counts remain in CSV. Unidentifiable tails are not extrapolated.")


def render(run_root,heading_root,output,analysis_subdir="analysis"):
    tables,heading,sources=load_audited_tables(run_root,heading_root,analysis_subdir)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"pdf.fonttype":42,"svg.fonttype":"none"})
    book=Supplement(output,sources)
    try:
        metric_pages(book,tables,heading)
        distribution_pages(book,tables,heading)
        survival_pages(book,tables,heading)
        return book.close()
    except Exception:
        book.pdf.close();raise


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root",required=True,type=Path)
    parser.add_argument("--heading-root",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--analysis-subdir",default="analysis")
    arguments=parser.parse_args()
    print(json.dumps(render(arguments.run_root,arguments.heading_root,arguments.output,arguments.analysis_subdir)["pdf"],indent=2))
