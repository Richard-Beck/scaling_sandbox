"""Render audited Figure 3 distance, geometry and direct-heading evidence."""
from pathlib import Path
import argparse
import json
import os
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

import signature_figure3_native_20261020 as native
import signature_figure3_synthesis_20261020 as synthesis

ROOT = native.ROOT
OUT = ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance"
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes 3", spring="Spring")
COLORS = ["#c98a25", "#298895", "#6e548f"]
read, sha, write = native.read, native.sha, native.write


def data():
    source = synthesis.DESTINATION
    generation = read(source/"generation.json")
    assert generation["A2_included"] and generation["B_included"]
    assert read(synthesis.B/"results/all_independent_native_audit.json")["complete"]
    names = dict(A1="A1_complete_region_evidence.csv", B="B_complete_region_evidence_and_geometry.csv",
        A2="A2_all_switching_brackets.csv", remap="A1_fixed_angle_geometry_remap.csv",
        heading="A1_direct_heading_check.csv", targets="A1_heading_target_time_brackets.csv",
        controls="A1_control_handover_convergence.csv", threshold="A1_threshold_handover_sensitivity.csv",
        Bcontrols="B_control_evidence_and_geometry.csv")
    return {key:pd.read_csv(source/name) for key,name in names.items()}, {name:sha(source/name) for name in names.values()}


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, pad=2)


def size_legend(fig, y=.95):
    handles = [Line2D([], [], color=c, lw=2, marker="o", ms=4, label=label) for c,label in zip(COLORS, ["reference size", "1.26× native size label", "1.59× native size label"])]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.52,y), ncol=3, frameon=False, fontsize=9)


def A1_curve(ax, frame, model):
    part = frame[(frame.model == model) & frame.control.isna() & (frame.classification != "internal")]
    for color,size in zip(COLORS, native.SIZES[model]):
        subset = part[part.size_parameter == size].sort_values("gap")
        external = subset[subset.classification == "external"]
        ax.plot(external.gap, external.handover_time, color=color, lw=1.2, marker="o", ms=3)
        borderline = subset[subset.classification == "borderline"]
        ax.scatter(borderline.gap, borderline.handover_time, edgecolors=color, facecolors="none", s=22)
        absent = subset[subset.handover_time.isna()]
        ax.scatter(absent.gap, np.full(len(absent), .035), transform=ax.get_xaxis_transform(), marker="x", color=color, s=16)
        censored = subset[subset.handover_lower_bound.notna()]
        ax.scatter(censored.gap, censored.handover_lower_bound, marker="^", color=color, s=25)
    ax.axvspan(-.13 if model == "spring" else -2., 0, color="#e9e9e9", zorder=-2)
    ax.set_yscale("log")
    ax.set_title(NAMES[model], fontsize=10, weight="bold")
    ax.set_xlabel("Initial edge gap g")
    ax.set_ylabel("Handover time")
    clean(ax)


def B_curve(ax, frame, model, ratio=1., rate=False):
    part = frame[(frame.model == model) & frame.control.isna() & (frame.target_ratio == ratio)]
    for color,size in zip(COLORS, native.SIZES[model]):
        subset = part[(part.size_parameter == size) & part.primary_claim_eligible].sort_values("initial_gap")
        x = subset.measured_gap.to_numpy()
        lo = subset.potential_angular_rate_lower.to_numpy() if rate else subset.resolution_lower.to_numpy()
        hi = subset.potential_angular_rate_upper.to_numpy() if rate else subset.resolution_upper.to_numpy()
        for a,b,c in zip(x,lo,hi):
            if not rate and np.isnan(c): ax.scatter(a,b,marker="^",color=color,s=25)
            elif rate and b == 0 and np.isfinite(c): ax.scatter(a,c,marker="v",color=color,s=25)
            elif rate and b > 0 and np.isnan(c): ax.scatter(a,b,marker="^",color=color,s=25)
            elif b > 0 and np.isfinite(c):
                middle = np.sqrt(b*c)
                ax.errorbar(a,middle,yerr=[[middle-b],[c-middle]],fmt="o",color=color,ms=3,lw=.8,capsize=1)
        # Connect only adjacent eligible, resolved original placements. Filtering
        # first would silently bridge a merged/unseeded/censored placement.
        ordered = part[part.size_parameter == size].sort_values("initial_gap")
        for left,right in zip(ordered.to_dict("records")[:-1], ordered.to_dict("records")[1:]):
            lowkey, highkey = ("potential_angular_rate_lower", "potential_angular_rate_upper") if rate else ("resolution_lower", "resolution_upper")
            if not all(r["primary_claim_eligible"] and np.isfinite(r["measured_gap"]) and r[lowkey] > 0 and np.isfinite(r[highkey]) for r in (left,right)):
                continue
            if right["measured_gap"] > left["measured_gap"]:
                ax.plot([left["measured_gap"],right["measured_gap"]],
                    [np.sqrt(left[lowkey]*left[highkey]),np.sqrt(right[lowkey]*right[highkey])],color=color,lw=.7,alpha=.65)
    ax.set_title(NAMES[model], fontsize=10, weight="bold")
    ax.set_xlabel("Measured withdrawal gap")
    ax.set_ylabel("Potential angular step / time" if rate else "Cue-free resolution time")
    ax.set_yscale("log")
    ax.text(.98,.04,f"{int((~part.primary_claim_eligible).sum())} alternate / unseeded",transform=ax.transAxes,ha="right",va="bottom",fontsize=6.5,
        bbox=dict(facecolor="white",edgecolor="none",alpha=.85,pad=1.))
    clean(ax)


def ring(ax, radius, width, support, C, gap, title):
    angle = np.linspace(0,2*np.pi,400)
    ax.plot(radius*np.cos(angle),radius*np.sin(angle),color="#bfc8ce",lw=3)
    incumbent = np.linspace(-np.pi*width/C,np.pi*width/C,100)
    theta = 2*np.pi*(gap+(width+support)/2)/C
    cue = np.linspace(theta-np.pi*support/C,theta+np.pi*support/C,40)
    ax.plot(radius*np.cos(incumbent),radius*np.sin(incumbent),color="#21689a",lw=6)
    ax.plot(radius*np.cos(cue),radius*np.sin(cue),color="#e58025",lw=6)
    start,end = np.pi*width/C,theta-np.pi*support/C
    arc = np.linspace(start,end,100)
    ax.plot(1.16*radius*np.cos(arc),1.16*radius*np.sin(arc),color="#333333",lw=1.2)
    middle=(start+end)/2
    ax.text(1.40*radius*np.cos(middle),1.40*radius*np.sin(middle),"g",ha="center",va="center",fontsize=12)
    ax.text(0,0,title,ha="center",va="center",fontsize=10)
    ax.set_aspect("equal");ax.set_xlim(-1.65*radius,1.65*radius);ax.set_ylim(-1.6*radius,1.6*radius);ax.axis("off")


def kymographs(axes, A1, model, gaps):
    selected, matrices = [], []
    for gap in gaps:
        row = A1[(A1.model == model) & (A1.size_parameter == 37.8) & A1.control.isna() & (abs(A1.gap-gap)<1e-8)].iloc[0]
        status = read(Path(row.native_status))
        tt,yy=[],[]
        for phase in status["phases"]:
            archive=synthesis.PRIMARY/phase["path"]
            assert sha(archive)==phase["sha256"]
            with np.load(archive) as z:t,y=z["time"].copy(),z["states"][:,:status["placement"]["cells"]].copy()
            use=slice(None) if not tt else slice(1,None)
            tt.extend((phase["absolute_start"]+t[use]).tolist());yy.extend(y[use])
        selected.append((row,status));matrices.append((np.asarray(tt),np.asarray(yy)))
    cap=np.ceil(1.4*max(row.handover_time for row,status in selected)/100)*100
    vmax=max(y[t<=cap].max() for t,y in matrices)
    for ax,(row,status),(t,y) in zip(axes,selected,matrices):
        pl=status["placement"];n=pl["cells"];old=pl["initial_10_percent_contour"]["origin"]/pl["basal_circumference"]
        x=((np.arange(n)+.5)/n-old+.25)%1*37.8;order=np.argsort(x);keep=t<=cap
        image=ax.pcolormesh(x[order],t[keep],y[keep][:,order],shading="nearest",vmin=0,vmax=vmax,cmap="viridis",rasterized=True)
        ax.axhline(row.handover_time,color="white",ls="--",lw=.8)
        ax.invert_yaxis();ax.set_title(f'{NAMES[model]} · g={row.gap:g}',fontsize=9)
        ax.set_xlabel("Native ring arc");ax.set_ylabel("Cue-on time")
        clean(ax)
    return image


def remap_curve(ax, frame, model, heading=None):
    part=frame[frame.model==model]
    for color,size in zip(COLORS,native.SIZES[model]):
        subset=part[part.size_parameter==size].sort_values("theta_degrees")
        ax.plot(subset.theta_degrees,subset.measured_time,color=color,lw=1.4,marker="o",ms=3)
        ax.plot(subset.theta_degrees,subset.geometry_only_time,color=color,lw=1.,ls="--")
    ax.set_yscale("log");ax.set_xlabel("Requested cue-site angle (°)");ax.set_ylabel("Handover time")
    ax.set_title(NAMES[model],fontsize=10,weight="bold");clean(ax)
    if heading is not None:
        inset=ax.inset_axes([.12,.58,.35,.30])
        for color,size in zip(COLORS,native.SIZES[model]):
            q=heading[(heading.model==model)&(heading.size_parameter==size)&heading.control.isna()&(heading.classification=="external")&(heading.coherence_at_handover>=.05)]
            inset.scatter(np.rad2deg(q.cue_axis_target_angle),q.handover_angle_fraction_of_requested,color=color,s=5)
        inset.axhline(1.,color="#888888",lw=.6);inset.axhline(.5,color="#aaaaaa",lw=.6,ls=":")
        inset.set_ylim(0,1.12);inset.set_xticks([0,180]);inset.set_yticks([.5,1]);inset.tick_params(labelsize=6,length=2)
        inset.set_title("Actual angle / target",fontsize=7,pad=2)


def main_plates(tables,out):
    plt.rcParams.update({"font.size":10,"axes.labelsize":9,"xtick.labelsize":8,"ytick.labelsize":8,
        "pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(10.5,13.5));fig.subplots_adjust(top=.955,bottom=.08,left=.075,right=.94)
    gs=fig.add_gridspec(4,12,height_ratios=[1.05,1.35,1.2,1.2],hspace=.78,wspace=1.1)
    fig.suptitle("Figure 3 · Separation affects native competition differently across models",fontsize=13,y=.995)
    ref=tables["A1"][(tables["A1"].model=="goryachev")&(tables["A1"].size_parameter==37.8)].iloc[0]
    ring(fig.add_subplot(gs[0,:3]),1.,ref.initial_10_percent_width,3.,37.8,2.,"Fixed C\nnear")
    ring(fig.add_subplot(gs[0,3:6]),1.,ref.initial_10_percent_width,3.,37.8,8.,"Fixed C\nfar")
    ax=fig.add_subplot(gs[0,6:]);ax.axis("off")
    ax.text(0,1,"A  Rotate one native cue at fixed size",weight="bold",fontsize=11,va="top")
    ax.text(0,.78,"g: initial incumbent 10% edge → cue near edge\nA1: continuous cue → complete-region handover\nA2: own-handover-scaled pulse → cue-free fate\nB: seed ratio 0.5 / 1 / 2 → cue-free peak resolution\n12 fresh qualified basals; native kinetics unchanged",fontsize=9.2,va="top",linespacing=1.65)
    kg=gs[1,:].subgridspec(1,8,width_ratios=[1,1,.045,.60,1,1,.045,.15],wspace=.35)
    axes=[fig.add_subplot(kg[0,i]) for i in [0,1,4,5]]
    image=kymographs(axes[:2],tables["A1"],"goryachev",[2.,8.]);fig.colorbar(image,cax=fig.add_subplot(kg[0,2]),label="Native active signal")
    image=kymographs(axes[2:],tables["A1"],"holmes_model3",[2.,8.]);fig.colorbar(image,cax=fig.add_subplot(kg[0,6]),label="Active Rac")
    axes[1].set_ylabel("");axes[3].set_ylabel("")
    fig.text(.04,.756,"B  Same-size native histories · white line: complete-region handover",fontsize=11,weight="bold")
    for i,model in enumerate(NAMES): A1_curve(fig.add_subplot(gs[2,i*3:(i+1)*3]),tables["A1"],model)
    for i,model in enumerate(NAMES): B_curve(fig.add_subplot(gs[3,i*3:(i+1)*3]),tables["B"],model)
    fig.text(.04,.491,"C1  Driven handover · open symbols: borderline placements",fontsize=11,weight="bold")
    fig.text(.04,.267,"C2  Cue-free resolution · ratio 1 · triangles retain finite-cap lower bounds",fontsize=11,weight="bold")
    fig.legend(handles=[Line2D([],[],color=c,marker="o",lw=1,label=f"size label ×{ratio}") for c,ratio in zip(COLORS,["1","1.26","1.59"])],loc="lower center",bbox_to_anchor=(.5,.018),ncol=3,frameon=False,fontsize=8)
    fig.savefig(out/"figure3_distance_main.png",dpi=220)
    fig.savefig(out/"supporting/vectors/figure3_distance_main.svg")
    second=plt.figure(figsize=(10.5,10.5));grid=second.add_gridspec(3,4,height_ratios=[1.,1.3,1.3],hspace=.7,wspace=.9)
    second.suptitle("Figure 3 · Geometry remaps challenges; activity direction is measured separately",fontsize=13,y=.99)
    for col,C in enumerate([30.,47.6]):
        base=tables["A1"][(tables["A1"].model=="wave_pinning")&(tables["A1"].size_parameter==C)].iloc[0]
        ax=second.add_subplot(grid[0,col])
        ring(ax,C/30.,base.initial_10_percent_width,3.,C,4.,f"C={C:g}")
        ax.set_xlim(-2.6,2.6);ax.set_ylim(-2.6,2.6)
    ax=second.add_subplot(grid[0,2:]);ax.axis("off")
    ax.text(0,1,"D  Physical gap and center angle",fontsize=11,weight="bold",va="top")
    ax.text(0,.74,r"$\theta=\frac{2\pi}{C}\left[g+\frac{w_{10}+s}{2}\right]$",fontsize=17,va="top")
    ax.text(0,.29,"Measured widths enter the center step.\nHandover and full-angle rotation are separate events.",fontsize=9,va="top",linespacing=1.5)
    for i,model in enumerate(["wave_pinning","goryachev"]):
        remap_curve(second.add_subplot(grid[1,i*2:(i+1)*2]),tables["remap"],model,tables["heading"])
        B_curve(second.add_subplot(grid[2,i*2:(i+1)*2]),tables["B"],model,rate=True)
    second.text(.05,.685,"E  Fixed-angle cue placement: measured handover (solid), small-cell remap (dashed)",fontsize=10.5,weight="bold")
    second.text(.05,.385,"F  Fixed physical scale: conditional geometric step / resolution time",fontsize=10.5,weight="bold")
    second.text(.5,.025,"Dashed remaps stop at observed small-cell support. Insets measure actual global angle at handover, R≥0.05.\n"
        "Potential step uses measured withdrawal centers and widths; direct heading-time tests and spatial outcomes are retained in the supplement.\n"
        "These are native polarity-activity readouts; observed migration and conditional OU persistence estimates remain separate.",ha="center",fontsize=8)
    second.subplots_adjust(top=.93,bottom=.13,left=.09,right=.96)
    size_legend(second,.955)
    second.savefig(out/"figure3_distance_translation.png",dpi=220)
    second.savefig(out/"supporting/vectors/figure3_distance_translation.svg")
    with PdfPages(out/"supporting/vectors/figure3_distance_main.pdf") as pdf:
        pdf.savefig(fig);pdf.savefig(second)
    plt.close(fig);plt.close(second)


def supplements(tables,out):
    destination=out/"supporting/vectors/figure3_distance_supplement.pdf"
    with PdfPages(destination) as pdf:
        # Every intended placement remains visible, including unavailable own
        # handover calibrations. Amount and field-stability gates stay separate.
        fig,axes=plt.subplots(4,3,figsize=(10.5,13.5))
        palette=dict(challenger="#e58025",incumbent="#21689a",unresolved="#888888",location_ambiguous="#aaaaaa")
        lookup=tables["A1"][tables["A1"].control.isna()].set_index("placement_id")
        for i,model in enumerate(NAMES):
            for j,size in enumerate(native.SIZES[model]):
                ax=axes[i,j];part=tables["A2"][(tables["A2"].model==model)&(tables["A2"].size_parameter==size)]
                for record in part.to_dict("records"):
                    gap=lookup.loc[record["placement_id"],"gap"]
                    if pd.isna(record.get("tested_trials")):
                        ax.scatter(np.full(5,gap),[.5,.75,1.,1.5,2.],marker="x",color="#aaaaaa",s=14)
                        continue
                    for trial in json.loads(record["tested_trials"]):
                        color=palette.get(trial["winner"],"#888888")
                        ax.scatter(gap,trial["multiplier"],marker="s",s=24,color=color)
                        if not trial["finite_tail_pass"]:ax.scatter(gap,trial["multiplier"],marker="x",s=12,color="black")
                ax.axvspan(-.13 if model=="spring" else -2,0,color="#eeeeee",zorder=-2)
                ax.set_title(f'{NAMES[model]} · size {size:g}',fontsize=9)
                ax.set_yticks([.5,1.,1.5,2.]);ax.set_xlabel("Initial edge gap");ax.set_ylabel("Pulse / own handover");clean(ax)
        fig.suptitle("S1 · All 545 cue-withdrawal pulses and 34 unavailable own-handover placements",fontsize=12)
        fig.text(.5,.018,"Blue: incumbent-associated final region; orange: challenger-associated; grey: unresolved. Grey ×: no true two-region handover.\nBlack × over a tested cell: finite native-field tail fails, independent of spatial fate. Shading: borderline placement.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.965));pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(4,3,figsize=(10.5,13.5))
        for i,model in enumerate(NAMES):
            for j,ratio in enumerate([.5,1.,2.]):
                B_curve(axes[i,j],tables["B"],model,ratio)
                axes[i,j].set_title(f'{NAMES[model]} · ratio {ratio:g}',fontsize=9)
        fig.suptitle("S2 · Cue-free resolution at all three complete-region withdrawal ratios",fontsize=12)
        size_legend(fig,.965)
        fig.text(.5,.012,"Primary gap: current global10%-contrast contours. ▲: resolution-time lower bound at completed horizon.\nAlternate, merged and unseeded cases remain in the complete evidence table and native-history reader.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.94));pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(10.5,9.5))
        for ax,model in zip(axes.flat,NAMES):remap_curve(ax,tables["remap"],model,tables["heading"])
        fig.suptitle("S3 · Fixed-angle handover remaps in all four models",fontsize=12)
        size_legend(fig,.955)
        fig.text(.5,.012,"Solid: measured native size; dashed: reference-size handover curve evaluated at the larger cell's measured-width gap.\nNo extrapolation or interpolation through unavailable/censored neighbors. Insets: actual direction at handover / requested cue angle.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.925));pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(4,2,figsize=(10.5,13.5))
        for i,model in enumerate(NAMES):
            ax,clock=axes[i]
            for color,size in zip(COLORS,native.SIZES[model]):
                part=tables["heading"][(tables["heading"].model==model)&(tables["heading"].size_parameter==size)&tables["heading"].control.isna()&(tables["heading"].classification=="external")]
                visible=part[part.coherence_at_handover>=.05]
                ax.scatter(np.rad2deg(visible.cue_axis_target_angle),np.rad2deg(visible.angle_at_handover),color=color,s=14)
                target=tables["targets"][(tables["targets"].model==model)&(tables["targets"].size_parameter==size)&tables["targets"].control.isna()&(tables["targets"].classification=="external")&(tables["targets"].target_fraction==1.)&(tables["targets"].coherence_floor==.05)]
                target=target.merge(part[["slot","native_handover_time"]],on="slot")
                for r in target.to_dict("records"):
                    if not np.isfinite(r["native_handover_time"]):continue
                    if r["status"]=="observed_saved_valid_segment_crossing":clock.scatter(r["native_handover_time"],r["upper"],color=color,s=14,marker="D" if r.get("prior_heading_gap") else "o")
                    elif r["status"]=="right_censored_saved_heading_at_completed_horizon":clock.scatter(r["native_handover_time"],r["lower"],marker="^",color=color,s=22)
                    else:clock.scatter(r["native_handover_time"],.04,transform=clock.get_xaxis_transform(),marker="x",color=color,s=14)
            ax.plot([0,180],[0,180],color="#888888",lw=.8);ax.plot([0,180],[0,90],color="#aaaaaa",lw=.8,ls=":")
            ax.set_xlim(0,180);ax.set_ylim(0,180);ax.set_title(NAMES[model],fontsize=10,weight="bold")
            ax.set_xlabel("Cue-site target angle (°)");ax.set_ylabel("Actual angle at handover (°)")
            clock.set_xscale("log");clock.set_yscale("log");clock.set_xlabel("Complete-region handover time");clock.set_ylabel("Full-angle target time")
            clean(ax);clean(clock)
        fig.suptitle("S4 · Handover typically reaches part of the requested angle; full-target timing is separate",fontsize=11)
        size_legend(fig,.965)
        fig.text(.5,.012,"R=|M₁|/M₀≥0.05 at usable samples and intervening conservative complex-moment chords. ▲: finite-cap timing lower bound.\n×: direction-quality/cancellation gap prevents a full timing bound. Exact-direction and0.1/0.2 gates plus50/90/99% targets are tabulated.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.94));pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(10.5,9.5))
        for ax,model in zip(axes.flat,NAMES):
            for color,size in zip(COLORS,native.SIZES[model]):
                part=tables["B"][(tables["B"].model==model)&(tables["B"].size_parameter==size)&tables["B"].control.isna()&tables["B"].primary_claim_eligible&tables["B"].resolution_upper.notna()]
                for winner,marker in [("challenger","o"),("incumbent","s")]:
                    selected=part[(part.winner==winner)&part.vector_turn_basal_to_resolution_or_cap.notna()]
                    ax.scatter(np.rad2deg(selected.potential_angular_step),np.rad2deg(abs(selected.vector_turn_basal_to_resolution_or_cap)),marker=marker,color=color,s=14,alpha=.7)
            ax.plot([0,180],[0,180],lw=.8,color="#888888");ax.set_xlim(0,180);ax.set_ylim(0,180)
            ax.set_title(NAMES[model],weight="bold",fontsize=10);ax.set_xlabel("Potential region-center step (°)");ax.set_ylabel("Actual net global-heading change (°)");clean(ax)
        fig.suptitle("S5 · Physical center steps and native global-heading changes can differ",fontsize=12)
        size_legend(fig,.955)
        fig.text(.5,.012,"Resolved primary B cases; basals and resolution endpoints require R≥0.05. ○ challenger-associated, □ incumbent-associated final region.\nActivity can move or coalesce to an intermediate direction. Endpoint angles are not an assumed cell-velocity law or an unwrapped path length.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.925));pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(4,3,figsize=(10.5,13.5))
        for i,model in enumerate(NAMES):
            for j,size in enumerate(native.SIZES[model]):
                ax=axes[i,j]
                base=tables["A1"][(tables["A1"].model==model)&(tables["A1"].size_parameter==size)&tables["A1"].control.isna()&(tables["A1"].classification=="external")].sort_values("gap")
                ax.plot(base.gap,base.handover_time,color="black",lw=1.2,label="Full region")
                for color,fraction in zip(["#4e8e79","#ba8649","#9c6b9d"],[.25,.5,.75]):
                    part=tables["threshold"][(tables["threshold"].model==model)&(tables["threshold"].size_parameter==size)&tables["threshold"].control.isna()&(tables["threshold"].classification=="external")&(tables["threshold"].contour_fraction==fraction)].sort_values("gap")
                    ax.plot(part.gap,part.upper,color=color,marker="o",ms=2,lw=.7,label=f"{fraction:g}")
                    absent=part[part.upper.isna()]
                    ax.scatter(absent.gap,np.full(len(absent),.035),transform=ax.get_xaxis_transform(),marker="x",s=12,color=color)
                ax.set_yscale("log");ax.set_title(f'{NAMES[model]} · {size:g}',fontsize=9);ax.set_xlabel("Initial edge gap");ax.set_ylabel("Handover time");clean(ax)
        axes[0,0].legend(fontsize=7)
        fig.suptitle("S6 · Amount sensitivity: local 25/50/75% peak contours versus complete-region excess",fontsize=11)
        fig.text(.5,.012,"Original excess above each region's local-peak contour; watershed basins and common current floor are unchanged.\nTimes are retained native-sample upper brackets. ×: unavailable continuous sampled two-region crossing; no cap substitution.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.96));pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(10.5,9.5))
        for ax,model in zip(axes[0],["wave_pinning","goryachev"]):
            for color,variant in [(COLORS[0],"fine_mesh"),(COLORS[1],"half_max_step_tighter_tolerances")]:
                part=tables["controls"][(tables["controls"].model==model)&(tables["controls"].variant==variant)]
                ax.plot(part.gap,100*part.relative_time_difference,color=color,marker="o",label=variant.replace("_"," "))
            ax.set_title(NAMES[model],fontsize=10,weight="bold");ax.set_xlabel("Matched initial gap");ax.set_ylabel("Relative handover difference (%)");clean(ax)
        axes[0,0].legend(fontsize=7)
        for ax in axes[1]:ax.axis("off")
        axes[1,0].text(0,1,"Numerical evidence",fontsize=12,weight="bold",va="top")
        axes[1,0].text(0,.85,"Fresh native basals: 12 / 12 qualified\nLong no-cue controls: 53,000 chemical / 5,300 Spring\nFine-mesh basals: 2 / 2 independently qualified\nA1 refinement controls: 12 / 12 native audits\nObserved control roots: 10 / 10 verified\nA2 native audit: 545 / 545\nB native audit: 381 / 381 (24 mirrors agree)\nB mesh/step controls: 36 / 36 native audits",fontsize=9.5,va="top",linespacing=1.6)
        axes[1,1].text(0,1,"Readout and source gates",fontsize=12,weight="bold",va="top")
        axes[1,1].text(0,.85,"Native positivity and conserved-pool gates pass.\nScalar fate and field stationarity remain separate.\nContinuous roots replay the original upper fields.\nB controls retain all seed statuses and winners.\n35 / 36 resolution brackets overlap.\nGoryachev g=4, ratio 1, fine mesh: 1040–1060\nversus 1080–1100 (about 4% timing sensitivity).\nFull native atlas: 1,150 histories; all sources linked.",fontsize=9.5,va="top",linespacing=1.6)
        fig.suptitle("S7 · Numerical, preparation and complete-source checks",fontsize=12)
        fig.tight_layout(rect=(0,.035,1,.96));pdf.savefig(fig);plt.close(fig)
    return destination


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=OUT);args=parser.parse_args()
    out=args.output.resolve()
    for category in ("vectors","metadata","data","viewers"):(out/"supporting"/category).mkdir(parents=True,exist_ok=True)
    tables,hashes=data()
    for name in hashes:shutil.copy2(synthesis.DESTINATION/name,out/"supporting/data"/name)
    main_plates(tables,out)
    supplement=supplements(tables,out)
    write(out/"supporting/metadata/main_generation.json",dict(generated_at=native.now(),native_integrations=0,
        source_data_sha256=hashes,renderer_sha256=sha(Path(__file__)),main_pages=2,
        outputs={str(p.relative_to(out)):sha(p) for p in [out/"figure3_distance_main.png",out/"figure3_distance_translation.png",supplement,
            *list((out/"supporting/vectors").glob("figure3_distance_main.*")),*list((out/"supporting/vectors").glob("figure3_distance_translation.*"))]}))
    print(out)


if __name__=="__main__":main()
