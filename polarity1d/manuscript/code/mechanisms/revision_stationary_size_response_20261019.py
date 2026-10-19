"""Read-only three-size stationary-cue response atlases from authoritative tables."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import html
import io
import json
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

REPO = Path(__file__).resolve().parents[4]
TOPIC = REPO / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/standard_angles_sizes"
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring")
NAMES = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "holmes_model3":"Holmes 3", "spring":"Spring ring"}
ANGLES = ("pi_8", "pi_4", "pi_2", "pi")
DEGREES = dict(pi_8=22.5, pi_4=45., pi_2=90., pi=180.)
SIZECOLORS = ("#397bab", "#40916c", "#cf7930")


def sha_bytes(data):return hashlib.sha256(data).hexdigest()
def sizes(model):return (1.,2.,3.) if model=="spring" else (20.,30.,40.)
def durations(model):return (.1,1.,10.,100.,1000.) if model=="spring" else (1.,10.,100.,300.,1000.)
def size_label(model,size):return ("C₀=" if model=="spring" else "C=")+f"{size:g}"
def number(row,key):
    try:return float(row.get(key,np.nan))
    except (TypeError,ValueError):return np.nan
def finite(value):return bool(np.isfinite(value))
def wrap(deg):return float(np.rad2deg(np.angle(np.exp(1j*np.deg2rad(deg)))))


def load(summary,manifest):
    with summary.open("rb") as stream:
        raw=stream.read();source_stat=os.fstat(stream.fileno())
    frame=pd.read_csv(io.BytesIO(raw)).copy();frame.attrs["source_snapshot_bytes"]=raw
    if "size_parameter" not in frame:
        frame["size_parameter"]=np.where(frame.model=="spring",pd.to_numeric(frame.resting_circumference),pd.to_numeric(frame.circumference))
    assert not frame.duplicated(["model","size_parameter","angle_label","cue_duration"]).any()
    expected={(m,s,a,t) for m in MODELS for s in sizes(m) for a in ANGLES for t in durations(m)}
    provenance=dict(summary_path=str(summary),summary_sha256=sha_bytes(raw),source_snapshot_modified_at=datetime.fromtimestamp(source_stat.st_mtime,timezone.utc).isoformat(),expected_grid="explicit user-requested sizes and inherited stationary angle/duration set")
    if manifest is not None:
        blob=manifest.read_bytes();frozen=json.loads(blob)
        keys=set()
        for c in frozen["cases"]:
            m=c["model"]
            size=c.get("size_parameter",c.get("resting_circumference") if m=="spring" else c.get("circumference"))
            keys.add((m,float(size),c["angle_label"],float(c["cue_duration"])))
        assert keys==expected,(len(keys),len(expected),keys-expected,expected-keys)
        provenance.update(manifest_path=str(manifest),manifest_sha256=sha_bytes(blob),expected_grid="frozen native manifest, checked against requested three-size grid")
    lookup={}
    for _,r in frame.iterrows():
        key=(r.model,float(r.size_parameter),r.angle_label,float(r.cue_duration))
        assert key in expected,key
        complete=str(r.get("status","complete"))=="complete" and finite(number(r,"R_final")) and finite(number(r,"M0_ratio_final"))
        lookup[key]=(r.to_dict(),complete)
    return frame,lookup,expected,provenance


def save(fig,stem,out):
    paths=[]
    for ext in ("png","pdf","svg"):
        p=(out if ext=="png" else out/"supporting/vectors")/(stem+"."+ext)
        fig.savefig(p,dpi=220,bbox_inches="tight");paths.append(p)
    plt.close(fig);return paths


def footer(fig,complete,generated):
    state="complete grid" if complete==240 else "pending conditions remain explicit"
    fig.text(.04,.017,f"Analyzed histories {complete}/240 · generated {generated} · {state} · native clocks differ across models",fontsize=7.5)


def matrix(ax,model,size,lookup,metric,title,cmap,vmin,vmax,heading=False,phase=None):
    ts=durations(model);z=np.full((4,5),np.nan);pending=np.ones((4,5),bool);weak=np.zeros((4,5),bool);rows={}
    for i,a in enumerate(ANGLES):
        for j,t in enumerate(ts):
            row,done=lookup.get((model,size,a,t),({},False));rows[i,j]=row
            if done:
                pending[i,j]=False
                weak[i,j]=heading and number(row,"R_"+phase)<.1
                if not weak[i,j]:z[i,j]=number(row,metric)
    cm=plt.get_cmap(cmap).copy();cm.set_bad("#eeeeee")
    image=ax.imshow(z,origin="lower",aspect="auto",cmap=cm,vmin=vmin,vmax=vmax,extent=(-.5,4.5,-.5,3.5))
    for i in range(4):
        for j in range(5):
            if pending[i,j]:
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor="white",edgecolor=".8",hatch="///",lw=.3))
                ax.text(j,i,"pending",ha="center",va="center",fontsize=6.3,color=".45")
            elif weak[i,j]:ax.text(j,i,"R<.1",ha="center",va="center",fontsize=7,color=".35")
            elif finite(z[i,j]):
                normalized=(z[i,j]-vmin)/(vmax-vmin) if vmax>vmin else .5
                rgba=cm(np.clip(normalized,0,1));luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                ink="white" if luminance<.45 else ".15"
                displayed="0" if abs(z[i,j])<1e-10 else (f"{z[i,j]:.1f}" if heading else f"{z[i,j]:.3g}")
                ax.text(j,i,displayed,ha="center",va="center",fontsize=7,color=ink)
            row=rows[i,j]
            if not pending[i,j] and phase=="final" and row.get("full_field_tail_status")=="fails_finite_tail_gate":
                ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,edgecolor="black",lw=1.5))
    ax.set_xticks(range(5),[f"{t:g}" for t in ts]);ax.set_yticks(range(4),[f"{DEGREES[a]:g}°" for a in ANGLES])
    ax.set_xlabel("Tested cue duration (native units)",fontsize=8);ax.set_ylabel("Prescribed cue angle",fontsize=8)
    ax.set_title(title,loc="left",fontsize=9.5)
    return image


def response_atlas(lookup,out,complete,generated):
    paths=[]
    for model in MODELS:
        fig,axes=plt.subplots(2,3,figsize=(11.8,6.7))
        for j,s in enumerate(sizes(model)):
            for i,(phase,field,label) in enumerate((("withdrawal","angle_withdrawal_deg","At withdrawal"),("final","angle_final_deg","After 13,000 cue-free units"))):
                im=matrix(axes[i,j],model,s,lookup,field,f"{size_label(model,s)} · {label}","twilight_shifted",-180,180,heading=True,phase=phase)
        fig.suptitle(f"{NAMES[model]} · all four stationary-cue angles across three sizes",fontsize=14,y=.99)
        fig.text(.04,.925,"Measured activity heading in P=M₁/M₁(0) coordinates. Grey R<.1 cells retain a valid weak vector but no interpreted heading.",fontsize=8.4)
        fig.text(.04,.897,"Durations are tested categories; cell spacing does not imply equal native time. Black outline: separate final full-field tail gate failed.",fontsize=8)
        fig.subplots_adjust(left=.075,right=.895,bottom=.10,top=.83,wspace=.28,hspace=.48)
        cb=fig.colorbar(im,cax=fig.add_axes((.92,.22,.014,.48)),ticks=(-180,-90,0,90,180));cb.set_label("Activity heading relative to initial M₁ (°)",fontsize=8)
        footer(fig,complete,generated);paths+=save(fig,"stationary_sizes_angles_"+model,out)
    return paths


def all_angles_overview(lookup,out,complete,generated):
    """One primary review composition; preserve every angle and native clock."""
    fig,axes=plt.subplots(4,4,figsize=(13.2,11.4))
    limits={}
    for a in ANGLES:
        values=[number(row,"angle_"+phase+"_deg") for key,(row,done) in lookup.items() if done and key[2]==a for phase in ("withdrawal","final") if number(row,"R_"+phase)>=.1 and finite(number(row,"angle_"+phase+"_deg"))]
        target=DEGREES[a];low=min(values+[0.]);high=max(values+[target]);pad=.055*max(target,high-low)
        limits[a]=(low-pad,high+pad)
    for i,model in enumerate(MODELS):
        for j,a in enumerate(ANGLES):
            ax=axes[i,j];count=0
            for size,color in zip(sizes(model),SIZECOLORS):
                x=np.asarray(durations(model));ys={p:[] for p in ("withdrawal","final")};weak={p:[] for p in ys};targets=[]
                for t in x:
                    row,done=lookup.get((model,size,a,float(t)),({},False));count+=int(done)
                    if done:targets.append(number(row,"target_angle_in_P_coordinates_deg"))
                    for phase in ys:
                        valid=done and number(row,"R_"+phase)>=.1
                        ys[phase].append(number(row,"angle_"+phase+"_deg") if valid else np.nan)
                        weak[phase].append(done and not valid)
                for phase,style,fill,marker in (("withdrawal","--","white","+"),("final","-",color,"x")):
                    ax.plot(x,ys[phase],color=color,ls=style,marker="o",ms=3.3,mfc=fill,lw=1.1,mew=.8)
                    selected=x[np.asarray(weak[phase])]
                    if len(selected):ax.plot(selected,np.full(len(selected),.035 if phase=="final" else .075),marker=marker,ls="none",color=color,ms=4,transform=ax.get_xaxis_transform(),clip_on=True)
                for k,t in enumerate(x):
                    row,done=lookup.get((model,size,a,float(t)),({},False))
                    if done and row.get("full_field_tail_status")=="fails_finite_tail_gate":
                        if number(row,"R_final")>=.1:ax.scatter(t,ys["final"][k],s=43,marker="s",facecolors="none",edgecolors="black",lw=.7,zorder=5)
                        else:ax.plot(t,.035,marker="s",ms=6,mfc="none",mec="black",mew=.7,ls="none",transform=ax.get_xaxis_transform(),zorder=5)
                if any(finite(t) for t in targets):ax.axhline(next(t for t in targets if finite(t)),color=".65",ls=":",lw=.5,zorder=0)
            ax.set_xscale("log");ax.set_xlim(min(durations(model))*.7,max(durations(model))*1.4);ax.set_ylim(*limits[a])
            ax.set_title(f"{NAMES[model]} · {DEGREES[a]:g}°",loc="left",fontsize=8.8,fontweight="bold")
            if complete<240:ax.text(.98,.96,f"n={count}/15",ha="right",va="top",transform=ax.transAxes,color=".4",fontsize=6.7)
            ax.set_xlabel("Cue duration (native units)",fontsize=7.5);ax.set_ylabel("Activity heading (°)",fontsize=7.5)
            ax.tick_params(labelsize=7)
    fig.suptitle("Three sizes · all four stationary-cue angles",fontsize=15,y=.995)
    fig.text(.5,.961,"Open/dashed: cue withdrawal · filled/solid: after 13,000 cue-free units · actual P-relative angles with a shared range in each column",ha="center",fontsize=8.7)
    handles=[Line2D([],[],color=c,label=l) for c,l in zip(SIZECOLORS,("Small: C20 / Spring C₀1","Middle: C30 / Spring C₀2","Large: C40 / Spring C₀3"))]
    handles += [Line2D([],[],color=".3",ls=s,marker="o",mfc=m,label=l) for s,m,l in (("--","white","At withdrawal"),("-",".3","After withdrawal"))]
    fig.legend(handles=handles,ncol=3,frameon=False,loc="lower center",bbox_to_anchor=(.5,.043),fontsize=8)
    fig.text(.5,.031,"Bottom +/×: weak-R withdrawal/final, no heading · □: final full-field tail failure · Five tested durations per curve; connecting segments guide the eye.",ha="center",fontsize=7.3)
    fig.tight_layout(rect=(0,.10,1,.94),h_pad=2.0,w_pad=2)
    footer(fig,complete,generated)
    return save(fig,"stationary_sizes_all_angles_overview",out)


def context_table(frame,out):
    records=[]
    fields=("circumference","resting_circumference","original_core_physical_width","full_physical_support","initial_angular_support_deg","physical_distance_at_onset","cue_centre_physical_at_onset","preparation_axis_fraction","cue_centre_fraction","requested_peak","initial_sampled_cue_peak","initial_cue_integral","overlap_fraction_cue_support","overlap_fraction_original_core","initial_unit_peak_cue_excess_fraction_total","M0_initial","M1_initial_real","M1_initial_imag","initial_M1_absolute_phase_deg","initial_M1_offset_from_preparation_deg")
    for model in MODELS:
        for size in sizes(model):
            for angle in ANGLES:
                group=frame[(frame.model==model)&(frame.size_parameter==size)&(frame.angle_label==angle)]
                record=dict(model=model,size_parameter=size,size_label=size_label(model,size),angle_label=angle,requested_angle_deg=DEGREES[angle],source_rows=len(group),context_available=bool(len(group)))
                for field in fields:
                    vals=pd.to_numeric(group[field],errors="coerce").dropna().values if field in group else np.array([])
                    if not len(vals) and field in ("M0_initial","M1_initial_real","M1_initial_imag"):
                        basal=frame[(frame.model==model)&(frame.size_parameter==size)]
                        if "input_sha256" in basal:assert len(basal.input_sha256.dropna().unique())<=1
                        vals=pd.to_numeric(basal[field],errors="coerce").dropna().values if field in basal else np.array([])
                    if len(vals):
                        assert np.allclose(vals,vals[0],rtol=1e-10,atol=1e-12),(model,size,angle,field,vals)
                        record[field]=float(vals[0])
                    else:record[field]=np.nan
                m1=np.hypot(record["M1_initial_real"],record["M1_initial_imag"])
                record.update(initial_abs_M1=m1,initial_q1=m1/record["M0_initial"] if record["M0_initial"]>0 else np.nan)
                record["initial_core_angular_width_deg"]=360*record["original_core_physical_width"]/record["circumference"]
                records.append(record)
    context=pd.DataFrame(records);context.to_csv(out/"supporting/data/stationary_sizes_basal_cue_context.csv",index=False)
    return context


def paired_size_contrasts(lookup,out):
    """Descriptive matched comparisons; no inferred threshold or uncertainty."""
    result=[]
    for model in MODELS:
        ss=sizes(model)
        for left,right in ((ss[0],ss[1]),(ss[1],ss[2]),(ss[0],ss[2])):
            for angle in ANGLES:
                for t in durations(model):
                    a,done_a=lookup.get((model,left,angle,t),({},False));b,done_b=lookup.get((model,right,angle,t),({},False))
                    row=dict(model=model,smaller_size=left,larger_size=right,angle_label=angle,cue_angle_deg=DEGREES[angle],cue_duration=t,matched_pair_available=bool(done_a and done_b))
                    for phase in ("withdrawal","final"):
                        valid=done_a and done_b and min(number(a,"R_"+phase),number(b,"R_"+phase))>=.1
                        row[phase+"_heading_pair_defined"]=bool(valid)
                        row[phase+"_larger_minus_smaller_circular_angle_deg"]=wrap(number(b,"angle_"+phase+"_deg")-number(a,"angle_"+phase+"_deg")) if valid else np.nan
                    for field in ("R_final","R_withdrawal","M0_ratio_final","M0_ratio_withdrawal","fixed_region_excess_final","fixed_region_rebound_from_sampled_minimum","initial_angular_support_deg","overlap_fraction_cue_support"):
                        row["larger_minus_smaller_"+field]=number(b,field)-number(a,field) if done_a and done_b else np.nan
                    result.append(row)
    pd.DataFrame(result).to_csv(out/"supporting/data/stationary_sizes_matched_contrasts.csv",index=False)


def plot_context(context,out,complete,generated):
    fig,axes=plt.subplots(4,2,figsize=(10.8,11.5))
    for i,model in enumerate(MODELS):
        ax=axes[i,0]
        g=context[(context.model==model)&(context.angle_label=="pi")]
        ax.plot(g.size_parameter,g.initial_angular_support_deg,"o-",color="#4a6272",lw=1.2,ms=4)
        ax.set_xticks(sizes(model));ax.set_xlabel("Resting C₀" if model=="spring" else "Circumference C")
        ax.set_ylabel("Initial full cue width (°)");ax.set_title(NAMES[model]+" · native footprint",loc="left",fontsize=10)
        z=np.full((3,4),np.nan)
        for j,s in enumerate(sizes(model)):
            for k,a in enumerate(ANGLES):
                r=context[(context.model==model)&(context.size_parameter==s)&(context.angle_label==a)].iloc[0];z[j,k]=r.overlap_fraction_cue_support
        ax=axes[i,1];im=ax.imshow(z,cmap="Blues",vmin=0,vmax=1,origin="lower",aspect="auto")
        for j in range(3):
            for k in range(4):ax.text(k,j,f"{z[j,k]:.1%}" if finite(z[j,k]) else "pending",ha="center",va="center",fontsize=8,color="white" if z[j,k]>.65 else ".25")
        ax.set_xticks(range(4),[f"{DEGREES[a]:g}°" for a in ANGLES]);ax.set_yticks(range(3),[size_label(model,s) for s in sizes(model)])
        ax.set_xlabel("Prescribed cue angle");ax.set_title("Fraction of cue support overlapping initial core",loc="left",fontsize=9)
    fig.suptitle("Size context · fixed physical support changes angular exposure",fontsize=14,y=.99)
    fig.text(.05,.95,"Linked context table retains basal physical M₀, |M₁|, q₁, original-core width, actual C and cue peak/integral for every size/angle.",fontsize=8.3)
    fig.tight_layout(rect=(0,.05,1,.925),h_pad=2.5,w_pad=3)
    footer(fig,complete,generated);return save(fig,"stationary_sizes_basal_cue_context",out)


def local_and_amount(lookup,out,complete,generated):
    paths=[]
    metrics=(("R_final","Final relative resultant R","viridis",0,1.2),
             ("M0_ratio_final","Final total active amount M₀/M₀(0)","viridis",None,None),
             ("fixed_region_excess_cue_on_sampled_minimum","Cue-on sampled original-region minimum","viridis",None,None),
             ("fixed_region_rebound_from_sampled_minimum","Final regional excess minus sampled minimum","viridis",None,None))
    for model in MODELS:
        fig,axes=plt.subplots(4,3,figsize=(11.7,11.8))
        for i,(field,label,cmap,lower,upper) in enumerate(metrics):
            v=[number(row,field) for key,(row,done) in lookup.items() if key[0]==model and done and finite(number(row,field))]
            vmin=lower if lower is not None else min(v,default=0)
            vmax=max(upper,max(v,default=upper)) if upper is not None else max(v,default=1)
            if vmax-vmin<1e-10:vmin-=.01;vmax+=.01
            for j,s in enumerate(sizes(model)):
                im=matrix(axes[i,j],model,s,lookup,field,size_label(model,s),cmap,vmin,vmax)
            cb=fig.colorbar(im,ax=list(axes[i,:]),fraction=.015,pad=.015);cb.ax.tick_params(labelsize=7)
        fig.suptitle(NAMES[model]+" · magnitude, amount and fixed-region response",fontsize=14,y=.99)
        fig.text(.04,.955,"All angles and sizes, same native histories. Extrema are sampled values; envelopes/cells are not independent-seed uncertainty.",fontsize=8.2)
        fig.text(.04,.93,"The original basal region is a fixed material interval. Depletion/rebound can reflect translation or redistribution and does not establish site lineage.",fontsize=8)
        fig.subplots_adjust(left=.075,right=.88,bottom=.085,top=.88,wspace=.28,hspace=.63)
        for i,(_,label,_,_,_) in enumerate(metrics):
            fig.text(.075,axes[i,0].get_position().y1+.034,label,fontsize=9.5,fontweight="bold")
        footer(fig,complete,generated);paths+=save(fig,"stationary_sizes_pattern_amount_"+model,out)
    return paths


def candidate(lookup,out,complete,generated):
    fig,axes=plt.subplots(4,2,figsize=(10.8,11.5))
    for i,model in enumerate(MODELS):
        for j,angle in enumerate(("pi_2","pi")):
            ax=axes[i,j]
            for s,color in zip(sizes(model),SIZECOLORS):
                x=np.asarray(durations(model));ys={p:[] for p in ("withdrawal","final")}
                for t in x:
                    r,done=lookup.get((model,s,angle,float(t)),({},False))
                    for p in ys:ys[p].append(number(r,"angle_"+p+"_deg") if done and number(r,"R_"+p)>=.1 else np.nan)
                for p,style,fill in (("withdrawal","--","white"),("final","-",color)):
                    ax.plot(x,ys[p],color=color,ls=style,marker="o",ms=3.5,mfc=fill,lw=1.1)
            ax.axhline(DEGREES[angle],color=".65",ls=":",lw=.7)
            ax.set(xscale="log",ylim=(-5,190),xlabel="Cue duration (native units)",ylabel="Activity heading (°)")
            ax.set_title(NAMES[model]+f" · {DEGREES[angle]:g}° cue",loc="left",fontsize=9.5)
    fig.suptitle("Optional focused candidate · 90° and 180° response",fontsize=14,y=.995)
    fig.text(.05,.953,"This is a candidate view, not a main-panel selection. All four angles remain in the primary response atlases.",fontsize=8.5)
    handles=[Line2D([],[],color=c,label=l) for c,l in zip(SIZECOLORS,("Small: C20 / C₀1","Middle: C30 / C₀2","Large: C40 / C₀3"))]
    handles += [Line2D([],[],color=".3",ls=s,marker="o",mfc=m,label=l) for s,m,l in (("--","white","At withdrawal"),("-",".3","After withdrawal"))]
    fig.legend(handles=handles,ncol=3,frameon=False,loc="lower center",bbox_to_anchor=(.5,.015),fontsize=8)
    fig.tight_layout(rect=(0,.10,1,.925),h_pad=2,w_pad=2)
    footer(fig,complete,generated);return save(fig,"stationary_sizes_optional_90_180",out)


def evidence(frame,out,generated):
    lines=["# Three-size stationary-cue evidence", "",f"Generated {generated}. This is a review checkpoint of the complete requested grid; the author-facing main composition remains open. No native integration is performed by this renderer.","",
        "## Matched native examples across every cue angle", "",
        "The entries below are retained activity headings in degrees, after 13,000 cue-free native units. Chemical examples use T=100; Spring uses T=10. They illustrate the complete response curves rather than define an effective-duration threshold. Sizes are C20/30/40 for chemical models and resting C₀1/2/3 for Spring.","",
        "| Model | Cue angle | Native T | Small | Middle | Large |","|---|---:|---:|---:|---:|---:|"]
    for model in MODELS:
        T=10. if model=="spring" else 100.
        for a in ANGLES:
            g=frame[(frame.model==model)&(frame.angle_label==a)&(frame.cue_duration==T)].sort_values("size_parameter")
            if len(g)!=3 or not (g.status=="complete").all():continue
            vals=[f"{number(r,'angle_final_deg'):.2f}" if number(r,"R_final")>=.1 else "masked" for _,r in g.iterrows()]
            lines.append(f"| {NAMES[model]} | {DEGREES[a]:g}° | {T:g} | "+" | ".join(vals)+" |")
    lines += ["", "Across these fixed-physical-cue protocols, larger rings generally retain smaller non-antipodal turns for the same native cue duration. The response curves retain all four angles rather than selecting only antipodal reversals. Antipodal conditions can change between return, reversal and weak-resultant endpoints on the tested grid; connecting segments do not identify a refined threshold or a distinct binary mechanism.","",
        "At Goryachev 90°/T100, withdrawal angles are 59.41°, 22.48° and 9.39° for C20/30/40, while retained angles are 89.09°, −0.014° and −0.003°. The smaller ring continues toward the cue after withdrawal; the other two return near their original axes. The frozen original-region excess changes from about .331/.686/.858 at withdrawal to approximately 0/1/1 at the final endpoint. This is measured redistribution and rebound in a fixed basal material region, not proof of a uniquely tracked site being destroyed or recreated.","",
        "## Magnitude and qualification",""]
    completed=frame[frame.status=="complete"]
    nonanti=completed[completed.angle_label!="pi"]
    if len(completed)==240:
        residual=float(nonanti.final_primary_rotated_basal_relative_contrast_error.max())
        lines.append(f"All {len(nonanti)} non-antipodal final profiles have one resolved peak and match a rotated copy of their own basal active profile to at most {residual:.3g} of its contrast. Their final R differs from 1 by at most {float(abs(nonanti.R_final-1).max()):.3g}. Thus a small retained turn is not generally a weak-activity endpoint.")
        lines.append("\nThe two final weak-resultant cases are Holmes C30/180°/T300 (R≈.0996, independently failed full-field tail gate) and C40/180°/T1000 (R≈.0399). Their headings remain masked. The separate field-tail flags are 228 passes, two failures and ten unavailable qualifications from compact retained histories. The other failure is Holmes C40/180°/T300. No unavailable flag is promoted to a pass from its angle alone.")
    lines += ["", "## What the size comparison controls", "",
        "The requested cue-centre angles and native durations are shared within each model; peak and physical support retain their native definitions. Fixed physical support means angular cue width shrinks with circumference, while the physical cue offset at a fixed angle grows. Initial core width, overlap, basal amount and first-moment strength also vary across sizes. The linked context table retains those quantities, including physical M₀, |M₁| and q₁, so P normalization does not make the basals appear quantitatively identical. These experiments measure the joint size response of the declared protocol and do not isolate a purely geometric contribution.","",
        "All initial overlap cases remain. Physical total activity is distinguished from resultant magnitude; the fixed original-region minimum and rebound are sampled observables, not lineage classifiers. Deterministic histories do not acquire independent-seed uncertainty from repeated time samples. Native clocks/doses remain uncalibrated across models, and activity orientation is not measured cell movement.","",
        "- [Single all-model, all-angle overview](stationary_sizes_all_angles_overview.png)",
        "- [Basal and cue context](supporting/data/stationary_sizes_basal_cue_context.csv)",
        "- [Matched descriptive size contrasts](supporting/data/stationary_sizes_matched_contrasts.csv)",
        "- [Retained exact source table snapshot](supporting/data/source_summary_snapshot.csv)",
        "- [Coverage and generation-time validation](supporting/metadata/validation.json)", ""]
    (out/"EVIDENCE.md").write_text("\n".join(lines))


def write_docs(frame,lookup,expected,provenance,artifacts,out,generated,producer):
    complete=sum(done for _,done in lookup.values())
    rows=[]
    for key in sorted(expected):
        r,done=lookup.get(key,({},False))
        rows.append(dict(model=key[0],size_parameter=key[1],angle_label=key[2],cue_duration=key[3],analysis_available=done,condition_id=r.get("condition_id"),native_status=r.get("status","not_in_analysis_snapshot")))
    pd.DataFrame(rows).to_csv(out/"supporting/data/stationary_sizes_coverage.csv",index=False)
    status=dict(generated_at=generated,producer_state=producer,analyzed_complete=complete,expected_cases=240,all_four_angles_primary=True,
        native_integrations_performed_by_this_script=0,provenance=provenance,renderer=str(Path(__file__)),renderer_sha256=sha_bytes(Path(__file__).read_bytes()),
        artifacts=[dict(path=str(p),sha256=sha_bytes(p.read_bytes())) for p in artifacts],visual_inspection="pending")
    (out/"supporting/metadata/validation.json").write_text(json.dumps(status,indent=2))
    main=[dict(title="All four models and all four angles",path=str(out/"stationary_sizes_all_angles_overview.png"))]
    support=[dict(title=NAMES[m]+": enlarged withdrawal/retention atlas",path=str(out/("stationary_sizes_angles_"+m+".png"))) for m in MODELS]
    support.append(dict(title="Basal size and cue geometry context",path=str(out/"stationary_sizes_basal_cue_context.png")))
    support += [dict(title=NAMES[m]+": magnitude, amount and region",path=str(out/("stationary_sizes_pattern_amount_"+m+".png"))) for m in MODELS]
    legend=["Three-size stationary-cue response review. All four prescribed cue angles are retained for each native model and each requested size: chemical C20/30/40 and Spring resting C₀1/2/3. Colour denotes size, open/dashed measurements are at cue withdrawal and filled/solid measurements follow 13,000 cue-free native units. Heading is arg P for physical P=M₁(t)/M₁(0), using instantaneous C(t) in Spring. Every angle column shares a degree range; duration clocks remain model-specific. Lines connect tested conditions and do not establish refined duration thresholds.",
        "Weak-resultant outcomes R<.1 have no interpreted heading and are marked at the panel bottom (+ withdrawal, × final). A square outline independently flags failure of the full-field final-tail gate; complete tail qualification and unavailable flags remain linked. Local fields and total active amount remain necessary to distinguish cancellation, redistribution and retained orientation.",
        "The context table restores basal physical M₀, |M₁| and q₁ removed from view by initial-moment normalization, plus basal/core dimensions and actual cue footprint/integral/overlap. Fixed physical support changes angular exposure and overlap across sizes, so this protocol does not isolate geometry from all other basal/input changes. All overlap cases remain; no site-lineage or movement law is assigned."]
    (out/"supporting/metadata/CAPTION.md").write_text("# Three-size response review caption\n\n"+"\n\n".join(legend)+"\n")
    (out/"supporting/metadata/content_checkpoint.json").write_text(json.dumps(dict(title="Three-size stationary-cue response review",role="review checkpoint; no final main-layout selection",assets=main,supplementary=support,legend=legend,generated_at=generated,analyzed_complete=complete,expected_cases=240),indent=2))
    links="\n".join(f"- [{NAMES[m]}: all four angles](stationary_sizes_angles_{m}.png)" for m in MODELS)
    readme=f"""# Three-size stationary-cue response review

Generated {generated}. The authoritative analysis snapshot contains {complete}/240 analyzed histories. Producer state: {producer}. {"The full requested grid is available." if complete==240 else "Missing conditions remain visibly pending and are never replaced by zeros."} This renderer performs no native integration and does not select a final main layout.

{links}
- [Single all-model, all-angle overview](stationary_sizes_all_angles_overview.png)
- [Basal and cue geometry context](stationary_sizes_basal_cue_context.png)
- [Physical basal magnitudes, cue footprint and overlap table](supporting/data/stationary_sizes_basal_cue_context.csv)
- [All240 expected-condition coverage](supporting/data/stationary_sizes_coverage.csv)
- [Matched larger-minus-smaller descriptive comparisons](supporting/data/stationary_sizes_matched_contrasts.csv)
- [Bounded quantitative evidence](EVIDENCE.md)

The primary review overview retains four model rows and all four cue-angle columns, with size in colour and withdrawal/retention in line and marker style. Every angle column has a shared degree range; every model keeps its native duration axis. Bottom +/× markers indicate weak-resultant withdrawal/final samples without assigning their displayed vertical position an angle. An enlarged per-model withdrawal/retention atlas is supporting evidence. Chemical sizes are C20/30/40; Spring uses resting C₀1/2/3 with its actual physical circumference. Durations in the enlarged cell atlases are tested categories; those cell spacings do not imply equal time. Actual measured heading is displayed only for R≥.1; weak vectors remain valid measurements in the linked magnitude tables. A black border in an enlarged final atlas flags failure of the separate final full-field tail gate. Unavailable tail qualification is retained in the authoritative table rather than invented from an angle endpoint.

P=M₁(t)/M₁(0) uses the initial complex physical activity moment, with M₁=C(t)mean[u exp(iθ)]. This rotates the initial axis and scales its initial magnitude once; it does not divide by current total activity. The context table retains initial physical M₀, |M₁| and q₁=|M₁|/M₀ so normalization does not hide basal magnitude differences. It also retains original-core physical width, actual basal circumference, native cue peak/integral, initial angular width and support/core overlap. Fixed physical support changes angular exposure as size changes; these tests do not isolate geometry from all changes in basal organization or cue overlap.

Separate pattern/amount sheets retain final R and physical M₀/M₀(0), the cue-on sampled minimum in the fixed original basal material region and rebound from the full-history sampled minimum. A regional minimum is not a continuous-time extremum, and depletion/rebound need not mean loss/rebirth of an independently tracked site. All overlap cases remain. These deterministic histories have no invented replicate uncertainty. Native clocks and doses are not calibrated across models, and activity orientation is not measured movement.

The optional90°/180° candidate, when requested, is a secondary view rather than a settled main-panel selection. All four angles remain primary review material. Native generation, exact reuses, physical-moment/core definitions, source hashes and censoring-aware crossings are owned by the authoritative campaign producer.
"""
    (out/"README.md").write_text(readme)
    document=f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Three-size stationary-cue review</title><style>body{{font:16px/1.5 system-ui;max-width:1100px;margin:2em auto;padding:0 1em}}img{{max-width:100%}}a{{color:#176181}}</style><h1>Three-size stationary-cue response review</h1><p>Input snapshot modified {html.escape(provenance["source_snapshot_modified_at"])} · rendered {html.escape(generated)} · analyzed {complete}/240 · producer {html.escape(producer)}. Opening or manually refreshing loads this snapshot; there is no timed browser refresh.</p>'
    for a in main:document+=f'<h2>{html.escape(a["title"])}</h2><img src="{Path(a["path"]).name}" alt="{html.escape(a["title"])}">'
    document+='<h2>Supporting views</h2><ul>'
    for a in support:document+=f'<li><a href="{Path(a["path"]).name}">{html.escape(a["title"])}</a></li>'
    document+='</ul>'
    document+='<p><a href="README.md">Definitions and limits</a> · <a href="supporting/data/stationary_sizes_basal_cue_context.csv">Basal/cue context</a> · <a href="supporting/data/stationary_sizes_coverage.csv">Complete expected coverage</a></p></html>'
    (out/"index.html").write_text(document)
    evidence(frame,out,generated)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path,default=TOPIC/"supporting/data/stationary_sizes_summary.csv")
    parser.add_argument("--manifest",type=Path)
    parser.add_argument("--out",type=Path,default=TOPIC/"response")
    parser.add_argument("--producer-state",choices=("preparing","running","complete","unknown"),default="unknown")
    parser.add_argument("--include-candidate",action="store_true")
    args=parser.parse_args()
    for folder in (args.out,args.out/"supporting/data",args.out/"supporting/metadata",args.out/"supporting/vectors"):folder.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    frame,lookup,expected,provenance=load(args.summary.resolve(),args.manifest.resolve() if args.manifest else None)
    snapshot=args.out/"supporting/data/source_summary_snapshot.csv";snapshot.write_bytes(frame.attrs.pop("source_snapshot_bytes"));provenance["retained_source_snapshot_path"]=str(snapshot)
    complete=sum(done for _,done in lookup.values());generated=datetime.now(timezone.utc).isoformat()
    context=context_table(frame,args.out)
    paired_size_contrasts(lookup,args.out)
    artifacts=all_angles_overview(lookup,args.out,complete,generated)
    artifacts+=response_atlas(lookup,args.out,complete,generated)
    artifacts+=plot_context(context,args.out,complete,generated)
    artifacts+=local_and_amount(lookup,args.out,complete,generated)
    if args.include_candidate:artifacts+=candidate(lookup,args.out,complete,generated)
    write_docs(frame,lookup,expected,provenance,artifacts,args.out,generated,args.producer_state)
    print(json.dumps(dict(output=str(args.out),analyzed=complete,expected=240,artifacts=len(artifacts))))


if __name__=="__main__":main()
