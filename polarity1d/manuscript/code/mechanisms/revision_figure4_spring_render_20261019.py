"""Two readable Figure 4 sheets with Spring; render-only scientific inputs."""
from pathlib import Path
from datetime import datetime,timezone
import csv
import json
import os
import math
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from revision_figure4_20261019 import ROOT,OUT,OLD,DATA,META,VECTORS,read_rows,select,sha

RUN=ROOT/"polarity1d/manuscript/data/runs/revision_spring_noise_main_20261019"
MODELS=("wave_pinning","holmes_model3","spring_ring")
NAMES={"wave_pinning":"Wave-pinning","holmes_model3":"Holmes 3","spring_ring":"Spring ring"}
COLORS=("#397bab","#40916c","#cf7930")


def values(rows,key):return np.asarray([float(r[key]) for r in rows])
def save(fig,stem):
    paths=[]
    for suffix in ("png","pdf","svg"):
        p=(OUT if suffix=="png" else VECTORS)/(stem+"."+suffix)
        fig.savefig(p,dpi=300,bbox_inches="tight");paths.append(p)
    plt.close(fig);return paths


def inputs():
    audit=json.loads((META/"spring_noise_analysis_audit.json").read_text());extra=json.loads((META/"spring_noise_supplement_audit.json").read_text())
    assert audit["cases"]==81 and audit["noisy_cases"]==72 and extra["cases"]==81
    assert extra["maximum_material_centroid_error"]<1e-10
    direction={}
    chemical=read_rows(OLD/"supporting/data/vision_noise_direction_correlations.csv",lambda r:r["condition"]=="noise" and r["window"]=="full" and float(r["gate"])==.1 and float(r["lag"])%40==0)
    spring=read_rows(DATA/"spring_noise_direction_correlations.csv",lambda r:r["condition"]=="noise" and r["window"]=="full" and float(r["gate"])==.1 and float(r["lag"])%4==0)
    endpoints={};observation={};representative={}
    ep=read_rows(DATA/"endpoint_statistics.csv",lambda r:r["condition"]=="noise" and float(r["fraction"])==.01 and r["gate"]=="strict_resolved_site" and float(r["lag"])==8000)
    sp=read_rows(DATA/"spring_noise_site_statistics.csv",lambda r:r["condition"]=="noise" and float(r["fraction"])==.01 and r["gate"]=="strict" and float(r["lag"])==800)
    ob=read_rows(DATA/"observation_statistics.csv",lambda r:r["condition"]=="noise" and float(r["fraction"])==.01)
    so=read_rows(DATA/"spring_noise_observation_statistics.csv",lambda r:r["condition"]=="noise" and float(r["fraction"])==.01)
    reps=read_rows(OLD/"supporting/data/vision_noise_representative_histories.csv",lambda r:r["role"]=="within_condition_descriptor_medoid" and float(r["fraction"])==.01)
    sr=read_rows(DATA/"spring_noise_representative_histories.csv",lambda r:float(r["fraction"])==.01)
    for model,C in (("wave_pinning",20),("holmes_model3",60),("spring_ring",2)):
        if model=="spring_ring":
            direction[model]=spring;endpoints[model]=sp;observation[model]=so
            representative[model],=select(sr,resting_circumference=C)
        else:
            direction[model]=select(chemical,model=model);endpoints[model]=select(ep,model=model);observation[model]=select(ob,model=model)
            representative[model],=select(reps,model=model,circumference=C)
    return direction,endpoints,observation,representative,audit,extra


def history(ax,model,row,label,common_vmax=None):
    with np.load(row["profile_archive"]) as z:
        t,u=z["time"],z["primary"]
        q=z["material_positions"] if model=="spring_ring" else z["x"]/float(row["circumference"])
    p=(q+.5)%1-.5;order=np.argsort(p)
    mesh=ax.pcolormesh(p[order],t/1000,u[:,order],vmin=0,vmax=common_vmax or float(u.max()),cmap="viridis",shading="nearest",rasterized=True)
    ax.set_ylim(t[-1]/1000,0);ax.set_xlim(-.5,.5)
    C=float(row["resting_circumference"] if model=="spring_ring" else row["circumference"])
    ax.set_title(f"{label}  {NAMES[model]}, {'C₀' if model=='spring_ring' else 'C'}={C:g}",loc="left",fontweight="bold",fontsize=10)
    ax.set_xlabel("Material boundary position / cycle")
    ax.set_ylabel("Time (10³ native units)")
    cb=plt.colorbar(mesh,ax=ax,fraction=.043,pad=.02);cb.set_label("Native active field",fontsize=8);cb.ax.tick_params(labelsize=7)
    ax.text(.02,.04,"seed "+row["seed"],transform=ax.transAxes,color="white",fontsize=7,bbox=dict(facecolor="black",edgecolor="none",alpha=.35,pad=2))


def sizes(model):return (1.,2.,3.) if model=="spring_ring" else (20.,40.,60.)
def size_key(model):return "resting_circumference" if model=="spring_ring" else "circumference"
def xlab(model):return "Resting circumference C₀ (Spring units)" if model=="spring_ring" else "Circumference C (native units)"


def main_sheets(direction,endpoints,observation,reps):
    artifacts=[]
    direction_lower=min(float(r["ci_lower"]) for model in MODELS for r in direction[model] if np.isfinite(float(r["ci_lower"])))
    fig,axes=plt.subplots(3,3,figsize=(12.3,11.1),gridspec_kw=dict(height_ratios=(1.35,1,1.05)))
    for j,model in enumerate(MODELS):
        history(axes[0,j],model,reps[model],"A"+str(j+1))
        ax=axes[1,j]
        for C,color in zip(sizes(model),COLORS):
            for f,style in ((.001,"--"),(.01,"-")):
                r=select(direction[model],fraction=f,**{size_key(model):C});x=values(r,"lag")
                divisor=1. if model=="spring_ring" else 1000.
                ax.plot(x/divisor,values(r,"direction_correlation"),color=color,ls=style)
                if f==.01:ax.fill_between(x/divisor,values(r,"ci_lower"),values(r,"ci_upper"),color=color,alpha=.13,lw=0)
        ax.set_ylim(min(-.08,direction_lower-.03),1.04);ax.axhline(0,color=".8",lw=.5)
        ax.set_xlabel("Lag (Spring native time)" if model=="spring_ring" else "Lag (10³ native time)")
        ax.set_ylabel("Activity-direction cosine")
        ax.set_title(f"B{j+1}  Directional memory",loc="left",fontweight="bold",fontsize=10)
        ax.set_xlim(0,800 if model=="spring_ring" else 8)
        ax=axes[2,j]
        metric="basal_arc_squared" if model=="spring_ring" else "displacement_squared"
        r=select(endpoints[model],metric=metric)
        unavailable=[]
        for row,color in zip(r,COLORS):
            C=float(row[size_key(model)]);y=np.sqrt(float(row["value"]));lo=np.sqrt(float(row["ci_lower"]));hi=np.sqrt(float(row["ci_upper"]))
            if not np.isfinite(y) or int(row["independent_seeds_with_pairs"])<6 or float(row["eligible_pair_fraction"])<.01:
                unavailable.append(f"{C:g}");continue
            ax.plot([C,C],[lo,hi],color=color,lw=1.2);ax.scatter(C,y,color=color,s=27)
            ax.annotate(f"{float(row['eligible_pair_fraction']):.0%}",(C,y),xytext=(4,6),textcoords="offset points",fontsize=7)
        ax.set_ylim(bottom=0)
        if model!="spring_ring":ax.set_ylim(0,7.8)
        else:
            basals=json.loads((RUN/"manifest.json").read_text())["basals"]
            ax.plot([b["resting_circumference"] for b in basals],[b["basal_circumference"]/np.sqrt(12) for b in basals],color=".5",ls=":",lw=1,label="Uniform-angle endpoints")
            ax.legend(frameon=False,fontsize=7.1,loc="upper left")
        ax.set_xticks(sizes(model));ax.set_xlabel(xlab(model));ax.set_ylabel("Folded basal-arc endpoint RMS\n(Spring native length)" if model=="spring_ring" else "Site RMS displacement\n(native length units)")
        ax.set_title(f"C{j+1}  {'Reference-arc shift' if model=='spring_ring' else 'Physical site displacement'}",loc="left",fontweight="bold",fontsize=9.5)
        if unavailable:ax.text(.03,.05,"Sparse endpoints: "+", ".join(unavailable),transform=ax.transAxes,fontsize=7)
    fig.suptitle("Figure 4 · Directional memory needs a physical site readout",fontsize=15,y=.995)
    fig.text(.5,.956,"Memory: both declared strengths; histories/site metrics: f=.01 · 12 independent seeds · native units differ across models",ha="center",fontsize=8.7)
    fig.legend(handles=[Line2D([],[],color=color,label=lab) for color,lab in zip(COLORS,("Small: C20 / Spring C₀1","Middle: C40 / Spring C₀2","Large: C60 / Spring C₀3"))]+
               [Line2D([],[],color=".3",ls=s,label=l) for s,l in (("--","f=.001"),("-","f=.01 + 95% seed interval"))],loc="lower center",bbox_to_anchor=(.5,-.005),ncol=3,frameon=False,fontsize=8.3)
    fig.tight_layout(rect=(0,.068,1,.92),h_pad=3.0,w_pad=2.0)
    artifacts+=save(fig,"figure4_with_spring_dynamics")
    fig,axes=plt.subplots(2,3,figsize=(12.2,7.3))
    largest_obs=max(float(r["ci_upper"]) for model in MODELS for r in observation[model] if r["metric"] in ("multiregion_fraction","weak_R_fraction"))
    for j,model in enumerate(MODELS):
        ax=axes[0,j]
        geometry_data=read_rows(DATA/"spring_noise_site_statistics.csv",lambda r:r["condition"]=="noise" and float(r["fraction"])==.01 and r["gate"]=="strict" and float(r["lag"])==4) if model=="spring_ring" else endpoints[model]
        for metric,style,label in (("native_cosine","-","Native ring"),("reference_cosine","--","Common reference ring")):
            r=[r for r in select(geometry_data,metric=metric) if int(r["independent_seeds_with_pairs"])>=6 and float(r["eligible_pair_fraction"])>=.01]
            x=values(r,size_key(model));y=values(r,"value")
            ax.plot(x,y,style,color="#6552a1",marker="o",ms=3.5,label=label)
            ax.fill_between(x,values(r,"ci_lower"),values(r,"ci_upper"),color="#6552a1",alpha=.09,lw=0)
        lower=min((float(r["ci_lower"]) for r in geometry_data if r["metric"] in ("native_cosine","reference_cosine") and np.isfinite(float(r["ci_lower"]))),default=0.)
        ax.set_xticks(sizes(model));ax.set_ylim(lower-.015 if model=="spring_ring" else min(0.,lower-.03),1.008 if model=="spring_ring" else 1.04);ax.set_xlabel(xlab(model));ax.set_ylabel("Eligible site-direction cosine")
        ax.set_title(f"D{j+1}  {NAMES[model]} · lag {4 if model=='spring_ring' else 8000}",loc="left",fontweight="bold",fontsize=9.2)
        ax.legend(frameon=False,fontsize=7.5,loc="lower right")
        ax=axes[1,j]
        for metric,color,style,label in (("multiregion_fraction","#327d9a","-","Multiple half-contrast regions"),("weak_R_fraction","#b24973","--","Weak resultant R < .1")):
            r=select(observation[model],metric=metric);x=values(r,size_key(model));y=values(r,"value")
            ax.plot(x,y,style,color=color,marker="o",ms=3.5,label=label)
            ax.fill_between(x,values(r,"ci_lower"),values(r,"ci_upper"),color=color,alpha=.12,lw=0)
        ax.set_xticks(sizes(model));ax.set_ylim(-.008,max(.1,1.08*largest_obs));ax.set_xlabel(xlab(model));ax.set_ylabel("Fraction of observed sample times")
        ax.set_title(f"E{j+1}  Local pattern and magnitude",loc="left",fontweight="bold",fontsize=9.5)
        ax.legend(frameon=False,fontsize=7.1,loc="upper left")
    fig.suptitle("Figure 4 continued · Geometry and passive local-pattern signatures",fontsize=14,y=.995)
    fig.text(.5,.94,"Same endpoints in each native/reference pair · chemical C*=60; Spring C*=3.04517 · Spring uses early lag4 before endpoint mixing",ha="center",fontsize=8.5)
    fig.tight_layout(rect=(0,.015,1,.91),h_pad=3,w_pad=2)
    artifacts+=save(fig,"figure4_with_spring_geometry_patterns")
    return artifacts


def spring_support():
    artifacts=[]
    reps=read_rows(DATA/"spring_noise_representative_histories.csv")
    vmax=max(float(np.load(r["profile_archive"])["primary"].max()) for r in reps)
    fig,axes=plt.subplots(2,3,figsize=(11.5,7.4),sharey=True)
    for i,f in enumerate((.001,.01)):
        for j,C in enumerate((1.,2.,3.)):
            row,=select(reps,fraction=f,resting_circumference=C);history(axes[i,j],"spring_ring",row,f"f={f:g}",common_vmax=vmax)
    fig.suptitle("Spring native histories: every prespecified condition medoid",fontsize=13)
    fig.tight_layout(rect=(0,0,1,.95));artifacts+=save(fig,"figure4_spring_all_histories")
    traces=read_rows(DATA/"spring_noise_ensemble_traces.csv")
    corr=read_rows(DATA/"spring_noise_correlations.csv",lambda r:float(r["lag"])%4==0)
    direction=read_rows(DATA/"spring_noise_direction_correlations.csv",lambda r:float(r["gate"])==.1 and float(r["lag"])%4==0)
    fig,axes=plt.subplots(3,2,figsize=(11,9.0),sharex=True)
    for i,(metric,label) in enumerate((("R","Relative first-moment magnitude R"),("initial_memory","Initial-axis memory Re P"),("circumference","Physical circumference C(t)"))):
        for j,f in enumerate((.001,.01)):
            ax=axes[i,j]
            for C,color in zip((1.,2.,3.),COLORS):
                for kind,style in (("noise","-"),("mean","--"),("basal",":")):
                    r=select(traces,resting_circumference=C,fraction=0. if kind=="basal" else f,condition=kind,metric=metric)
                    ax.plot(values(r,"time"),values(r,"value"),ls=style,color=color,lw=1.15)
                    if kind=="noise":ax.fill_between(values(r,"time"),values(r,"ci_lower"),values(r,"ci_upper"),color=color,alpha=.10,lw=0)
            ax.set_title(f"f={f:g}",loc="left");ax.set_ylabel(label)
            if i==2:ax.set_xlabel("Spring native time")
    fig.suptitle("Spring magnitude, initial memory and dynamic geometry with all controls",fontsize=13)
    fig.legend(handles=[Line2D([],[],color=c,label=f"C₀={C}") for C,c in zip((1,2,3),COLORS)]+[Line2D([],[],color=".3",ls=s,label=l) for s,l in (("-","Noise"),("--","Mean input"),(":","No added input"))],loc="lower center",bbox_to_anchor=(.5,-.02),ncol=6,frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,.03,1,.96));artifacts+=save(fig,"figure4_spring_magnitude_geometry_controls")
    fig,axes=plt.subplots(2,3,figsize=(12,7.4),sharex=True,sharey=True)
    lower=min(float(r["ci_lower"]) for r in corr if r["condition"]=="noise" and np.isfinite(float(r["ci_lower"])))
    for i,f in enumerate((.001,.01)):
        for j,C in enumerate((1.,2.,3.)):
            ax=axes[i,j]
            for window,style in (("full","-"),("first_half","--"),("last_half",":")):
                r=select(corr,resting_circumference=C,fraction=f,condition="noise",window=window)
                ax.plot(values(r,"lag"),values(r,"correlation"),color=COLORS[j],ls=style,lw=1.2)
                if window=="full":ax.fill_between(values(r,"lag"),values(r,"ci_lower"),values(r,"ci_upper"),color=COLORS[j],alpha=.13,lw=0)
            d=select(direction,resting_circumference=C,fraction=f,condition="noise",window="full")
            ax.plot(values(d,"lag"),values(d,"direction_correlation"),color=".3",lw=.8,ls="-.")
            ax.set_title(f"C₀={C:g}, f={f:g}",loc="left");ax.set_ylim(min(-.1,lower-.03),1.05);ax.set_ylabel("Lag correlation")
            if i==1:ax.set_xlabel("Lag (Spring native time)")
    fig.suptitle("Spring finite-window memory; no stationarity assumption or exponential fit",fontsize=13)
    fig.legend(handles=[Line2D([],[],color=".3",ls=s,label=l) for s,l in (("-","Full vector: full window"),("--","First half"),(":","Last half"),("-.","R≥.1 direction: full window"))],loc="lower center",bbox_to_anchor=(.5,-.015),ncol=4,frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,.045,1,.96));artifacts+=save(fig,"figure4_spring_memory_windows")
    cdf=read_rows(DATA/"spring_noise_increment_cdf.csv",lambda r:float(r["lag"])==.5 and float(r["gate"])==.1)
    fig,axes=plt.subplots(1,3,figsize=(11.5,3.8),sharex=True,sharey=True)
    for j,C in enumerate((1.,2.,3.)):
        ax=axes[j]
        for f,style,color in ((.001,"--","#337fab"),(.01,"-","#bd6036")):
            r=select(cdf,resting_circumference=C,fraction=f)
            ax.plot(values(r,"threshold_radians"),1-values(r,"cdf"),ls=style,color=color,label=f"f={f:g}")
            ax.fill_between(values(r,"threshold_radians"),1-values(r,"ci_upper"),1-values(r,"ci_lower"),color=color,alpha=.12,lw=0)
        ax.set_xscale("log");ax.set_yscale("log");ax.set_xlim(1e-6,np.pi);ax.set_ylim(1e-5,1.05)
        ax.set_title(f"Spring C₀={C:g}",loc="left");ax.set_xlabel("Absolute .5-unit activity-angle increment (rad)")
    axes[0].set_ylabel("Fraction above increment");axes[0].legend(frameon=False,fontsize=8)
    fig.suptitle("Spring descriptive increment tails, both endpoint R≥.1; no event classifier",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.94));artifacts+=save(fig,"figure4_spring_increment_tails")
    site=read_rows(DATA/"spring_noise_site_statistics.csv",lambda r:r["condition"]=="noise" and r["gate"]=="strict" and float(r["lag"])%4==0)
    fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True)
    for j,C in enumerate((1.,2.,3.)):
        for f,color in ((.001,"#397bab"),(.01,"#c16b31")):
            r=select(site,resting_circumference=C,fraction=f,metric="basal_arc_squared")
            axes[0,j].plot(values(r,"lag"),np.sqrt(values(r,"value")),color=color)
            axes[0,j].fill_between(values(r,"lag"),np.sqrt(values(r,"ci_lower")),np.sqrt(values(r,"ci_upper")),color=color,alpha=.12,lw=0)
            basal=float(r[0]["basal_circumference"])
            axes[0,j].axhline(basal/np.sqrt(12),color=".55",ls=":",lw=.8)
            for metric,style in (("native_cosine","-"),("reference_cosine","--")):
                rr=select(site,resting_circumference=C,fraction=f,metric=metric)
                axes[1,j].plot(values(rr,"lag"),values(rr,"value"),color=color,ls=style,lw=1.2)
            axes[2,j].plot(values(r,"lag"),values(r,"eligible_pair_fraction"),color=color)
        for i in range(3):
            axes[i,j].set_title(f"Spring C₀={C:g}",loc="left",fontsize=9)
            axes[i,j].axvline(4,color=".75",lw=.6);axes[i,j].axvline(50,color=".75",lw=.6)
            axes[i,j].set_xlim(0,800)
        axes[0,j].set_ylabel("Folded basal-arc endpoint RMS")
        axes[1,j].set_ylabel("Site cosine: native / remapped")
        axes[2,j].set_ylabel("Eligible endpoint fraction");axes[2,j].set_ylim(0,1.04);axes[2,j].set_xlabel("Lag (Spring native time)")
    fig.suptitle("Spring endpoint geometry across lag: short shifts and late circular mixing",fontsize=13)
    fig.legend(handles=[Line2D([],[],color=c,label=l) for c,l in (("#397bab","f=.001"),("#c16b31","f=.01"))]+
               [Line2D([],[],color=".3",ls=s,label=l) for s,l in (("-","Native site cosine"),("--","Common-C* site cosine"),(":","Uniform-angle RMS reference"))],
               loc="lower center",bbox_to_anchor=(.5,-.01),ncol=3,frameon=False,fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.96));artifacts+=save(fig,"figure4_spring_endpoint_lag_sensitivity")
    return artifacts


def endpoint_window_summary():
    site=read_rows(DATA/"spring_noise_site_statistics.csv",lambda r:r["condition"]=="noise" and r["gate"]=="strict" and float(r["lag"]) in (4.,50.,800.))
    manifest=json.loads((RUN/"manifest.json").read_text());summary=[]
    for C in (1.,2.,3.):
        for f in (.001,.01):
            jobs=[j for j in manifest["jobs"] if j["kind"]=="noise" and j["resting_circumference"]==C and j["fraction"]==f]
            traces=[]
            for j in jobs:
                folder=RUN/"workers"/f"{j['slot']:04d}";status=json.loads((folder/"status.json").read_text())
                assert sha(folder/"trace.npz")==status["output_sha256"]["trace.npz"]
                with np.load(folder/"trace.npz") as z:
                    P=z["P_real"]+1j*z["P_imag"]
                    strict=(abs(P)>=.1)&(z["current_half_contrast_regions"]==1)&(z["peaks_prominence_ge_005"]<=1)&(z["contrast"]>=.05)
                    traces.append((z["dominant_region_fraction"].copy(),strict))
            for lag in (4.,50.,800.):
                rr=select(site,resting_circumference=C,fraction=f,lag=lag);by={r["metric"]:r for r in rr};base=by["basal_arc_squared"]
                r0=float(base["basal_circumference"]);rms=np.sqrt(float(base["value"]));uniform=r0/np.sqrt(12)
                k=int(lag/.5);samples=[]
                for q,v in traces:
                    mask=v[k:]&v[:-k];samples.append(abs((q[k:]-q[:-k]+.5)%1-.5)[mask])
                samples=np.sort(np.concatenate(samples));n=len(samples)
                distance=float(max(np.max(np.arange(1,n+1)/n-2*samples),np.max(2*samples-np.arange(n)/n))) if n else np.nan
                summary.append(dict(resting_C=C,fraction=f,lag=lag,basal_C=r0,basal_arc_RMS=rms,
                    RMS_ci_lower=np.sqrt(float(base["ci_lower"])),RMS_ci_upper=np.sqrt(float(base["ci_upper"])),uniform_angle_RMS=uniform,RMS_over_uniform_reference=rms/uniform,
                    midpoint_arc_RMS=np.sqrt(float(by["midpoint_arc_squared"]["value"])),
                    native_site_cosine=float(by["native_cosine"]["value"]),reference_C_site_cosine=float(by["reference_cosine"]["value"]),
                    eligible_pair_fraction=float(base["eligible_pair_fraction"]),independent_seeds_with_pairs=int(base["independent_seeds_with_pairs"]),
                    descriptive_CDF_distance_from_uniform=distance,uniform_CDF_comparison="descriptive only; temporal origins are correlated; no p-value or uniformity claim"))
    p=DATA/"spring_noise_short_long_endpoint_summary.csv"
    with p.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    return summary


def increment_native_records():
    """Expose actual native context for sampled large global-angle changes."""
    manifest=json.loads((RUN/"manifest.json").read_text());records=[];examples={}
    for job in (j for j in manifest["jobs"] if j["kind"]=="noise"):
        folder=RUN/"workers"/f"{job['slot']:04d}";status=json.loads((folder/"status.json").read_text())
        assert sha(folder/"trace.npz")==status["output_sha256"]["trace.npz"]
        with np.load(folder/"trace.npz") as z:
            p=z["P_real"]+1j*z["P_imag"];R=abs(p)
            a=abs(np.angle(p[1:]*np.conj(p[:-1])));valid=(R[1:]>=.1)&(R[:-1]>=.1)
            indices=np.where(valid)[0];indices=indices[np.argsort(a[indices])[-5:][::-1]]
            for k in indices:
                row=dict(resting_C=job["resting_circumference"],fraction=job["fraction"],seed=job["seed"],slot=job["slot"],
                    t0=float(z["time"][k]),t1=float(z["time"][k+1]),absolute_angle_radians=float(a[k]),
                    R0=float(R[k]),R1=float(R[k+1]),regions0=int(z["current_half_contrast_regions"][k]),regions1=int(z["current_half_contrast_regions"][k+1]),
                    nearest_saved_profile_before=4*math.floor(float(z["time"][k])/4),nearest_saved_profile_after=4*math.ceil(float(z["time"][k+1])/4),
                    trace_archive=str(folder/"trace.npz"),profile_archive=str(folder/"primary.npz"),
                    selection="five largest .5-unit increments per seed conditional on endpoint R>=.1; no event classifier")
                records.append(row)
                C=job["resting_circumference"]
                if job["fraction"]==.01 and (C not in examples or row["absolute_angle_radians"]>examples[C]["absolute_angle_radians"]):examples[C]=row
    path=DATA/"spring_noise_large_increment_native_context.csv"
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    fig,axes=plt.subplots(2,3,figsize=(12,7.3))
    for j,C in enumerate((1.,2.,3.)):
        row=examples[C];middle=(row["t0"]+row["t1"])/2;folder=Path(row["trace_archive"]).parent
        with np.load(folder/"primary.npz") as z:
            t=z["time"].copy();mask=abs(t-middle)<=12;u=z["primary"][mask];q=z["material_positions"].copy()
        x=(q+.5)%1-.5;order=np.argsort(x);ax=axes[0,j]
        mesh=ax.pcolormesh(x[order],t[mask]-middle,u[:,order],vmin=0,vmax=float(u.max()),cmap="viridis",shading="nearest",rasterized=True)
        ax.axhline(0,color="white",ls=":",lw=1);ax.set_ylim(12,-12);ax.set_xlim(-.5,.5)
        ax.set_title(f"C₀={C:g}, seed {row['seed']}",loc="left",fontsize=9)
        ax.set_xlabel("Material boundary position / cycle");ax.set_ylabel("Time from selected interval (native)")
        cb=fig.colorbar(mesh,ax=ax,fraction=.04,pad=.02);cb.set_label("Native active field",fontsize=7)
        with np.load(folder/"trace.npz") as z:
            t=z["time"];mask=abs(t-middle)<=12;p=z["P_real"][mask]+1j*z["P_imag"][mask];t=t[mask]-middle
        ax=axes[1,j]
        for y,color,style,label in ((p.real,".25","-","Re P"),(p.imag,".6","--","Im P"),(abs(p),"#b24972","-","R=|P|")):
            ax.plot(t,y,color=color,ls=style,lw=1.1,label=label)
        ax.axhline(.1,color="#b24972",ls=":",lw=.7);ax.axvspan(row["t0"]-middle,row["t1"]-middle,color="#cf7930",alpha=.2)
        ax.set_xlim(-12,12);ax.set_xlabel("Time from selected interval (native)");ax.set_ylabel("Physical activity vector / initial M₁")
        ax.set_title(f"Selected |Δangle|={row['absolute_angle_radians']:.3f} rad; t={row['t0']:g}–{row['t1']:g}",loc="left",fontsize=8.3)
        if j==0:ax.legend(frameon=False,fontsize=7,ncol=3,loc="lower left")
    fig.suptitle("Spring largest sampled angle increments: native context, not an event classification",fontsize=12.3)
    fig.text(.5,.932,"Selected maxima at f=.01 with both endpoint R≥.1 · scalar cadence .5, saved fields every4 · narrow interval shaded below",ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.01,1,.90),h_pad=2.5,w_pad=2)
    return save(fig,"figure4_spring_large_increment_context")


def evidence_report():
    direction=read_rows(DATA/"spring_noise_direction_correlations.csv",lambda r:r["condition"]=="noise" and r["window"]=="full" and float(r["gate"])==.1 and float(r["lag"])==800)
    point=read_rows(DATA/"spring_noise_short_long_endpoint_summary.csv",lambda r:float(r["lag"])==800)
    obs=read_rows(DATA/"spring_noise_observation_statistics.csv",lambda r:r["condition"]=="noise")
    seeds=read_rows(DATA/"spring_noise_seed_summary.csv",lambda r:r["condition"]=="noise")
    increments=read_rows(DATA/"spring_noise_seed_increments.csv",lambda r:r["condition"]=="noise" and float(r["lag"])==.5 and float(r["gate"])==.1)
    text=["# Figure 4: final four-model evidence and bounded claims", "",
        "The two main sheets compare Wave-pinning, Holmes 3 and Spring in their separate native units. Complete Goryachev results and its forced-state controls remain linked supporting evidence. All 243 original chemical cases and all 81 new Spring cases are retained, including both declared strengths, every original replicate and matched mean/no-input controls.","",
        "## Two claims", "",
        "1. Angular retention and physical endpoint displacement can respond differently to size. Wave-pinning's stronger-input angular retention rises while its site RMS decreases; Holmes' physical endpoint RMS is nonmonotone despite rising angular retention. Spring's weaker-input angular retention rises with size even though basal-reference arc RMS is nonmonotone. Stronger Spring input removes long-lag direction retention at all three sizes, so these results do not support a universal enlargement advantage.",
        "2. Pairwise remapping and passive local-pattern/magnitude measurements add distinct information. Remapping identical eligible site endpoints changes their angular interpretation without changing their physical shifts. Multiple local regions, weak global resultants and sampled increment tails describe observed organization; they do not establish functional competition, an abrupt-event mechanism, site ancestry or cell movement.","",
        "## Spring quantitative comparison", "",
        "All entries below use 12 independent seeds and lag 800 native units. Parentheses are pointwise 95% whole-seed bootstrap intervals. Site coverage is the fraction of candidate endpoint pairs passing the strict gate, not the fraction of time with proven site ancestry.","",
        "| f | Resting C₀ | Direction cosine | Basal-reference endpoint RMS | Eligible pairs | Mean R | Multiple-region time | Weak-R time |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for f in (.001,.01):
        for C in (1.,2.,3.):
            d,=select(direction,resting_circumference=C,fraction=f);s,=select(point,resting_C=C,fraction=f)
            o={r["metric"]:r for r in select(obs,resting_circumference=C,fraction=f)}
            text.append(f"| {f:g} | {C:g} | {float(d['direction_correlation']):.3f} ({float(d['ci_lower']):.3f}, {float(d['ci_upper']):.3f}) | {float(s['basal_arc_RMS']):.3f} ({float(s['RMS_ci_lower']):.3f}, {float(s['RMS_ci_upper']):.3f}) | {float(s['eligible_pair_fraction']):.1%} | {float(o['mean_R']['value']):.3f} | {float(o['multiregion_fraction']['value']):.3%} | {float(o['weak_R_fraction']['value']):.3%} |")
    text += ["", "At stronger input, lag-800 RMS values 0.298/0.592/0.863 approach the uniform wrapped-angle references 0.301/0.590/0.879. These folded endpoint measurements discard winding and cannot measure travel. In the main remapping panel Spring instead uses early lag 4: the small and middle rings have similar basal-arc RMS (0.0607/0.0657), native site cosines 0.940/0.980 and common-reference site cosines 0.9923/0.9910. The full curves and lag 4/50/800 summaries remain accessible. The CDF distance from uniform is descriptive only, with correlated temporal origins and no p-value.", "",
        "The largest observed relative circumference excursion across all noisy trajectories is "+f"{max(float(r['maximum_relative_C_excursion']) for r in seeds):.3%}. "+"Basal-reference arc is the primary site metric; midpoint-current-circumference arc and circumference increments are retained as metric sensitivities. A fixed material site changes neither primary coordinate under uniform stretch.","",
        "## Sampled increment evidence", "",
        "The original chemical diagnostic at cadence 5 contains no >π/4 increment conditional on both R≥.1; its maximum is 0.4074 rad. This does not exclude faster changes. Spring has rare large sampled global-angle increments at cadence .5. They are displayed with actual native field and scalar-vector context, selected as diagnostics rather than medoids or event classes.", "",
        "| f | Resting C₀ | Eligible .5-unit pairs | Samples >π/4 | Largest sampled increment (rad) |", "|---|---:|---:|---:|---:|"]
    for f in (.001,.01):
        for C in (1.,2.,3.):
            r=select(increments,resting_circumference=C,fraction=f)
            text.append(f"| {f:g} | {C:g} | {sum(int(x['pairs_eligible']) for x in r)} | {sum(int(x['count_above_pi_over_4']) for x in r)} | {max(float(x['max_absolute_angle']) for x in r):.4f} |")
    text += ["", "The largest stronger-input examples have endpoint R values 0.474→0.218 (C₀1), 0.140→0.107 (C₀2) and 0.148→0.145 (C₀3); their current-region counts are 1→2, 2→2 and 2→2. Thus a large global-angle increment alone does not identify uninterrupted relocation of a single activity site. Saved primary fields are spaced by 4 units and cannot resolve the detailed subinterval field evolution at scalar cadence .5.","",
        "## Goryachev forced-state qualification", "",
        "At stronger input, the C20 deterministic mean control loses polarity. Its noisy trajectories therefore cannot be presented as clean spontaneous wandering of the original polarized state. Complete native histories, both strengths, vector/magnitude controls, eligibility and local-pattern statistics are preserved in the chemical interim sheet and original report. Strict lag-8000 site endpoints are too sparse for an inferential physical-site curve (1/0/0 eligible pairs for C20/40/60).", "",
        "## Validation and source access", "",
        "The main audit verifies all 81 Spring cases and 64,881 saved profiles against physical moments, counts, source hashes and exact time grids. The supplemental replay independently checks all saved dominant-region centroids (maximum material-coordinate error 2.23×10⁻¹⁶), regenerates all 72 exact noisy forcing paths, and retains first/last-half site and whole-seed increment-CDF analyses. Independent scientific review separately reproduces the core estimates and whole-seed intervals. Frozen pilot failures and later prospective estimator sensitivity gates remain explicit.","",
        "- [Two-sheet caption and integration metadata](supporting/metadata/content.json)",
        "- [Spring methods and limits](SPRING_METHODS_EVIDENCE_LIMITS.md)",
        "- [Core source and statistics audit](supporting/metadata/spring_noise_analysis_audit.json)",
        "- [Centroid, half-window, CDF and forcing audit](supporting/metadata/spring_noise_supplement_audit.json)",
        "- [Early/late endpoint measurements](supporting/data/spring_noise_short_long_endpoint_summary.csv)",
        "- [Selected scalar increment records and original native archive paths](supporting/data/spring_noise_large_increment_native_context.csv)",
        "- [Preserved original chemical interim](figure4_stochastic_revision.png)", ""]
    (OUT/"RESULTS_WITH_SPRING.md").write_text("\n".join(text))


def docs(artifacts,audit,extra,endpoints,observation):
    main1="figure4_with_spring_dynamics";main2="figure4_with_spring_geometry_patterns"
    legend=[
        "Figure 4. Directional memory needs a physical site readout. Memory panels B show both predeclared strengths: dashed f=.001, solid f=.01. Histories and site/pattern metrics use f=.01. Every noisy group retains all 12 independent seeds; shading and bars are pointwise 95% whole-seed intervals (B shades the stronger input). Native units differ across models.",
        "A1–A3: native active-field histories from the predeclared condition medoids: Wave-pinning C20, Holmes C60 and Spring C₀2. Each model uses a zero-based native activity scale; boundary position is material position around the ring.",
        "B1–B3: whole-window activity-direction cosine, conditional on both endpoint relative first-moment magnitudes R≥.1. The physical activity moment M₁(t) is divided by its initial complex value. Activity direction is distinct from cell movement.",
        "C1–C3: dominant half-contrast-site endpoint RMS shift at lag 8000 for chemical rings or 800 for Spring; percentages are eligible endpoint pairs. Spring uses basal-reference arc shift Cbasal×wrapped material-position change, with the dotted uniform-angle reference Cbasal/√12. Long-lag folded shifts can approach circular mixing and are not travel distances. Spring has its own native length axis.",
        "D1–D3: native site-angle cosine and remapping of the same eligible shortest shifts to common circumference C*=60 for chemical rings or 3.04517 for Spring. Chemical panels use lag 8000; Spring uses early lag 4 before late endpoint mixing. Full lag curves and 4/50/800 comparisons remain supporting evidence. This is a geometric coordinate comparison, not a causal decomposition of global activity-vector memory.",
        "E1–E3: observed time fractions with multiple current half-contrast activity regions and weak relative resultant R<.1. Region multiplicity does not establish functional competition. All 6401 scalar samples per seed enter each fraction. Supporting increments retain rare >π/4 Spring changes at cadence .5 and none in chemical traces at cadence 5, conditional on both endpoint R≥.1; these sampled diagnostics do not define an event class or exclude faster changes.",
        "Complete Goryachev results remain supporting evidence: its high-input C20 mean control also loses polarity, precluding an intact-state wandering interpretation. Both declared strengths, all mean/no-input controls and medoids, magnitude and increment distributions, eligibility sensitivities and full/half-window analyses are retained. Methods define Spring's approved stretch/dilution/diffusion carryover, mean-force global-extension closure and finite transported-material noise law."]
    main_claims=["Angular-memory and physical site-shift measurements can show different size trends. Spring makes the forcing dependence explicit: weak-input lag-800 direction cosine rises .584/.832/.958 across sizes despite nonmonotone basal-arc RMS, whereas stronger input gives .019/−.017/.028 with substantial mean relative first-moment magnitude remaining and late folded shifts near the uniform-angle reference. All three main models retain their native units and controls.",
                 "Pairwise coordinate remapping exposes the geometric part of site-angle retention, while local-region and magnitude readouts distinguish observable organizations without classifying abrupt events or inferring competition."]
    sources=[("Original chemical-ring stochastic report",OLD/"RESULTS.md"),("Original noise design",OLD/"DESIGN.md"),
        ("Complete Goryachev/chemical interim figure",OUT/"figure4_stochastic_revision.png"),("Frozen Spring main design",RUN/"manifest.json"),
        ("Paired finite-law estimator sensitivity",META/"spring_noise_bridge_audit.json"),("Spring full analysis audit",META/"spring_noise_analysis_audit.json"),
        ("Spring centroid/half-window/CDF/forcing audit",META/"spring_noise_supplement_audit.json"),("Spring full seed coverage",DATA/"spring_noise_seed_summary.csv"),
        ("Spring methods and limits",OUT/"SPRING_METHODS_EVIDENCE_LIMITS.md"),("Final four-model quantitative evidence",OUT/"RESULTS_WITH_SPRING.md"),("Spring raw full/half-window memory",DATA/"spring_noise_correlations.csv"),
        ("Spring site and stretch metrics",DATA/"spring_noise_site_statistics.csv"),("Spring first/last-half site metrics",DATA/"spring_noise_site_half_windows.csv"),
        ("Spring magnitude-gated increment distributions",DATA/"spring_noise_increment_cdf.csv"),("Spring actual forcing audits",DATA/"spring_noise_realized_input_audits.csv")]
    sources.append(("Spring early/late endpoint summary",DATA/"spring_noise_short_long_endpoint_summary.csv"))
    sources.extend([("Independent scientific review",OUT.parent/"supporting/metadata/INDEPENDENT_REVIEW.md"),
                    ("Independent requirement and evidence audit",OUT.parent/"supporting/metadata/REQUIREMENT_EVIDENCE_AUDIT.json")])
    if (META/"spring_amount_spatial_balance_validation.json").exists():
        sources.extend([("Spring active amount and spatial-balance methods",OUT/"SPRING_AMOUNT_SPATIAL_BALANCE.md"),
                        ("Spring scalar amount/profile audit",META/"spring_amount_spatial_balance_validation.json"),
                        ("Spring active amount and balance condition summaries",DATA/"spring_amount_spatial_balance/condition_window_summaries.csv"),
                        ("All81 Spring scalar case exports",DATA/"spring_amount_spatial_balance/case_index.csv")])
    supplementary=[("Spring: every native medoid","figure4_spring_all_histories"),("Spring magnitude, geometry and controls","figure4_spring_magnitude_geometry_controls"),
        ("Spring full/half-window memory","figure4_spring_memory_windows"),("Spring increment tails","figure4_spring_increment_tails"),
        ("Spring short/long-lag geometry and circular mixing","figure4_spring_endpoint_lag_sensitivity"),
        ("Spring largest sampled increment native context","figure4_spring_large_increment_context"),
        ("All chemical models including Goryachev","figure4_stochastic_revision"),("Chemical memory: both strengths","figure4_both_strengths_memory"),
        ("Chemical magnitude and controls","figure4_magnitude_controls"),("Chemical increment tails","figure4_increment_tails"),("Chemical site-gate sensitivity","figure4_site_gate_sensitivity")]
    if (OUT/"figure4_spring_amount_spatial_balance.png").exists():
        supplementary.append(("Spring active amount and spatial first harmonic","figure4_spring_amount_spatial_balance"))
    content=dict(title="Figure 4 · Directional memory and physical site dynamics",role="stochastic dynamics",claim=main_claims[0],narrative=main_claims,
        assets=[dict(kind="main",path=str(OUT/(stem+".png")),pdf=str(VECTORS/(stem+".pdf")),svg=str(VECTORS/(stem+".svg")),caption=cap)
                for stem,cap in ((main1,"Native histories, activity-direction memory and physical/reference-arc site shifts."),(main2,"Paired geometric remapping and passive local-pattern/magnitude readouts."))],
        main_asset=str(OUT/(main1+".png")),legend=legend,
        measurements=["Activity-direction correlation and full-vector memory","Physical/reference-arc site RMS with endpoint coverage","Native/common-reference site-angle retention on identical pairs","Local-region and weak-resultant occupancy","Gated increment distributions and full/half-window sensitivity"],
        limitations=["Native units are not calibrated across models.","Activity direction does not assign movement.","Endpoint pairs do not establish site ancestry.","Spring reference-arc shifts separate material change from uniform stretch.","Remapping is a coordinate comparison, not a causal decomposition.","Finite-knot/noise-law and finite-window limitations remain explicit."],
        sources=[dict(label=label,path=str(path.relative_to(ROOT))) for label,path in sources],
        supplementary=[dict(title=label,path=str(OUT/(stem+".png")),caption=label) for label,stem in supplementary]+
            [dict(title=label,path=str(OLD/file),caption=label) for label,file in
             (("Original chemical medoids, all18 conditions","vision_noise_selected_histories.png"),
              ("Goryachev lowest-resultant diagnostic","vision_noise_lowest_resultant_diagnostic.png"))])
    (META/"content.json").write_text(json.dumps(content,indent=2))
    (META/"CAPTION_WITH_SPRING.md").write_text("# Figure 4 revised caption\n\n"+"\n\n".join(legend)+"\n")
    report=dict(generated_at=datetime.now(timezone.utc).isoformat(),renderer=str(Path(__file__)),renderer_sha256=sha(__file__),analysis_only=True,
        original_chemical_cases=243,new_spring_cases=81,total_noisy_cases=216+72,
        spring_full_audit_sha256=sha(META/"spring_noise_analysis_audit.json"),spring_supplement_audit_sha256=sha(META/"spring_noise_supplement_audit.json"),
        output_hashes=[dict(path=str(p),sha256=sha(p)) for p in artifacts],rendered_inspection="pending visual inspection")
    (META/"figure4_with_spring_render_validation.json").write_text(json.dumps(report,indent=2))


def main():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.65,
                         "lines.linewidth":1.4,"xtick.labelsize":8,"ytick.labelsize":8,"svg.fonttype":"none","pdf.fonttype":42})
    d,e,o,r,a,x=inputs();endpoint_window_summary();artifacts=main_sheets(d,e,o,r);artifacts+=spring_support();artifacts+=increment_native_records();evidence_report();docs(artifacts,a,x,e,o)
    print(json.dumps(dict(artifacts=len(artifacts),output=str(OUT),spring_cases=a["cases"])))


if __name__=="__main__":main()
