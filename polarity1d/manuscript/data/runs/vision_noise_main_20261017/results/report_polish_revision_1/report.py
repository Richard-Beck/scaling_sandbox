"""Render/report the complete frozen replicated noisy-ring experiment.

This source selects no input, launches no simulation and changes no estimator.
Full-resolution exports remain authoritative; render-only thinning reduces
plotting memory while preserving every curve's prescribed endpoint.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import os
import subprocess
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MODELS=("wave_pinning","goryachev","holmes_model3")
NAMES=dict(zip(MODELS,("Wave-pinning","Goryachev","Holmes 3")))
COLORS={20.:"#2f73ad",40.:"#42865c",60.:"#c66b25"}
PREFIX="vision_noise_"


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(2**20),b""):digest.update(block)
    return digest.hexdigest()


def read_rows(data,name,predicate=None):
    with (data/(PREFIX+name+".csv")).open(newline="") as stream:
        return [row for row in csv.DictReader(stream) if predicate is None or predicate(row)]


def select(rows,**query):
    return [r for r in rows if all(r[k]==v if isinstance(v,str) else float(r[k])==float(v) for k,v in query.items())]


def array(rows,key):return np.asarray([float(r[key]) for r in rows])

def figure_save(fig,out,stem,artifacts):
    vector=out/"supporting/vectors";vector.mkdir(parents=True,exist_ok=True)
    for extension in ("png","pdf","svg"):
        path=(out if extension=="png" else vector)/(PREFIX+stem+"."+extension)
        fig.savefig(path,dpi=180 if extension=="png" else None,bbox_inches="tight");artifacts.append(path)
    plt.close(fig)


def common_limits(values,include=(0.,1.),padding=.05):
    values=np.asarray(list(values)+list(include),float);values=values[np.isfinite(values)]
    lower=float(values.min());upper=float(values.max());margin=padding*max(upper-lower,.2)
    return lower-margin,upper+margin


def legend(fig,extra=(),bottom=-.01):
    handles=[Line2D([],[],color=COLORS[C],label=f"C = {C:g}") for C in COLORS]+list(extra)
    fig.legend(handles=handles,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.5,bottom))


def trace_plot(rows,metric,ylabel,out,stem,artifacts):
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True,sharey=True)
    displayed=select(rows,metric=metric)
    limits=common_limits([float(r[k]) for r in displayed for k in ("mean","ci_lower","ci_upper")])
    for i,model in enumerate(MODELS):
        for j,f in enumerate((.001,.01)):
            ax=axes[i,j]
            for C,color in COLORS.items():
                noise=select(displayed,model=model,circumference=C,condition="noise",fraction=f)
                mean=select(displayed,model=model,circumference=C,condition="mean",fraction=f)
                t=array(noise,"time")
                ax.fill_between(t,array(noise,"ci_lower"),array(noise,"ci_upper"),color=color,alpha=.15,lw=0)
                ax.plot(t,array(noise,"mean"),color=color,lw=1.3)
                ax.plot(array(mean,"time"),array(mean,"mean"),color=color,lw=.85,ls="--")
            ax.set_ylim(*limits);ax.axhline(0,color=".75",lw=.5);ax.axhline(1,color=".85",lw=.5)
            ax.set_title(f"{NAMES[model]}  |  f = {f:g}",loc="left")
            if j==0:ax.set_ylabel(ylabel)
            if i==2:ax.set_xlabel("Native time")
    legend(fig,[Line2D([],[],color=".2",label="Noise: mean + 95% seed-bootstrap interval"),Line2D([],[],color=".2",ls="--",label="Deterministic mean control")])
    fig.suptitle("12 independent realizations per model, circumference and input strength",y=1.005)
    fig.tight_layout(rect=(0,.065,1,1));figure_save(fig,out,stem,artifacts)


def correlation_plot(corr,direction,out,artifacts):
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True,sharey=True)
    limits=common_limits([float(r[k]) for r in corr+direction for k in ("correlation","direction_correlation","ci_lower","ci_upper") if k in r],include=(-.05,1.))
    for i,model in enumerate(MODELS):
        for j,f in enumerate((.001,.01)):
            ax=axes[i,j]
            for C,color in COLORS.items():
                for source,key,style in ((corr,"correlation","-"),(direction,"direction_correlation","--")):
                    rows=select(source,model=model,circumference=C,condition="noise",fraction=f)
                    t=array(rows,"lag")
                    ax.plot(t,array(rows,key),color=color,ls=style,lw=1.2 if style=="-" else .9)
                    if style=="-":ax.fill_between(t,array(rows,"ci_lower"),array(rows,"ci_upper"),color=color,alpha=.12,lw=0)
            ax.set_ylim(*limits);ax.axhline(0,color=".8",lw=.5);ax.axhline(1/np.e,color=".7",lw=.5,ls=":")
            ax.set_title(f"{NAMES[model]}  |  f = {f:g}",loc="left")
            if j==0:ax.set_ylabel("Lag correlation")
            if i==2:ax.set_xlabel("Lag (native time)")
    legend(fig,[Line2D([],[],color=".2",label="Full vector + 95% seed interval"),Line2D([],[],color=".2",ls="--",label="Direction, both R ≥ 0.1; coverage below")])
    fig.suptitle("Whole-window raw memory; weak vectors remain in the full-vector statistic",y=1.005)
    fig.tight_layout(rect=(0,.065,1,1));figure_save(fig,out,"correlations",artifacts)


def coverage_plot(direction,site,out,artifacts):
    fig,axes=plt.subplots(3,3,figsize=(11,8),sharex=True)
    physical=select(site,metric="shortest_physical_displacement_squared")
    ylim=common_limits([np.sqrt(float(r[k])) for r in physical if r["model"]!="goryachev" for k in ("value","ci_lower","ci_upper") if np.isfinite(float(r[k]))],include=(0.,.1))
    for i,model in enumerate(MODELS):
        for C,color in COLORS.items():
            for f,style in ((.001,"--"),(.01,"-")):
                rows=select(direction,model=model,circumference=C,condition="noise",fraction=f)
                axes[i,0].plot(array(rows,"lag"),array(rows,"eligible_pair_fraction"),color=color,ls=style,lw=1.1)
                rows=select(physical,model=model,circumference=C,condition="noise",fraction=f)
                if model!="goryachev":
                    axes[i,1].plot(array(rows,"lag"),np.sqrt(array(rows,"value")),color=color,ls=style,lw=1.1)
                    axes[i,1].fill_between(array(rows,"lag"),np.sqrt(array(rows,"ci_lower")),np.sqrt(array(rows,"ci_upper")),color=color,alpha=.10,lw=0)
                axes[i,2].plot(array(rows,"lag"),array(rows,"eligible_pair_fraction"),color=color,ls=style,lw=1.1)
        axes[i,0].set_ylabel(NAMES[model]+"\nEligible fraction")
        axes[i,1].set_ylabel("Physical RMS endpoint displacement")
        axes[i,2].set_ylabel("Eligible fraction")
        axes[i,0].set_ylim(-.02,1.04);axes[i,2].set_ylim(-.02,1.04);axes[i,1].set_ylim(*ylim)
        for ax in axes[i]:
            ax.set_xlim(0,8000)
            if i==2:ax.set_xlabel("Lag (native time)")
    axes[1,1].text(.05,.75,"Conservative peak gate leaves sparse\nor no eligible Goryachev endpoints.\nHigh input at lag 8000: C20 has one\npair; C40/C60 have none.\nNo population drift curve shown;\nraw summaries remain in the CSV.",transform=axes[1,1].transAxes,fontsize=7.5,va="top")
    for ax,title in zip(axes[0],("Direction pairs, R ≥ 0.1","Resolved-region endpoint displacement","Resolved-region endpoint pairs")):ax.set_title(title,fontsize=9)
    legend(fig,[Line2D([],[],color=".2",ls="--",label="f = 0.001"),Line2D([],[],color=".2",label="f = 0.01")])
    fig.suptitle("Scalar cadence 5; finite-cadence endpoint comparisons do not prove site identity",y=1.005)
    fig.tight_layout(rect=(0,.065,1,1));figure_save(fig,out,"coverage_geometry",artifacts)


def integral_plot(rows,out,artifacts):
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True,sharey=True)
    for i,model in enumerate(MODELS):
        for j,f in enumerate((.001,.01)):
            ax=axes[i,j]
            for window,marker,offset in (("full","o",-1.),("first_half","^",0.),("last_half","s",1.)):
                chosen=sorted(select(rows,model=model,condition="noise",fraction=f,window=window),key=lambda r:float(r["circumference"]))
                x=array(chosen,"circumference")+offset;y=array(chosen,"restricted_integral")
                error=np.vstack([y-array(chosen,"ci_lower"),array(chosen,"ci_upper")-y])
                # The original estimate need not lie inside a finite-bootstrap percentile interval.
                for k,row in enumerate(chosen):
                    C=float(row["circumference"])
                    ax.plot([x[k],x[k]],[float(row["ci_lower"]),float(row["ci_upper"])],color=COLORS[C],lw=.9)
                    ax.scatter([x[k]],[y[k]],color=COLORS[C],marker=marker,s=22,zorder=3,facecolors="none" if window=="first_half" else COLORS[C])
            ax.set_title(f"{NAMES[model]}  |  f = {f:g}",loc="left");ax.axhline(0,color=".8",lw=.5)
            if j==0:ax.set_ylabel("Signed restricted area, lag 0–8000\n(native time units)")
            if i==2:ax.set_xlabel("Circumference C")
            ax.set_xticks((20,40,60))
    legend(fig,[Line2D([],[],color=".2",marker=m,ls="",label=w) for m,w in (("o","Full window"),("^","First half"),("s","Last half"))])
    fig.suptitle("Window dependence is retained; no burn-in, forced exponential fit or extrapolation",y=1.005)
    fig.tight_layout(rect=(0,.065,1,1));figure_save(fig,out,"restricted_integrals",artifacts)


def history_plot(selected,out,artifacts):
    rows=[r for r in selected if r["role"]=="within_condition_descriptor_medoid"]
    assert len(rows)==18
    fig,axes=plt.subplots(3,6,figsize=(15,8),sharey=True)
    for i,model in enumerate(MODELS):
        cases=[]
        for f in (.001,.01):
            for C in COLORS:
                row,=select(rows,model=model,circumference=C,fraction=f)
                with np.load(row["profile_archive"]) as z:
                    time=z["time"].copy();primary=z["primary"].copy();x=z["x"].copy()
                position=(x/C+.5)%1-.5;order=np.argsort(position)
                cases.append((row,time,position[order],primary[:,order]))
        vmin=min(c[3].min() for c in cases);vmax=max(c[3].max() for c in cases)
        for j,(row,time,position,primary) in enumerate(cases):
            ax=axes[i,j]
            mesh=ax.pcolormesh(position,time,primary,shading="nearest",vmin=vmin,vmax=vmax,cmap="viridis",rasterized=True)
            ax.set_title(f"C{float(row['circumference']):g}, f={float(row['fraction']):g}\nseed {row['seed']}",fontsize=7)
            ax.set_xlim(-.5,.5);ax.set_ylim(32000,0)
            if j==0:ax.set_ylabel(NAMES[model]+"\nNative time")
            if i==2:ax.set_xlabel("Position / C")
        fig.colorbar(mesh,ax=axes[i,:].tolist(),fraction=.02,pad=.01,label="Native active field")
    fig.subplots_adjust(left=.065,right=.87,bottom=.07,top=.92,wspace=.17,hspace=.3)
    fig.suptitle("Every prespecified condition medoid; one raw native colour scale per model",y=.985)
    figure_save(fig,out,"selected_histories",artifacts)


def weakest_diagnostic(selected,run,manifest,out,artifacts):
    row,=select(selected,role="lowest_sampled_resultant_diagnostic")
    with np.load(row["archive"]) as trace:
        t=trace["time"].copy();P=trace["P_real"]+1j*trace["P_imag"];R=abs(P)
        q1=trace["q1"].copy();q2=trace["q2"].copy();contrast=trace["contrast"].copy()
    j=int(np.argmin(R));event=float(t[j]);C=float(row["circumference"])
    with np.load(row["profile_archive"]) as z:pt=z["time"].copy();primary=z["primary"].copy();x=z["x"].copy()
    k=int(np.argmin(abs(pt-event)))
    mean_job,=[job for job in manifest["jobs"] if job["model"]==row["model"] and job["circumference"]==C and job["mean_fraction"]==float(row["fraction"]) and job["kind"]=="mean"]
    with np.load(run/"workers"/f"{mean_job['slot']:04d}"/"primary.npz") as z:mean_primary=z["primary"][k].copy()
    position=(x/C+.5)%1-.5;order=np.argsort(position)
    fig,axes=plt.subplots(1,3,figsize=(12,3.4))
    span=(t>=max(0,event-2000))&(t<=min(32000,event+2000))
    axes[0].plot(t[span],R[span],label="R");axes[0].plot(t[span],q1[span],label="q₁");axes[0].plot(t[span],q2[span],label="q₂")
    axes[0].axvline(event,color=".5",ls=":",lw=.8);axes[0].set_xlabel("Native time");axes[0].set_ylabel("Resultant / harmonic ratio");axes[0].legend(frameon=False,fontsize=8)
    axes[1].plot(position[order],primary[k,order],label="Noisy state");axes[1].plot(position[order],mean_primary[order],ls="--",label="Mean control")
    axes[1].set_xlabel("Position / C");axes[1].set_ylabel("Native active field");axes[1].legend(frameon=False,fontsize=8)
    window=(pt>=max(0,event-2000))&(pt<=min(32000,event+2000))
    mesh=axes[2].pcolormesh(position[order],pt[window],primary[window][:,order],shading="nearest",cmap="viridis",rasterized=True)
    axes[2].invert_yaxis();axes[2].set_xlabel("Position / C");axes[2].set_ylabel("Native time")
    fig.colorbar(mesh,ax=axes[2],fraction=.045,label="Native active field")
    fig.suptitle(f"Lowest sampled noisy R: {NAMES[row['model']]} C{C:g}, f={float(row['fraction']):g}, seed {row['seed']}\nscalar event t={event:g}; nearest saved profile t={pt[k]:g}; diagnostic, not a typical realization",fontsize=10)
    fig.tight_layout(rect=(0,0,1,.88));figure_save(fig,out,"lowest_resultant_diagnostic",artifacts)
    return dict(model=row["model"],circumference=C,fraction=float(row["fraction"]),seed=row["seed"],
                event_time=event,R=float(R[j]),q1=float(q1[j]),q2=float(q2[j]),contrast=float(contrast[j]),profile_time=float(pt[k]))


def input_audit_export(run,manifest,data):
    rows=[];sources=[]
    for job in manifest["jobs"]:
        if job["kind"]!="noise":continue
        folder=run/"workers"/f"{job['slot']:04d}";path=folder/"noise_audit.json";mode=folder/"noise_modes.csv"
        audit=json.loads(path.read_text());mu=float(job["mean"])
        rows.append(dict(model=job["model"],circumference=job["circumference"],fraction=job["mean_fraction"],seed=job["seed"],
                         local_mean_over_mu=audit["empirical_input_mean"]/mu,
                         local_second_moment_over_mu_squared=audit["empirical_input_second_moment"]/mu**2,
                         local_variance_over_mu_squared=(audit["empirical_input_second_moment"]-audit["empirical_input_mean"]**2)/mu**2,
                         expected_local_variance_over_mu_squared=np.e-1,
                         gaussian_spatial_covariance_ell=audit["empirical_spatial_covariance_at_ell"],
                         expected_gaussian_spatial_covariance_ell=audit["theoretical_spatial_covariance_at_ell"],
                         gaussian_temporal_covariance_tau=audit["empirical_temporal_covariance_at_tau"],
                         expected_gaussian_temporal_covariance_tau=audit["theoretical_temporal_covariance_at_tau"],
                         input_spatial_covariance_ell_over_mu_squared=audit["empirical_input_spatial_covariance_at_ell"]/mu**2,
                         expected_input_spatial_covariance_ell_over_mu_squared=audit["theoretical_input_spatial_covariance_at_ell"]/mu**2,
                         input_temporal_covariance_tau_over_mu_squared=audit["empirical_input_temporal_covariance_at_tau"]/mu**2,
                         expected_input_temporal_covariance_tau_over_mu_squared=audit["theoretical_input_temporal_covariance_at_tau"]/mu**2,
                         integrated_input_mean=audit["empirical_integrated_input_mean"],
                         expected_integrated_input_mean=audit["theoretical_integrated_input_mean"],
                         integrated_input_variance_over_mu_squared=audit["empirical_integrated_input_variance"]/mu**2,
                         expected_integrated_input_variance_over_mu_squared=audit["theoretical_integrated_input_variance"]/mu**2,
                         first_input_mode_second_moment_over_mu_squared=audit["empirical_first_angular_input_mode_second_moment"]/mu**2,
                         expected_first_input_mode_second_moment_over_mu_squared=audit["theoretical_first_angular_input_mode_second_moment"]/mu**2,
                         noise_audit_path=str(path),noise_audit_sha256=sha(path),noise_modes_path=str(mode),noise_modes_sha256=sha(mode)))
        sources.extend((path,mode))
    assert len(rows)==216
    path=data/(PREFIX+"input_audits.csv")
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    lines=[]
    for model in MODELS:
        for f in (.001,.01):
            for C in COLORS:
                group=select(rows,model=model,circumference=C,fraction=f);assert len(group)==12
                mean=np.mean([float(r["local_mean_over_mu"]) for r in group]);variance=np.mean([float(r["local_variance_over_mu_squared"]) for r in group])
                zero=np.mean([float(r["integrated_input_variance_over_mu_squared"]) for r in group]);first=np.mean([float(r["first_input_mode_second_moment_over_mu_squared"]) for r in group])
                lines.append(f"| {NAMES[model]} | {C:g} | {f:g} | {mean:.4g} | {variance:.4g} | {zero:.4g} / {float(group[0]['expected_integrated_input_variance_over_mu_squared']):.4g} | {first:.4g} / {float(group[0]['expected_first_input_mode_second_moment_over_mu_squared']):.4g} |")
    return path,sources,"\n".join(lines)


def fmt(value):return "unavailable" if not np.isfinite(float(value)) else f"{float(value):.4g}"


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--run",type=Path,required=True);parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--reuse-figures",action="store_true");parser.add_argument("--message",type=Path)
    args=parser.parse_args();run=args.run.resolve();out=args.output.resolve();root=run.parents[2];repo=root.parents[1]
    data=out/"supporting/data";metadata=out/"supporting/metadata"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    audit=json.loads((metadata/(PREFIX+"analysis_audit.json")).read_text());manifest=json.loads((run/"manifest.json").read_text())
    assert audit["cases"]==243 and audit["groups"]==45 and audit["maxlag_native_time"]==8000 and audit["site_descriptor_cadence"]=="scalar"
    assert audit["bootstrap_draws"]==1000 and audit["bootstrap_seed"]==920731
    traces=read_rows(data,"ensemble_traces",lambda r:r["metric"] in ("initial_memory","R") and float(r["time"])%20==0)
    corr=read_rows(data,"lag_correlations",lambda r:r["window"]=="full" and r["condition"] in ("noise","mean") and float(r["lag"])%20==0)
    direction=read_rows(data,"direction_correlations",lambda r:r["window"]=="full" and r["condition"]=="noise" and float(r["gate"])==.1 and float(r["lag"])%20==0)
    persistence=read_rows(data,"persistence")
    site=read_rows(data,"site_displacement",lambda r:r["window"]=="full" and r["condition"]=="noise" and float(r["lag"])%20==0)
    selected=read_rows(data,"representative_histories");summaries=read_rows(data,"seed_summary")
    assert len(summaries)==243 and all(int(r["n_seeds"])==12 for r in persistence if r["condition"]=="noise")
    artifacts=[]
    input_path,input_sources,input_table=input_audit_export(run,manifest,data);artifacts.append(input_path)
    if args.reuse_figures:
        # Preserve the reviewed native-history choices and all unchanged panels.
        previous=json.loads((metadata/(PREFIX+"main_report_manifest.json")).read_text())
        diagnostic=previous["diagnostic"]
        for stem in ("initial_memory","magnitude","restricted_integrals","selected_histories","lowest_resultant_diagnostic"):
            paths=[out/(PREFIX+stem+".png")]+[out/"supporting/vectors"/(PREFIX+stem+"."+ext) for ext in ("pdf","svg")]
            assert all(path.is_file() for path in paths);artifacts.extend(paths)
        coverage_plot(direction,site,out,artifacts)
        correlation_plot(corr,direction,out,artifacts)
    else:
        trace_plot(traces,"initial_memory","Initial-axis memory: mean Re P",out,"initial_memory",artifacts)
        trace_plot(traces,"R","Mean relative magnitude |P|",out,"magnitude",artifacts)
        correlation_plot(corr,direction,out,artifacts);coverage_plot(direction,site,out,artifacts)
        integral_plot(persistence,out,artifacts);history_plot(selected,out,artifacts)
        diagnostic=weakest_diagnostic(selected,run,manifest,out,artifacts)
    table=[];mean_table=[]
    for model in MODELS:
        for f in (.001,.01):
            for C in COLORS:
                initial,=select(traces,model=model,circumference=C,condition="noise",fraction=f,metric="initial_memory",time=32000)
                pers,=select(persistence,model=model,circumference=C,condition="noise",fraction=f,window="full")
                d,=select(direction,model=model,circumference=C,fraction=f,lag=8000)
                table.append(f"| {NAMES[model]} | {C:g} | {f:g} | {fmt(initial['mean'])} [{fmt(initial['ci_lower'])}, {fmt(initial['ci_upper'])}] | {fmt(pers['restricted_integral'])} [{fmt(pers['ci_lower'])}, {fmt(pers['ci_upper'])}] | {fmt(d['direction_correlation'])} | {fmt(d['eligible_pair_fraction'])} |")
                control,=select(summaries,model=model,circumference=C,condition="mean",fraction=f)
                mean_table.append(f"| {NAMES[model]} | {C:g} | {f:g} | {fmt(control['final_R'])} | {fmt(control['final_total_active_ratio'])} |")
    status=audit["manifest_audit"]["statuses"]
    sizes={str(path):int(subprocess.check_output(["du","-s","-B1",str(path)],text=True).split()[0]) for path in (run,out)}
    actual_message=args.message.read_text().strip() if args.message else "Actual numerical summaries appear below; final interpretive prose awaits direct review of the completed curves and uncertainty."
    text=f'''# Replicated noisy rings: polarity memory, magnitude and local geometry

**All 243 prescribed integrations completed and all 216 noisy realizations are
included.** The experiment contains 12 independent realizations in each of the
18 model/circumference/strength groups, plus 18 deterministic mean and nine
zero-input controls. Every realization reaches 32000 native time units. This
report supplies the noise component of the manuscript mechanism-and-measurement
checkpoint; it does not assign a cell velocity or calibrate migration persistence.
Generated {datetime.now(timezone.utc).isoformat()}.

**Effective executed protocol:** $T=32000$; scalar and site-descriptor cadence 5;
saved-primary-profile cadence 40; shared maximum lag 8000. The operative
[frozen configuration](../../../../data/runs/vision_noise_main_20261017/configuration.json)
and [frozen analysis plan](../../../../data/runs/vision_noise_main_20261017/plan_freeze.json)
supersede prospective timing/cadence wording in the design-stage notes.
The [frozen scientific design](../../../../data/runs/vision_noise_main_20261017/protocol/scientific_design.md)
records the prospective rationale; it is not a live producer-status page.

{actual_message}

## Question and declared comparison

Does enlargement change measurable whole-boundary polarity memory under local
fluctuating activation? Wave-pinning, Goryachev and Holmes begin from their exact
qualifying C20/C40/C60 basals at spatial step 0.05. The positive input law,
preparation, kinetic parameters, two strengths and 12 independent seeds per
condition were frozen before the main producer launched. Both strengths are
retained, including conditions in which the deterministic mean alters the
polarized state.

The native external activation input is

$$
I(s,t)=fA_m\\exp[Z(s,t)-1/2],\\quad
f\\in\\{{0.001,0.01\\}},\\quad
A_m\\in\\{{0.1,6.4,1\\}}.\\tag{{1}}
$$

The Gaussian field has unit local ensemble variance, physical spatial
covariance width 1.5 and OU knot correlation time 50. Knots are separated by
5 and use analytically variance-normalized linear interpolation. No
per-realization mean removal, amplitude normalization or state clipping is
introduced. This is the declared interpolated driving law, not an exact
continuous-time OU path. Local noise statistics are fixed within each model;
total input and whole-ring angular input modes change with circumference.
The three native input/time scales are not calibrated equal biological exposures.

## Realized forcing and whole-ring modes

The [216 realized-input audits](supporting/data/{PREFIX}input_audits.csv)
link every worker noise-audit JSON and per-knot mode trace with exact hashes.
They retain realized local means/variances, Gaussian and lognormal spatial and
temporal covariances, integrated zeroth-mode variance and first-angular-mode
second moment alongside finite-ring expectations. The table averages the 12
realized audit summaries within each declared condition; it never rescales
any input realization. Expected normalized local mean and variance are 1 and
$e-1$, respectively. Zeroth-mode variance and first-mode second moment below
are divided by $\\mu^2$; each entry is observed / theoretical.

| Model | C | f | Local mean / μ | Local variance / μ² | Zeroth-mode variance / μ² | First-mode second moment / μ² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{input_table}

Thus equal local covariance is not equal whole-ring angular forcing or equal
total input at different circumferences. The complete technical archive verifier
regenerated every realization from its frozen source/seed; finite realized
moment deviations remain observations rather than normalization targets.

## Initial memory and magnitude

![Replicated initial-axis memory]({PREFIX}initial_memory.png)

Solid curves are the independent-realization mean; shading is the 2.5–97.5%
percentile interval from 1000 complete-realization bootstrap resamples. Dashed
curves are deterministic mean-input controls. The same full range is used
across panels and contains every displayed mean and interval. No individual
seed or strength is removed because it gives an inconvenient result.

![Replicated vector magnitude]({PREFIX}magnitude.png)

For the native active field $u(s,t)$, the measured moments and aligned vector are

$$
M_k(t)=\\int_0^C u(s,t)e^{{2\\pi i k s/C}}\\,\\mathrm ds,
\\qquad P(t)=\\frac{{M_1(t)}}{{M_1(0)}},\\qquad
q_k(t)=\\frac{{|M_k(t)|}}{{M_0(t)}}.\\tag{{2}}
$$

The activity vector is initially $(1,0)$ and its relative magnitude is $R=|P|$;
R can exceed one because its denominator is the initial moment. Initial-axis memory is $\\langle\\operatorname{{Re}}P(t)\\rangle$;
it can change through direction, magnitude or both. A small resultant can
coexist with active opposing sites and cannot by itself establish homogeneous
local activity. Total activity, $q_1$, $q_2$, and peak/region diagnostics remain
in the seedwise tables. Uniform mean input can act on a spatially nonuniform
inactive field, so its direct contribution to the active moment need not vanish.

## Lag memory and direction

![Replicated vector and direction lag correlations]({PREFIX}correlations.png)

The solid statistic is the raw uncentered complex-vector lag correlation: the
ratio of ensemble-averaged within-realization endpoint cross-products to the
square root of both endpoint energies. It retains weak resultants and subtracts
no temporal or ensemble mean vector. Each dashed curve is the conditional
unit-direction cosine with both endpoint magnitudes at least 0.1. These answer
different questions; an amplitude-memory change is not automatically slower
turning. Confidence bands use whole seeds, never time origins as independent
replicates. Complete direction intervals and magnitude gates 0.05/0.2 are in
[the direction table](supporting/data/{PREFIX}direction_correlations.csv).

![Direction eligibility and physical geometry]({PREFIX}coverage_geometry.png)

Direction-pair coverage is shown explicitly. Physical displacement uses producer
scalar descriptors every 5 units, independently validated against all saved
40-unit primary profiles. Both endpoints must have one current half-contrast
region, at most one peak with native prominence at least 0.05, contrast at
least 0.05 and R at least 0.1. Missing eligible pairs remain unavailable. The
same physical-displacement scale is used throughout; angular persistence can
increase merely because the same physical displacement subtends a smaller
angle on a larger ring. These finite-cadence endpoint comparisons do not prove
continued site identity, define a diffusion coefficient or describe cell travel.
The tables preserve all eligibility fractions and bootstrap-defined fractions.
The absolute prominence gate is deliberately conservative: small secondary
peaks can fail it while a clear dominant half-contrast incumbent remains.
Missing eligible pairs do not establish instability, and every detected
0.05-prominence peak is not labelled a functional competitor or protrusion.
Raw dominant-site descriptors and region counts remain available. Goryachev's
physical RMS fragments are not displayed as a population curve: at high input
and lag 8000, C20 has only one eligible pair among 57612, while C40/C60 have
none. The collapsed one-pair interval is not a precision estimate. Complete
computed values and bootstrap-defined fractions remain in the unmodified CSV;
this restriction does not prevent observing its strong dominant incumbent.

## Restricted memory and adaptation

![Full and half-window restricted integrals]({PREFIX}restricted_integrals.png)

The integral is the raw signed area through lag 8000, including negative
correlations. It is an observation-window descriptor in native time units,
not a fitted exponential persistence time. Full/first-half/last-half windows
use the same maximum lag and no chosen burn-in. Differences between them expose
continuing adaptation; agreement alone does not prove stationarity. A curve
that does not cross $1/e$ within this horizon does not supply a resolved
exponential time. Weak endpoint-energy statuses remain explicit in
[the complete integral table](supporting/data/{PREFIX}persistence.csv).

| Model | C | f | Final mean Re P [95% interval] | Restricted area [95% interval] | Direction cosine at lag 8000 | Eligible direction pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
'''+"\n".join(table)+f'''

## Deterministic mean regimes

The mean controls determine whether noise is probing the original kind of
polarized state or a different forced-state regime. In particular the pilot's
high-input Goryachev C20 mean destroyed its polarized resultant. The complete
main controls are retained here rather than counting that regime as a clean
noise-specific directional-persistence comparison.

| Model | C | f | Mean-control final R | Mean-control final total activity / initial |
| --- | ---: | ---: | ---: | ---: |
'''+"\n".join(mean_table)+f'''

## Local histories selected by the frozen rule

![Every condition medoid]({PREFIX}selected_histories.png)

One medoid is shown for every noisy model/C/strength condition. The frozen
descriptors are final initial-axis memory, minimum magnitude, time-mean
magnitude and the profile fraction with multiple current half-contrast regions.
The sum of within-group median deviations, scaled by descriptor IQR with a
unit zero-IQR fallback, determines the medoid; ties use case identifier.
These are descriptive representatives, not the histories most favourable to a
size effect. Each model row has one raw native activity colour scale. The
same angular axis spans different physical lengths at different C. Exact seeds,
archives and hashes are retained in [the selection table](supporting/data/{PREFIX}representative_histories.csv).

![Lowest-resultant diagnostic]({PREFIX}lowest_resultant_diagnostic.png)

The additional diagnostic is the globally lowest sampled noisy resultant,
selected independently of which size trend it suggests. It occurs in
{NAMES[diagnostic['model']]} at C={diagnostic['circumference']:g},
f={diagnostic['fraction']:g}, seed {diagnostic['seed']},
$t={diagnostic['event_time']:g}$: $R={diagnostic['R']:.5g}$,
$q_1={diagnostic['q1']:.5g}$, $q_2={diagnostic['q2']:.5g}$ and native
contrast {diagnostic['contrast']:.5g}. The profile snapshot is the nearest
saved time, $t={diagnostic['profile_time']:g}$; it is not passed off as the
exact scalar-event state. The matched mean control and raw local profile
separate weak whole-ring asymmetry from local activity loss or cancellation.
No dominant-site switch alone is labelled incumbent replacement.

## Corresponding experimental measurement

A whole-boundary activity reporter could provide the same first angular moment,
its magnitude and its temporal correlations. Local site locations, intensities,
region counts and physical displacement would distinguish persistent asymmetry
from local wandering, multisite cancellation or loss of contrast. This is a
polarity-memory readout motivated by the enlargement/persistence hypothesis;
it is distinct from measured migration trajectories and conditional OU velocity
persistence estimates in the manuscript's empirical component.

## Complete execution and numerical evidence

The independent analysis verified {audit['manifest_audit']['frozen_source_count']}
frozen source files, {audit['manifest_audit']['external_count']} original
provenance/input files and {audit['manifest_audit']['worker_output_hashes_verified']}
worker-output hashes. Every exact basal start is checked. At each saved primary profile, raw moments
and periodic half-contrast region/peak descriptors are independently reproduced;
scalar descriptors supply the intervening five-unit site series.
Maximum retained raw-moment discrepancy is
{max(r['moment_error'] for r in audit['sources']):.3g}; maximum conserved-pool
relative drift is {max(r['max_pool_relative_error'] for r in status):.3g}.
All 243 complete statuses and the prescribed 12 independent seeds per group
are required; a failed or missing trajectory cannot be silently censored.
The signed integrals, full raw curves, independent-seed intervals and weak
vector/pair statuses are retained in [the analysis audit](supporting/metadata/{PREFIX}analysis_audit.json)
and [the seed summaries](supporting/data/{PREFIX}seed_summary.csv).

A [separate numerical amendment](../../../../data/runs/vision_noise_main_20261017/results/analysis_numerical_amendment/AMENDMENT.md)
preserves the original frozen estimator while replacing cancellation-prone tail
subtraction and weak FFT products by positive suffix sums and direct endpoint
products. It leaves the producer, all trajectories and scientific definitions
unchanged; no normalized correlation is clipped. All 1601 lags in all three
windows of all 27 actual mean/zero controls pass independent direct-product
regression, with maximum absolute correlation error below $1.281\\times10^{{-12}}$.
Extreme tiny-tail/large-range fixtures also pass, with genuinely unresolved
endpoint energies left undefined. A separate
[center-observability adapter amendment](../../../../data/runs/vision_noise_main_20261017/results/statistics_centroid_amendment_20261017/AMENDMENT.md)
retains all field/count checks but compares positions only when a current
component exists. A flat control otherwise carries its previous scalar- versus
profile-cadence position, which is not an observed center. Such carried labels are
flagged as unobserved current centers in the validation audit and cannot meet
the unchanged direction/site gates.
All 243 cases and 194643 saved profiles pass the actual strict adapter validation;
no center tolerance or mathematical readout is changed. Seed/group tables retain numerical fallback
and sensitivity fields alongside weak-state statuses.

The frozen pilot numerical checks distinguish fixed-forcing ODE accuracy from
sensitivity to the interpolated input path. Tightening the ODE solver changes
P by at most $7.71\\times10^{{-8}}$ over the three checked cases. Coupling the
5-to-2.5-unit interpolation changes maximum P by 0.00546/0.00176/0.01807 in
Wave-pinning/Goryachev/Holmes, with raw-profile differences of
2.31%/0.219%/12.74% of basal contrast. Those are finite-pilot path sensitivity
checks, not a theorem of long-run stochastic convergence or a claim that
local histories are numerically identical. The [reviewed pilot](../vision_noise_pilot_20261017/RESULTS.md)
retains every control, check and interpretation limit.

Allocated storage measured after generation: native main group
{sizes[str(run)]:,} bytes and this report group {sizes[str(out)]:,} bytes before
this final prose/manifest. The root data-manifest owner receives exact measured
sizes; scientific importance is the complete replicated stochastic comparison.

A fresh-context actual-result critique will assess the curves, local histories,
uncertainty and scientific claims before this report is folded into the final
three-figure checkpoint. Output generation alone does not establish an accepted
manuscript figure or a general enlargement effect.
'''
    path=out/"RESULTS.md";path.write_text(text);artifacts.append(path)
    source_paths=[Path(__file__).resolve(),run/"manifest.json",run/"plan_freeze.json",run/"protocol/analysis_source.py",metadata/(PREFIX+"analysis_audit.json")]
        # Numerical amendment is separate from the original prospective protocol.
    source_paths += [Path(audit["generator_path"]),run/"results/analysis_numerical_amendment/REGRESSION.json",run/"results/analysis_numerical_amendment/AMENDMENT.md",run/"results/statistics_centroid_amendment_20261017/REGRESSION.json",run/"results/statistics_centroid_amendment_20261017/AMENDMENT.md"]
    source_paths += [run/"results/statistics_workflow/report.py",metadata/"FRESH_NUMERICAL_AMENDMENT_REVIEW.md",metadata/"FRESH_FLAT_CENTER_AMENDMENT_REVIEW.md"]
    if args.message:source_paths.append(args.message.resolve())
    source_paths += input_sources
    source_paths += list(data.glob(PREFIX+"*.csv"))
    manifest_path=metadata/(PREFIX+"main_report_manifest.json")
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root="../..",source_sha256={str(p.relative_to(repo)):sha(p) for p in source_paths},
                sha256={str(p.relative_to(out)):sha(p) for p in artifacts},complete_cases=243,noise_replicates=216,noise_groups=18,
                no_integration=True,raw_prescribed_statistics=True,figure_thinning="plot-only times/lags every20; authoritative exports remain5",
                diagnostic=diagnostic,allocated_bytes_before_final_prose=sizes)
    manifest_path.write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(dict(report=str(path),artifacts=len(artifacts),allocated_bytes=sizes),indent=2))

if __name__=="__main__":main()
