"""Render the completed bounded noise pilot; no simulations or tuning."""
from pathlib import Path
from datetime import datetime,timezone
import csv
import hashlib
import json
import os
import sys
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
RUN=ROOT/"data/runs/vision_noise_pilot_20261017"
OUT=ROOT/"figures/length_mechanisms_20261008/dev/vision_noise_pilot_20261017"
DATA=OUT/"supporting/data"
META=OUT/"supporting/metadata"
PREFIX="vision_noise_pilot_"
MODELS=("wave_pinning","goryachev","holmes_model3")
NAMES={"wave_pinning":"Wave-pinning","goryachev":"Goryachev","holmes_model3":"Holmes 3"}
COLORS={20.:"#2f73ad",60.:"#c66b25"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(name):
    with (DATA/(PREFIX+name+".csv")).open() as stream:
        return list(csv.DictReader(stream))


def subset(source,**conditions):
    return [r for r in source if all(str(r[key])==str(value) if isinstance(value,str)
                                    else float(r[key])==float(value) for key,value in conditions.items())]


def save(fig,stem,artifacts):
    fig.savefig(OUT/(PREFIX+stem+".png"),dpi=180,bbox_inches="tight")
    (OUT/"supporting/vectors").mkdir(parents=True,exist_ok=True)
    for extension in ("pdf","svg"):
        path=OUT/"supporting/vectors"/(PREFIX+stem+"."+extension)
        fig.savefig(path,bbox_inches="tight");artifacts.append(path)
    artifacts.append(OUT/(PREFIX+stem+".png"));plt.close(fig)


def setup():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.labelsize":9,
                         "xtick.labelsize":8,"ytick.labelsize":8,"axes.spines.top":False,
                         "axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})


def dynamics(traces,metric,title,ylim,stem,artifacts):
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True,sharey=True)
    for i,model in enumerate(MODELS):
        for j,fraction in enumerate((.001,.01)):
            ax=axes[i,j]
            for C in (20.,60.):
                for condition,style in (("noise","-"),("mean","--"),("zero",":")):
                    chosen=subset(traces,model=model,circumference=C,condition=condition,
                                  fraction=0. if condition=="zero" else fraction,metric=metric)
                    ax.plot([float(r["time"]) for r in chosen],[float(r["mean"]) for r in chosen],
                            style,color=COLORS[C],lw=1.25 if condition=="noise" else .85,
                            alpha=1 if condition!="zero" else .6)
            ax.set_title(f"{NAMES[model]}  |  f = {fraction:g}",loc="left")
            ax.set_ylim(*ylim);ax.axhline(1,color=".8",lw=.5,zorder=0)
            ax.axhline(0,color=".75",lw=.5,zorder=0)
            if j==0:ax.set_ylabel(title)
            if i==2:ax.set_xlabel("Native time")
    handles=[Line2D([],[],color=COLORS[C],label=f"C = {C:g}") for C in (20.,60.)]
    handles += [Line2D([],[],color=".25",ls=style,label=label) for style,label in
                (("-","One noisy realization"),("--","Deterministic mean"),(":","Zero input"))]
    fig.legend(handles=handles,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.5,-.012))
    fig.suptitle("Bounded pilot: one realization per condition; no replication interval",y=1.005)
    fig.tight_layout(rect=(0,.065,1,1));save(fig,stem,artifacts)


def correlation_plot(correlations,directions,artifacts):
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True,sharey=True)
    for i,model in enumerate(MODELS):
        for j,fraction in enumerate((.001,.01)):
            ax=axes[i,j]
            for C in (20.,60.):
                chosen=subset(correlations,model=model,circumference=C,condition="noise",fraction=fraction,window="full")
                ax.plot([float(r["lag"]) for r in chosen],[float(r["correlation"]) for r in chosen],color=COLORS[C],lw=1.3)
                chosen=subset(directions,model=model,circumference=C,condition="noise",fraction=fraction,window="full",gate=.1)
                ax.plot([float(r["lag"]) for r in chosen],[float(r["direction_correlation"]) for r in chosen],color=COLORS[C],ls="--",lw=.9)
                chosen=subset(correlations,model=model,circumference=C,condition="mean",fraction=fraction,window="full")
                values=[float(r["correlation"]) if r["status"]=="available" else np.nan for r in chosen]
                ax.plot([float(r["lag"]) for r in chosen],values,color=COLORS[C],ls=":",lw=.9)
            ax.set_title(f"{NAMES[model]}  |  f = {fraction:g}",loc="left")
            ax.set_ylim(-.5,1.04);ax.axhline(0,color=".8",lw=.5);ax.axhline(1/np.e,color=".7",lw=.6,ls=":")
            if model=="goryachev" and fraction==.01:
                ax.text(.03,.06,"C20 mean: weak endpoint energies;\nnormalized curve withheld here",transform=ax.transAxes,fontsize=7.5)
            if j==0:ax.set_ylabel("Uncentered lag correlation")
            if i==2:ax.set_xlabel("Lag (native time)")
    handles=[Line2D([],[],color=COLORS[C],label=f"C = {C:g}") for C in (20.,60.)]
    handles += [Line2D([],[],color=".25",ls=s,label=l) for s,l in
                (("-","Noisy vector"),("--","Noisy direction, R ≥ 0.1"),(":","Mean-input vector; weak values withheld"))]
    fig.legend(handles=handles,loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(.5,-.023))
    fig.suptitle("Full-window pilot curves; direction eligibility is reported separately",y=1.005)
    fig.tight_layout(rect=(0,.075,1,1));save(fig,"correlations",artifacts)


def kymographs(summaries,artifacts):
    fig,axes=plt.subplots(3,4,figsize=(12,8),sharey=True)
    for i,model in enumerate(MODELS):
        cases=[]
        for f,C in ((.001,20.),(.001,60.),(.01,20.),(.01,60.)):
            row,=subset(summaries,model=model,circumference=C,condition="noise",fraction=f)
            path=Path(row["archive"]).parent/"primary.npz"
            with np.load(path) as z:
                time=z["time"].copy();primary=z["primary"].copy();x=z["x"].copy()
            cut=((x/C+.5)%1)-.5;order=np.argsort(cut)
            cases.append((f,C,time,cut[order],primary[:,order]))
        vmin=min(case[4].min() for case in cases);vmax=max(case[4].max() for case in cases)
        for j,(f,C,time,position,primary) in enumerate(cases):
            ax=axes[i,j]
            mesh=ax.pcolormesh(position,time,primary,shading="nearest",vmin=vmin,vmax=vmax,rasterized=True,cmap="viridis")
            ax.set_title(f"C{C:g}, f={f:g}",fontsize=9)
            ax.set_xlim(-.5,.5);ax.set_ylim(4000,0)
            ax.axvline(0,color="white",lw=.4,ls=":",alpha=.6)
            if j==0:ax.set_ylabel(NAMES[model]+"\nNative time")
            if i==2:ax.set_xlabel("Position / C (periodic cut)")
        fig.colorbar(mesh,ax=axes[i,:].tolist(),pad=.012,fraction=.025,label="Native active field")
    fig.suptitle("Every noisy pilot history: shared native colour range within each model",y=.98)
    fig.subplots_adjust(left=.08,right=.87,bottom=.07,top=.93,wspace=.2,hspace=.27)
    save(fig,"all_noise_histories",artifacts)


def geometry_diagnostic(directions,site,artifacts):
    fig,axes=plt.subplots(3,3,figsize=(11,8),sharex=True)
    for i,model in enumerate(MODELS):
        for C in (20.,60.):
            for fraction,style in ((.001,"--"),(.01,"-")):
                chosen=subset(directions,model=model,circumference=C,condition="noise",fraction=fraction,window="full",gate=.1)
                axes[i,0].plot([float(r["lag"]) for r in chosen],[float(r["eligible_pair_fraction"]) for r in chosen],
                               color=COLORS[C],ls=style,lw=1.1)
                chosen=subset(site,model=model,circumference=C,condition="noise",fraction=fraction,window="full",metric="shortest_physical_displacement_squared")
                axes[i,1].plot([float(r["lag"]) for r in chosen],[np.sqrt(float(r["value"])) for r in chosen],
                               color=COLORS[C],ls=style,lw=1.1)
                axes[i,2].plot([float(r["lag"]) for r in chosen],[float(r["eligible_pair_fraction"]) for r in chosen],
                               color=COLORS[C],ls=style,lw=1.1)
        axes[i,0].set_ylabel(NAMES[model]+"\nEligible fraction")
        axes[i,1].set_ylabel("Physical RMS displacement")
        axes[i,2].set_ylabel("Eligible fraction")
        axes[i,0].set_ylim(-.02,1.04);axes[i,2].set_ylim(-.02,1.04)
        axes[i,1].set_ylim(-.03,1.85)
        for ax in axes[i]:
            ax.set_xlim(0,1000)
            if i==2:ax.set_xlabel("Lag (native time)")
    for ax,title in zip(axes[0],("Direction-pair coverage (R ≥ 0.1)","Single-region endpoint displacement","Single-region pair coverage")):
        ax.set_title(title,fontsize=9)
    axes[1,1].text(.05,.72,"High-input G: no eligible pairs\nat lag 1000; low input coverage<6%",transform=axes[1,1].transAxes,fontsize=7.5)
    handles=[Line2D([],[],color=COLORS[C],label=f"C = {C:g}") for C in (20.,60.)]
    handles += [Line2D([],[],color=".2",ls=style,label=f"f = {f:g}") for f,style in ((.001,"--"),(.01,"-"))]
    fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,bbox_to_anchor=(.5,-.01))
    fig.suptitle("Coverage and physical geometry: endpoint comparisons do not establish site ancestry",y=1.005)
    fig.tight_layout(rect=(0,.04,1,1));save(fig,"coverage_geometry",artifacts)


def input_audits(manifest):
    result=[]
    for job in manifest["jobs"]:
        if job["kind"]!="noise":continue
        directory=RUN/"workers"/f"{job['slot']:04d}"
        audit=json.loads((directory/"noise_audit.json").read_text())
        result.append(dict(model=job["model"],circumference=job["circumference"],fraction=job["mean_fraction"],seed=job["seed"],
                           local_mean_over_declared=audit["empirical_input_mean"]/job["mean"],
                           local_second_moment_over_mean_squared=audit["empirical_input_second_moment"]/job["mean"]**2,
                           gaussian_spatial_covariance_ell=audit["empirical_spatial_covariance_at_ell"],
                           gaussian_temporal_covariance_tau=audit["empirical_temporal_covariance_at_tau"],
                           integrated_input_variance=audit["empirical_integrated_input_variance"],
                           expected_integrated_input_variance=audit["theoretical_integrated_input_variance"],
                           first_input_mode_second_moment=audit["empirical_first_angular_input_mode_second_moment"],
                           expected_first_input_mode_second_moment=audit["theoretical_first_angular_input_mode_second_moment"],
                           source=str(directory/"noise_audit.json"),source_sha256=sha(directory/"noise_audit.json")))
    path=DATA/(PREFIX+"input_audits.csv")
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=result[0]);writer.writeheader();writer.writerows(result)
    return result,path


def main():
    setup();artifacts=[]
    traces=rows("ensemble_traces");correlations=rows("lag_correlations");directions=rows("direction_correlations")
    summaries=rows("seed_summary");persistence=rows("persistence")
    audit=json.loads((META/(PREFIX+"analysis_audit.json")).read_text())
    manifest=json.loads((RUN/"manifest.json").read_text())
    numerical_path=ROOT/"data/runs/vision_noise_numerical_audit_20261017/NUMERICAL_COMPARISON.json"
    numerical=json.loads(numerical_path.read_text())
    assert numerical["complete"]
    assert len(summaries)==30 and len({r["case_id"] for r in summaries})==30
    assert all(int(r["n_seeds"])==1 for r in persistence)
    if "--report-only" in sys.argv:
        for stem in ("initial_memory","magnitude","correlations","coverage_geometry","all_noise_histories"):
            paths=[OUT/(PREFIX+stem+".png")]+[OUT/"supporting/vectors"/(PREFIX+stem+"."+ext) for ext in ("pdf","svg")]
            assert all(path.is_file() for path in paths)
            artifacts.extend(paths)
    else:
        dynamics(traces,"initial_memory","Initial-axis memory  Re P",(-.5,1.2),"initial_memory",artifacts)
        dynamics(traces,"R","Relative magnitude  |P|",(-.02,1.2),"magnitude",artifacts)
        correlation_plot(correlations,directions,artifacts)
        geometry_diagnostic(directions,rows("site_displacement"),artifacts)
        kymographs(summaries,artifacts)
    input_rows,input_path=input_audits(manifest);artifacts.append(input_path)
    status=audit["manifest_audit"]["statuses"]
    table=[]
    for model in MODELS:
        for f in (.001,.01):
            for C in (20.,60.):
                row,=subset(summaries,model=model,circumference=C,condition="noise",fraction=f)
                mean,=subset(summaries,model=model,circumference=C,condition="mean",fraction=f)
                lag,=subset(persistence,model=model,circumference=C,condition="noise",fraction=f,window="full")
                table.append(f"| {NAMES[model]} | {C:g} | {f:g} | {float(row['final_memory']):.4f} | {float(row['minimum_R']):.4f} | {float(row['final_R']):.4f} | {float(mean['final_R']):.4g} | {float(lag['correlation_at_maxlag']):.4f} |")
    zero=subset(summaries,condition="zero")
    zero_drift=max(abs(float(r["final_memory"])-1) for r in zero)
    text=f'''# Noisy-ring pilot: actual response regimes and numerical feasibility

**All 30 prescribed pilot cases completed.** This is one noisy realization at
each model/size/strength, plus its deterministic mean and zero-input controls;
it is not a replicated persistence comparison. Native time/input units are
uncalibrated across models. Main-ensemble launch remains a separate decision.
Generated {datetime.now(timezone.utc).isoformat()}.

## Message and decision supported

The frozen input produces a useful range without changing its parameters to
obtain a preferred size trend. Wave-pinning changes weakly over this horizon;
Holmes exhibits axis changes while retaining a strong resultant; Goryachev at
C20 and the higher mean reaches weak/multiregion states. Crucially, **the
Goryachev C20 deterministic high-mean control itself loses polarity**. That
condition belongs to a different forced-state regime; its noisy response is
not a clean comparison of directional memory around the original polarized
branch. C60 retains a substantial resultant under the same local mean.

These are actual one-seed regimes for assessing the predeclared ensemble, not
proof of a size trend. Both strengths and every control remain visible. The
local noise mean and physical correlation length/time are held fixed within
model; total input and its whole-ring angular modes vary with circumference.

## Input and readouts

The native positive activation input is

$$
I(s,t)=fA_m\\exp[Z(s,t)-1/2],\\qquad
f\\in\\{{0.001,0.01\\}},\\qquad A_m\\in\\{{0.1,6.4,1\\}}.\\tag{{1}}
$$

The three $A_m$ values belong to Wave-pinning, Goryachev and Holmes respectively.
They are native cue references, not equal biological doses. The Gaussian field
has unit ensemble pointwise variance, periodized Gaussian spatial covariance
with physical Gaussian width 1.5 and exact OU knots with correlation time 50.
Knots are separated by 5; linear interpolation is divided by its analytic
standard deviation. This is a declared interpolated stochastic law, not an
exact continuous OU path. The ensemble input mean is $fA_m$; no realization is
renormalized, demeaned or clipped. Uniform mean-input controls can contribute
to the active first moment through multiplication by the nonuniform inactive
field; a uniform input harmonic alone does not imply zero stimulus moment.

The primary state vector is $P(t)=M_1(t)/M_1(0)$, initially $(1,0)$.
Initial-axis memory $\\operatorname{{Re}}P$ is separate from the uncentered
lag correlation. Magnitude, total activity, both harmonic ratios and native
profiles accompany it. Correlations use every 5-unit lag through 1000 and are
computed separately over full/first/last-half windows without selecting burn-in.
The signed restricted integral retains negative correlations and has no forced
exponential fit. A one-seed pilot has **no replication confidence interval**.
The [frozen statistical definitions](../vision_noise_20261017/DESIGN.md) give
all equations, segmentation rules and interpretation limits.

## Initial memory and magnitude, including controls

![Every pilot memory trace]({PREFIX}initial_memory.png)

Solid lines are individual noisy histories, dashed lines their deterministic
means and dotted lines zero input. Blue/orange denote C20/C60. The common axis
keeps weak changes visibly small. Initial memory can change because magnitude
changes, direction changes, or both; it is not automatically a turning rate.

![Every pilot magnitude trace]({PREFIX}magnitude.png)

The lower high-input Goryachev C20 mean curve establishes the distinct
forced-state regime. A small first moment alone does not identify homogeneous
activity: the local profiles and second harmonic remain necessary.

## Lag memory and direction are kept distinct

![Vector and direction correlations]({PREFIX}correlations.png)

Solid curves are full-window vector correlations; dashed curves are
magnitude-gated direction cosines. Direction results require both endpoint
magnitudes to exceed 0.1. Their usable-pair fractions and 0.05/0.2 sensitivities
are retained in the [direction table](supporting/data/{PREFIX}direction_correlations.csv).
The deterministic mean is dotted where endpoint RMS magnitudes exceed 0.05.
Its weak Goryachev C20 values are withheld from this display, explicitly flagged
and still retained in the raw table. Normalizing a nearly vanished vector must
not produce a claim of persistent heading.

Most conditions have not resolved strong decorrelation by lag 1000. The
restricted integral is an observation-window descriptor, not an extrapolated
persistence time. The [full/half-window table](supporting/data/{PREFIX}persistence.csv)
exposes adaptation and weak endpoint energies; it should be read beside these
raw curves rather than used as a replacement for them.

| Model | C | f | Noisy final Re P | Noisy min R | Noisy final R | Mean final R | Noisy C(1000) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
'''+"\n".join(table)+f'''

## Pair coverage and physical versus angular motion

![Coverage and physical displacement]({PREFIX}coverage_geometry.png)

The single-region endpoint diagnostic additionally requires one current
half-contrast region, at most one peak with native prominence at least 0.05,
contrast at least 0.05 and R at least 0.1 at both endpoints. Missing eligible
pairs remain unavailable. In particular, high-input Goryachev has zero eligible
single-region pairs at lag 1000, and low-input coverage is below 6% there; this
pilot supplies no resolved Goryachev physical-drift comparison under that
criterion. The native histories remain the appropriate evidence about its
local pattern dynamics.

The high-input Wave-pinning example illustrates the geometric distinction:
C20/C60 direction cosines at lag 1000 are 0.9699/0.9949, whereas eligible physical
RMS displacements are 0.7884/0.9853 native length units (pair coverage
98.7%/90.7%). Here the larger ring has stronger angular memory while its
observed endpoint physical displacement is larger. This one-seed comparison
is not a population trend; it shows why angular memory must not be relabelled
slower physical site motion. All panels keep the same physical-displacement
axis; the small Goryachev values have very low eligibility, rather than proving
exceptionally stable same-site motion.

## Every noisy native history

![All twelve noisy histories]({PREFIX}all_noise_histories.png)

All twelve stochastic histories are shown, so no pilot example is selected to
support enlargement increasing persistence. Each row uses one native activity
colour scale across both sizes and strengths; the periodic cut keeps the
original incumbent at position zero. The same angular axis spans different
physical distances at C20 and C60. Quantitative
[endpoint-conditioned physical displacement](supporting/data/{PREFIX}site_displacement.csv)
therefore accompanies angular correlations. These 20-unit profile-cadence endpoints do
not prove continued site identity. A dominant-site change alone is not called
replacement. The [profile table](supporting/data/{PREFIX}profile_diagnostics.csv)
retains current/frozen half-contrast region counts, prominent peaks and tracking
ambiguities.

The lowest noisy Goryachev C20 high-input resultant occurs at $t=3770$:
$R=0.00721$, while contrast is 3.281, $q_2=0.1970$ and two current half-contrast
regions with three prominent peaks remain. Its deterministic mean control ends
with contrast approximately 2.95×10⁻¹² and R approximately 3.01×10⁻¹³. Thus the
noisy low first moment is not equivalent to homogeneous activity. These exact
scalar events are retained in the corresponding worker trace archives.

The medoid candidate table is [retained](supporting/data/{PREFIX}representative_histories.csv)
for the later replicated figure. With one seed per condition its within-group
medoids are simply those individual histories, not representative population
examples. A separate lowest-resultant history is labelled a diagnostic.

## Actual input statistics and numerical accounting

The [twelve input audits](supporting/data/{PREFIX}input_audits.csv) preserve
realized local means/second moments, spatial and temporal covariance checks,
and the zeroth/first whole-ring input modes alongside their finite-ring
expectations. Finite-realization deviations are reported, not corrected by
rescaling. Whole-ring mode amplitudes must not be mistaken for equal forcing
at different circumferences merely because local covariance is fixed.

The independent read-only analysis verified {audit['manifest_audit']['frozen_source_count']}
frozen executable source hashes, {audit['manifest_audit']['external_count']} original
input/provenance hashes and {audit['manifest_audit']['worker_output_hashes_verified']}
worker-output hashes. All 30 initial primary profiles equal their exact copied
basals. All scalar and profile grids reach 4000 at their prescribed 5/20-unit
cadences. Every retained profile independently reproduces the raw zeroth/first/
second moments and native region/peak counts. Maximum raw-moment discrepancy
is {max(r['moment_error'] for r in audit['sources']):.3g}; maximum conserved-pool
relative error is {max(r['max_pool_relative_error'] for r in status):.3g}.
The maximum zero-control final-memory deviation is {zero_drift:.3g}.
Worker elapsed time ranges from {min(r['elapsed_seconds'] for r in status):.2f}
to {max(r['elapsed_seconds'] for r in status):.2f} seconds.

The [analysis audit](supporting/metadata/{PREFIX}analysis_audit.json) records
complete exact sources and outputs. The analysis also passes direct versus FFT
lag-product checks, an analytic rotating-vector correlation, periodic harmonic
and peak tests, and magnitude-gate pair-count tests. The [complete numerical comparison](../../../../data/runs/vision_noise_numerical_audit_20261017/NUMERICAL_COMPARISON.json)
contains six independently verified C60 high-input checks. Tightening ODE
tolerances for the identical forcing changes $P$ by at most
$4.53\\times10^{{-8}}$, $1.40\\times10^{{-9}}$ and $7.70\\times10^{{-8}}$
in Wave-pinning, Goryachev and Holmes respectively. The coupled 5-to-2.5-unit
interpolation comparison changes maximum $P$ by 0.00546, 0.00176 and 0.01807;
maximum raw-profile changes are 2.31%, 0.219% and 12.74% of basal contrast, and
maximum dominant-region position differences are 0.0566, 0.0156 and 0.1866
native length units. Native profiles are therefore not numerically identical.
This comparison changes the declared interpolated driving path and measures
sensitivity to that construction; it is not additional ODE error or proof of
continuous-OU convergence.

Independent direct recomputation of full-window vector correlations finds
maximum coarse/fine curve differences 0.000110, 0.0000176 and 0.000563. Their
restricted integrals change 997.469→997.534, 998.324→998.323 and
992.855→992.469 native time units. Thus these global pilot summaries are less
sensitive than individual local profiles over this check. Prominent-peak counts
can change under the interpolation refinement, especially in Goryachev; this
reinforces the explicit coverage limits on its local-site displacement summary.
A different uncoupled random realization would not be a convergence check.

## Why these displays and what could be measured

Memory and magnitude traces distinguish maintained directional asymmetry from
loss of its resultant; mean controls expose a change of the forced state.
Vector/direction curves prevent amplitude fluctuations being relabelled as
turning. The native histories explain whether vector changes accompany a moving
single site, multiple regions or loss of local contrast. The same ingredients
could be obtained from a whole-boundary activity reporter, with physical site
positions retained alongside the angular vector. This proposes an observable
polarity-memory measurement; it supplies no cell velocity law or calibrated
migration persistence estimate.

## Reviewed pilot and bounded main continuation

The [fresh actual-result critique](supporting/metadata/FRESH_RESULT_CRITIQUE.md)
finds no remaining essential blocker after inspecting all five actual plots,
controls, raw archives and complete numerical comparisons. The parent has
reviewed this evidence and authorized the bounded main ensemble. No extra
amplitude has been selected from these outcomes, and no missing or failed
realization is silently censored.

The main design retains both declared strengths and explicitly labels the
high-mean Goryachev C20 regime. It uses 12 independent seeds per
model/size/strength at C20/40/60 over 32000 native units, with a shared maximum
lag of 8000. The longer prospective horizon addresses unresolved weak-pilot
decay without selecting a favourable size trend; it is not yet an ensemble
result. Main site diagnostics use the scalar cadence of 5, independently
checked at every saved profile snapshot; main saved profiles use cadence 40.
Pilot physical-displacement diagnostics remain at their declared cadence of 20.
This report records the completed pilot, not the outcome of that continuation.
'''
    path=OUT/"RESULTS.md";path.write_text(text);artifacts.append(path)
    manifest_path=META/(PREFIX+"report_manifest.json")
    sources={str(p.relative_to(REPO)):sha(p) for p in [Path(__file__),RUN/"manifest.json",META/(PREFIX+"analysis_audit.json"),numerical_path]}
    for file in DATA.glob(PREFIX+"*.csv"):
        sources[str(file.relative_to(REPO))]=sha(file)
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root="../..",source_sha256=sources,
                sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},
                no_integration=True,cases=30,replicated_inference=False,
                selection="All twelve noisy pilot histories; all mean and zero controls retained in scalar figures")
    manifest_path.write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(dict(report=str(path),artifacts=len(artifacts)),indent=2))


if __name__=="__main__":main()
