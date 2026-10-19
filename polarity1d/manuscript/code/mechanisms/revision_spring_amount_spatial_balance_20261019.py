"""Export active-only Spring amount and normalized first harmonic; no integration."""
from pathlib import Path
from datetime import datetime,timezone
import csv
import gzip
import hashlib
import json
import os
import time
for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[name]="1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[4]
RUN=ROOT/"polarity1d/manuscript/data/runs/revision_spring_noise_main_20261019"
OUT=ROOT/"polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4"
DATA=OUT/"supporting/data/spring_amount_spatial_balance"
META=OUT/"supporting/metadata"
VECTORS=OUT/"supporting/vectors"
B=2000
SEED=2026101911
METRICS=("M0_ratio","q1_ratio","R","paired_product")


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_csv(path,rows):
    op=gzip.open if str(path).endswith(".gz") else open
    with op(path,"wt",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def load():
    manifest=json.loads((RUN/"manifest.json").read_text());cases=[];audits=[];index=[]
    for job in manifest["jobs"]:
        folder=RUN/"workers"/f"{job['slot']:04d}";status=json.loads((folder/"status.json").read_text())
        assert status["status"]=="complete" and status["manifest_sha256"]==sha(RUN/"manifest.json")
        for name in ("trace.npz","primary.npz"):assert sha(folder/name)==status["output_sha256"][name]
        with np.load(folder/"trace.npz") as z:
            t=z["time"].copy();C=z["circumference"].copy();M0=z["M0"].copy();M1=z["M1_real"]+1j*z["M1_imag"]
            P=z["P_real"]+1j*z["P_imag"];q1=z["q1"].copy();total=z["amount"].copy()
        np.testing.assert_array_equal(t,np.arange(6401)*.5)
        assert np.all(M0>0) and q1[0]>0
        A=M0/M0[0];Q=q1/q1[0];R=abs(P);product=A*Q
        with np.load(folder/"primary.npz") as z:
            pt=z["time"];u=z["primary"];v=z["inactive"];pc=z["circumference"];q=z["material_positions"]
            np.testing.assert_array_equal(pt,np.arange(801)*4.)
            ii=np.searchsorted(t,pt);np.testing.assert_array_equal(C[ii],pc)
            fm0=pc*u.mean(axis=1);fm1=pc*np.mean(u*np.exp(2j*np.pi*q),axis=1);ftotal=pc*(u+v).mean(axis=1)
            audit=dict(slot=job["slot"],saved_profiles=len(pt),scalar_samples=len(t),trace_sha256=sha(folder/"trace.npz"),primary_sha256=sha(folder/"primary.npz"),
                maximum_active_M0_field_error=float(abs(fm0-M0[ii]).max()),maximum_M1_field_error=float(abs(fm1-M1[ii]).max()),
                maximum_q1_field_error=float(abs(abs(fm1)/fm0-q1[ii]).max()),maximum_conserved_total_field_error=float(abs(ftotal-total[ii]).max()),
                maximum_scalar_q1_definition_error=float(abs(abs(M1)/M0-q1).max()),maximum_scalar_P_definition_error=float(abs(M1/M1[0]-P).max()),
                maximum_factorization_error=float(abs(R-product).max()))
        assert max(audit[k] for k in audit if k.startswith("maximum_"))<1e-10
        meta=dict(slot=job["slot"],resting_C=job["resting_circumference"],condition=job["kind"],fraction=job["fraction"],seed=job["seed"],basal_C=float(C[0]))
        case=dict(**meta,time=t,M0_ratio=A,q1_ratio=Q,R=R,paired_product=product)
        path=DATA/"cases"/f"slot_{job['slot']:04d}.npz"
        np.savez_compressed(path,time=t,M0=M0,M1_real=M1.real,M1_imag=M1.imag,q1=q1,M0_ratio=A,q1_ratio=Q,R=R,paired_product=product,
            factorization_residual=R-product,circumference=C,conserved_active_inactive_amount=total)
        index.append(dict(**meta,scalar_cadence=.5,field_audit_cadence=4.,scalar_samples=len(t),horizon=3200.,initial_M0=float(M0[0]),initial_abs_M1=float(abs(M1[0])),
            initial_q1=float(q1[0]),initial_conserved_total=float(total[0]),trace_export=str(path),trace_export_sha256=sha(path),
            native_trace=str(folder/"trace.npz"),native_trace_sha256=audit["trace_sha256"],native_primary=str(folder/"primary.npz"),native_primary_sha256=audit["primary_sha256"]))
        cases.append(case);audits.append(audit)
    write_csv(DATA/"case_index.csv",index)
    return cases,audits,index


def bounds(samples):return np.percentile(samples,(2.5,97.5),axis=0)


def summaries(cases):
    rng=np.random.default_rng(SEED);time_rows=[];case_rows=[];group_rows=[];plot={}
    groups=sorted(set((c["resting_C"],c["condition"],c["fraction"]) for c in cases))
    for C,condition,f in groups:
        group=[c for c in cases if (c["resting_C"],c["condition"],c["fraction"])==(C,condition,f)]
        n=len(group);t=group[0]["time"]
        assert n==(12 if condition=="noise" else 1)
        arrays={k:np.asarray([c[k] for c in group]) for k in METRICS}
        weights=rng.multinomial(n,np.ones(n)/n,size=B)/n if n>1 else None
        base=dict(resting_C=C,basal_C=group[0]["basal_C"],condition=condition,fraction=f,independent_cases=n,scalar_cadence=.5)
        aggregate={k:arrays[k].mean(axis=0) for k in METRICS};interval={}
        for k in METRICS:
            interval[k]=bounds(weights@arrays[k]) if n>1 else (np.full(len(t),np.nan),np.full(len(t),np.nan))
        mean_product_of_means=aggregate["M0_ratio"]*aggregate["q1_ratio"]
        for j,tt in enumerate(t):
            row=dict(**base,time=float(tt))
            for k in METRICS:row.update({k+"_mean":float(aggregate[k][j]),k+"_ci_lower":float(interval[k][0][j]),k+"_ci_upper":float(interval[k][1][j])})
            row.update(product_of_component_means=float(mean_product_of_means[j]),cross_case_component_covariance=float(aggregate["paired_product"][j]-mean_product_of_means[j]),
                ensemble_factorization_residual=float(aggregate["R"][j]-aggregate["paired_product"][j]))
            time_rows.append(row)
        plot[C,condition,f]=dict(time=t,mean=aggregate,interval=interval)
        for window,start,stop in (("full",0.,3200.),("first_half",0.,1600.),("last_half",1600.,3200.)):
            keep=(t>=start)&(t<=stop);means={k:arrays[k][:,keep].mean(axis=1) for k in METRICS}
            for i,c in enumerate(group):
                row=dict(**{k:c[k] for k in ("slot","resting_C","basal_C","condition","fraction","seed")},window=window,start=start,stop=stop,samples=int(keep.sum()))
                for k in METRICS:row.update({k+"_sample_mean":float(means[k][i]),k+"_sample_minimum":float(arrays[k][i,keep].min()),k+"_sample_maximum":float(arrays[k][i,keep].max())})
                row.update(product_of_component_sample_means=float(means["M0_ratio"][i]*means["q1_ratio"][i]),within_case_time_component_covariance=float(means["paired_product"][i]-means["M0_ratio"][i]*means["q1_ratio"][i]))
                case_rows.append(row)
            row=dict(**base,window=window,start=start,stop=stop,samples_per_case=int(keep.sum()))
            samples={k:weights@means[k] for k in METRICS} if n>1 else {}
            for k in METRICS:
                lo,hi=bounds(samples[k]) if n>1 else (np.nan,np.nan)
                row.update({k+"_mean":float(means[k].mean()),k+"_ci_lower":float(lo),k+"_ci_upper":float(hi)})
            prod=means["M0_ratio"].mean()*means["q1_ratio"].mean();cov=means["paired_product"].mean()-prod
            cov_ci=bounds(samples["paired_product"]-samples["M0_ratio"]*samples["q1_ratio"]) if n>1 else (np.nan,np.nan)
            row.update(product_of_component_means=float(prod),joint_case_time_component_covariance=float(cov),component_covariance_ci_lower=float(cov_ci[0]),component_covariance_ci_upper=float(cov_ci[1]),
                ensemble_factorization_residual=float(means["R"].mean()-means["paired_product"].mean()))
            group_rows.append(row)
    write_csv(DATA/"ensemble_scalar_traces.csv.gz",time_rows);write_csv(DATA/"case_window_summaries.csv",case_rows);write_csv(DATA/"condition_window_summaries.csv",group_rows)
    return plot,group_rows


def render(plot):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig,axes=plt.subplots(3,2,figsize=(11.5,9),sharex=True,sharey="row")
    colors=("#397bab","#40916c","#cf7930")
    for i,(metric,label) in enumerate((("M0_ratio","Active amount M₀/M₀(0)"),("q1_ratio","Relative first harmonic q₁/q₁(0)"),("R","Relative first moment R"))):
        for j,f in enumerate((.001,.01)):
            ax=axes[i,j]
            for C,color in zip((1.,2.,3.),colors):
                for condition,style in (("noise","-"),("mean","--"),("basal",":")):
                    g=plot[C,condition,0. if condition=="basal" else f]
                    ax.plot(g["time"],g["mean"][metric],color=color,ls=style,lw=1.1)
                    if condition=="noise":ax.fill_between(g["time"],*g["interval"][metric],color=color,alpha=.11,lw=0)
            ax.set_title(f"f={f:g}",loc="left",fontsize=10);ax.set_ylabel(label);ax.axhline(1,color=".65",lw=.5)
            ax.set_xlim(0,3200)
            if i==2:ax.set_xlabel("Spring native time")
    fig.suptitle("Spring activity amount and spatial first harmonic",fontsize=15,y=.99)
    fig.text(.5,.941,r"$q_1=|M_1|/M_0,\qquad R=(M_0/M_0(0))(q_1/q_1(0))$",ha="center",fontsize=13)
    fig.text(.5,.909,"Exact per-case/time identity; population mean R is the mean paired product, not the product of separate means.",ha="center",fontsize=8.5)
    handles=[Line2D([],[],color=c,label=f"C₀={int(C)}") for C,c in zip((1,2,3),colors)]
    handles += [Line2D([],[],color=".3",ls=style,label=name) for style,name in (("-","Noise"),("--","Mean input"),(":","No added input"))]
    fig.legend(handles=handles,ncol=6,frameon=False,loc="lower center",bbox_to_anchor=(.5,.018),fontsize=8)
    fig.text(.5,.007,"All 81 native cases; scalar cadence 0.5, field checks every 4; shared row scales; pointwise 95% seed intervals for noisy groups only.",ha="center",fontsize=7.8)
    fig.tight_layout(rect=(0,.065,1,.88),h_pad=2.1,w_pad=2.2)
    files=[]
    for ext in ("png","pdf","svg"):
        p=(OUT if ext=="png" else VECTORS)/("figure4_spring_amount_spatial_balance."+ext);fig.savefig(p,dpi=250,bbox_inches="tight");files.append(p)
    plt.close(fig);return files


def docs(audits,index,rows,artifacts,elapsed):
    methods=r"""# Spring active amount and spatial first harmonic

This adds no integration or forcing change. All81 original cases are retained:72 noisy histories (three sizes × two strengths ×12 seeds), six deterministic mean-input controls and three no-added-input controls. The scalar trace stores active-only M₀ separately from the conserved active-plus-inactive `amount`.

$$
M_0(t)=C(t)\,\operatorname{mean}_j u_j(t),\qquad
M_1(t)=C(t)\,\operatorname{mean}_j[u_j(t)e^{2\pi i q_j}],
\qquad q_1(t)=\frac{|M_1(t)|}{M_0(t)}.
\tag{1}
$$

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad
R(t)=|P(t)|=\frac{M_0(t)}{M_0(0)}\frac{q_1(t)}{q_1(0)}.
\tag{2}
$$

Every scalar record is retained at its actual cadence0.5 over0–3200 native units (6401 samples), without an R gate or adaptation exclusion. Every saved active/inactive profile is replayed at its actual cadence4 (801 profiles per case). Physical M₀ and M₁ are reconstructed with instantaneous circumference; no profile interpolation is introduced. These cadence4 checks support the stored scalar readouts, rather than claim an independent full-field reconstruction at every cadence0.5 frame. Scalar q₁/P definitions and the paired factorization are checked independently. The exported conserved total is separately labelled and is never substituted for active M₀.

`case_index.csv` identifies all81 native trace/profile hashes, exact initial physical quantities and exported case NPZ hashes. Each case NPZ contains physical M₀/M₁, q₁, M₀/M₀(0), q₁/q₁(0), R, the paired product/residual, C(t) and separately labelled conserved total. `ensemble_scalar_traces.csv.gz` supplies every scalar time point for all15 conditions. Controls have one deterministic trajectory and no confidence intervals. Every noisy condition uses all12 independent seeds.

Intervals use2000 whole-seed multinomial bootstrap bags (RNG seed2026101911). The same bags are applied jointly to every component and time; time samples are not independent replicates. Full/first-half/last-half windows are closed0–3200,0–1600 and1600–3200, with6401/3201/3201 samples. The midpoint belongs to both half-windows; they are paired sensitivities rather than independent replicates. Their summaries are arithmetic sample means, not stationarity estimates or new correlation fits. Case minima/maxima are explicitly sampled extrema.

Mean R equals the mean of paired products. It generally differs from the product of separate component means. Both exports retain that latter product, the appropriate component covariance and the factorization residual; bootstrap window covariance is recomputed from the same jointly resampled cases. Covariance here is the exact mean-product difference (population-moment convention), without a Bessel correction. Cross-case covariance at a fixed time is distinguished from covariance over the joint case/time population within a window. No additive causal shares are inferred.

q₁ is the normalized first harmonic of a nonnegative activity field, not a proof of profile identity, a site-lineage label or a general coherence measure. R may exceed1. Amount and first-harmonic changes are measured here rather than inferred from mean R alone. Native Spring units/forcing differ from the chemical campaign; their common M₀/P/q₁ definitions permit a methods comparison, not a pooled exposure, timescale or model ranking. The chemical exports use cadence5 and a32000-unit horizon; Spring uses cadence0.5 and3200. The supporting sheet changes no current main layout.
"""
    # Numeric prose is written with spaces without changing formulas or filenames.
    for a,b in (("All81","All 81"),("all81","all 81"),("cases:72","cases: 72"),("×12","× 12"),("cadence0.5","cadence 0.5"),("over0–3200","over 0–3200"),("cadence4","cadence 4"),("all15","all 15"),("all12","all 12"),("use2000","use 2000"),("seed2026101911","seed 2026101911"),("closed0–3200,0–1600","closed 0–3200, 0–1600"),("with6401","with 6401"),("exceed1","exceed 1"),("cadence5","cadence 5"),("a32000","a 32000"),("and3200","and 3200")):
        methods=methods.replace(a,b)
    methods += "\n## Measured full-window means\n\nAll entries use the same 6401 scalar samples per case and all 12 noisy seeds. Brackets are pointwise 95% whole-seed intervals for the full-window sample mean. R and the paired product agree to roundoff.\n\n| f | Resting C₀ | M₀/M₀(0) | q₁/q₁(0) | R |\n|---|---:|---:|---:|---:|\n"
    for row in rows:
        if row["condition"]!="noise" or row["window"]!="full":continue
        fields=[]
        for key in ("M0_ratio","q1_ratio","R"):
            fields.append(f"{row[key+'_mean']:.6f} [{row[key+'_ci_lower']:.6f}, {row[key+'_ci_upper']:.6f}]")
        methods += f"| {row['fraction']:g} | {row['resting_C']:g} | "+" | ".join(fields)+" |\n"
    methods += "\nFor stronger input at C₀3, the measured full-window mean active-amount ratio is 1.00350, while the first-harmonic ratio is about 0.90458 and mean R is 0.90744. The lower mean relative first moment therefore accompanies lower normalized first harmonic rather than reduced mean active amount in this window. This does not assert constant amount at each frame, rigid translation, a stationary distribution or a causal decomposition of the size response. All deterministic controls and half-window values remain in the tables.\n"
    (OUT/"SPRING_AMOUNT_SPATIAL_BALANCE.md").write_text(methods)
    caption="Spring amount/spatial-balance supplement. Rows show active-only M₀/M₀(0), normalized first-harmonic ratio q₁/q₁(0), and R, with shared vertical scales across strengths within each row. All81 original cases, both strengths and all matched controls remain. Curves use scalar cadence0.5; every stored profile independently checks physical moments at cadence4. Shading is pointwise95% whole-seed bootstrap uncertainty across12 noisy seeds, not temporal variability; deterministic controls carry no uncertainty interval. The exact product identity holds within each case/time; mean R is the mean paired product, not the product of separately averaged components. Native axes and declared forcing remain Spring-specific. Neither R nor q₁ proves rigid profile translation."
    for a,b in (("All81","All 81"),("cadence0.5","cadence 0.5"),("cadence4","cadence 4"),("pointwise95","pointwise 95"),("across12","across 12")):caption=caption.replace(a,b)
    (META/"SPRING_AMOUNT_SPATIAL_BALANCE_CAPTION.md").write_text(caption+"\n")
    outputs=[p for p in DATA.rglob('*') if p.is_file()]+artifacts+[OUT/"SPRING_AMOUNT_SPATIAL_BALANCE.md",META/"SPRING_AMOUNT_SPATIAL_BALANCE_CAPTION.md"]
    report=dict(generated_at=datetime.now(timezone.utc).isoformat(),new_native_integrations=0,source_script=str(Path(__file__)),source_script_sha256=sha(__file__),
        native_manifest_sha256=sha(RUN/"manifest.json"),cases=81,noisy_cases=72,scalar_samples=81*6401,profiles_verified=81*801,scalar_cadence=.5,profile_audit_cadence=4.,
        active_M0_field="M0",conserved_total_field="amount",bootstrap_draws=B,bootstrap_seed=SEED,uncertainty_unit="whole independent seed; same joint bags across components/times/windows",
        audits=audits,maximum_errors={k:max(a[k] for a in audits) for k in audits[0] if k.startswith("maximum_")},elapsed_seconds=elapsed,
        outputs=[dict(path=str(p),sha256=sha(p)) for p in outputs],visual_inspection="pending",independent_review="pending")
    (META/"spring_amount_spatial_balance_validation.json").write_text(json.dumps(report,indent=2))


def main():
    started=time.monotonic();(DATA/"cases").mkdir(parents=True,exist_ok=True)
    cases,audits,index=load();plot,rows=summaries(cases);artifacts=render(plot);docs(audits,index,rows,artifacts,time.monotonic()-started)
    print(json.dumps(dict(cases=81,output=str(DATA),elapsed_seconds=time.monotonic()-started)))


if __name__=="__main__":main()
