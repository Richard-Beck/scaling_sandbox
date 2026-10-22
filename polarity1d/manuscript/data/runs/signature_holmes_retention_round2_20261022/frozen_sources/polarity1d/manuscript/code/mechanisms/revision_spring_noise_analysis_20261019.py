"""Independent read-only Spring ring noise audits, statistics and render inputs."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
from scipy.signal import find_peaks
from vision_noise_analysis_stable_v2_20261017 import correlation_parts, direction_parts

ROOT=Path(__file__).resolve().parents[4]
FIG=ROOT/"polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4"
META=FIG/"supporting/metadata"
DATA=FIG/"supporting/data"
DRAW_COUNT=2000
BOOTSTRAP_SEED=2026101905


def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,x):Path(path).write_text(json.dumps(x,indent=2,allow_nan=True)+"\n")
def rows(path):return list(csv.DictReader(Path(path).open()))


def csvwrite(name,values):
    p=DATA/(name+".csv")
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for v in values for k in v)))
        w.writeheader();w.writerows(values)
    return p


def ratio(x,y):
    x,y=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float))
    return np.divide(x,y,out=np.full(x.shape,np.nan),where=y>0)


def load(run,job):
    folder=run/"workers"/f"{job['slot']:04d}"
    s=read(folder/"status.json")
    assert s["status"]=="complete" and s["producer"]=="terminal"
    assert s["manifest_sha256"]==sha(run/"manifest.json")
    for name,h in s["output_sha256"].items():assert sha(folder/name)==h
    with np.load(folder/"trace.npz") as z:r={k:z[k].copy() for k in z.files}
    r["P"]=r["P_real"]+1j*r["P_imag"]
    r["R"]=abs(r["P"])
    r["strict"]=(r["current_half_contrast_regions"]==1)&(r["peaks_prominence_ge_005"]<=1)&(r["contrast"]>=.05)&(r["R"]>=.1)
    r["one_region"]=(r["current_half_contrast_regions"]==1)&(r["contrast"]>=.05)&(r["R"]>=.1)
    r["job"],r["status"],r["folder"]=job,s,folder
    return r


def direct_statistics(r,maxlag):
    P,R,q,C=r["P"],r["R"],r["dominant_region_fraction"],r["circumference"]
    raw,direction,coverage,site=[],[],[],[]
    for k in range(maxlag+1):
        a,b=(P,P) if k==0 else (P[k:],P[:-k])
        raw.append(float(np.real(np.vdot(b,a))/np.sqrt(np.sum(abs(a)**2)*np.sum(abs(b)**2))))
        gate=(abs(a)>=.1)&(abs(b)>=.1)
        direction.append(float(np.mean(np.cos(np.angle(a[gate]*b[gate].conj())))) if gate.any() else np.nan)
        v=r["strict"] if k==0 else r["strict"][k:]&r["strict"][:-k]
        dq=np.zeros(len(q)) if k==0 else (q[k:]-q[:-k]+.5)%1-.5
        coverage.append(float(np.mean(v)))
        site.append(float(np.sqrt(np.mean((C[0]*dq[v])**2))) if v.any() else np.nan)
    return dict(raw=np.array(raw),direction=np.array(direction),coverage=np.array(coverage),site_RMS=np.array(site))


def refinement_audit(run,original):
    m=read(run/"manifest.json")
    plan=read(META/"spring_noise_refinement_plan.json")
    assert len([j for j in m["jobs"] if j["variant"]=="fine_knots"])==4
    records=[]
    for j in [j for j in m["jobs"] if j["variant"]!="standard"]:
        base,=[b for b in m["jobs"] if b["variant"]=="standard" and b["kind"]=="noise" and b["seed"]==j["seed"]]
        a,b=load(run,base),load(run,j)
        aa=read(a["folder"]/"noise_audit.json");ab=read(b["folder"]/"noise_audit.json")
        assert aa["latent_fine_knots_sha256"]==ab["latent_fine_knots_sha256"]
        old=[v for v in read(original/"manifest.json")["jobs"] if v["kind"]=="noise" and v["seed"]==j["seed"] and v["variant"]=="standard"]
        if old:
            previous=read(original/"workers"/f"{old[0]['slot']:04d}"/"noise_audit.json")
            assert aa["preserved_base_knots_sha256"]==previous["latent_fine_knots_sha256"]
        ta,tb=direct_statistics(a,200),direct_statistics(b,200)
        with np.load(a["folder"]/"primary.npz") as z:ua=z["primary"].copy()
        with np.load(b["folder"]/"primary.npz") as z:ub=z["primary"].copy()
        site_gate=a["strict"]&b["strict"]
        dq=(a["dominant_region_fraction"]-b["dominant_region_fraction"]+.5)%1-.5
        diffs={name:float(np.nanmax(abs(ta[name]-tb[name]))) for name in ta}
        morph=dict(weak_R_fraction_difference=float(abs(np.mean(a["R"]<.1)-np.mean(b["R"]<.1))),
                   multiregion_fraction_difference=float(abs(np.mean(a["current_half_contrast_regions"]>1)-np.mean(b["current_half_contrast_regions"]>1))))
        physical_tolerance=np.maximum(.01,.1*ta["site_RMS"])
        site_defined=np.isfinite(ta["site_RMS"])&np.isfinite(tb["site_RMS"])
        gates=dict(raw_correlation=diffs["raw"]<=.05,direction_correlation=diffs["direction"]<=.05,
            weak_R=morph["weak_R_fraction_difference"]<=.02,multiregion=morph["multiregion_fraction_difference"]<=.02,
            site_coverage=diffs["coverage"]<=.05,
            physical_RMS=bool(np.all(abs(ta["site_RMS"][site_defined]-tb["site_RMS"][site_defined])<=physical_tolerance[site_defined])) if site_defined.any() else None)
        records.append(dict(seed=j["seed"],fraction=j["fraction"],variant=j["variant"],
            baseline_slot=base["slot"],variant_slot=j["slot"],preserved_original_quarter_unit_knots=bool(old),
            identical_latent_fine_knots=True,statistic_max_differences=diffs,morphology=morph,gates=gates,
            every_available_gate_passes=all(v for v in gates.values() if v is not None),
            max_pointpath_P_difference=float(np.max(abs(a["P"]-b["P"]))),
            max_primary_difference_over_initial_contrast=float(np.max(abs(ua-ub))/np.ptp(ua[0])),
            max_primary_relative_L1=float(np.max(np.mean(abs(ua-ub),axis=1)/np.mean(ua,axis=1))),
            max_circumference_difference=float(np.max(abs(a["circumference"]-b["circumference"]))),
            max_gated_site_angle_difference=float(np.max(abs(2*np.pi*dq[site_gate]))) if site_gate.any() else np.nan,
            gated_site_comparison_fraction=float(site_gate.mean()),
            scalar_samples=len(a["time"]),profiles=len(ua)))
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),run=str(run),analysis_source=str(Path(__file__)),
        source_sha256=sha(__file__),plan_sha256=sha(META/"spring_noise_refinement_plan.json"),
        original_pilot=str(original),original_source_unchanged=True,paired_cases=records,
        all_estimator_gates_pass=all(r["every_available_gate_passes"] for r in records),
        actual_worker_elapsed_seconds=[r["elapsed_seconds"] for r in read(run/"results/status.json")["statuses"]],
        scope="Finite-knot forcing-law sensitivity, distinct from same-law ODE error; no continuous-OU or long-horizon convergence proof")
    write(META/"spring_noise_bridge_audit.json",result)
    print(json.dumps(result,indent=2))
    return result


def bounds(x):
    x=np.asarray(x)
    out=np.full((2,)+x.shape[1:],np.nan)
    valid=np.isfinite(x).any(axis=0)
    if valid.any():out[:,valid]=np.nanpercentile(x[:,valid],(2.5,97.5),axis=0)
    return out


def summarize(numerator,denominator,w):
    n=np.asarray(numerator,float);d=np.asarray(denominator,float)
    value=ratio(n.sum(axis=0),d.sum(axis=0))
    draws=ratio(w@n,w@d) if len(w) else np.empty((0,)+value.shape)
    lo,hi=bounds(draws) if len(w) else (np.full(value.shape,np.nan),np.full(value.shape,np.nan))
    return value,lo,hi,np.mean(np.isfinite(draws),axis=0) if len(w) else np.full(value.shape,np.nan)


def validate_profiles(r):
    with np.load(r["folder"]/"primary.npz") as z:
        t,u,C,q=z["time"],z["primary"],z["circumference"],z["material_positions"]
    ix=np.rint(t/.5).astype(int)
    assert np.array_equal(r["time"][ix],t)
    M1=C*(u@np.exp(2j*np.pi*q))/len(q)
    old=r["M1_real"][ix]+1j*r["M1_imag"][ix]
    e1=float(np.max(abs(M1-old)))
    e0=float(np.max(abs(C*np.mean(u,axis=1)-r["M0"][ix])))
    assert e1<1e-11 and e0<1e-11
    contrast=np.ptp(u,axis=1)
    mask=u>=u.min(axis=1)[:,None]+.5*contrast[:,None]
    count=np.sum(mask&~np.roll(mask,1,axis=1),axis=1)
    count[mask.all(axis=1)]=1;count[contrast<=1e-12]=0
    np.testing.assert_array_equal(count,r["current_half_contrast_regions"][ix])
    pe=[]
    for field in u:
        pk,_=find_peaks(np.tile(field,3),prominence=.05)
        pe.append(np.sum((pk>=len(q))&(pk<2*len(q))))
    np.testing.assert_array_equal(pe,r["peaks_prominence_ge_005"][ix])
    return dict(slot=r["job"]["slot"],profiles=len(t),max_M1_error=e1,max_M0_error=e0,regions_match=True,peaks_match=True)


def analyze(run):
    m=read(run/"manifest.json")
    assert m["mode"]=="main" and len(m["jobs"])==81
    for relative,h in m["source_sha256"].items():assert sha(run/"frozen_sources"/relative)==h
    records=[load(run,j) for j in m["jobs"]]
    profile_audit=[validate_profiles(r) for r in records]
    rng=np.random.default_rng(BOOTSTRAP_SEED)
    correlation_rows,direction_rows,trace_rows,site_rows,obs_rows,seed_rows,seed_site_rows,increment_rows,persistence_rows=[],[],[],[],[],[],[],[],[]
    selection=[]
    reference_C=max(b["basal_circumference"] for b in m["basals"])
    lags=np.arange(1601)*.5
    for C0 in (1.,2.,3.):
        for condition,fraction in [("noise",.001),("noise",.01),("mean",.001),("mean",.01),("basal",0.)]:
            rr=[r for r in records if r["job"]["resting_circumference"]==C0 and r["job"]["kind"]==condition and r["job"]["fraction"]==fraction]
            n=len(rr);assert n==(12 if condition=="noise" else 1)
            w=rng.multinomial(n,np.full(n,1/n),DRAW_COUNT)/n if n>1 else np.empty((0,n))
            base=dict(model="spring_ring",resting_circumference=C0,basal_circumference=rr[0]["circumference"][0],
                      condition=condition,fraction=fraction,n_seeds=n,native_time_calibration="none")
            descriptors=[]
            for r in rr:
                seed=dict(**base,seed=r["job"]["seed"],slot=r["job"]["slot"],observations=len(r["time"]),
                    final_memory=float(r["P"][-1].real),minimum_R=float(r["R"].min()),mean_R=float(r["R"].mean()),final_R=float(r["R"][-1]),
                    multiregion_fraction=float(np.mean(r["current_half_contrast_regions"]>1)),weak_R_fraction=float(np.mean(r["R"]<.1)),
                    strict_site_fraction=float(r["strict"].mean()),one_region_fraction=float(r["one_region"].mean()),
                    circumference_minimum=float(r["circumference"].min()),circumference_maximum=float(r["circumference"].max()),
                    maximum_relative_C_excursion=float(np.max(abs(r["circumference"]/r["circumference"][0]-1))),
                    archive=str(r["folder"]/"trace.npz"),archive_sha256=sha(r["folder"]/"trace.npz"))
                seed_rows.append(seed);descriptors.append([seed[k] for k in ("final_memory","minimum_R","mean_R","multiregion_fraction")])
                for lag in (.5,4.,10.,50.,200.,800.):
                    k=int(lag/.5);angle=abs(np.angle(r["P"][k:]*r["P"][:-k].conj()))
                    for gate in (.05,.1,.2):
                        valid=(r["R"][k:]>=gate)&(r["R"][:-k]>=gate)
                        a=angle[valid]
                        increment_rows.append(dict(**base,seed=seed["seed"],slot=seed["slot"],lag=lag,gate=gate,
                            eligible_pair_fraction=float(valid.mean()),pairs_eligible=int(valid.sum()),
                            angle_q50=float(np.quantile(a,.5)) if len(a) else np.nan,
                            angle_q95=float(np.quantile(a,.95)) if len(a) else np.nan,
                            angle_q99=float(np.quantile(a,.99)) if len(a) else np.nan,
                            max_absolute_angle=float(a.max()) if len(a) else np.nan,
                            count_above_pi_over_4=int(np.sum(a>np.pi/4)),classifier="none"))
            if condition=="noise":
                d=np.array(descriptors);iqr=np.quantile(d,.75,axis=0)-np.quantile(d,.25,axis=0);iqr[iqr==0]=1
                distances=np.sum(abs(d-np.median(d,axis=0))/iqr,axis=1)
                chosen=min(range(n),key=lambda i:(distances[i],rr[i]["job"]["slot"]))
                selection.append(dict(**base,slot=rr[chosen]["job"]["slot"],seed=rr[chosen]["job"]["seed"],
                    profile_archive=str(rr[chosen]["folder"]/"primary.npz"),medoid_distance=float(distances[chosen]),
                    role="within_condition_descriptor_medoid"))
            groupseeds=seed_rows[-n:]
            for metric in ("mean_R","final_R","multiregion_fraction","weak_R_fraction","strict_site_fraction","one_region_fraction","maximum_relative_C_excursion"):
                vals=np.array([s[metric] for s in groupseeds]);mean=float(vals.mean())
                lo,hi=np.percentile(w@vals,(2.5,97.5)) if len(w) else (np.nan,np.nan)
                obs_rows.append(dict(**base,metric=metric,value=mean,ci_lower=lo,ci_upper=hi))
            for metric in ("R","initial_memory","circumference"):
                data=np.stack([r["P"].real if metric=="initial_memory" else r[metric] for r in rr])[:,::8]
                lo,hi=np.percentile(w@data,(2.5,97.5),axis=0) if len(w) else (np.full(data.shape[1],np.nan),np.full(data.shape[1],np.nan))
                for k,t in enumerate(rr[0]["time"][::8]):trace_rows.append(dict(**base,time=t,metric=metric,value=float(data[:,k].mean()),ci_lower=lo[k],ci_upper=hi[k]))
            for window,start,stop in (("full",0.,3200.),("first_half",0.,1600.),("last_half",1600.,3200.)):
                parts=[]
                for r in rr:
                    take=(r["time"]>=start)&(r["time"]<=stop)
                    parts.append(correlation_parts(r["P"][take],1600))
                nu=np.stack([p[0] for p in parts]);eo=np.stack([p[1] for p in parts]);el=np.stack([p[2] for p in parts])
                value=ratio(nu.mean(axis=0),np.sqrt(eo.mean(axis=0)*el.mean(axis=0)))
                draws=ratio(w@nu,np.sqrt((w@eo)*(w@el))) if len(w) else np.empty((0,1601))
                lo,hi=bounds(draws) if len(w) else (np.full(1601,np.nan),np.full(1601,np.nan))
                restricted=float(np.trapezoid(value,lags))
                area_draws=np.trapezoid(draws,lags,axis=1) if len(w) else np.array([])
                area_lo,area_hi=np.percentile(area_draws,(2.5,97.5)) if len(w) else (np.nan,np.nan)
                below=np.flatnonzero(value<=1/np.e)
                persistence_rows.append(dict(**base,window=window,restricted_integral=restricted,ci_lower=area_lo,ci_upper=area_hi,
                    restricted_horizon=800.,first_sampled_1_over_e_crossing=float(lags[below[0]]) if len(below) else np.nan,
                    minimum_endpoint_rms=float(np.sqrt(min(eo.mean(axis=0).min(),el.mean(axis=0).min()))),
                    exponential_fit=False,burn_in_selected=False,negative_correlations_retained=True))
                for k,lag in enumerate(lags):correlation_rows.append(dict(**base,window=window,lag=lag,correlation=value[k],ci_lower=lo[k],ci_upper=hi[k],
                    endpoint_rms_origin=float(np.sqrt(eo[:,k].mean())),endpoint_rms_later=float(np.sqrt(el[:,k].mean())),
                    status="weak_endpoint_rms" if min(np.sqrt(eo[:,k].mean()),np.sqrt(el[:,k].mean()))<.05 else ("defined" if np.isfinite(value[k]) else "undefined")))
                for gate in (.05,.1,.2):
                    parts=[direction_parts(r["P"][(r["time"]>=start)&(r["time"]<=stop)],1600,gate) for r in rr]
                    nu=np.stack([p[0] for p in parts]);den=np.stack([p[1] for p in parts])
                    value,lo,hi,defined=summarize(nu,den,w)
                    for k,lag in enumerate(lags):direction_rows.append(dict(**base,window=window,gate=gate,lag=lag,
                        direction_correlation=value[k],ci_lower=lo[k],ci_upper=hi[k],bootstrap_defined_fraction=defined[k],
                        eligible_pair_fraction=float(den[:,k].sum()/sum(p[2][k] for p in parts))))
            # Actual displacement distributions, kept distinct from global-vector angles.
            for gate_name in ("strict","one_region"):
                ns={k:[] for k in ("basal_arc_squared","midpoint_arc_squared","native_cosine","reference_cosine","delta_C_squared")};dens=[]
                for r in rr:
                    C,q=r["circumference"],r["dominant_region_fraction"]
                    sums={k:[] for k in ns};counts=[]
                    for k,lag in enumerate(lags):
                        v=r[gate_name] if k==0 else r[gate_name][k:]&r[gate_name][:-k]
                        dq=np.zeros(len(q)) if k==0 else (q[k:]-q[:-k]+.5)%1-.5
                        midpoint=C if k==0 else .5*(C[k:]+C[:-k])
                        dC=np.zeros(len(C)) if k==0 else C[k:]-C[:-k]
                        ds=C[0]*dq
                        vals=(ds**2,(midpoint*dq)**2,np.cos(2*np.pi*dq),np.cos(2*np.pi*ds/reference_C),dC**2)
                        counts.append(int(v.sum()))
                        for name,val in zip(sums,vals):sums[name].append(float(val[v].sum()))
                        if lag in (.5,4.,10.,50.,200.,800.):
                            seed_site_rows.append(dict(**base,seed=r["job"]["seed"],slot=r["job"]["slot"],gate=gate_name,lag=lag,
                                pairs_eligible=int(v.sum()),eligible_pair_fraction=float(v.mean()),reference_C=reference_C,
                                **{key+"_sum":float(val[v].sum()) for key,val in zip(sums,vals)}))
                    dens.append(counts)
                    for name in ns:ns[name].append(sums[name])
                den=np.array(dens)
                for metric,nu in ns.items():
                    value,lo,hi,defined=summarize(nu,den,w)
                    for k,lag in enumerate(lags):site_rows.append(dict(**base,gate=gate_name,lag=lag,metric=metric,value=value[k],ci_lower=lo[k],ci_upper=hi[k],
                        eligible_pair_fraction=float(den[:,k].sum()/(n*(6401-k))),independent_seeds_with_pairs=int(np.sum(den[:,k]>0)),
                        eligible_pairs=int(den[:,k].sum()),bootstrap_defined_fraction=defined[k],reference_C=reference_C,same_site_ancestry_assumed=False))
    paths=[csvwrite(name,values) for name,values in (("spring_noise_correlations",correlation_rows),("spring_noise_direction_correlations",direction_rows),
        ("spring_noise_ensemble_traces",trace_rows),("spring_noise_site_statistics",site_rows),("spring_noise_observation_statistics",obs_rows),
        ("spring_noise_seed_summary",seed_rows),("spring_noise_seed_site_statistics",seed_site_rows),
        ("spring_noise_seed_increments",increment_rows),("spring_noise_representative_histories",selection),("spring_noise_persistence",persistence_rows))]
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),run=str(run),manifest_sha256=sha(run/"manifest.json"),
        cases=len(records),noisy_cases=sum(r["job"]["kind"]=="noise" for r in records),all_hashes_verified=True,
        profile_audit=profile_audit,total_profiles=sum(a["profiles"] for a in profile_audit),
        maximum_profile_M1_error=max(a["max_M1_error"] for a in profile_audit),
        estimator_source="vision_noise_analysis_stable_v2_20261017.py",estimator_sha256=sha(Path(__file__).with_name("vision_noise_analysis_stable_v2_20261017.py")),
        new_metrics_source=sha(__file__),bootstrap_draws=DRAW_COUNT,bootstrap_seed=BOOTSTRAP_SEED,unit="whole independent seed",
        scalar_cadence=.5,profile_cadence=4.,horizon=3200.,maximum_lag=800.,
        primary_physical_quantity="Basal-reference arc shift Cbasal*wrapped_delta_q; not motion of material or cell",reference_C=reference_C,
        outputs=[dict(path=str(p),sha256=sha(p)) for p in paths])
    write(META/"spring_noise_analysis_audit.json",audit)
    print(json.dumps({k:v for k,v in audit.items() if k not in ("profile_audit","outputs")},indent=2))
    return audit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("command",choices=("refinement","analyze"));p.add_argument("--run",type=Path,required=True)
    p.add_argument("--original",type=Path,default=ROOT/"polarity1d/manuscript/data/runs/revision_spring_noise_pilot_20261019")
    a=p.parse_args();META.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    if a.command=="refinement":refinement_audit(a.run.resolve(),a.original.resolve())
    else:analyze(a.run.resolve())


if __name__=="__main__":main()
