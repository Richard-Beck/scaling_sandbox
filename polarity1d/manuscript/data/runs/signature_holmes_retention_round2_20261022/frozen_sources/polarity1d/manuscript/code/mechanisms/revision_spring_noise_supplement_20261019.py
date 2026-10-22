"""Additional read-only Spring coverage: centroids, half windows and increments.

Adds diagnostics to the frozen main analysis without replacing its outputs.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import importlib.util
import json
import os
import sys
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
from revision_spring_noise_analysis_20261019 import load,read,sha,write,csvwrite,ratio,bounds,ROOT,FIG,META,DATA


def centroid_replay(r):
    with np.load(r["folder"]/"primary.npz") as z:t,u,q=z["time"],z["primary"],z["material_positions"]
    n=len(q);errors=[];flat=0
    for time,field in zip(t,u):
        if np.ptp(field)<=1e-12:flat+=1;continue
        m=field>=field.min()+.5*np.ptp(field)
        components=[]
        if m.all():components=[np.arange(n)]
        else:
            starts=np.flatnonzero(m&~np.roll(m,1))
            for start in starts:
                run=[]
                for step in range(n):
                    k=(int(start)+step)%n
                    if not m[k]:break
                    run.append(k)
                components.append(np.asarray(run))
        masses=[np.sum(field[ix]-field.min()) for ix in components]
        ix=components[int(np.argmax(masses))]
        numerator=np.dot(field[ix]-field.min(),np.exp(2j*np.pi*q[ix]))
        if abs(numerator)==0:continue
        center=np.angle(numerator)/(2*np.pi)%1
        old=r["dominant_region_fraction"][int(round(time/.5))]
        errors.append(abs((center-old+.5)%1-.5))
    maximum=float(max(errors,default=0.))
    assert maximum<1e-10
    return dict(slot=r["job"]["slot"],saved_profiles=len(t),observed_current_centroids=len(errors),flat_unobserved_profiles=flat,
                maximum_material_fraction_error=maximum,maximum_basal_arc_error=maximum*r["circumference"][0])


def main(run):
    m=read(run/"manifest.json");assert len(m["jobs"])==81
    records=[load(run,j) for j in m["jobs"]]
    centroid=[centroid_replay(r) for r in records]
    rng=np.random.default_rng(2026101906);draws=2000;reference_C=max(b["basal_circumference"] for b in m["basals"])
    cdf_rows,half_rows,input_rows=[],[],[]
    thresholds=np.r_[0.,np.geomspace(1e-6,np.pi,181)]
    # Regenerate declared forcing byte-for-byte from frozen producer inputs.
    producer=run/"frozen_sources/code/mechanisms/revision_spring_noise_20261019.py"
    sys.path.insert(0,str(producer.parent));spec=importlib.util.spec_from_file_location("frozen_spring_noise_producer",producer)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    for r in records:
        j=r["job"];a=read(r["folder"]/"noise_audit.json")
        base=dict(slot=j["slot"],condition=j["kind"],resting_circumference=j["resting_circumference"],fraction=j["fraction"],seed=j["seed"],
                  source_path=str(r["folder"]/"noise_audit.json"),source_sha256=sha(r["folder"]/"noise_audit.json"),
                  basal_circumference=r["circumference"][0],actual_C_minimum=float(r["circumference"].min()),actual_C_maximum=float(r["circumference"].max()),
                  physical_ell_minimum=float(.15*r["circumference"].min()/r["circumference"][0]),physical_ell_maximum=float(.15*r["circumference"].max()/r["circumference"][0]))
        if j["kind"]=="noise":
            noise=module.MaterialNoise(j["cells"],r["circumference"][0],m["horizon"],j["seed"],4*j["fraction"],knot=.125,coupling=.125)
            assert noise.latent_sha==a["latent_fine_knots_sha256"]
            I=noise.mean*np.exp(noise.knots-.5)
            M0=noise.basal_C*I.mean(axis=1)
            M1=noise.basal_C*(I@np.exp(2j*np.pi*(np.arange(noise.n)+.5)/noise.n))/noise.n
            base.update(local_mean_over_mu=a["empirical_input_mean_over_mu"],local_variance_over_mu2=a["empirical_input_variance_over_mu2"],
                expected_local_variance_over_mu2=a["expected_input_variance_over_mu2"],
                basal_input_integral_mean=float(M0.mean()),basal_input_integral_variance_over_mu2=float(M0.var()/noise.mean**2),
                expected_basal_input_integral_variance_over_mu2=a["theoretical_integrated_input_variance_over_mu2"],
                basal_input_first_angular_second_moment_over_mu2=float(np.mean(abs(M1)**2)/noise.mean**2),
                expected_basal_input_first_angular_second_moment_over_mu2=a["expected_first_angular_input_second_moment_over_mu2"],
                knots_regenerated_exact=True,regenerated_fine_knots_sha256=noise.latent_sha,
                global_mode_metric="Basal physical quadrature of imposed material field; dynamic C changes are separately reported")
        input_rows.append(base)
    for C0 in (1.,2.,3.):
        for kind,fraction in (("noise",.001),("noise",.01),("mean",.001),("mean",.01),("basal",0.)):
            rr=[r for r in records if r["job"]["resting_circumference"]==C0 and r["job"]["kind"]==kind and r["job"]["fraction"]==fraction]
            n=len(rr);w=rng.multinomial(n,np.full(n,1/n),draws)/n if n>1 else np.empty((0,n))
            base=dict(model="spring_ring",resting_circumference=C0,basal_circumference=rr[0]["circumference"][0],condition=kind,fraction=fraction,n_seeds=n)
            if kind=="noise":
                for lag in (.5,4.,10.,50.):
                    k=int(lag/.5)
                    for gate in (.05,.1,.2):
                        sums=[];counts=[];available=[]
                        for r in rr:
                            valid=(r["R"][k:]>=gate)&(r["R"][:-k]>=gate)
                            a=abs(np.angle(r["P"][k:]*r["P"][:-k].conj()))[valid]
                            sums.append(np.searchsorted(np.sort(a),thresholds,side="right"));counts.append(len(a));available.append(len(valid))
                        sums=np.array(sums);counts=np.array(counts)
                        estimate=ratio(sums.sum(axis=0),counts.sum());d=ratio(w@sums,(w@counts)[:,None]);lo,hi=bounds(d)
                        for i,th in enumerate(thresholds):cdf_rows.append(dict(**base,lag=lag,gate=gate,threshold_radians=th,cdf=estimate[i],ci_lower=lo[i],ci_upper=hi[i],
                            eligible_pair_fraction=float(counts.sum()/sum(available)),bootstrap_draws=draws,uncertainty_unit="whole independent seed"))
            for window,start,stop in (("first_half",0.,1600.),("last_half",1600.,3200.)):
                for gate_name in ("strict","one_region"):
                    names=("basal_arc_squared","midpoint_arc_squared","native_cosine","reference_cosine","delta_C_squared")
                    numerators={name:[] for name in names};denominators=[]
                    for r in rr:
                        take=(r["time"]>=start)&(r["time"]<=stop)
                        q,C,v=r["dominant_region_fraction"][take],r["circumference"][take],r[gate_name][take]
                        ns={name:[] for name in names};counts=[]
                        for k in range(1601):
                            mask=v if k==0 else v[k:]&v[:-k]
                            dq=np.zeros(len(q)) if k==0 else (q[k:]-q[:-k]+.5)%1-.5
                            midpoint=C if k==0 else .5*(C[k:]+C[:-k])
                            dC=np.zeros(len(C)) if k==0 else C[k:]-C[:-k]
                            ds=r["circumference"][0]*dq
                            vals=(ds**2,(midpoint*dq)**2,np.cos(2*np.pi*dq),np.cos(2*np.pi*ds/reference_C),dC**2)
                            for name,values in zip(names,vals):ns[name].append(float(values[mask].sum()))
                            counts.append(int(mask.sum()))
                        denominators.append(counts)
                        for name in names:numerators[name].append(ns[name])
                    den=np.array(denominators)
                    for metric in names:
                        nu=np.array(numerators[metric]);value=ratio(nu.sum(axis=0),den.sum(axis=0));d=ratio(w@nu,w@den) if len(w) else np.empty((0,1601))
                        lo,hi=bounds(d) if len(w) else (np.full(1601,np.nan),np.full(1601,np.nan))
                        for k in range(1601):half_rows.append(dict(**base,window=window,window_start=start,window_stop=stop,gate=gate_name,metric=metric,lag=k*.5,
                            value=value[k],ci_lower=lo[k],ci_upper=hi[k],eligible_pair_fraction=float(den[:,k].sum()/(n*(3201-k))),
                            independent_seeds_with_pairs=int(np.sum(den[:,k]>0)),reference_C=reference_C,same_site_ancestry_assumed=False))
    files=[csvwrite("spring_noise_increment_cdf",cdf_rows),csvwrite("spring_noise_site_half_windows",half_rows),csvwrite("spring_noise_realized_input_audits",input_rows)]
    report=dict(generated_at=datetime.now(timezone.utc).isoformat(),analysis_only=True,new_native_integrations=0,cases=len(records),
        source_sha256=sha(__file__),centroid_audit=centroid,maximum_material_centroid_error=max(c["maximum_material_fraction_error"] for c in centroid),
        exact_noise_regenerations=sum(r.get("knots_regenerated_exact",False) for r in input_rows),
        half_window_site_maxlag=800.,increment_lags=[.5,4.,10.,50.],bootstrap_draws=draws,bootstrap_seed=2026101906,
        outputs=[dict(path=str(p),sha256=sha(p)) for p in files])
    write(META/"spring_noise_supplement_audit.json",report)
    print(json.dumps({k:v for k,v in report.items() if k not in ("centroid_audit","outputs")},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--run",type=Path,required=True);a=p.parse_args();main(a.run.resolve())
