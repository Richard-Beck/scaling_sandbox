"""Incremental canonical-data collection and manual-refresh broad-strength review.

This independent read-only collector never changes frozen native sources or
inputs, submits native cases, or silently counts numerical companions as seeds.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import html
import json
import os
import re
import signal
import subprocess
import sys
import time

for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
sys.dont_write_bytecode=True
import numpy as np
from scipy.signal import find_peaks

ROOT=Path(__file__).resolve().parents[4]
DEFAULT=ROOT/"polarity1d/manuscript/data/runs/revision_noise_strength_broad_20261019"
NAMES=dict(wave_pinning="Wave-pinning",goryachev="Goryachev",holmes_model3="Holmes 3",spring_ring="Spring ring")
STOP_REASON=None


def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
    return h.hexdigest()
def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+".tmp")
    q.write_text(json.dumps(clean(v),indent=2,allow_nan=False)+"\n");q.replace(p)
def csvwrite(p,rr):
    p=Path(p);q=p.with_suffix(p.suffix+".tmp")
    with q.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rr for k in r)))
        w.writeheader();w.writerows(rr)
    q.replace(p)
def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f))


def scheduler(run):
    launch=read(run/"submission.json");handles=[r["job_id"] for r in launch["arrays"]]
    states={};error=None
    try:
        command=["sacct","-j",",".join(handles),"--format=JobID,State,ExitCode,ElapsedRaw","--noheader","--parsable2"]
        result=subprocess.run(command,text=True,capture_output=True,timeout=30,check=True)
        for line in result.stdout.splitlines():
            parts=line.split("|");match=re.fullmatch(r"(\d+)_(\d+)",parts[0])
            if match:states[int(match[2])]=dict(array=match[1],state=parts[1].split()[0],exit_code=parts[2],elapsed_seconds=int(parts[3]))
    except Exception as e:error=f"{type(e).__name__}: {e}"
    return dict(queried_at=now(),handles=handles,states=states,error=error)


def trace(path):
    with np.load(path) as z:
        r={k:z[k].copy() for k in ("time","contrast","current_half_contrast_regions","peaks_prominence_ge_005","M0","q1","q2")}
        r["P"]=z["P_real"]+1j*z["P_imag"]
    r["R"]=abs(r["P"]);r["A"]=r["M0"]/r["M0"][0];r["Q"]=r["q1"]/r["q1"][0]
    r["contrast_ratio"]=r["contrast"]/r["contrast"][0]
    return r


def sample_runs(mask,t):
    d=np.diff(np.r_[False,mask,False].astype(int));a=np.flatnonzero(d==1);b=np.flatnonzero(d==-1)-1
    spans=t[b]-t[a]
    return dict(sample_runs=len(a),run_q90=float(np.quantile(spans,.9)) if len(a) else np.nan,
                longest_observed_span=float(spans.max()) if len(a) else np.nan,
                left_window_censored=bool(len(a) and a[0]==0),right_window_censored=bool(len(b) and b[-1]==len(t)-1))


def describe(r,window,lo,hi,lags):
    take=(r["time"]>=lo)&(r["time"]<=hi);t=r["time"][take];dt=float(t[1]-t[0])
    P=r["P"][take];R=r["R"][take];contrast=r["contrast_ratio"][take]
    multi=r["current_half_contrast_regions"][take]>1;qualified=multi&(contrast>=.1)
    base=dict(window=window,first_time=float(t[0]),last_time=float(t[-1]),observations=len(t),scalar_cadence=dt,
        multiregion_fraction=float(multi.mean()),qualified_multiregion_fraction=float(qualified.mean()),
        qualified_basal001_multiregion_fraction=float(np.mean(multi&(contrast>=.01))),
        weak_R_fraction=float(np.mean(R<.1)),low_contrast_fraction=float(np.mean(contrast<.1)),
        low_contrast_001_fraction=float(np.mean(contrast<.01)),
        qualified_absolute005_multiregion_fraction=float(np.mean(multi&(r["contrast"][take]>=.05))),
        multiple_native005_peaks_fraction=float(np.mean(r["peaks_prominence_ge_005"][take]>1)),
        mean_R=float(R.mean()),mean_A=float(r["A"][take].mean()),mean_Q=float(r["Q"][take].mean()),
        mean_raw_q1=float(r["q1"][take].mean()),mean_q2=float(r["q2"][take].mean()),mean_contrast_ratio=float(contrast.mean()),
        mean_paired_AQ=float(np.mean(r["A"][take]*r["Q"][take])),
        factorization_maximum_error=float(np.max(abs(R-r["A"][take]*r["Q"][take]))))
    for label,mask in (("multiregion",multi),("qualified_multiregion",qualified)):
        base.update({label+"_"+k:v for k,v in sample_runs(mask,t).items()})
    memory=[]
    for lag in lags:
        k=int(round(lag/dt));assert abs(k*dt-lag)<1e-7
        a,b=P[k:],P[:-k];gate=(abs(a)>=.1)&(abs(b)>=.1)
        n=int(gate.sum());available=len(gate)
        cosine=float(np.sum(np.cos(np.angle(a[gate]*b[gate].conj())))) if n else 0.
        memory.append(dict(window=window,lag=lag,pairs_available=available,pairs_eligible=n,
            eligible_fraction=float(n/available),direction_cosine_sum=cosine,
            direction_correlation=cosine/n if n else np.nan,
            support=n>=100 and n/available>=.01,
            raw_product_sum=float(np.real(np.vdot(b,a))),energy_a=float(np.sum(abs(a)**2)),energy_b=float(np.sum(abs(b)**2))))
    return base,memory


def profile_audit(path,r):
    with np.load(path) as z:
        t,u=z["time"],z["primary"];C=np.broadcast_to(z["circumference"],len(t))
        q=z["material_positions"] if "material_positions" in z else z["x"]/C[0]
    ix=np.searchsorted(r["time"],t);np.testing.assert_allclose(r["time"][ix],t,rtol=0,atol=1e-9)
    contrast=np.ptp(u,axis=1);above=u>=u.min(axis=1)[:,None]+.5*contrast[:,None]
    count=np.sum(above&~np.roll(above,1,axis=1),axis=1);count[above.all(axis=1)]=1;count[contrast<=1e-12]=0
    np.testing.assert_array_equal(count,r["current_half_contrast_regions"][ix])
    M0=C*np.mean(u,axis=1);M1=C*(u@np.exp(2j*np.pi*q))/len(q)
    eP=float(np.max(abs(M1/M1[0]-r["P"][ix])));e0=float(np.max(abs(M0-r["M0"][ix])))
    assert max(eP,e0)<1e-9
    relative_peaks=[]
    for field in u:
        peaks,_=find_peaks(np.tile(field,3),prominence=.05*contrast[0]);relative_peaks.append(np.count_nonzero((peaks>=len(q))&(peaks<2*len(q))))
    return dict(profiles=len(t),profile_cadence=float(t[1]-t[0]),max_physical_P_error=eP,max_physical_M0_error=e0,
                relative_basal005_multiple_peak_profile_fraction=float(np.mean(np.array(relative_peaks)>1)))


def summarize_case(base,tracepath,profilepath,existing_profile_audit=None):
    r=trace(tracepath);assert len(r["time"])==6401
    assert np.isfinite(r["P"]).all() and np.isfinite(r["M0"]).all()
    platform=base["platform"];T=float(r["time"][-1]);lags=[5.,40.,100.,500.,2000.,8000.] if platform=="chemical" else [.5,4.,10.,50.,200.,800.]
    descriptor=[];memory=[]
    for name,lo,hi in (("full",0,T),("first_half",0,T/2),("last_half",T/2,T)):
        row,mm=describe(r,name,lo,hi,lags);descriptor.append(dict(**base,**row));memory.extend(dict(**base,**v) for v in mm)
    assert max(v["factorization_maximum_error"] for v in descriptor)<1e-9
    return dict(descriptors=descriptor,memory=memory,profile_audit=existing_profile_audit or profile_audit(profilepath,r))


def numeric_checks(run,m):
    results=[];pending=[]
    means=("mean_R","mean_A","mean_Q","mean_contrast_ratio")
    occupancy=("multiregion_fraction","qualified_multiregion_fraction","weak_R_fraction","low_contrast_fraction")
    for tight in [j for j in m["jobs"] if j["variant"]=="tight_ode"]:
        standard=m["jobs"][tight["paired_standard_slot"]];folders=[run/"workers"/f"{j['slot']:04d}" for j in (standard,tight)]
        if not all((p/"status.json").exists() and read(p/"status.json").get("complete") for p in folders):
            pending.append(tight["model"]);continue
        audits=[read(p/"noise_audit.json") for p in folders]
        key="gaussian_knots_sha256" if tight["platform"]=="chemical" else "latent_fine_knots_sha256"
        assert audits[0][key]==audits[1][key]
        standard_coarse,standard_dense,tight_coarse,tight_dense=[trace(p/name) for p,name in
            ((folders[0],"trace.npz"),(folders[0],"dense_trace.npz"),(folders[1],"trace.npz"),(folders[1],"dense_trace.npz"))]
        for p in folders:
            with np.load(p/"trace.npz") as a,np.load(p/"dense_trace.npz") as b:
                for key2 in ("time","P_real","P_imag","M0","contrast","current_half_contrast_regions"):
                    np.testing.assert_allclose(a[key2],b[key2][::10],rtol=0,atol=1e-9 if key2=="time" else 0)
        comparisons=(("ODE_canonical",standard_coarse,tight_coarse),("observation_standard",standard_coarse,standard_dense),("observation_tight",tight_coarse,tight_dense))
        T=m["native_protocols"][tight["platform"]]["horizon"];lags=m["numerical_check"]["lag_checkpoints"][tight["platform"]]
        for label,a,b in comparisons:
            for name,lo,hi in (("full",0,T),("first_half",0,T/2),("last_half",T/2,T)):
                ar,am=describe(a,name,lo,hi,lags);br,bm=describe(b,name,lo,hi,lags)
                base=dict(model=tight["model"],size=tight["size"],seed=tight["seed"],check=label,window=name,latent_sha256=audits[0][key])
                for metric in means+occupancy:
                    difference=abs(br[metric]-ar[metric]);tol=.02 if metric in occupancy else max(.02,.02*abs(ar[metric]))
                    results.append(dict(**base,metric=metric,lag=None,reference=ar[metric],comparison=br[metric],absolute_difference=difference,tolerance=tol,status="pass" if difference<=tol else "fail"))
                for metric in ("multiregion_run_q90","qualified_multiregion_run_q90"):
                    av,bv=ar[metric],br[metric];both=np.isfinite(av) and np.isfinite(bv)
                    difference=abs(bv-av) if both else np.nan
                    tol=max(2*m["native_protocols"][tight["platform"]]["scalar_interval"],.2*max(av,bv)) if both else np.nan
                    state=("pass" if difference<=tol else "fail") if both else "unavailable_both" if not np.isfinite(av) and not np.isfinite(bv) else "availability_changed"
                    results.append(dict(**base,metric=metric,lag=None,reference=av,comparison=bv,absolute_difference=difference,tolerance=tol,status=state))
                for aa,bb in zip(am,bm):
                    difference=abs(aa["eligible_fraction"]-bb["eligible_fraction"])
                    results.append(dict(**base,metric="direction_eligible_fraction",lag=aa["lag"],reference=aa["eligible_fraction"],comparison=bb["eligible_fraction"],absolute_difference=difference,tolerance=.02,status="pass" if difference<=.02 else "fail"))
                    supported=aa["support"] and bb["support"];difference=abs(aa["direction_correlation"]-bb["direction_correlation"]) if supported else np.nan
                    state=("pass" if difference<=.05 else "fail") if supported else "unavailable_both" if not aa["support"] and not bb["support"] else "availability_changed"
                    results.append(dict(**base,metric="direction_correlation",lag=aa["lag"],reference=aa["direction_correlation"],comparison=bb["direction_correlation"],absolute_difference=difference,tolerance=.05,status=state))
    return results,pending


def collect(run):
    run=Path(run).resolve();m=read(run/"manifest.json");out=run/"results/partial";out.mkdir(parents=True,exist_ok=True)
    cachepath=out/"case_cache.json";cache=read(cachepath) if cachepath.exists() else {}
    sched=scheduler(run);sources=[];worker_inventory=[]
    for old in rows(run/"results/exact_reuse_cases.csv"):
        base=dict(case_id="legacy_"+old["platform"]+"_"+old["legacy_slot"],model=old["model"],platform=old["platform"],size=float(old["size"]),fraction=float(old["fraction"]),kind=old["kind"],seed=int(old["seed"]),legacy=True)
        sources.append((base,Path(old["trace"]),Path(old["primary"]),old["trace_sha256"],old["primary_sha256"]))
    for j in m["jobs"]:
        folder=run/"workers"/f"{j['slot']:04d}";p=folder/"status.json";s=read(p) if p.exists() else {}
        ss=sched["states"].get(j["slot"],{})
        worker_inventory.append(dict(slot=j["slot"],model=j["model"],size=j["size"],fraction=j["fraction"],kind=j["kind"],variant=j["variant"],
            native_status=s.get("status","not_started"),producer=s.get("producer","not_started"),complete=s.get("complete",False),
            scheduler_state=ss.get("state","not_reported"),scheduler_exit_code=ss.get("exit_code"),elapsed_seconds=s.get("elapsed_seconds",ss.get("elapsed_seconds")),error=s.get("error",s.get("message"))))
        if s.get("complete") and j["variant"]=="standard":
            base=dict(case_id=f"new_{j['slot']:04d}",model=j["model"],platform=j["platform"],size=j["size"],fraction=j["fraction"],kind=j["kind"],seed=j["seed"],legacy=False)
            sources.append((base,folder/"trace.npz",folder/"primary.npz",s["output_sha256"]["trace.npz"],s["output_sha256"]["primary.npz"]))
    changed=False
    for base,tp,pp,th,ph in sources:
        key=base["case_id"]
        if key not in cache:
            assert sha(tp)==th and sha(pp)==ph
            item=summarize_case(base,tp,pp);item.update(trace=str(tp),trace_sha256=th,profile=str(pp),profile_sha256=ph)
            cache[key]=clean(item);changed=True
        else:
            assert cache[key]["trace_sha256"]==th and cache[key]["profile_sha256"]==ph
            if "qualified_basal001_multiregion_fraction" not in cache[key]["descriptors"][0]:
                item=summarize_case(base,tp,pp,cache[key]["profile_audit"])
                item.update(trace=str(tp),trace_sha256=th,profile=str(pp),profile_sha256=ph)
                cache[key]=clean(item);changed=True
    if changed:write(cachepath,cache)
    descriptors=[v for c in cache.values() for v in c["descriptors"]];memory=[v for c in cache.values() for v in c["memory"]]
    conditions=[]
    exclude={"size","fraction","seed","legacy","observations","scalar_cadence","first_time","last_time"}
    metrics=sorted({k for r in descriptors for k,v in r.items() if isinstance(v,(float,int)) and not isinstance(v,bool) and k not in exclude}
                   |{"multiregion_run_q90","qualified_multiregion_run_q90","multiregion_longest_observed_span","qualified_multiregion_longest_observed_span"})
    for model in NAMES:
        sizes=(1.,2.,3.) if model=="spring_ring" else (20.,40.,60.)
        for size in sizes:
            for kind,f in [(kind,f) for f in m["plan"]["grid"] for kind in ("noise","mean")]+[("basal",0.)]:
                target=12 if f in (.001,.01) and kind=="noise" else 2 if kind=="noise" else 1
                for window in ("full","first_half","last_half"):
                    rr=[r for r in descriptors if r["model"]==model and r["size"]==size and r["kind"]==kind and r["fraction"]==f and r["window"]==window]
                    row=dict(model=model,size=size,fraction=f,kind=kind,window=window,available_cases=len(rr),target_cases=target,
                        coverage="complete" if len(rr)==target else "partial" if rr else "not_available",
                        inference="exploratory2seedtarget" if kind=="noise" and target==2 else "retained12seeds" if kind=="noise" else "deterministic_control")
                    for metric in metrics:
                        vals=[r[metric] for r in rr if r.get(metric) is not None]
                        row[metric]=float(np.mean(vals)) if vals else np.nan
                        row[metric+"_seed_min"]=float(min(vals)) if vals else np.nan;row[metric+"_seed_max"]=float(max(vals)) if vals else np.nan
                    conditions.append(row)
    numerical,pending=numeric_checks(run,m)
    csvwrite(out/"case_window_descriptors.csv",descriptors);csvwrite(out/"case_window_memory.csv",memory)
    profile_readouts=[dict(**{k:c["descriptors"][0][k] for k in ("case_id","model","size","fraction","kind","seed")},
                           **c["profile_audit"],source=c["profile"],source_sha256=c["profile_sha256"]) for c in cache.values()]
    csvwrite(out/"case_profile_readouts.csv",profile_readouts)
    csvwrite(out/"condition_window_summaries.csv",conditions);csvwrite(out/"worker_inventory.csv",worker_inventory)
    if numerical:csvwrite(out/"upper_strength_numerical_checks.csv",numerical)
    failures=[r for r in worker_inventory if r["native_status"]=="failed" or r["scheduler_state"] in ("FAILED","TIMEOUT","OUT_OF_MEMORY","CANCELLED","NODE_FAIL")]
    terminal=all(r["producer"]=="terminal" for r in worker_inventory) and not any(r["scheduler_state"] in ("RUNNING","PENDING","COMPLETING","CONFIGURING") for r in worker_inventory)
    status=dict(generated_at=now(),data_observed_at=sched["queried_at"],producer="terminal_with_failures" if terminal and failures else "terminal" if terminal else "running",
        collector_state="running" if os.environ.get("SLURM_JOB_ID") else "one_shot",collector_job_id=os.environ.get("SLURM_JOB_ID"),
        collector_observation="Native state is observed at the stated time; this static page cannot claim current state if the collector stops",
        complete=terminal and not failures,new_completed_outputs=sum(bool(r["complete"]) for r in worker_inventory),requested_new_outputs=256,
        scientific_cases_available=len(cache),legacy_cases=324,new_scientific_cases=len(cache)-324,numerical_companions_complete=sum(r["complete"] and r["variant"]=="tight_ode" for r in worker_inventory),
        failures=len(failures),scheduler=sched,manifest_sha256=sha(run/"manifest.json"),collector_sha256=sha(__file__),
        upper_checks_pending=pending,upper_defined_failed=sum(r["status"]=="fail" for r in numerical),
        upper_availability_changes=sum(r["status"]=="availability_changed" for r in numerical),
        numerical_qualification="pending_root_and_independent_review",replication_started=False,
        profile_replay=dict(cases=len(cache),profiles=sum(c["profile_audit"]["profiles"] for c in cache.values()),maximum_P_error=max(c["profile_audit"]["max_physical_P_error"] for c in cache.values())))
    write(out/"status.json",status)
    render_html(run,out,status,conditions)
    return status


def render_html(run,out,status,conditions):
    esc=html.escape
    def value(r,k,pct=False):
        v=r.get(k);return "—" if v is None or not np.isfinite(v) else f"{100*v:.2f}%" if pct else f"{v:.3f}"
    parts=["<!doctype html><html><head><meta charset='utf-8'><title>Broad-strength pilot: partial observations</title><style>body{font:16px system-ui;margin:2rem auto;max-width:1180px;padding:0 1rem;color:#23303a}table{border-collapse:collapse;width:100%;font-size:14px;margin:1rem 0 2rem}td,th{border-bottom:1px solid #d9dfe3;padding:.45rem;text-align:right}th:first-child,td:first-child{text-align:left}.notice{background:#fff1cd;padding:1rem;line-height:1.5}.meta{color:#53616d}h2{margin-top:2rem}</style></head><body>",
        "<h1>Broad-strength noisy-ring pilot</h1>",
        f"<p class='meta'>Data generated {esc(status['generated_at'])} · last observed native state: {esc(status['producer'])} at {esc(status['data_observed_at'])} · collector: {esc(status.get('collector_state','unknown'))} · native arrays {', '.join(status['scheduler']['handles'])}</p>",
        f"<p>{status['new_completed_outputs']}/256 new outputs complete; {status['scientific_cases_available']} scientific cases available, including all324 retained cases. Failures: {status['failures']}. Upper-strength numerical cases pending: {esc(', '.join(status['upper_checks_pending']) or 'none')}.</p>",
        "<p class='notice'>New strengths target only two independent seeds; the original0.001/0.01 levels retain12. These are incomplete or exploratory observations, not final frequencies or discrete phase boundaries. Increasing f scales positive input mean and fluctuation amplitude together. Tables use the full native window and canonical cadence. Reload manually for newer data while the collector runs; if it stops, this page remains a timestamped observation, not live native status. No timed refresh.</p>",
        "<p><a href='condition_window_summaries.csv'>All windows and seed ranges</a> · <a href='case_window_descriptors.csv'>Every scalar case</a> · <a href='case_profile_readouts.csv'>Profile-cadence prominence sensitivity</a> · <a href='worker_inventory.csv'>Native/scheduler status</a> · <a href='status.json'>Exact status</a></p>"]
    for model,name in NAMES.items():
        for size in ((1.,2.,3.) if model=="spring_ring" else (20.,40.,60.)):
            parts.append(f"<h2>{name}, {'C₀' if model=='spring_ring' else 'C'}={size:g}</h2><table><tr><th>f</th><th>Noisy n/target</th><th>Multiple regions</th><th>Mean R</th><th>Mean M₀ ratio</th><th>Mean q₁ ratio</th><th>Low contrast</th><th>Matched mean R</th></tr>")
            for f in [1e-4,3e-4,.001,.003,.01,.03,.1,.3,1.]:
                r,=[r for r in conditions if r['model']==model and r['size']==size and r['fraction']==f and r['kind']=='noise' and r['window']=='full']
                c,=[r for r in conditions if r['model']==model and r['size']==size and r['fraction']==f and r['kind']=='mean' and r['window']=='full']
                parts.append(f"<tr><td>{f:g}</td><td>{r['available_cases']}/{r['target_cases']}</td><td>{value(r,'multiregion_fraction',True)}</td><td>{value(r,'mean_R')}</td><td>{value(r,'mean_A')}</td><td>{value(r,'mean_Q')}</td><td>{value(r,'low_contrast_fraction',True)}</td><td>{value(c,'mean_R')}</td></tr>")
            parts.append("</table>")
    parts.append("<p>Multiple regions use the current half-contrast profile; low contrast means below10% of the initial basal contrast. Amplitude-qualified region counts and native/relative prominence sensitivities are exported separately. R is relative first-moment magnitude, distinct from amount and coherence. Sample-run spans are sampled descriptors, not event durations. Native units and amplitudes are not calibrated across models.</p></body></html>")
    tmp=out/"index.html.tmp";tmp.write_text("\n".join(parts));tmp.replace(out/"index.html")


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("command",choices=("collect","monitor"));p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--interval",type=int,default=60)
    a=p.parse_args()
    if a.command=="monitor":assert os.environ.get("SLURM_JOB_ID"),"Long-running monitor uses SLURM"
    def stop(signum,frame):globals()["STOP_REASON"]="signal_"+str(signum)
    if a.command=="monitor":
        for sig in (signal.SIGUSR1,signal.SIGTERM,signal.SIGINT):signal.signal(sig,stop)
    s=None;error=None
    try:
        while STOP_REASON is None:
            s=collect(a.run);print(json.dumps({k:v for k,v in s.items() if k not in ("scheduler",)},indent=2),flush=True)
            if a.command=="collect" or s["producer"].startswith("terminal"):break
            for _ in range(a.interval):
                if STOP_REASON is not None:break
                time.sleep(1)
    except Exception as exc:
        error=f"{type(exc).__name__}: {exc}";raise
    finally:
        if a.command=="monitor":
            out=Path(a.run)/"results/partial"
            if s is None and (out/"status.json").exists():s=read(out/"status.json")
            if s is not None:
                s["collector_state"]="stopped_native_terminal" if s["producer"].startswith("terminal") else "stopped_error" if error else "stopped_"+(STOP_REASON or "requested")
                s["collector_stopped_at"]=now();s["collector_stop_reason"]=error or STOP_REASON or "Native arrays terminal"
                write(out/"status.json",s)
                # Preserve the already rendered data and exact observation time.
                page=out/"index.html"
                if page.exists():
                    text=page.read_text().replace("· collector: running ·","· collector: "+html.escape(s["collector_state"])+" ·")
                    notice="<p class='notice'>Collector stopped at "+html.escape(s["collector_stopped_at"])+". Native state above is the last timestamped observation; "+html.escape(s["collector_stop_reason"])+".</p>"
                    temp=out/"index.html.tmp";temp.write_text(text.replace("</body>",notice+"</body>"));temp.replace(page)


if __name__=="__main__":main()
