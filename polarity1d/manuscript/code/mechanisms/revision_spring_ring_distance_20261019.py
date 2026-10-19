"""Bounded Spring-ring distance x finite-duration experiment.

Exact shared basals, fixed physical input support, complete coupled histories.
Resting circumference, basal circumference and evolving physical distances are
separate observables. Every integration is a frozen-source SLURM worker.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import time
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from spring_ring_model_20261019 import Closure, physical_raised_cosine

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"data/runs/revision_spring_distance_duration_20261019"
BASALS=ROOT/"data/runs/revision_ring_spring_20261019"
PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"


def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,z):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_suffix(p.suffix+".tmp");temp.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");temp.replace(p)


def core_weights(q,u):
    peak=int(u.argmax());origin=q[peak];s=(q-origin+.5)%1-.5;order=np.argsort(s);x,p=s[order],u[order]
    level=p.min()+.5*np.ptp(p);inside=np.flatnonzero(p>=level);assert np.all(np.diff(inside)==1)
    a,b=int(inside[0]),int(inside[-1]);assert a>0 and b<len(q)-1
    left=x[a-1]+(level-p[a-1])/(p[a]-p[a-1])*(x[a]-x[a-1])
    right=x[b]+(level-p[b])/(p[b+1]-p[b])*(x[b+1]-x[b])
    dq=1/len(q);weights=np.maximum(0,np.minimum(s+dq/2,right)-np.maximum(s-dq/2,left))
    return dict(weights=weights,left=float(left),right=float(right),origin=float(origin),peak_cell=peak,
                fraction_width=float(weights.sum()),background=float(p.min()))


def launch(run):
    if run.exists():raise FileExistsError("Fresh run required; retained output immutable")
    summary=read(BASALS/"results/preparation_summary.json");assert summary["complete"] and summary["qualified"]==3
    independent=BASALS/"results/root_independent_preparation_holdout.json";assert read(independent)["passed"]
    for directory in ("inputs","frozen_sources","workers","logs","results"):(run/directory).mkdir(parents=True)
    paths=[Path(__file__),Path(__file__).with_name("spring_ring_model_20261019.py"),Path(__file__).with_name("spring_ring_closure_20261017.py"),ROOT/"code/core/models/spring.py",ROOT/"code/core/models/__init__.py"]
    sources={}
    for path in paths:
        relative=path.relative_to(ROOT);target=run/"frozen_sources"/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);sources[str(relative)]=sha(target)
    jobs=[];exclusions=[];inputs={}
    for C0 in (1.,2.,3.):
        status=read(BASALS/"preparations"/f"C{C0:g}"/"status.json");assert status["status"]=="qualified"
        original=Path(status["checkpoint"]);assert sha(original)==status["checkpoint_sha256"]
        target=run/"inputs"/f"C{C0:g}.npz";shutil.copy2(original,target);inputs[str(target.relative_to(run))]=sha(target)
        with np.load(target) as z:state=z["state"];q=z["material_positions"];n=len(q)
        Cb=float(state[-1]);core=core_weights(q,state[:n])
        for fraction in (.30,.35,.40,.45,.50):
            center=(core["origin"]+fraction)%1;gap=fraction*Cb-core["right"]*Cb-.1
            record=dict(model="spring",resting_circumference=C0,basal_circumference=Cb,cells=n,fraction=fraction,
                cue_center_fraction=center,physical_distance_at_onset=fraction*Cb,normal_angle_deg=360*fraction,
                basal_core_physical_width=core["fraction_width"]*Cb,initial_short_edge_gap=gap,
                input=str(target.relative_to(run)),input_sha256=inputs[str(target.relative_to(run))])
            for duration in (.1,1.,10.,100.,1000.):
                if gap<=0:
                    exclusions.append(dict(**record,cue_duration=duration,status="initial_geometric_overlap_exclusion"));continue
                jobs.append(dict(**record,slot=len(jobs),cue_duration=duration))
    manifest=dict(generated_at=now(),experiment="Three-size Spring native distance-duration pilot",jobs=jobs,exclusions=exclusions,
        fractions=[.30,.35,.40,.45,.50],durations=[.1,1.,10.,100.,1000.],cue_peak=4.,cue_full_physical_support=.2,
        recovery=13000.,solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),source_sha256=sources,input_sha256=inputs,
        shared_basals_root=str(BASALS),independent_preparation_holdout=str(independent),independent_preparation_holdout_sha256=sha(independent),
        interpretation="New mean-force global-extension ring analogue; prescribed material placement, time-varying physical arc distances; finite sampled outcomes, not unique minimum duration or native2Dmechanics",
        requested_conditions=len(jobs)+len(exclusions),native_jobs=len(jobs))
    write(run/"manifest.json",manifest)
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job=subprocess.check_output(["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:30:00",f"--array=0-{len(jobs)-1}",
        "--job-name=Spring-distance-duration","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command],text=True).strip()
    write(run/"launch.json",dict(generated_at=now(),array_job_id=job,manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(dict(array_job_id=job,native_jobs=len(jobs),exclusions=len(exclusions),run=str(run)),indent=2))


def samples(end):
    return np.unique(np.r_[0.,np.geomspace(1e-4,min(.1,end),25),np.arange(.1,min(10.,end),.1),
        np.arange(10.,min(100.,end),1.),np.arange(100.,min(1000.,end),5.),np.arange(1000.,end,50.),end])


def crossing(t,values,threshold):
    hit=np.flatnonzero(values<=threshold)
    if not len(hit):return dict(status="right_censored",lower=float(t[-1]),upper=None)
    i=int(hit[0]);return dict(status="observed_saved_sample_crossing",lower=float(t[max(0,i-1)]),upper=float(t[i]))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Native integration requires SLURM")
    manifest=read(run/"manifest.json");job=manifest["jobs"][slot]
    for relative,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/relative)==digest
    for relative,digest in manifest["input_sha256"].items():assert sha(run/relative)==digest
    directory=run/"workers"/f"{slot:03d}";directory.mkdir()
    record=dict(generated_at=now(),producer="running",status="challenge",complete=False,job=job,slurm_job_id=os.environ["SLURM_JOB_ID"],phases=[])
    write(directory/"status.json",record);began=time.monotonic()
    with np.load(run/job["input"]) as z:state=z["state"].copy();q=z["material_positions"].copy()
    initial=state.copy();n=job["cells"];model=Closure(job["resting_circumference"],n);amount0=model.amount(state)
    core=core_weights(q,state[:n]);parts=[]
    for name,duration,peak,offset in (("challenge",job["cue_duration"],4.,0.),("recovery",13000.,0.,job["cue_duration"])):
        record.update(status=name,generated_at=now());write(directory/"status.json",record)
        def cue(y):return physical_raised_cosine(model,y,job["cue_center_fraction"],peak,.2) if peak else (0.,0.)
        def fun(t,y):return model.rhs(y,cue(y)[0])
        def jac(t,y):return model.jacobian(y,*cue(y))
        solution=solve_ivp(fun,(0.,duration),state,method="BDF",rtol=1e-9,atol=1e-11,max_step=10.,t_eval=samples(duration),jac=jac)
        if not solution.success:raise RuntimeError(solution.message)
        states=solution.y.T;amounts=states[:,-1]*np.mean(states[:,:n]+states[:,n:-1],axis=1)
        mass_error=float(abs(amounts/amount0-1).max());minimum=float(states[:,:-1].min())
        assert np.isfinite(states).all() and mass_error<=1e-7 and minimum>=-5e-10
        archive=directory/(name+".npz")
        np.savez_compressed(archive,time=solution.t,states=states,material_positions=q,cue_duration=job["cue_duration"],
            cue_center_fraction=job["cue_center_fraction"],cue_peak=peak,physical_amount=amounts)
        record["phases"].append(dict(kind=name,path=str(archive.relative_to(run)),sha256=sha(archive),samples=len(solution.t),
            max_amount_relative_error=mass_error,minimum_concentration=minimum,nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu))
        parts.append((solution.t+offset,states));state=states[-1].copy()
    np.testing.assert_array_equal(parts[0][1][-1],parts[1][1][0])
    t=np.r_[parts[0][0],parts[1][0][1:]];y=np.vstack((parts[0][1],parts[1][1][1:]));C=y[:,-1];u=y[:,:n]
    oldraw=C*(u@core["weights"]);background_amount0=initial[-1]*core["background"]*core["weights"].sum()
    E0=oldraw[0]-background_amount0;excess=(oldraw-background_amount0)/E0
    frozen_concentration_excess=C*((u-core["background"])@core["weights"])/E0
    theta=2*np.pi*(q-core["origin"]);M0=C*u.mean(axis=1);M1=C*(u@np.exp(1j*theta))/n;M2=C*(u@np.exp(2j*theta))/n
    P=M1/M1[0];R=abs(P);angles=np.rad2deg(np.angle(P));idxoff=len(parts[0][0])-1
    cue_difference=np.rad2deg(np.angle(np.exp(1j*np.deg2rad(angles-360*job["fraction"]))))
    initial_patch=np.maximum(0,np.minimum((q-job["cue_center_fraction"]+.5)%1-.5+1/(2*n),.1/initial[-1])-np.maximum((q-job["cue_center_fraction"]+.5)%1-.5-1/(2*n),-.1/initial[-1]))
    cue_fixed_material=C*(u@initial_patch);cue_excess=(cue_fixed_material-initial[-1]*core["background"]*initial_patch.sum())/E0
    scales=np.r_[np.full(2*n,2.),job["resting_circumference"]];tail=t>=t[-1]-1000
    tail_range=float(np.max(np.ptp(y[tail],axis=0)/scales));rhs=float(np.max(abs(model.rhs(y[-1],0.))/scales)*1000)
    final_u=u[-1];peaks,_=find_peaks(np.tile(final_u,3),prominence=.05);peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
    half=final_u>=final_u.min()+.5*np.ptp(final_u);components=int(np.sum(half&~np.roll(half,1)))
    output={**job,"status":"complete","cue_free_time":13000.,"physical_C_min":float(C.min()),"physical_C_max":float(C.max()),"physical_C_final":float(C[-1]),
        "physical_distance_min":float(job["fraction"]*C.min()),"physical_distance_max":float(job["fraction"]*C.max()),
        "old_core_excess_at_withdrawal":float(excess[idxoff]),"old_core_excess_final":float(excess[-1]),"old_core_raw_final":float(oldraw[-1]/oldraw[0]),
        "old_core_frozen_concentration_excess_final":float(frozen_concentration_excess[-1]),"cue_fixed_material_excess_final":float(cue_excess[-1]),
        "old_core_excess_minimum":float(excess.min()),"old_core_excess_rebound":float(excess[-1]-excess[idxoff]),
        "R_at_withdrawal":float(R[idxoff]),"R_final":float(R[-1]),"q1_final":float(abs(M1[-1])/M0[-1]),"q2_final":float(abs(M2[-1])/M0[-1]),
        "angle_at_withdrawal":float(angles[idxoff]) if R[idxoff]>=.1 else None,"angle_final":float(angles[-1]) if R[-1]>=.1 else None,
        "cue_angle_error_final":float(cue_difference[-1]) if R[-1]>=.1 else None,"cue_aligned_final_15deg":bool(abs(cue_difference[-1])<=15 and R[-1]>=.1),
        "old_aligned_final_15deg":bool(abs(angles[-1])<=15 and R[-1]>=.1),"resolved_peaks_final":len(peaks),"half_contrast_components_final":components,
        "tail_1000_scaled_range":tail_range,"rhs_projected_1000":rhs,"tail_stationary":bool(tail_range<=1e-4 and rhs<=1e-4),
        "maximum_amount_relative_error":max(p["max_amount_relative_error"] for p in record["phases"])}
    for threshold in (.9,.5,.1):
        tag=str(threshold).replace(".","p");event=crossing(t,excess,threshold);during=crossing(t[:idxoff+1],excess[:idxoff+1],threshold)
        output.update({f"loss_{tag}_{k}":v for k,v in event.items()});output.update({f"cue_on_loss_{tag}_{k}":v for k,v in during.items()})
    np.savez_compressed(directory/"readout_traces.npz",time=t,physical_C=C,old_core_raw_amount=oldraw,
        old_core_transported_background_excess=excess,old_core_frozen_concentration_excess=frozen_concentration_excess,
        cue_fixed_material_excess=cue_excess,M0=M0,M1_real=M1.real,M1_imag=M1.imag,P_real=P.real,P_imag=P.imag,R=R,q1=abs(M1)/M0,q2=abs(M2)/M0,
        angle_deg_unmasked=angles,physical_distance=job["fraction"]*C)
    write(directory/"readouts.json",output)
    record.update(generated_at=now(),producer="finished",status="complete",complete=True,readouts=str(directory/"readouts.json"),
        readouts_sha256=sha(directory/"readouts.json"),traces_sha256=sha(directory/"readout_traces.npz"),elapsed_seconds=time.monotonic()-began)
    write(directory/"status.json",record);print(json.dumps(output,indent=2))


def collect(run):
    manifest=read(run/"manifest.json");rows=list(manifest["exclusions"]);statuses=[]
    for job in manifest["jobs"]:
        directory=run/"workers"/f'{job["slot"]:03d}';path=directory/"status.json"
        record=read(path) if path.exists() else dict(status="not_started",complete=False)
        statuses.append(record)
        if record.get("complete") and record["status"]=="complete":
            assert sha(directory/"readouts.json")==record["readouts_sha256"]
            row=read(directory/"readouts.json")
            for phase in record["phases"]:
                prefix="actual_cue" if phase["kind"]=="challenge" else "actual_recovery"
                row[prefix+"_path"]=str(run/phase["path"])
                row[prefix+"_sha256"]=phase["sha256"]
            rows.append(row)
        else:rows.append(dict(**job,status=record["status"]))
    frame=pd.DataFrame(rows);frame.to_csv(run/"results/distance_duration_readouts.csv",index=False)
    complete=all(s.get("complete") and s["status"]=="complete" for s in statuses)
    write(run/"results/collection.json",dict(generated_at=now(),complete=complete,native_complete=sum(s.get("complete",False) for s in statuses),
        requested=len(rows),excluded=len(manifest["exclusions"]),states={key:sum(s["status"]==key for s in statuses) for key in set(s["status"] for s in statuses)}))
    print(json.dumps(read(run/"results/collection.json"),indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=["launch","worker","collect"])
    parser.add_argument("--run",type=Path,default=DEFAULT);parser.add_argument("--slot",type=int);a=parser.parse_args()
    if a.action=="launch":launch(a.run.resolve())
    elif a.action=="collect":collect(a.run.resolve())
    else:worker(a.run.resolve(),a.slot)
