"""Frozen first stage: qualify five newly used intermediate-size ring basals.

Three C35 chemical basals retain their existing exact native preparation;
two Spring basals are prepared from fresh homogeneous native states via SLURM.
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
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[name]="1"
import numpy as np
from scipy.signal import find_peaks
from spring_ring_model_20261019 import Closure, CuePhase, run_protocol

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"data/runs/revision_arc_narrow_basals_20261019"
PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS=("wave_pinning","goryachev","holmes_model3")


def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+".tmp")
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");tmp.replace(p)


def freeze(run):
    if run.exists():raise FileExistsError("Fresh immutable basal stage required")
    for folder in ("inputs","frozen_sources","workers","basals","logs","results"):(run/folder).mkdir(parents=True)
    reference=ROOT/"data/runs/revision_standard_angles_sizes_20261019"
    accepted=read(reference/"manifest.json")
    sources={}
    for rel,digest in accepted["source_sha256"].items():
        original=reference/"frozen_sources"/rel;assert sha(original)==digest
        target=run/"frozen_sources"/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(original,target);sources[rel]=digest
    relative=Path(__file__).relative_to(ROOT);target=run/"frozen_sources"/relative;shutil.copy2(__file__,target);sources[str(relative)]=sha(target)
    for name in ("spring_ring_model_20261019.py","spring_ring_closure_20261017.py"):
        assert sha(Path(__file__).with_name(name))==sources["code/mechanisms/"+name]
    old=ROOT/"data/runs/chemical_ring_functional_v2_20261017";oldmanifest=read(old/"manifest.json")
    reused=[];inputs={}
    for model in MODELS:
        job=next(j for j in oldmanifest["jobs"] if j["model"]==model and j["circumference"]==35)
        preparation=old/job["preparation"];statuspath=preparation/"workers/00/status.json";status=read(statuspath)
        assert status["complete"] and status["status"]=="qualifying_single_site"
        assert any(a.get("verified_settled") for a in status["attempts"])
        original=preparation/status["checkpoint"];assert sha(original)==status["checkpoint_sha256"]
        target=run/"basals"/f"{model}_C35.npz";shutil.copy2(original,target)
        reused.append(dict(model=model,size_parameter=35.,execution="exact_existing_preparation_reuse",input=str(target.relative_to(run)),
                           input_sha256=sha(target),original_path=str(original),original_sha256=sha(original),
                           qualification_status_path=str(statuspath),qualification_status_sha256=sha(statuspath),
                           preparation_manifest=str(preparation/"manifest.json"),preparation_manifest_sha256=sha(preparation/"manifest.json"),
                           fixed_scales=job["fixed_scales"]))
        inputs[str(target.relative_to(run))]=sha(target)
    jobs=[]
    for size in (2.5,2.75):
        n=int(round(100*size));model=Closure(size,n);initial=model.initial()
        path=run/"inputs"/f"spring_C{size:g}_homogeneous.npz"
        np.savez_compressed(path,state=initial,material_positions=model.material_positions,resting_circumference=size)
        inputs[str(path.relative_to(run))]=sha(path)
        jobs.append(dict(slot=len(jobs),model="spring",size_parameter=size,resting_circumference=size,cells=n,
                         preparation_axis_fraction=float(model.material_positions[0]),input=str(path.relative_to(run)),input_sha256=sha(path)))
    manifest=dict(generated_at=now(),experiment="New intermediate-size basals for the authorized narrow-range stationary-cue extension",
                  requested_basals=5,exact_existing_chemical_preparations=reused,jobs=jobs,source_sha256=sources,input_sha256=inputs,
                  reference_native_manifest=str(reference/"manifest.json"),reference_native_manifest_sha256=sha(reference/"manifest.json"),
                  native_protocol=dict(preparation_duration=10.,cue_peak=4.,cue_full_physical_support=.2,cue_free_release=13000.,output_interval=5.),
                  solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
                  qualification=dict(two_cue_free_windows=[[1000.,2000.],[2000.,3000.]],window_tolerance=1e-4,
                                     holdout=[3000.,13000.],holdout_tolerance=1e-3,last_window_tolerance=1e-4,
                                     rhs_projected1000_tolerance=1e-4,minimum_contrast=.05,one_peak=True,one_halfcontrast_region=True,
                                     amount_relative_tolerance=1e-7,minimum_concentration=-5e-10),
                  unchanged_native_parameters=True,new_condition_grid="Pending actual basal edge feasibility; no challenge jobs in this stage")
    write(run/"manifest.json",manifest);write(run/"freeze.json",dict(generated_at=now(),manifest_sha256=sha(run/"manifest.json"),fresh_preparation_jobs=2,existing_preparations_to_independently_qualify=3))
    collect(run);print(json.dumps(read(run/"freeze.json"),indent=2))


def submit(run):
    assert not(run/"launch.json").exists();manifest=read(run/"manifest.json")
    assert sha(run/"manifest.json")==read(run/"freeze.json")["manifest_sha256"]
    for rel,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(run/rel)==digest
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:05:00","--array=0-1",
          "--job-name=arc-narrow-basals","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command]
    job=subprocess.check_output(args,text=True).strip();write(run/"launch.json",dict(generated_at=now(),array_job_id=job,sbatch_arguments=args,manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(read(run/"launch.json"),indent=2))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Native preparation requires SLURM")
    plan=read(run/"manifest.json");job=plan["jobs"][slot]
    for rel,digest in plan["source_sha256"].items():assert sha(run/"frozen_sources"/rel)==digest
    assert sha(run/job["input"])==job["input_sha256"]
    directory=run/"workers"/f"{slot:03d}";directory.mkdir();began=time.monotonic()
    status=dict(generated_at=now(),producer="running",status="preparation",complete=False,job=job,slurm_job_id=os.environ["SLURM_JOB_ID"])
    write(directory/"status.json",status)
    try:
        with np.load(run/job["input"]) as z:initial=z["state"].copy()
        size,n=job["size_parameter"],job["cells"];model=Closure(size,n)
        result=run_protocol(model,[CuePhase(10.,4.,job["preparation_axis_fraction"],.2,.1),CuePhase(13000.,output_interval=5.)],initial_state=initial)
        tt,yy=result["time"],result["states"];state=yy[-1].copy();scales=np.r_[np.full(2*n,2.),size]
        history=directory/"preparation_and_release.npz"
        np.savez_compressed(history,time=tt,states=yy,material_positions=model.material_positions,physical_amount=result["physical_amount"],resting_circumference=size)
        def window(a,b,tolerance):
            keep=(tt>=10+a)&(tt<=10+b);t,y=tt[keep],yy[keep];assert len(t)>2
            excursion=float(np.max(np.ptp(y,axis=0)/scales))
            drift=float(np.max(abs(y[-1]-y[-2])/scales)/(t[-1]-t[-2])*1000)
            return dict(cue_free_window=[a,b],samples=len(t),scaled_range=excursion,projected_last_saved_rate=drift,passed=max(excursion,drift)<=tolerance)
        windows=[window(1000,2000,1e-4),window(2000,3000,1e-4)]
        hold=window(3000,13000,1e-3);tail=window(12000,13000,1e-4)
        rhs=float(np.max(abs(model.rhs(state,0.))/scales)*1000)
        u=state[:n];contrast=float(np.ptp(u));mask=u>=u.min()+contrast/2
        peaks,_=find_peaks(np.tile(u,3),prominence=.05);peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
        regions=int(np.sum(mask&~np.roll(mask,1)))
        error=float(result["max_saved_amount_relative_error"]);minimum=float(yy[:,:-1].min())
        passed=all(w["passed"] for w in windows+[hold,tail]) and rhs<=1e-4 and contrast>=.05 and len(peaks)==regions==1 and error<=1e-7 and minimum>=-5e-10
        target=run/"basals"/f"spring_C{size:g}_{'qualified' if passed else 'unqualified'}.npz"
        np.savez_compressed(target,state=state,material_positions=model.material_positions,resting_circumference=size,physical_circumference=state[-1],fixed_scales=scales)
        status.update(generated_at=now(),producer="finished",status="qualified" if passed else "unqualified",complete=True,
                      checkpoint=str(target),checkpoint_sha256=sha(target),history_path=str(history),history_sha256=sha(history),
                      fixed_scales=scales.tolist(),two_windows=windows,holdout=hold,final_window=tail,rhs_projected1000=rhs,
                      features=dict(contrast=contrast,half_contrast_regions=regions,resolved_peaks=len(peaks),peak_cells=peaks.tolist(),physical_circumference=float(state[-1])),
                      amount_relative_error=error,minimum_concentration=minimum,elapsed_seconds=time.monotonic()-began)
    except Exception as error:
        status.update(generated_at=now(),producer="finished",status="error",complete=True,error=f"{type(error).__name__}: {error}");raise
    finally:write(directory/"status.json",status)
    print(json.dumps(status,indent=2))


def collect(run):
    plan=read(run/"manifest.json");records=[]
    for job in plan["jobs"]:
        path=run/"workers"/f'{job["slot"]:03d}'/"status.json"
        records.append(read(path) if path.exists() else dict(job=job,status="not_started",complete=False))
    write(run/"results/preparation_summary.json",dict(generated_at=now(),complete=all(r["complete"] for r in records),
          fresh_qualified=sum(r["status"]=="qualified" for r in records),fresh_planned=2,existing_chemical=plan["exact_existing_chemical_preparations"],records=records))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("action",choices=["freeze","submit","worker","collect"]);p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--slot",type=int)
    a=p.parse_args();run=a.run.resolve();globals()[a.action](run,a.slot) if a.action=="worker" else globals()[a.action](run)
