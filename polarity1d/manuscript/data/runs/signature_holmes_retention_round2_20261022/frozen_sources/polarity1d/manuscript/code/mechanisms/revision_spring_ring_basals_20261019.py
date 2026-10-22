"""Shared exact Spring-ring preparations for revised Figures 3 and 4.

C0=2 preserves the independently validated exact basal. C0=1/3 use fresh
uniform initial states, the same native physical cue and a full release.
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
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from scipy.signal import find_peaks
from spring_ring_model_20261019 import Closure, CuePhase, run_protocol

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"data/runs/revision_ring_spring_20261019"
PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"
RETAINED=ROOT/"data/runs/spring_ring_preparation_pilot_20261017/workers/qualified_basal.npz"


def now():return datetime.now(timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,record):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n");temporary.replace(path)
def read(path):return json.loads(Path(path).read_text())


def features(state,n):
    u=state[:n];contrast=float(np.ptp(u));mask=u>=u.min()+.5*contrast
    regions=int(np.sum(mask&~np.roll(mask,1))) if contrast>1e-12 else 0
    peaks,_=find_peaks(np.tile(u,3),prominence=.05)
    peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
    return dict(contrast=contrast,half_contrast_regions=regions,resolved_peaks=len(peaks),
        peak_cells=peaks.tolist(),peak=float(u.max()),background=float(u.min()),physical_circumference=float(state[-1]))


def launch(run):
    if run.exists():raise FileExistsError("Preserve completed output; fresh run required")
    for name in ("frozen_sources","preparations","logs","results","inputs"):(run/name).mkdir(parents=True)
    paths=[Path(__file__),Path(__file__).with_name("spring_ring_model_20261019.py"),
        Path(__file__).with_name("spring_ring_closure_20261017.py"),ROOT/"code/core/models/spring.py",
        ROOT/"code/core/models/__init__.py"]
    hashes={}
    for p in paths:
        relative=p.relative_to(ROOT);target=run/"frozen_sources"/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);hashes[str(relative)]=sha(target)
    shutil.copy2(RETAINED,run/"inputs/retained_C2.npz")
    jobs=[dict(slot=i,resting_circumference=C,cells=int(100*C),reused=(C==2.)) for i,C in enumerate((1.,2.,3.))]
    manifest=dict(generated_at=now(),model="Spring mean-force global-extension circular analogue",jobs=jobs,
        initial="Fresh native uniform concentrations; no reflected/remapped checkpoint at new sizes",
        preparation_cue_duration=10.,basal_release=13000.,cue_peak=4.,cue_full_physical_support=.2,
        source_sha256=hashes,retained_input_sha256=sha(run/"inputs/retained_C2.npz"),
        retained_original=str(RETAINED),retained_original_audit=str(ROOT/"data/runs/spring_ring_preparation_pilot_20261017/results/audit.json"),
        native_time_and_length_uncalibrated=True,solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
        qualification=dict(tail_window=1000.,scaled_range=1e-4,rhs_projected_to_1000=1e-4,minimum_contrast=.05,
            one_half_contrast_region=True,one_resolved_peak=True,physical_amount_relative_tolerance=1e-7))
    write(run/"manifest.json",manifest)
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job=subprocess.check_output(["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:30:00","--array=0-2",
        "--job-name=revision-Spring-basals","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command],text=True).strip()
    write(run/"launch.json",dict(generated_at=now(),array_job_id=job,manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(dict(array_job_id=job,run=str(run),jobs=jobs),indent=2))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("All native integration uses SLURM")
    manifest=read(run/"manifest.json");job=manifest["jobs"][slot]
    for p,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/p)==digest
    assert sha(run/"inputs/retained_C2.npz")==manifest["retained_input_sha256"]
    C0=job["resting_circumference"];n=job["cells"];model=Closure(C0,n)
    directory=run/"preparations"/f"C{C0:g}";directory.mkdir()
    record=dict(generated_at=now(),producer="running",status="preparing",job=job,slurm_job_id=os.environ["SLURM_JOB_ID"])
    write(directory/"status.json",record);began=time.monotonic()
    if job["reused"]:
        with np.load(run/"inputs/retained_C2.npz") as z:state=z["state"].copy()
        source_audit=read(Path(manifest["retained_original_audit"]));assert source_audit["passed"]
        amount_error=0.;tail_range=source_audit["holdout"]["excursion"]
        record["retained_qualification_audit"]=manifest["retained_original_audit"]
        record["retained_qualification_audit_sha256"]=sha(Path(manifest["retained_original_audit"]))
    else:
        result=run_protocol(model,[CuePhase(10.,4.,model.material_positions[0],.2,.1),CuePhase(13000.,output_interval=5.)])
        state=result["states"][-1]
        archive=directory/"preparation_and_release.npz"
        np.savez_compressed(archive,time=result["time"],states=result["states"],material_positions=model.material_positions,
            physical_amount=result["physical_amount"],resting_circumference=C0)
        record["history_path"]=str(archive);record["history_sha256"]=sha(archive)
        amount_error=result["max_saved_amount_relative_error"]
        scales=np.r_[np.full(2*n,2.),C0]
        tail=result["time"]>=result["time"][-1]-1000
        tail_range=float(np.max(np.ptp(result["states"][tail],axis=0)/scales))
    scales=np.r_[np.full(2*n,2.),C0]
    rhs=float(np.max(abs(model.rhs(state,0.))/scales)*1000.)
    f=features(state,n)
    passed=(tail_range<=1e-4 and rhs<=1e-4 and f["contrast"]>=.05 and f["half_contrast_regions"]==1 and f["resolved_peaks"]==1 and amount_error<=1e-7)
    checkpoint=directory/("qualified_basal.npz" if passed else "unqualified_final.npz")
    np.savez_compressed(checkpoint,state=state,material_positions=model.material_positions,resting_circumference=C0,
        fixed_scales=scales,physical_circumference=state[-1])
    record.update(generated_at=now(),producer="finished",status="qualified" if passed else "unqualified",complete=True,
        checkpoint=str(checkpoint),checkpoint_sha256=sha(checkpoint),fixed_scales=scales.tolist(),
        tail_1000_scaled_range=float(tail_range),rhs_projected_1000=rhs,physical_amount=model.amount(state),
        amount_relative_error=float(amount_error),minimum_concentration=float(state[:-1].min()),features=f,
        elapsed_seconds=time.monotonic()-began)
    write(directory/"status.json",record);print(json.dumps(record,indent=2))


def collect(run):
    manifest=read(run/"manifest.json");rows=[]
    for job in manifest["jobs"]:
        path=run/"preparations"/f'C{job["resting_circumference"]:g}'/"status.json"
        rows.append(read(path) if path.exists() else dict(job=job,status="not_started",complete=False))
    complete=all(r.get("complete",False) for r in rows)
    write(run/"results/preparation_summary.json",dict(generated_at=now(),complete=complete,
        qualified=sum(r["status"]=="qualified" for r in rows),planned=len(rows),records=rows))
    print(json.dumps(dict(complete=complete,statuses=[(r["job"]["resting_circumference"],r["status"]) for r in rows]),indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=["launch","worker","collect"])
    parser.add_argument("--run",type=Path,default=DEFAULT);parser.add_argument("--slot",type=int)
    args=parser.parse_args()
    if args.action=="launch":launch(args.run.resolve())
    elif args.action=="collect":collect(args.run.resolve())
    else:worker(args.run.resolve(),args.slot)
