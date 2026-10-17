"""Bounded, frozen SLURM distance pilot from exact N800 chemical ring states."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import scipy
from scipy.integrate import BDF
from scipy.signal import find_peaks
import ring_distance_prototype_20261017 as ring
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/chemical_ring_distance_pilot_20261017"
BASAL = ROOT / "data/runs/ring_preparation_pilot_N800_20261017"
FIG = ROOT / "figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
PEAKS = {"wave_pinning":0.1,"goryachev":6.4,"holmes_model3":1.0}
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+".tmp");temp.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");temp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def weights_about(grid,center,left,right):
    arc=grid.arc(grid.x,center)
    return np.maximum(0,np.minimum(arc+grid.dx/2,right)-np.maximum(arc-grid.dx/2,left))


def original_core(grid,primary):
    peak=int(np.argmax(primary));origin=float(grid.x[peak]);arc=grid.arc(grid.x,origin)
    order=np.argsort(arc);x,u=arc[order],primary[order]
    level=float(u.min()+0.5*np.ptp(u));inside=np.flatnonzero(u>=level)
    assert np.all(np.diff(inside)==1)
    a,b=int(inside[0]),int(inside[-1]);assert a>0 and b<len(x)-1
    left=x[a-1]+(level-u[a-1])/(u[a]-u[a-1])*(x[a]-x[a-1])
    right=x[b]+(level-u[b])/(u[b+1]-u[b])*(x[b+1]-x[b])
    weights=weights_about(grid,origin,left,right)
    np.testing.assert_allclose(weights.sum(),right-left,atol=1e-13)
    return dict(peak_cell=peak,origin=origin,level=level,left=float(left),right=float(right),
                width=float(right-left),weights=weights)


def init(run):
    if run.exists():raise ValueError("Preserve existing output; use a fresh run")
    if "data/runs" not in str(run):raise ValueError("Use a fresh data/runs directory")
    for name in ("inputs","frozen_sources","workers","logs","results"):(run/name).mkdir(parents=True)
    paths=[Path(__file__),Path(__file__).with_name("ring_distance_prototype_20261017.py"),
           Path(__file__).with_name("length_dev_output_layout.py")]
    for folder in (ROOT/"code/core/models",ROOT/"code/core/src/polarity1d"):paths+=sorted(folder.glob("*.py"))
    hashes={}
    for path in paths:
        rel=path.relative_to(ROOT);dest=run/"frozen_sources"/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,dest);hashes[str(rel)]=sha(dest)
    manifest=read(BASAL/"manifest.json");audit=read(BASAL/"results/independent_audit.json");assert audit["passed"]
    grid=ring.RingGrid(40,800);models={};jobs=[];excluded=[]
    for base in manifest["jobs"]:
        key=base["model"];status=read(BASAL/"workers"/f'{base["slot"]:02d}'/"status.json")
        assert status["status"]=="qualifying_single_site"
        original=BASAL/status["checkpoint"];assert sha(original)==status["checkpoint_sha256"]
        path=run/"inputs"/f"{key}.npz";shutil.copy2(original,path);assert sha(path)==status["checkpoint_sha256"]
        with np.load(path) as z:state=z["state"].copy()
        core=original_core(grid,state[0]);reflect=(2*core["peak_cell"]-np.arange(grid.cells))%grid.cells
        scales=np.asarray([base["fixed_field_scales"][name] for name in base["state_names"]])
        models[key]=dict(checkpoint=str(path.relative_to(run)),checkpoint_sha256=sha(path),
            original_checkpoint=str(original.relative_to(ROOT)),state_names=base["state_names"],
            fixed_field_scales=base["fixed_field_scales"],diffusivities=base["diffusivities"],
            conserved_groups=base["conserved_groups"],incumbent_peak_cell=core["peak_cell"],
            fixed_original_core={k:v for k,v in core.items() if k!="weights"},
            initial_full_field_reflection_error=float(np.max(np.abs(state-state[:,reflect])/scales[:,None])))
        for distance,control in [(8,False),(12,False),(16,False),(20,False),(-12,True)]:
            offset=int(round(distance/grid.dx));center=(core["peak_cell"]+offset)%grid.cells
            cue=ring.ring_cue(grid,center,3,PEAKS[key]);overlap=bool(np.any((core["weights"]>0)&(cue>0)))
            if overlap:
                excluded.append(dict(model=key,signed_distance=distance,reason="Cue support overlaps fixed original half-contrast core"));continue
            jobs.append(dict(slot=len(jobs),model=key,signed_distance=distance,
                shortest_arc=abs(distance),cue_center_cell=int(center),cue_peak=PEAKS[key],
                role="reflected cue numerical control" if control else "distance case"))
    assert len(jobs)==14 and len(excluded)==1 and excluded[0]["model"]=="holmes_model3"
    result=dict(created_at=now(),experiment="First bounded chemical ring distance pilot",models=models,jobs=jobs,
        excluded=excluded,circumference=40,cells=800,cue_total_physical_support=3,
        cue_duration=1000.0,basal_duration=13000.0,basal_external_input=0.0,
        sampling="0.1 native unit through10;1 through100;5 thereafter, restarted at withdrawal",
        solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.0),
        limits=dict(allowed_negative_roundoff=5e-10,pool_relative_tolerance=1e-7,worker_soft_seconds=3300),
        tail=dict(windows=[[11000,12000],[12000,13000]],tolerance=1e-4,
                  interpretation="Two late-window full-field stationarity diagnostics only; not a new front-establishment or takeover gate."),
        readout_definition=dict(original_core="Frozen basal half-contrast interval; fractional finite-volume edge cells; raw primary concentration mean and integral.",
            cue_region="Fixed total-width3 geometric support, fractional edge cells; raw primary mean/integral even after withdrawal.",
            spatial="Complete native fields, every raw local maximum and prominence/valley bases; frozen basal half-contrast region count is descriptive only.",
            identity="A global maximum or one moving cap is not an identified new front; full spatial histories distinguish motion, coexistence, competition and merger.",
            controls="Reflected cues restart the same unmodified basal array; numerical symmetry controls, not biological replicates."),
        source_sha256=hashes,basal_manifest_sha256=sha(BASAL/"manifest.json"),
        basal_audit_sha256=sha(BASAL/"results/independent_audit.json"),
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,executable=PYTHON),
        scope="No retuning, duration optimization or exhaustive threshold search; native reactions/diffusion/pools unchanged.")
    write(run/"manifest.json",result);write(run/"status.json",dict(generated_at=now(),producer="not started",complete=False))
    return result


def launch(run):
    manifest=init(run);runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command="exec "+shlex.join(["env","PYTHONDONTWRITEBYTECODE=1","OPENBLAS_NUM_THREADS=1","OMP_NUM_THREADS=1","MKL_NUM_THREADS=1",PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=01:00:00",
          "--signal=B:USR1@60",f'--array=0-{len(manifest["jobs"])-1}',"--job-name=ring-distance",
          "--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap",command]
    job=subprocess.check_output(args,text=True).strip();write(run/"launch.json",dict(generated_at=now(),array_job_id=job,command=args,manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(dict(array_job_id=job,workers=len(manifest["jobs"]),run=str(run)),indent=2))


def sampling(end):
    result=np.r_[np.arange(0,10.00001,.1),np.arange(11,101,1),np.arange(105,end+0.01,5)]
    return np.unique(np.r_[result[result<=end],end])


def phase(model,grid,state,cue,end,manifest,deadline,initial_pools):
    species=len(model.state_names);wanted=sampling(end);times=[0.0];states=[state.copy()];index=1
    solver=BDF(lambda t,y:ring.native_rhs(model,grid,y.reshape(species,grid.cells),cue).ravel(),0,state.ravel(),end,
        rtol=manifest["solver"]["rtol"],atol=manifest["solver"]["atol"],max_step=manifest["solver"]["max_step"],
        jac_sparsity=ring.periodic_sparsity(species,grid.cells))
    status,message="complete","";worst_pool=0.;minimum=float(state.min())
    while solver.status=="running":
        if STOP or time.monotonic()>=deadline:status,message="timeout","SLURM advance signal or soft wall limit";break
        message=solver.step() or ""
        if solver.status=="failed":status="numerical_error";break
        accepted=solver.y.reshape(species,grid.cells);minimum=min(minimum,float(accepted.min()))
        worst_pool=max(worst_pool,float(np.max(np.abs(ring.pool_amounts(model,grid,accepted)/initial_pools-1))))
        dense=solver.dense_output()
        while index<len(wanted) and wanted[index]<=solver.t:
            times.append(float(wanted[index]));states.append(dense(wanted[index]).reshape(species,grid.cells).copy());index+=1
        if minimum < -manifest["limits"]["allowed_negative_roundoff"] or worst_pool>manifest["limits"]["pool_relative_tolerance"]:
            status,message="numerical_error","Concentration or conserved-pool check failed";break
    if times[-1]!=solver.t:times.append(float(solver.t));states.append(solver.y.reshape(species,grid.cells).copy())
    return dict(time=np.asarray(times),states=np.asarray(states),status=status,message=message,
        max_pool_relative_error=worst_pool,minimum_concentration=minimum,nfev=int(solver.nfev),njev=int(solver.njev),nlu=int(solver.nlu))


def readouts(grid,states,times,original,core,cue_center,model):
    oldw=core["weights"];cuew=weights_about(grid,float(grid.x[cue_center]),-1.5,1.5)
    primary=states[:,0];oldamount=primary@oldw;cueamount=primary@cuew
    angle=2*np.pi*(grid.x-core["origin"])/grid.circumference
    harmonic=(primary@np.exp(1j*angle))/primary.sum(axis=1)
    peak_frames=[];peak_cells=[];peak_values=[];prominences=[];left_bases=[];right_bases=[];left_values=[];right_values=[]
    counts=[];resolved=[];global_cells=[]
    for frame,values in enumerate(primary):
        maxima,properties=find_peaks(np.tile(values,3),prominence=0)
        keep=(maxima>=grid.cells)&(maxima<2*grid.cells);maxima=maxima[keep]-grid.cells
        left=properties["left_bases"][keep]%grid.cells;right=properties["right_bases"][keep]%grid.cells
        peak_frames.extend([frame]*len(maxima));peak_cells.extend(maxima);peak_values.extend(values[maxima]);prominences.extend(properties["prominences"][keep])
        left_bases.extend(left);right_bases.extend(right);left_values.extend(values[left]);right_values.extend(values[right])
        counts.append(len(ring.circular_components(values>=core["level"])))
        resolved.append(int(np.sum(properties["prominences"][keep]>=0.05)));global_cells.append(int(np.argmax(values)))
    pools=np.column_stack([states[:,list(group),:].sum(axis=(1,2))*grid.dx for group in model.conserved_groups])
    return dict(old_core_amount=oldamount,old_core_mean=oldamount/oldw.sum(),cue_region_amount=cueamount,
        cue_region_mean=cueamount/cuew.sum(),primary_total=primary.sum(axis=1)*grid.dx,pools=pools,
        harmonic_real=harmonic.real,harmonic_imag=harmonic.imag,harmonic_magnitude=np.abs(harmonic),
        frozen_basal_half_contrast_components=np.asarray(counts),raw_local_peak_count=np.asarray([np.sum(np.asarray(peak_frames)==i) for i in range(len(times))]),
        peaks_over_native_prominence_reference=np.asarray(resolved),global_max_cell=np.asarray(global_cells),
        peak_frame=np.asarray(peak_frames,dtype=int),peak_cell=np.asarray(peak_cells,dtype=int),peak_value=np.asarray(peak_values),
        peak_prominence=np.asarray(prominences),peak_left_valley_cell=np.asarray(left_bases,dtype=int),peak_right_valley_cell=np.asarray(right_bases,dtype=int),
        peak_left_valley_value=np.asarray(left_values),peak_right_valley_value=np.asarray(right_values))


def tail_diagnostic(times,states,scales):
    records=[]
    for a,b in ((11000,12000),(12000,13000)):
        keep=(times>=a)&(times<=b)
        if keep.sum()<2:records.append(dict(window=[a,b],available=False));continue
        t,y=times[keep],states[keep];variation=np.ptp(y,axis=0).max(axis=1)/scales
        drift=np.abs(y[-1]-y[-2]).max(axis=1)/(t[-1]-t[-2])*1000/scales
        records.append(dict(window=[a,b],available=True,steady=bool(np.max(np.r_[variation,drift])<=1e-4),relative_range=variation.tolist(),projected_relative_drift=drift.tolist()))
    return dict(windows=records,both_late_windows_steady=all(row.get("steady",False) for row in records))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("All integration requires SLURM")
    signal.signal(signal.SIGUSR1,lambda s,f:globals().__setitem__("STOP",True));signal.signal(signal.SIGTERM,lambda s,f:globals().__setitem__("STOP",True))
    manifest=read(run/"manifest.json");job=manifest["jobs"][slot];base=manifest["models"][job["model"]]
    for path,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/path)==digest,path
    assert sha(run/base["checkpoint"])==base["checkpoint_sha256"]
    with np.load(run/base["checkpoint"]) as z:original=z["state"].copy()
    model=ring.native_model(job["model"]);grid=ring.RingGrid(manifest["circumference"],manifest["cells"])
    core=original_core(grid,original[0]);scales=np.asarray([base["fixed_field_scales"][name] for name in base["state_names"]])
    initial_pools=ring.pool_amounts(model,grid,original);cue=ring.ring_cue(grid,job["cue_center_cell"],3,job["cue_peak"])
    directory=run/"workers"/f"{slot:02d}";directory.mkdir();record=dict(generated_at=now(),producer="running",complete=False,job=job,status="starting",phases=[])
    status_path=directory/"status.json";write(status_path,record);deadline=time.monotonic()+manifest["limits"]["worker_soft_seconds"]
    state=original.copy()
    try:
        for name,end,signal_values in (("cue",1000.,cue),("recovery",13000.,np.zeros(grid.cells))):
            record.update(generated_at=now(),status=name);write(status_path,record)
            result=phase(model,grid,state,signal_values,end,manifest,deadline,initial_pools)
            metrics=readouts(grid,result["states"],result["time"],original,core,job["cue_center_cell"],model)
            path=directory/f"{name}.npz";np.savez_compressed(path,time=result["time"],states=result["states"],x=grid.x,
                stimulus=signal_values,original_core_weights=core["weights"],**metrics)
            info={key:value for key,value in result.items() if key not in ("time","states")}
            info.update(phase=name,path=str(path.relative_to(run)),sha256=sha(path),samples=len(result["time"]))
            record["phases"].append(info);state=result["states"][-1].copy()
            if name=="recovery":
                record["tail"]=tail_diagnostic(result["time"],result["states"],scales)
                record["final_alignment_to_original"]=ring.translation_compare(state,original,scales)
                record["final_raw_old_core_mean"]=float(metrics["old_core_mean"][-1]);record["final_raw_cue_mean"]=float(metrics["cue_region_mean"][-1])
            if result["status"]!="complete":record["status"]=result["status"];break
        else:record["status"]="complete_steady_tail" if record["tail"]["both_late_windows_steady"] else "complete_unsettled_tail"
        record.update(generated_at=now(),producer="finished",complete=True)
    except Exception as error:
        record.update(generated_at=now(),producer="finished",complete=True,status="error",error=f"{type(error).__name__}: {error}");write(status_path,record);raise
    write(status_path,record);print(json.dumps(dict(job=job,status=record["status"]),indent=2))


def status(run):
    manifest=read(run/"manifest.json");workers=[]
    for job in manifest["jobs"]:
        path=run/"workers"/f'{job["slot"]:02d}'/"status.json"
        workers.append(read(path) if path.exists() else dict(job=job,complete=False,status="not_started"))
    result=dict(generated_at=now(),complete=all(row["complete"] for row in workers),workers=workers)
    result["producer"]="finished" if result["complete"] else "running";write(run/"status.json",result)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=["launch","worker","status"])
    parser.add_argument("--run",type=Path,default=DEFAULT);parser.add_argument("--slot",type=int)
    args=parser.parse_args();run=args.run.resolve()
    if args.action=="launch":launch(run)
    elif args.action=="worker":worker(run,args.slot)
    else:
        result=status(run);print(json.dumps(dict(complete=result["complete"],workers=[dict(slot=row["job"]["slot"],model=row["job"]["model"],distance=row["job"]["signed_distance"],status=row["status"]) for row in result["workers"]]),indent=2))


if __name__=="__main__":main()
