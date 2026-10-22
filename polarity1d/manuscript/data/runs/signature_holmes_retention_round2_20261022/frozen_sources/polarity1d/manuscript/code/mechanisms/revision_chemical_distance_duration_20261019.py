"""Targeted C40 distance-duration extension motivated by retained trajectories.

Five placements/model x four fresh durations, with exact retained1000 references.
No monotonicity assumption, binary search or invented threshold.
"""
from pathlib import Path
import argparse
import json
import os
import shlex
import shutil
import subprocess
import time
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import chemical_ring_functional_v2_20261017 as campaign
import ring_distance_prototype_20261017 as ring
import chemical_ring_distance_pilot_20261017 as distance
import chemical_ring_functional_collect_20261017 as oldcollect

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/runs/chemical_ring_functional_v2_20261017"
DEFAULT=ROOT/"data/runs/revision_chemical_distance_duration_20261019"
PYTHON=campaign.PYTHON
read,write,sha,now=campaign.read,campaign.write,campaign.sha,campaign.now
FRACTIONS={"wave_pinning":[.20,.25,.30,.40,.50],"goryachev":[.10,.20,.30,.40,.50],"holmes_model3":[.25,.30,.35,.40,.50]}


def samples(end):
    return np.unique(np.r_[0.,np.geomspace(1e-4,min(.1,end),25),np.arange(.1,min(10.,end),.1),
        np.arange(10.,min(100.,end),1.),np.arange(100.,min(1000.,end),5.),np.arange(1000.,end,50.),end])


def launch(run):
    if run.exists():raise FileExistsError("Fresh output directory required")
    for name in ("inputs","frozen_sources","workers","logs","results"):(run/name).mkdir(parents=True)
    paths=[Path(__file__)]+[Path(__file__).with_name(name) for name in ("chemical_ring_functional_v2_20261017.py",
        "chemical_ring_distance_pilot_20261017.py","ring_distance_prototype_20261017.py","ring_preparation_pilot_20261017.py",
        "length_dev_output_layout.py","chemical_ring_functional_collect_20261017.py")]
    for folder in (ROOT/"code/core/models",ROOT/"code/core/src/polarity1d"):paths+=list(folder.glob("*.py"))
    frozen={}
    for path in paths:
        rel=path.relative_to(ROOT);target=run/"frozen_sources"/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);frozen[str(rel)]=sha(target)
    matrix=pd.read_csv(SOURCE/"results/matrix.csv");oldmanifest=read(SOURCE/"manifest.json");jobs=[];references=[];inputs={}
    for model,fractions in FRACTIONS.items():
        base=next(j for j in oldmanifest["jobs"] if j["model"]==model and j["circumference"]==40)
        for fraction in fractions:
            row=matrix[(matrix.model==model)&(matrix.circumference==40)&np.isclose(matrix.fraction,fraction)].iloc[0]
            assert row.status=="complete" and row.short_edge_to_cue_gap>0
            oldstatus=read(Path(row.status_path));source=Path(oldstatus["checkpoint"]);assert sha(source)==row.checkpoint_sha256
            target=run/"inputs"/(model+"_C40.npz")
            if not target.exists():shutil.copy2(source,target)
            inputs[str(target.relative_to(run))]=sha(target)
            common=dict(model=model,circumference=40,cells=int(row.cells),fraction=fraction,short_route=float(row.short_route),
                cue_center_cell=int(row.cue_center_cell),cue_peak=campaign.previous.PEAKS[model],cue_full_width=3.,recovery=13000.,
                fixed_scales=base["fixed_scales"],input=str(target.relative_to(run)),input_sha256=sha(target),
                source_cue=str(row.cue_path),source_cue_sha256=sha(row.cue_path),source_recovery=str(row.recovery_path),source_recovery_sha256=sha(row.recovery_path))
            references.append(dict(**common,duration=1000.,status="exact_retained_reference"))
            for duration in (1.,10.,100.,300.):jobs.append(dict(**common,slot=len(jobs),duration=duration))
    manifest=dict(generated_at=now(),experiment="Trajectory-motivated C40 chemical distance-duration pilot",jobs=jobs,reused_references=references,
        fractions=FRACTIONS,fresh_durations=[1.,10.,100.,300.],reused_duration=1000.,recovery=13000.,source_sha256=frozen,input_sha256=inputs,
        source_matrix=str(SOURCE/"results/matrix.csv"),source_matrix_sha256=sha(SOURCE/"results/matrix.csv"),
        solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
        rationale="C40 retained histories show195-460 Gory half-loss, nonmonotone245-830 WP half-loss and no Holmes half-loss by1000; use new shorter cues to inspect retained angle and regional activity separately",
        interpretation="Sixty fresh fixed-size native branches; preserve successful islands and shortest tested passing duration; no unique minimum or general distance law")
    write(run/"manifest.json",manifest)
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job=subprocess.check_output(["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:30:00",f"--array=0-{len(jobs)-1}",
        "--job-name=chemical-distance-duration","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command],text=True).strip()
    write(run/"launch.json",dict(generated_at=now(),array_job_id=job,manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(dict(array_job_id=job,native_jobs=len(jobs),retained_references=len(references),run=str(run)),indent=2))


def bracket(t,values,threshold):
    hit=np.flatnonzero(values<=threshold)
    if not len(hit):return dict(status="right_censored",lower=float(t[-1]),upper=None)
    i=int(hit[0]);return dict(status="observed_saved_sample_crossing",lower=float(t[max(0,i-1)]),upper=float(t[i]))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Native integration uses SLURM")
    manifest=read(run/"manifest.json");job=manifest["jobs"][slot]
    for rel,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(run/rel)==digest
    directory=run/"workers"/f"{slot:03d}";directory.mkdir()
    record=dict(generated_at=now(),producer="running",status="challenge",complete=False,job=job,slurm_job_id=os.environ["SLURM_JOB_ID"],phases=[])
    write(directory/"status.json",record);began=time.monotonic()
    model=ring.native_model(job["model"]);grid=ring.RingGrid(40,job["cells"])
    with np.load(run/job["input"]) as z:basal=z["state"].copy()
    state=basal.copy();pools=ring.pool_amounts(model,grid,basal);core=distance.original_core(grid,basal[0]);parts=[]
    cue=ring.ring_cue(grid,job["cue_center_cell"],3.,job["cue_peak"])
    for name,duration,signal,offset in (("challenge",job["duration"],cue,0.),("recovery",13000.,np.zeros(grid.cells),job["duration"])):
        record.update(status=name,generated_at=now());write(directory/"status.json",record)
        t,y,stats=campaign.integrate(model,grid,state,signal,duration,samples(duration),pools,began+1700.)
        assert t[-1]==duration
        archive=directory/(name+".npz");np.savez_compressed(archive,time=t,states=y,x=grid.x,stimulus=signal)
        record["phases"].append(dict(kind=name,path=str(archive.relative_to(run)),sha256=sha(archive),samples=len(t),**stats))
        parts.append((t+offset,y));state=y[-1].copy()
    np.testing.assert_array_equal(parts[0][1][-1],parts[1][1][0])
    t=np.r_[parts[0][0],parts[1][0][1:]];y=np.concatenate((parts[0][1],parts[1][1][1:]));p=y[:,0]
    bg=float(basal[0].min());oldraw=p@core["weights"];background=bg*core["weights"].sum();E0=oldraw[0]-background;excess=(oldraw-background)/E0
    patch=distance.weights_about(grid,float(grid.x[job["cue_center_cell"]]),-1.5,1.5)
    patch_excess=((p@patch)-bg*patch.sum())/E0
    theta=2*np.pi*(grid.x-core["origin"])/40;M0=grid.dx*p.sum(axis=1);M1=grid.dx*(p@np.exp(1j*theta));M2=grid.dx*(p@np.exp(2j*theta));P=M1/M1[0];R=abs(P)
    angles=np.rad2deg(np.angle(P));angles=np.where(np.isclose(abs(angles),180.,atol=1e-7,rtol=0),180.,angles)
    error=np.rad2deg(np.angle(np.exp(1j*np.deg2rad(angles-360*job["fraction"]))))
    idxoff=len(parts[0][0])-1;tail=t>=t[-1]-1000;scales=np.asarray(job["fixed_scales"])[:,None]
    tail_range=float(np.max(np.ptp(y[tail],axis=0)/scales));rhs=float(np.max(abs(ring.native_rhs(model,grid,y[-1],np.zeros(grid.cells)))/scales)*1000)
    spatial=oldcollect.identity(p[-1],core,job["cue_center_cell"],campaign.spatial(grid,p[-1],core,job["cue_center_cell"]))
    output=dict(**job,status="complete",cue_duration=job["duration"],physical_distance=job["short_route"],angular_separation_deg=360*job["fraction"],
        old_core_excess_at_withdrawal=float(excess[idxoff]),old_core_excess_final=float(excess[-1]),old_core_raw_final=float(oldraw[-1]/oldraw[0]),
        old_core_excess_minimum=float(excess.min()),old_core_excess_minimum_time=float(t[excess.argmin()]),old_core_excess_rebound=float(excess[-1]-excess[idxoff]),
        cue_patch_excess_final=float(patch_excess[-1]),R_at_withdrawal=float(R[idxoff]),R_final=float(R[-1]),q1_final=float(abs(M1[-1])/M0[-1]),q2_final=float(abs(M2[-1])/M0[-1]),
        angle_at_withdrawal=float(angles[idxoff]) if R[idxoff]>=.1 else None,angle_final=float(angles[-1]) if R[-1]>=.1 else None,
        cue_angle_error_final=float(error[-1]) if R[-1]>=.1 else None,cue_aligned_final_15deg=bool(abs(error[-1])<=15 and R[-1]>=.1),
        old_aligned_final_15deg=bool(abs(angles[-1])<=15 and R[-1]>=.1),resolved_peaks_final=spatial["resolved_peak_count"],outcome=spatial["location_association"],
        tail_1000_scaled_range=tail_range,rhs_projected_1000=rhs,tail_stationary=bool(tail_range<=1e-4 and rhs<=1e-4),
        maximum_pool_relative_error=max(phase["max_pool_error"] for phase in record["phases"]))
    for threshold in (.9,.5,.1):
        tag=str(threshold).replace(".","p");event=bracket(t,excess,threshold);during=bracket(t[:idxoff+1],excess[:idxoff+1],threshold)
        output.update({f"loss_{tag}_{k}":v for k,v in event.items()});output.update({f"cue_on_loss_{tag}_{k}":v for k,v in during.items()})
    np.savez_compressed(directory/"readout_traces.npz",time=t,old_core_raw_amount=oldraw,old_core_excess=excess,cue_patch_excess=patch_excess,
        M0=M0,M1_real=M1.real,M1_imag=M1.imag,P_real=P.real,P_imag=P.imag,R=R,q1=abs(M1)/M0,q2=abs(M2)/M0,angle_deg_unmasked=angles)
    write(directory/"readouts.json",output)
    record.update(generated_at=now(),producer="finished",status="complete",complete=True,readouts_sha256=sha(directory/"readouts.json"),
        traces_sha256=sha(directory/"readout_traces.npz"),elapsed_seconds=time.monotonic()-began)
    write(directory/"status.json",record);print(json.dumps(output,indent=2))


def collect(run):
    manifest=read(run/"manifest.json");rows=[];statuses=[]
    for job in manifest["jobs"]:
        directory=run/"workers"/f'{job["slot"]:03d}';path=directory/"status.json"
        record=read(path) if path.exists() else dict(status="not_started",complete=False);statuses.append(record)
        if record.get("complete") and record["status"]=="complete":
            assert sha(directory/"readouts.json")==record["readouts_sha256"]
            row=read(directory/"readouts.json")
            for phase in record["phases"]:
                prefix="actual_cue" if phase["kind"]=="challenge" else "actual_recovery"
                row[prefix+"_path"]=str(run/phase["path"])
                row[prefix+"_sha256"]=phase["sha256"]
            row["lineage_reference_1000_cue_path"]=row["source_cue"]
            row["lineage_reference_1000_recovery_path"]=row["source_recovery"]
            rows.append(row)
        else:rows.append(dict(**job,status=record["status"]))
    pd.DataFrame(rows).to_csv(run/"results/distance_duration_readouts.csv",index=False)
    write(run/"results/collection.json",dict(generated_at=now(),complete=all(r.get("complete") and r["status"]=="complete" for r in statuses),
        native_complete=sum(r.get("complete",False) for r in statuses),planned=len(statuses),reused_references=len(manifest["reused_references"]),
        states={key:sum(r["status"]==key for r in statuses) for key in set(r["status"] for r in statuses)}))
    print(json.dumps(read(run/"results/collection.json"),indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("action",choices=["launch","worker","collect"]);p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--slot",type=int);a=p.parse_args()
    if a.action=="launch":launch(a.run.resolve())
    elif a.action=="collect":collect(a.run.resolve())
    else:worker(a.run.resolve(),a.slot)
