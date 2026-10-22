"""Bounded 500/2000-duration checks around retained 1000-unit ring challenges."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
import chemical_ring_functional_v2_20261017 as campaign
import ring_distance_prototype_20261017 as ring
import chemical_ring_distance_pilot_20261017 as distance
import chemical_ring_functional_collect_20261017 as oldcollect

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"data/runs/vision_ring_duration_20261017"
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/vision_ring_20261017"
SOURCE=ROOT/"data/runs/chemical_ring_functional_v2_20261017"
PYTHON=campaign.PYTHON
CASES=[("wave_pinning",40,.4),("goryachev",40,.4),("holmes_model3",40,.4),
       ("holmes_model3",40,.3),("holmes_model3",30,.4),("goryachev",40,.5),
       ("goryachev",80,.5),("holmes_model3",40,.5),("wave_pinning",80,.5)]
read,write,sha,now=campaign.read,campaign.write,campaign.sha,campaign.now


def initialize(run):
    if run.exists():raise ValueError("Fresh output required")
    for folder in ("inputs","frozen_sources","workers","logs","results"):(run/folder).mkdir(parents=True)
    paths=[Path(__file__)]+[Path(__file__).with_name(name) for name in (
        "chemical_ring_functional_v2_20261017.py","chemical_ring_distance_pilot_20261017.py",
        "ring_distance_prototype_20261017.py","ring_preparation_pilot_20261017.py","length_dev_output_layout.py",
        "chemical_ring_functional_collect_20261017.py")]
    for folder in (ROOT/"code/core/models",ROOT/"code/core/src/polarity1d"):paths+=list(folder.glob("*.py"))
    frozen={}
    for p in paths:
        rel=p.relative_to(ROOT);dest=run/"frozen_sources"/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,dest);frozen[str(rel)]=sha(dest)
    matrix=list(csv.DictReader((SOURCE/"results/matrix.csv").open()));jobs=[];inputs={}
    for model,C,fraction in CASES:
        row=next(r for r in matrix if r["model"]==model and float(r["circumference"])==C and float(r["fraction"])==fraction)
        assert row["status"]=="complete"
        s=read(Path(row["status_path"]));original=Path(s["checkpoint"])
        assert sha(original)==row["checkpoint_sha256"]
        dest=run/"inputs"/f"{model}_C{C}.npz"
        if not dest.exists():shutil.copy2(original,dest)
        inputs[str(dest.relative_to(run))]=sha(dest)
        for duration in (500.,2000.):
            jobs.append(dict(slot=len(jobs),model=model,circumference=C,cells=int(float(row["cells"])),
                fraction=fraction,short_route=float(row["short_route"]),cue_center_cell=int(float(row["cue_center_cell"])),
                cue_peak=campaign.previous.PEAKS[model],cue_full_width=3.,duration=duration,recovery=13000.,
                input=str(dest.relative_to(run)),input_sha256=sha(dest),
                source_cue=row["cue_path"],source_cue_sha256=sha(row["cue_path"]),
                source_recovery=row["recovery_path"],source_recovery_sha256=sha(row["recovery_path"]),
                source_worker_status=row["status_path"],source_worker_status_sha256=sha(row["status_path"])))
    m=dict(created_at=now(),experiment="Bounded finite ring cue-duration/readout sensitivity",
        jobs=jobs,durations=[500.,1000.,2000.],reused_duration=1000.,recovery=13000.,
        recovery_readouts=[3000.,8000.,13000.],solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
        frozen_source_sha256=frozen,input_sha256=inputs,plan_sha256=sha(TOPIC/"PLAN.md"),
        interpretation="Nine coverage anchors, not a duration threshold search. Same exact basal/cue; native units not cross-model calibration.")
    write(run/"manifest.json",m);write(run/"provenance.json",dict(manifest_sha256=sha(run/"manifest.json")))
    return m


def verify(run):
    m=read(run/"manifest.json");assert sha(run/"manifest.json")==read(run/"provenance.json")["manifest_sha256"]
    for k,d in m["frozen_source_sha256"].items():assert sha(run/"frozen_sources"/k)==d,k
    for k,d in m["input_sha256"].items():assert sha(run/k)==d
    return m


def load_phase(path):
    with np.load(path) as z:
        result={"time":z["time"].copy(),"x":z["x"].copy(),"stimulus":z["stimulus"].copy()}
        if "states" in z.files:
            result.update(primary=z["states"][:,0].copy(),states=z["states"].copy())
        else:
            result.update(primary=z["primary"].copy(),snapshot_times=z["snapshot_times"].copy(),snapshot_states=z["snapshot_states"].copy())
        return result


def exact_state(phase,t):
    key="time" if "states" in phase else "snapshot_times"
    found=np.flatnonzero(np.isclose(phase[key],t,rtol=0,atol=1e-10))
    if not len(found):return None
    return phase["states" if "states" in phase else "snapshot_states"][found[0]].copy()


def store(path,grid,times,states,stimulus,snapshot_times):
    chosen=[]
    for t in snapshot_times:
        idx=np.flatnonzero(np.isclose(times,t,rtol=0,atol=1e-10))
        if len(idx):chosen.append(int(idx[0]))
    chosen=np.unique(chosen)
    np.savez_compressed(path,time=times,primary=states[:,0],x=grid.x,stimulus=stimulus,
        snapshot_times=times[chosen],snapshot_states=states[chosen])


def worker(run,slot):
    assert os.environ.get("SLURM_JOB_ID"),"All simulation requires SLURM"
    m=verify(run);job=m["jobs"][slot];out=run/"workers"/f"{slot:02d}";out.mkdir()
    status=dict(created_at=now(),producer="running",complete=False,job=job,phases=[])
    write(out/"status.json",status);began=time.monotonic()
    try:
        with np.load(run/job["input"]) as z:basal=z["state"].copy()
        source=Path(job["source_cue"]);assert sha(source)==job["source_cue_sha256"]
        previous=load_phase(source);np.testing.assert_array_equal(exact_state(previous,0.),basal)
        model=ring.native_model(job["model"]);grid=ring.RingGrid(job["circumference"],job["cells"])
        cue=ring.ring_cue(grid,job["cue_center_cell"],3.,job["cue_peak"]);zero=np.zeros(grid.cells)
        np.testing.assert_array_equal(previous["stimulus"],cue)
        pools=ring.pool_amounts(model,grid,basal);core=distance.original_core(grid,basal[0])
        if job["duration"]==500 and exact_state(previous,500.) is not None:
            state=exact_state(previous,500.);keep=previous["time"]<=500
            status["cue_history"]=dict(path=str(source),sha256=sha(source),source_stop=500.,
                kind="exact retained cue prefix",samples=int(keep.sum()))
            cue_times=previous["time"][keep];cue_primary=previous["primary"][keep]
        else:
            start=1000. if job["duration"]==2000 else 0.
            state=exact_state(previous,start);assert state is not None
            end=job["duration"]-start
            samples=np.unique(np.r_[0.,np.arange(.1,10.0001,.1),np.arange(11,min(100,end)+.1),
                np.arange(105,end+.1,5),end])
            times,states,stats=campaign.integrate(model,grid,state,cue,end,samples,pools,began+3300)
            path=out/"cue_new.npz";store(path,grid,times,states,cue,[0.,1.,10.,500.,end])
            status["phases"].append(dict(phase="cue_new",path=str(path),sha256=sha(path),start_from_onset=start,**stats))
            if start==1000:
                cue_times=np.r_[previous["time"],times[1:]+1000]
                cue_primary=np.concatenate([previous["primary"],states[1:,0]])
                status["cue_history"]=dict(kind="exact retained 1000 prefix plus native continuation",prefix_path=str(source),
                    prefix_sha256=sha(source),extension_path=str(path),extension_sha256=sha(path))
                np.testing.assert_array_equal(states[0],exact_state(previous,1000.))
            else:
                cue_times,cue_primary=times,states[:,0]
                errors=[]
                for t in (0.,1.,10.):
                    oldstate=exact_state(previous,t);i=int(np.argmin(abs(times-t)))
                    errors.append(float(np.max(abs(states[i]-oldstate))))
                status["replayed_prefix_native_maxabs_at_0_1_10"]=errors
                assert max(errors)<1e-6
                prefix_errors=[]
                for i,t in enumerate(times):
                    idx=np.flatnonzero(np.isclose(previous["time"],t,rtol=0,atol=1e-10))
                    if len(idx):prefix_errors.append(float(np.max(abs(states[i,0]-previous["primary"][idx[0]]))))
                assert len(prefix_errors)==len(times)
                prefix_error=max(prefix_errors);prefix_scale=max(1.,float(basal[0].max()))
                assert prefix_error<1e-6 and prefix_error/prefix_scale<1e-7
                status["replayed_primary_prefix_check"]=dict(samples=len(prefix_errors),
                    end_time=float(times[-1]),max_absolute_error=prefix_error,
                    fixed_scale=prefix_scale,max_scaled_error=prefix_error/prefix_scale)
                status["cue_history"]=dict(kind="native 500 cue from exact basal; missing full retained state at500",
                    path=str(path),sha256=sha(path),source_primary_comparison_maxabs=max(errors))
            state=states[-1].copy()
        times,states,stats=campaign.integrate(model,grid,state,zero,13000.,np.r_[0.,np.arange(50,13001,50)],pools,began+3300)
        path=out/"recovery.npz";store(path,grid,times,states,zero,[0.,500.,1000.,3000.,8000.,13000.])
        status["phases"].append(dict(phase="recovery",path=str(path),sha256=sha(path),**stats))
        readouts=[]
        for t in (3000.,8000.,13000.):
            i=int(np.flatnonzero(times==t)[0]);spatial=campaign.spatial(grid,states[i,0],core,job["cue_center_cell"])
            spatial=oldcollect.identity(states[i,0],core,job["cue_center_cell"],spatial)
            theta=2*np.pi*(grid.x-core["origin"])/grid.circumference
            initial_moment=grid.dx*(basal[0]@np.exp(1j*theta));moment=grid.dx*(states[i,0]@np.exp(1j*theta))/initial_moment
            readouts.append(dict(cue_free_time=t,time_from_onset=job["duration"]+t,
                original_core_amount=float(states[i,0]@core["weights"]),
                original_core_amount_fraction=float(states[i,0]@core["weights"]/(basal[0]@core["weights"])),
                Px=float(moment.real),Py=float(moment.imag),R=float(abs(moment)),
                angle_deg=float(np.angle(moment,deg=True)) if abs(moment)>=.1 else None,
                spatial=spatial))
        combined_t=np.r_[cue_times,times[1:]+job["duration"]]
        combined_primary=np.concatenate([cue_primary,states[1:,0]])
        np.savez_compressed(out/"primary_history.npz",time_from_onset=combined_t,primary=combined_primary,
            x=grid.x,basal_state=basal,cue=cue,core_weights=core["weights"],cue_duration=job["duration"])
        status.update(complete=True,producer="finished",status="complete",readouts=readouts,
            combined_primary_path=str(out/"primary_history.npz"),combined_primary_sha256=sha(out/"primary_history.npz"),
            core={k:v for k,v in core.items() if k!="weights"},wall_seconds=time.monotonic()-began)
    except Exception as e:
        status.update(complete=True,producer="finished",status="error",error=f"{type(e).__name__}: {e}")
        write(out/"status.json",status);raise
    write(out/"status.json",status)
    print(json.dumps(dict(slot=slot,status=status["status"],wall_seconds=status["wall_seconds"])))


def submit(run,slots,label):
    verify(run);runner=run/"frozen_sources/code/mechanisms"/Path(__file__).name
    command="exec "+shlex.join(["env","OPENBLAS_NUM_THREADS=1","OMP_NUM_THREADS=1","MKL_NUM_THREADS=1",PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=01:00:00",
        "--job-name=vision-ring-duration","--array="+slots,"--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap",command]
    job=subprocess.check_output(args,text=True).strip().split(";")[0]
    write(run/f"launch_{label}.json",dict(job_id=job,command=args,manifest_sha256=sha(run/"manifest.json")))
    print(job)


def collect(run):
    m=verify(run);rows=[];statuses=[]
    for job in m["jobs"]:
        s=read(run/"workers"/f'{job["slot"]:02d}'/"status.json");assert s["status"]=="complete"
        assert sha(s["combined_primary_path"])==s["combined_primary_sha256"]
        statuses.append(s)
        for r in s["readouts"]:
            rows.append(dict(model=job["model"],circumference=job["circumference"],fraction=job["fraction"],
                cue_duration=job["duration"],recovery=job["recovery"],cue_free_time=r["cue_free_time"],
                original_core_amount=r["original_core_amount"],original_core_amount_fraction=r["original_core_amount_fraction"],
                Px=r["Px"],Py=r["Py"],R=r["R"],angle_deg=r["angle_deg"],
                resolved_peak_count=r["spatial"]["resolved_peak_count"],location_association=r["spatial"]["location_association"],
                history_path=s["combined_primary_path"],history_sha256=s["combined_primary_sha256"]))
    with (run/"results/duration_readouts.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write(run/"results/verification.json",dict(passed=True,complete=True,workers=len(statuses),
        max_pool_error=max(p["max_pool_error"] for s in statuses for p in s["phases"]),
        max_replayed_prefix_absolute_error=max(max(max(s.get("replayed_prefix_native_maxabs_at_0_1_10",[0.])),
            s.get("replayed_primary_prefix_check",{}).get("max_absolute_error",0.)) for s in statuses),
        max_replayed_native_0_1_10_error=max(max(s.get("replayed_prefix_native_maxabs_at_0_1_10",[0.])) for s in statuses),
        collector_sha256=sha(Path(__file__)),
        frozen_worker_sha256=m["frozen_source_sha256"]["code/mechanisms/"+Path(__file__).name],
        statuses=statuses))
    print(json.dumps(dict(complete=True,conditions=18,readouts=len(rows))))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("action",choices=["initialize","submit","worker","collect"])
    p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--slots",default="0-17");p.add_argument("--label",default="main");p.add_argument("--slot",type=int)
    a=p.parse_args();r=a.run.resolve()
    if a.action=="initialize":print(json.dumps(dict(jobs=len(initialize(r)["jobs"]))))
    elif a.action=="submit":submit(r,a.slots,a.label)
    elif a.action=="worker":worker(r,a.slot)
    else:collect(r)
