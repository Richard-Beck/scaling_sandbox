"""Frozen four-model experiment at pi, pi/2, pi/4 and pi/8.

All initial overlap cases are retained. The nominal preparation axis defines
cue placement; the measured M1 offset defines exact targets in P coordinates.
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
import chemical_ring_functional_v2_20261017 as campaign
import ring_distance_prototype_20261017 as ring
from spring_ring_model_20261019 import Closure,physical_raised_cosine

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"data/runs/revision_standard_angles_duration_20261019"
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/standard_angles"
PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS=("wave_pinning","goryachev","holmes_model3","spring")
ANGLES=(("pi",.5),("pi_2",.25),("pi_4",.125),("pi_8",.0625))
PEAKS=dict(wave_pinning=.1,goryachev=6.4,holmes_model3=1.,spring=4.)


def now():return datetime.now(timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,z):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");tmp.replace(path)


def samples(end):
    return np.unique(np.r_[0.,np.geomspace(1e-4,min(.1,end),25),np.arange(.1,min(10.,end),.1),
        np.arange(10.,min(100.,end),1.),np.arange(100.,min(1000.,end),5.),np.arange(1000.,end,50.),end])


def core_info(q,u,origin):
    arc=(q-origin+.5)%1-.5;order=np.argsort(arc);x,p=arc[order],u[order]
    level=float(p.min()+.5*np.ptp(p));inside=np.flatnonzero(p>=level);a,b=int(inside[0]),int(inside[-1])
    assert np.all(np.diff(inside)==1) and a>0 and b<len(q)-1
    left=x[a-1]+(level-p[a-1])/(p[a]-p[a-1])*(x[a]-x[a-1])
    right=x[b]+(level-p[b])/(p[b+1]-p[b])*(x[b+1]-x[b])
    dq=1/len(q);weights=np.maximum(0,np.minimum(arc+dq/2,right)-np.maximum(arc-dq/2,left))
    return dict(left=float(left),right=float(right),weights=weights,width=float(right-left),background=float(p.min()),level=level)


def interval_weights(q,centre,halfwidth):
    arc=(q-centre+.5)%1-.5;dq=1/len(q)
    return np.maximum(0,np.minimum(arc+dq/2,halfwidth)-np.maximum(arc-dq/2,-halfwidth))


def chemical_cue(grid,fraction,peak,width):
    # All four requested fractions are exact native cell offsets on N800.
    # This uses an exact cyclic translation, never rounding an arbitrary centre.
    offset=fraction*grid.cells
    if offset==int(offset):
        return ring.ring_cue(grid,int(offset),width,peak),"exact native cyclic translation at requested on-grid centre"
    centre=grid.x[0]+fraction*grid.circumference
    distance=grid.distance(grid.x,centre)
    signal=np.where(distance<width/2,peak*.5*(1+np.cos(2*np.pi*distance/width)),0.)
    return signal,"continuous periodic raised cosine at requested centre; no cell snap"


def full_state(path,t):
    with np.load(path) as z:
        tt=z["time"] if "states" in z else z["snapshot_times"]
        found=np.flatnonzero(np.isclose(tt,t,rtol=0,atol=1e-10));assert len(found)==1,(path,t)
        return z["states" if "states" in z else "snapshot_states"][int(found[0])].copy()


def source_run(path):
    for parent in Path(path).parents:
        if (parent/"manifest.json").exists():return parent
    raise FileNotFoundError(path)


def reuse_proof(case,candidate,run,dynamics):
    cue_path=Path(candidate["actual_cue_path"]);recovery_path=Path(candidate["actual_recovery_path"])
    assert sha(cue_path)==candidate["actual_cue_sha256"] and sha(recovery_path)==candidate["actual_recovery_sha256"]
    with np.load(run/case["input"]) as z:initial=z["state"].copy()
    np.testing.assert_array_equal(full_state(cue_path,0.),initial)
    np.testing.assert_array_equal(full_state(cue_path,case["cue_duration"]),full_state(recovery_path,0.))
    with np.load(cue_path) as z:
        assert z["time"][-1]==case["cue_duration"]
        if case["model"]!="spring":
            with np.load(run/case["cue_input"]) as expected:np.testing.assert_array_equal(z["stimulus"],expected["initial_cue"])
        else:
            assert float(z["cue_peak"])==4. and float(z["cue_center_fraction"])==case["cue_centre_fraction"]
    with np.load(recovery_path) as z:
        assert z["time"][-1]==13000.
        if "stimulus" in z:assert np.count_nonzero(z["stimulus"])==0
        elif "cue_peak" in z:assert float(z["cue_peak"])==0.
    origin=source_run(cue_path)
    for relative,digest in dynamics.items():
        frozen=origin/"frozen_sources"/relative
        assert frozen.exists() and sha(frozen)==digest,(origin,relative)
    if case["model"]=="spring":
        old=read(origin/"manifest.json");assert old["cue_full_physical_support"]==.2 and old["cue_peak"]==4.
    return dict(status="exact_archive_reuse",actual_cue_path=str(cue_path),actual_cue_sha256=sha(cue_path),
        actual_recovery_path=str(recovery_path),actual_recovery_sha256=sha(recovery_path),source_run=str(origin),
        exact_initial_full_state=True,exact_cue_recovery_full_state_continuity=True,exact_cue_profile_or_continuous_parameters=True,
        exact_duration_and_recovery_horizon=True,unchanged_native_dynamics_verified=True)


def launch(run):
    if run.exists():raise FileExistsError("Fresh campaign required; existing evidence immutable")
    for folder in ("inputs","frozen_sources","workers","logs","results"):(run/folder).mkdir(parents=True)
    paths=[Path(__file__)]+[Path(__file__).with_name(name) for name in ("chemical_ring_functional_v2_20261017.py",
        "ring_distance_prototype_20261017.py","ring_preparation_pilot_20261017.py","chemical_ring_distance_pilot_20261017.py",
        "length_dev_output_layout.py","spring_ring_model_20261019.py","spring_ring_closure_20261017.py")]
    for folder in (ROOT/"code/core/models",ROOT/"code/core/src/polarity1d"):paths+=list(folder.glob("*.py"))
    frozen={}
    for path in paths:
        rel=path.relative_to(ROOT);target=run/"frozen_sources"/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);frozen[str(rel)]=sha(target)
    oldtable=pd.read_csv(ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/distance_response/supporting/data/distance_duration_all_models.csv")
    oldmatrix=pd.read_csv(ROOT/"data/runs/chemical_ring_functional_v2_20261017/results/matrix.csv")
    oldplan=read(ROOT/"data/runs/chemical_ring_functional_v2_20261017/manifest.json")
    cases=[];jobs=[];contexts=[];inputs={};rejected=[]
    for model in MODELS:
        spring=model=="spring"
        original=ROOT/("data/runs/revision_ring_spring_20261019/preparations/C2/qualified_basal.npz" if spring else f"data/runs/revision_chemical_distance_duration_20261019/inputs/{model}_C40.npz")
        target=run/"inputs"/(model+"_basal.npz");shutil.copy2(original,target);inputs[str(target.relative_to(run))]=sha(target)
        with np.load(target) as z:state=z["state"].copy()
        if spring:
            n=(len(state)-1)//2;q=(np.arange(n)+.5)/n;Cb=float(state[-1]);active=state[:n];inactive=state[n:-1];scales=np.r_[np.full(2*n,2.),2.].tolist()
            relevant=["code/mechanisms/spring_ring_model_20261019.py","code/mechanisms/spring_ring_closure_20261017.py","code/core/models/spring.py"]
        else:
            n=state.shape[1];grid=ring.RingGrid(40,n);q=grid.x/40;Cb=40.;active,inactive=state[0],state[1]
            scales=next(j["fixed_scales"] for j in oldplan["jobs"] if j["model"]==model and j["circumference"]==40)
            relevant=["code/mechanisms/ring_distance_prototype_20261017.py",f"code/core/models/{model}.py"]
            relevant += [p for p in frozen if p.startswith("code/core/src/polarity1d/")]
        dynamics={key:frozen[key] for key in relevant};q0=float(q[0]);core=core_info(q,active,q0)
        M1=Cb*np.mean(active*np.exp(2j*np.pi*q));offset=float(np.rad2deg(np.angle(M1))-360*q0)
        for label,fraction in ANGLES:
            centre=(q0+fraction)%1;width=.2 if spring else 3.;peak=PEAKS[model]
            if spring:
                signal,dC=physical_raised_cosine(Closure(2.,n),state,centre,peak,width);evaluation="continuous physical raised cosine at exact material centre; no cell snap"
            else:signal,evaluation=chemical_cue(grid,fraction,peak,width);dC=np.zeros(n)
            cuefile=run/"inputs"/f"{model}_{label}_cue.npz"
            np.savez_compressed(cuefile,material_positions=q,initial_cue=signal,initial_cue_dC=dC,centre_fraction=centre,
                prescribed_angle_radians=2*np.pi*fraction,initial_state=state)
            inputs[str(cuefile.relative_to(run))]=sha(cuefile)
            left,right=core["left"]*Cb,core["right"]*Cb;d=fraction*Cb
            overlap=sum(max(0.,min(right,d+width/2+k*Cb)-max(left,d-width/2+k*Cb)) for k in (-1,0,1))
            patch=interval_weights(q,centre,width/(2*Cb));total_excess=Cb*np.mean(active-core["background"])
            weighted_excess=Cb*np.mean((signal/peak)*(active-core["background"]))
            context=dict(model=model,circumference=Cb,resting_circumference=2. if spring else None,cells=n,
                angle_label=label,fraction=fraction,requested_angle_deg=360*fraction,requested_angle_radians=2*np.pi*fraction,
                preparation_axis_fraction=q0,preparation_axis_absolute_deg=360*q0,initial_M1_absolute_phase_deg=float(np.rad2deg(np.angle(M1))),
                initial_M1_offset_from_preparation_deg=offset,target_angle_in_P_coordinates_deg=360*fraction-offset,
                cue_centre_fraction=centre,cue_centre_physical_at_onset=centre*Cb,physical_distance_at_onset=d,
                full_physical_support=width,initial_angular_support_deg=360*width/Cb,requested_peak=peak,
                initial_sampled_cue_peak=float(signal.max()),initial_cue_integral=float(Cb*signal.mean()),
                initial_cue_shape_sha256=hashlib.sha256(np.asarray(signal,dtype='<f8').tobytes()).hexdigest(),cue_evaluation=evaluation,
                original_core_left_physical=left,original_core_right_physical=right,original_core_physical_width=(right-left),
                cue_support_core_intersection=overlap,overlap_fraction_cue_support=overlap/width,overlap_fraction_original_core=overlap/(right-left),
                initial_overlap=bool(overlap>1e-12),separate_challenger_interpretation_eligible=bool(overlap<=1e-12),
                initial_unit_peak_cue_excess_fraction_total=float(weighted_excess/total_excess),
                initial_geometric_support_excess_fraction_total=float(Cb*((active-core["background"])@patch)/total_excess),
                initial_cue_weighted_active=float(signal@active/signal.sum()),initial_cue_weighted_inactive=float(signal@inactive/signal.sum()),
                input=str(target.relative_to(run)),input_sha256=sha(target),cue_input=str(cuefile.relative_to(run)),cue_input_sha256=sha(cuefile),fixed_scales=scales)
            contexts.append(context)
            for duration in ((.1,1.,10.,100.,1000.) if spring else (1.,10.,100.,300.,1000.)):
                case=dict(**context,condition_id=len(cases),cue_duration=duration,recovery=13000.,native_units_uncalibrated=True)
                candidate=oldtable[(oldtable.model==model)&np.isclose(oldtable.fraction,fraction)&np.isclose(oldtable.cue_duration,duration)&(oldtable.status=="complete")]
                if spring:candidate=candidate[np.isclose(candidate.resting_circumference,2.)]
                reuse=None
                if len(candidate):reuse=candidate.iloc[0].to_dict()
                elif not spring and duration==1000.:
                    subset=oldmatrix[(oldmatrix.model==model)&(oldmatrix.circumference==40)&np.isclose(oldmatrix.fraction,fraction)&(oldmatrix.status=="complete")]
                    if len(subset):
                        r=subset.iloc[0];reuse=dict(actual_cue_path=r.cue_path,actual_cue_sha256=sha(r.cue_path),actual_recovery_path=r.recovery_path,actual_recovery_sha256=sha(r.recovery_path))
                if reuse is not None:
                    try:proof=reuse_proof(case,reuse,run,dynamics)
                    except Exception as error:rejected.append(dict(condition_id=case["condition_id"],reason=f"{type(error).__name__}: {error}"));proof=None
                else:proof=None
                if proof is not None:case.update(execution="exact_archive_reuse",reuse=proof)
                else:
                    case.update(execution="new_native_branch",worker_slot=len(jobs));jobs.append(case)
                cases.append(case)
    manifest=dict(generated_at=now(),experiment="Four-model standardized-angle finite-duration experiment",cases=cases,jobs=jobs,
        requested_conditions=len(cases),overlap_cases_retained=True,excluded_for_overlap=0,angles_radians=[np.pi,np.pi/2,np.pi/4,np.pi/8],
        source_sha256=frozen,input_sha256=inputs,rejected_reuse_candidates=rejected,
        angular_reference="Declared nominal preparation axis, never a fresh argmax. P angles are relative to measured physical M1(0); exact target offsets are stored.",
        source_original_basals="Exact current C40 chemical and C0=2 Spring basals; no reflection or remapping",
        solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
        interpretation="Shared centre angles with original model-specific physical widths/peaks/native clocks. Initial overlap remains measured; fixed-region readouts are not distinct challenger/site lineage outcomes.")
    assert len(cases)==80
    write(run/"manifest.json",manifest);pd.DataFrame(contexts).to_csv(run/"results/initial_context.csv",index=False)
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job=subprocess.check_output(["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:30:00",f"--array=0-{len(jobs)-1}",
        "--job-name=standard-angle-duration","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command],text=True).strip()
    write(run/"launch.json",dict(generated_at=now(),array_job_id=job,requested=80,new_native=len(jobs),exact_reuse=80-len(jobs),manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(read(run/"launch.json"),indent=2))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Native integration requires SLURM")
    manifest=read(run/"manifest.json");case=manifest["jobs"][slot]
    for rel,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(run/rel)==digest
    directory=run/"workers"/f"{slot:03d}";directory.mkdir();began=time.monotonic()
    record=dict(generated_at=now(),producer="running",complete=False,status="challenge",condition_id=case["condition_id"],phases=[],slurm_job_id=os.environ["SLURM_JOB_ID"])
    write(directory/"status.json",record)
    try:
        with np.load(run/case["input"]) as z:state=z["state"].copy()
        with np.load(run/case["cue_input"]) as z:q=z["material_positions"].copy();fixedcue=z["initial_cue"].copy()
        spring=case["model"]=="spring";n=case["cells"]
        if spring:model=Closure(2.,n);amount0=model.amount(state)
        else:model=ring.native_model(case["model"]);grid=ring.RingGrid(40,n);pools=ring.pool_amounts(model,grid,state)
        for name,duration,peak in (("challenge",case["cue_duration"],case["requested_peak"]),("recovery",13000.,0.)):
            record.update(status=name,generated_at=now());write(directory/"status.json",record)
            if spring:
                def cue(y):return physical_raised_cosine(model,y,case["cue_centre_fraction"],peak,.2) if peak else (0.,0.)
                solution=solve_ivp(lambda t,y:model.rhs(y,cue(y)[0]),(0.,duration),state,method="BDF",jac=lambda t,y:model.jacobian(y,*cue(y)),
                    rtol=1e-9,atol=1e-11,max_step=10.,t_eval=samples(duration))
                if not solution.success:raise RuntimeError(solution.message)
                tt=solution.t;yy=solution.y.T;amount=yy[:,-1]*(yy[:,:n]+yy[:,n:-1]).mean(axis=1)
                stats=dict(max_pool_error=float(abs(amount/amount0-1).max()),minimum=float(yy[:,:-1].min()),nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu)
            else:tt,yy,stats=campaign.integrate(model,grid,state,fixedcue if peak else np.zeros(n),duration,samples(duration),pools,began+1700.)
            assert tt[-1]==duration and stats["max_pool_error"]<1e-7 and stats["minimum"]>=-5e-10 and np.isfinite(yy).all()
            archive=directory/(name+".npz")
            if spring:np.savez_compressed(archive,time=tt,states=yy,material_positions=q,cue_peak=peak,cue_centre_fraction=case["cue_centre_fraction"],full_physical_support=.2)
            else:np.savez_compressed(archive,time=tt,states=yy,x=grid.x,stimulus=fixedcue if peak else np.zeros(n))
            record["phases"].append(dict(kind=name,path=str(archive.relative_to(run)),sha256=sha(archive),**stats));state=yy[-1].copy()
        record.update(status="complete",complete=True,producer="finished",elapsed_seconds=time.monotonic()-began)
    except Exception as error:
        record.update(status="error",complete=True,producer="finished",error=f"{type(error).__name__}: {error}");raise
    finally:record["generated_at"]=now();write(directory/"status.json",record)
    print(json.dumps(record,indent=2))


def collect(run):
    manifest=read(run/"manifest.json");rows=[]
    for case in manifest["cases"]:
        row={**case};row.pop("reuse",None)
        if case["execution"]=="exact_archive_reuse":
            row.update(status="complete",**{k:v for k,v in case["reuse"].items() if k.startswith("actual_")})
        else:
            directory=run/"workers"/f'{case["worker_slot"]:03d}';path=directory/"status.json"
            record=read(path) if path.exists() else dict(status="not_started",complete=False)
            row["status"]=record["status"]
            if record["status"]=="complete":
                for phase in record["phases"]:
                    prefix="actual_cue" if phase["kind"]=="challenge" else "actual_recovery"
                    row[prefix+"_path"]=str(run/phase["path"]);row[prefix+"_sha256"]=phase["sha256"]
        rows.append(row)
    pd.DataFrame(rows).to_csv(run/"results/condition_inventory.csv",index=False)
    record=dict(generated_at=now(),requested=80,complete=sum(r["status"]=="complete" for r in rows),
        exact_reuse=sum(r["execution"]=="exact_archive_reuse" for r in rows),new_native=len(manifest["jobs"]),
        initial_overlap_conditions=sum(r["initial_overlap"] for r in rows),excluded_for_overlap=0,
        states={key:sum(r["status"]==key for r in rows) for key in set(r["status"] for r in rows)})
    write(run/"results/collection.json",record);print(json.dumps(record,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("action",choices=["launch","worker","collect"]);p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--slot",type=int);a=p.parse_args()
    if a.action=="launch":launch(a.run.resolve())
    elif a.action=="collect":collect(a.run.resolve())
    else:worker(a.run.resolve(),a.slot)
