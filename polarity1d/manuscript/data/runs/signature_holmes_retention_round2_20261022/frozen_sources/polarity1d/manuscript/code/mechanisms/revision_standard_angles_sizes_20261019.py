"""Immutable three-size extension of the requested stationary-cue experiment."""
from pathlib import Path
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
import revision_standard_angles_duration_20261019 as base
from spring_ring_model_20261019 import Closure,physical_raised_cosine

ROOT=base.ROOT
DEFAULT=ROOT/"data/runs/revision_standard_angles_sizes_20261019"
OUT=base.TOPIC.parent/"standard_angles_sizes"
PYTHON=base.PYTHON
sha=base.sha
read=base.read
write=base.write
now=base.now


def freeze(run):
    if run.exists():raise FileExistsError("A fresh immutable campaign is required")
    for folder in ("inputs","frozen_sources","workers","logs","results"):(run/folder).mkdir(parents=True)
    oldrun=ROOT/"data/runs/chemical_ring_functional_v2_20261017";oldplan=read(oldrun/"manifest.json")
    paths=[Path(__file__),Path(base.__file__)]+[Path(__file__).with_name(name) for name in (
        "chemical_ring_functional_v2_20261017.py","ring_distance_prototype_20261017.py","ring_preparation_pilot_20261017.py",
        "chemical_ring_distance_pilot_20261017.py","length_dev_output_layout.py","spring_ring_model_20261019.py","spring_ring_closure_20261017.py")]
    for folder in (ROOT/"code/core/models",ROOT/"code/core/src/polarity1d"):paths+=list(folder.glob("*.py"))
    sources={}
    for path in paths:
        relative=path.relative_to(ROOT);target=run/"frozen_sources"/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
        sources[str(relative)]=sha(target)
    standard=pd.read_csv(base.TOPIC/"supporting/data/standard_angles_summary.csv")
    springtable=pd.read_csv(base.TOPIC.parent/"spring_distance/supporting/data/spring_distance_duration_summary.csv")
    matrix=pd.read_csv(oldrun/"results/matrix.csv")
    cases=[];jobs=[];contexts=[];inputs={};basals=[];rejected=[]
    for model in base.MODELS:
        spring=model=="spring"
        for size in ((1.,2.,3.) if spring else (20.,30.,40.)):
            if spring:
                original=ROOT/f"data/runs/revision_ring_spring_20261019/preparations/C{size:g}/qualified_basal.npz"
                statuspath=original.with_name("status.json");status=read(statuspath)
                assert status["status"]=="qualified" and status["complete"] and sha(original)==status["checkpoint_sha256"]
                qualification=dict(status="qualified",status_path=str(statuspath),status_sha256=sha(statuspath),
                    independent_holdout_path=str(ROOT/"data/runs/revision_ring_spring_20261019/results/root_independent_preparation_holdout.json"))
                qualification["independent_holdout_sha256"]=sha(qualification["independent_holdout_path"])
                scales=status["fixed_scales"]
            else:
                oldjob=next(j for j in oldplan["jobs"] if j["model"]==model and j["circumference"]==size)
                original=oldrun/oldjob["checkpoint"];assert sha(original)==oldjob["checkpoint_sha256"]
                if size in (20.,30.):
                    statuspath=Path(oldjob["original_checkpoint"]).with_name("status.json");status=read(statuspath)
                    assert status["status"]=="qualifying_single_site" and status["complete"] and status["checkpoint_sha256"]==sha(original)
                    assert any(a.get("verified_settled",False) for a in status["attempts"])
                else:
                    # C40 is exactly the accepted original80 state; qualification remains linked through its native lineage.
                    statuspath=Path(oldjob["existing_functional_status"]);status=read(statuspath);assert status["complete"]
                qualification=dict(status="qualified",status_path=str(statuspath),status_sha256=sha(statuspath),
                    native_grid_manifest=str(oldrun/"manifest.json"),native_grid_manifest_sha256=sha(oldrun/"manifest.json"),
                    original_checkpoint=oldjob["original_checkpoint"])
                scales=oldjob["fixed_scales"]
            target=run/"inputs"/f"{model}_C{size:g}_basal.npz";shutil.copy2(original,target);inputs[str(target.relative_to(run))]=sha(target)
            with np.load(target) as z:state=z["state"].copy()
            if spring:
                n=(len(state)-1)//2;q=(np.arange(n)+.5)/n;Cb=float(state[-1]);u,v=state[:n],state[n:-1]
                relevant=["code/mechanisms/spring_ring_model_20261019.py","code/mechanisms/spring_ring_closure_20261017.py","code/core/models/spring.py"]
            else:
                n=state.shape[1];Cb=size;grid=base.ring.RingGrid(Cb,n);q=grid.x/Cb;u,v=state[0],state[1]
                relevant=["code/mechanisms/ring_distance_prototype_20261017.py",f"code/core/models/{model}.py"]
                relevant += [p for p in sources if p.startswith("code/core/src/polarity1d/")]
            assert n==int(100*size if spring else 20*size)
            dynamics={rel:sources[rel] for rel in relevant};q0=float(q[0]);core=base.core_info(q,u,q0)
            M1=Cb*np.mean(u*np.exp(2j*np.pi*q));offset=float(np.rad2deg(np.angle(M1))-360*q0)
            basals.append(dict(model=model,size_parameter=size,cells=n,basal_circumference=Cb,input=str(target.relative_to(run)),
                input_sha256=sha(target),original_path=str(original),original_sha256=sha(original),qualification=qualification,
                preparation_axis_fraction=q0,initial_M1_absolute_phase_deg=float(np.rad2deg(np.angle(M1))),initial_M1_offset_from_preparation_deg=offset))
            for label,fraction in base.ANGLES:
                centre=(q0+fraction)%1;width=.2 if spring else 3.;peak=base.PEAKS[model]
                if spring:
                    signal,dC=physical_raised_cosine(Closure(size,n),state,centre,peak,width);evaluation="Continuous physical raised cosine at requested material centre; no cell snap"
                else:
                    signal,evaluation=base.chemical_cue(grid,fraction,peak,width);dC=np.zeros(n)
                cuefile=run/"inputs"/f"{model}_C{size:g}_{label}_cue.npz"
                np.savez_compressed(cuefile,material_positions=q,initial_cue=signal,initial_cue_dC=dC,centre_fraction=centre,initial_state=state,
                    prescribed_angle_radians=2*np.pi*fraction)
                inputs[str(cuefile.relative_to(run))]=sha(cuefile)
                left,right=core["left"]*Cb,core["right"]*Cb;d=fraction*Cb
                overlap=sum(max(0.,min(right,d+width/2+k*Cb)-max(left,d-width/2+k*Cb)) for k in (-1,0,1))
                patch=base.interval_weights(q,centre,width/(2*Cb));total_excess=Cb*np.mean(u-core["background"])
                context=dict(model=model,size_parameter=size,size_label=(f"C0={size:g}" if spring else f"C={size:g}"),circumference=Cb,
                    resting_circumference=size if spring else None,cells=n,angle_label=label,fraction=fraction,requested_angle_deg=360*fraction,
                    requested_angle_radians=2*np.pi*fraction,preparation_axis_fraction=q0,preparation_axis_absolute_deg=360*q0,
                    initial_M1_absolute_phase_deg=float(np.rad2deg(np.angle(M1))),initial_M1_offset_from_preparation_deg=offset,
                    target_angle_in_P_coordinates_deg=360*fraction-offset,cue_centre_fraction=centre,cue_centre_physical_at_onset=centre*Cb,
                    physical_distance_at_onset=d,full_physical_support=width,initial_angular_support_deg=360*width/Cb,requested_peak=peak,
                    initial_sampled_cue_peak=float(signal.max()),initial_cue_integral=float(Cb*signal.mean()),
                    initial_cue_shape_sha256=hashlib.sha256(np.asarray(signal,dtype='<f8').tobytes()).hexdigest(),cue_evaluation=evaluation,
                    original_core_left_physical=left,original_core_right_physical=right,original_core_physical_width=right-left,
                    cue_support_core_intersection=overlap,overlap_fraction_cue_support=overlap/width,overlap_fraction_original_core=overlap/(right-left),
                    initial_overlap=bool(overlap>1e-12),separate_challenger_interpretation_eligible=bool(overlap<=1e-12),
                    cue_support_contains_preparation_axis=bool(d<=width/2),
                    initial_unit_peak_cue_excess_fraction_total=float(Cb*np.mean((signal/peak)*(u-core["background"]))/total_excess),
                    initial_geometric_support_excess_fraction_total=float(Cb*((u-core["background"])@patch)/total_excess),
                    initial_cue_weighted_active=float(signal@u/signal.sum()),initial_cue_weighted_inactive=float(signal@v/signal.sum()),
                    input=str(target.relative_to(run)),input_sha256=sha(target),cue_input=str(cuefile.relative_to(run)),cue_input_sha256=sha(cuefile),fixed_scales=scales)
                contexts.append(context)
                for duration in ((.1,1.,10.,100.,1000.) if spring else (1.,10.,100.,300.,1000.)):
                    case=dict(**context,condition_id=len(cases),cue_duration=duration,recovery=13000.,native_units_uncalibrated=True)
                    candidate=standard[(standard.model==model)&np.isclose(standard.fraction,fraction)&np.isclose(standard.cue_duration,duration)]
                    candidate=candidate[np.isclose(candidate.resting_circumference,size)] if spring else candidate[np.isclose(candidate.circumference,size)]
                    reuse=None;originid=None
                    if len(candidate):reuse=candidate.iloc[0].to_dict();originid=int(candidate.iloc[0].condition_id)
                    elif spring:
                        subset=springtable[np.isclose(springtable.resting_circumference,size)&np.isclose(springtable.fraction,fraction)&np.isclose(springtable.cue_duration,duration)]
                        if len(subset):reuse=subset.iloc[0].to_dict()
                    elif duration==1000.:
                        subset=matrix[(matrix.model==model)&np.isclose(matrix.circumference,size)&np.isclose(matrix.fraction,fraction)&(matrix.status=="complete")]
                        if len(subset):
                            r=subset.iloc[0];reuse=dict(actual_cue_path=r.cue_path,actual_cue_sha256=sha(r.cue_path),actual_recovery_path=r.recovery_path,actual_recovery_sha256=sha(r.recovery_path))
                    proof=None
                    if reuse is not None:
                        try:
                            # Existing helper checks actual archive/state/cue/duration/native dynamics exactly.
                            proof=reuse_proof(case,reuse,run,dynamics)
                        except Exception as error:rejected.append(dict(condition_id=case["condition_id"],reason=f"{type(error).__name__}: {error}"))
                    if proof is not None:case.update(execution="exact_archive_reuse",reuse=proof,source_condition_id=originid)
                    else:case.update(execution="new_native_branch",worker_slot=len(jobs));jobs.append(case)
                    cases.append(case)
    assert len(cases)==240 and len(contexts)==48 and len(basals)==12
    manifest=dict(generated_at=now(),experiment="Three-size stationary-cue angle-duration extension",requested_conditions=240,
        cases=cases,jobs=jobs,basals=basals,source_sha256=sources,input_sha256=inputs,rejected_reuse_candidates=rejected,
        fixed_chemical_sizes=[20,30,40],resting_Spring_sizes=[1,2,3],angles_radians=[np.pi,np.pi/2,np.pi/4,np.pi/8],
        overlap_cases_retained=True,excluded_for_overlap=0,solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
        cue_full_physical_support=dict(chemical=3.,spring=.2),
        angular_reference="Declared preparation axis; exact measured initial physical M1 phase offset retained. Continuous centres are not snapped.",
        interpretation="Same native stationary physical-support cue protocol; all initial overlaps retained, without separate-challenger inference. Native clocks, widths and total material differ across sizes/models.")
    write(run/"manifest.json",manifest);pd.DataFrame(contexts).to_csv(run/"results/initial_context.csv",index=False)
    write(run/"results/basal_inventory.json",dict(generated_at=now(),basals=basals))
    summary=dict(generated_at=now(),requested=240,exact_reuse=240-len(jobs),new_native=len(jobs),initial_overlap_conditions=sum(c["initial_overlap"] for c in cases),
        excluded_for_overlap=0,qualified_basals=12,manifest_sha256=sha(run/"manifest.json"))
    write(run/"freeze.json",summary);print(json.dumps(summary,indent=2))


def reuse_proof(case,candidate,run,dynamics):
    # Generalized only for archive centre-key spelling and Spring resting size.
    cue=Path(candidate["actual_cue_path"]);release=Path(candidate["actual_recovery_path"])
    assert sha(cue)==candidate["actual_cue_sha256"] and sha(release)==candidate["actual_recovery_sha256"]
    with np.load(run/case["input"]) as z:initial=z["state"].copy()
    np.testing.assert_array_equal(base.full_state(cue,0.),initial)
    np.testing.assert_array_equal(base.full_state(cue,case["cue_duration"]),base.full_state(release,0.))
    with np.load(cue) as z:
        assert z["time"][-1]==case["cue_duration"]
        if case["model"]=="spring":
            key="cue_centre_fraction" if "cue_centre_fraction" in z else "cue_center_fraction"
            assert float(z[key])==case["cue_centre_fraction"] and float(z["cue_peak"])==4.
            if "full_physical_support" in z:assert float(z["full_physical_support"])==.2
        else:
            with np.load(run/case["cue_input"]) as expected:np.testing.assert_array_equal(z["stimulus"],expected["initial_cue"])
    with np.load(release) as z:
        assert z["time"][-1]==13000.
        assert (float(z["cue_peak"])==0.) if "cue_peak" in z else not np.count_nonzero(z["stimulus"])
    origin=base.source_run(cue)
    for rel,digest in dynamics.items():assert sha(origin/"frozen_sources"/rel)==digest,(origin,rel)
    if case["model"]=="spring":
        old=read(origin/"manifest.json")
        # Original standard-angle manifest stores physical support per case; earlier campaign stores it globally.
        if "cue_full_physical_support" in old:assert old["cue_full_physical_support"]==.2
        else:assert any(c["model"]=="spring" and c["full_physical_support"]==.2 and c["cue_centre_fraction"]==case["cue_centre_fraction"] for c in old["cases"])
    return dict(status="exact_archive_reuse",actual_cue_path=str(cue),actual_cue_sha256=sha(cue),actual_recovery_path=str(release),actual_recovery_sha256=sha(release),
        source_run=str(origin),exact_initial_full_state=True,exact_phase_restart=True,exact_cue_profile_or_continuous_parameters=True,
        exact_duration_and_recovery_horizon=True,unchanged_native_dynamics_verified=True)


def submit(run):
    assert not (run/"launch.json").exists()
    manifest=read(run/"manifest.json");assert sha(run/"manifest.json")==read(run/"freeze.json")["manifest_sha256"]
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job=subprocess.check_output(["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:30:00",f'--array=0-{len(manifest["jobs"])-1}',
        "--job-name=stationary-three-sizes","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command],text=True).strip()
    write(run/"launch.json",dict(**read(run/"freeze.json"),array_job_id=job,submitted_at=now()));print(json.dumps(read(run/"launch.json"),indent=2))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Native integration requires SLURM")
    manifest=read(run/"manifest.json");case=manifest["jobs"][slot]
    for rel,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/rel)==digest
    for rel in (case["input"],case["cue_input"]):assert sha(run/rel)==manifest["input_sha256"][rel]
    directory=run/"workers"/f"{slot:03d}";directory.mkdir();began=time.monotonic()
    status=dict(generated_at=now(),producer="running",status="challenge",complete=False,condition_id=case["condition_id"],phases=[],slurm_job_id=os.environ["SLURM_JOB_ID"])
    write(directory/"status.json",status)
    try:
        with np.load(run/case["input"]) as z:state=z["state"].copy()
        with np.load(run/case["cue_input"]) as z:q=z["material_positions"].copy();fixedcue=z["initial_cue"].copy()
        spring=case["model"]=="spring";n=case["cells"]
        if spring:model=Closure(case["resting_circumference"],n);amount0=model.amount(state)
        else:model=base.ring.native_model(case["model"]);grid=base.ring.RingGrid(case["circumference"],n);pools=base.ring.pool_amounts(model,grid,state)
        for name,duration,peak in (("challenge",case["cue_duration"],case["requested_peak"]),("recovery",13000.,0.)):
            status.update(status=name,generated_at=now());write(directory/"status.json",status)
            if spring:
                def cue(y):return physical_raised_cosine(model,y,case["cue_centre_fraction"],peak,.2) if peak else (0.,0.)
                solution=solve_ivp(lambda t,y:model.rhs(y,cue(y)[0]),(0.,duration),state,method="BDF",jac=lambda t,y:model.jacobian(y,*cue(y)),
                    rtol=1e-9,atol=1e-11,max_step=10.,t_eval=base.samples(duration))
                if not solution.success:raise RuntimeError(solution.message)
                tt=solution.t;yy=solution.y.T;amount=yy[:,-1]*(yy[:,:n]+yy[:,n:-1]).mean(axis=1)
                stats=dict(max_pool_error=float(abs(amount/amount0-1).max()),minimum=float(yy[:,:-1].min()),nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu)
            else:tt,yy,stats=base.campaign.integrate(model,grid,state,fixedcue if peak else np.zeros(n),duration,base.samples(duration),pools,began+1700.)
            assert tt[-1]==duration and stats["max_pool_error"]<1e-7 and stats["minimum"]>=-5e-10 and np.isfinite(yy).all()
            path=directory/(name+".npz")
            if spring:np.savez_compressed(path,time=tt,states=yy,material_positions=q,cue_peak=peak,cue_centre_fraction=case["cue_centre_fraction"],full_physical_support=.2,resting_circumference=case["resting_circumference"])
            else:np.savez_compressed(path,time=tt,states=yy,x=grid.x,stimulus=fixedcue if peak else np.zeros(n))
            status["phases"].append(dict(kind=name,path=str(path.relative_to(run)),sha256=sha(path),**stats));state=yy[-1].copy()
        status.update(status="complete",complete=True,producer="finished",elapsed_seconds=time.monotonic()-began)
    except Exception as error:status.update(status="error",complete=True,producer="finished",error=f"{type(error).__name__}: {error}");raise
    finally:status["generated_at"]=now();write(directory/"status.json",status)


def collect(run):
    manifest=read(run/"manifest.json");rows=[]
    for case in manifest["cases"]:
        row={k:v for k,v in case.items() if k!="reuse"}
        if case["execution"]=="exact_archive_reuse":row.update(status="complete",**{k:v for k,v in case["reuse"].items() if k.startswith("actual_")})
        else:
            directory=run/"workers"/f'{case["worker_slot"]:03d}';path=directory/"status.json";status=read(path) if path.exists() else dict(status="not_started")
            row["status"]=status["status"]
            if status["status"]=="complete":
                for phase in status["phases"]:
                    prefix="actual_cue" if phase["kind"]=="challenge" else "actual_recovery";row[prefix+"_path"]=str(run/phase["path"]);row[prefix+"_sha256"]=phase["sha256"]
        rows.append(row)
    pd.DataFrame(rows).to_csv(run/"results/condition_inventory.csv",index=False)
    counts=dict(generated_at=now(),requested=240,complete=sum(r["status"]=="complete" for r in rows),exact_reuse=240-len(manifest["jobs"]),new_native=len(manifest["jobs"]),
        initial_overlap_conditions=sum(r["initial_overlap"] for r in rows),excluded_for_overlap=0,states=pd.Series([r["status"] for r in rows]).value_counts().to_dict())
    write(run/"results/collection.json",counts);print(json.dumps(counts,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("action",choices=["freeze","submit","worker","collect"]);p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--slot",type=int);a=p.parse_args();run=a.run.resolve()
    if a.action=="freeze":freeze(run)
    elif a.action=="submit":submit(run)
    elif a.action=="collect":collect(run)
    else:worker(run,a.slot)
