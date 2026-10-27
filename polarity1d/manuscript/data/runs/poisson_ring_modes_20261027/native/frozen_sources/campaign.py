"""Frozen broad-strength pilot around the unchanged four-model F4 noise laws.

prepare/check perform no integrations. Native work requires worker under SLURM;
submit requires an independent manifest-specific preflight acceptance receipt.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import copy
import csv
import hashlib
import importlib
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time
import traceback

for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
MS = ROOT / "polarity1d/manuscript"
TOPIC = MS / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/broad_strength_campaign_20261019"
DEFAULT = MS / "data/runs/revision_noise_strength_broad_20261019"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
GRID = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, .1, .3, 1.)
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring_ring")
REFERENCES = dict(wave_pinning=.1, goryachev=6.4, holmes_model3=1., spring_ring=4.)
LEGACY = dict(chemical=MS / "data/runs/vision_noise_main_20261017",
              spring=MS / "data/runs/revision_spring_noise_main_20261019")
STOP = False


def now(): return datetime.now(timezone.utc).isoformat()
def read(p): return json.loads(Path(p).read_text())
def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""): h.update(block)
    return h.hexdigest()
def ah(a): return hashlib.sha256(np.asarray(a,dtype="<f8",order="C").tobytes()).hexdigest()
def write(p, value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+".tmp");tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");tmp.replace(p)
def csvwrite(p, records):
    with Path(p).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in records for k in r)))
        w.writeheader();w.writerows(records)


def prepare(run):
    run=Path(run).resolve()
    assert not run.exists(), "Choose a fresh immutable preparation path"
    assert run.is_relative_to(MS/"data/runs")
    plan=TOPIC/"PLAN.md";assert plan.is_file()
    for name in ("inputs","frozen_sources","protocol","workers","logs","results"):
        (run/name).mkdir(parents=True)
    shutil.copy2(__file__,run/"frozen_sources/campaign.py")
    shutil.copy2(plan,run/"protocol/PLAN.md")
    manifests={platform:read(p/"manifest.json") for platform,p in LEGACY.items()}
    sources={"frozen_sources/campaign.py":sha(__file__)}
    for platform,old in LEGACY.items():
        for rel,digest in manifests[platform]["source_sha256"].items():
            p=old/"frozen_sources"/rel;assert sha(p)==digest
            target=run/"frozen_sources"/platform/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
            sources[str(target.relative_to(run))]=digest
    basals=[];reuse=[];runtimes=[]
    for platform,old in LEGACY.items():
        m=manifests[platform]
        for j in m["jobs"]:
            model=j.get("model","spring_ring");size=float(j.get("circumference",j.get("resting_circumference")))
            inputrel=j["input"] if platform=="chemical" else j["path"]
            digest=j["input_sha256"] if platform=="chemical" else j["sha256"]
            identity=(model,size)
            found=[b for b in basals if (b["model"],b["size"])==identity]
            if not found:
                source=old/inputrel;assert sha(source)==digest
                target=run/"inputs"/f"{model}_C{size:g}.npz";shutil.copy2(source,target)
                with np.load(source) as z: state=z["state"].copy()
                b=dict(model=model,size=size,platform=platform,cells=j["cells"],input=str(target.relative_to(run)),
                       input_sha256=digest,state_sha256=ah(state),legacy_input=str(source),
                       basal_circumference=float(state[-1]) if platform=="spring" else size)
                basals.append(b)
            else: assert found[0]["input_sha256"]==digest
            folder=old/"workers"/f"{j['slot']:04d}";status=read(folder/"status.json")
            assert status["status"]=="complete" and status["producer"]=="terminal"
            assert status["job"]==j and status["manifest_sha256"]==sha(old/"manifest.json")
            for name,digest in status["output_sha256"].items(): assert sha(folder/name)==digest
            kind=j["kind"];fraction=float(j.get("mean_fraction",j.get("fraction")))
            reuse.append(dict(model=model,size=size,platform=platform,kind=kind,fraction=fraction,seed=j["seed"],
                replicate=j.get("replicate",0),legacy_slot=j["slot"],legacy_run=str(old),status_sha256=sha(folder/"status.json"),
                input_sha256=found[0]["input_sha256"] if found else basals[-1]["input_sha256"],
                trace=str(folder/"trace.npz"),trace_sha256=status["output_sha256"]["trace.npz"],
                primary=str(folder/"primary.npz"),primary_sha256=status["output_sha256"]["primary.npz"]))
            allocated=sum(p.stat().st_blocks*512 for p in folder.iterdir() if p.is_file())
            runtimes.append(dict(model=model,size=size,kind=kind,fraction=fraction,seed=j["seed"],
                elapsed_seconds=status["elapsed_seconds"],allocated_bytes=allocated,status_path=str(folder/"status.json")))
    assert len(basals)==12 and len(reuse)==324
    jobs=[]
    for mi,model in enumerate(MODELS):
        sizes=(1.,2.,3.) if model=="spring_ring" else (20.,40.,60.)
        for si,size in enumerate(sizes):
            b,=[b for b in basals if b["model"]==model and b["size"]==size]
            for fi,f in enumerate(GRID):
                if f in (.001,.01):continue
                for rep in range(2):
                    seed=202610192000000+mi*100000+si*10000+fi*100+rep
                    dense=bool(size==sizes[-1] and f==1. and rep==0)
                    jobs.append(dict(**b,kind="noise",fraction=f,mean=REFERENCES[model]*f,seed=seed,replicate=rep,
                        variant="standard",dense_observation=dense,paired_standard_slot=None))
                jobs.append(dict(**b,kind="mean",fraction=f,mean=REFERENCES[model]*f,seed=0,replicate=0,
                    variant="standard",dense_observation=False,paired_standard_slot=None))
    assert len(jobs)==252
    for i,j in enumerate(jobs):j["slot"]=i
    for j in list(jobs):
        if j["dense_observation"]:
            jobs.append(dict(j,slot=len(jobs),variant="tight_ode",paired_standard_slot=j["slot"]))
    assert len(jobs)==256
    planvalues=dict(grid=list(GRID),sigma=1.,new_level_pilot_seeds=2,anchors_retained_seeds=12,
        final_uniform_replication_target=12,scientific_pilot_jobs=252,tight_ode_companions=4,
        reused_cases=324,available_scientific_cases_after_pilot=576,uniform_grid_cases_at_12_seeds=1416,
        uniform_grid_new_jobs_after_reuse=1092,
        scope="Broad positive activation strength only; mean and SD scale together; no sigma/correlation arm",
        controls="Uniform input equals declared ensemble mean; no path-mean normalization; all no-added-input controls retained",
        optional_range_extensions=[1e-5,10.],automatic_range_extension=False,
        replication="Pilot locates sampled change intervals; root/review freeze replication allocation before expansion. Uniform12 across all base strengths is a complete option, not yet submitted.")
    manifest=dict(generated_at=now(),experiment="Four-model broad input-strength exploratory pilot",status="frozen_not_submitted",
        plan=planvalues,basals=basals,jobs=jobs,source_sha256=sources,
        legacy={k:dict(run=str(p),manifest_sha256=sha(p/"manifest.json"),source_sha256=manifests[k]["source_sha256"]) for k,p in LEGACY.items()},
        native_protocols={
            "chemical":dict(noise=manifests["chemical"]["noise"],horizon=32000.,scalar_interval=5.,profile_interval=40.,full_interval=1000.,solver=manifests["chemical"]["solver"],
                            diagnostics=manifests["chemical"]["diagnostics"]),
            "spring":dict(noise=manifests["spring"]["noise"],horizon=3200.,scalar_interval=.5,profile_interval=4.,full_interval=100.,solver=manifests["spring"]["solver"])},
        numerical_check=dict(selected="Largest size, f=1, first new seed, each model",dense_factor=10,
            forcing="Same frozen generator, horizon, seed, knots, coupling, ell, tau and sigma; compare latent SHA exactly",
            tight_solver=dict(rtol=1e-10,atol=1e-12,max_step_factor=.5),
            gates=dict(occupancy_absolute_difference=.02,normalized_mean_absolute_difference=.02,
                       direction_correlation_absolute_difference=.05,direction_pair_coverage_absolute_difference=.02,
                       sampled_run_q90_relative_difference=.2,sampled_run_q90_absolute_floor_cadences=2),
            normalized_mean_tolerance="max(0.02,0.02*abs(reference_mean)); reference=standard for ODE, canonical coarse for observation check",
            windows="Inclusive full[0,T], first[0,T/2], last[T/2,T]; shared midpoint retained, windows not independent",
            lag_checkpoints=dict(chemical=[5.,40.,100.,500.,2000.,8000.],spring=[.5,4.,10.,50.,200.,800.]),
            direction_support="Both endpoint R>=0.1; at least100 eligible pairs and eligible fraction>=0.01 in both comparisons; otherwise unavailable",
            occupancy_gates=["raw current half-contrast multiregion","contrast>=0.1*basal_contrast AND multiregion","R<0.1","contrast<0.1*basal_contrast"],
            normalized_means=["R","M0/M0initial","q1/q1initial","contrast/basal_contrast"],
            morphology_sensitivity="Also export contrast floors0.01*basal and absolute0.05; prominence0.05 native and0.05*basal_contrast",
            pointwise_paths="Export as diagnostic; not sole acceptance criterion",undefined="Unavailable; no invented pass"),
        resources=dict(ordinary=dict(cpus=1,memory="3G",walltime="02:00:00",soft_seconds=6900),
                       dense=dict(cpus=1,memory="6G",walltime="04:00:00",soft_seconds=14100),qos="small",array_throttle=None),
        protocol_sha256=sha(run/"protocol/PLAN.md"),environment=dict(executable=PYTHON,python=sys.version,numpy=np.__version__))
    write(run/"manifest.json",manifest)
    csvwrite(run/"results/exact_reuse_cases.csv",reuse)
    csvwrite(run/"results/existing_runtime_storage.csv",runtimes)
    baseline_hours=0.;baseline_bytes=0.
    for b in basals:
        noise=[r for r in runtimes if r["model"]==b["model"] and r["size"]==b["size"] and r["kind"]=="noise" and r["fraction"]==.01]
        mean,=[r for r in runtimes if r["model"]==b["model"] and r["size"]==b["size"] and r["kind"]=="mean" and r["fraction"]==.01]
        baseline_hours+=7*(2*np.mean([r["elapsed_seconds"] for r in noise])+mean["elapsed_seconds"])/3600
        baseline_bytes+=7*3*np.mean([r["allocated_bytes"] for r in noise])
    write(run/"results/preparation.json",dict(generated_at=now(),manifest_sha256=sha(run/"manifest.json"),
        native_jobs_submitted=0,reused_cases=324,new_scientific_jobs=252,numerical_companions=4,
        baseline_252_worker_hours=float(baseline_hours),baseline_252_allocated_bytes=float(baseline_bytes),
        evidence="Extrapolation using f=.01 same-size old worker elapsed/allocated storage, not a measured prediction at larger f",
        exact_reuse_case_index_sha256=sha(run/"results/exact_reuse_cases.csv"),runtime_table_sha256=sha(run/"results/existing_runtime_storage.csv")))
    write(run/"status.json",dict(generated_at=now(),producer="frozen_not_submitted",complete=False,
        scientific_new_cases=252,numerical_companions=4,legacy_cases_available=324))
    return dict(run=str(run),manifest_sha256=sha(run/"manifest.json"),jobs=256,new_scientific_jobs=252,native_jobs_submitted=0)


def native_module(run,platform):
    if platform=="spring":
        # A single check process visits both frozen roots. Do not let its earlier
        # chemical models package choose the Spring native module by accident.
        for key in list(sys.modules):
            if key=="models" or key.startswith("models."):sys.modules.pop(key)
    path=run/"frozen_sources"/platform/"code/mechanisms"
    sys.path.insert(0,str(path))
    return importlib.import_module("vision_noise_20261017" if platform=="chemical" else "revision_spring_noise_20261019")


def check(run):
    run=Path(run).resolve();m=read(run/"manifest.json")
    for rel,digest in m["source_sha256"].items():assert sha(run/rel)==digest
    assert sha(run/"protocol/PLAN.md")==m["protocol_sha256"]
    for b in m["basals"]:
        assert sha(run/b["input"])==b["input_sha256"]
        with np.load(run/b["input"]) as z:assert ah(z["state"])==b["state_sha256"]
    assert len(m["jobs"])==256 and len({j["seed"] for j in m["jobs"] if j["kind"]=="noise" and j["variant"]=="standard"})==168
    for j in m["jobs"]:
        assert j["mean"]==REFERENCES[j["model"]]*j["fraction"] and j["fraction"] in GRID
        if j["variant"]=="tight_ode":
            p=m["jobs"][j["paired_standard_slot"]]
            assert all(j[k]==p[k] for k in ("seed","fraction","model","size","input_sha256","dense_observation"))
    # Algebraic input-channel and conservation checks, with NO ODE integration.
    errors=[]
    for platform in ("chemical","spring"):
        native=native_module(run,platform)
        for b in [v for v in m["basals"] if v["platform"]==platform]:
            with np.load(run/b["input"]) as z:y=z["state"].copy()
            if platform=="chemical":
                grid=native.ring.RingGrid(b["size"],b["cells"]);model=native.ring.native_model(b["model"])
                for f in GRID:
                    I=np.full(b["cells"],REFERENCES[b["model"]]*f)
                    diff=model.reaction(0.,grid.x,y,I)-model.reaction(0.,grid.x,y,np.zeros_like(I))
                    expected=I*y[1]/(7.5 if b["model"]=="holmes_model3" else 1.)
                    error=float(np.max(abs(diff[0]-expected)));assert error<1e-10
                    assert np.max(abs(diff[0]+diff[1]))<1e-10
                    errors.append(error)
            else:
                model=native.Closure(b["size"],b["cells"]);n=b["cells"]
                for f in GRID:
                    I=np.full(n,4*f);diff=model.rhs(y,I)-model.rhs(y,np.zeros(n))
                    expected=I*y[n:-1];error=float(np.max(abs(diff[:n]-expected)));assert error<1e-10
                    assert np.max(abs(diff[:n]+diff[n:-1]))<1e-10 and diff[-1]==0
                    errors.append(error)
    result=dict(generated_at=now(),passed=True,manifest_sha256=sha(run/"manifest.json"),native_integrations=0,
                frozen_source_files_checked=len(m["source_sha256"]),qualified_basals_checked=12,
                algebraic_model_size_strength_checks=len(errors),maximum_input_channel_error=max(errors),
                reuse_scope="All324 old cases retain original source/input/law/solver/cadence and output SHA proofs; old two strengths are not rerun",
                next="Independent protocol/source preflight receipt and root launch authorization required before submit")
    write(run/"results/local_preflight.json",result);return result


def canonicalize_dense(directory,factor,scalar_dt,profile_dt):
    for name,dt in (("trace",scalar_dt),("primary",profile_dt)):
        source=directory/(name+".npz");target=directory/("dense_"+name+".npz");source.rename(target)
        with np.load(target) as z:
            n=len(z["time"]);data={k:(z[k][::factor].copy() if z[k].ndim and z[k].shape[0]==n else z[k].copy()) for k in z.files}
        expected=np.arange(len(data["time"]))*dt
        np.testing.assert_allclose(data["time"],expected,rtol=0,atol=1e-9);data["time"]=expected
        np.savez_compressed(source,**data)


def worker(run,slot):
    assert os.environ.get("SLURM_JOB_ID"),"Native integration requires SLURM"
    run=Path(run).resolve();m=read(run/"manifest.json");j=m["jobs"][slot]
    for rel,digest in m["source_sha256"].items():assert sha(run/rel)==digest
    assert sha(run/"protocol/PLAN.md")==m["protocol_sha256"]
    assert sha(run/j["input"])==j["input_sha256"]
    directory=run/"workers"/f"{slot:04d}";directory.mkdir(exist_ok=False)
    status=dict(generated_at=now(),started_at=now(),producer="running",status="running",complete=False,
                job=j,slurm_job_id=os.environ["SLURM_JOB_ID"],manifest_sha256=sha(run/"manifest.json"))
    write(directory/"status.json",status);start=time.monotonic()
    native=None
    def stop(*args):
        globals()["STOP"]=True
        if native is not None:setattr(native,"STOP",True)
    signal.signal(signal.SIGUSR1,stop);signal.signal(signal.SIGTERM,stop)
    spec=m["native_protocols"][j["platform"]];dense=bool(j["dense_observation"]);factor=10 if dense else 1
    resource=m["resources"]["dense" if dense else "ordinary"];deadline=start+resource["soft_seconds"]
    try:
        native=native_module(run,j["platform"])
        with np.load(run/j["input"]) as z:initial=z["state"].copy()
        tight=j["variant"]=="tight_ode";solver=copy.deepcopy(spec["solver"])
        if tight:solver.update(rtol=1e-10,atol=1e-12,max_step=solver["max_step"]*.5)
        if j["platform"]=="chemical":
            grid=native.ring.RingGrid(j["size"],j["cells"]);model=native.ring.native_model(j["model"])
            forcing=native.NoiseInput(grid,spec["horizon"],j["seed"],j["mean"],kind=j["kind"],**spec["noise"])
            write(directory/"noise_audit.json",forcing.audit())
            modes=forcing.mode_records()
            with (directory/"noise_modes.csv").open("w",newline="") as f:
                w=csv.writer(f);w.writerow(modes);w.writerows(zip(*modes.values()))
            nm=dict(horizon=spec["horizon"],output_interval=spec["scalar_interval"]/factor,
                primary_interval=spec["profile_interval"]/factor,full_interval=spec["full_interval"],
                solver=solver,diagnostics=spec["diagnostics"])
            result=native.integrate(model,grid,initial,forcing,nm,deadline,directory,status)
            status.update(result)
            if not result["complete"]:raise RuntimeError(result["message"])
        else:
            from scipy.integrate import solve_ivp
            n=j["cells"];closure=native.Closure(j["size"],n);noise=spec["noise"]
            forcing=native.MaterialNoise(n,initial[-1],spec["horizon"],j["seed"],j["mean"],kind=j["kind"],
                knot=noise["knot"],coupling=noise["coupling"],tau=noise["tau"],ell=noise["ell"],sigma=noise["sigma"],bridge_base=noise.get("bridge_base"))
            write(directory/"noise_audit.json",forcing.audit())
            t=np.arange(int(round(spec["horizon"]/(spec["scalar_interval"]/factor)))+1)*(spec["scalar_interval"]/factor)
            def rhs(t,y):
                if STOP or time.monotonic()>=deadline:raise TimeoutError("SLURM signal or soft walltime limit")
                return closure.rhs(y,forcing(t))
            solution=solve_ivp(rhs,(0.,spec["horizon"]),initial,method="BDF",t_eval=t,
                rtol=solver["rtol"],atol=solver["atol"],max_step=solver["max_step"],
                jac=lambda t,y:closure.jacobian(y,forcing(t)))
            assert solution.success and len(solution.t)==len(t),solution.message
            y=solution.y.T;u,v,C=y[:,:n],y[:,n:-1],y[:,-1]
            assert np.isfinite(y).all() and y[:,:-1].min()>=-5e-10 and C.min()>0
            amount=C*np.mean(u+v,axis=1);error=float(np.max(abs(amount/amount[0]-1)));assert error<=1e-7
            q=closure.material_positions;M0=C*np.mean(u,axis=1);M1=C*(u@np.exp(2j*np.pi*q))/n;M2=C*(u@np.exp(4j*np.pi*q))/n;P=M1/M1[0]
            desc=native.diagnostics(u,q)
            np.savez_compressed(directory/"trace.npz",time=t,circumference=C,P_real=P.real,P_imag=P.imag,
                M1_real=M1.real,M1_imag=M1.imag,M0=M0,q1=abs(M1)/M0,q2=abs(M2)/M0,
                amount=amount,local_input_mean=np.array([forcing(tt).mean() for tt in t]),material_positions=q,**desc)
            pi=np.arange(0,len(t),8);fi=np.flatnonzero(np.isclose(t/spec["full_interval"],np.rint(t/spec["full_interval"])))
            np.savez_compressed(directory/"primary.npz",time=t[pi],primary=u[pi],inactive=v[pi],circumference=C[pi],material_positions=q)
            np.savez_compressed(directory/"full_states.npz",time=t[fi],states=y[fi],material_positions=q)
            np.savez_compressed(directory/"final_state.npz",state=y[-1],time=spec["horizon"])
            status.update(status="complete",native_time=spec["horizon"],all_finite=True,minimum_concentration=float(y[:,:-1].min()),
                max_amount_relative_error=error,circumference_range=[float(C.min()),float(C.max())],
                nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu)
        if dense:canonicalize_dense(directory,factor,spec["scalar_interval"],spec["profile_interval"])
        with np.load(directory/"trace.npz") as z:assert len(z["time"])==6401
        with np.load(directory/"primary.npz") as z:assert len(z["time"])==801
        status.update(status="complete",complete=True,canonical_scalar_samples=6401,canonical_profile_samples=801)
    except Exception as e:
        status.update(status="failed",complete=False,error=f"{type(e).__name__}: {e}",traceback=traceback.format_exc())
        raise
    finally:
        status.update(generated_at=now(),finished_at=now(),producer="terminal",elapsed_seconds=time.monotonic()-start,
                      output_sha256={p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!="status.json"})
        write(directory/"status.json",status)
    return dict(slot=slot,status=status["status"],elapsed_seconds=status["elapsed_seconds"])


def submit(run,review):
    run=Path(run).resolve();m=read(run/"manifest.json");r=read(review)
    assert r.get("passed") is True and r.get("manifest_sha256")==sha(run/"manifest.json"),"Manifest-specific independent preflight required"
    local=read(run/"results/local_preflight.json")
    assert local["passed"] and local["manifest_sha256"]==sha(run/"manifest.json") and not (run/"submission.json").exists()
    handles=[]
    for group,dense in (("ordinary",False),("dense",True)):
        prior=run/("submission_"+group+".json")
        if prior.exists():
            handles.append(read(prior));continue
        slots=[str(j["slot"]) for j in m["jobs"] if bool(j["dense_observation"])==dense]
        res=m["resources"][group]
        command="exec "+shlex.join([PYTHON,str(run/"frozen_sources/campaign.py"),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
        args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem="+res["memory"],"--time="+res["walltime"],
              "--signal=B:USR1@60","--array="+",".join(slots),"--job-name=f4-strength-"+group,
              "--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap",command]
        job=subprocess.check_output(args,text=True).strip().split(";")[0]
        record=dict(group=group,job_id=job,slots=len(slots),command=args,submitted_at=now())
        write(run/("submission_"+group+".json"),record);handles.append(record)
    result=dict(submitted_at=now(),manifest_sha256=sha(run/"manifest.json"),preflight_review=str(Path(review).resolve()),
                preflight_review_sha256=sha(review),arrays=handles)
    write(run/"submission.json",result);return result


def status(run):
    run=Path(run).resolve();m=read(run/"manifest.json");counts={};complete=0
    records=[]
    for j in m["jobs"]:
        p=run/"workers"/f"{j['slot']:04d}"/"status.json"
        s=read(p) if p.exists() else {}
        state=s.get("status","not_started");counts[state]=counts.get(state,0)+1
        complete+=bool(s.get("complete"))
        records.append(dict(slot=j["slot"],model=j["model"],size=j["size"],fraction=j["fraction"],kind=j["kind"],
                            variant=j["variant"],state=state,producer=s.get("producer","not_started"),
                            native_time=s.get("native_time"),elapsed_seconds=s.get("elapsed_seconds")))
    terminal=all(r["producer"]=="terminal" for r in records)
    z=dict(generated_at=now(),producer="terminal" if terminal else "submitted" if (run/"submission.json").exists() else "frozen_not_submitted",
           complete=complete==len(m["jobs"]),counts=counts,completed_cases=complete,requested_jobs=len(m["jobs"]),
           original_cases_retained=324,numerical_companions=4,manifest_sha256=sha(run/"manifest.json"))
    csvwrite(run/"results/live_worker_inventory.csv",records);write(run/"status.json",z);return z


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("command",choices=("prepare","check","worker","submit","status"))
    parser.add_argument("--run",type=Path,default=DEFAULT);parser.add_argument("--slot",type=int);parser.add_argument("--review",type=Path)
    a=parser.parse_args()
    if a.command=="prepare":result=prepare(a.run)
    elif a.command=="check":result=check(a.run)
    elif a.command=="worker":result=worker(a.run,a.slot)
    elif a.command=="status":result=status(a.run)
    else:assert a.review is not None;result=submit(a.run,a.review)
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
