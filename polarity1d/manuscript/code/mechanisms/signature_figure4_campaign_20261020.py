"""Prospective native Figure 4 ring experiments following the E0 interim.

The retained lognormal ring forcing and native solvers are copied unchanged.
This is separate from Figure 2's bounded, linear gap-noise experiment. Planning
and freezing perform no integration; production requires a delivered E0 interim
and explicit experiment choices, with SLURM workers only.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shlex
import shutil
import signal
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/"polarity1d/manuscript").is_dir())
MS=ROOT/"polarity1d/manuscript"
SOURCE=MS/"data/runs/revision_noise_strength_replication_scoped_20261019"
CANONICAL_CENSUS=ROOT/"docs/tasks/briefs/20261020_figure4_source_census.json"
PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS=("wave_pinning","goryachev","holmes_model3","spring_ring")
REFERENCES=dict(wave_pinning=.1,goryachev=6.4,holmes_model3=1.,spring_ring=4.)
CENTRAL_FIVE=(.001,.003,.01,.03,.1)
STOP=False
ODE_PREFIX_RULE=dict(version=1,relative_tolerance=2e-7,absolute_tolerance=2e-9,
    gaussian_prefix="Exact SHA256 of all original latent Gaussian knots.",
    initial_state="Exact qualified native initial-state/input hashes.",
    scalar_and_profile="Every original numeric field at every original saved time, with no field interpolation; all original fields must remain present.",
    coupled_state="Every species and Spring circumference at the original horizon, read from the exact shared full-state checkpoint time.",
    spatial_grid="Exact original x/material-position arrays.",
    interpretation="Numerical evolution equivalence, not bitwise identity at the original stopping cap. Missing or failed checks exclude long-arm qualification and never substitute the original shorter run.")


def now():return datetime.now(timezone.utc).isoformat()
def read(path):return json.loads(Path(path).read_text())
def sha(path):
    result=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""):result.update(block)
    return result.hexdigest()
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");temporary.replace(path)


def environment():
    import platform
    import numpy as np
    import scipy
    return dict(python=sys.version,executable=sys.executable,numpy=np.__version__,scipy=scipy.__version__,
        platform=platform.platform(),threads={name:os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")})


def sizes(model,intermediate=False):
    return ((1.5,2.5) if intermediate else (1.,2.,3.)) if model=="spring_ring" else (
        (30.,50.) if intermediate else (20.,40.,60.))


def protocol(platform,stage,native_protocols):
    value=copy.deepcopy(native_protocols[platform])
    if stage in ("E2","E4"):value["profile_interval"]=value["scalar_interval"]
    if stage=="E4":value["horizon"]*=4
    return value


def grid(stage,*,strengths=None,triggers=None,pilot_replicates=None):
    """Count noisy trajectories separately from new deterministic controls."""
    if stage=="E1_pilot":
        if set(strengths or {})!=set(MODELS):raise ValueError("Specify the separate resource/schema pilot strength per model")
        pilot=[]
        for model in MODELS:
            f=float(strengths[model]);replicate=(pilot_replicates or {}).get(model,0)
            if not math.isfinite(f) or f<=0 or not 0<=replicate<12:raise ValueError("Invalid prospective resource-pilot choice")
            for sigma in (1.,2.):
                pilot.append(dict(stage=stage,model=model,size=max(sizes(model)),platform="spring" if model=="spring_ring" else "chemical",
                    fraction=f,mean=REFERENCES[model]*f,sigma=sigma,kind="noise",replicate=replicate,
                    seed=202610204900000+MODELS.index(model)*100+replicate,slot=len(pilot)))
        return pilot
    if stage not in ("E1","E2","E3","E4"):raise ValueError("Unknown Figure 4 experiment")
    if stage=="E1":
        if set(strengths or {})!=set(MODELS):raise ValueError("E1 requires two explicit E0-informed strengths per model")
        for model,levels in strengths.items():
            if len(levels)!=2 or len(set(levels))!=2 or any(not math.isfinite(f) or f<=0 for f in levels):
                raise ValueError("E1 requires two distinct finite positive strengths; exact sigma=1 companion coverage is checked when freezing")
    groups=[]
    if stage=="E4":
        if not triggers:raise ValueError("E4 requires the actual E0 longest-usable-lag triggers")
        seen=set()
        for trigger in triggers:
            model=trigger["model"];size=float(trigger["size"]);f=float(trigger["fraction"])
            if model not in MODELS or size not in sizes(model):raise ValueError("E4 trigger has an unknown model/size")
            if (not 0<f or not math.isfinite(f) or not math.isfinite(trigger["heading_cosine"])
                    or not math.isfinite(trigger["longest_usable_lag"]) or trigger["heading_cosine"]<=.8 or trigger["longest_usable_lag"]<=0):
                raise ValueError("E4 requires positive f and cosine>0.8 at the longest usable lag")
            if (model,size,f) in seen:raise ValueError("Repeated E4 trigger")
            seen.add((model,size,f));groups.append((model,size,f,1.,trigger))
    else:
        for model in MODELS:
            levels=strengths[model] if stage=="E1" else (.005,.02,.05) if stage=="E2" else CENTRAL_FIVE
            for size in sizes(model,stage=="E3"):
                for f in levels:
                    for sigma in ((.25,.5,2.) if stage=="E1" else (1.,)):
                        groups.append((model,size,float(f),sigma,None))
    rows=[]
    for group_index,(model,size,f,sigma,trigger) in enumerate(groups):
        for replicate in range(12):
            rows.append(dict(stage=stage,model=model,size=size,platform="spring" if model=="spring_ring" else "chemical",
                fraction=f,mean=REFERENCES[model]*f,sigma=sigma,kind="noise",replicate=replicate,
                seed=202610204000000+int(stage[1])*100000000+(group_index//3 if stage=="E1" else group_index)*100+replicate,
                latent_pair_group=group_index//3 if stage=="E1" else None,
                trigger=trigger))
    if stage in ("E2","E3","E4"):
        for model,size,f in sorted({(m,s,f) for m,s,f,_,_ in groups}):
            rows.append(dict(stage=stage,model=model,size=size,platform="spring" if model=="spring_ring" else "chemical",
                fraction=f,mean=REFERENCES[model]*f,sigma=0.,kind="mean",replicate=0,seed=0))
    if stage in ("E3","E4"):
        for model,size in sorted({(m,s) for m,s,_,_,_ in groups}):
            rows.append(dict(stage=stage,model=model,size=size,platform="spring" if model=="spring_ring" else "chemical",
                fraction=0.,mean=0.,sigma=0.,kind="basal",replicate=0,seed=0))
    for slot,row in enumerate(rows):row["slot"]=slot
    return rows


def retained_cases(source):
    """Canonical case cache; duplicate physical outputs are not new seeds."""
    cache=read(source/"results/partial/case_cache.json");result=[];seen={}
    for case in cache.values():
        descriptor=case["descriptors"][0]
        key=tuple(descriptor[k] for k in ("model","size","fraction","kind","seed"))
        identity=(case["trace_sha256"],case["profile_sha256"])
        if key in seen:
            if seen[key]!=identity:raise ValueError("Conflicting retained outputs for one scientific seed")
            continue
        seen[key]=identity
        result.append({k:descriptor[k] for k in ("case_id","model","size","fraction","kind","seed")}|
            {k:case[k] for k in ("trace","trace_sha256","profile","profile_sha256")})
    return result


def select_conditions(jobs,conditions):
    """Split declared conditions only, preserving every sigma and seed cell."""
    if conditions is None:return jobs,[]
    keys=[(row["model"],float(row["size"]),float(row["fraction"])) for row in conditions]
    available={(j["model"],j["size"],j["fraction"]) for j in jobs}
    if not keys or len(set(keys))!=len(keys) or not set(keys)<=available:
        raise ValueError("A condition subset must contain unique declared model/size/f conditions")
    selected=[];deferred=[]
    for job in jobs:
        row=copy.deepcopy(job);row["full_stage_slot"]=job["slot"]
        (selected if (job["model"],job["size"],job["fraction"]) in keys else deferred).append(row)
    for slot,row in enumerate(selected):row["slot"]=slot
    return selected,deferred


def reused_controls_and_companions(stage,jobs,retained):
    needed=[]
    if stage in ("E1","E1_pilot"):
        for model,size,f in sorted({(j["model"],j["size"],j["fraction"]) for j in jobs}):
            needed.append((model,size,f,"mean"))
            companions=[r for r in retained if (r["model"],r["size"],r["fraction"],r["kind"])==(model,size,f,"noise")]
            if len(companions)!=12:raise ValueError("E1 lacks twelve retained sigma=1 companions")
            needed.extend((model,size,f,"noise",r["seed"]) for r in companions)
    if stage in ("E1","E1_pilot","E2"):
        needed.extend((model,size,0.,"basal") for model,size in sorted({(j["model"],j["size"]) for j in jobs}))
    records=[]
    for requirement in needed:
        selected=[r for r in retained if tuple(r[k] for k in ("model","size","fraction","kind"))==requirement[:4]
                  and (len(requirement)==4 or r["seed"]==requirement[4])]
        if len(selected)!=1:raise ValueError("Missing or duplicated exact retained control/companion")
        record=copy.deepcopy(selected[0])
        for field in ("trace","profile"):
            if sha(record[field])!=record[field+"_sha256"]:raise ValueError("Retained evidence changed")
        records.append(record)
    return records


def prepare(out,stage,gate,*,strengths=None,triggers=None,extra_basals=None,companion_census=None,condition_subset=None,source=SOURCE):
    out=Path(out).resolve();source=Path(source).resolve();gate=Path(gate).resolve()
    if out.exists() or not out.is_relative_to((MS/"data/runs").resolve()):raise ValueError("Choose a fresh manuscript data/runs directory")
    release_digest=sha(gate);release=read(gate)
    if release.get("e0_interim_delivered") is not True and not (stage=="E1_pilot" and release.get("pilot_preparation_authorized") is True):
        raise ValueError("Deliver E0 A1–A4/A7 interim before freezing new native experiments")
    if release.get("e0_interim_delivered") is True:
        if not release.get("interim_artifacts"):raise ValueError("Identify concrete delivered E0 interim artifacts before freezing native production")
        for artifact,digest in release["interim_artifacts"].items():
            if sha(artifact)!=digest:raise ValueError("Delivered E0 interim artifact changed")
    if stage=="E1_pilot" and release.get("pilot_strengths")!=strengths:raise ValueError("Resource/schema pilot strengths must match their distinct preparation record")
    if stage=="E1" and (release.get("e1_strengths_fixed") is not True or release.get("e1_strengths")!=strengths):
        raise ValueError("E1 strengths must match the fixed post-interim design receipt")
    if stage in ("E3","E4") and release.get("e1_e2_complete") is not True:raise ValueError("E3/E4 follow completed E1/E2")
    jobs=grid(stage,strengths=strengths,triggers=triggers,pilot_replicates=release.get("pilot_replicates"))
    full_stage_jobs=len(jobs)
    if condition_subset is not None and stage!="E1":raise ValueError("Only E1 may split pending sigma1 companion conditions")
    jobs,deferred=select_conditions(jobs,condition_subset);parent=read(source/"manifest.json")
    basals={(b["model"],float(b["size"])):dict(b,source=str(source/b["input"])) for b in parent["basals"]}
    for record in extra_basals or []:
        if record.get("qualified") is not True:raise ValueError("An intermediate basal is unqualified")
        basals[(record["model"],float(record["size"]))]=record
    required={(j["model"],j["size"]) for j in jobs}
    if not required<=set(basals):raise ValueError("Fresh qualified intermediate basals are required before E3")
    retained=retained_cases(source)
    if companion_census is None and source==SOURCE.resolve() and CANONICAL_CENSUS.exists():
        companion_census=CANONICAL_CENSUS
    if companion_census:
        companion_digest=sha(companion_census);supplied=read(companion_census)
        if isinstance(supplied,dict) and "source_manifest_sha256" in supplied and supplied["source_manifest_sha256"]!=sha(source/"manifest.json"):
            raise ValueError("Companion census points to another retained source campaign")
        extra=supplied.get("cases",[]) if isinstance(supplied,dict) else supplied
        identities={tuple(r[k] for k in ("model","size","fraction","kind","seed")):(r["trace_sha256"],r["profile_sha256"]) for r in retained}
        for record in extra:
            if record.get("sigma",1.)!=1. and record["kind"]=="noise":continue
            if record.get("horizon",3200. if record["model"]=="spring_ring" else 32000.)!=(3200. if record["model"]=="spring_ring" else 32000.):continue
            key=tuple(record[k] for k in ("model","size","fraction","kind","seed"));identity=(record["trace_sha256"],record["profile_sha256"])
            if key in identities:
                if identities[key]!=identity and not record.get("observation_companion_of"):raise ValueError("Conflicting sigma=1 companion identities")
            else:retained.append(record);identities[key]=identity
    reuses=reused_controls_and_companions(stage,jobs,retained)
    for job in jobs:
        job["independent_noise_seed"]=job["kind"]=="noise" and stage!="E1_pilot"
        job["new_gaussian_path"]=job["kind"]=="noise"
        if job["kind"] in ("mean","basal"):
            old=[r for r in retained if (r["model"],r["size"],r["fraction"],r["kind"])==
                 (job["model"],job["size"],job["fraction"],job["kind"])]
            if len(old)==1:job["observation_companion_of"]=dict(old[0],
                role="extended_horizon_control" if stage=="E4" else "denser_saved_profile_control")
    if stage in ("E1","E1_pilot","E2"):
        for job in jobs:
            if job["kind"]!="noise":continue
            existing=sorted([r for r in retained if (r["model"],r["size"],r["fraction"],r["kind"])==
                (job["model"],job["size"],job["fraction"],"noise")],key=lambda r:r["seed"])
            if stage=="E2" and job["replicate"]>=len(existing):continue
            old=existing[job["replicate"]];audit=Path(old["trace"]).parent/"noise_audit.json"
            status_path=Path(old["trace"]).parent/"status.json";old_status=read(status_path)
            old_basal_sha=old_status["job"].get("input_sha256") or old_status["job"].get("sha256")
            if old_status.get("complete") is not True or old_basal_sha!=basals[(job["model"],job["size"])]["input_sha256"]:
                raise ValueError("The paired sigma=1 path must use the identical qualified basal")
            job.update(seed=old["seed"],new_gaussian_path=False,paired_sigma1_source=dict(noise_audit=str(audit),noise_audit_sha256=sha(audit),
                gaussian_knots_sha256=read(audit)["gaussian_knots_sha256"],original_case=old,
                original_status=str(status_path),original_status_sha256=sha(status_path),exact_basal_sha256=old_basal_sha))
            if stage=="E2":job["observation_companion_of"]=dict(old,role="denser_saved_profile_noise_replay")
    if stage=="E4":
        for model,size,f in sorted({(j["model"],j["size"],j["fraction"]) for j in jobs if j["kind"]=="noise"}):
            existing=sorted([r for r in retained if (r["model"],r["size"],r["fraction"],r["kind"])==(model,size,f,"noise")],key=lambda r:r["seed"])
            subset=[j for j in jobs if (j["model"],j["size"],j["fraction"],j["kind"])==(model,size,f,"noise")]
            if len(existing)>12:raise ValueError("More than twelve retained canonical seeds for an E4 condition")
            for job,old in zip(subset,existing):
                job.update(seed=old["seed"],new_gaussian_path=False)
                audit=Path(old["trace"]).parent/"noise_audit.json"
                status_path=Path(old["trace"]).parent/"status.json";old_status=read(status_path)
                old_basal_sha=old_status["job"].get("input_sha256") or old_status["job"].get("sha256")
                if old_status.get("complete") is not True or old_basal_sha!=basals[(model,size)]["input_sha256"]:
                    raise ValueError("An E4 extension must retain the exact original basal")
                job["retained_prefix"]=dict(noise_audit=str(audit),noise_audit_sha256=sha(audit),
                    gaussian_knots_sha256=read(audit)["gaussian_knots_sha256"],
                    horizon=3200. if model=="spring_ring" else 32000.,original_case=old,
                    original_status=str(status_path),original_status_sha256=sha(status_path),exact_basal_sha256=old_basal_sha)
    for directory in ("inputs","frozen_sources","workers","logs","results"):(out/directory).mkdir(parents=True)
    sources={};inputs={};frozen_basals=[]
    for relative,digest in parent["source_sha256"].items():
        original=source/relative
        if sha(original)!=digest:raise ValueError("Retained native source changed")
        target=out/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(original,target)
        sources[relative]=digest
    target=out/"frozen_sources/figure4_campaign_20261020.py";shutil.copy2(__file__,target)
    sources[str(target.relative_to(out))]=sha(target)
    for key in sorted(required):
        basal=copy.deepcopy(basals[key]);original=Path(basal["source"])
        if sha(original)!=basal["input_sha256"]:raise ValueError("Exact qualified basal changed")
        relative=f"inputs/{basal['model']}_C{basal['size']:g}.npz";shutil.copy2(original,out/relative)
        if sha(out/relative)!=basal["input_sha256"]:raise ValueError("Basal changed while freezing")
        basal.update(input=relative);inputs[relative]=sha(out/relative);frozen_basals.append(basal)
    for job in jobs:
        basal=next(b for b in frozen_basals if (b["model"],float(b["size"]))==(job["model"],job["size"]))
        job.update({k:basal[k] for k in ("cells","input","input_sha256","state_sha256","basal_circumference")})
        job["protocol"]=protocol(job["platform"],stage,parent["native_protocols"])
    shutil.copy2(gate,out/"inputs/e0_interim_gate.json");inputs["inputs/e0_interim_gate.json"]=sha(out/"inputs/e0_interim_gate.json")
    if inputs["inputs/e0_interim_gate.json"]!=release_digest or sha(gate)!=release_digest:raise ValueError("Interim/design receipt changed while freezing")
    if companion_census:
        shutil.copy2(companion_census,out/"inputs/companion_census.json");inputs["inputs/companion_census.json"]=sha(out/"inputs/companion_census.json")
        if inputs["inputs/companion_census.json"]!=companion_digest or sha(companion_census)!=companion_digest:raise ValueError("Companion census changed while freezing")
    manifest=dict(generated_utc=now(),stage=stage,status="frozen_not_submitted",jobs=jobs,basals=frozen_basals,
        native_protocols=parent["native_protocols"],source_sha256=sources,input_sha256=inputs,
        parent_source=str(source),parent_manifest_sha256=sha(source/"manifest.json"),reused_cases=reuses,
        full_stage_planned_jobs=full_stage_jobs,selected_condition_subset=condition_subset,
        deferred_planned_cells=deferred,
        subset_scope="Prospectively declared E1 model/size/f conditions, retaining all sigma and seed cells; deferred cells require a separate exact companion-backed run and are not failures or censoring.",
        counts={kind:sum(j["kind"]==kind for j in jobs) for kind in ("noise","mean","basal")},
        central_five=list(CENTRAL_FIVE),strength_choices=strengths,e4_triggers=triggers,
        retained_native_prefix_equivalence_rule=ODE_PREFIX_RULE if stage=="E4" else None,
        forcing="Retained lognormal mean*exp(sigma*Z-sigma^2/2); no clipping or realized normalization.",
        physical_arc="Chemical x; Spring actual circumference(t)*material_positions, never resting circumference.",
        spring_forcing_frame="Retained material field with basal physical ell=.15; instantaneous physical ell scales as C(t)/Cbasal.",
        resources=dict(cpus=1,qos="small",memory="8G" if stage=="E4" else "4G",walltime="08:00:00" if stage=="E4" else "02:00:00"),
        replication_unit="Twelve independent seed blocks per condition; E1 sigma values share each retained sigma=1 Gaussian field and must be compared with paired whole-seed resampling.",
        native_production_held_for_E0_interim=release.get("e0_interim_delivered") is not True,
        environment=environment(),
        full_brief_complete=False,automatic_submission=False)
    write(out/"manifest.json",manifest);write(out/"freeze.json",dict(manifest_sha256=sha(out/"manifest.json")))
    verify(out)
    native_preflight(out)
    worker=out/"worker.sh";worker.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec "+shlex.join([PYTHON,"-B",str(target),"worker","--out",str(out)])+" --slot \"$SLURM_ARRAY_TASK_ID\"\n");worker.chmod(0o750)
    args=["sbatch","--parsable","--partition=red","--qos=small","--cpus-per-task=1","--mem="+manifest["resources"]["memory"],"--time="+manifest["resources"]["walltime"],f"--array=0-{len(jobs)-1}","--job-name=F4-"+stage,"--output="+str(out/"logs/%A_%a.out"),"--error="+str(out/"logs/%A_%a.err"),str(worker)]
    (out/"submission_command.txt").write_text(shlex.join(args)+"\n")
    return dict(output=str(out),manifest_sha256=sha(out/"manifest.json"),counts=manifest["counts"],submitted=False)


def verify(out):
    manifest=read(out/"manifest.json")
    if read(out/"freeze.json")["manifest_sha256"]!=sha(out/"manifest.json"):raise ValueError("Frozen manifest changed")
    for mapping in (manifest["source_sha256"],manifest["input_sha256"]):
        for relative,digest in mapping.items():
            if sha(out/relative)!=digest:raise ValueError("Frozen source/input changed: "+relative)
    return manifest


def native_module(out,platform):
    # A preflight can visit several frozen roots and both native platforms in
    # one process. Import only this payload's source, including the core package
    # named polarity1d, which otherwise collides with the repository namespace.
    names=("models","polarity1d","vision_noise_20261017","ring_distance_prototype_20261017",
        "revision_spring_noise_20261019","spring_ring_model_20261019","length_dev_output_layout")
    for name in list(sys.modules):
        if any(name==prefix or name.startswith(prefix+".") for prefix in names):sys.modules.pop(name)
    path=out/"frozen_sources/campaign.py"
    spec=importlib.util.spec_from_file_location("unchanged_f4_native_dispatch",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.native_module(out,platform)


def native_preflight(out):
    """Check exact input-channel algebra without evolving any model."""
    import numpy as np
    manifest=verify(out);rows=[]
    for platform in ("chemical","spring"):
        records=[b for b in manifest["basals"] if b["model"] in MODELS and (b["model"]=="spring_ring")== (platform=="spring")]
        if not records:continue
        native=native_module(out,platform)
        for basal in records:
            with np.load(out/basal["input"]) as data:state=data["state"].copy()
            state_digest=hashlib.sha256(np.asarray(state,dtype="<f8",order="C").tobytes()).hexdigest()
            if state_digest!=basal["state_sha256"] or not np.isfinite(state).all():raise ValueError("Frozen basal state changed")
            n=basal["cells"];mean=REFERENCES[basal["model"]]*.01
            if platform=="chemical":
                grid=native.ring.RingGrid(basal["size"],n);model=native.ring.native_model(basal["model"])
                command=mean*(1+.2*np.sin(2*np.pi*grid.x/basal["size"]))
                difference=model.reaction(0.,grid.x,state,command)-model.reaction(0.,grid.x,state,np.zeros(n))
                expected=command*state[1]/(7.5 if basal["model"]=="holmes_model3" else 1.)
                active=state[0];rhs=native.ring.native_rhs(model,grid,state,np.zeros(n))
                regions=len(native.ring.circular_components(active>=active.min()+.5*np.ptp(active)))
                error=float(np.max(abs(difference[0]-expected)));pool_error=float(np.max(abs(difference[0]+difference[1])))
            else:
                closure=native.Closure(basal["size"],n);active=state[:n]
                command=mean*(1+.2*np.sin(2*np.pi*closure.material_positions))
                difference=closure.rhs(state,command)-closure.rhs(state,np.zeros(n))
                error=float(np.max(abs(difference[:n]-command*state[n:-1])))
                pool_error=float(np.max(abs(difference[:n]+difference[n:-1])))
                if difference[-1]!=0 or not np.isclose(state[-1],basal["basal_circumference"],atol=1e-12,rtol=0):raise ValueError("Spring geometry/input channel changed")
                rhs=closure.rhs(state,np.zeros(n));regions=int(native.diagnostics(active[None,:],closure.material_positions)["current_half_contrast_regions"][0])
            if error>1e-10 or pool_error>1e-10 or np.ptp(active)<.05 or regions!=1:raise ValueError("Native basal/input-channel preflight failed")
            rows.append(dict(model=basal["model"],size=basal["size"],input_sha256=sha(out/basal["input"]),
                state_sha256=state_digest,input_channel_error=error,input_pool_error=pool_error,regions=regions,
                no_input_rhs_maximum_absolute=float(np.max(abs(rhs))),loaded_native_module=str(native.__file__)))
    report=dict(generated_utc=now(),passed=True,manifest_sha256=sha(out/"manifest.json"),native_integrations=0,basals=rows,
        scope="Exact source/input identity, single-region basal topology and nonuniform native input-channel algebra; original qualification/holdouts remain separate evidence.")
    write(out/"results/native_preflight.json",report);return report


class QueriedInput:
    """Observe actual input evaluations without changing the retained field."""
    def __init__(self,forcing):
        self.forcing=forcing;self.calls=0;self.count=0;self.total=0.
        self.minimum=float("inf");self.maximum=float("-inf")
    def __call__(self,t):
        import numpy as np
        values=self.forcing(t)
        if not np.isfinite(values).all():raise FloatingPointError("Nonfinite retained lognormal input")
        self.calls+=1;self.count+=values.size;self.total+=float(values.sum())
        self.minimum=min(self.minimum,float(values.min()));self.maximum=max(self.maximum,float(values.max()))
        return values
    def report(self,reference):
        return dict(evaluations=self.calls,values=self.count,
            minimum=self.minimum if self.count else None,maximum=self.maximum if self.count else None,
            maximum_over_native_reference=self.maximum/reference if self.count else None,
            evaluation_weighted_mean=self.total/self.count if self.count else None,
            scope="All native/diagnostic queries, including repeated implicit-solver queries; this mean is not a time average.",
            clipping_or_pathwise_normalization=False)


def input_statistics(forcing,job,spec):
    """Retain finite-realization input moments, including the response window."""
    import numpy as np
    audit=forcing.audit();knots=getattr(forcing,"knots",None)
    audit.update(sigma=job["sigma"],theoretical_input_mean=job["mean"],
        theoretical_input_variance=job["mean"]**2*math.expm1(job["sigma"]**2) if job["kind"]=="noise" else 0.,
        theoretical_native_reference_fraction=job["fraction"],finite_run_mean_is_not_normalized=True)
    windows={"full":(0.,spec["horizon"]),"postburn":(.5*spec["horizon"],spec["horizon"]),
        "postburn_first_half":(.5*spec["horizon"],.75*spec["horizon"]),
        "postburn_second_half":(.75*spec["horizon"],spec["horizon"]),
        "postburn_late_quarter":(.875*spec["horizon"],spec["horizon"])}
    sampled={};knot=spec["noise"]["knot"]
    for name,(start,stop) in windows.items():
        if knots is None:
            value=0. if job["kind"]=="basal" else job["mean"]
            sampled[name]=dict(time_start=start,time_stop=stop,mean=value,variance=0.,minimum=value,maximum=value,
                quantiles={str(q):value for q in (.01,.5,.9,.99,.999)})
        else:
            first=int(math.ceil(start/knot));last=int(math.floor(stop/knot))+1
            values=job["mean"]*np.exp(job["sigma"]*knots[first:last]-.5*job["sigma"]**2)
            if not np.isfinite(values).all():raise FloatingPointError("Nonfinite sampled lognormal input")
            sampled[name]=dict(time_start=first*knot,time_stop=(last-1)*knot,mean=float(values.mean()),
                variance=float(values.var()),minimum=float(values.min()),maximum=float(values.max()))
            if name in ("full","postburn"):
                sampled[name]["quantiles"]={str(q):float(v) for q,v in zip((.01,.5,.9,.99,.999),np.quantile(values,[.01,.5,.9,.99,.999]))}
    audit["realized_input_windows"]=sampled
    audit["realized_input_sampling_scope"]="Uniform latent-knot times and native cells; knot extrema/quantiles are sampled, not a continuous-time supremum. Actual solver-query extrema are reported separately."
    return audit


def check_paired_sigma1_field(forcing,job):
    if "paired_sigma1_source" not in job:return
    import numpy as np
    source=job["paired_sigma1_source"]
    if sha(source["noise_audit"])!=source["noise_audit_sha256"]:raise ValueError("Paired sigma=1 forcing receipt changed")
    if sha(source["original_status"])!=source["original_status_sha256"]:raise ValueError("Paired sigma=1 native receipt changed")
    measured=hashlib.sha256(np.asarray(forcing.knots,dtype="<f8",order="C").tobytes()).hexdigest()
    if measured!=source["gaussian_knots_sha256"]:raise ValueError("E1 did not retain its sigma=1 latent Gaussian field")


def forcing_preflight(out,shard=None,shards=None):
    """Regenerate frozen inputs and their distribution receipts, without ODEs."""
    import numpy as np
    start=time.monotonic();manifest=verify(out);rows=[];paired={}
    if (shard is None)!=(shards is None):raise ValueError("Declare both forcing-preflight shard and shard count")
    if shards is not None and (not 1<=shards<=len(manifest["jobs"]) or not 0<=shard<shards):raise ValueError("Invalid forcing-preflight shard")
    selected=manifest["jobs"] if shard is None else manifest["jobs"][shard::shards]
    for platform in ("chemical","spring"):
        jobs=[j for j in selected if j["platform"]==platform]
        if not jobs:continue
        native=native_module(out,platform)
        for job in jobs:
            spec=job["protocol"];noise=dict(spec["noise"],sigma=job["sigma"])
            references=check_retained_reference_artifacts(job)
            if platform=="chemical":
                field=native.NoiseInput(native.ring.RingGrid(job["size"],job["cells"]),
                    spec["horizon"],job["seed"],job["mean"],kind=job["kind"],**noise)
            else:
                field=native.MaterialNoise(job["cells"],job["basal_circumference"],spec["horizon"],job["seed"],job["mean"],
                    kind=job["kind"],knot=noise["knot"],coupling=noise["coupling"],tau=noise["tau"],ell=noise["ell"],
                    sigma=noise["sigma"],bridge_base=noise.get("bridge_base"))
            check_paired_sigma1_field(field,job);check_retained_prefix(field,job,noise["knot"])
            audit=input_statistics(field,job,spec)
            if getattr(field,"knots",None) is not None:
                digest=hashlib.sha256(np.asarray(field.knots,dtype="<f8",order="C").tobytes()).hexdigest()
                key=(job["model"],job["size"],job["fraction"],job["seed"])
                if key in paired and paired[key]!=digest:raise ValueError("Sigma companions have different latent paths")
                paired[key]=digest
            row=dict(slot=job["slot"],model=job["model"],size=job["size"],fraction=job["fraction"],sigma=job["sigma"],
                seed=job["seed"],input_sha256=job["input_sha256"],native_input_audit=audit,retained_native_reference_artifacts=references)
            path=out/"results/forcing_preflight"/f'{job["slot"]:04d}.json';write(path,row)
            rows.append({k:row[k] for k in ("slot","model","size","fraction","sigma","seed")}|dict(
                receipt=str(path),receipt_sha256=sha(path),sampled_mean=audit["realized_input_windows"]["full"]["mean"],
                sampled_maximum=audit["realized_input_windows"]["full"]["maximum"]))
            del field
    report=dict(generated_utc=now(),manifest_sha256=sha(out/"manifest.json"),passed=True,native_integrations=0,
        production_held_for_E0_interim=manifest.get("native_production_held_for_E0_interim",False),workers=rows,
        selected_tasks=[j["slot"] for j in selected],shard=shard,shards=shards,full_manifest_workers=len(manifest["jobs"]),
        wall_seconds=time.monotonic()-start,peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope="Exact frozen native forcing and latent-path checks only; no evolution or response/hypothesis analysis.")
    path=out/"results/forcing_preflight.json" if shard is None else out/"results/forcing_preflight_shards"/f"shard_{shard:03d}.json"
    write(path,report);return report


def collect_forcing_preflight(out,shards):
    manifest=verify(out);covered=[];rows=[];missing=[];failed=[];proofs=[]
    if not 1<=shards<=len(manifest["jobs"]):raise ValueError("Invalid forcing-preflight shard count")
    for shard in range(shards):
        path=out/"results/forcing_preflight_shards"/f"shard_{shard:03d}.json"
        if not path.exists():missing.append(shard);continue
        record=read(path);expected=[j["slot"] for j in manifest["jobs"][shard::shards]]
        if record.get("passed") is not True:failed.append(shard);continue
        if record["manifest_sha256"]!=sha(out/"manifest.json") or record["selected_tasks"]!=expected or record["shards"]!=shards or record["shard"]!=shard:
            raise ValueError("Forcing-preflight shard coverage/identity changed")
        if len(record["workers"])!=len(expected) or {r["slot"] for r in record["workers"]}!=set(expected):raise ValueError("Forcing-preflight shard lost workers")
        for row in record["workers"]:
            proof=Path(row["receipt"])
            if sha(proof)!=row["receipt_sha256"]:raise ValueError("Forcing-preflight worker receipt changed")
            detailed=read(proof);job=manifest["jobs"][row["slot"]]
            if any(detailed[key]!=job[key] for key in ("slot","model","size","fraction","sigma","seed","input_sha256")):
                raise ValueError("Forcing-preflight worker identity changed")
        covered.extend(expected);rows.extend(record["workers"]);proofs.append(dict(path=str(path),sha256=sha(path)))
    complete=not missing and not failed and len(covered)==len(set(covered))==len(manifest["jobs"])
    report=dict(generated_utc=now(),manifest_sha256=sha(out/"manifest.json"),passed=complete,native_integrations=0,
        missing_shards=missing,failed_shards=failed,shards=shards,workers=sorted(rows,key=lambda r:r["slot"]),shard_receipts=proofs,
        selected_tasks=sorted(covered),full_manifest_workers=len(manifest["jobs"]),
        scope="All frozen initial inputs and native Gaussian fields/prefixes qualified before evolution; no response or scientific hypothesis is evaluated.")
    write(out/"results/forcing_preflight.json",report);return report


def compare_sigma1_evolution(directory,job):
    """Independent σ=1 pilot compares every common scalar/profile/state field."""
    if job["sigma"]!=1. or "paired_sigma1_source" not in job:return
    import numpy as np
    old=job["paired_sigma1_source"]["original_case"]
    checks=[]
    for name,source in (("trace.npz",Path(old["trace"])),("primary.npz",Path(old["profile"])),
                        ("final_state.npz",Path(old["trace"]).parent/"final_state.npz")):
        with np.load(source) as original,np.load(directory/name) as fresh:
            source_index=fresh_index=None
            if name=="primary.npz" and not np.array_equal(original["time"],fresh["time"]):
                # Dense E2 observations replay the same native evolution. Check
                # every coarser saved frame exactly, without interpolating fields.
                coarse,fine=(original["time"],fresh["time"]) if len(original["time"])<=len(fresh["time"]) else (fresh["time"],original["time"])
                indices=np.searchsorted(fine,coarse)
                if np.any(indices>=len(fine)):raise ValueError("Sigma=1 profile clocks have different coverage")
                np.testing.assert_allclose(fine[indices],coarse,rtol=0,atol=1e-9)
                if len(original["time"])<=len(fresh["time"]):fresh_index=indices
                else:source_index=indices
            required={field for field in original.files if original[field].dtype.kind in "biufc"}
            if not required<=set(fresh.files):raise ValueError("Sigma=1 companion dropped saved scientific fields: "+name)
            for field in sorted(set(original.files)&set(fresh.files)):
                a,b=original[field],fresh[field]
                if a.dtype.kind not in "biufc" or b.dtype.kind not in "biufc":continue
                if field not in ("x","material_positions"):
                    if source_index is not None and a.ndim and a.shape[0]==len(original["time"]):a=a[source_index]
                    if fresh_index is not None and b.ndim and b.shape[0]==len(fresh["time"]):b=b[fresh_index]
                np.testing.assert_allclose(a,b,rtol=2e-7,atol=2e-9,equal_nan=True,err_msg=name+":"+field)
                checks.append(dict(archive=name,field=field,shape=list(a.shape),bitwise_equal=bool(np.array_equal(a,b,equal_nan=True)),
                    maximum_absolute_difference=float(np.max(abs(a-b))) if a.dtype.kind not in "b" and a.size else 0.))
    write(directory/"sigma1_evolution_equivalence.json",dict(passed=True,checks=checks,
        original_case=old,rtol=2e-7,atol=2e-9,scope="Full saved scalar history, profile history and final coupled state; exact original Gaussian field already checked."))


def check_retained_evolution_prefix(directory,job):
    """Audit a long solve against all original saved native observations."""
    if "retained_prefix" not in job:return
    import numpy as np
    prefix=job["retained_prefix"];old=prefix["original_case"];checks=[]
    if sha(prefix["original_status"])!=prefix["original_status_sha256"]:raise ValueError("Original native-prefix status changed")
    rtol=ODE_PREFIX_RULE["relative_tolerance"];atol=ODE_PREFIX_RULE["absolute_tolerance"]
    for archive,key in (("trace.npz","trace"),("primary.npz","profile")):
        path=Path(old[key])
        if sha(path)!=old[key+"_sha256"]:raise ValueError("Original native-prefix observation changed")
        with np.load(path) as source,np.load(directory/archive) as fresh:
            indices=np.searchsorted(fresh["time"],source["time"])
            if np.any(indices>=len(fresh["time"])):raise ValueError("Extended observations do not cover their complete original prefix")
            np.testing.assert_array_equal(fresh["time"][indices],source["time"])
            required={field for field in source.files if source[field].dtype.kind in "biufc"}
            if not required<=set(fresh.files):raise ValueError("Extended run dropped original numeric fields")
            for field in sorted(required):
                a,b=source[field],fresh[field]
                if field in ("x","material_positions"):np.testing.assert_array_equal(a,b)
                else:
                    if a.ndim and a.shape[0]==len(source["time"]):
                        if not b.ndim or b.shape[0]!=len(fresh["time"]):raise ValueError("Extended native observation layout changed")
                        b=b[indices]
                    np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,equal_nan=True,err_msg=archive+":"+field)
                checks.append(dict(archive=archive,field=field,shape=list(a.shape),bitwise_equal=bool(np.array_equal(a,b,equal_nan=True)),
                    maximum_absolute_difference=float(np.max(abs(a-b))) if a.dtype.kind!="b" and a.size else 0.))
    original_final=Path(old["trace"]).parent/"final_state.npz";status=read(prefix["original_status"])
    if sha(original_final)!=status["output_sha256"]["final_state.npz"]:raise ValueError("Original final coupled state changed")
    with np.load(original_final) as original,np.load(directory/"full_states.npz") as fresh:
        if float(original["time"])!=prefix["horizon"]:raise ValueError("Original native checkpoint does not reach its declared cap")
        indices=np.flatnonzero(fresh["time"]==prefix["horizon"])
        if len(indices)!=1:raise ValueError("Extended run needs one exact original-cap full-state checkpoint")
        a,b=original["state"],fresh["states"][indices[0]]
        np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,equal_nan=True)
        checks.append(dict(archive="full_states.npz",field="state_at_original_horizon",shape=list(a.shape),
            bitwise_equal=bool(np.array_equal(a,b,equal_nan=True)),maximum_absolute_difference=float(np.max(abs(a-b)))))
    write(directory/"retained_native_prefix_equivalence.json",dict(passed=True,rule=ODE_PREFIX_RULE,
        original_horizon=prefix["horizon"],original_case=old,checks=checks))


def spring_native_solve(rhs,jacobian,initial,times,solver,checkpoint=None):
    """Keep the native BDF settings, optionally restarting at an original cap.

    A longer t_span can change dense-output values at the original terminal
    step. A separately declared restart reproduces that terminal step before
    continuing with the same global-time RHS and forcing. It changes neither
    the retained tolerance nor the native equations; the original attempt is
    retained separately and every saved prefix observation is checked again.
    """
    import numpy as np
    from scipy.integrate import solve_ivp
    times=np.asarray(times)
    kwargs=dict(method="BDF",rtol=solver["rtol"],atol=solver["atol"],max_step=solver["max_step"],jac=jacobian)
    if checkpoint is None:
        return solve_ivp(rhs,(float(times[0]),float(times[-1])),initial,t_eval=times,**kwargs)
    checkpoint=float(checkpoint)
    matches=np.flatnonzero(times==checkpoint)
    if len(matches)!=1 or not times[0]<checkpoint<times[-1]:
        raise ValueError("Native Spring restart needs one exact internal observation checkpoint")
    cut=int(matches[0]);first=solve_ivp(rhs,(float(times[0]),checkpoint),initial,t_eval=times[:cut+1],**kwargs)
    if not first.success or len(first.t)!=cut+1:return first
    last=solve_ivp(rhs,(checkpoint,float(times[-1])),first.y[:,-1],t_eval=times[cut+1:],**kwargs)
    if not last.success:return last
    from types import SimpleNamespace
    return SimpleNamespace(success=True,message="Native BDF restarted at the declared original horizon",
        t=np.concatenate((first.t,last.t)),y=np.concatenate((first.y,last.y),axis=1),
        nfev=first.nfev+last.nfev,njev=first.njev+last.njev,nlu=first.nlu+last.nlu)


class OriginalCapBDF:
    """A native BDF proxy retaining one observation/normalization closure.

    The first solver terminates exactly at the original cap. Its final dense
    step is consumed before the next step initializes a fresh native BDF with
    the same global RHS and settings. No saved fields or readouts are rebased.
    """
    def __init__(self,factory,checkpoint,fun,t0,y0,t_bound,**kwargs):
        if not t0<checkpoint<t_bound:raise ValueError("Native BDF restart must be internal to the integration interval")
        self._factory=factory;self._fun=fun;self._bound=t_bound;self._checkpoint=checkpoint;self._kwargs=kwargs
        self._solver=factory(fun,t0,y0,checkpoint,**kwargs);self._last=False
        self._totals=dict(nfev=0,njev=0,nlu=0)
    def __getattr__(self,name):
        if name=="status" and not self._last and self._solver.status=="finished":return "running"
        if name in self._totals:return self._totals[name]+getattr(self._solver,name)
        return getattr(self._solver,name)
    def step(self):
        if not self._last and self._solver.status=="finished":
            for name in self._totals:self._totals[name]+=getattr(self._solver,name)
            self._solver=self._factory(self._fun,self._checkpoint,self._solver.y.copy(),self._bound,**self._kwargs)
            self._last=True
        return self._solver.step()
    def dense_output(self):return self._solver.dense_output()


def worker(out,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Figure 4 native integration requires SLURM")
    import numpy as np
    start=time.monotonic();manifest=verify(out);launch_gate=check_launch_gate(out,manifest)
    job=manifest["jobs"][slot];directory=out/"workers"/f"{slot:04d}"
    directory.mkdir(exist_ok=False);status=dict(generated_at=now(),producer="running",status="running",complete=False,
        job=job,slurm_job_id=os.environ["SLURM_JOB_ID"],manifest_sha256=sha(out/"manifest.json"))
    write(directory/"status.json",status)
    native=None;forcing=None;law_audit=None
    def stop(*args):
        global STOP
        STOP=True
        if native is not None:setattr(native,"STOP",True)
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGUSR1,stop)
    try:
        native=native_module(out,job["platform"]);spec=job["protocol"];noise=dict(spec["noise"],sigma=job["sigma"])
        with np.load(out/job["input"]) as data:initial=data["state"].copy()
        deadline=start+(28000. if job["stage"]=="E4" else 6900.)
        if job["platform"]=="chemical":
            grid=native.ring.RingGrid(job["size"],job["cells"]);model=native.ring.native_model(job["model"])
            forcing=native.NoiseInput(grid,spec["horizon"],job["seed"],job["mean"],kind=job["kind"],**noise)
            check_retained_prefix(forcing,job,noise["knot"])
            check_paired_sigma1_field(forcing,job);law_audit=input_statistics(forcing,job,spec)
            write(directory/"noise_audit.json",law_audit);forcing=QueriedInput(forcing)
            checkpoint=job.get("native_solver_restart_at_original_horizon")
            if checkpoint is not None and ("retained_prefix" not in job or checkpoint!=job["retained_prefix"]["horizon"]):
                raise ValueError("Native restart must match the exact declared retained original horizon")
            original_bdf=native.BDF
            if checkpoint is not None:
                native.BDF=lambda fun,t0,y0,t_bound,**kwargs:OriginalCapBDF(original_bdf,checkpoint,fun,t0,y0,t_bound,**kwargs)
            try:
                status.update(native.integrate(model,grid,initial,forcing,dict(horizon=spec["horizon"],output_interval=spec["scalar_interval"],
                    primary_interval=spec["profile_interval"],full_interval=spec["full_interval"],solver=spec["solver"],diagnostics=spec["diagnostics"]),deadline,directory,status))
            finally:native.BDF=original_bdf
            if not status["complete"]:raise RuntimeError(status.get("message","Native chemical solve incomplete"))
        else:
            n=job["cells"];closure=native.Closure(job["size"],n)
            forcing=native.MaterialNoise(n,initial[-1],spec["horizon"],job["seed"],job["mean"],kind=job["kind"],
                knot=noise["knot"],coupling=noise["coupling"],tau=noise["tau"],ell=noise["ell"],sigma=noise["sigma"],bridge_base=noise.get("bridge_base"))
            check_retained_prefix(forcing,job,noise["knot"])
            check_paired_sigma1_field(forcing,job);law_audit=input_statistics(forcing,job,spec)
            write(directory/"noise_audit.json",law_audit);forcing=QueriedInput(forcing)
            t=np.arange(int(round(spec["horizon"]/spec["scalar_interval"]))+1)*spec["scalar_interval"]
            def rhs(t,y):
                if STOP or time.monotonic()>=deadline:raise TimeoutError("SLURM signal or native worker soft deadline")
                return closure.rhs(y,forcing(t))
            checkpoint=job.get("native_solver_restart_at_original_horizon")
            if checkpoint is not None and ("retained_prefix" not in job or checkpoint!=job["retained_prefix"]["horizon"]):
                raise ValueError("Native restart must match the exact declared retained original horizon")
            solution=spring_native_solve(rhs,lambda t,y:closure.jacobian(y,forcing(t)),initial,t,spec["solver"],checkpoint)
            if not solution.success or len(solution.t)!=len(t):raise RuntimeError(solution.message)
            y=solution.y.T;u,v,C=y[:,:n],y[:,n:-1],y[:,-1]
            if not np.isfinite(y).all() or y[:,:-1].min() < -5e-10 or C.min()<=0:raise ValueError("Native Spring field validity failed")
            amount=C*np.mean(u+v,axis=1);error=float(np.max(abs(amount/amount[0]-1)))
            if error>1e-7:raise ValueError("Native Spring physical amount drift")
            q=closure.material_positions;M0=C*np.mean(u,axis=1);M1=C*(u@np.exp(2j*np.pi*q))/n;M2=C*(u@np.exp(4j*np.pi*q))/n;P=M1/M1[0]
            np.savez_compressed(directory/"trace.npz",time=t,circumference=C,P_real=P.real,P_imag=P.imag,
                M1_real=M1.real,M1_imag=M1.imag,M0=M0,q1=abs(M1)/M0,q2=abs(M2)/M0,amount=amount,
                local_input_mean=np.array([forcing(tt).mean() for tt in t]),material_positions=q,**native.diagnostics(u,q))
            stride=int(round(spec["profile_interval"]/spec["scalar_interval"]))
            if not np.isclose(stride*spec["scalar_interval"],spec["profile_interval"]):raise ValueError("Profile clock must divide scalar clock")
            pi=np.arange(0,len(t),stride);fi=np.flatnonzero(np.isclose(t/spec["full_interval"],np.rint(t/spec["full_interval"])))
            np.savez_compressed(directory/"primary.npz",time=t[pi],primary=u[pi],inactive=v[pi],circumference=C[pi],material_positions=q)
            np.savez_compressed(directory/"full_states.npz",time=t[fi],states=y[fi],material_positions=q)
            np.savez_compressed(directory/"final_state.npz",state=y[-1],time=spec["horizon"])
            status.update(complete=True,status="complete",native_time=spec["horizon"],all_finite=True,
                minimum_concentration=float(y[:,:-1].min()),max_amount_relative_error=error,circumference_range=[float(C.min()),float(C.max())],nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu)
        status.update(native_integration_complete=True,validation_phase="observation_clocks_and_native_evolution")
        for name,dt in (("trace.npz",spec["scalar_interval"]),("primary.npz",spec["profile_interval"])):
            with np.load(directory/name) as data:
                expected=np.arange(int(round(spec["horizon"]/dt))+1)*dt
                np.testing.assert_allclose(data["time"],expected,rtol=0,atol=1e-8)
        law_audit["actual_native_input_queries"]=forcing.report(REFERENCES[job["model"]])
        write(directory/"noise_audit.json",law_audit)
        compare_sigma1_evolution(directory,job)
        check_retained_evolution_prefix(directory,job)
        if "retained_prefix" in job:status["retained_native_prefix_equivalence_passed"]=True
        artifacts={path.name:sha(path) for path in directory.iterdir() if path.is_file() and path.name not in ("result.json","status.json")}
        write(directory/"result.json",dict(generated_utc=now(),job=job,manifest_sha256=sha(out/"manifest.json"),
            artifact_sha256=artifacts,environment=environment(),launch_gate=launch_gate,wall_seconds=time.monotonic()-start,peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    except Exception as error:
        if law_audit is not None and isinstance(forcing,QueriedInput):
            law_audit["actual_native_input_queries"]=forcing.report(REFERENCES[job["model"]]);write(directory/"noise_audit.json",law_audit)
        status.update(complete=False,status="failed",failure_scope="native_validation" if status.get("native_integration_complete") else "native_integration_or_input",
            error=f"{type(error).__name__}: {error}",traceback=traceback.format_exc())
        raise
    finally:
        status.update(generated_at=now(),producer="terminal",elapsed_seconds=time.monotonic()-start,
            output_sha256={p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!="status.json"})
        write(directory/"status.json",status)
    return dict(slot=slot,status="complete",wall_seconds=time.monotonic()-start)


def check_launch_gate(out,manifest):
    """Pilot preparation may precede E0 delivery; its native execution may not."""
    if not manifest.get("native_production_held_for_E0_interim"):return dict(e0_delivery_verified_in_frozen_inputs=True)
    path=out/"results/e0_interim_launch_release.json"
    if not path.is_file():raise RuntimeError("Native pilot remains held until concrete E0 interim delivery")
    release=read(path)
    if release.get("e0_interim_delivered") is not True or release.get("manifest_sha256")!=sha(out/"manifest.json") or not release.get("interim_artifacts"):
        raise RuntimeError("E0 delivery release does not identify this payload and its concrete interim artifacts")
    for artifact,digest in release["interim_artifacts"].items():
        if sha(artifact)!=digest:raise ValueError("Delivered E0 artifact changed")
    return dict(e0_delivery_release_path=str(path),e0_delivery_release_sha256=sha(path))


def check_retained_prefix(forcing,job,knot):
    """An extended existing seed must preserve its complete Gaussian prefix."""
    if "retained_prefix" not in job:return
    import numpy as np
    old=job["retained_prefix"]
    if sha(old["noise_audit"])!=old["noise_audit_sha256"]:raise ValueError("Retained forcing audit changed")
    if sha(old["original_status"])!=old["original_status_sha256"]:raise ValueError("Retained native-prefix receipt changed")
    count=int(round(old["horizon"]/knot))+1
    measured=hashlib.sha256(np.asarray(forcing.knots[:count],dtype="<f8",order="C").tobytes()).hexdigest()
    if measured!=old["gaussian_knots_sha256"]:raise ValueError("E4 changed a retained seed's Gaussian input prefix")


def check_retained_reference_artifacts(job):
    """Verify all old native-prefix reference files before a costly long solve."""
    if "retained_prefix" not in job:return {}
    prefix=job["retained_prefix"];old=prefix["original_case"]
    if sha(prefix["original_status"])!=prefix["original_status_sha256"]:raise ValueError("Retained native reference status changed")
    status=read(prefix["original_status"]);paths={Path(old[key]):old[key+"_sha256"] for key in ("trace","profile")}
    paths[Path(old["trace"]).parent/"final_state.npz"]=status["output_sha256"]["final_state.npz"]
    for path,digest in paths.items():
        if sha(path)!=digest:raise ValueError("Retained native reference artifact changed: "+str(path))
    return {str(path):digest for path,digest in paths.items()}


def collect(out):
    manifest=verify(out);missing=[];failed=[];receipts=[];cases=[]
    for job in manifest["jobs"]:
        folder=out/"workers"/f"{job['slot']:04d}";result_path=folder/"result.json";status_path=folder/"status.json"
        if not result_path.exists():
            if status_path.exists() and read(status_path).get("status")=="failed":failed.append(job["slot"])
            else:missing.append(job["slot"])
            continue
        result=read(result_path);status=read(status_path)
        if result["job"]!=job or result["manifest_sha256"]!=sha(out/"manifest.json") or status.get("complete") is not True:
            raise ValueError("Native worker result does not match its frozen experiment")
        actual={p.name for p in folder.iterdir() if p.is_file() and p.name not in ("status.json","result.json")}
        if actual!=set(result["artifact_sha256"]):raise ValueError("Scientific artifact inventory is incomplete")
        for name,digest in result["artifact_sha256"].items():
            if sha(folder/name)!=digest:raise ValueError("Native artifact changed: "+str(folder/name))
        if status["output_sha256"]["result.json"]!=sha(result_path):raise ValueError("Native result changed after terminal receipt")
        receipts.append(dict(slot=job["slot"],result_sha256=sha(result_path),scientific_artifacts=len(actual),
            native_complete=True,all_artifact_hashes_verified=True,wall_seconds=result["wall_seconds"],
            peak_resident_memory_KiB=result["peak_resident_memory_KiB"]))
        cases.append({key:job[key] for key in ("model","size","fraction","kind","seed","sigma","independent_noise_seed","new_gaussian_path")}|dict(
            case_id=f"{manifest['stage']}_{job['slot']:04d}",horizon=job["protocol"]["horizon"],
            scalar_interval=job["protocol"]["scalar_interval"],profile_interval=job["protocol"]["profile_interval"],
            trace=str(folder/"trace.npz"),trace_sha256=result["artifact_sha256"]["trace.npz"],
            profile=str(folder/"primary.npz"),profile_sha256=result["artifact_sha256"]["primary.npz"],
            basal_input_sha256=job["input_sha256"],native_manifest_sha256=sha(out/"manifest.json"),
            native_result_sha256=sha(result_path),observation_companion_of=job.get("observation_companion_of"),
            retained_prefix=job.get("retained_prefix"),paired_sigma1_source=job.get("paired_sigma1_source"),
            retained_native_prefix_equivalence=dict(path=str(folder/"retained_native_prefix_equivalence.json"),
                sha256=result["artifact_sha256"]["retained_native_prefix_equivalence.json"])
                if "retained_native_prefix_equivalence.json" in result["artifact_sha256"] else None))
    report=dict(generated_utc=now(),stage=manifest["stage"],manifest_sha256=sha(out/"manifest.json"),
        completed=len(receipts),total=len(manifest["jobs"]),missing_slots=missing,failed_slots=failed,
        complete=len(receipts)==len(manifest["jobs"]),all_observed_artifact_hashes_verified=bool(receipts),
        counts=manifest["counts"],scientific_analysis_complete=False,
        scope="Native execution and exact output identities only; missing work is not censoring and no scientific hypothesis is decided here.")
    write(out/"results/artifact_verification.json",dict(generated_utc=now(),workers=receipts))
    write(out/"results/case_census.json",dict(generated_utc=now(),source_manifest_sha256=manifest["parent_manifest_sha256"],
        native_run_manifest_sha256=sha(out/"manifest.json"),cases=cases,scientific_analysis_complete=False))
    write(out/"results/status.json",report);return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest="action",required=True)
    prep=sub.add_parser("prepare");prep.add_argument("--out",type=Path,required=True);prep.add_argument("--stage",choices=("E1_pilot","E1","E2","E3","E4"),required=True)
    prep.add_argument("--gate",type=Path,required=True);prep.add_argument("--strengths",type=Path);prep.add_argument("--triggers",type=Path);prep.add_argument("--basals",type=Path)
    prep.add_argument("--companion-census",type=Path,help="Additional audited sigma=1/control records, for strengths beyond the fully replicated original grid.")
    prep.add_argument("--condition-subset",type=Path,help="E1 only: complete model/size/f conditions for a disjoint companion-ready batch.")
    check=sub.add_parser("verify");check.add_argument("--out",type=Path,required=True)
    preflight=sub.add_parser("preflight");preflight.add_argument("--out",type=Path,required=True)
    forcing=sub.add_parser("forcing-preflight");forcing.add_argument("--out",type=Path,required=True);forcing.add_argument("--shard",type=int);forcing.add_argument("--shards",type=int)
    forcing_collect=sub.add_parser("collect-forcing-preflight");forcing_collect.add_argument("--out",type=Path,required=True);forcing_collect.add_argument("--shards",type=int,required=True)
    collection=sub.add_parser("collect");collection.add_argument("--out",type=Path,required=True)
    work=sub.add_parser("worker");work.add_argument("--out",type=Path,required=True);work.add_argument("--slot",type=int,required=True)
    args=parser.parse_args();out=args.out.resolve()
    if args.action=="prepare":result=prepare(out,args.stage,args.gate,strengths=read(args.strengths) if args.strengths else None,triggers=read(args.triggers) if args.triggers else None,extra_basals=read(args.basals) if args.basals else None,companion_census=args.companion_census,condition_subset=read(args.condition_subset) if args.condition_subset else None)
    elif args.action=="verify":result=dict(verified=True,jobs=len(verify(out)["jobs"]),manifest_sha256=sha(out/"manifest.json"))
    elif args.action=="preflight":result=native_preflight(out)
    elif args.action=="forcing-preflight":result=forcing_preflight(out,args.shard,args.shards)
    elif args.action=="collect-forcing-preflight":result=collect_forcing_preflight(out,args.shards)
    elif args.action=="collect":result=collect(out)
    else:result=worker(out,args.slot)
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
