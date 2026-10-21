"""Independent native-field/forcing audits for the prospective Figure 4 runs.

This reconstructs saved global observations and physical invariants. It does
not estimate persistence, assign regimes or decide any Figure 4 hypothesis.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import time
import traceback

import numpy as np

import signature_figure4_campaign_20261020 as C


def array_sha(values):
    return hashlib.sha256(np.asarray(values,dtype="<f8",order="C").tobytes()).hexdigest()


def file_snapshot(path):
    """Local unchanged-file guard after full hashing; NFS device IDs may vary."""
    info=Path(path).stat()
    return dict(inode=info.st_ino,size=info.st_size,mtime_ns=info.st_mtime_ns,ctime_ns=info.st_ctime_ns)


def check_snapshots(rows):
    for path,record in rows.items():
        if file_snapshot(path)!=record:raise ValueError("Hashed native scientific file changed: "+path)


def profile_observations(active,positions,circumference):
    """Physical circular amounts, normalized initial moment and raw regions."""
    active=np.asarray(active);positions=np.asarray(positions);circumference=np.asarray(circumference)
    if active.ndim!=2 or positions.shape!=(active.shape[1],):raise ValueError("Native profile layout changed")
    physical=np.broadcast_to(circumference,(len(active),))/active.shape[1]
    amount=physical*np.sum(active,axis=1)
    first=physical*(active@np.exp(2j*np.pi*positions))
    second=physical*(active@np.exp(4j*np.pi*positions))
    if abs(first[0])<=1e-12*amount[0]:raise ValueError("Basal moment is not a direction reference")
    polarity=first/first[0]
    q1=np.divide(abs(first),amount,out=np.zeros_like(amount),where=amount>0)
    q2=np.divide(abs(second),amount,out=np.zeros_like(amount),where=amount>0)
    contrast=np.ptp(active,axis=1)
    mask=active>=active.min(axis=1)[:,None]+.5*contrast[:,None]
    regions=np.sum(mask&~np.roll(mask,1,axis=1),axis=1)
    regions=np.where(contrast>1e-12,regions,0)
    return dict(M0=amount,M1_real=first.real,M1_imag=first.imag,P_real=polarity.real,P_imag=polarity.imag,
        q1=q1,q2=q2,contrast=contrast,current_half_contrast_regions=regions)


def _compare(observed,reconstructed,field):
    observed=np.asarray(observed);reconstructed=np.asarray(reconstructed)
    if observed.shape!=reconstructed.shape or not np.isfinite(observed).all() or not np.isfinite(reconstructed).all():
        raise ValueError("Nonfinite or incompatible native observation: "+field)
    if field=="current_half_contrast_regions":np.testing.assert_array_equal(observed,reconstructed)
    else:np.testing.assert_allclose(observed,reconstructed,rtol=2e-9,atol=2e-10,err_msg=field)
    return float(np.max(abs(observed-reconstructed))) if observed.size else 0.


def check_native_artifact_reuse(run,job,manifest):
    """Rehash the actual native execution behind explicitly adapted metadata."""
    reuse=job.get("native_artifact_reuse")
    if reuse is None:return None
    original=Path(reuse["run"]).resolve();source=C.verify(original)
    if C.sha(original/"manifest.json")!=reuse["native_manifest_sha256"]:
        raise ValueError("Reused native source manifest changed")
    source_job=source["jobs"][reuse["slot"]];folder=original/"workers"/f"{reuse['slot']:04d}"
    if C.sha(folder/"status.json")!=reuse["status_sha256"] or C.sha(folder/"result.json")!=reuse["result_sha256"]:
        raise ValueError("Reused actual native status/result changed")
    status=C.read(folder/"status.json");result=C.read(folder/"result.json")
    if not status.get("complete") or result["job"]!=source_job or result["manifest_sha256"]!=reuse["native_manifest_sha256"]:
        raise ValueError("Only complete actually executed native workers can be reused")
    if C.sha(folder/"result.json")!=status["output_sha256"]["result.json"]:
        raise ValueError("Reused native terminal result binding changed")
    proof=Path(reuse["worker_audit"])
    if C.sha(proof)!=reuse["worker_audit_sha256"]:raise ValueError("Reused independent native worker audit changed")
    source_audit=original/"results/independent_native_audit_20261020"
    if (Path(reuse["source_audit_manifest"]).resolve()!=source_audit/"manifest.json"
            or Path(reuse["source_audit_status"]).resolve()!=source_audit/"status.json"
            or C.sha(source_audit/"manifest.json")!=reuse["source_audit_manifest_sha256"]
            or C.sha(source_audit/"status.json")!=reuse["source_audit_status_sha256"]):
        raise ValueError("Reused full source-audit manifest/status changed")
    audit_manifest,_=verify_audit(source_audit);aggregate=C.read(source_audit/"status.json")
    links=[row for row in aggregate["receipts"] if row["slot"]==reuse["slot"]]
    if (not aggregate["full_native_run_covered"] or aggregate["native_manifest_sha256"]!=reuse["native_manifest_sha256"]
            or aggregate["audit_manifest_sha256"]!=reuse["source_audit_manifest_sha256"] or reuse["slot"] not in audit_manifest["tasks"]
            or len(links)!=1 or links[0].get("passed") is not True or links[0]["sha256"]!=reuse["worker_audit_sha256"]
            or Path(links[0]["path"]).resolve()!=proof.resolve()):
        raise ValueError("Reused actual worker is not bound to its full independent source audit")
    prior=C.read(proof)
    if prior.get("passed") is not True or prior["native_manifest_sha256"]!=reuse["native_manifest_sha256"] or prior["native_result_sha256"]!=reuse["result_sha256"] or prior["native_status_sha256"]!=reuse["status_sha256"]:
        raise ValueError("Reused actual native worker lacks its exact independent audit")
    check_snapshots(prior["hashed_native_file_snapshots"])
    ignored={"slot","full_stage_slot","native_artifact_reuse","native_recovery_source","native_solver_restart_at_original_horizon"}
    if {k:v for k,v in job.items() if k not in ignored}!={k:v for k,v in source_job.items() if k not in ignored}:
        raise ValueError("Reused native worker changes a scientific condition or Gaussian seed")
    if source_job.get("full_stage_slot",source_job["slot"])!=job.get("full_stage_slot",job["slot"]):
        raise ValueError("Reused native worker changes the original parent slot")
    if job.get("native_solver_restart_at_original_horizon")!=source_job.get("native_solver_restart_at_original_horizon"):
        raise ValueError("Reused native solver restart declaration changed")
    for relative,digest in source["source_sha256"].items():
        if relative.startswith(("frozen_sources/chemical/","frozen_sources/spring/")) or relative in (
                "frozen_sources/campaign.py","frozen_sources/replication_core.py","frozen_sources/replication_driver.py"):
            if manifest["source_sha256"].get(relative)!=digest:
                raise ValueError("Reused native model/solver/forcing source differs")
    expected={p.name for p in folder.iterdir() if p.is_file() and p.name not in ("result.json","status.json")}
    if set(result["artifact_sha256"])!=expected or result["artifact_sha256"]!=reuse["artifact_sha256"]:
        raise ValueError("Reused native artifact inventory changed")
    selected=Path(run)/"workers"/f"{job['slot']:04d}"
    for name,digest in result["artifact_sha256"].items():
        if (selected/name).resolve()!=(folder/name).resolve() or C.sha(folder/name)!=digest:
            raise ValueError("Selected reused artifact is not the exact original native file")
    provenance_paths=[original/"manifest.json",folder/"result.json",folder/"status.json",proof,
        source_audit/"manifest.json",source_audit/"freeze.json",source_audit/"status.json"]
    provenance_paths.extend(Path(path) for path in audit_manifest["source_sha256"])
    return dict(actual_native_run=str(original),actual_native_slot=reuse["slot"],
        native_manifest_sha256=reuse["native_manifest_sha256"],native_result_sha256=reuse["result_sha256"],
        native_status_sha256=reuse["status_sha256"],exact_native_artifacts=len(expected),
        source_provenance_snapshots={str(path):file_snapshot(path) for path in provenance_paths},
        metadata_only_adapter=True,native_reintegration_performed=False)


def audit_worker(run,slot,*,manifest=None):
    started=time.monotonic();run=Path(run).resolve();manifest=C.verify(run) if manifest is None else manifest
    job=manifest["jobs"][slot];folder=run/"workers"/f"{slot:04d}"
    status=C.read(folder/"status.json");result=C.read(folder/"result.json")
    if not status.get("complete") or status["job"]!=job or status["manifest_sha256"]!=C.sha(run/"manifest.json") or result["job"]!=job or result["manifest_sha256"]!=C.sha(run/"manifest.json"):
        raise ValueError("Native worker is not complete on its frozen experiment")
    reused=check_native_artifact_reuse(run,job,manifest)
    artifacts={p.name for p in folder.iterdir() if p.is_file() and p.name not in ("status.json","result.json")}
    if artifacts!=set(result["artifact_sha256"]):raise ValueError("Native artifact inventory changed")
    snapshots={} if reused is None else dict(reused["source_provenance_snapshots"])
    for name,digest in result["artifact_sha256"].items():
        snapshots[str(folder/name)]=file_snapshot(folder/name)
        if C.sha(folder/name)!=digest:raise ValueError("Native scientific artifact changed: "+name)
        check_snapshots({str(folder/name):snapshots[str(folder/name)]})
    if C.sha(folder/"result.json")!=status["output_sha256"]["result.json"]:raise ValueError("Native terminal receipt changed")
    for name in ("result.json","status.json"):snapshots[str(folder/name)]=file_snapshot(folder/name)
    with np.load(run/job["input"]) as data:
        initial=data["state"].copy()
        initial_positions=data["material_positions"].copy() if job["platform"]=="spring" else data["x"].copy()
    if array_sha(initial)!=job["state_sha256"]:raise ValueError("Exact native initial state changed")
    n=job["cells"];spring=job["platform"]=="spring";spec=job["protocol"]
    with np.load(folder/"primary.npz") as data:
        t=data["time"].copy();active=data["primary"].copy()
        positions=data["material_positions"].copy() if spring else data["x"].copy()/job["size"]
        np.testing.assert_array_equal(data["material_positions"] if spring else data["x"],initial_positions)
        circumference=data["circumference"].copy()
        if spring:
            inactive=data["inactive"]
            if inactive.shape!=active.shape or not np.isfinite(inactive).all() or inactive.min() < -5e-10:raise ValueError("Invalid inactive Spring profiles")
    np.testing.assert_array_equal(t,np.arange(int(round(spec["horizon"]/spec["profile_interval"]))+1)*spec["profile_interval"])
    expected_active=initial[:n] if spring else initial[0]
    np.testing.assert_array_equal(active[0],expected_active)
    np.testing.assert_allclose(positions,(np.arange(n)+.5)/n,rtol=0,atol=2e-15)
    if active.shape!=(len(t),n) or not np.isfinite(active).all() or active.min() < -5e-10 or np.any(circumference<=0):
        raise ValueError("Invalid native active profiles or physical geometry")
    if not spring and float(circumference)!=job["size"]:raise ValueError("Chemical physical circumference changed")
    if spring and float(circumference[0])!=initial[-1]:raise ValueError("Spring actual initial circumference changed")
    reconstructed=profile_observations(active,positions,circumference);differences={}
    with np.load(folder/"trace.npz") as trace:
        np.testing.assert_array_equal(trace["time"],np.arange(int(round(spec["horizon"]/spec["scalar_interval"]))+1)*spec["scalar_interval"])
        indices=np.searchsorted(trace["time"],t);np.testing.assert_array_equal(trace["time"][indices],t)
        for field,prediction in reconstructed.items():differences[field]=_compare(trace[field][indices],prediction,field)
        if spring:np.testing.assert_array_equal(trace["circumference"][indices],circumference)
    del active
    with np.load(folder/"full_states.npz") as data:
        full_time=data["time"].copy();fields=data["states"].copy()
    np.testing.assert_array_equal(full_time,np.arange(int(round(spec["horizon"]/spec["full_interval"]))+1)*spec["full_interval"])
    np.testing.assert_array_equal(fields[0],initial)
    expected_shape=(len(full_time),2*n+1) if spring else (len(full_time),*initial.shape)
    if fields.shape!=expected_shape or not np.isfinite(fields).all():raise ValueError("Full native coupled state changed")
    if spring:
        if fields[:,:-1].min() < -5e-10 or np.any(fields[:,-1]<=0):raise ValueError("Spring native field positivity failed")
        pools=(fields[:,-1]*np.mean(fields[:,:n]+fields[:,n:-1],axis=1))[:,None]
    else:
        if fields.min() < -5e-10:raise ValueError("Native field positivity failed")
        pools=np.stack([job["size"]*np.mean(fields[:,k]+fields[:,k+1],axis=1) for k in range(0,initial.shape[0],2)],axis=1)
    pool_drift=float(np.max(abs(pools/pools[0]-1)))
    if pool_drift>1e-7:raise ValueError("Saved native physical pool conservation failed")
    with np.load(folder/"final_state.npz") as final:
        if float(final["time"])!=spec["horizon"]:raise ValueError("Final native horizon changed")
        np.testing.assert_allclose(final["state"],fields[-1],rtol=2e-9,atol=2e-10)
    del fields
    native=C.native_module(run,job["platform"]);noise=dict(spec["noise"],sigma=job["sigma"])
    if spring:
        forcing=native.MaterialNoise(n,job["basal_circumference"],spec["horizon"],job["seed"],job["mean"],kind=job["kind"],
            knot=noise["knot"],coupling=noise["coupling"],tau=noise["tau"],ell=noise["ell"],sigma=noise["sigma"],bridge_base=noise.get("bridge_base"))
    else:forcing=native.NoiseInput(native.ring.RingGrid(job["size"],n),spec["horizon"],job["seed"],job["mean"],kind=job["kind"],**noise)
    law=C.read(folder/"noise_audit.json")
    if forcing.knots is not None:
        measured=array_sha(forcing.knots)
        if measured!=law["gaussian_knots_sha256"]:raise ValueError("Native forcing cannot be reproduced from the declared seed/law")
        C.check_paired_sigma1_field(forcing,job);C.check_retained_prefix(forcing,job,noise["knot"])
    else:measured=None
    if law["theoretical_input_mean"]!=job["mean"] or law["sigma"]!=job["sigma"] or not law["finite_run_mean_is_not_normalized"]:
        raise ValueError("Native fixed-mean forcing declaration changed")
    queries=law["actual_native_input_queries"]
    if queries["clipping_or_pathwise_normalization"] or queries["minimum"]<0 or not np.isfinite(queries["maximum"]):raise ValueError("Native input-query guard failed")
    if "retained_prefix" in job:
        prefix=C.read(folder/"retained_native_prefix_equivalence.json")
        if not prefix["passed"] or prefix["rule"]!=C.ODE_PREFIX_RULE:raise ValueError("E4 retained native evolution prefix is unqualified")
    check_snapshots(snapshots)
    return dict(generated_utc=C.now(),passed=True,slot=slot,model=job["model"],size=job["size"],kind=job["kind"],sigma=job["sigma"],seed=job["seed"],
        native_manifest_sha256=C.sha(run/"manifest.json"),native_result_sha256=C.sha(folder/"result.json"),native_status_sha256=C.sha(folder/"status.json"),
        native_artifacts_hashed=len(artifacts),hashed_native_file_snapshots=snapshots,exact_initial_native_fields=True,full_coupled_fields_positive_and_finite=True,
        full_physical_pool_relative_drift=pool_drift,global_observations_recomputed_from_profiles=differences,
        profile_samples_audited=len(t),scalar_samples_audited_at_saved_profiles=len(t),
        saved_profile_cadence=spec["profile_interval"],native_scalar_cadence=spec["scalar_interval"],
        unobserved_scalar_field_scope="Only shared native-profile/scalar times can be reconstructed; full scalar cadence is covered when profile_interval==scalar_interval.",
        exact_forcing_gaussian_sha256=measured,forcing_independently_regenerated=True,
        actual_native_artifact_reuse=reused,
        retained_native_prefix_qualified="retained_prefix" in job,wall_seconds=time.monotonic()-started,
        scientific_hypotheses_evaluated=False)


def prepare(run,out,shards=24,tasks=None):
    run=Path(run).resolve();out=Path(out).resolve();manifest=C.verify(run)
    if out.exists() or not out.is_relative_to(run/"results"):raise ValueError("Choose a fresh native audit directory within this run's results")
    selected=list(range(len(manifest["jobs"]))) if tasks is None else list(tasks)
    if not selected or len(set(selected))!=len(selected) or not set(selected)<=set(range(len(manifest["jobs"]))):raise ValueError("Invalid native audit task subset")
    for folder in ("frozen_sources","workers","logs"):(out/folder).mkdir(parents=True)
    sources={}
    for path in (Path(__file__).resolve(),Path(C.__file__).resolve()):
        target=out/"frozen_sources"/path.name;shutil.copy2(path,target);sources[str(target)]=C.sha(target)
    record=dict(generated_utc=C.now(),run=str(run),native_manifest_sha256=C.sha(run/"manifest.json"),tasks=selected,
        full_run_workers=len(manifest["jobs"]),shards=min(shards,len(selected)),source_sha256=sources,
        scope="Independent source/artifact/field/forcing integrity only; no persistence or hypothesis inference.")
    C.write(out/"manifest.json",record);C.write(out/"freeze.json",dict(manifest_sha256=C.sha(out/"manifest.json")))
    return dict(out=str(out),tasks=len(selected),shards=record["shards"],manifest_sha256=C.sha(out/"manifest.json"))


def verify_audit(out):
    if C.sha(out/"manifest.json")!=C.read(out/"freeze.json")["manifest_sha256"]:raise ValueError("Frozen native audit manifest changed")
    manifest=C.read(out/"manifest.json")
    for path,digest in manifest["source_sha256"].items():
        if C.sha(path)!=digest:raise ValueError("Frozen native auditor changed")
    run=Path(manifest["run"])
    if C.sha(run/"manifest.json")!=manifest["native_manifest_sha256"]:raise ValueError("Frozen native run changed")
    return manifest,C.verify(run)


def worker(out,shard):
    manifest,native=verify_audit(out);run=Path(manifest["run"])
    selected=manifest["tasks"][shard::manifest["shards"]];rows=[]
    for slot in selected:
        try:row=audit_worker(run,slot,manifest=native)
        except Exception as error:row=dict(generated_utc=C.now(),slot=slot,passed=False,error=f"{type(error).__name__}: {error}",traceback=traceback.format_exc())
        C.write(out/"workers"/f"{slot:04d}.json",row);rows.append(row)
    report=dict(generated_utc=C.now(),shard=shard,tasks=selected,passed=all(r["passed"] for r in rows),failed=[r["slot"] for r in rows if not r["passed"]])
    C.write(out/f"shard_{shard:03d}.json",report);return report


def collect(out):
    manifest,native=verify_audit(out);rows=[];missing=[];failed=[]
    for slot in manifest["tasks"]:
        path=out/"workers"/f"{slot:04d}.json"
        if not path.exists():missing.append(slot);continue
        row=C.read(path)
        if not row["passed"]:failed.append(slot)
        elif row["native_manifest_sha256"]!=manifest["native_manifest_sha256"]:raise ValueError("Native audit receipt refers to another run")
        elif "hashed_native_file_snapshots" in row:check_snapshots(row["hashed_native_file_snapshots"])
        rows.append(dict(slot=slot,path=str(path),sha256=C.sha(path),passed=row["passed"]))
    report=dict(generated_utc=C.now(),passed=not missing and not failed,selected_workers=len(manifest["tasks"]),audited_workers=len(rows),
        full_native_run_workers=len(native["jobs"]),full_native_run_covered=len(manifest["tasks"])==len(native["jobs"]),
        missing_tasks=missing,failed_tasks=failed,native_manifest_sha256=manifest["native_manifest_sha256"],
        audit_manifest_sha256=C.sha(out/"manifest.json"),receipts=rows,scientific_hypotheses_evaluated=False)
    C.write(out/"status.json",report);return report


def merge(run,parts,out):
    """Bind disjoint scoped audits to every selected native worker exactly once."""
    run=Path(run).resolve();out=Path(out).resolve();native=C.verify(run);receipts=[];coverage=[];all_rows=[]
    if out.exists() or not out.is_relative_to(run/"results"):raise ValueError("Choose a fresh merged receipt within this run's results")
    for part in parts:
        part=Path(part).resolve();manifest,_=verify_audit(part)
        if Path(manifest["run"])!=run:raise ValueError("Mixed native run in audit aggregation")
        status=collect(part)
        if not status["passed"]:raise ValueError("Scoped native audit incomplete or failed")
        overlap=set(coverage)&set(manifest["tasks"])
        if overlap:raise ValueError("Native audit scopes overlap: "+str(sorted(overlap)))
        coverage.extend(manifest["tasks"]);all_rows.extend(status["receipts"])
        receipts.append(dict(path=str(part/"status.json"),sha256=C.sha(part/"status.json"),manifest_sha256=C.sha(part/"manifest.json"),tasks=manifest["tasks"]))
    if set(coverage)!=set(range(len(native["jobs"]))):raise ValueError("Native audit aggregation does not cover the complete frozen run")
    report=dict(generated_utc=C.now(),passed=True,native_manifest_sha256=C.sha(run/"manifest.json"),
        selected_workers=len(coverage),full_native_run_covered=True,missing_tasks=[],failed_tasks=[],
        scoped_audits=receipts,receipts=sorted(all_rows,key=lambda row:row["slot"]),scientific_hypotheses_evaluated=False)
    C.write(out,report);return report


def e4_population(job):
    if not isinstance(job.get("new_gaussian_path"),bool):raise ValueError("E4 requires an explicit Gaussian population identity")
    original="retained_prefix" in job
    if original==job["new_gaussian_path"]:raise ValueError("E4 original/top-up identity conflicts with retained native prefix")
    return "original_seed_continuation" if original else "additional_independent_seed"


def selected_native_source_gates(run,native):
    """Independently enforce exact preferred-original/failed-prefix selection."""
    run=Path(run).resolve();authority=run/"results/complete_original_slot_selection_authority.json";record=C.read(authority)
    original=Path(native["original_attempt"]["run"]).resolve();source=C.verify(original)
    original_digest=C.sha(original/"manifest.json")
    recipe_path=Path(native["recovery_recipe"]["path"]);recipe=C.read(recipe_path)
    if (C.sha(recipe_path)!=native["recovery_recipe"]["sha256"] or record["recipe_sha256"]!=C.sha(recipe_path)
            or original_digest!=native["original_attempt"]["manifest_sha256"] or record["original_manifest_sha256"]!=original_digest
            or recipe["original_manifest_sha256"]!=original_digest or recipe["prefix_rule"]!=C.ODE_PREFIX_RULE
            or recipe["native_adapter_sha256"]!=C.sha(run/"frozen_sources/figure4_campaign_20261020.py")):
        raise ValueError("E4 selected native-source recipe/original manifest binding changed")
    if (record["native_manifest_sha256"]!=C.sha(run/"manifest.json") or record["all_original_slots"]!=len(native["jobs"])
            or len(source["jobs"])!=len(native["jobs"]) or record["response_outcomes_used_for_selection"]
            or record["selection"]!=native["actual_native_selection"]):
        raise ValueError("E4 selected native-source authority differs from the frozen manifest")
    ignored={"slot","full_stage_slot","native_artifact_reuse","native_recovery_source","native_solver_restart_at_original_horizon"}
    batches={}
    for selection in record["selection"]:
        slot=selection["original_slot"]
        if slot not in range(len(native["jobs"])) or slot in batches:raise ValueError("E4 original native-source slots overlap or are unknown")
        job=native["jobs"][slot];old_job=source["jobs"][slot];reuse=job.get("native_artifact_reuse")
        if selection["actual"]!=reuse or {k:v for k,v in job.items() if k not in ignored}!={k:v for k,v in old_job.items() if k not in ignored}:
            raise ValueError("E4 selected native condition/source differs from the original slot")
        original_status=original/"workers"/f"{slot:04d}"/"status.json"
        if Path(selection["original_status"]["path"]).resolve()!=original_status or C.sha(original_status)!=selection["original_status"]["sha256"]:
            raise ValueError("Preserved original native attempt path/hash changed")
        status=C.read(original_status);role=selection["role"]
        if status.get("complete"):
            if role!="passing_original" or Path(reuse["run"]).resolve()!=original or reuse["slot"]!=slot:
                raise ValueError("A passing original native outcome must be preferred exactly")
        elif (role!="qualified_old_cap_restart" or status.get("status")!="failed" or not status.get("native_integration_complete")
                or "check_retained_evolution_prefix" not in status.get("traceback","") or "retained_prefix" not in old_job
                or job.get("native_solver_restart_at_original_horizon")!=old_job["retained_prefix"]["horizon"]
                or Path(reuse["run"]).resolve()==original):
            raise ValueError("A restart may replace only an actual failed native-prefix outcome")
        actual_status=Path(reuse["run"])/"workers"/f"{reuse['slot']:04d}"/"status.json"
        batches[slot]=dict(native_source_selection_authority=str(authority),native_source_selection_authority_sha256=C.sha(authority),
            role=role,actual_native_run=reuse["run"],actual_native_slot=reuse["slot"],
            actual_native_slurm_job_id=C.read(actual_status)["slurm_job_id"],
            actual_native_worker_audit=reuse["worker_audit"],actual_native_worker_audit_sha256=reuse["worker_audit_sha256"])
    if set(batches)!=set(range(len(native["jobs"]))):raise ValueError("E4 original native-source slots are incomplete")
    return batches


def e4_gates(run,audit,out):
    """Link every long noisy path to numerical, native-prefix and control gates."""
    run=Path(run).resolve();audit=Path(audit).resolve();out=Path(out).resolve();native=C.verify(run)
    if native["stage"]!="E4" or out.exists() or not out.is_relative_to(run/"results"):raise ValueError("Choose a fresh E4 numerical gate receipt")
    master=C.read(audit)
    if not master["passed"] or not master["full_native_run_covered"] or master["native_manifest_sha256"]!=C.sha(run/"manifest.json"):
        raise ValueError("E4 requires a complete independent native audit")
    detailed={}
    for record in master["receipts"]:
        path=Path(record["path"])
        if C.sha(path)!=record["sha256"]:raise ValueError("Independent E4 worker audit changed")
        row=C.read(path)
        if not row["passed"]:raise ValueError("Failed independent E4 worker audit")
        check_snapshots(row["hashed_native_file_snapshots"])
        detailed[record["slot"]]=dict(path=str(path),sha256=record["sha256"])
    expected=set(range(len(native["jobs"])))
    if set(detailed)!=expected or len(master["receipts"])!=len(expected):raise ValueError("E4 native worker audit coverage is incomplete")
    batches={}
    if native.get("native_execution_performed_in_immutable_source_runs"):
        batches=selected_native_source_gates(run,native)
    for path in (() if native.get("native_execution_performed_in_immutable_source_runs") else sorted((run/"results/submission_batches").glob("*.json"))):
        record=C.read(path)
        if record["manifest_sha256"]!=C.sha(run/"manifest.json"):raise ValueError("E4 native submission refers to another manifest")
        gate=Path(record["engineering_gate"])
        if C.sha(gate)!=record["engineering_gate_sha256"] or C.read(gate).get("passed") is not True:raise ValueError("E4 native engineering gate is unqualified")
        for slot in record["task_indices"]:
            if slot in batches:raise ValueError("E4 native submission task scopes overlap")
            batches[slot]=dict(submission=str(path),submission_sha256=C.sha(path),array_job_id=record["array_job_id"],
                engineering_gate=str(gate),engineering_gate_sha256=record["engineering_gate_sha256"])
    if set(batches)!=expected:raise ValueError("E4 engineering release coverage is incomplete")
    means={(j["model"],j["size"],j["fraction"]):j for j in native["jobs"] if j["kind"]=="mean"}
    zeros={(j["model"],j["size"]):j for j in native["jobs"] if j["kind"]=="basal"}
    rows=[]
    for job in native["jobs"]:
        if job["kind"]!="noise":continue
        population=e4_population(job);prefix=None;folder=run/"workers"/f"{job['slot']:04d}"
        if population=="original_seed_continuation":
            path=folder/"retained_native_prefix_equivalence.json";record=C.read(path)
            if not record["passed"] or record["rule"]!=C.ODE_PREFIX_RULE:raise ValueError("E4 original native evolution prefix is unqualified")
            prefix=dict(path=str(path),sha256=C.sha(path),passed=True)
        mean=means[(job["model"],job["size"],job["fraction"])];zero=zeros[(job["model"],job["size"])]
        if mean["protocol"]["horizon"]!=job["protocol"]["horizon"] or zero["protocol"]["horizon"]!=job["protocol"]["horizon"]:
            raise ValueError("E4 controls do not cover the complete long window")
        rows.append(dict(slot=job["slot"],model=job["model"],size=job["size"],fraction=job["fraction"],seed=job["seed"],
            population=population,new_gaussian_path=job["new_gaussian_path"],retained_prefix=job.get("retained_prefix"),
            native_audit=detailed[job["slot"]],native_engineering_release=batches[job["slot"]],retained_native_prefix_equivalence=prefix,
            matched_mean=dict(slot=mean["slot"],native_audit=detailed[mean["slot"]],native_engineering_release=batches[mean["slot"]]),
            no_input=dict(slot=zero["slot"],native_audit=detailed[zero["slot"]],native_engineering_release=batches[zero["slot"]])))
    declared={(t["model"],t["size"],t["fraction"]) for t in native["e4_triggers"]}
    observed={(r["model"],r["size"],r["fraction"]) for r in rows}
    if len(declared)!=52 or observed!=declared or len(rows)!=624 or len(means)!=52 or len(zeros)!=12:
        raise ValueError("E4 numerical projection requires all52 literal triggers and all688 native tasks")
    qualified=[]
    for model,size,fraction in sorted(declared):
        group=[r for r in rows if (r["model"],r["size"],r["fraction"])==(model,size,fraction)]
        if len(group)!=12 or len({r["seed"] for r in group})!=12:raise ValueError("E4 condition lost a declared seed outcome")
        qualified.append(dict(model=model,C_label=size,f=fraction,native_seed_outcomes=12,
            original_seed_continuations=sum(r["population"]=="original_seed_continuation" for r in group),
            mean_slot=group[0]["matched_mean"]["slot"],zero_slot=group[0]["no_input"]["slot"],
            numerical_integrity_only=True,scientific_long_arm_qualification=False))
    report=dict(generated_utc=C.now(),passed=True,complete=True,all_preselected_arms_complete=True,
        all_required_prefix_checks_passed=True,all_required_controls_passed=True,qualified_conditions=qualified,
        numerical_integrity_only=True,scientific_long_arm_qualification=False,
        scientific_scope="Complete native evolution/control/prefix integrity; stationarity, forcing, persistence and hypothesis qualification remain separate statistical gates.",
        native_manifest_sha256=C.sha(run/"manifest.json"),independent_native_master=dict(path=str(audit),sha256=C.sha(audit)),
        noise_workers=len(rows),mean_controls=len(means),no_input_controls=len(zeros),cases=rows)
    C.write(out,report);return report


def complete_e1_e2(ready,late,e2,out):
    """Release the next declared native arms after every E1/E2 numerical gate."""
    ready=Path(ready).resolve();late=Path(late).resolve();e2=Path(e2).resolve();out=Path(out).resolve()
    if out.exists() or not out.is_relative_to(e2/"results"):raise ValueError("Choose a fresh E1/E2 completion receipt")
    stages=[];identities=[];gates={};basals={};native_sources=None
    for run,stage,audit_path in ((ready,"E1",ready/"results/complete_native_independent_audit_20261020.json"),
            (late,"E1",late/"results/independent_native_audit_20261020/status.json"),
            (e2,"E2",e2/"results/complete_native_independent_audit_20261020.json")):
        manifest=C.verify(run);status=C.read(run/"results/status.json");census=C.read(run/"results/case_census.json");audit=C.read(audit_path)
        if manifest["stage"]!=stage or not status["complete"] or status["missing_slots"] or status["failed_slots"]:
            raise ValueError("Every declared E1/E2 native task must complete")
        if status["manifest_sha256"]!=C.sha(run/"manifest.json") or census["native_run_manifest_sha256"]!=C.sha(run/"manifest.json"):
            raise ValueError("Native E1/E2 collection refers to another manifest")
        if not audit["passed"] or audit.get("full_native_run_covered") is not True or audit["native_manifest_sha256"]!=C.sha(run/"manifest.json"):
            raise ValueError("Every E1/E2 native task needs its independent field/forcing audit")
        if len(census["cases"])!=len(manifest["jobs"]) or {r["slot"] for r in audit["receipts"]}!=set(range(len(manifest["jobs"]))):
            raise ValueError("E1/E2 native census/audit coverage is incomplete")
        for job,case in zip(manifest["jobs"],census["cases"]):
            folder=run/"workers"/f"{job['slot']:04d}";result=C.read(folder/"result.json")
            if any(case[key]!=job[key] for key in ("model","size","fraction","kind","seed","sigma","new_gaussian_path")):
                raise ValueError("Native E1/E2 census changed a trajectory identity")
            for key,name in (("trace","trace.npz"),("profile","primary.npz")):
                if case[key]!=str(folder/name) or case[key+"_sha256"]!=result["artifact_sha256"][name]:raise ValueError("Native E1/E2 census changed a scientific observation")
        for record in audit["receipts"]:
            if C.sha(record["path"])!=record["sha256"]:raise ValueError("Independent E1/E2 native receipt changed")
            detailed=C.read(record["path"])
            if not detailed["passed"]:raise ValueError("Independent E1/E2 native receipt failed")
            check_snapshots(detailed["hashed_native_file_snapshots"])
        copied={k:v for k,v in manifest["source_sha256"].items() if k!="frozen_sources/figure4_campaign_20261020.py"}
        if native_sources is None:native_sources=copied
        elif copied!=native_sources:raise ValueError("E1/E2 native model/solver source closures differ")
        for basal in manifest["basals"]:
            key=(basal["model"],basal["size"])
            if key in basals and basals[key]!=basal["input_sha256"]:raise ValueError("E1/E2 native basals differ")
            basals[key]=basal["input_sha256"]
        identities.append(manifest["jobs"])
        for path in (run/"results/status.json",run/"results/case_census.json",audit_path):gates[str(path)]=C.sha(path)
        stages.append(dict(run=str(run),stage=stage,workers=len(manifest["jobs"]),manifest_sha256=C.sha(run/"manifest.json"),
            native_status=str(run/"results/status.json"),native_census=str(run/"results/case_census.json"),independent_native_audit=str(audit_path)))
    base_gate=ready/"inputs/e0_interim_gate.json";release=C.read(base_gate)
    if release.get("e0_interim_delivered") is not True or release.get("e1_strengths_fixed") is not True:raise ValueError("Declared E0 delivery and E1 selections are required")
    for path,digest in release["interim_artifacts"].items():
        if C.sha(path)!=digest:raise ValueError("Delivered corrected E0 interim changed")
    key=lambda j:tuple(j[field] for field in ("model","size","fraction","sigma","replicate","kind"))
    actual=identities[0]+identities[1];expected=C.grid("E1",strengths=release["e1_strengths"])
    if len(actual)!=864 or len({key(j) for j in actual})!=864 or {key(j) for j in actual}!={key(j) for j in expected}:
        raise ValueError("E1 roots do not form the exact disjoint 864-cell declared campaign")
    if len(identities[2])!=468 or {key(j) for j in identities[2]}!={key(j) for j in C.grid("E2")}:
        raise ValueError("E2 does not cover its complete declared transition/control campaign")
    release.update(generated_utc=C.now(),e1_e2_complete=True,native_stage_gates=gates,native_stages=stages,
        completed_E1_native_workers=864,completed_E2_native_workers=468,numerical_integrity_only=True,scientific_analysis_complete=False,
        scope="All prospectively declared E1/E2 native tasks, exact basals, model sources, saved fields and seed laws qualified. This releases the already declared E3/E4 native arms; scientific statistics remain separate.")
    C.write(out,release);return release


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest="action",required=True)
    prep=sub.add_parser("prepare");prep.add_argument("--run",type=Path,required=True);prep.add_argument("--out",type=Path,required=True);prep.add_argument("--shards",type=int,default=24);prep.add_argument("--tasks",type=int,nargs="+")
    work=sub.add_parser("worker");work.add_argument("--out",type=Path,required=True);work.add_argument("--shard",type=int,required=True)
    collect_command=sub.add_parser("collect");collect_command.add_argument("--out",type=Path,required=True)
    merge_command=sub.add_parser("merge");merge_command.add_argument("--run",type=Path,required=True);merge_command.add_argument("--parts",type=Path,nargs="+",required=True);merge_command.add_argument("--out",type=Path,required=True)
    gates_command=sub.add_parser("e4-gates");gates_command.add_argument("--run",type=Path,required=True);gates_command.add_argument("--audit",type=Path,required=True);gates_command.add_argument("--out",type=Path,required=True)
    stages_command=sub.add_parser("complete-e1-e2");stages_command.add_argument("--ready",type=Path,required=True);stages_command.add_argument("--late",type=Path,required=True);stages_command.add_argument("--e2",type=Path,required=True);stages_command.add_argument("--out",type=Path,required=True)
    one=sub.add_parser("one");one.add_argument("--run",type=Path,required=True);one.add_argument("--slot",type=int,required=True)
    args=parser.parse_args()
    if args.action=="prepare":result=prepare(args.run,args.out,args.shards,args.tasks)
    elif args.action=="worker":result=worker(args.out.resolve(),args.shard)
    elif args.action=="collect":result=collect(args.out.resolve())
    elif args.action=="merge":result=merge(args.run,args.parts,args.out)
    elif args.action=="e4-gates":result=e4_gates(args.run,args.audit,args.out)
    elif args.action=="complete-e1-e2":result=complete_e1_e2(args.ready,args.late,args.e2,args.out)
    else:result=audit_worker(args.run,args.slot)
    print(__import__("json").dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
