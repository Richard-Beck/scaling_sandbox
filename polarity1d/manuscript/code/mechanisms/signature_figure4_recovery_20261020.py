"""Immutable engineering recovery of explicitly failed Figure 4 native gates.

Successful native arrays are reused exactly. Failed attempts and all original
receipts remain in their original run; a recovery run has its own complete
manifest, adapted reuse receipts, independent native audit and census.
"""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
from pathlib import Path
import shlex
import shutil
import subprocess

import signature_figure4_campaign_20261020 as C


def platform_adapter_sha(path,platform):
    """Ignore an unrelated platform adapter added after a platform pilot froze."""
    tree=ast.parse(Path(path).read_text());worker=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="worker")
    branch=next(n for n in ast.walk(worker) if isinstance(n,ast.If) and isinstance(n.test,ast.Compare)
        and any(isinstance(v,ast.Constant) and v.value=="chemical" for v in n.test.comparators))
    name="OriginalCapBDF" if platform=="chemical" else "spring_native_solve"
    helper=next(n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)
    selected=branch.body if platform=="chemical" else branch.orelse
    return hashlib.sha256((ast.dump(helper,include_attributes=False)+repr([ast.dump(n,include_attributes=False) for n in selected])).encode()).hexdigest()


def freeze_base(original,out):
    original=Path(original).resolve();out=Path(out).resolve();manifest=copy.deepcopy(C.verify(original))
    if out.exists() or not out.is_relative_to(C.MS/"data/runs"):
        raise ValueError("Choose a fresh project-local immutable recovery root")
    for name in ("inputs","frozen_sources","workers","results","logs"):(out/name).mkdir(parents=True)
    for mapping in (manifest["source_sha256"],manifest["input_sha256"]):
        for relative,digest in mapping.items():
            target=out/relative;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(original/relative,target)
            if C.sha(target)!=digest:raise ValueError("Original recovery input/source changed while copying")
    old_driver=out/"frozen_sources/original_figure4_campaign_20261020.py"
    shutil.copy2(original/"frozen_sources/figure4_campaign_20261020.py",old_driver)
    manifest["source_sha256"][str(old_driver.relative_to(out))]=C.sha(old_driver)
    for source,name in ((Path(C.__file__),"figure4_campaign_20261020.py"),(Path(__file__),Path(__file__).name)):
        target=out/"frozen_sources"/name;shutil.copy2(source,target)
        manifest["source_sha256"][str(target.relative_to(out))]=C.sha(target)
    # The recovery CLI imports these canonical names when its frozen copy runs
    # on SLURM; the native worker continues using the separately named driver.
    for source in (Path(C.__file__),Path(__file__).with_name("signature_figure4_native_audit_20261020.py")):
        target=out/"frozen_sources"/source.name;shutil.copy2(source,target)
        manifest["source_sha256"][str(target.relative_to(out))]=C.sha(target)
    manifest.update(generated_utc=C.now(),status="frozen_not_submitted",
        original_attempt=dict(run=str(original),manifest_sha256=C.sha(original/"manifest.json")),
        recovery_scope="Native numerical/source gates only; no response outcome selection or hypothesis changes.",
        native_restart_normalization="Original t0 basal/M1 reference and one complete global-time forcing/observation closure; no checkpoint rebasing.")
    return manifest


def finish_freeze(out,manifest):
    out=Path(out).resolve();C.write(out/"manifest.json",manifest)
    C.write(out/"freeze.json",dict(manifest_sha256=C.sha(out/"manifest.json")))
    C.native_preflight(out)
    script=out/"worker.sh";driver=out/"frozen_sources/figure4_campaign_20261020.py"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec "+
        shlex.join([C.PYTHON,"-B",str(driver),"worker","--out",str(out)])+' --slot "$SLURM_ARRAY_TASK_ID"\n')
    script.chmod(0o750)
    return dict(run=str(out),manifest_sha256=C.sha(out/"manifest.json"),workers=len(manifest["jobs"]))


def prepare_pilot(original,out,slots,*,metadata=None):
    original=Path(original).resolve();out=Path(out).resolve();source=C.verify(original)
    if len(set(slots))!=len(slots) or not slots:raise ValueError("Select distinct actual failed native workers")
    rows=[]
    for slot in slots:
        job=copy.deepcopy(source["jobs"][slot]);status_path=original/"workers"/f"{slot:04d}"/"status.json"
        status=C.read(status_path)
        if (status.get("status")!="failed" or not status.get("native_integration_complete") or "retained_prefix" not in job
                or "check_retained_evolution_prefix" not in status.get("traceback","")):
            raise ValueError("This recovery pilot is limited to a preserved native-prefix validation failure")
        job.update(slot=len(rows),full_stage_slot=slot,
            native_solver_restart_at_original_horizon=job["retained_prefix"]["horizon"],
            native_recovery_source=dict(run=str(original),slot=slot,status=str(status_path),status_sha256=C.sha(status_path),
                native_manifest_sha256=C.sha(original/"manifest.json")))
        rows.append(job)
    manifest=freeze_base(original,out)
    manifest.update(jobs=rows,recovery_pilot=True,counts={k:sum(j["kind"]==k for j in rows) for k in ("noise","mean","basal")},
        resources=dict(cpus=1,memory="4G",qos="small",walltime="01:00:00"))
    manifest.update(metadata or {})
    return finish_freeze(out,manifest)


def prepare_retry(original,out,slots,pilot,recipe):
    original=Path(original).resolve();out=Path(out).resolve();pilot=Path(pilot).resolve();recipe=Path(recipe).resolve()
    source=C.verify(original);qualified=C.verify(pilot);authority=C.read(recipe)
    if authority["native_adapter_sha256"]!=C.sha(C.__file__) or authority["original_manifest_sha256"]!=C.sha(original/"manifest.json"):
        raise ValueError("The frozen recovery adapter/original cohort changed")
    audit_status=pilot/"results/independent_native_audit_20261020/status.json";audit=C.read(audit_status)
    if audit.get("passed") is not True or len(qualified["jobs"])!=1:raise ValueError("Native retry requires its complete platform engineering pilot")
    pilot_job=qualified["jobs"][0];platform=pilot_job["platform"]
    if pilot_job["full_stage_slot"]!={"spring":361,"chemical":228}[platform]:raise ValueError("Unexpected platform recovery qualification case")
    if any(source["jobs"][slot]["platform"]!=platform for slot in slots):raise ValueError("A platform pilot cannot release another native platform")
    expected=platform_adapter_sha(C.__file__,platform)
    if platform_adapter_sha(pilot/"frozen_sources/figure4_campaign_20261020.py",platform)!=expected:
        raise ValueError("Native platform adapter differs from its qualified pilot")
    candidate=qualified_candidate(pilot,0)
    result=prepare_pilot(original,out,slots,metadata=dict(recovery_pilot=False,actual_failure_retry_batch=True,
        resources=dict(cpus=1,memory="8G",qos="small",walltime="08:00:00"),
        prospective_recovery_recipe=dict(path=str(recipe),sha256=C.sha(recipe)),
        platform_engineering_release=dict(platform=platform,native_adapter_platform_sha256=expected,
            pilot=str(pilot),pilot_manifest_sha256=C.sha(pilot/"manifest.json"),
            pilot_full_independent_audit=dict(path=str(audit_status),sha256=C.sha(audit_status)),candidate=candidate)))
    C.write(out/"results/native_engineering_release_20261020.json",dict(generated_utc=C.now(),passed=True,
        manifest_sha256=C.sha(out/"manifest.json"),platform=platform,exact_native_adapter_fingerprint=expected,
        qualified_pilot=candidate,recipe=dict(path=str(recipe),sha256=C.sha(recipe)),
        selected_failed_original_slots=list(slots),scientific_conditions_unchanged=True,prefix_tolerance_unchanged=True))
    return dict(result,manifest_sha256=C.sha(out/"manifest.json"),platform=platform)


def qualified_candidate(run,slot):
    """Require an actually passed independent engineering receipt for reuse."""
    run=Path(run).resolve();manifest=C.verify(run);job=manifest["jobs"][slot]
    folder=run/"workers"/f"{slot:04d}";status=C.read(folder/"status.json");result=C.read(folder/"result.json")
    audit_root=run/"results/independent_native_audit_20261020"
    receipt=audit_root/"workers"/f"{slot:04d}.json"
    audit=C.read(receipt)
    if not status.get("complete") or audit.get("passed") is not True:
        raise ValueError("Selected recovery candidate has no passed actual native engineering outcome")
    if result["job"]!=job or result["manifest_sha256"]!=C.sha(run/"manifest.json") or audit["native_manifest_sha256"]!=C.sha(run/"manifest.json"):
        raise ValueError("Selected recovery candidate manifest binding changed")
    if audit["native_result_sha256"]!=C.sha(folder/"result.json") or audit["native_status_sha256"]!=C.sha(folder/"status.json"):
        raise ValueError("Selected recovery candidate terminal evidence changed")
    from signature_figure4_native_audit_20261020 import check_snapshots,verify_audit
    audit_manifest,_=verify_audit(audit_root);aggregate=C.read(audit_root/"status.json")
    links=[row for row in aggregate["receipts"] if row["slot"]==slot]
    if (not aggregate["full_native_run_covered"] or aggregate["native_manifest_sha256"]!=C.sha(run/"manifest.json")
            or aggregate["audit_manifest_sha256"]!=C.sha(audit_root/"manifest.json")
            or slot not in audit_manifest["tasks"] or len(links)!=1 or links[0].get("passed") is not True
            or Path(links[0]["path"]).resolve()!=receipt or links[0]["sha256"]!=C.sha(receipt)):
        raise ValueError("Selected actual worker is not bound to its unchanged full independent source audit")
    check_snapshots(audit["hashed_native_file_snapshots"])
    if "retained_prefix" in job:
        prefix=C.read(folder/"retained_native_prefix_equivalence.json")
        if prefix.get("passed") is not True or prefix["rule"]!=C.ODE_PREFIX_RULE:raise ValueError("Selected original-seed continuation lost its unchanged prefix gate")
    return dict(run=str(run),slot=slot,native_manifest_sha256=C.sha(run/"manifest.json"),
        result_sha256=C.sha(folder/"result.json"),status_sha256=C.sha(folder/"status.json"),
        artifact_sha256=result["artifact_sha256"],worker_audit=str(receipt),worker_audit_sha256=C.sha(receipt),
        source_audit_manifest=str(audit_root/"manifest.json"),source_audit_manifest_sha256=C.sha(audit_root/"manifest.json"),
        source_audit_status=str(audit_root/"status.json"),source_audit_status_sha256=C.sha(audit_root/"status.json"))


def prepare_selection(original,out,recovery_roots,recipe):
    """Bind all original688 slots to exact passing original/repaired artifacts."""
    original=Path(original).resolve();out=Path(out).resolve();recipe=Path(recipe).resolve()
    source=C.verify(original);authority=C.read(recipe)
    if source["stage"]!="E4" or len(source["jobs"])!=688 or authority["original_manifest_sha256"]!=C.sha(original/"manifest.json"):
        raise ValueError("Recovery selection requires the complete frozen original E4 cohort")
    if authority["native_adapter_sha256"]!=C.sha(C.__file__) or authority["prefix_rule"]!=C.ODE_PREFIX_RULE:
        raise ValueError("The prospectively declared native recovery adapter/rule changed")
    candidates={}
    for root in recovery_roots:
        root=Path(root).resolve();manifest=C.verify(root)
        if manifest["original_attempt"]["manifest_sha256"]!=C.sha(original/"manifest.json"):
            raise ValueError("Recovery candidate belongs to another original cohort")
        for job in manifest["jobs"]:
            path=root/"workers"/f"{job['slot']:04d}"/"status.json"
            if not path.exists() or not C.read(path).get("complete"):continue
            parent=job["full_stage_slot"]
            if parent in candidates:raise ValueError("More than one successful retry exists for an original slot")
            candidates[parent]=qualified_candidate(root,job["slot"])
    jobs=[];selection=[]
    for job in source["jobs"]:
        slot=job["slot"];status_path=original/"workers"/f"{slot:04d}"/"status.json"
        if not status_path.exists():raise ValueError("An original planned native slot has not run")
        status=C.read(status_path)
        if status.get("status") not in ("complete","failed"):raise ValueError("Do not select while an original native slot remains live")
        if status.get("complete"):
            candidate=qualified_candidate(original,slot);role="passing_original"
        else:
            if slot not in candidates:raise ValueError("Failed original slot has no qualified native repair: "+str(slot))
            candidate=candidates[slot];role="qualified_old_cap_restart"
        selected=copy.deepcopy(job);selected["native_artifact_reuse"]=candidate
        actual=C.verify(Path(candidate["run"]))["jobs"][candidate["slot"]]
        for name in ("native_solver_restart_at_original_horizon","native_recovery_source"):
            if name in actual:selected[name]=copy.deepcopy(actual[name])
        jobs.append(selected);selection.append(dict(original_slot=slot,role=role,actual=candidate,
            original_status=dict(path=str(status_path),sha256=C.sha(status_path))))
    manifest=freeze_base(original,out);manifest.update(jobs=jobs,recovery_pilot=False,
        recovery_recipe=dict(path=str(recipe),sha256=C.sha(recipe)),
        actual_native_selection=selection,
        counts={k:sum(j["kind"]==k for j in jobs) for k in ("noise","mean","basal")},
        native_execution_performed_in_immutable_source_runs=True)
    finish_freeze(out,manifest)
    write_reuse_workers(out,jobs,selection)
    C.write(out/"results/complete_original_slot_selection_authority.json",dict(generated_utc=C.now(),
        native_manifest_sha256=C.sha(out/"manifest.json"),all_original_slots=688,
        original_manifest_sha256=C.sha(original/"manifest.json"),recipe_sha256=C.sha(recipe),
        passing_originals=sum(r["role"]=="passing_original" for r in selection),
        actual_restarts=sum(r["role"]=="qualified_old_cap_restart" for r in selection),selection=selection,
        response_outcomes_used_for_selection=False,all_original_attempts_retained=True))
    return dict(run=str(out),manifest_sha256=C.sha(out/"manifest.json"),selected_workers=688)


def write_reuse_workers(out,jobs,selection):
    """Adapt metadata only, after exact source audit snapshots were validated."""
    out=Path(out).resolve()
    for job,row in zip(jobs,selection):
        candidate=row["actual"];source_folder=Path(candidate["run"])/"workers"/f"{candidate['slot']:04d}"
        folder=out/"workers"/f"{job['slot']:04d}";folder.mkdir()
        for name in candidate["artifact_sha256"]:(folder/name).symlink_to(source_folder/name)
        provenance=dict(generated_utc=C.now(),metadata_only_adapter=True,original_selection=row,
            native_reintegration_performed=False,raw_scientific_artifacts_reused_exactly=True)
        C.write(folder/"native_reuse_provenance.json",provenance)
        result=copy.deepcopy(C.read(source_folder/"result.json"));result.update(generated_utc=C.now(),job=job,
            manifest_sha256=C.sha(out/"manifest.json"),metadata_only_reuse_adapter=provenance)
        result["artifact_sha256"]["native_reuse_provenance.json"]=C.sha(folder/"native_reuse_provenance.json")
        C.write(folder/"result.json",result)
        # Raw byte hashes already bind the unchanged actual files and their
        # validated snapshots. The complete selected audit/collector rehashes
        # those bytes independently; metadata preparation need not read them
        # again solely to copy the same values into a new terminal envelope.
        output_hashes=dict(result["artifact_sha256"]);output_hashes["result.json"]=C.sha(folder/"result.json")
        status=copy.deepcopy(C.read(source_folder/"status.json"));status.update(generated_at=C.now(),producer="terminal_immutable_reuse_adapter",job=job,
            manifest_sha256=C.sha(out/"manifest.json"),metadata_only_reuse_adapter=provenance,
            output_sha256=output_hashes)
        C.write(folder/"status.json",status)


def collect_selection(out):
    out=Path(out).resolve();manifest=C.verify(out)
    if not manifest.get("native_execution_performed_in_immutable_source_runs"):
        raise ValueError("This collector is restricted to the declared selected native cohort")
    if (out/"results/status.json").exists() or (out/"results/case_census.json").exists():
        raise ValueError("Never regenerate an existing selected native census/acceptance")
    report=C.collect(out)
    if not report["complete"]:raise ValueError("Selected native cohort is incomplete")
    census=C.read(out/"results/case_census.json")
    for job,case,row in zip(manifest["jobs"],census["cases"],manifest["actual_native_selection"]):
        case.update(original_slot=job.get("full_stage_slot",job["slot"]),native_artifact_reuse=job["native_artifact_reuse"],
            native_recovery_source=job.get("native_recovery_source"),
            native_solver_restart_at_original_horizon=job.get("native_solver_restart_at_original_horizon"),
            execution_selection_role=row["role"],recovery_recipe=manifest["recovery_recipe"],
            scientific_conditions_unchanged=True)
    census["native_execution_selection_authority"]=dict(path=str(out/"results/complete_original_slot_selection_authority.json"),
        sha256=C.sha(out/"results/complete_original_slot_selection_authority.json"))
    C.write(out/"results/case_census.json",census)
    return report


def submit_retry(out):
    """Submit a frozen actual-failure batch and its independent engineering gate."""
    from signature_figure4_native_audit_20261020 import prepare
    out=Path(out).resolve();manifest=C.verify(out)
    if not manifest.get("actual_failure_retry_batch") or (out/"results/submission.json").exists():
        raise ValueError("Submit each declared actual-failure batch exactly once")
    common=["sbatch","--parsable","--partition=red","--qos=small","--cpus-per-task=1"]
    driver=out/"frozen_sources/figure4_campaign_20261020.py"
    preflight=out/"results/forcing_preflight.sh"
    preflight.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"+
        shlex.join([C.PYTHON,"-B",str(driver),"forcing-preflight","--out",str(out)])+"\n"+
        shlex.join([C.PYTHON,"-c","import json; assert json.load(open("+repr(str(out/"results/forcing_preflight.json"))+
            "))['passed'] is True"])+"\n");preflight.chmod(0o750)
    before=subprocess.check_output(common+["--mem=4G","--time=00:15:00","--job-name=F4-E4-retry-preflight",
        "--output="+str(out/"logs/preflight_%j.out"),"--error="+str(out/"logs/preflight_%j.err"),str(preflight)],text=True).strip()
    native=subprocess.check_output(common+["--mem=8G","--time=08:00:00","--dependency=afterok:"+before,
        f"--array=0-{len(manifest['jobs'])-1}","--job-name=F4-E4-actual-prefix-retry",
        "--output="+str(out/"logs/native_%A_%a.out"),"--error="+str(out/"logs/native_%A_%a.err"),str(out/"worker.sh")],text=True).strip()
    audit=out/"results/independent_native_audit_20261020";record=prepare(out,audit,shards=24)
    worker=audit/"worker.sh";auditor=audit/"frozen_sources/signature_figure4_native_audit_20261020.py"
    worker.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec "+
        shlex.join([C.PYTHON,"-B",str(auditor),"worker","--out",str(audit)])+' --shard "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
    after=subprocess.check_output(common+["--mem=8G","--time=00:30:00","--dependency=afterany:"+native,
        f"--array=0-{record['shards']-1}","--job-name=F4-E4-retry-native-audit",
        "--output="+str(out/"logs/audit_%A_%a.out"),"--error="+str(out/"logs/audit_%A_%a.err"),str(worker)],text=True).strip()
    collector=audit/"collect.sh";collector.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"+
        shlex.join([C.PYTHON,"-B",str(auditor),"collect","--out",str(audit)])+"\n"+
        shlex.join([C.PYTHON,"-c","import json; assert json.load(open("+repr(str(audit/"status.json"))+
            "))['passed'] is True"])+"\n");collector.chmod(0o750)
    done=subprocess.check_output(common+["--mem=4G","--time=00:15:00","--dependency=afterany:"+after,
        "--job-name=F4-E4-retry-audit-collect","--output="+str(out/"logs/collect_%j.out"),
        "--error="+str(out/"logs/collect_%j.err"),str(collector)],text=True).strip()
    report=dict(generated_utc=C.now(),manifest_sha256=C.sha(out/"manifest.json"),forcing_preflight_job_id=before,
        native_array_job_id=native,independent_audit_array_job_id=after,independent_audit_collect_job_id=done,
        original_slots=[job["full_stage_slot"] for job in manifest["jobs"]],
        actual_failed_cases_only=True,successful_original_runs_reexecuted=0,native_resources=manifest["resources"])
    C.write(out/"results/submission.json",report);return report


def main():
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest="command",required=True)
    pilot=sub.add_parser("prepare-pilot");pilot.add_argument("--original",type=Path,required=True)
    pilot.add_argument("--out",type=Path,required=True);pilot.add_argument("--slots",type=int,nargs="+",required=True)
    retry=sub.add_parser("prepare-retry");retry.add_argument("--original",type=Path,required=True)
    retry.add_argument("--out",type=Path,required=True);retry.add_argument("--slots",type=int,nargs="+",required=True)
    retry.add_argument("--pilot",type=Path,required=True);retry.add_argument("--recipe",type=Path,required=True)
    selection=sub.add_parser("prepare-selection");selection.add_argument("--original",type=Path,required=True)
    selection.add_argument("--out",type=Path,required=True);selection.add_argument("--recovery-roots",type=Path,nargs="+",required=True)
    selection.add_argument("--recipe",type=Path,required=True)
    collection=sub.add_parser("collect-selection");collection.add_argument("--out",type=Path,required=True)
    submission=sub.add_parser("submit-retry");submission.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if args.command=="prepare-pilot":print(prepare_pilot(args.original,args.out,args.slots))
    elif args.command=="prepare-retry":print(prepare_retry(args.original,args.out,args.slots,args.pilot,args.recipe))
    elif args.command=="prepare-selection":print(prepare_selection(args.original,args.out,args.recovery_roots,args.recipe))
    elif args.command=="collect-selection":print(collect_selection(args.out))
    elif args.command=="submit-retry":print(submit_retry(args.out))


if __name__=="__main__":main()
