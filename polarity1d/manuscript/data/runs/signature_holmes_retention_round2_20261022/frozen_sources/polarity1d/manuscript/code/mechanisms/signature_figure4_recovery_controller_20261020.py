"""Release the complete selected E4 cohort after every actual source audit.

This controller performs metadata preparation and read-only validation only.
It never integrates a native model or promotes scientific long-arm eligibility.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shlex
import shutil
import subprocess

import signature_figure4_campaign_20261020 as C
import signature_figure4_recovery_20261020 as R
import signature_figure4_native_audit_20261020 as A


def script(path,commands):
    path.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"+
        "\n".join(shlex.join(command) for command in commands)+"\n")
    path.chmod(0o750)


def submit(path,name,logs,*,memory="4G",wall="00:15:00",dependency=None,array=None):
    command=["sbatch","--parsable","--partition=red","--qos=small","--cpus-per-task=1","--mem="+memory,"--time="+wall,
        "--job-name="+name,"--output="+str(logs/(name+"_%A_%a.out")),"--error="+str(logs/(name+"_%A_%a.err"))]
    if dependency:command.append("--dependency=afterok:"+":".join(dependency))
    if array:command.append("--array="+array)
    return subprocess.check_output(command+[str(path)],text=True).strip()


def run(plan_path):
    plan_path=Path(plan_path).resolve();home=plan_path.parent
    if C.sha(plan_path)!=C.read(home/"freeze.json")["plan_sha256"]:raise ValueError("Frozen selected-source controller plan changed")
    plan=C.read(plan_path)
    for path,digest in plan["source_sha256"].items():
        if C.sha(path)!=digest:raise ValueError("Frozen selected-source controller code changed")
    original=Path(plan["original_run"]);out=Path(plan["selected_run"])
    if C.sha(original/"manifest.json")!=plan["original_manifest_sha256"] or C.sha(plan["original_attempt_census"]["path"])!=plan["original_attempt_census"]["sha256"]:
        raise ValueError("Original native attempt authority changed")
    original_audit=plan["original_independent_attempt_audit"]
    if C.sha(original_audit["path"])!=original_audit["sha256"]:raise ValueError("Original independent attempt audit was regenerated")
    if C.sha(plan["recipe"]["path"])!=plan["recipe"]["sha256"]:raise ValueError("Prospective recovery recipe changed")
    for source in plan["actual_source_runs"]:
        root=Path(source["run"])
        if C.sha(root/"manifest.json")!=source["manifest_sha256"]:raise ValueError("An actual recovery source manifest changed")
        audit_root=root/"results/independent_native_audit_20261020";A.verify_audit(audit_root)
        if C.sha(audit_root/"manifest.json")!=source["independent_audit_manifest_sha256"]:
            raise ValueError("The prospectively frozen actual source audit manifest changed")
        pinned=source.get("accepted_independent_native_gate")
        if pinned is not None and C.sha(pinned["path"])!=pinned["sha256"]:raise ValueError("An already accepted actual source audit was regenerated")
        audit=C.read(audit_root/"status.json")
        if audit.get("passed") is not True or not audit["full_native_run_covered"]:raise ValueError("Every actual repair source needs its full independent native gate")
    prepared=R.prepare_selection(original,out,[Path(row["run"]) for row in plan["actual_source_runs"]],Path(plan["recipe"]["path"]))
    audit=out/"results/independent_native_audit_20261020";record=A.prepare(out,audit,shards=24)
    native_auditor=audit/"frozen_sources/signature_figure4_native_audit_20261020.py"
    worker=audit/"worker.sh";script(worker,[[C.PYTHON,"-B",str(native_auditor),"worker","--out",str(audit),"--shard","__SHARD__"]])
    worker.write_text(worker.read_text().replace("__SHARD__",'"$SLURM_ARRAY_TASK_ID"'))
    audit_id=submit(worker,"F4-E4-selected-native-audit",out/"logs",memory="8G",wall="00:30:00",array=f"0-{record['shards']-1}")
    audit_collect=audit/"collect.sh";script(audit_collect,[[C.PYTHON,"-B",str(native_auditor),"collect","--out",str(audit)],
        [C.PYTHON,"-c","import json; assert json.load(open("+repr(str(audit/"status.json"))+"))['passed'] is True"]])
    audit_collect_id=submit(audit_collect,"F4-E4-selected-audit-collect",out/"logs",dependency=[audit_id])
    collector=out/"results/collect_selected_native.sh";script(collector,[[C.PYTHON,"-B",str(out/"frozen_sources/signature_figure4_recovery_20261020.py"),"collect-selection","--out",str(out)]])
    collector_id=submit(collector,"F4-E4-selected-native-collect",out/"logs",memory="4G",wall="00:30:00")
    auditor=out/"frozen_sources/signature_figure4_native_audit_20261020.py"
    master=out/"results/complete_native_independent_audit_20261020.json";merger=out/"results/merge_selected_native.sh"
    script(merger,[[C.PYTHON,"-B",str(auditor),"merge","--run",str(out),"--parts",str(audit),"--out",str(master)]])
    master_id=submit(merger,"F4-E4-selected-native-master",out/"logs",dependency=[audit_collect_id,collector_id])
    projection=out/"results/complete_native_numerical_control_prefix_gate_20261020.json";gate=out/"results/project_selected_native_gates.sh"
    script(gate,[[C.PYTHON,"-B",str(auditor),"e4-gates","--run",str(out),"--audit",str(master),"--out",str(projection)]])
    projection_id=submit(gate,"F4-E4-selected-control-prefix",out/"logs",dependency=[master_id])
    report=dict(generated_utc=C.now(),prepared=prepared,controller_plan=dict(path=str(plan_path),sha256=C.sha(plan_path)),
        source_selection_authority=dict(path=str(out/"results/complete_original_slot_selection_authority.json"),sha256=C.sha(out/"results/complete_original_slot_selection_authority.json")),
        native_integrations=0,independent_audit_array_job_id=audit_id,independent_audit_collect_job_id=audit_collect_id,
        native_collector_job_id=collector_id,complete_native_master_job_id=master_id,numerical_control_prefix_gate_job_id=projection_id,
        complete_native_master_path=str(master),numerical_control_prefix_gate_path=str(projection),scientific_qualification=False)
    C.write(out/"results/selected_native_gate_submission_20261020.json",report)
    C.write(home/"submitted_selected_native_gates.json",report)
    return report


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--plan",type=Path,required=True)
    args=parser.parse_args();print(run(args.plan))


if __name__=="__main__":main()
