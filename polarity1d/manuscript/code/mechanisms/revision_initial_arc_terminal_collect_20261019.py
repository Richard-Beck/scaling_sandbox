"""One guarded terminal inventory/readout refresh; no integrations or rendering."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import collections
import csv
import hashlib
import json
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)


def scheduler_state(array_id, expected):
    command=["sacct","-j",str(array_id),"--format=JobID%40,State%40,ExitCode","--parsable2","--noheader"]
    try:
        answer=subprocess.run(command,capture_output=True,text=True,timeout=45,check=True)
        tasks={}
        for line in answer.stdout.splitlines():
            fields=line.split("|")
            if len(fields)<3:continue
            job,state,exitcode=(x.strip() for x in fields[:3])
            if job.startswith(str(array_id)+"_") and job[len(str(array_id))+1:].isdigit():
                tasks[job]=dict(state=state,exit_code=exitcode)
        counts=dict(collections.Counter(v["state"] for v in tasks.values()))
        succeeded=sum(v["state"]=="COMPLETED" and v["exit_code"]=="0:0" for v in tasks.values())
        return dict(array_job_id=array_id,queried_at=stamp(),expected_tasks=expected,observed_tasks=len(tasks),state_counts=counts,
                    successful_tasks=succeeded,all_expected_tasks_completed_successfully=len(tasks)==expected and succeeded==expected,
                    task_records=tasks,query_error=None)
    except Exception as error:
        return dict(array_job_id=array_id,queried_at=stamp(),expected_tasks=expected,observed_tasks=0,state_counts={},
                    successful_tasks=0,all_expected_tasks_completed_successfully=False,task_records={},query_error=f"{type(error).__name__}: {error}")


def main(guardpath):
    guardpath=Path(guardpath).resolve();guard=json.loads(guardpath.read_text())
    directory=guardpath.parent;out=Path(guard["derived_output"]);native=Path(guard["native_run"])
    status=dict(started_at=stamp(),source_guard_path=str(guardpath),source_guard_sha256=sha(guardpath),
                frozen_wrapper_sha256=sha(__file__),stage="source_guard",steps=[],error=None,source_guard_passed=False,
                native_inventory_refresh_succeeded=False,derived_readout_refresh_succeeded=False)
    try:
        assert sha(__file__)==guard["frozen_wrapper_sha256"],"Frozen terminal wrapper changed"
        failures=[]
        for name,digest in guard["live_source_sha256"].items():
            path=Path(name)
            if not path.is_file() or sha(path)!=digest:failures.append(name)
        if failures:raise RuntimeError("Live source/input hash guard rejected: "+", ".join(failures))
        status["source_guard_passed"]=True
        for stage,command in (("native_inventory_refresh",[guard["python"],guard["producer_script"],"collect","--run",str(native)]),
                              ("derived_readout_refresh",[guard["python"],guard["collector_script"],"--run",str(native),"--out",str(out),"--terminal"])):
            status["stage"]=stage
            answer=subprocess.run(command,capture_output=True,text=True,timeout=60 if stage=="native_inventory_refresh" else 480)
            status["steps"].append(dict(stage=stage,command=command,return_code=answer.returncode,stdout=answer.stdout,stderr=answer.stderr))
            if answer.returncode:raise RuntimeError(f"{stage} exited {answer.returncode}; see captured output")
            status[stage+"_succeeded"]=True
    except Exception as error:
        status["error"]=f"{type(error).__name__}: {error}"
    status["scheduler"]=scheduler_state(guard["native_array_job_id"],guard["expected_native_tasks"])
    inventory=native/"results/condition_inventory.csv"
    try:
        with inventory.open() as stream:rows=list(csv.DictReader(stream))
        status["terminal_native_inventory_counts"]=dict(collections.Counter(r["status"] for r in rows))
        status["terminal_native_inventory_rows"]=len(rows)
        native_complete=len(rows)==240 and all(r["status"]=="complete" for r in rows)
    except Exception as error:
        status["terminal_native_inventory_read_error"]=f"{type(error).__name__}: {error}"
        native_complete=False
    if status["error"]:
        state="terminal_collection_error"
    elif not native_complete:
        state="terminal_native_outputs_incomplete_or_failed"
    elif not status["scheduler"]["all_expected_tasks_completed_successfully"]:
        state="terminal_readouts_complete_accounting_unresolved_or_failed"
    else:
        state="finished"
    failed_stage=status["stage"] if status["error"] else None
    status.update(finished_at=stamp(),stage="terminal",failed_stage=failed_stage,producer=state,passed=state=="finished")
    atomic_json(directory/"terminal_status.json",status)
    latestpath=out/"supporting/metadata/latest_status.json"
    latest=json.loads(latestpath.read_text()) if latestpath.exists() else dict(requested=240,complete=None)
    latest.update(producer=state,terminal_collection=True,terminal_checked_at=status["finished_at"],
                  terminal_status_path=str(directory/"terminal_status.json"),terminal_status_sha256=sha(directory/"terminal_status.json"),
                  terminal_collection_error=status["error"],terminal_native_inventory_counts=status.get("terminal_native_inventory_counts",{}),
                  terminal_native_scheduler_state_counts=status["scheduler"]["state_counts"],
                  terminal_native_scheduler_successful_tasks=status["scheduler"]["successful_tasks"],
                  derived_readout_refresh_succeeded=status["derived_readout_refresh_succeeded"])
    if status["error"]:
        latest["compact_data_may_contain_partial_export"]=failed_stage=="derived_readout_refresh"
        latest["artifact_sha256"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()
                                  and p.name!="latest_status.json" and not p.name.endswith((".tmp",".lock"))}
    atomic_json(latestpath,latest)
    print(json.dumps({k:v for k,v in status.items() if k not in ("steps","scheduler")},indent=2))
    return 0 if status["passed"] else 1


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--guard",required=True,type=Path)
    sys.exit(main(parser.parse_args().guard))
