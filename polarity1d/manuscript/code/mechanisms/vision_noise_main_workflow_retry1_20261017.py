"""Bounded SLURM statistics/report workflow for the frozen complete ensemble.

No simulation, failed-case omission, or estimator modification is permitted.
The prospective estimator source is executed byte-for-byte from the run's
protocol. A separate explicit index records every prescribed case and status.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import traceback

PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"


def now():return datetime.now(timezone.utc).isoformat()
def read(path):return json.loads(Path(path).read_text())
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(2**20),b""):digest.update(block)
    return digest.hexdigest()
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");temporary.replace(path)


def verify_plan(run):
    freeze=read(run/"plan_freeze.json");manifest=read(run/"manifest.json")
    assert len(manifest["jobs"])==243
    counts={kind:sum(j["kind"]==kind for j in manifest["jobs"]) for kind in ("noise","mean","basal")}
    assert counts==dict(noise=216,mean=18,basal=9)
    assert manifest["horizon"]==32000 and manifest["output_interval"]==5 and manifest["primary_interval"]==40
    plan=freeze["analysis"]
    assert plan["max_lag_native"]==8000
    assert plan["bootstrap"]==dict(draws=1000,seed=920731,unit="whole independent complete realization",interval="percentile2.5–97.5; recompute full ratio")
    assert plan["site_descriptors"]["cadence_native"]==5 and plan["site_descriptors"]["source"]=="scalar"
    for record in freeze["files"].values():assert sha(run/record["copy"])==record["copy_sha256"]
    assert sha(run/"configuration.json")==manifest["configuration_sha256"]==freeze["configuration_sha256"]
    for relative,digest in manifest["source_sha256"].items():assert sha(run/"frozen_sources"/relative)==digest
    return freeze,manifest


def complete_index(run,workflow):
    workflow.mkdir(parents=True,exist_ok=True)
    freeze,manifest=verify_plan(run);rows=[]
    for job in manifest["jobs"]:
        folder=run/"workers"/f"{job['slot']:04d}";path=folder/"status.json"
        assert path.exists(),f"Missing prescribed worker {job['slot']}"
        status=read(path)
        assert status.get("producer")=="terminal" and status.get("complete") and status.get("status")=="complete",(job["slot"],status.get("status"))
        assert status["job"]==job and status["native_time"]==32000
        assert status["samples"]==6401 and status["full_samples"]==33
        assert status["manifest_sha256"]==sha(run/"manifest.json")
        rows.append(dict(slot=job["slot"],model=job["model"],circumference=job["circumference"],
                         condition=job["condition"],fraction=job["mean_fraction"],seed=job["seed"],replicate=job["replicate"],
                         archive=str(folder/"trace.npz"),profile_archive=str(folder/"primary.npz"),
                         status="complete",status_path=str(path),status_sha256=sha(path),
                         input=str(run/job["input"]),input_sha256=job["input_sha256"],
                         trace_sha256=status["output_sha256"]["trace.npz"],profile_sha256=status["output_sha256"]["primary.npz"]))
    for model in ("wave_pinning","goryachev","holmes_model3"):
        for C in (20,40,60):
            for fraction in (.001,.01):
                group=[r for r in rows if r["model"]==model and r["circumference"]==C and r["fraction"]==fraction and r["condition"]=="noise"]
                assert len(group)==12 and len({r["seed"] for r in group})==12
                assert {r["replicate"] for r in group}==set(range(12))
    path=workflow/"all_prescribed_cases.csv"
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    return path


def prepare(run,output,report):
    run=Path(run).resolve();output=Path(output).resolve();report=Path(report).resolve()
    verify_plan(run);folder=run/"results/statistics_workflow_retry1"
    review=output/"supporting/metadata/FRESH_FLAT_CENTER_AMENDMENT_REVIEW.md"
    assert review.is_file(),"Fresh numerical amendment review required before workflow freeze"
    amendment=run/"results/flat_center_validation_amendment"
    regression=read(amendment/"REGRESSION.json")
    assert regression["passed"] and regression["validated_cases"]==243 and regression["profile_snapshots"]==194643
    assert not regression["statistical_or_scalar_eligibility_change"] and not regression["tolerance_loosened"]
    assert sha(amendment/"corrected_analysis_source.py")==regression["corrected_sha256"]
    assert sha(run/"protocol/analysis_source.py")==regression["original_frozen_analysis_sha256"]
    folder.mkdir(parents=True,exist_ok=True)
    assert not(folder/"submission.json").exists(),"Workflow already submitted; inspect exact handle instead of launching again"
    for source,name in ((Path(__file__).resolve(),"workflow.py"),(report,"report.py"),(amendment/"corrected_analysis_source.py","amended_analysis.py")):
        target=folder/name
        if target.exists():assert sha(target)==sha(source),f"Prepared source changed; preserve its record before replacement: {target}"
        else:shutil.copy2(source,target)
    record=dict(prepared_at=now(),run=str(run),output=str(output),resources=dict(qos="small",cpus=1,memory="16G",walltime="01:30:00"),
                memory_reason="Frozen estimator retains millions of scalar/CSV dictionaries; 16GB avoids the borderline 8GB peak without changing scientific output",
                sources={name:sha(folder/name) for name in ("workflow.py","report.py","amended_analysis.py")},
                analysis_source=str(folder/"amended_analysis.py"),analysis_sha256=sha(folder/"amended_analysis.py"),
                original_analysis_source=str(run/"protocol/analysis_source.py"),original_analysis_sha256=sha(run/"protocol/analysis_source.py"),
                numerical_regression_sha256=sha(run/"results/analysis_numerical_amendment/REGRESSION.json"),
                flat_center_regression_sha256=sha(amendment/"REGRESSION.json"),
                preserved_failed_job="21126773",
                numerical_review_path=str(review),numerical_review_sha256=sha(review),
                manifest_sha256=sha(run/"manifest.json"),plan_sha256=sha(run/"plan_freeze.json"))
    write(folder/"preparation.json",record)
    return record


def submit(run,dependency=None):
    run=Path(run).resolve();folder=run/"results/statistics_workflow_retry1";prepared=read(folder/"preparation.json")
    assert not(folder/"submission.json").exists(),"Already submitted; inspect exact job"
    if dependency is None:complete_index(run,folder)
    command="exec "+shlex.join(["env","PYTHONDONTWRITEBYTECODE=1","OPENBLAS_NUM_THREADS=1","OMP_NUM_THREADS=1","MKL_NUM_THREADS=1",PYTHON,str(folder/"workflow.py"),"worker","--run",str(run)])
    args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=16G","--time=01:30:00","--job-name=vision-noise-stats","--output="+str(folder/"slurm_%j.out"),"--error="+str(folder/"slurm_%j.err")]
    if dependency is not None:args.append("--dependency=afterok:"+str(dependency))
    args += ["--wrap",command]
    job=subprocess.check_output(args,text=True).strip().split(";")[0]
    record=dict(submitted_at=now(),job_id=job,dependency=dependency,command=args,preparation_sha256=sha(folder/"preparation.json"))
    write(folder/"submission.json",record);return record


def worker(run):
    assert os.environ.get("SLURM_JOB_ID"),"Heavy analysis must use SLURM"
    run=Path(run).resolve();folder=run/"results/statistics_workflow_retry1";prepared=read(folder/"preparation.json")
    for name,digest in prepared["sources"].items():assert sha(folder/name)==digest
    assert sha(prepared["analysis_source"])==prepared["analysis_sha256"]
    assert sha(run/"protocol/analysis_source.py")==prepared["original_analysis_sha256"]
    assert sha(prepared["numerical_review_path"])==prepared["numerical_review_sha256"]
    state=dict(started_at=now(),producer="running",stage="complete-case validation",slurm_job_id=os.environ["SLURM_JOB_ID"],complete=False)
    write(folder/"status.json",state)
    try:
        index=complete_index(run,folder)
        state.update(stage="frozen prospective analysis",case_index_sha256=sha(index));write(folder/"status.json",state)
        command=[PYTHON,prepared["analysis_source"],"--index",str(run/"manifest.json"),"--output",prepared["output"],"--maxlag","8000","--bootstrap","1000","--site-cadence","scalar","--prefix","vision_noise_"]
        subprocess.run(command,check=True)
        state.update(stage="actual-result report");write(folder/"status.json",state)
        subprocess.run([PYTHON,str(folder/"report.py"),"--run",str(run),"--output",prepared["output"]],check=True)
        state.update(complete=True,status="complete",stage="terminal")
    except Exception as error:
        state.update(status="failed",error=f"{type(error).__name__}: {error}",traceback=traceback.format_exc())
        raise
    finally:
        state.update(producer="terminal",finished_at=now());write(folder/"status.json",state)
    return state


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("prepare","submit","index","worker"));parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--output",type=Path);parser.add_argument("--report",type=Path);parser.add_argument("--dependency")
    args=parser.parse_args()
    if args.command=="prepare":result=prepare(args.run,args.output,args.report)
    elif args.command=="submit":result=submit(args.run,args.dependency)
    elif args.command=="index":result=dict(path=str(complete_index(args.run.resolve(),args.run.resolve()/"results/statistics_workflow_retry1")))
    else:result=worker(args.run)
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
