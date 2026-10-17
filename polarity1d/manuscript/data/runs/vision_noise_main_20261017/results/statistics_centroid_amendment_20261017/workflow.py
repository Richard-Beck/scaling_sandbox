"""Statistics/report-only retry for verified native ensemble; no reintegration."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, importlib.util, json, os, shlex, subprocess, traceback
FOLDER=Path(__file__).resolve().parent
RUN=FOLDER.parents[1]
PYTHON="/home/4473331/.conda/envs/scaling_sandbox/bin/python"
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,j):
 p=Path(p);tmp=p.with_suffix(p.suffix+".tmp");tmp.write_text(json.dumps(j,indent=2,allow_nan=False)+"\n");tmp.replace(p)
def prior():
 spec=importlib.util.spec_from_file_location("prior_workflow",FOLDER/"previous_workflow_source.py")
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def validate():
 prepared=read(FOLDER/"preparation.json")
 for name,digest in prepared["source_sha256"].items():assert sha(FOLDER/name)==digest
 assert sha(RUN/"manifest.json")==prepared["manifest_sha256"]
 assert sha(RUN/"protocol/analysis_source.py")==prepared["original_prospective_source_sha256"]
 assert sha(RUN/"results/statistics_workflow/status.json")==prepared["failed_status_sha256"]
 assert sha(RUN/"results/execution_verification.json")==prepared["execution_verification_sha256"]
 assert read(RUN/"results/execution_verification.json")["passed"]
 prior().verify_plan(RUN)
 return prepared

def worker():
 assert os.environ.get("SLURM_JOB_ID"),"Statistics require SLURM"
 prepared=validate();state=dict(started_at=now(),producer="running",stage="complete-case validation",complete=False,job_id=os.environ["SLURM_JOB_ID"])
 write(FOLDER/"status.json",state)
 try:
  index=prior().complete_index(RUN,FOLDER);state.update(stage="unchanged prospective statistics",case_index_sha256=sha(index));write(FOLDER/"status.json",state)
  subprocess.run([PYTHON,str(FOLDER/"analysis.py"),"--index",str(RUN/"manifest.json"),"--output",prepared["output"],"--maxlag","8000","--bootstrap","1000","--site-cadence","scalar","--prefix","vision_noise_"],check=True)
  state.update(stage="actual-result report");write(FOLDER/"status.json",state)
  subprocess.run([PYTHON,str(FOLDER/"report.py"),"--run",str(RUN),"--output",prepared["output"]],check=True)
  state.update(complete=True,status="complete",stage="terminal")
 except Exception as error:
  state.update(status="failed",error=str(error),traceback=traceback.format_exc());raise
 finally:
  state.update(producer="terminal",finished_at=now());write(FOLDER/"status.json",state)
 return state

def submit():
 validate();assert not(FOLDER/"submission.json").exists(),"Already submitted; inspect exact handle"
 prior().complete_index(RUN,FOLDER)
 command="exec "+shlex.join(["env","PYTHONDONTWRITEBYTECODE=1","OPENBLAS_NUM_THREADS=1","OMP_NUM_THREADS=1","MKL_NUM_THREADS=1",PYTHON,str(FOLDER/"workflow.py"),"worker"])
 args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=16G","--time=01:30:00","--job-name=vision-noise-stats-centroid","--output="+str(FOLDER/"slurm_%j.out"),"--error="+str(FOLDER/"slurm_%j.err"),"--wrap",command]
 job=subprocess.check_output(args,text=True).strip().split(";")[0]
 record=dict(submitted_at=now(),job_id=job,command=args,preparation_sha256=sha(FOLDER/"preparation.json"),verified_native_cases=243,native_integrations=0)
 write(FOLDER/"submission.json",record);return record
if __name__=="__main__":
 parser=argparse.ArgumentParser();parser.add_argument("command",choices=("submit","worker"));args=parser.parse_args()
 print(json.dumps(submit() if args.command=="submit" else worker(),indent=2))
