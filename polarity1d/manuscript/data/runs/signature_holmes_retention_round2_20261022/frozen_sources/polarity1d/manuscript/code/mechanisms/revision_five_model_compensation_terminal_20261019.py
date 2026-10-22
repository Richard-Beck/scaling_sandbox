"""Guarded terminal derived collection; no native integration or submission."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import os
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'data/runs/revision_five_model_compensation_baseline_v2_20261019'
REPORT=ROOT/'code/mechanisms/revision_five_model_compensation_report_20261019.py'
OUTPUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/five_model_compensation_20261019/baseline'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def now():return datetime.now(timezone.utc).isoformat()
def write(path,data):
    temp=path.with_name(path.name+'.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)


def prepare(run,out):
    folder=run/'derived_terminal_collection';assert not folder.exists();folder.mkdir()
    p=read(run/'provenance.json');paths=[run/'manifest.json',run/'provenance.json',run/'launch.json',REPORT,Path(__file__)]
    paths += [run/'frozen_sources'/r for r in p['frozen_sources']]
    paths += [run/r for r in p['frozen_inputs']]
    guards={str(q):sha(q) for q in paths}
    write(folder/'guard.json',dict(created_at=now(),run=str(run),out=str(out),report=str(REPORT),expected_hashes=guards,
        native_job_id=read(run/'launch.json')['array_job_id'],purpose='Read-only cached terminal native-results collection; no integration.'))
    frozen=folder/'terminal_collector.py';frozen.write_bytes(Path(__file__).read_bytes())
    shell=folder/'worker.sh';shell.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(frozen)+' execute --run '+str(run)+'\n');shell.chmod(0o750)
    job=read(run/'launch.json')['array_job_id']
    command='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:30:00 --dependency=afterany:'+job+' --job-name=F2-baseline-terminal --output='+str(folder/'%j.out')+' --error='+str(folder/'%j.err')+' '+str(shell)
    (folder/'submission_command.txt').write_text(command+'\n')
    print(json.dumps(dict(guard=str(folder/'guard.json'),guard_sha256=sha(folder/'guard.json'),guards=len(guards),
        wrapper=str(shell),wrapper_sha256=sha(shell),command=command),indent=2))


def execute(run):
    folder=run/'derived_terminal_collection';guard=read(folder/'guard.json')
    for rel,digest in guard['expected_hashes'].items():assert sha(Path(rel))==digest,rel
    command=['/home/4473331/.conda/envs/scaling_sandbox/bin/python',guard['report'],'--run',str(run),'--out',guard['out']]
    start=now();result=subprocess.run(command,capture_output=True,text=True)
    (folder/'collector.stdout').write_text(result.stdout);(folder/'collector.stderr').write_text(result.stderr)
    current=Path(guard['out'])/'current'
    status=read(current/'metadata/latest_status.json') if current.exists() else None
    write(folder/'execution.json',dict(started_at=start,finished_at=now(),guard_sha256=sha(folder/'guard.json'),
        scheduler_job_id=os.environ.get('SLURM_JOB_ID'),command=command,returncode=result.returncode,
        actual_native_producer_status=status,actual_generation=str(current.resolve()) if current.exists() else None))
    assert result.returncode==0,result.stderr
    assert status['terminal_scheduler_observed'],'Terminal scheduler state unavailable: result retained, no automatic restart.'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','execute']);parser.add_argument('--run',type=Path,default=DEFAULT);parser.add_argument('--out',type=Path,default=OUTPUT)
    a=parser.parse_args()
    if a.action=='prepare':prepare(a.run.resolve(),a.out.resolve())
    else:execute(a.run.resolve())
