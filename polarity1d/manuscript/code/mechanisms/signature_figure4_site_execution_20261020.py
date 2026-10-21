"""Engineering-only DAG for additive literal A3 observations of late stages."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
import signature_figure4_site_20261020 as S


def checked(controller,plan_path):
    freeze=S.E.read_json(controller/'freeze.json')
    for path,digest in freeze['source_sha256'].items():S.E.require_hash(path,digest)
    if str(Path(__file__).resolve()) not in freeze['source_sha256']:raise ValueError('Execute the frozen literal-site DAG source')
    S.E.require_hash(plan_path,freeze['plan_sha256'][str(plan_path)])
    plan=S.E.read_json(plan_path);S.E.require_hash(plan['rule_path'],plan['rule_sha256'])
    return plan


def submit(script,role,dependency=None):
    cmd=['sbatch','--parsable','--job-name='+role]
    if dependency:cmd.append('--dependency=afterok:'+str(dependency))
    cmd.append(str(script));r=subprocess.run(cmd,check=True,capture_output=True,text=True);job=r.stdout.strip().split(';')[0]
    if not job.isdigit():raise ValueError('Unambiguous SLURM job ID required')
    return dict(job_id=job,command=cmd,script_sha256=S.E.digest(script))


def script(path,command,memory,array=None):
    lines=['#!/bin/bash','#SBATCH --partition=red','#SBATCH --qos=small','#SBATCH --cpus-per-task=1',f'#SBATCH --mem={memory}G','#SBATCH --time=01:00:00',f'#SBATCH --output={path.parent}/logs/{path.stem}_%A_%a.log']
    if array is not None:lines.append(f'#SBATCH --array=0-{array-1}')
    path.write_text('\n'.join(lines+['set -euo pipefail','export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1',command])+'\n')


def chain(controller,plan_path):
    plan=checked(controller,plan_path);parent=Path(plan['parent']);pipeline=parent/'results/pipeline_submission.json';p=S.E.read_json(pipeline)
    if p['stage']!=plan['stage'] or not p['all_arms_submitted'] or p['plan_sha256']!=plan['parent_execution_plan_sha256'] or p['inventory_sha256']!=S.E.digest(parent/'inventory.json'):
        raise ValueError('The complete prepared normalization DAG is required')
    receipt=controller/f"results/{plan['stage']}_chain.json"
    if receipt.exists():raise FileExistsError('Do not submit a duplicate literal-site chain')
    py=shlex.quote(sys.executable);q=shlex.quote;target=controller/f"launch_{plan['stage']}.slurm"
    script(target,f'{py} -B {q(str(Path(__file__).resolve()))} launch --controller {q(str(controller))} --plan {q(str(plan_path))}',4)
    job=submit(target,'F4site-'+plan['stage']+'-prepare',p['final_collect_job_id'])
    S.E.write_json(receipt,dict(stage=plan['stage'],parent_pipeline_sha256=S.E.digest(pipeline),parent_final_collect_job_id=p['final_collect_job_id'],submission=job))
    return job


def launch(controller,plan_path):
    plan=checked(controller,plan_path);out=Path(plan['out']);S.prepare(out,Path(plan['parent']),Path(plan['rule_path']))
    m=S.E.read_json(out/'manifest.json');n=len(m['groups']);py=shlex.quote(sys.executable);q=shlex.quote;sources=out/'frozen_sources';jobs=[]
    for role,memory,command in [('measure',8,f'{py} -B {q(str(sources/"signature_figure4_site_20261020.py"))} worker --out {q(str(out))} --task "$SLURM_ARRAY_TASK_ID"'),
            ('audit',8,f'{py} -B {q(str(sources/"signature_figure4_site_audit_20261020.py"))} --out {q(str(out))} --task "$SLURM_ARRAY_TASK_ID"'),
            ('collect',4,f'{py} -B {q(str(sources/"signature_figure4_site_20261020.py"))} collect --out {q(str(out))}')]:
        path=out/(role+'.slurm');script(path,command,memory,n if role!='collect' else None)
        record=submit(path,'F4site-'+plan['stage']+'-'+role,jobs[-1]['job_id'] if jobs else None);jobs.append(dict(role=role,**record))
        S.E.write_json(out/'results/submission.json',dict(stage=plan['stage'],manifest_sha256=S.E.digest(out/'manifest.json'),source_execution_freeze_sha256=S.E.digest(controller/'freeze.json'),jobs=jobs,all_submitted=role=='collect',no_native_integrations=True))
    return dict(prepared=True,groups=n,jobs=jobs)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=('chain','launch'));p.add_argument('--controller',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    result=(chain if a.command=='chain' else launch)(a.controller.resolve(),a.plan.resolve());print(json.dumps(result,indent=2))
