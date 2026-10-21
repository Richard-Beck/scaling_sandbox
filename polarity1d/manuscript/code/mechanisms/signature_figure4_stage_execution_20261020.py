"""Launch declared read-only stage analyses after complete native qualification."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
import signature_figure4_stage_analysis_20261020 as S


def launch(controller,plan_path):
    controller=controller.resolve();plan_path=plan_path.resolve();freeze=S.E.read_json(controller/'freeze.json')
    for path,digest in freeze['source_sha256'].items():S.E.require_hash(path,digest)
    actual=str(Path(__file__).resolve())
    if actual not in freeze['source_sha256']:raise ValueError('Use the actual frozen stage-execution source')
    if str(plan_path) not in freeze['plan_sha256']:raise ValueError('Only a predeclared frozen execution plan may launch')
    S.E.require_hash(plan_path,freeze['plan_sha256'][str(plan_path)]);plan=S.E.read_json(plan_path);out=Path(plan['out']).resolve()
    if plan.get('mode')=='original_seeds':
        S.prepare_original_seed_sensitivity(out,Path(plan['parent']))
    else:
        for run in plan['native_roots']:S.E.require_hash(Path(run['path'])/'manifest.json',run['manifest_sha256'])
        S.E.require_hash(plan['gate']['path'],plan['gate']['sha256'])
        rule=plan.get('forced_rule')
        if rule:S.E.require_hash(rule['path'],rule['sha256'])
        S.prepare(out,plan['stage'],Path(plan['e0']),[Path(r['path']).resolve() for r in plan['native_roots']],Path(plan['gate']['path']),S.E.read_json(rule['path']) if rule else None)
    S.E.write_json(out/'results/prepared_intake.json',dict(passed=True,stage=plan['stage'],inventory_sha256=S.E.digest(out/'inventory.json'),
        freeze_sha256=S.E.digest(out/'freeze.json'),execution_plan_path=str(plan_path),execution_plan_sha256=S.E.digest(plan_path),
        execution_source_freeze_sha256=S.E.digest(controller/'freeze.json'),no_native_integrations=True))
    S.verify(out)
    arms=S.E.read_json(out/'results/arms.json')['arms'];submitted=[]
    receipt=out/'results/pipeline_submission.json'
    def submit(role,script,dependency=None,memory=None):
        command=['sbatch','--parsable','--job-name=F4'+plan['stage']+'-'+role]
        if dependency:command+=['--dependency=afterok:'+':'.join(dependency)]
        if memory:command+=['--mem='+str(memory)+'G']
        command.append(str(script));result=subprocess.run(command,check=True,capture_output=True,text=True)
        job=result.stdout.strip().split(';')[0]
        if not job.isdigit():raise ValueError('SLURM did not return an unambiguous job ID')
        submitted.append(dict(role=role,job_id=job,script=str(script),script_sha256=S.E.digest(script),command=command))
        S.E.write_json(receipt,dict(stage=plan['stage'],plan_sha256=S.E.digest(plan_path),inventory_sha256=S.E.digest(out/'inventory.json'),submissions=submitted,all_arms_submitted=False))
        return job
    def script(path,commands,memory):
        path.write_text('\n'.join(['#!/bin/bash','#SBATCH --partition=red','#SBATCH --qos=small','#SBATCH --cpus-per-task=1',f'#SBATCH --mem={memory}G','#SBATCH --time=01:00:00',f'#SBATCH --output={out}/logs/{path.stem}_%j.log','set -euo pipefail','export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1']+commands)+'\n')
    (out/'logs').mkdir(exist_ok=True);exports=[];python=shlex.quote(sys.executable);sources=out/'frozen_sources'
    resources=plan.get('resources',{});collector_memory=resources.get('collector_gb',16)
    for arm in arms:
        child=Path(arm['path']);label='s'+str(arm['sigma']).replace('.','p')
        analysis=submit('analysis-'+label,child/'submit_analysis.slurm',memory=resources.get('analysis_gb'))
        audit=submit('audit-'+label,child/'submit_audit.slurm',[analysis],memory=resources.get('audit_gb'))
        path=child/'submit_export_and_audit_aggregate.slurm'
        script(path,[f'{python} -B {shlex.quote(str(sources/"signature_figure4_e0_20261020.py"))} collect --out {shlex.quote(str(child))}',
                     f'{python} -B {shlex.quote(str(sources/"signature_figure4_audit_20261020.py"))} collect --out {shlex.quote(str(child))}'],collector_memory)
        exports.append(submit('export-'+label,path,[audit]))
    path=out/'submit_stage_collect.slurm'
    script(path,[f'{python} -B {shlex.quote(str(sources/"signature_figure4_stage_analysis_20261020.py"))} collect --out {shlex.quote(str(out))}'],collector_memory)
    final=submit('collect',path,exports)
    S.E.write_json(receipt,dict(stage=plan['stage'],plan_sha256=S.E.digest(plan_path),inventory_sha256=S.E.digest(out/'inventory.json'),submissions=submitted,all_arms_submitted=True,final_collect_job_id=final,
        declared_condition_groups={str(arm['sigma']):arm['condition_groups'] for arm in arms},no_native_integrations=True))
    return dict(prepared=True,stage=plan['stage'],out=str(out),final_collect_job_id=final,submissions=submitted)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--controller',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);args=p.parse_args()
    print(json.dumps(launch(args.controller,args.plan),indent=2))
