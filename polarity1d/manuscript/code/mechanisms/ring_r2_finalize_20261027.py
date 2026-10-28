"""Dependency-driven final R2 refresh, without polling or simulation retries.

Follow the two immutable campaigns' declared controller/finalizer job handles.
Read their status only after a scheduler dependency has completed. A terminal
execution failure is recorded for attention and never treated as censoring.
"""
from __future__ import annotations
import argparse,hashlib,json,os,shlex,shutil,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript'
EXISTING=MS/'data/runs/poisson_ring_modes_20261027'
PER_CELL=MS/'data/runs/ring_r2_per_cell_20261027'
FOLDER=PER_CELL/'release_finalization'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
REBUILD=MS/'code/mechanisms/rebuild_figureR2_20261027.py'
CADENCE_AUDIT=MS/'code/mechanisms/audit_r2_cadence_20261027.py'
ACTIVE={'PENDING','RUNNING','CONFIGURING','COMPLETING','SUSPENDED','RESIZING','REQUEUED'}
def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def status(state,**fields):
    result=dict(generated_utc=now(),state=state,slurm_job_id=os.environ.get('SLURM_JOB_ID'),entrypoint_sha256=sha(__file__),**fields)
    write(FOLDER/'status.json',result);return result
def scheduler(jobs):
    """Record both current scheduler state and accounting exit status."""
    jobs=sorted(set(str(j) for j in jobs if j));joined=','.join(jobs)
    if not jobs:return dict(jobs={},squeue={},sacct={})
    query=subprocess.run(['squeue','-h','-j',joined,'-o','%i|%T|%r'],capture_output=True,text=True)
    live={}
    for line in query.stdout.splitlines():
        parts=line.strip().split('|',2)
        if len(parts)==3:live[parts[0]]=dict(state=parts[1],reason=parts[2])
    acct=subprocess.run(['sacct','-n','-P','-X','-j',joined,'--format=JobIDRaw,State,ExitCode'],capture_output=True,text=True)
    accounts={}
    for line in acct.stdout.splitlines():
        parts=line.strip().split('|')
        if len(parts)>=3 and parts[0] in jobs:accounts[parts[0]]=dict(state=parts[1].split()[0].rstrip('+'),exit_code=parts[2])
    return dict(jobs={j:dict(current=live.get(j),accounting=accounts.get(j)) for j in jobs},
      squeue=dict(returncode=query.returncode,stderr=query.stderr),sacct=dict(returncode=acct.returncode,stderr=acct.stderr))
def submit(dependencies,initial=False):
    deps=sorted(set(str(x) for x in dependencies if x));FOLDER.mkdir(parents=True,exist_ok=True)
    if not deps:raise ValueError('A successor watcher requires live dependency handles')
    command=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time=04:00:00',
      '--job-name=R2_finalize','--dependency=afterany:'+':'.join(deps),
      '--output='+str(FOLDER/'%j.out'),'--error='+str(FOLDER/'%j.err'),
      '--wrap',shlex.join([PY,'-B',str(Path(__file__).resolve()),'run','--predecessors',','.join(deps)])]
    result=subprocess.run(command,capture_output=True,text=True)
    receipt=dict(generated_utc=now(),command=command,dependencies=deps,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
      entrypoint=str(Path(__file__).resolve().relative_to(ROOT)),entrypoint_sha256_at_submission=sha(__file__),
      source_loading='Entrypoint and rebuild are loaded from the latest workspace code only when this dependency job starts; frozen simulation workers are unchanged')
    if result.returncode:
        write(FOLDER/'submission_failure.json',receipt);status('NEEDS_ATTENTION',reason='SLURM submission failed',receipt=receipt);return None
    jid=result.stdout.strip().split(';')[0];receipt['job_id']=jid
    write(FOLDER/f'submission_{jid}.json',receipt);write(FOLDER/'latest_submission.json',receipt)
    if initial:write(FOLDER/'initial_submission.json',receipt)
    print(json.dumps(receipt));return jid
def producer(run):
    p=run/'results/status.json'
    if not p.exists():return dict(path=str(run.relative_to(ROOT)),state='MISSING_STATUS',handles=[])
    obj=read(p);state=obj.get('state','MISSING_STATE')
    if state=='ANALYZING':handles=[obj.get('finalizer_job_id')]
    elif state=='EXTENDING':handles=[obj.get('controller_job_id')]
    elif state=='RUNNING':
        receipt=read(run/'submission.json') if (run/'submission.json').exists() else {}
        handles=[receipt.get('controller_job_id')]
    else:handles=[]
    return dict(path=str(run.relative_to(ROOT)),state=state,status=obj,handles=[str(x) for x in handles if x])
def run(predecessors):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Run the finalizer through SLURM')
    producers=[producer(EXISTING),producer(PER_CELL)]
    handles=[j for p in producers for j in p['handles']]
    receipt=scheduler(handles+predecessors)
    if all(p['state']=='COMPLETE' for p in producers):
        status('REFRESHING',producers=producers,scheduler=receipt)
        # Complete-grid validation supplements the immutable prelaunch choice.
        # It does not change cadence, estimator, inputs, or simulation records.
        audit_folder=PER_CELL/'validation'/('post_launch_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%f'))
        validation=subprocess.run([PY,'-B',str(CADENCE_AUDIT),'--output',str(audit_folder)],capture_output=True,text=True)
        (FOLDER/'cadence_stdout.log').write_text(validation.stdout)
        (FOLDER/'cadence_stderr.log').write_text(validation.stderr)
        summary_path=audit_folder/'summary.json'
        summary=read(summary_path) if summary_path.exists() else {}
        expected_cohorts=sum(c['model'] in ('wave_pinning','spring_ring') for c in read(EXISTING/'manifest.json')['conditions'])
        cadence_pass=bool(validation.returncode==0 and summary.get('complete') and summary.get('cohorts_checked')==expected_cohorts and summary.get('frozen_selection_unchanged') and not summary.get('failures') and not summary.get('pending_cohorts'))
        figure=MS/'figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/figureR2/supporting'
        write(figure/'metadata/record_cadence_finalization.json',dict(generated_utc=now(),passed=cadence_pass,returncode=validation.returncode,expected_cohorts=expected_cohorts,audit=str(audit_folder.relative_to(ROOT)),stderr=validation.stderr[-4000:],policy='Render the complete available figure even when this check fails; record the discrepancy without retuning or silently changing cadence.'))
        for source,target in (('summary.json',figure/'metadata/record_cadence_postlaunch_summary.json'),('snapshot.json',figure/'metadata/record_cadence_postlaunch_snapshot.json'),('cadence.csv',figure/'data/record_cadence_postlaunch_validation.csv')):
            if (audit_folder/source).exists():
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(audit_folder/source,target)
        result=subprocess.run([PY,'-B',str(REBUILD),'all'],capture_output=True,text=True)
        (FOLDER/'rebuild_stdout.log').write_text(result.stdout);(FOLDER/'rebuild_stderr.log').write_text(result.stderr)
        if result.returncode:
            status('NEEDS_ATTENTION',reason='Compact extraction/rebuild failed',returncode=result.returncode,producers=producers,scheduler=receipt)
            raise RuntimeError('R2 final refresh failed; see release_finalization/rebuild_stderr.log')
        # Revision 1: check the revised event-excess algebra and every display
        # encoding on the final data, then refresh scoped links and inventories.
        checks=[
          [PY,'-B',str(MS/'code/mechanisms/audit_ring_revision1_20261028.py'),'--require-complete'],
          [PY,'-B',str(MS/'code/mechanisms/organize_length_dev.py'),'--group','ring_r1r2_20261027/figureR2'],
          [PY,'-B',str(MS/'code/mechanisms/refresh_ring_release_20261028.py')]]
        validation_receipts=[]
        for command in checks:
            checked=subprocess.run(command,capture_output=True,text=True,cwd=ROOT)
            validation_receipts.append(dict(command=command,returncode=checked.returncode,stdout=checked.stdout,stderr=checked.stderr))
            write(FOLDER/'revision1_validation.json',validation_receipts)
            if checked.returncode:
                status('NEEDS_ATTENTION',reason='Revision-1 final-data validation failed',command=command,returncode=checked.returncode,producers=producers,scheduler=receipt)
                raise RuntimeError('Revision-1 final-data validation failed; see revision1_validation.json')
        status('REFRESHED_AWAITING_REVIEW',producers=producers,scheduler=receipt,rebuild_code_sha256=sha(REBUILD),cadence_validation_passed=cadence_pass,cadence_audit=str(audit_folder.relative_to(ROOT)),
          provenance=str((MS/'figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/figureR2/supporting/metadata/provenance.json').relative_to(ROOT)),
          next_step='Independent final-data audit and visual review, followed by root commit/push; this watcher does not commit or publish')
        return
    blockers=[];pending=[]
    for p in producers:
        if p['state']=='COMPLETE':continue
        if p['state'] not in {'RUNNING','EXTENDING','ANALYZING'}:
            blockers.append(dict(producer=p['path'],reason='Terminal or missing producer status: '+p['state']));continue
        if not p['handles']:
            blockers.append(dict(producer=p['path'],reason='No declared controller/finalizer handle'));continue
        for job in p['handles']:
            observed=receipt['jobs'].get(job,{})
            current=observed.get('current');accounting=observed.get('accounting')
            if current and current['state'] in ACTIVE and current.get('reason')!='DependencyNeverSatisfied':pending.append(job)
            elif not current and accounting and accounting['state'] in ACTIVE:pending.append(job)
            else:blockers.append(dict(producer=p['path'],job_id=job,reason='Controller/finalizer is terminal, missing, or dependency-blocked before producer completion',observed=observed))
    if blockers:
        status('NEEDS_ATTENTION',reason='Campaign execution/finalization requires attention; no automatic rerun',blockers=blockers,producers=producers,scheduler=receipt)
        return
    successor=submit(pending)
    if successor:status('WAITING_DEPENDENCY',successor_job_id=successor,dependencies=pending,producers=producers,scheduler=receipt)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=('submit','run'));parser.add_argument('--predecessors',default='21384387,21388801');a=parser.parse_args()
    dependencies=[s for s in a.predecessors.split(',') if s]
    if a.command=='submit':
        if (FOLDER/'initial_submission.json').exists():print(json.dumps(read(FOLDER/'initial_submission.json')))
        else:submit(dependencies,initial=True)
    else:run(dependencies)
