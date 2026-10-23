"""Only authorized new noise campaign: >=.003, five sizes, 24 -> 96 seeds.

Native equations, forcing, solver, burn-in and horizons are imported unchanged
from the completed round-three campaign. Only durable output differs.
"""
from common import *
import argparse,subprocess,shlex,time,traceback
import numpy as np
import figure4_round3_campaign_20261022 as native
from mode_analysis import sufficient,aggregate_modes

def prepare():
    RUN.mkdir(parents=True,exist_ok=True)
    for d in ('workers','logs','submissions','precision'):(RUN/d).mkdir(exist_ok=True)
    m=read(SOURCE/'manifest.json');conditions=[c for c in m['conditions'] if c['fraction']>=.003]
    frames={p:int(round(v['horizon']/(2*v['scalar_interval'])))+1 for p,v in m['native_protocols'].items()}
    per_seed={p:n*(33*8+4+8+4)+4096 for p,n in frames.items()}
    initial=sum(per_seed[c['platform']]*24 for c in conditions)
    maximum=sum(per_seed[c['platform']]*96 for c in conditions)
    fs=os.statvfs(RUN);free=fs.f_bavail*fs.f_frsize
    plan=dict(source_manifest=str(SOURCE/'manifest.json'),source_sha256=digest(SOURCE/'manifest.json'),
      conditions=conditions,initial_seeds=24,maximum_seeds=96,SE_log_gamma_target=.1,
      modes=33,mode_dtype='complex64',analysis_frames=frames,analysis_cadence={p:v['scalar_interval'] for p,v in m['native_protocols'].items()},
      initial_estimated_bytes=initial,maximum_estimated_bytes=maximum,reserved_bytes=2*maximum,
      filesystem_free_bytes=free,quota_headroom='unknown; shared filesystem free space is not a personal quota',
      storage='Uncompressed byte-payload planning for modes, |P|, time and circumference, not a strict allocated-storage bound; block and metadata overhead require measured allocation. No native profiles saved.',
      generator='unchanged round3 solve; new seeds 202610230600000 + original_condition_index*1000 + replicate',
      native_protocols=m['native_protocols'],source_code_sha256=digest(CODE/'figure4_round3_campaign_20261022.py'))
    if free<2*maximum:raise RuntimeError('Insufficient twofold storage reserve: reduce cadence or mode count before seeds')
    write(RUN/'manifest.json',plan);write(META/'modes_storage_prelaunch.json',plan)
    return plan

def worker(slot,topup=False):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Native integration requires SLURM')
    m=read(RUN/'manifest.json');ci,within=divmod(slot,72 if topup else 24)
    rep=within+24 if topup else within;c=m['conditions'][ci]
    if digest(SOURCE/'manifest.json')!=m['source_sha256']:raise RuntimeError('Source manifest changed')
    folder=RUN/'workers'/f'c{ci:03d}'/f'seed{rep:03d}';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'status.json').exists() and read(folder/'status.json').get('complete'):return
    started=time.monotonic();status=dict(condition=c,replicate=rep,complete=False,job_id=os.environ['SLURM_JOB_ID'])
    try:
        seed=202610230600000+c['condition_index']*1000+rep
        t,p,C,basal,checks,audit=native.solve(SOURCE,c,seed)
        n=p.shape[1];theta=2*np.pi*(np.arange(n)+.5)/n
        # Positive Fourier convention matches existing M1 and heading exactly.
        modes=(p@np.exp(1j*np.outer(theta,np.arange(min(33,n//2+1))))/n).astype('complex64')
        P=np.abs(modes[:,1]/modes[:,0]).astype('float32')
        temp=folder/'modes.tmp.npz'
        np.savez_compressed(temp,time=t.astype('float64'),circumference=C.astype('float32'),modes=modes,P_abs=P)
        temp.replace(folder/'modes.npz')
        audit['mode_only_output']=True;write(folder/'noise_audit.json',audit)
        status.update(complete=True,seed=seed,validation=checks,sha256=digest(folder/'modes.npz'),
                      bytes=(folder/'modes.npz').stat().st_size)
    except Exception as e:status.update(error=str(e),traceback=traceback.format_exc());raise
    finally:status['wall_seconds']=time.monotonic()-started;write(folder/'status.json',status)

def submit():
    m=prepare();receipts=[]
    for model in MODELS:
        slots=[i for i,c in enumerate(m['conditions']) if c['model']==model]
        ranges=','.join(f'{i*24}-{i*24+23}' for i in slots)
        cmd=['sbatch','--parsable','--partition=red','--qos=xxlarge','--cpus-per-task=1','--mem=2G','--time=02:00:00',
             '--array='+ranges,'--job-name=F4r6-mode-'+model[:3],f'--output={RUN}/logs/%A_%a.out',f'--error={RUN}/logs/%A_%a.err',
             '--wrap',shlex.join([PYTHON,'-B',str(Path(__file__).resolve()),'worker','--slot'])+' "$SLURM_ARRAY_TASK_ID"']
        p=subprocess.run(cmd,text=True,capture_output=True)
        receipts.append(dict(command=shlex.join(cmd),returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,
             qos_reason='4800 independent one-core native integrations; xxlarge permits burst throughput, same as round3; no concurrency throttle'))
        write(RUN/'submissions'/f'{model}.json',receipts[-1])
    write(META/'modes_submission.json',receipts)
    ids=[r['stdout'].strip().split(';')[0] for r in receipts if r['returncode']==0]
    if len(ids)==4:
        cmd=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time=02:00:00',
             '--dependency=afterany:'+':'.join(ids),'--job-name=F4r6-precision',f'--output={RUN}/logs/precision_%j.out',
             '--wrap',shlex.join([PYTHON,'-B',str(Path(__file__).resolve()),'precision','--submit-topup'])]
        p=subprocess.run(cmd,text=True,capture_output=True);write(RUN/'submissions/precision.json',dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
    return receipts

def precision(submit_topup=False,early=False):
    # An early precision pass may overlap the dependent final controller in
    # wall time. Serialize allocation decisions, not native worker execution.
    import fcntl
    with (RUN/'precision.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return precision_locked(submit_topup,early)

def precision_locked(submit_topup=False,early=False):
    m=read(RUN/'manifest.json');rows=[];topup=[]
    ledger_path=RUN/'allocation_targets.json'
    ledger=read(ledger_path) if ledger_path.exists() else dict(targets={str(i):24 for i in range(len(m['conditions']))},topup_job_ids=[])
    for ci,c in enumerate(m['conditions']):
        paths=sorted((RUN/'workers'/f'c{ci:03d}').glob('seed*/modes.npz'))
        target=ledger['targets'][str(ci)]
        if len(paths)<target:
            rows.append(dict(condition_id=c['condition_id'],n_seeds=len(paths),target_seeds=target,status='awaiting_allocated_cohort'))
            continue
        ss=[];circumferences=[]
        for path in paths:
            with np.load(path) as z:
                ss.append(sufficient(z['time'],z['modes'],c['model'],z['P_abs']))
                circumferences.append(float(z['circumference'].mean()))
        if not ss:rows.append(dict(condition_id=c['condition_id'],n_seeds=0,status='missing'));continue
        result,curves=aggregate_modes(ss,dict(c,C_mean=float(np.mean(circumferences))),identity=('targeted',ci))
        first=result[1];se=first['SE_log_gamma']
        row=dict(condition_id=c['condition_id'],n_seeds=len(paths),target_seeds=target,SE_log_gamma=se,status=first['status'],gamma=first['gamma'])
        # Unresolved/invalid estimates get the same permitted topup, never a made-up zero rate.
        if target==24 and len(paths)==24 and (not np.isfinite(se) or se>.1):topup.append(ci)
        rows.append(row);csvwrite(RUN/'precision'/f'c{ci:03d}_rates.csv',result)
        csvwrite(RUN/'precision'/f'c{ci:03d}_curves.csv',curves)
    csvwrite(META/'modes_precision.csv',rows)
    write(RUN/'precision/status.json',dict(conditions=rows,topup_conditions=topup,early=early,SE_rule='SE ln gamma >0.1 or unresolved at 24 ->96; cap96'))
    if topup and submit_topup:
        array=','.join(f'{i*72}-{i*72+71}' for i in topup)
        cmd=['sbatch','--parsable','--partition=red','--qos=xxlarge','--cpus-per-task=1','--mem=2G','--time=02:00:00',
             '--array='+array,'--job-name=F4r6-topup',f'--output={RUN}/logs/%A_%a.out',f'--error={RUN}/logs/%A_%a.err',
             '--wrap',shlex.join([PYTHON,'-B',str(Path(__file__).resolve()),'worker','--topup','--slot'])+' "$SLURM_ARRAY_TASK_ID"']
        p=subprocess.run(cmd,text=True,capture_output=True)
        receipt=dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,early=early,conditions=topup,additional_seeds=72*len(topup),allocation_rule='24 ->96 only; seed slots24–95; no duplicate allocated condition')
        write(RUN/'submissions'/f'topup_batch{len(ledger["topup_job_ids"]):02d}.json',receipt)
        if p.returncode==0:
            dep=p.stdout.strip().split(';')[0]
            ledger['topup_job_ids'].append(dep)
            for ci in topup:ledger['targets'][str(ci)]=96
    write(ledger_path,ledger);write(META/'modes_allocation_targets.json',ledger)
    dependencies=[]
    if not early and submit_topup:
        for jid in ledger['topup_job_ids']:
            active=subprocess.run(['squeue','-h','-j',jid,'-o','%A'],capture_output=True,text=True)
            if active.stdout.strip():dependencies.append(jid)
    if dependencies:
        finish=[PYTHON,'-B',str(Path(__file__).resolve()),'precision']
        q=subprocess.run(['sbatch','--parsable','--qos=small','--cpus-per-task=1','--mem=4G','--time=02:00:00','--dependency=afterany:'+':'.join(dependencies),
            '--job-name=F4r6-final',f'--output={RUN}/logs/final_%j.out','--wrap',shlex.join(finish)],text=True,capture_output=True)
        write(RUN/'submissions/final.json',dict(returncode=q.returncode,stdout=q.stdout,stderr=q.stderr,dependencies=dependencies))
    elif not early or submit_topup:
        subprocess.run([PYTHON,'-B',str(OUT/'rebuild_figure4.py')],check=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','submit','worker','precision']);p.add_argument('--slot',type=int);p.add_argument('--topup',action='store_true');p.add_argument('--submit-topup',action='store_true');p.add_argument('--early',action='store_true');a=p.parse_args()
    if a.action=='worker':worker(a.slot,a.topup)
    elif a.action=='precision':precision(a.submit_topup,a.early)
    else:print(globals()[a.action]())
