"""R2 reduced per-cell event arm; native dynamics and frozen rate estimator reused.

Only the arrival hazard and retained output schema differ from Richard's frozen
campaign. Existing per-perimeter records are read only. No parameter/seed search.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, inspect, json, math, os
from pathlib import Path
import shlex, shutil, subprocess, sys
from datetime import datetime, timezone
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript'
SOURCE=MS/'data/runs/poisson_ring_modes_20261027'
RUN=MS/'data/runs/ring_r2_per_cell_20261027'
OUT=MS/'figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/figureR2'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
def now(): return datetime.now(timezone.utc).isoformat()
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');q.replace(p)
def table(p,rows):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def load_module(p,name):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def legacy():
    o=load_module(SOURCE/'frozen_sources/production.py','r2_original_campaign');o.RUN=SOURCE;return o
def segment(run,slot,level):return run/'workers'/f'{slot:03d}'/'segments'/f'{level:02d}'
def available(run,c,final=False):
    levels=read(run/'expected_levels.json');ll=[int(levels[str(s)]) for s in c['slots']]
    if len(set(ll))!=1:return None
    level=ll[0]
    if all((segment(run,s,level)/'status.json').exists() and read(segment(run,s,level)/'status.json').get('complete') for s in c['slots']):return level
    return None
def series(run,slot,level,mode=1):
    return np.concatenate([np.load(segment(run,slot,k)/'u_modes.npy',mmap_mode='r')[:,mode].astype(complex) for k in range(level+1)])

def validate():
    """Check each complete target-model cohort against full cadence and its CI."""
    o=legacy();e=o.engine();m=read(SOURCE/'manifest.json');rows=[]
    for c in m['conditions']:
        if c['model'] not in ('wave_pinning','spring_ring'):continue
        level=available(SOURCE,c)
        if level is None:continue
        job=m['jobs'][c['slots'][0]];xs=[series(SOURCE,s,level) for s in c['slots']]
        raw=np.array([e.acf_raw(x) for x in xs]);g,kind,b=e.half_rate(raw.mean(0)/raw[:,0].mean(),job['frame_interval'])
        lo,hi,los,his,p=o.bootstrap_rates(raw,job['frame_interval'])
        for factor in (1,2,4,8,16,32,64,128):
            thin=np.array([e.acf_raw(x[::factor]) for x in xs]);gg,kk,bb=e.half_rate(thin.mean(0)/thin[:,0].mean(),factor*job['frame_interval'])
            val=gg if gg is not None else bb
            passed=(kk==kind and ((lo<=val<=hi) if kind=='resolved' else np.isclose(val,b,rtol=.002)))
            rows.append(dict(condition=c['index'],model=c['model'],size=c['size'],rate=c['rate_level'],amplitude=c['amplitude_level'],extension_level=level,
              factor=factor,cadence=factor*job['frame_interval'],mode_cutoff=1,full_gamma=g,full_status=kind,full_bound=b,
              full_ci_low=lo,full_ci_high=hi,reduced_gamma=gg,reduced_status=kk,reduced_bound=bb,passes=bool(passed)))
    if not rows:raise RuntimeError('No complete source cohorts')
    choice={}
    for model in ('wave_pinning','spring_ring'):
        q=[r for r in rows if r['model']==model]
        good=[f for f in (1,2,4,8,16,32,64,128) if all(r['passes'] for r in q if r['factor']==f)]
        choice[model]=max(good)
    # The event-window observable is specified in existing native frames. Keep
    # that cadence in a tiny P-only stream while modes use the selected cadence.
    table(RUN/'validation/cadence.csv',rows)
    write(RUN/'validation/selection.json',dict(generated_utc=now(),factors=choice,mode_cutoff=1,
      tested_cohorts=len(rows)//8,source_manifest_sha256=sha(SOURCE/'manifest.json'),
      criterion='Every reduced pooled mode-1 rate inside its full-record exact whole-seed 95% bootstrap CI, same censor status; noncrossings preserve horizon bound',
      P_cadence='Native .5 chemical / .05 Spring, retained for 2/5/20-native-frame event windows',
      cutoff_justification='Mode-1 estimator reads only active mode1; mode0 retained solely to form P. Higher modes and inactive modes never enter gamma1.'))
    print(json.dumps(read(RUN/'validation/selection.json')))

def transformed_engine():
    """Apply audited literal changes to an exact frozen evolution source."""
    original=(SOURCE/'frozen_sources/engine.py').read_text()
    changes=[
      ('return np.r_[delta, job["rho"] * C]','return np.r_[delta, job["cell_rate"]]'),
      ('bottom = sparse.csc_matrix(([job["rho"]], ([0], [physical_dimension-1])), shape=(1, physical_dimension))','bottom = sparse.csc_matrix((1, physical_dimension))'),
      ('for species in ("u", "v")}','for species in ("u",)}'),
      ('saved["v"][next_sample] = coefficients[1]',''),
      ('forcing = Bursts(job, n)','forcing = Bursts(job, n)\n    forcing.log_path = folder / "events.csv"'),
      ('self.events.append((t, peak, center))','self.events.append((t, peak, center))\n        if hasattr(self, "log_path"):\n            with self.log_path.open("a") as event_file:\n                event_file.write(f"{t:.17g},{peak:.17g},{center:.17g}\\n")'),
      ('nm = job["maximum_mode"] + 1','nm = 2'),
    ]
    for a,b in changes:
        if original.count(a)!=1:raise ValueError('Frozen evolution transformation is no longer exact: '+a)
        original=original.replace(a,b)
    return original,changes

def prepare():
    if (RUN/'manifest.json').exists():return
    check=read(RUN/'validation/selection.json');m=read(SOURCE/'manifest.json')
    for part in ('frozen_sources','workers','logs','results','submissions','inputs'):(RUN/part).mkdir(parents=True,exist_ok=True)
    # Basal and native source references are shared read-only; no copy of output.
    (RUN/'native').symlink_to(SOURCE/'native',target_is_directory=True)
    basals=[]
    for b in m['basals']:
        if b['model'] in ('wave_pinning','spring_ring'):
            p=RUN/b['input'];p.symlink_to(SOURCE/b['input']);basals.append(b)
    engine,changes=transformed_engine();(RUN/'frozen_sources/engine.py').write_text(engine)
    shutil.copy2(__file__,RUN/'frozen_sources/campaign.py')
    # The native-frame P stream is generated directly by dense output in the
    # same solve; low modes are reduced immediately afterwards. Its initial
    # schema temporarily holds only active modes 0/1 at native cadence.
    jobs=[];conditions=[]
    for c in m['conditions']:
        if c['model'] not in ('wave_pinning','spring_ring') or c['size'] in (40.,2.):continue
        slots=[]
        for old in c['slots']:
            j=dict(m['jobs'][old]);j.update(source_slot=old,slot=len(jobs),condition_index=len(conditions),maximum_mode=1,supply_arm='per_cell')
            ref=40.
            if j['model']=='spring_ring':
                basal=next(b for b in m['basals'] if b['model']=='spring_ring' and b['size']==2.)
                with np.load(SOURCE/basal['input']) as z:ref=float(z['state'][-1])
            j.update(cell_rate=j['rho']*ref,reference_circumference=ref,mode_stride=check['factors'][j['model']])
            slots.append(j['slot']);jobs.append(j)
        conditions.append(dict(c,index=len(conditions),source_condition=c['index'],slots=slots))
    initial=0;events=0
    for j in jobs:
        frames=int(j['sample_duration']/j['frame_interval'])+1
        # P complex64 + time float64 + optional C float64, low modes complex64.
        initial+=frames*(16+(8 if j['model']=='spring_ring' else 0))+math.ceil(frames/j['mode_stride'])*16
        events+=(j['burn']+j['sample_duration'])*j['cell_rate']
    cap=4*initial+4*events*70+len(jobs)*200000
    manifest=dict(created_utc=now(),jobs=jobs,conditions=conditions,basals=basals,maximum_extension_level=2,
      source_campaign=str(SOURCE.relative_to(ROOT)),source_manifest_sha256=sha(SOURCE/'manifest.json'),
      source_sha256={'frozen_sources/engine.py':sha(RUN/'frozen_sources/engine.py'),'frozen_sources/campaign.py':sha(RUN/'frozen_sources/campaign.py')},
      literal_engine_changes=changes,validation=check,initial_estimated_bytes=int(initial+events*70),
      maximum_estimated_bytes=int(cap),existing_initial_bytes=m['initial_logical_coefficient_bytes'],
      estimated_fraction_of_existing_initial=cap/m['initial_logical_coefficient_bytes'],
      shared_center='All C40 / Spring size2 records are existing-arm records. Spring shared center retains exact original dynamic hazard; outer per-cell sizes use the fixed exact basal size2 physical circumference.',
      seeds='Exactly the corresponding existing outer-size source seeds, unchanged',records='P/time/native cadence, active modes0/1 at validated cadence, Spring C, event log, bounded native restart state; no v or full profile records')
    write(RUN/'manifest.json',manifest);write(RUN/'freeze.json',dict(manifest_sha256=sha(RUN/'manifest.json')))
    write(RUN/'expected_levels.json',{str(j['slot']):0 for j in jobs})
    for name,command in (('worker.sh','worker'),('analysis.sh','analyze')):
        p=RUN/name;p.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+shlex.join([PY,'-B',str(RUN/'frozen_sources/campaign.py'),command,'--slot'])+' "$SLURM_ARRAY_TASK_ID" --level "${R2_LEVEL:-0}"\n');p.chmod(0o750)
    write(RUN/'results/storage_estimate.json',{k:manifest[k] for k in ('initial_estimated_bytes','maximum_estimated_bytes','existing_initial_bytes','estimated_fraction_of_existing_initial')})
    print(json.dumps(read(RUN/'results/storage_estimate.json')))

def engine():
    e=load_module(RUN/'frozen_sources/engine.py','r2_per_cell_engine');e.RUN=RUN;return e
def verify():
    m=read(RUN/'manifest.json');assert sha(RUN/'manifest.json')==read(RUN/'freeze.json')['manifest_sha256']
    for p,h in m['source_sha256'].items():assert sha(RUN/p)==h
    return m
def worker(slot,level):
    m=verify();e=engine();e.worker(slot,level)
    p=segment(RUN,slot,level);status=read(p/'status.json')
    if (p/'P.npy').exists():return
    u=np.load(p/'u_modes.npy',mmap_mode='r');P=(u[:,1].conj()/u[:,0]).astype('complex64')
    np.save(p/'P.npy',P);np.save(p/'heading_time.npy',np.load(p/'time.npy'))
    stride=m['jobs'][slot]['mode_stride'];native_time=np.load(p/'heading_time.npy')
    tick=np.rint((native_time-m['jobs'][slot]['burn'])/m['jobs'][slot]['frame_interval']).astype('int64')
    keep=(tick%stride)==0;thin=np.array(u[keep]);del u
    # Keep native time stream explicitly; mode timestamps are its stride subset.
    np.save(p/'u_modes.npy',thin);np.save(p/'time.npy',native_time[keep])
    status['reduced_record']=dict(mode_stride=stride,mode_cutoff=1,P_native_cadence=m['jobs'][slot]['frame_interval'])
    write(p/'status.json',status)
def sbatch(name,command=None,array=None,dependency=None,level=0,hours='24:00:00'):
    cmd=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time='+hours,'--job-name='+name,
      '--output='+str(RUN/'logs/%A_%a.out'),'--error='+str(RUN/'logs/%A_%a.err')]
    if dependency:cmd+=['--dependency='+dependency]
    if array:cmd+=['--array='+array,'--export=ALL,R2_LEVEL='+str(level),str(RUN/'worker.sh')]
    else:cmd+=['--wrap',shlex.join([PY,'-B',str(RUN/'frozen_sources/campaign.py'),command])]
    r=subprocess.run(cmd,text=True,capture_output=True)
    receipt=dict(generated_utc=now(),command=cmd,returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
    if r.returncode:write(RUN/'submissions/failed.json',receipt);raise RuntimeError(r.stderr)
    jid=r.stdout.strip().split(';')[0];write(RUN/'submissions'/f'{name}_{jid}.json',receipt);return jid
def submit():
    m=verify()
    if (RUN/'submission.json').exists():return
    jid=sbatch('R2_percell',array='0-'+str(len(m['jobs'])-1));control=sbatch('R2_cellcontrol','control',dependency='afterany:'+jid,hours='02:00:00')
    write(RUN/'submission.json',dict(generated_utc=now(),array_job_id=jid,controller_job_id=control,jobs=len(m['jobs']),array_throttle=None));print(read(RUN/'submission.json'))
def control():
    m=verify();e=engine();levels=read(RUN/'expected_levels.json');extend=[];rows=[]
    for c in m['conditions']:
        level=available(RUN,c)
        if level is None:raise RuntimeError('Incomplete per-cell execution: '+str(c['index']))
        j=m['jobs'][c['slots'][0]];raw=np.array([e.acf_raw(series(RUN,s,level)) for s in c['slots']])
        g,k,b=e.half_rate(raw.mean(0)/raw[:,0].mean(),j['frame_interval']*j['mode_stride']);T=j['sample_duration']*2**level
        want=k=='beyond_analysis_lag' or (k=='resolved' and T<8*math.log(2)/g)
        if want and level<2:extend+=c['slots']
        rows.append(dict(condition=c['index'],gamma=g,status=k,bound=b,extension_level=level,record_duration=T,needs_extension=want,at_cap=level==2))
    table(RUN/'results/horizon_screen.csv',rows)
    if extend:
        lev={int(levels[str(s)])+1 for s in extend};assert len(lev)==1;level=lev.pop()
        for s in extend:levels[str(s)]=level
        write(RUN/'expected_levels.json',levels)
        jid=sbatch('R2_cellextend'+str(level),array=','.join(map(str,extend)),level=level)
        cj=sbatch('R2_cellcontrol','control',dependency='afterany:'+jid,hours='02:00:00')
        write(RUN/'results/status.json',dict(state='EXTENDING',generated_utc=now(),array_job_id=jid,controller_job_id=cj,extended_histories=len(extend)))
    else:
        write(RUN/'results/status.json',dict(state='COMPLETE',generated_utc=now(),histories=len(m['jobs'])))
        # Rebuild compact records and figure after execution. Existing producer
        # may still extend; extraction records its precise current status.
        command=shlex.join([PY,'-B',str(MS/'code/mechanisms/rebuild_figureR2_20261027.py'),'all'])
        jid=subprocess.run(['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time=04:00:00','--job-name=R2_rebuild','--output='+str(RUN/'logs/rebuild_%j.out'),'--error='+str(RUN/'logs/rebuild_%j.err'),'--wrap',command],capture_output=True,text=True)
        write(RUN/'results/rebuild_submission.json',dict(returncode=jid.returncode,stdout=jid.stdout,stderr=jid.stderr))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=('validate','prepare','submit','worker','control'));a.add_argument('--slot',type=int);a.add_argument('--level',type=int,default=0);args=a.parse_args()
    if args.command=='worker':worker(args.slot,args.level)
    else:globals()[args.command]()
