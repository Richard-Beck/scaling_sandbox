"""Paired weak-input fast-cytosol campaign for ring Figure R1.

Uses round7.integrate unchanged and the round3 descriptor/exemplar extractors.
Every seed and latent Gaussian field is paired to its retained native run.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, os, shlex, shutil, subprocess, sys, time
from datetime import datetime, timezone
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
CODE=ROOT/'polarity1d/manuscript/code/mechanisms'
BASE=ROOT/'polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022'
RUN=ROOT/'polarity1d/manuscript/data/runs/ring_r1_fast_cytosol_20261027'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
def now(): return datetime.now(timezone.utc).isoformat()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp')
    tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def module(name,p):
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def prepare():
    if (RUN/'manifest.json').exists(): raise ValueError('Campaign already frozen')
    idx=pd.read_csv(BASE/'results/compact_index.csv')
    q=idx[(idx.model=='wave_pinning')&idx.fraction.isin([.003,.005,.01])].sort_values(['fraction','size','replicate'])
    if len(q)!=5760 or not (q.groupby(['size','fraction']).size()==384).all(): raise ValueError('Unexpected native allocation')
    for part in ('frozen_sources','inputs','results','logs'): (RUN/part).mkdir(parents=True,exist_ok=True)
    files=[Path(__file__),CODE/'round7_noise_20261026.py',CODE/'figure4_round3_campaign_20261022.py',CODE/'signature_figure4_regions_20261020.py']
    sources={}
    for p in files:
        target=RUN/'frozen_sources'/p.name;shutil.copy2(p,target);sources[str(target.relative_to(RUN))]=sha(target)
    inputs={}
    for C in (20,30,40,50,60):
        p=BASE/f'inputs/wave_pinning_C{C}.npz';target=RUN/'inputs'/p.name;shutil.copy2(p,target);inputs[str(target.relative_to(RUN))]=sha(target)
    jobs=[]
    for r in q.itertuples(index=False):
        audit=Path(r.compact).with_name('noise_audit.json')
        common=dict(model='wave_pinning',C=float(r.size),f=float(r.fraction),replicate=int(r.replicate),seed=int(r.seed),
                    end=float(r.horizon),burn=float(r.horizon)/2,dt=.05,native_compact=r.compact,native_compact_sha256=r.compact_sha256,
                    gaussian_knots_sha256=read(audit)['gaussian_knots_sha256'])
        jobs.append(dict(slot=len(jobs),D=1000.,solver_check=False,**common))
        if r.fraction==.005 and r.size in (20,60) and r.replicate<8:
            jobs.append(dict(slot=len(jobs),D=10.,solver_check=True,**common))
    # Uncompressed descriptor estimate: exact fixed-ring one-region records.
    descriptor=3201*(18*4+8+2+4+4+1+1)+32768
    exemplar=30*2000*400*2
    estimated=len(jobs)*descriptor+exemplar
    fs=os.statvfs(RUN);available=fs.f_bavail*fs.f_frsize
    if available<2*estimated: raise RuntimeError('Insufficient shared-filesystem margin')
    write(RUN/'manifest.json',dict(created_utc=now(),jobs=jobs,source_sha256=sources,input_sha256=inputs,
        native_manifest=str(BASE/'manifest.json'),native_manifest_sha256=sha(BASE/'manifest.json'),
        native_index_sha256=sha(BASE/'results/compact_index.csv'),
        solver='Unchanged round7 integrate: IMEX Euler dt=.05, native dx=.05, cyclic Thomas implicit diffusion',
        pairing='Exact native per-condition seed, initial state, end/burn and Gaussian knot SHA256; only inactive diffusivity changes',
        extraction='Unchanged round3.extract and round3.exemplar; native histories remain in worker memory',
        storage=dict(planning_bytes_per_seed=descriptor,estimated_bytes=estimated,reservation_bytes=2*estimated,
                     shared_filesystem_available_bytes=available,personal_quota_available_bytes=None,
                     quota_status='quota -sv returned no limit; shared filesystem free space is not personal quota'),
        n_fast=len(q),n_solver_check=16))
    script=RUN/'worker.sh';script.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+shlex.join([PY,'-B',str(RUN/'frozen_sources'/Path(__file__).name),'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"\n');script.chmod(0o750)
    print(json.dumps(dict(jobs=len(jobs),estimated_GiB=estimated/2**30,manifest=str(RUN/'manifest.json'))))
def worker(slot,benchmark=False):
    manifest=read(RUN/'manifest.json');j=manifest['jobs'][slot]
    target=RUN/('benchmark' if benchmark else 'results')/f'{slot:04d}'
    if (target/'status.json').exists() and read(target/'status.json').get('complete'):return
    target.mkdir(parents=True,exist_ok=True)
    for relative,digest in {**manifest['source_sha256'],**manifest['input_sha256']}.items():
        if sha(RUN/relative)!=digest: raise ValueError('Frozen source/input changed: '+relative)
    r7=module('r1_round7',RUN/'frozen_sources/round7_noise_20261026.py')
    r3=module('r1_round3',RUN/'frozen_sources/figure4_round3_campaign_20261022.py')
    helper=module('r1_native_helper',BASE/'frozen_sources/signature_figure4_campaign_20261020.py')
    native=helper.native_module(BASE,'chemical')
    with np.load(RUN/f"inputs/wave_pinning_C{j['C']:g}.npz") as z: state=z['state'].copy()
    basal=state[0].copy();end=100. if benchmark else j['end'];burn=0. if benchmark else j['burn']
    forcing=native.NoiseInput(native.ring.RingGrid(j['C'],state.shape[1]),end,j['seed'],j['f']*.1,
                             **read(BASE/'manifest.json')['native_protocols']['chemical']['noise'])
    knot_sha=hashlib.sha256(np.asarray(forcing.knots,dtype='<f8').tobytes()).hexdigest()
    if not benchmark and knot_sha!=j['gaussian_knots_sha256']:raise ValueError('Noise pairing failed')
    started=time.time()
    profiles,merr,minimum=r7.integrate(state,np.asarray(forcing.knots),j['f']*.1,0,np.array([.1,j['D']]),j['dt'],end,burn,.05)
    if not np.isfinite(profiles).all() or minimum < -1e-8 or merr>1e-7: raise ValueError('Numerical guard failed')
    times=burn+5*np.arange(len(profiles));fields,protocol,summary=r3.extract(times,profiles,j['C'],basal)
    np.savez_compressed(target/'compact.npz',**fields)
    if j['replicate']<2:
        np.savez_compressed(target/'exemplar.npz',**r3.exemplar(fields,profiles))
    write(target/'status.json',dict(job=j,complete=True,generated_utc=now(),elapsed_seconds=time.time()-started,
        gaussian_knots_sha256=knot_sha,paired_noise_verified=not benchmark,minimum_concentration=minimum,
        maximum_pool_relative_error=merr,compact_sha256=sha(target/'compact.npz'),region_protocol=protocol,region_summary=summary,
        slurm_job_id=os.environ.get('SLURM_JOB_ID'),benchmark=benchmark))
    print(json.dumps(dict(slot=slot,seconds=time.time()-started,mass_error=merr,minimum=minimum)),flush=True)
def submit():
    m=read(RUN/'manifest.json')
    cmd=['sbatch','--parsable','--partition=red','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--mem=2G','--time=01:00:00',
         '--array=0-'+str(len(m['jobs'])-1),'--job-name=ring_R1','--output='+str(RUN/'logs/%A_%a.out'),
         '--error='+str(RUN/'logs/%A_%a.err'),str(RUN/'worker.sh')]
    p=subprocess.run(cmd,text=True,capture_output=True)
    write(RUN/'submission.json',dict(generated_utc=now(),command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,
        one_line_command=shlex.join(cmd),qos_justification='Thousands of independent one-core runs; scheduler controls concurrency; no array throttle'))
    print(p.stdout or p.stderr,flush=True)
def collect():
    m=read(RUN/'manifest.json');rows=[];missing=[]
    for j in m['jobs']:
        folder=RUN/'results'/f"{j['slot']:04d}";p=folder/'status.json'
        if not p.exists() or not read(p).get('complete'):missing.append(j['slot']);continue
        status=read(p)
        if sha(folder/'compact.npz')!=status['compact_sha256']:raise ValueError('Compact digest failed')
        rows.append(dict(**j,compact=str(folder/'compact.npz'),compact_sha256=status['compact_sha256'],
                         exemplar=str(folder/'exemplar.npz') if (folder/'exemplar.npz').exists() else '',
                         paired_noise_verified=status['paired_noise_verified']))
    pd.DataFrame(rows).to_csv(RUN/'results/compact_index.csv',index=False)
    write(RUN/'status.json',dict(generated_utc=now(),complete=not missing,completed=len(rows),expected=len(m['jobs']),missing_slots=missing))
    print(json.dumps(dict(completed=len(rows),missing=len(missing))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','worker','submit','collect']);p.add_argument('--slot',type=int);p.add_argument('--benchmark',action='store_true');a=p.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='worker':worker(a.slot,a.benchmark)
    elif a.command=='submit':submit()
    else:collect()
