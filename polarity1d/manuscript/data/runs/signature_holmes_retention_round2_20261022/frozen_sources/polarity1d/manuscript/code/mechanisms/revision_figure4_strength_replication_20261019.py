"""Prepare/check/dispatch bounded replication of three existing strength levels.

No native integration in prepare/check. Frozen pilot scientific functions are
copied byte-for-byte. Submit requires manifest-specific independent acceptance
and all four prospective numerical/observation receipts; no automatic launch.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
MS=ROOT/'polarity1d/manuscript'
PILOT=MS/'data/runs/revision_noise_strength_broad_20261019'
DEFAULT=MS/'data/runs/revision_noise_strength_replication_20261019'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
PILOT_MANIFEST='795023e22468754dd481e00d09a97bd6d371b15f88b8ea0b75c3a4c55e9c4dce'
SCIENTIFIC_DRIVER='4676276ed5eba25c14f214cf02f2f84fa137f708d06d981ebabad2af4a3ad59f'
MODELS=('wave_pinning','goryachev','holmes_model3','spring_ring')
LEVELS=(.003,.03,.1)
REF=dict(wave_pinning=.1,goryachev=6.4,holmes_model3=1.,spring_ring=4.)
REVIEW=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/supporting/metadata'


def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp')
    q.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');q.replace(p)
def csvread(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def csvwrite(p,rr):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rr for k in r)));w.writeheader();w.writerows(rr)
def engine(run):
    path=Path(run)/'frozen_sources/campaign.py';assert sha(path)==SCIENTIFIC_DRIVER
    spec=importlib.util.spec_from_file_location('unchanged_pilot_engine',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


PLAN='''# Bounded common-strength replication: frozen, not submitted

Replicate f=0.003,0.03,0.1 from two to twelve independent seeds at every
original F4 model/size: WP/Goryachev/Holmes C20/40/60; Spring C0=1/2/3.
This adds360 noisy workers in36 groups, with no new mean/no-input controls.
All576 scientific pilot/retained cases remain accessible on exact hash proofs.
After successful execution there are936 scientific cases; the other four
new broad endpoints remain exploratory two-seed levels. No new strength,
correlation, sigma, horizon, basal, native solver or observation law is introduced.

The positive lognormal input mean and SD continue scaling together. Means are
the declared ensemble means (model reference times f), without path-mean
normalization. References are WP0.1, Gory6.4, Holmes1, Spring4. Chemical noise
retains ell1.5,tau50,knots5,coupling2.5; Spring retains basal-physical ell0.15,
tau5,finite knots/coupling0.125,sigma1 and normalized material interpolation.
Spring physical M0/M1 use actual C(t); model clocks are not calibrated.

Chemical horizon32000/scalar5/profile40/full1000 and Spring3200/0.5/4/100
are unchanged. Every job uses the original standard solver from the accepted
pilot; no dense/tight companion is counted as an independent replicate.
Actual functions are copied unchanged from the frozen first-stage campaign.
Proposed seeds are pilot replicate indices2..11, preserving existing0/1.

Selection is common across all models/sizes:0.003 tests memory loss preceding
frequent multiple regions,0.03 resolves the observed organization changes and
large seed/window spreads,0.1 checks higher occupancy versus low-contrast
responses. Original0.001/0.01 already have12 seeds. Grid changes are not phases;
region counts do not establish competition, site lineage or event rates.

Retain M0,q1,R=AQ per frame, contrast/amplitude qualification, native and
basal-relative prominence, raw/qualified multiple-region occupancy, sampled
positive-run spans/censoring, vector/direction memory and eligible coverage.
Analyze full/first/last halves at the same lag checkpoints; the midpoint is
shared. Resample whole independent seeds for intervals after execution.
Memory direction pairs require both endpoint R>=0.1; sparse support remains
unavailable. More seeds estimate finite-window variability, not stationarity.

Resource proposal: QoSsmall,1CPU,3GiB,2h cap per ordinary worker, no throttle.
The measured first two paths at these36 conditions predict63.1CPUh and3.83GB
allocated native output;30% storage reserve4.98GB. These are estimates rather
than upper bounds. Guards retain finite fields, nonnegative pools/conservation,
positive physical C and immutable source/input checks. Failures are errors,
not scientific absence or censoring; do not automatically resubmit.

Preparation is not submission. Submit requires local source/reuse checks,
independent preflight for this exact manifest and all four first-stage
upper-strength numerical/observation acceptance receipts, including Spring.
If Spring observation gates require a protocol change, this payload cannot
be silently changed or launched: prepare a fresh reviewed version instead.
Root reviews the concrete payload before any launch. No native job has been
authorized for automatic submission by this file.
'''


def prepare(run):
    run=Path(run).resolve();assert not run.exists() and run.is_relative_to(MS/'data/runs')
    assert sha(PILOT/'manifest.json')==PILOT_MANIFEST
    pm=read(PILOT/'manifest.json')
    for sub in ('inputs','frozen_sources','protocol','workers','logs','results'):(run/sub).mkdir(parents=True)
    sources={}
    for rel,digest in pm['source_sha256'].items():
        source=PILOT/rel;assert sha(source)==digest;target=run/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);sources[rel]=digest
    shutil.copy2(__file__,run/'frozen_sources/replication_driver.py');sources['frozen_sources/replication_driver.py']=sha(__file__)
    (run/'protocol/PLAN.md').write_text(PLAN)
    shutil.copy2(PILOT/'protocol/PLAN.md',run/'protocol/PILOT_PLAN.md')
    for b in pm['basals']:
        assert sha(PILOT/b['input'])==b['input_sha256'];shutil.copy2(PILOT/b['input'],run/b['input'])
    retained=[]
    for old in csvread(PILOT/'results/exact_reuse_cases.csv'):
        root=Path(old['legacy_run']);status=root/'workers'/f"{int(old['legacy_slot']):04d}"/'status.json'
        retained.append(dict(model=old['model'],size=float(old['size']),kind=old['kind'],fraction=float(old['fraction']),seed=int(old['seed']),
            source_run=str(root),source_slot=int(old['legacy_slot']),source_manifest_sha256=sha(root/'manifest.json'),status=str(status),status_sha256=old['status_sha256'],
            input_sha256=old['input_sha256'],trace=old['trace'],trace_sha256=old['trace_sha256'],primary=old['primary'],primary_sha256=old['primary_sha256'],origin='original_retained'))
    for j in pm['jobs']:
        if j['variant']!='standard':continue
        folder=PILOT/'workers'/f"{j['slot']:04d}";s=read(folder/'status.json')
        assert s['complete'] and s['status']=='complete' and s['producer']=='terminal' and s['job']==j and s['manifest_sha256']==PILOT_MANIFEST
        retained.append(dict(model=j['model'],size=j['size'],kind=j['kind'],fraction=j['fraction'],seed=j['seed'],source_run=str(PILOT),source_slot=j['slot'],
            source_manifest_sha256=PILOT_MANIFEST,status=str(folder/'status.json'),status_sha256=sha(folder/'status.json'),input_sha256=j['input_sha256'],
            trace=str(folder/'trace.npz'),trace_sha256=s['output_sha256']['trace.npz'],primary=str(folder/'primary.npz'),primary_sha256=s['output_sha256']['primary.npz'],origin='broad_pilot'))
    assert len(retained)==576
    jobs=[]
    for model in MODELS:
        for size in ((1.,2.,3.) if model=='spring_ring' else (20.,40.,60.)):
            for f in LEVELS:
                old=sorted([j for j in pm['jobs'] if j['model']==model and j['size']==size and j['fraction']==f and j['kind']=='noise' and j['variant']=='standard'],key=lambda j:j['replicate'])
                assert len(old)==2 and old[1]['seed']==old[0]['seed']+1 and not any(j['dense_observation'] for j in old)
                for rep in range(2,12):jobs.append(dict(old[0],slot=len(jobs),seed=old[0]['seed']+rep,replicate=rep))
    assert len(jobs)==360
    m=dict(generated_at=now(),experiment='Bounded common-strength replication',status='frozen_not_submitted',
        plan=dict(common_strengths=list(LEVELS),models=4,sizes_per_model=3,additional_seeds_per_group=10,target_seeds=12,new_noisy_cases=360,retained_cases=576,total_scientific_cases_after=936,
                  sigma=1.,new_controls=0,new_basals=0,new_correlations=False,new_horizons=False,automatic_extension=False,automatic_submission=False),
        basals=copy.deepcopy(pm['basals']),jobs=jobs,source_sha256=sources,native_protocols=copy.deepcopy(pm['native_protocols']),resources=copy.deepcopy(pm['resources']),
        numerical_check=copy.deepcopy(pm['numerical_check']),protocol_sha256=sha(run/'protocol/PLAN.md'),
        reuse=dict(pilot_run=str(PILOT),pilot_manifest_sha256=PILOT_MANIFEST,index='results/exact_reuse_cases.csv',count=576),
        numerical_acceptance_required=dict(manifest_sha256=PILOT_MANIFEST,receipts={model:str(REVIEW/f'FIGURE4_BROAD_NUMERICAL_{model.upper()}_INDEPENDENT_CHECK.json') for model in MODELS}),
        environment=copy.deepcopy(pm['environment']))
    csvwrite(run/'results/exact_reuse_cases.csv',retained);m['reuse']['index_sha256']=sha(run/'results/exact_reuse_cases.csv')
    write(run/'manifest.json',m)
    write(run/'status.json',dict(generated_at=now(),producer='frozen_not_submitted',native_jobs_submitted=0,requested_jobs=360,retained_scientific_cases=576,complete=False))
    return dict(run=str(run),manifest_sha256=sha(run/'manifest.json'),jobs=360,submitted=False)


def numerical(run):
    m=read(Path(run)/'manifest.json');rr=[]
    for model,p in m['numerical_acceptance_required']['receipts'].items():
        path=Path(p);d=read(path) if path.exists() else {}
        passed=d.get('passed') is True and d.get('model')==model and d.get('manifest_sha256')==PILOT_MANIFEST
        rr.append(dict(model=model,path=p,exists=path.exists(),passed=passed,sha256=sha(path) if path.exists() else None))
    return dict(ready=all(r['passed'] for r in rr),receipts=rr)


def check(run):
    run=Path(run).resolve();m=read(run/'manifest.json');e=engine(run)
    for rel,digest in m['source_sha256'].items():assert sha(run/rel)==digest
    assert sha(run/'protocol/PLAN.md')==m['protocol_sha256'] and sha(run/m['reuse']['index'])==m['reuse']['index_sha256']
    for b in m['basals']:
        assert sha(run/b['input'])==b['input_sha256']
        with np.load(run/b['input']) as z:assert e.ah(z['state'])==b['state_sha256']
    retained=csvread(run/m['reuse']['index']);old_seeds={int(r['seed']) for r in retained if r['kind']=='noise'}
    assert len(retained)==576 and len(m['jobs'])==360 and len({j['seed'] for j in m['jobs']})==360
    for r in retained:
        assert sha(r['status'])==r['status_sha256'];s=read(r['status'])
        assert s['status']=='complete' and s['producer']=='terminal' and s['manifest_sha256']==r['source_manifest_sha256']
        assert s['job']['seed']==int(r['seed']) and s['job']['kind']==r['kind']
        assert sha(r['trace'])==r['trace_sha256'] and sha(r['primary'])==r['primary_sha256']
        b,=[b for b in m['basals'] if b['model']==r['model'] and b['size']==float(r['size'])];assert b['input_sha256']==r['input_sha256']
    for j in m['jobs']:
        assert j['seed'] not in old_seeds and j['kind']=='noise' and j['variant']=='standard' and not j['dense_observation']
        assert j['fraction'] in LEVELS and j['mean']==REF[j['model']]*j['fraction'] and 2<=j['replicate']<=11
    errors=[]
    for platform in ('chemical','spring'):
        native=e.native_module(run,platform)
        assert Path(native.__file__).resolve().is_relative_to(run/'frozen_sources'/platform)
        for b in [b for b in m['basals'] if b['platform']==platform]:
            with np.load(run/b['input']) as z:y=z['state'].copy()
            n=b['cells'];q=(np.arange(n)+.5)/n
            for f in LEVELS:
                for shape in (np.ones(n),1+.2*np.cos(2*np.pi*q)):
                    I=REF[b['model']]*f*shape
                    if platform=='chemical':
                        grid=native.ring.RingGrid(b['size'],n);model=native.ring.native_model(b['model'])
                        d=model.reaction(0.,grid.x,y,I)-model.reaction(0.,grid.x,y,np.zeros(n));expected=I*y[1]/(7.5 if b['model']=='holmes_model3' else 1.)
                        err=float(np.max(abs(d[0]-expected)));assert np.max(abs(d[0]+d[1]))<1e-10
                    else:
                        model=native.Closure(b['size'],n);d=model.rhs(y,I)-model.rhs(y,np.zeros(n));expected=I*y[n:-1]
                        err=float(np.max(abs(d[:n]-expected)));assert np.max(abs(d[:n]+d[n:-1]))<1e-10 and d[-1]==0
                    assert err<1e-10;errors.append(err)
    gates=numerical(run)
    result=dict(generated_at=now(),passed=True,manifest_sha256=sha(run/'manifest.json'),native_integrations=0,source_files=len(m['source_sha256']),basals=12,
        exact_retained_cases=576,new_noisy_cases=360,algebraic_uniform_nonuniform_checks=len(errors),maximum_input_channel_error=max(errors),
        numerical_qualification=gates,ready_to_submit=False,reason='Requires exact-manifest independent preflight and root review; numerical gate state retained explicitly')
    write(run/'results/local_preflight.json',result);return result


def submit(run,review):
    run=Path(run).resolve();m=read(run/'manifest.json');r=read(review)
    assert r.get('passed') is True and r.get('manifest_sha256')==sha(run/'manifest.json'),'Independent exact-manifest preflight required'
    local=read(run/'results/local_preflight.json');assert local['passed'] and local['manifest_sha256']==sha(run/'manifest.json')
    gates=numerical(run);assert gates['ready'],'All four upper-strength numerical/observation acceptances required before submission'
    assert not (run/'submission.json').exists(),'No automatic resubmission'
    for rel,digest in m['source_sha256'].items():assert sha(run/rel)==digest
    args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=02:00:00','--signal=B:USR1@60',
          '--array=0-359','--job-name=f4-strength-repl','--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err'),
          '--wrap','exec '+shlex.join([PYTHON,str(run/'frozen_sources/replication_driver.py'),'worker','--run',str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"']
    handle=subprocess.check_output(args,text=True).strip().split(';')[0]
    receipt=dict(submitted_at=now(),job_id=handle,manifest_sha256=sha(run/'manifest.json'),preflight_review=str(Path(review).resolve()),preflight_sha256=sha(review),numerical_acceptance=gates,command=args)
    write(run/'submission.json',receipt);return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('prepare','check','worker','submit'));p.add_argument('--run',type=Path,default=DEFAULT);p.add_argument('--slot',type=int);p.add_argument('--review',type=Path);a=p.parse_args()
    if a.command=='prepare':r=prepare(a.run)
    elif a.command=='check':r=check(a.run)
    elif a.command=='worker':r=engine(a.run).worker(a.run,a.slot)
    else:assert a.review is not None;r=submit(a.run,a.review)
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
