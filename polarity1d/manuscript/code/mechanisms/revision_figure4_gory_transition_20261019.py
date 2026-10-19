"""Bounded Gory C60 strength refinement; scientific worker is unchanged.

Native submit requires exact-manifest independent preflight. Parent authorized
this26-case scope; completed936 main evidence is never modified.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,copy,hashlib,importlib.util,json,os,shlex,shutil,subprocess,sys
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
sys.dont_write_bytecode=True
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript';PILOT=MS/'data/runs/revision_noise_strength_broad_20261019'
REPL=MS/'data/runs/revision_noise_strength_replication_scoped_20261019'
DEFAULT=MS/'data/runs/revision_gory_strength_transition_20261019'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
PILOT_SHA='795023e22468754dd481e00d09a97bd6d371b15f88b8ea0b75c3a4c55e9c4dce'
LEVELS=(.015,.02)
PLAN='''# Bounded Goryachev C60 transition refinement

The completed936-case common-grid experiment leaves (.01,.03) unresolved at
C60: the .01 condition retains contrast, while .03 is often low contrast. It
does not exclude a narrower contrast-retaining multiple-region interval.

Add exactly f=.015 and .02, each with12 independent noise seeds and one smooth
ensemble-mean control:26 cases. Use the exact qualified GoryC60 basal, native
chemical parameters/mesh, positive whole-signal lognormal law (sigma1), spatial
length1.5, OU time50, applied knots5/coupling2.5, horizon32000, scalar cadence5,
native profiles40, full state1000 and standard solver tolerances unchanged.
The input reference remains6.4, so means are .096 and .128. No path-mean
normalization, forcing retuning, variance-only arm or new event classifier.

Preserve all936 native cases, sources, strict Spring failure and accepted
ten-metric scope. This refinement is reported separately, with exact old
.01/.03 anchors and C60 no-input context. Primary readouts use full/half
contrast-qualified region occupancy, contrast, R/A/Q, conditional direction
with all-pair coverage; no exclusive phases or lifetime claim. Saved native
profiles and fixed-time/declared medoid examples explain observable changes.

Gory's selected upper-strength standard/tight/dense screen already passed for
the retained numerical scope; no forcing or solver is changed here. Finite,
pool and source/input checks remain per worker. A finite two-point refinement
can reveal a sampled intermediate organization; it cannot exclude narrower
windows between tested values or establish a phase boundary.

New seed namespace202610195800000 +100*level_index +replicate; disjoint from
all completed936 seeds. Distinct strengths use independent realizations.
Controls receive the declared ensemble mean and use no stochastic seed.

Resources:26 one-CPU/3GB small-QoS tasks,1h cap,3300s soft limit, no throttle.
Measured neighboring C60 medians459–715s and max719s imply about4.5CPUh and
0.5GiB native output. No native submission without a manifest-specific
independent passed:true receipt. Parent explicitly authorized this bounded
scope after review; no further author confirmation is required.
'''

def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def core(run):
    path=Path(run)/'frozen_sources/campaign.py';assert sha(path)==read(run/'manifest.json')['source_sha256']['frozen_sources/campaign.py']
    spec=importlib.util.spec_from_file_location('unchanged_gory_strength_worker',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def prepare(run):
    run=Path(run).resolve();assert not run.exists() and run.is_relative_to(MS/'data/runs')
    assert sha(PILOT/'manifest.json')==PILOT_SHA
    parent=read(PILOT/'manifest.json')
    for name in ('inputs','frozen_sources','protocol','workers','logs','results'):(run/name).mkdir(parents=True)
    for rel,digest in parent['source_sha256'].items():
        assert sha(PILOT/rel)==digest;target=run/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(PILOT/rel,target)
    shutil.copy2(__file__,run/'frozen_sources/transition.py')
    sources=dict(parent['source_sha256']);sources['frozen_sources/transition.py']=sha(__file__)
    basal,=[b for b in parent['basals'] if b['model']=='goryachev' and b['size']==60]
    basal=copy.deepcopy(basal);assert sha(PILOT/basal['input'])==basal['input_sha256'];shutil.copy2(PILOT/basal['input'],run/basal['input'])
    cache=read(REPL/'results/partial/case_cache.json');assert len(cache)==936
    oldseeds={d['descriptors'][0]['seed'] for d in cache.values() if d['descriptors'][0]['kind']=='noise'}
    retained=[]
    for case_id,c in cache.items():
        d=next(x for x in c['descriptors'] if x['window']=='last_half')
        if d['model']=='goryachev' and d['size']==60 and (d['fraction'] in (.01,.03) or d['kind']=='basal'):
            retained.append(dict(case_id=case_id,descriptor=d,**{k:c[k] for k in ('trace','profile','trace_sha256','profile_sha256')}))
    assert len(retained)==27
    jobs=[]
    for level,f in enumerate(LEVELS):
        for replicate in range(12):
            seed=202610195800000+100*level+replicate;assert seed not in oldseeds
            jobs.append(dict(**basal,kind='noise',fraction=f,mean=6.4*f,seed=seed,replicate=replicate,variant='standard',dense_observation=False,slot=len(jobs)))
        jobs.append(dict(**basal,kind='mean',fraction=f,mean=6.4*f,seed=0,replicate=0,variant='standard',dense_observation=False,slot=len(jobs)))
    (run/'protocol/PLAN.md').write_text(PLAN)
    manifest=dict(generated_at=now(),experiment='Separate bounded GoryC60 transition refinement',jobs=jobs,basals=[basal],retained_context=retained,
        parent_pilot_manifest_sha256=PILOT_SHA,completed936_manifest_sha256=sha(REPL/'manifest.json'),
        completed936_status_sha256=sha(REPL/'results/partial/status.json'),native_protocols=copy.deepcopy(parent['native_protocols']),
        source_sha256=sources,protocol_sha256=sha(run/'protocol/PLAN.md'),
        resources=dict(ordinary=dict(memory='3G',walltime='01:00:00',soft_seconds=3300),dense=copy.deepcopy(parent['resources']['dense'])),
        root_authorization='Parent explicitly authorized preparation and submission of26 cases after exact independent preflight; no change to936main closure',
        interpretation='Two sampled intermediate strengths; no absence, common-dose or phase-boundary claim')
    write(run/'manifest.json',manifest);return check(run)

def check(run):
    run=Path(run).resolve();m=read(run/'manifest.json')
    for rel,digest in m['source_sha256'].items():assert sha(run/rel)==digest
    assert sha(run/'protocol/PLAN.md')==m['protocol_sha256']
    assert len(m['jobs'])==26 and len(m['basals'])==1
    parent=read(PILOT/'manifest.json');assert m['native_protocols']==parent['native_protocols']
    b=m['basals'][0];assert sha(run/b['input'])==b['input_sha256']
    with np.load(run/b['input']) as z:y=z['state'].copy()
    native=core(run).native_module(run,'chemical');grid=native.ring.RingGrid(60,1200);model=native.ring.native_model('goryachev')
    errors=[]
    for f in LEVELS:
        input=6.4*f*(1+.2*np.sin(2*np.pi*grid.x/60))
        diff=model.reaction(0.,grid.x,y,input)-model.reaction(0.,grid.x,y,np.zeros_like(input))
        errors.append(float(np.max(abs(diff[0]-input*y[1]))));assert errors[-1]<1e-10
    for c in m['retained_context']:
        for key in ('trace','profile'):assert sha(c[key])==c[key+'_sha256']
    noise=[j for j in m['jobs'] if j['kind']=='noise'];assert len(noise)==24 and len({j['seed'] for j in noise})==24
    for level,f in enumerate(LEVELS):
        subset=[j for j in m['jobs'] if j['fraction']==f];assert len(subset)==13
        for j in subset:
            assert j['mean']==6.4*f and not j['dense_observation'] and j['variant']=='standard'
            if j['kind']=='noise':assert j['seed']==202610195800000+100*level+j['replicate']
    result=dict(generated_at=now(),passed=True,manifest_sha256=sha(run/'manifest.json'),native_integrations=0,
        unchanged_worker_sha256=m['source_sha256']['frozen_sources/campaign.py'],frozen_sources=len(m['source_sha256']),
        new_jobs=26,new_noise_seeds=24,matched_mean_controls=2,retained_anchor_context_cases=27,
        maximum_nonuniform_input_channel_error=max(errors),source_basal_and_native_protocol_identity=True,
        independent_preflight_required_before_submit=True)
    write(run/'results/local_preflight.json',result);return result

def submit(run,review):
    run=Path(run).resolve();check(run);r=read(review);assert r.get('passed') is True and r['manifest_sha256']==sha(run/'manifest.json')
    assert not (run/'submission.json').exists()
    command='exec '+shlex.join([PYTHON,str(run/'frozen_sources/campaign.py'),'worker','--run',str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=['sbatch','--parsable','--qos=small','--cpus-per-task=1','--mem=3G','--time=01:00:00','--signal=B:USR1@60','--array=0-25','--job-name=f4-gory-transition','--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err'),'--wrap',command]
    job=subprocess.check_output(args,text=True).strip().split(';')[0]
    result=dict(submitted_at=now(),job_id=job,command=args,manifest_sha256=sha(run/'manifest.json'),preflight_review=str(Path(review).resolve()),preflight_sha256=sha(review),native_jobs=26,root_authorization=read(run/'manifest.json')['root_authorization'])
    write(run/'submission.json',result);return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','check','submit'));p.add_argument('--run',type=Path,default=DEFAULT);p.add_argument('--review',type=Path);a=p.parse_args()
    print(json.dumps(prepare(a.run) if a.action=='prepare' else check(a.run) if a.action=='check' else submit(a.run,a.review),indent=2))

if __name__=='__main__':main()
