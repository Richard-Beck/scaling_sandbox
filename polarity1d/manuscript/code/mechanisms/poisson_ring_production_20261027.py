"""Frozen 432-history Poisson ring campaign, bounded extensions and inference."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p/"polarity1d/manuscript").is_dir())
MS = ROOT/"polarity1d/manuscript"
TIMING = MS/"data/runs/poisson_ring_modes_timing_20261027"
PRODUCTION = MS/"data/runs/poisson_ring_modes_20261027"
SMOKE = MS/"data/runs/poisson_ring_modes_preflight_20261027"
OUT = MS/"figures/length_mechanisms_20261008/dev/poisson_ring_modes_20261027"
RUN = PRODUCTION
PY = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
ENGINE = None


def now(): return datetime.now(timezone.utc).isoformat()
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p, v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+".tmp");q.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n");q.replace(p)
def table(p, rows):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def engine():
    global ENGINE
    if ENGINE is None:
        p=RUN/"frozen_sources/engine.py"
        spec=importlib.util.spec_from_file_location("poisson_frozen_engine",p)
        ENGINE=importlib.util.module_from_spec(spec);spec.loader.exec_module(ENGINE)
        ENGINE.RUN=RUN
    return ENGINE


def verify():
    m=read(RUN/"manifest.json")
    if sha(RUN/"manifest.json")!=read(RUN/"freeze.json")["manifest_sha256"]:
        raise ValueError("Manifest changed")
    for p,d in m["source_sha256"].items():
        if sha(RUN/p)!=d: raise ValueError("Frozen source changed: "+p)
    seal=RUN/'postprocessing/freeze.json'
    if seal.exists() and sha(RUN/'postprocessing/production.py')!=read(seal)['source_sha256']:
        raise ValueError('Postprocessing source changed')
    archive_seal=RUN/'postprocessing/archive_v1.freeze.json'
    if archive_seal.exists():
        for name,digest in read(archive_seal)['files'].items():
            if sha(RUN/'postprocessing'/name)!=digest:raise ValueError('Compact archive source changed: '+name)
    return m


def compact_module():
    p=RUN/'postprocessing/compact_statistics.py'
    spec=importlib.util.spec_from_file_location('poisson_compact_statistics',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def prepare(smoke=False):
    if RUN.exists(): raise ValueError("Fresh immutable run required")
    for name in ("inputs","frozen_sources","workers","logs","results","submissions","analysis"):
        (RUN/name).mkdir(parents=True)
    shutil.copytree(TIMING/"native",RUN/"native")
    shutil.copy2(Path(__file__).with_name("poisson_ring_modes_20261027.py"),RUN/"frozen_sources/engine.py")
    shutil.copy2(__file__,RUN/"frozen_sources/production.py")
    E=engine();basals=[]
    for b in read(TIMING/"manifest.json")["basals"]:
        shutil.copy2(TIMING/b["input"],RUN/b["input"]);basals.append(b)
    jobs=[];conditions=[]
    if smoke:
        for full in (False,True):
            for mi,model in enumerate(E.MODELS):
                spring=model=="spring_ring"
                j=E.configuration(model,3. if spring else 60.,2,2)
                j.update(burn=2. if spring else 20.,sample_duration=(10. if spring else 100.)*(2 if full else 1),
                         seed=202610271000+mi,wall_seconds=1100.,slot=len(jobs),condition_index=mi,
                         smoke_kind="uninterrupted" if full else "resumed")
                jobs.append(j)
    else:
        for model in E.MODELS:
            spring=model=="spring_ring"
            for size in ((1.,2.,3.) if spring else (20.,40.,60.)):
                for ir in range(3):
                    for ia in range(3):
                        ci=len(conditions);slots=[]
                        for rep in range(4):
                            j=E.configuration(model,size,ir,ia,rep)
                            j.update(burn=1600. if spring else 16000.,sample_duration=12800. if spring else 128000.,
                                     slot=len(jobs),condition_index=ci,seed=202610271000000+len(jobs),wall_seconds=23.5*3600)
                            slots.append(j["slot"]);jobs.append(j)
                        conditions.append(dict(index=ci,model=model,size=size,rate_level=E.LEVELS[ir],
                                               amplitude_level=E.LEVELS[ia],slots=slots))
        assert len(jobs)==432 and len(conditions)==108 and len({j['seed'] for j in jobs})==432
    sources={str(p.relative_to(RUN)):sha(p) for part in ("native","frozen_sources") for p in sorted((RUN/part).rglob("*.py"))}
    initial_bytes=0
    for j in jobs:
        frames=int(round(j['sample_duration']/j['frame_interval']))+1
        initial_bytes+=frames*(16*(j['maximum_mode']+1)+8+(8 if j['model']=='spring_ring' else 0))
    m=dict(created_utc=now(),kind="preflight" if smoke else "production",jobs=jobs,conditions=conditions,
           basals=basals,source_sha256=sources,maximum_extension_level=2,repeats=4,
           extension_rule="All four seeds continue together until pooled mode-1 half-lag <= record/8, or level2 cap; first-frame decay does not trigger lengthening",
           maximum_sample_duration=dict(chemical=512000.,spring=51200.),
           estimator="Uncentered, actual-pair-count lag means; pool raw numerator and power over four seeds before normalization; ln(2)/linearly interpolated first half-height lag",
           bootstrap="Exact whole-seed nonparametric bootstrap: all35 count vectors for4 draws from4 seeds, with multinomial probabilities; no frame bootstrap",
           forcing=read(TIMING/'manifest.json')['forcing'],
           spring_intensity=read(TIMING/'manifest.json')['spring_intensity'],
           retained="Separate complex64 u/v Fourier coefficients and frame times; Spring circumference; one bounded full coupled restart state per segment, no profile histories",
           initial_logical_coefficient_bytes=initial_bytes,
           initial_logical_coefficient_gib=initial_bytes/2**30,
           maximum_logical_coefficient_gib=4*initial_bytes/2**30,
           resources=dict(partition="red",qos="small",cpus=1,memory="2G",walltime="24:00:00",array_throttle=None),
           authorization="User approved production launch after reviewing timing estimates, 27 October 2026")
    write(RUN/"manifest.json",m);write(RUN/"freeze.json",dict(manifest_sha256=sha(RUN/"manifest.json")))
    write(RUN/"expected_levels.json",{str(j['slot']):0 for j in jobs})
    entry=RUN/"frozen_sources/production.py"
    for name,command in (("worker.sh","worker"),("analysis.sh","analyze-condition")):
        script=RUN/name
        cmd=shlex.join([PY,"-B",str(entry),command,"--run",str(RUN),"--slot"])
        script.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+cmd+' "$SLURM_ARRAY_TASK_ID"'+(' --level "${POISSON_LEVEL:-0}"' if command=='worker' else '')+'\n')
        script.chmod(0o750)
    if not smoke:
        for part in ("data","metadata","vectors","viewers"):(OUT/"supporting"/part).mkdir(parents=True,exist_ok=True)
        text=(TIMING/'README.md').read_text().split('## Completed timing measurements')[0]
        text=text[text.index('The stochastic input is'):]
        (RUN/'README.md').write_text('# Poisson ring Fourier production campaign\n\n'
            '432 independent histories: 9 rate/amplitude conditions × 4 models × 3 sizes × 4 repeats.\n\n'
            'Chemical rings discard 16,000 time units, then initially record 128,000 at interval 0.5; Spring discards 1,600, then records 12,800 at interval 0.05. Modes are 0 through floor(C/ell), with ell=1.5 for chemical rings and .15 for Spring. Native local parameters and exact ring preparations are unchanged.\n\n'
            'Arrival densities are .0001/.0003/.001 (chemical) and .01/.03/.1 (Spring); exponential amplitude means are .3/1/3 times the established model references: .1, 6.4, 1 and 4. The two axes are varied independently.\n\n'
            'All four repeats of a condition continue together when the pooled mode-1 half-height lag exceeds one eighth of the record, up to 512,000 chemical or 51,200 Spring sampled units. Extension restarts retain the exact random-generator states, pending hazard threshold, live bursts and one final coupled state. Complete prefixes are preserved. Missing execution is an error; noncrossing at the cap remains an upper rate bound.\n\n'
            'Only u/v Fourier histories, frame time and Spring circumference are retained, plus small restart states. Segments concatenate in chronological order with no duplicate boundary frame. The valid extent of a running segment is recorded in progress.json; preallocated rows after that count are not data.\n\n'
            'The time average at each lag uses its actual pair count. Raw numerators and powers are pooled across four equal-duration seeds before normalization. The point estimate uses every frame lag through half the record. Exact four-seed bootstrap intervals enumerate all 35 count compositions with multinomial probabilities; censor statuses are preserved. No per-seed rate averaging or coherence exclusion is applied.\n\n'
            'For Spring plots, wavelength uses the sampled mean physical circumference across all four seeds; the material mode index is unchanged as circumference varies. Biological migration, conditional OU estimates and simulated activity-mode memory remain separate readouts.\n\n'+text)
        status("PREPARED",completed=0)
    print(json.dumps(dict(run=str(RUN),jobs=len(jobs),initial_gib=initial_bytes/2**30)))


def sbatch(name, *, array=None, command=None, script=None, dependency=None, level=0, hours="24:00:00",memory="2G"):
    cmd=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem='+memory,
         '--time='+hours,'--job-name='+name,'--output='+str(RUN/'logs'/('%A_%a.out' if array else '%j.out')),
         '--error='+str(RUN/'logs'/('%A_%a.err' if array else '%j.err'))]
    if array:cmd+=['--array='+array,'--export=ALL,POISSON_LEVEL='+str(level)]
    if dependency:cmd+=['--dependency='+dependency]
    if script:cmd+=[str(script)]
    else:
        entry=RUN/'postprocessing/archive_v1.py'
        if not entry.exists():entry=RUN/'postprocessing/production.py'
        if not entry.exists():entry=RUN/'frozen_sources/production.py'
        cmd+=['--wrap',shlex.join([PY,'-B',str(entry),command,'--run',str(RUN)])]
    r=subprocess.run(cmd,capture_output=True,text=True,check=True)
    jid=r.stdout.strip().split(';')[0]
    write(RUN/'submissions'/f'{name}_{jid}.json',dict(created_utc=now(),job_id=jid,command=cmd))
    return jid


def submit(smoke=False):
    m=verify()
    if (RUN/'submission.json').exists():raise ValueError('Already submitted')
    if not smoke:
        check=read(SMOKE/'results/validation.json')
        if not check.get('pass'):raise ValueError('Restart preflight has not passed')
    jid=sbatch('poisson_preflight' if smoke else 'poisson_modes',array='0-'+str(len(m['jobs'])-1),script=RUN/'worker.sh',
               hours='00:20:00' if smoke else '24:00:00')
    receipt=dict(created_utc=now(),job_id=jid,jobs=len(m['jobs']),array_throttle=None)
    if not smoke:
        receipt['controller_job_id']=sbatch('poisson_control',command='control',dependency='afterany:'+jid,hours='04:00:00',memory='4G')
        status('RUNNING',completed=0,submission=receipt)
    write(RUN/'submission.json',receipt)
    storage();print(json.dumps(receipt))


def worker(slot,level):
    m=verify();E=engine();E.worker(slot,level)
    if m['kind']=='preflight' and slot<4:E.worker(slot,1)


def segment(slot,level):return RUN/'workers'/f'{slot:03d}'/'segments'/f'{level:02d}'


def series(slot,mode,species='u'):
    level=int(read(RUN/'expected_levels.json')[str(slot)])
    arrays=[np.load(segment(slot,k)/f'{species}_modes.npy',mmap_mode='r')[:,mode].astype(np.complex128) for k in range(level+1)]
    return np.concatenate(arrays)


def validate_smoke():
    m=verify();comparisons=[]
    for slot in range(4):
        a,b=segment(slot,1),segment(slot+4,0)
        assert read(a/'status.json')['complete'] and read(b/'status.json')['complete']
        with np.load(a/'checkpoint.npz') as z:ya=z['state'][:-1]
        with np.load(b/'checkpoint.npz') as z:yb=z['state'][:-1]
        err=float(np.max(abs(ya-yb)))
        scale=float(np.max(abs(yb)))
        assert err < 3e-5*max(1.,scale),(slot,err,scale)
        ta=np.r_[np.load(segment(slot,0)/'time.npy'),np.load(a/'time.npy')]
        tb=np.load(b/'time.npy');assert np.allclose(ta,tb,rtol=0,atol=1e-9)
        for species in ('u','v'):
            xa=np.concatenate([np.load(segment(slot,0)/f'{species}_modes.npy'),np.load(a/f'{species}_modes.npy')])
            xb=np.load(b/f'{species}_modes.npy')
            assert np.max(abs(xa-xb)) < 3e-5*max(1.,float(np.max(abs(xb))))
        ca,cb=read(a/'checkpoint.json'),read(b/'checkpoint.json')
        assert ca['rng']==cb['rng'] and ca['event_count']==cb['event_count']
        assert abs(ca['next_hazard_threshold']-cb['next_hazard_threshold'])<1e-12
        comparisons.append(dict(model=m['jobs'][slot]['model'],max_final_state_difference=err,scale=scale,
                                matching_rng=True,matching_frame_times=True,events=ca['event_count']))
    result=dict(generated_utc=now(),**{'pass':True},comparisons=comparisons)
    write(RUN/'results/validation.json',result);print(json.dumps(result))


def status(state,**fields):
    r=dict(generated_utc=now(),state=state,run=str(RUN.relative_to(ROOT)),**fields)
    write(RUN/'results/status.json',r)
    if RUN==PRODUCTION:
        write(OUT/'supporting/metadata/status.json',r)
        (OUT/'README.md').write_text('# Poisson ring Fourier decorrelation campaign\n\n'
            f"Generated {r['generated_utc']}. Producer state: **{state}**.\n\n"
            '[Frozen design and mathematical definitions](../../../../data/runs/poisson_ring_modes_20261027/README.md). '
            '432 histories, 108 conditions, four independent seeds per condition. Only Fourier histories and bounded restart states are retained.\n\n'
            'The initial sampled duration is 128,000 chemical or 12,800 Spring native time units after burn-in; insufficient mode-1 records continue together to the declared fourfold cap. Whole-seed bootstrap intervals use the four repeats. Final rate tables and figures are produced after execution and extension checks finish.\n')


def control():
    m=verify();levels=read(RUN/'expected_levels.json');missing=[]
    for j in m['jobs']:
        p=segment(j['slot'],int(levels[str(j['slot'])]))/'status.json'
        if not p.exists() or not read(p).get('complete'):missing.append(j['slot'])
    if missing:
        status('EXECUTION_INCOMPLETE',missing=missing);storage()
        raise RuntimeError('Missing/failed executions, not censoring: '+str(missing))
    E=engine();rows=[];extend=[]
    for c in m['conditions']:
        ll={int(levels[str(s)]) for s in c['slots']};assert len(ll)==1
        level=ll.pop();j=m['jobs'][c['slots'][0]];T=j['sample_duration']*2**level
        raw=np.array([E.acf_raw(series(s,1)) for s in c['slots']]);curve=raw.mean(0)/raw[:,0].mean()
        gamma,kind,bound=E.half_rate(curve,j['frame_interval'])
        want=(kind=='beyond_analysis_lag' or (kind=='resolved' and T<8*math.log(2)/gamma))
        if want and level<m['maximum_extension_level']:extend+=c['slots']
        rows.append(dict(condition=c['index'],model=c['model'],size=c['size'],rate_level=c['rate_level'],
                         amplitude_level=c['amplitude_level'],extension_level=level,record_duration=T,
                         gamma1=gamma,status=kind,bound=bound,adequate_record=not want,
                         at_cap=level==m['maximum_extension_level']))
    table(RUN/'results/mode1_horizon_screen.csv',rows)
    if extend:
        nexts={int(levels[str(s)])+1 for s in extend};assert len(nexts)==1
        nextlevel=nexts.pop()
        for s in extend:levels[str(s)]=nextlevel
        write(RUN/'expected_levels.json',levels)
        jid=sbatch('poisson_extend'+str(nextlevel),array=','.join(map(str,extend)),script=RUN/'worker.sh',level=nextlevel)
        controller=sbatch('poisson_control',command='control',dependency='afterany:'+jid,hours='04:00:00',memory='4G')
        status('EXTENDING',completed=432,extension_level=nextlevel,extended_histories=len(extend),array_job_id=jid,controller_job_id=controller)
    else:
        script=RUN/'postprocessing/archive_analysis.sh'
        if not script.exists():script=RUN/'analysis.sh'
        jid=sbatch('poisson_analysis',array='0-107',script=script,hours='02:00:00',memory='4G')
        final=sbatch('poisson_finalize',command='finalize',dependency='afterany:'+jid,hours='02:00:00',memory='4G')
        status('ANALYZING',completed=432,analysis_job_id=jid,finalizer_job_id=final)
    storage()


def bootstrap_counts():
    counts=np.array([c for c in itertools.product(range(5),repeat=4) if sum(c)==4],float)
    weights=np.array([math.factorial(4)/math.prod(math.factorial(int(x)) for x in c)/4**4 for c in counts])
    assert len(counts)==35 and np.isclose(weights.sum(),1)
    return counts,weights


def weighted_endpoint(values,kinds,weights,q):
    order=np.argsort(values,kind='stable');i=order[min(np.searchsorted(np.cumsum(weights[order]),q),len(order)-1)]
    return float(values[i]),str(kinds[i])


def bootstrap_rates(raw,cadence):
    counts,weights=bootstrap_counts();curves=counts@raw;curves/=curves[:,0,None]
    hit=curves[:,1:]<=.5;has=hit.any(1);first=hit.argmax(1)+1
    rates=np.full(35,math.log(2)/((curves.shape[1]-1)*cadence));kinds=np.full(35,'upper_bound',dtype=object)
    for i in np.flatnonzero(has):
        j=int(first[i])
        if j==1:rates[i]=math.log(2)/cadence;kinds[i]='lower_bound'
        else:
            lag=cadence*(j-1+(curves[i,j-1]-.5)/(curves[i,j-1]-curves[i,j]))
            rates[i]=math.log(2)/lag;kinds[i]='resolved'
    lo,los=weighted_endpoint(rates,kinds,weights,.025);hi,his=weighted_endpoint(rates,kinds,weights,.975)
    return lo,hi,los,his,float(weights[kinds=='resolved'].sum())


def analyze_condition(ci):
    m=verify();c=m['conditions'][ci];E=engine();j=m['jobs'][c['slots'][0]]
    folder=RUN/'analysis'/f'{ci:03d}';folder.mkdir(exist_ok=True)
    level=int(read(RUN/'expected_levels.json')[str(c['slots'][0])]);T=j['sample_duration']*2**level
    C=c['size']
    if c['model']=='spring_ring':
        C=float(np.mean([np.concatenate([np.load(segment(s,k)/'circumference.npy',mmap_mode='r') for k in range(level+1)]).mean() for s in c['slots']]))
    rows=[];display=[];lag_index=None
    for n in range(j['maximum_mode']+1):
        raw=np.array([E.acf_raw(series(s,n)) for s in c['slots']]);power=raw[:,0].mean()
        if not np.isfinite(raw).all() or power<=0:raise ValueError('Invalid raw mode power')
        curve=raw.mean(0)/power
        gamma,kind,bound=E.half_rate(curve,j['frame_interval'])
        lo,hi,los,his,resolved_probability=bootstrap_rates(raw,j['frame_interval'])
        rows.append(dict(condition=ci,model=c['model'],size=c['size'],mean_physical_C=C,
                         rate_level=c['rate_level'],amplitude_level=c['amplitude_level'],rho=j['rho'],
                         amplitude_mean=j['amplitude_mean'],seeds=4,mode=n,wavelength=C/n if n else None,
                         gamma=gamma,status=kind,bound=bound,analysis_lag=(len(curve)-1)*j['frame_interval'],
                         crossing_half_lag=math.log(2)/gamma if gamma else None,
                         gamma_ci_low=lo,gamma_ci_high=hi,ci_low_status=los,ci_high_status=his,
                         bootstrap_resolved_probability=resolved_probability,record_duration=T,
                         extension_level=level,mean_mode_power=power))
        if lag_index is None:
            lag_index=np.unique(np.r_[np.arange(min(1025,len(curve))),np.geomspace(1,len(curve)-1,min(2048,len(curve))).astype(int)])
        display.append(curve[lag_index].astype('float32'))
    table(folder/'rates.csv',rows)
    np.savez_compressed(folder/'pooled_acf.npz',lag=lag_index*j['frame_interval'],mode=np.arange(len(rows)),R=np.array(display))
    write(folder/'status.json',dict(complete=True,generated_utc=now(),condition=ci,source_manifest_sha256=sha(RUN/'manifest.json'),
                                    extension_level=level,rows=len(rows),bootstrap_count_vectors=35))
    if (RUN/'postprocessing/compact_statistics.py').exists():
        compact_module().archive_condition(sys.modules[__name__],ci)


def storage():
    result=subprocess.run(['du','-sx','-B1',str(RUN)],capture_output=True,text=True,check=True)
    size=int(result.stdout.split()[0]);r=dict(generated_utc=now(),allocated_bytes=size,allocated_gib=size/2**30,
                                           path=str(RUN.relative_to(ROOT)),status=read(RUN/'results/status.json').get('state') if (RUN/'results/status.json').exists() else 'PREFLIGHT')
    write(RUN/'results/storage.json',r)
    if RUN!=PRODUCTION:return
    marker='<!-- poisson-ring-modes-20261027 begin -->';end='<!-- poisson-ring-modes-20261027 end -->'
    text=(f'{marker}\n\n## Poisson ring Fourier campaign, 27 October 2026\n\n'
          f"Measured **{r['generated_utc']}** with `du -sx -B1`: **{size:,} allocated bytes ({size/2**30:.6f} GiB)** at `polarity1d/manuscript/data/runs/poisson_ring_modes_20261027/`. Producer state: **{r['status']}**. This non-atomic measurement covers current allocation; it is not a final size forecast.\n\n"
          'Scientifically authoritative new ring assay: 432 independent histories, 108 conditions and four seeds per condition. Preserve per-seed complex active/inactive Fourier coefficients, exact inputs/frozen sources, bounded continuation checkpoints and scheduler/inference records. No full activity-profile histories are stored. Incomplete execution is distinct from finite-horizon censoring. Inventory does not authorize deletion.\n\n')
    m=read(RUN/'manifest.json')
    text+=f"The first-wave logical Fourier/time allocation is estimated at **{m['initial_logical_coefficient_gib']:.3f} GiB**, and the all-conditions extension cap at **{m['maximum_logical_coefficient_gib']:.3f} GiB** before small metadata/analysis additions. Actual allocated storage is measured above and updated by the campaign controller and finalizer. [Measurement receipt](polarity1d/manuscript/data/runs/poisson_ring_modes_20261027/results/storage.json).\n\n{end}\n"
    p=ROOT/'DATA_MANIFEST.md';old=p.read_text()
    if marker in old:
        start=old.index(marker);stop=old.index(end,start)+len(end);new=old[:start]+text.rstrip()+old[stop:]
    else:new=old.rstrip()+'\n\n'+text
    tmp=p.with_suffix('.md.poisson.tmp');tmp.write_text(new);tmp.replace(p)


def finalize():
    m=verify();rows=[];missing=[]
    if (RUN/'postprocessing/compact_statistics.py').exists():
        try:compact_module().check_archive(sys.modules[__name__])
        except Exception as exc:
            status('COMPACT_ARCHIVE_INCOMPLETE',error=repr(exc));storage();raise
    for c in m['conditions']:
        folder=RUN/'analysis'/f"{c['index']:03d}"
        if not (folder/'status.json').exists() or not read(folder/'status.json').get('complete'):
            missing.append(c['index']);continue
        with (folder/'rates.csv').open() as f:rows+=list(csv.DictReader(f))
    if missing:
        status('ANALYSIS_INCOMPLETE',missing_conditions=missing);storage();raise RuntimeError(str(missing))
    table(OUT/'supporting/data/mode_rates.csv',rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    E=engine();artifacts=[]
    for model in E.MODELS:
        fig,axes=plt.subplots(3,3,figsize=(10,9),sharex=True,sharey=True,constrained_layout=True)
        sizes=sorted({float(r['size']) for r in rows if r['model']==model});colors=dict(zip(sizes,('#76a9cf','#2878ad','#123b60')))
        for ir,rate in enumerate(E.LEVELS):
            for ia,amp in enumerate(E.LEVELS):
                ax=axes[ir,ia]
                for size in sizes:
                    selected=[r for r in rows if r['model']==model and r['rate_level']==rate and r['amplitude_level']==amp and float(r['size'])==size and int(r['mode'])>0]
                    for k,r in enumerate(selected):
                        x=float(r['wavelength']);resolved=r['status']=='resolved';y=float(r['gamma'] if resolved else r['bound'])
                        marker='o' if resolved else ('v' if r['status']=='beyond_analysis_lag' else '^')
                        ax.plot(x,y,marker,ms=3.5,color=colors[size],mfc=colors[size] if resolved else 'white',label=f'C={size:g}' if k==0 else None)
                        if r['ci_low_status']==r['ci_high_status']=='resolved':ax.vlines(x,float(r['gamma_ci_low']),float(r['gamma_ci_high']),color=colors[size],alpha=.4,lw=.65)
                ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.15)
                ax.set_title(f'{rate} rate; {amp} amplitude',fontsize=9)
                if ir==2:ax.set_xlabel('Wavelength λ (native length)')
                if ia==0:ax.set_ylabel('Decorrelation rate γ (native time⁻¹)')
        axes[0,0].legend(fontsize=7);fig.suptitle(E.NAMES[model]+' — Poisson bursts, four seeds per condition')
        png=OUT/f'poisson_ring_modes_{model}.png';pdf=OUT/'supporting/vectors'/f'poisson_ring_modes_{model}.pdf'
        fig.savefig(png,dpi=170);fig.savefig(pdf);plt.close(fig);artifacts+=[png,pdf]
    status('COMPLETE',completed=432,conditions=108,mode_rows=len(rows))
    p=OUT/'README.md';text=p.read_text()+'\n'
    for model in E.MODELS:text+=f"\n![{E.NAMES[model]}](poisson_ring_modes_{model}.png)\n"
    text+='\n[Mode-rate CSV](supporting/data/mode_rates.csv). Filled circles are resolved first half-height crossings; open downward/upward triangles show upper/lower rate bounds. Lines show pointwise 95% exact four-seed bootstrap intervals when both endpoints are resolved; endpoint censor statuses remain in the CSV. Mode zero is retained in the CSV and histories but is omitted from wavelength plots.\n'
    p.write_text(text)
    write(OUT/'supporting/metadata/artifacts.json',dict(generated_utc=now(),files={str(p.relative_to(OUT)):sha(p) for p in artifacts+[OUT/'supporting/data/mode_rates.csv']}))
    subprocess.run([PY,'-B',str(MS/'code/mechanisms/organize_length_dev.py'),'--group','poisson_ring_modes_20261027'],check=True,capture_output=True,text=True)
    storage()


def checks():
    counts,weights=bootstrap_counts()
    assert np.isclose(weights[(counts==1).all(1)][0],24/256)
    # Exact bootstrap must include the unresampled cohort and respect power pooling.
    raw=np.array([[1.,.7,.4],[2.,1.8,1.2],[3.,1.2,.3],[4.,3.2,1.6]])
    assert not np.allclose(raw.mean(0)/raw[:,0].mean(),(raw/raw[:,0,None]).mean(0))
    result=bootstrap_rates(raw,1.)
    assert all(np.isfinite(result[k]) for k in (0,1,4))
    print(json.dumps(dict(exact_bootstrap='pass',raw_power_pooling='pass',count_vectors=len(counts))))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=('prepare','prepare-smoke','submit','submit-smoke','worker','validate-smoke','control','analyze-condition','finalize','checks','storage'))
    ap.add_argument('--run',type=Path);ap.add_argument('--slot',type=int);ap.add_argument('--level',type=int,default=0);a=ap.parse_args()
    RUN=a.run.resolve() if a.run else SMOKE if a.command in ('prepare-smoke','submit-smoke','validate-smoke') else PRODUCTION
    if a.command in ('prepare','prepare-smoke'):prepare(a.command=='prepare-smoke')
    elif a.command in ('submit','submit-smoke'):submit(a.command=='submit-smoke')
    elif a.command=='worker':worker(a.slot,a.level)
    elif a.command=='analyze-condition':analyze_condition(a.slot)
    elif a.command=='validate-smoke':validate_smoke()
    else:globals()[a.command]()
