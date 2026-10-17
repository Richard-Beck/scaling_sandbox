"""Normalized two-assay pipeline: gradient reversal and finite countercue."""
from __future__ import annotations
from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace
from datetime import datetime,timezone
from collections import Counter
import argparse,csv,gzip,hashlib,json,os,sys,subprocess,shlex,math,shutil
import numpy as np
CODE=Path(__file__).resolve().parent
ROOT=CODE.parents[1]
sys.path.insert(0,str(ROOT/'polarity1d/src'))
from models import wave_pinning,otsuji,goryachev,legi,holmes_model3 as ho,debelly as db,spring as sp
from polarity1d.contract import Grid1D
from polarity1d.simulation import simulate
from polarity1d.solvers import laplacian_no_flux
PROTOCOL=json.loads((CODE/'protocol.json').read_text())
MODELS=['wave_pinning','otsuji','goryachev','legi','debelly','holmes_model3','spring']
RD_MODELS=set(MODELS[:4]+['holmes_model3'])
def make_model(key):
    return {'wave_pinning':wave_pinning,'otsuji':otsuji,'goryachev':goryachev,'legi':legi}[key].make_model()
def _sim(model, state, length, cells, stimulus, end, interval=2.5, max_step=10.0):
    state = np.asarray(state, dtype=float)
    if model.nonnegative:
        minimum = float(np.min(state))
        if minimum < -1e-8:
            raise ValueError(f"restart state violates nonnegative invariant by {minimum:g}")
        state = np.maximum(state, 0.0)
    return simulate(model, initial_state=state, length=float(length), cells=int(cells),
                    protocol=stimulus, end_time=float(end), output_interval=float(interval),
                    max_step=float(max_step), backend="bdf")


def _active(tr):
    if "active" in tr.fields:
        return tr.fields["active"]
    return tr.observable


def _orientation(tr):
    u = np.maximum(_active(tr), 0.0)
    center = 2.0 * tr.x / (tr.x[0] + tr.x[-1]) - 1.0
    denom = u.sum(axis=1)
    return np.divide(u @ center, denom, out=np.full(len(denom), np.nan), where=denom > 1e-12)


rd=SimpleNamespace(make_model=make_model,Grid1D=Grid1D,_sim=_sim,_active=_active,_orientation=_orientation,RD_MODELS=RD_MODELS)
ds=db
CAP=float(PROTOCOL["shared"]["preparation_cap"])
TOL=PROTOCOL["shared"]["relative_tolerance"]
WINDOW=float(PROTOCOL["shared"]["confirmation_window"])
DB_NAMES=['rac', 'rho', 'mca_bound', 'mca_unbound', 'displacement', 'velocity', 'tension', 'raw_left_displacement', 'raw_right_displacement', 'left_boundary_displacement', 'right_boundary_displacement']

class Cue:
    """Identical linear spatial cues within a model; assays differ in timing."""
    def __init__(self,kind,model,stop=None):
        self.kind={'train':'left','counter':'right'}.get(kind,kind)
        self.model=model;self.stop=stop
        self.event_times=(() if stop is None else (float(stop),)) + ((20.0,) if self.kind=='crossfade' else ())
    def __call__(self,t,x,L):
        q=np.asarray(x)/L
        amp=PROTOCOL['normalization']['amplitude_by_model'][self.model]
        if self.stop is not None and t>=self.stop:return np.zeros_like(q)
        if self.kind=='zero':return np.zeros_like(q)
        if self.kind=='left':return amp*(1-q)
        if self.kind=='right':return amp*q
        if self.kind=='crossfade':
            w=np.clip(t/20,0,1);return amp*((1-w)*(1-q)+w*q)
        raise ValueError(self.kind)

def model_for(key):
    if key=='holmes_model3':return replace(ho.HolmesModel3Eq56(),rac_activation=ho.ACTIVATION)
    if key=='debelly':return ds.DeBelly(ds.DeBellyParameters(**db._parameters().__dict__))
    if key=='spring':return sp
    return rd.make_model(key)

def initial(key,L,n=201):
    m=model_for(key)
    if key=='debelly':return None
    if key=='spring':return m.initial(n,float(L),perturb=1e-3)
    if key=='holmes_model3':return ho._low_state(m,L,n)
    grid=rd.Grid1D(L,n)
    return m.initial_state(grid)

def run(key,L,state,cue,end,interval=5,n=201,max_step=10):
    m=model_for(key)
    if key=='debelly':
        tr=ds.simulate_debelly(m,protocol=cue,length=L,cells=n,end_time=end,output_interval=interval,step=.5,initial_checkpoint=state)
        fields=tr.fields;active=fields['rac'];x=tr.x;t=tr.time
        final={name:fields[name][-1].copy() for name in DB_NAMES}
        ori=db._orientation(active)
    elif key=='spring':
        tr=sp._run(L,state,cue,end,sample=interval,cells=n)
        active=tr['active'];x=tr['x_fraction']*L;t=tr['time'];ori=tr['orientation']
        fields={name:tr[name] for name in ['active','inactive','length']};final=sp._end(tr)
    else:
        sim=ho._sim if key=='holmes_model3' else rd._sim
        tr=sim(m,state,L,n,cue,end,interval=interval,max_step=max_step)
        fields={name:tr.fields[name] for name in m.state_names};active=rd._active(tr)
        x=tr.x;t=tr.time;ori=ho._orientation(tr) if key=='holmes_model3' else rd._orientation(tr)
        final=np.vstack([fields[name][-1] for name in m.state_names])
    active=np.maximum(active,0)
    split=max(2,n//3);base=np.min(active,axis=1)
    return SimpleNamespace(time=np.asarray(t),fields=fields,active=active,x=np.asarray(x),state=final,
        orientation=np.asarray(ori),contrast=np.ptp(active,axis=1),
        left=np.max(active[:,:split],axis=1)-base,right=np.max(active[:,-split:],axis=1)-base)

def concatenate(a,b,offset):
    if a is None:
        b.time=b.time+offset;return b
    a.time=np.concatenate([a.time,offset+b.time[1:]])
    for name in ['active','orientation','contrast','left','right']:
        setattr(a,name,np.concatenate([getattr(a,name),getattr(b,name)[1:]],axis=0))
    a.fields={name:np.concatenate([a.fields[name],b.fields[name][1:]],axis=0) for name in a.fields}
    a.state=b.state;return a

def rates(key,L,tr,cue):
    if key=='debelly':return {}
    m=model_for(key);state=tr.state
    if key=='spring':
        rate=m.rhs(0,state,L0=L,d1=80,b=4,cue=0)
        physical=state[-1]-state[-2];q=np.linspace(0,1,len(tr.x))
        profile=cue(0,q*physical,physical)
        u=state[:-2:2];v=state[1:-2:2];delta=3.0+80.0*(physical-L)
        transfer=m.auth['f'](u,v,4.0+profile,5,6,2,delta)-m.auth['f'](u,v,4.0,5,6,2,delta)
        rate[:-2:2]+=transfer;rate[1:-2:2]-=transfer
        return dict(active=rate[:-2:2],inactive=rate[1:-2:2],length=rate[-1]-rate[-2])
    grid=rd.Grid1D(L,len(tr.x))
    rate=np.asarray(m.diffusivities)[:,None]*laplacian_no_flux(state,grid.dx)+m.reaction(0,grid.x,state,cue(0,grid.x,L))
    return {name:rate[i] for i,name in enumerate(m.state_names)}

def assessment(key,L,tr,cue,window):
    keep=tr.time>=tr.time[-1]-window
    t=tr.time[keep];fields={name:np.asarray(y)[keep] for name,y in tr.fields.items()}
    override={}
    if key=='debelly':
        fields={name:fields[name] for name in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','boundary_separation']}
        fields['displacement']=fields['displacement']-np.mean(fields['displacement'],axis=1)[:,None]
        fields['physical_length']=L+fields.pop('boundary_separation')
        override={'displacement':L,'physical_length':L}
    if key=='spring':override={'length':L}
    rr=rates(key,L,tr,cue);metrics={}
    for name,y in fields.items():
        scale=override.get(name,max(float(np.max(abs(y))),1e-8))
        change=float(np.max(np.ptp(y,axis=0)))/scale
        derivative=rr.get(name,(y[-1]-y[-2])/(t[-1]-t[-2]))
        drift=float(np.max(abs(derivative)))*window/scale
        metrics[name]=dict(relative_range=change,projected_relative_drift=drift,steady=bool(change<=TOL and drift<=TOL))
    return dict(steady=all(v['steady'] for v in metrics.values()),metrics=metrics,window=window,
                maximum_relative_range=max(v['relative_range'] for v in metrics.values()),
                maximum_projected_drift=max(v['projected_relative_drift'] for v in metrics.values()))

def settle(key,L,state,kind,window=WINDOW,n=201):
    cue=Cue(kind,key);trace=None;elapsed=0.;checks=[]
    while elapsed<CAP:
        duration=window if trace is None else min(100.,CAP-elapsed)
        chunk=run(key,L,state,cue,duration,interval=5,n=n)
        trace=concatenate(trace,chunk,elapsed);elapsed+=duration;state=chunk.state
        result=assessment(key,L,trace,cue,window)
        checks.append(dict(time=elapsed,**result));print(f'{key} L={L:g} {kind} t={elapsed:g} steady={result["steady"]}',flush=True)
        if result['steady']:break
        # Retain enough recent fields for confirmation; histories are saved as scalar checks.
        keep=trace.time>=elapsed-window
        trace.time=trace.time[keep]
        for name in ['active','orientation','contrast','left','right']:setattr(trace,name,getattr(trace,name)[keep])
        trace.fields={name:y[keep] for name,y in trace.fields.items()}
    return state,trace,dict(cue=kind,settled=result['steady'],time_to_settle=elapsed if result['steady'] else None,
                           elapsed=elapsed,cap=CAP,confirmation_window=window,checks=checks)

def left_eligible(tr,key,finite=False):
    keep=tr.time>=tr.time[-1]-WINDOW
    o=tr.orientation[keep];c=tr.contrast[keep]
    good=bool(np.max(o)<-.05 and np.ptp(o)<.005 and np.min(c)>=.05)
    extra={}
    if key=='debelly' and finite:
        keep=tr.time>=tr.time[-1]-1000
        extra={name:float(np.ptp(db._x50(tr.fields[name][keep]))) for name in ['rac','mca_bound']}
        good=good and all(v<=.002 for v in extra.values())
    return good,dict(orientation_max=float(np.max(o)),orientation_range=float(np.ptp(o)),minimum_contrast=float(np.min(c)),x50_ranges=extra)

def prepare(task):
    key,L,assay=task['model'],task['length'],task['assay'];n=task['cells']
    if assay not in PROTOCOL['active_assays']:raise ValueError('Assay is not active: '+assay)
    state=initial(key,L,n);phases=[]
    tr=None
    if not phases or phases[-1]['settled']:
        if assay=='gradient_reversal':kind='left'
        elif key=='spring':kind='train'
        else:kind='train' if key in ['holmes_model3','debelly'] else 'left'
        if kind is not None:
            state,tr,phase=settle(key,L,state,kind,n=n);phases.append(phase)
        if assay=='finite_countercue' and phases[-1]['settled']:
            state,tr,phase=settle(key,L,state,'zero',window=1000 if key=='debelly' else WINDOW,n=n);phases.append(phase)
    settled=all(p['settled'] for p in phases)
    criteria={}
    if settled:
        eligible,criteria=left_eligible(tr,key,finite=assay=='finite_countercue')
    else:eligible=False
    result=dict(**task,settled=settled,eligible=eligible,phases=phases,eligibility_diagnostics=criteria,
                status='qualified' if eligible else 'ineligible' if settled else 'unsettled',
                total_preparation_time=sum(p['elapsed'] for p in phases),cap_per_phase=CAP,relative_tolerance=TOL)
    return result,state,tr

def save_checkpoint(path,key,state,tr):
    contents=dict(time=tr.time,x=tr.x,active=tr.active,orientation=tr.orientation,contrast=tr.contrast,left=tr.left,right=tr.right)
    contents.update({'field_'+name:y for name,y in tr.fields.items()})
    if key=='debelly':contents.update({'state_'+name:value for name,value in state.items()})
    else:contents['state']=state
    np.savez_compressed(path,**contents)

def read_checkpoint(path,key):
    with np.load(path) as z:
        state={name:z['state_'+name].copy() for name in DB_NAMES} if key=='debelly' else z['state'].copy()
        tr=SimpleNamespace(**{name:z[name].copy() for name in ['time','x','active','orientation','contrast','left','right']},
                           fields={name[6:]:z[name].copy() for name in z.files if name.startswith('field_')},state=state)
    return state,tr

def sustained(time,good,dwell,start=0):
    for i in np.flatnonzero((time>=start)&good):
        j=np.searchsorted(time,time[i]+dwell,side='left')
        if j<len(time) and np.all(good[i:j+1]):return float(time[i]-start)
    return None


def save_trace(path,tr):
    fields={name:tr.fields[name] for name in ['length','boundary_separation','mca_bound','rho','tension'] if name in tr.fields}
    np.savez_compressed(path,time=tr.time,x=tr.x,active=tr.active,orientation=tr.orientation,
                        contrast=tr.contrast,left=tr.left,right=tr.right,**fields)


def gradient(task,state):
    key,L=task['model'],task['length'];n=task['cells']
    cue=Cue('crossfade',key)
    tr=run(key,L,state,cue,20000,interval=2.5,n=n)
    strict=key in ['spring','debelly']
    good=((tr.orientation>.05)&(tr.contrast>.05)) if strict else ((tr.orientation>=.05)&(tr.contrast>=.05))
    endpoint=sustained(tr.time,good,20)
    return dict(status='success' if endpoint is not None else 'censored',eligible=True,endpoint=endpoint,
                horizon=20000.,crossfade_duration=20.,endpoint_definition='From crossfade onset: right orientation and contrast 0.05 sustained 20 native units'),tr


def finite(task,state,durations):
    key,L=task['model'],task['length'];n=task['cells'];tests=[];pulse_profiles=[];final_profiles=[];tail_orientations=[]
    recovery_time=PROTOCOL["shared"]["recovery"]
    for D in durations:
        D=float(D)
        interval=min(2.5,D/20)
        pulse=run(key,L,state,Cue('counter',key,stop=D),D,interval=interval,n=n)
        recovery=run(key,L,pulse.state,Cue('zero',key),recovery_time,interval=2.5,n=n)
        keep=recovery.time>=recovery.time[-1]-500
        o=recovery.orientation[keep];c=recovery.contrast[keep]
        stable=bool(np.ptp(o)<.005 and np.min(c)>=.05)
        success=bool(stable and np.min(o)>.05)
        tests.append(dict(duration=D,success=success,settled=stable,final_orientation=float(o[-1]),
                          orientation_range_last_500=float(np.ptp(o)),minimum_contrast_last_500=float(np.min(c))))
        pulse_profiles.append(pulse.active[-1]);final_profiles.append(recovery.active[-1]);tail_orientations.append(o)
    success=[t['duration'] for t in tests if t['success']]
    result=dict(status='success' if success else 'censored',eligible=True,endpoint=min(success) if success else None,
                tests=tests,recovery_horizon=recovery_time,
                endpoint_definition='Shortest TESTED duration producing sustained autonomous right polarity after fixed recovery')
    diagnostics=dict(durations=np.asarray(durations),pulse_final_active=np.asarray(pulse_profiles),
                     recovery_final_active=np.asarray(final_profiles),recovery_tail_orientation=np.asarray(tail_orientations))
    return result,diagnostics

def now():return datetime.now(timezone.utc).isoformat()


def convert(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)


def write(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n');temporary.replace(path)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify():
    frozen=json.loads((HERE/'frozen_provenance.json').read_text())
    if sources()!=frozen['source_hashes']:raise RuntimeError('Runtime inputs differ from frozen campaign')
    return frozen['provenance_id']


def run_preparation(slot):
    if 'steady_state' in protocol():
        import settling
        return settling.prepare_pair(sys.modules[__name__],slot)
    pid=verify();task=[t for t in tasks() if t['action']=='steady_state_rerun'][slot];idx=task['task_index']
    record,state,tr=prepare(task)
    checkpoint=HERE/'checkpoints'/f'task_{idx:04d}.npz'
    save_checkpoint(checkpoint,task['model'],state,tr)
    record.update(provenance_id=pid,checkpoint=str(checkpoint.relative_to(ROOT)),checkpoint_sha256=sha(checkpoint),
                  generated_at=now(),slurm_job=os.environ.get('SLURM_ARRAY_JOB_ID'),slurm_slot=os.environ.get('SLURM_ARRAY_TASK_ID'))
    verify();write(HERE/'preparations'/f'task_{idx:04d}.json',record)
    print(json.dumps({k:record[k] for k in ['task_index','model','assay','length','status','total_preparation_time']}),flush=True)


def run_response(slot):
    pid=verify();job=json.loads((HERE/'response_manifest.json').read_text())[slot]
    task=tasks()[job['task_index']];idx=task['task_index'];assay=task['assay'];key=task['model']
    prep=json.loads((HERE/'preparations'/f'task_{idx:04d}.json').read_text())
    if not prep['eligible']:raise RuntimeError('Response manifest contains unqualified preparation')
    checkpoint=ROOT/prep['checkpoint']
    if sha(checkpoint)!=job['checkpoint_sha256']:raise RuntimeError('Prepared checkpoint changed')
    state,tr=read_checkpoint(checkpoint,key)
    if assay=='gradient_reversal':row,trace=gradient(task,state)
    elif assay=='finite_countercue':
        row,diagnostics=finite(task,state,job['durations'])
        diag=HERE/'trials'/f'chunk_{slot:05d}.npz';np.savez_compressed(diag,**diagnostics)
        row.update(durations=job['durations'],chunk_slot=slot,diagnostic_trace=str(diag.relative_to(ROOT)))
    else:raise ValueError(assay)
    if assay=='gradient_reversal':
        diag=HERE/'results'/f'trace_{idx:04d}.npz';save_trace(diag,trace)
        row['diagnostic_trace']=str(diag.relative_to(ROOT))
    row.update(**task,provenance_id=pid,generated_at=now(),preparation_file=str((HERE/'preparations'/f'task_{idx:04d}.json').relative_to(ROOT)),
               checkpoint_sha256=job['checkpoint_sha256'])
    out=HERE/'trials'/f'chunk_{slot:05d}.json' if assay=='finite_countercue' else HERE/'data'/f'task_{idx:04d}.json'
    verify();write(out,row)
    print(json.dumps(dict(slot=slot,task_index=idx,model=key,assay=assay,status=row['status'])),flush=True)


def response_plan():
    pid=verify();manifest=tasks();grids=protocol()['duration_grids'];jobs=[];counts={}
    for task in manifest:
        if task['action']=='reuse_reference':continue
        idx=task['task_index']
        p=HERE/'preparations'/f'task_{idx:04d}.json'
        if not p.exists():
            write(HERE/'data'/f'task_{idx:04d}.json',dict(**task,status='preparation_execution_error',eligible=False,endpoint=None,reason='Preparation did not finish; distinct from scientific settling timeout'))
            continue
        prep=json.loads(p.read_text());counts[prep['status']]=counts.get(prep['status'],0)+1
        if not prep['eligible']:
            write(HERE/'data'/f'task_{idx:04d}.json',dict(**task,status=prep['status'],eligible=False,endpoint=None,preparation=prep))
            continue
        base=dict(task_index=idx,checkpoint_sha256=prep['checkpoint_sha256'])
        if task['assay']=='finite_countercue':
            durations=grids[task['model']]['refined_durations']
            for i in range(0,len(durations),10):jobs.append(dict(**base,durations=durations[i:i+10]))
        else:jobs.append(base)
    for i,j in enumerate(jobs):j['slot']=i
    write(HERE/'response_manifest.json',jobs)
    write(HERE/'response_manifest_sha256.json',dict(sha256=sha(HERE/'response_manifest.json')))
    write(HERE/'preparation_summary.json',dict(generated_at=now(),counts=counts,response_jobs=len(jobs)))
    return jobs


# One immutable protocol and one scratch workspace per campaign.
def protocol():return PROTOCOL
def tasks():return protocol()['tasks']
def sources():
    paths=[CODE/'runner.py',CODE/'plot.py',CODE/'protocol.json',CODE/'environment.json']+list((CODE/'models').glob('*.py'))
    paths+=list((ROOT/'polarity1d/src/polarity1d').glob('*.py'))
    paths.append(CODE/'baseline/records.json.gz')
    if 'reuse_reference' in protocol():paths.append(CODE/protocol()['reuse_reference'])
    if 'steady_state' in protocol():paths.extend([CODE/'settling.py', CODE/protocol()['steady_state']['fixed_scales']])
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}

def setup():
    if (HERE/'frozen_provenance.json').exists():raise RuntimeError('Run already initialized; use a new --work directory')
    for name in ['data','preparations','checkpoints','trials','results','diagnostics','logs','histories']:(HERE/name).mkdir(parents=True,exist_ok=True)
    hashes=sources();pid=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    write(HERE/'frozen_provenance.json',dict(created_at=now(),provenance_id=pid,source_hashes=hashes,protocol=protocol()['shared']))

def submit(action, *, array=None, dependency=None):
    python=ROOT/'polarity1d/four_assay/environment.json'
    executable=json.loads(python.read_text())['python']
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([executable,str(CODE/'runner.py'),action,'--work',str(HERE),'--results',str(args.results.resolve())])
    if array:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    sbatch_args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=02:00:00',f'--job-name=polarity-{action}',f'--chdir={ROOT}',f'--output={HERE}/logs/{action}_%A_%a.out',f'--error={HERE}/logs/{action}_%A_%a.err']
    if array:sbatch_args.append('--array='+array)
    if dependency:sbatch_args.append('--dependency=afterok:'+dependency)
    return subprocess.check_output(sbatch_args+['--wrap',command],text=True).strip().split(';')[0]

def launch():
    if (args.results.resolve()/'records.json.gz').exists():raise RuntimeError('Results already exist; preserve them or choose another --results')
    setup();running_status('preparations_running');n=sum(t['action']=='steady_state_rerun' for t in tasks())
    if 'steady_state' in protocol():n=sum(t['assay']=='gradient_reversal' for t in tasks())
    preparation=submit('prepare',array=f'0-{n-1}')
    controller=submit('stage2',dependency=preparation)
    write(HERE/'launch.json',dict(stage='preparations_running',preparation_job_id=preparation,controller_job_id=controller,launched_at=now()))
    print(json.dumps(json.loads((HERE/'launch.json').read_text())))

def stage2():
    jobs=response_plan()
    running_status('responses_running')
    if any(not (HERE/'preparations'/f'task_{t["task_index"]:04d}.json').exists() for t in tasks() if t['action']=='steady_state_rerun'):raise RuntimeError('Missing preparation records')
    response=submit('respond',array=f'0-{len(jobs)-1}')
    final=submit('finalize',dependency=response)
    launch=json.loads((HERE/'launch.json').read_text());launch.update(stage='responses_running',response_job_id=response,finalize_job_id=final,response_jobs=len(jobs));write(HERE/'launch.json',launch)
    print(json.dumps(launch))

def load_reference():
    with gzip.open(CODE/protocol()['reuse_reference'],'rt') as f:return json.load(f)


def collect():
    cases=[];preps=[];jobs=json.loads((HERE/'response_manifest.json').read_text())
    reference=load_reference()
    reference_cases={r['task_index']:r for r in reference['cases']}
    reference_preps={r['task_index']:r for r in reference['preparations']}
    for task in tasks():
        if task['action']=='reuse_reference':
            cases.append(reference_cases[task['task_index']]);preps.append(reference_preps[task['task_index']]);continue
        idx=task['task_index'];p=HERE/'preparations'/f'task_{idx:04d}.json'
        prep=json.loads(p.read_text()) if p.exists() else None
        if prep:preps.append(prep)
        if task['assay']=='finite_countercue' and prep and prep['eligible']:
            planned=[j for j in jobs if j['task_index']==idx]
            chunks=[json.loads((HERE/'trials'/f'chunk_{j["slot"]:05d}.json').read_text()) for j in planned]
            tests=sorted([t for c in chunks for t in c['tests']],key=lambda t:t['duration'])
            expected=protocol()['duration_grids'][task['model']]['refined_durations']
            if [t['duration'] for t in tests]!=expected:raise RuntimeError(f'Incomplete duration coverage: task {idx}')
            passing=[t['duration'] for t in tests if t['success']]
            row=dict(**task,status='success' if passing else 'censored',eligible=True,endpoint=min(passing) if passing else None,tests=tests,tested_count=len(tests),planned_trial_count=len(expected),complete=True,passing_durations=passing,nonmonotonic_success_pattern=any(tests[i]['success'] and not tests[i+1]['success'] for i in range(len(tests)-1)))
        else:row=json.loads((HERE/'data'/f'task_{idx:04d}.json').read_text())
        cases.append(row)
    if len(cases)!=len(tasks()) or len(preps)!=len(tasks()):raise RuntimeError('Incomplete campaign')
    return dict(cases=cases,preparations=preps,generated_at=now(),provenance=json.loads((HERE/'frozen_provenance.json').read_text()),
                reuse=dict(reference=protocol().get('reuse_reference'),sha256=sha(CODE/protocol()['reuse_reference']),reused_cases=sum(t['action']=='reuse_reference' for t in tasks()),fresh_cases=sum(t['action']=='steady_state_rerun' for t in tasks())))

def csvwrite(path,rows):
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def compare_chained(out):
    """Audit the preparation dependency graph and compare intentionally new inputs."""
    import settling
    with gzip.open(out/'records.json.gz','rt') as f:fresh=json.load(f)
    previous=load_reference()
    def key(t):return t['model'],t['assay'],t['length']
    before={key(t):t for t in previous['cases']}
    preps={key(t):t for t in fresh['preparations']}
    expected={key(t) for t in tasks()}
    errors=[];rows=[]
    if {key(t) for t in fresh['cases']}!=expected or set(preps)!=expected:errors.append('Incomplete case/preparation coverage')
    audited=set();maximum_history_error=0.
    for t in tasks():
        k=key(t);p=preps[k];g=preps[(t['model'],'gradient_reversal',t['length'])]
        if p['phases'][0]!=g['phases'][0]:errors.append(f'Conditioning differs between assays: {k}')
        if p['phases'][0]['reference_checkpoint_sha256']!=g['checkpoint_sha256']:errors.append(f'Wrong conditioning reference: {k}')
        if t['assay']=='finite_countercue' and len(p['phases'])==2:
            if p['phases'][1]['input_checkpoint_sha256']!=g['checkpoint_sha256']:errors.append(f'Broken conditioning -> release dependency: {k}')
            if p['phases'][1]['reference_checkpoint_sha256']!=p['checkpoint_sha256']:errors.append(f'Wrong release reference: {k}')
        checkpoint=ROOT/p['checkpoint']
        if sha(checkpoint)!=p['checkpoint_sha256']:errors.append(f'Checkpoint hash changed: {k}')
        state,tr=read_checkpoint(checkpoint,t['model'])
        eligible,_=left_eligible(tr,t['model'],finite=t['assay']=='finite_countercue') if p['settled'] else (False,{})
        if eligible!=p['eligible']:errors.append(f'Eligibility mismatch: {k}')
        for phase in p['phases']:
            if phase['history'] in audited:continue
            audited.add(phase['history'])
            if not phase['settled']:errors.append(f'Unverified profile: {k} {phase["cue"]}');continue
            with np.load(ROOT/phase['history']) as z:
                times=z['time'];distances=np.zeros(len(times))
                if not np.all(np.diff(times)==protocol()['steady_state']['sample_interval']):errors.append(f'History sampling mismatch: {k}')
                for name,scale in phase['fixed_scales'].items():
                    y=z['field_'+name]
                    if not np.isfinite(y).all():errors.append(f'Nonfinite history: {k} {name}')
                    distances=np.maximum(distances,abs(y-y[-1]).reshape(len(times),-1).max(axis=1)/scale)
                maximum_history_error=max(maximum_history_error,float(abs(distances-z['maximum_scaled_distance']).max()))
                outside=np.flatnonzero(distances>phase['arrival_tolerance'])
                arrival=float(times[outside[-1]+1]) if len(outside) and outside[-1]+1<len(times) else 0. if not len(outside) else None
                if arrival!=phase['time_to_settle'] or arrival is None or times[-1]-arrival<phase['confirmation_window']:errors.append(f'Arrival mismatch: {k}')
                if not phase['checks'][-1].get('verified'):errors.append(f'Holdout not accepted: {k}')
    for row in fresh['cases']:
        old=before[key(row)]
        rows.append(dict(model=row['model'],assay=row['assay'],length=row['length'],
                         previous_status=old['status'],status=row['status'],previous_endpoint=old.get('endpoint'),endpoint=row.get('endpoint'),
                         lost_readout=old['status']=='success' and row['status']!='success',
                         gained_readout=old['status']!='success' and row['status']=='success'))
    losses=[p for p in rows if p['lost_readout']];gains=[p for p in rows if p['gained_readout']]
    csvwrite(out/'comparison_to_previous.csv',rows)
    result=dict(generated_at=now(),passed=not errors and not losses,execution_valid=not errors,
                comparison_kind='Revision 11 versus revision 10: all preparation states intentionally refined; no exact response controls',
                source_revision=11,coverage_complete={key(t) for t in fresh['cases']}==expected,
                cases=len(fresh['cases']),preparations=len(preps),distinct_settling_events=len(audited),
                finite_duration_tests=sum(len(t.get('tests',[])) for t in fresh['cases']),
                lost_readouts=losses,gained_readouts=gains,audit_errors=errors,
                maximum_history_distance_error=maximum_history_error,
                previous_success=sum(t['status']=='success' for t in previous['cases']),
                success=sum(t['status']=='success' for t in fresh['cases']))
    write(out/'comparison_to_previous.json',result)
    # Retain an explicit v8 comparison for established output consumers.
    with gzip.open(CODE/'baseline/records.json.gz','rt') as f:base=json.load(f)
    original={key(t):t for t in base['cases']}
    baseline_rows=[dict(model=t['model'],assay=t['assay'],length=t['length'],
                        baseline_status=original[key(t)]['status'],status=t['status'],
                        baseline_endpoint=original[key(t)].get('endpoint'),endpoint=t.get('endpoint')) for t in fresh['cases']]
    csvwrite(out/'comparison_to_baseline.csv',baseline_rows)
    write(out/'comparison_to_baseline.json',dict(generated_at=now(),source_revision=11,
          comparison_kind='Descriptive v8 comparison; all settled response inputs intentionally changed',
          execution_valid=not errors,passed=not errors and all(t['status']=='success' for t in fresh['cases'] if original[key(t)]['status']=='success'),
          lost_readouts=[p for p in baseline_rows if p['baseline_status']=='success' and p['status']!='success']))
    if errors:raise RuntimeError('Chained preparation audit failed: '+str(errors[:5]))
    print(json.dumps(result));return result

def export_tables(out):
    with gzip.open(out/'records.json.gz','rt') as f:records=json.load(f)
    endpoints=[];settling=[]
    for r in records['cases']:
        row={k:r.get(k) for k in ['task_index','baseline_task_index','model','assay','length','status','eligible','endpoint','tested_count','planned_trial_count']}
        row.update(original_status=row['status'],original_endpoint=row['endpoint'],mode=row['status'])
        endpoints.append(row)
    for p in records['preparations']:
        for i,phase in enumerate(p['phases']):
            settling.append(dict(task_index=p['task_index'],model=p['model'],assay=p['assay'],length=p['length'],phase_index=i,cue=phase['cue'],settled=phase['settled'],time_to_settle=phase['time_to_settle'],elapsed=phase['elapsed'],cap=phase['cap'],confirmation_window=phase['confirmation_window'],preparation_status=p['status'],eligible=p['eligible'],time_definition=phase.get('time_definition'),reference_checkpoint=phase.get('reference_checkpoint'),input_checkpoint_sha256=phase.get('input_checkpoint_sha256')))
    csvwrite(out/'endpoints.csv',endpoints);csvwrite(out/'settling_times.csv',settling)
    def read(name):return list(csv.DictReader((out/name).open()))
    return read('endpoints.csv'),read('settling_times.csv')


def compare(out):
    """Validate paired coverage and exact controls, allowing intended cue effects."""
    if 'steady_state' in protocol():return compare_chained(out)
    with gzip.open(CODE/'baseline/records.json.gz','rt') as f:base=json.load(f)
    with gzip.open(out/'records.json.gz','rt') as f:fresh=json.load(f)
    def key(row):return row['model'],row['assay'],row['length']
    expected={key(t) for t in tasks()}
    actual={key(t) for t in fresh['cases']}
    coverage_complete=actual==expected and len(fresh['cases'])==len(expected) and len(fresh['preparations'])==len(expected)
    before={key(t):t for t in base['cases']};base_preps={key(t):t for t in base['preparations']}
    reference=load_reference();previous_cases={key(t):t for t in reference['cases']}
    previous_preps={key(t):t for t in reference['preparations']}
    conditioning_differences=[];previous_losses=[];previous_control_tests=0
    after_preps={key(t):t for t in fresh['preparations']}
    mismatches=[];comparisons=[];finite_tests=0;control_tests=0
    def same(a,b,path):
        if isinstance(a,dict):
            if not isinstance(b,dict):mismatches.append(dict(path=path,reason='missing dictionary'));return
            for k,v in a.items():same(v,b.get(k),path+'/'+k)
        elif isinstance(a,list):
            if not isinstance(b,list) or len(a)!=len(b):mismatches.append(dict(path=path,reason='list length differs'));return
            for i,(x,y) in enumerate(zip(a,b)):same(x,y,path+'/'+str(i))
        elif isinstance(a,(int,float)) and not isinstance(a,bool):
            if b is None or not np.isclose(a,b,rtol=1e-9,atol=1e-10):mismatches.append(dict(path=path,baseline=a,fresh=b))
        elif a!=b:mismatches.append(dict(path=path,baseline=a,fresh=b))
    for row in fresh['cases']:
        k=key(row);old=before[k];prep=after_preps[k];bp=base_preps[k]
        unchanged=(row['assay']=='gradient_reversal' and row['model']!='spring') or row['model'] not in ['holmes_model3','debelly','spring']
        previous=previous_cases[k];previous_prep=previous_preps[k]
        if not (row['model']=='spring' and row['assay']=='gradient_reversal'):
            for field in ['status','eligible','endpoint']:same(previous.get(field),row.get(field),'normalized_v9/'+str(k)+'/'+field)
            if 'tests' in previous:same(previous['tests'],row.get('tests',[]),'normalized_v9/'+str(k)+'/tests');previous_control_tests+=len(row.get('tests',[]))
            same(previous_prep['phases'],prep['phases'],'normalized_v9/'+str(k)+'/phases')
            same(previous_prep['eligibility_diagnostics'],prep['eligibility_diagnostics'],'normalized_v9/'+str(k)+'/eligibility')
        if previous['status']=='success' and row['status']!='success':previous_losses.append(dict(model=row['model'],assay=row['assay'],length=row['length'],status=row['status']))
        tests=row.get('tests',[]);finite_tests+=len(tests)
        if row['assay']=='finite_countercue' and prep['eligible']:
            if [t['duration'] for t in tests]!=protocol()['duration_grids'][row['model']]['refined_durations']:
                coverage_complete=False
        if unchanged:
            for field in ['status','eligible','endpoint']:same(old.get(field),row.get(field),str(k)+'/'+field)
            if 'tests' in old:same(old['tests'],tests,str(k)+'/tests');control_tests+=len(tests)
            same(bp['phases'],prep['phases'],str(k)+'/phases')
            same(bp['eligibility_diagnostics'],prep['eligibility_diagnostics'],str(k)+'/eligibility')
        lost=old['status']=='success' and row['status']!='success'
        gained=old['status']!='success' and row['status']=='success'
        comparisons.append(dict(model=row['model'],assay=row['assay'],length=row['length'],
            task_index=row['task_index'],baseline_task_index=old['task_index'],unchanged_control=unchanged,
            baseline_status=old['status'],normalized_status=row['status'],lost_readout=lost,gained_readout=gained,
            baseline_endpoint=old.get('endpoint'),normalized_endpoint=row.get('endpoint'),
            baseline_preparation_time=bp['total_preparation_time'],normalized_preparation_time=prep['total_preparation_time']))
    for model in MODELS:
        lengths=sorted({t['length'] for t in tasks() if t['model']==model})
        for L in lengths:
            a=after_preps[(model,'gradient_reversal',L)]['phases'][0]
            b=after_preps[(model,'finite_countercue',L)]['phases'][0]
            if a['settled']!=b['settled'] or a['time_to_settle']!=b['time_to_settle']:
                conditioning_differences.append(dict(model=model,length=L,gradient=a['time_to_settle'],countercue=b['time_to_settle']))
    losses=[r for r in comparisons if r['lost_readout']];gains=[r for r in comparisons if r['gained_readout']]
    csvwrite(out/'comparison_to_baseline.csv',comparisons)
    groups=[]
    for model in MODELS:
        for assay in protocol()['active_assays']:
            rows=[r for r in comparisons if r['model']==model and r['assay']==assay]
            groups.append(dict(model=model,assay=assay,cases=len(rows),
                baseline_success=sum(r['baseline_status']=='success' for r in rows),
                normalized_success=sum(r['normalized_status']=='success' for r in rows),
                losses=sum(r['lost_readout'] for r in rows),gains=sum(r['gained_readout'] for r in rows)))
    result=dict(generated_at=now(),passed=coverage_complete and not mismatches and not losses and not previous_losses and not conditioning_differences,
        comparison_kind='Readout preservation versus v8/v9; exact unchanged controls; shared conditioning verified; Spring gradient preparation changed',
        coverage_complete=coverage_complete,cases=len(fresh['cases']),preparations=len(fresh['preparations']),
        baseline_success=sum(r['baseline_status']=='success' for r in comparisons),
        normalized_success=sum(r['normalized_status']=='success' for r in comparisons),
        finite_duration_tests=finite_tests,unchanged_control_pulse_tests=control_tests,
        control_mismatch_count=len(mismatches),control_mismatches=mismatches[:100],
        lost_readouts=losses,gained_readouts=gains,groups=groups,source_revision=protocol().get('revision',9),
        normalized_v9_lost_readouts=previous_losses,normalized_v9_control_pulse_tests=previous_control_tests,
        shared_conditioning_differences=conditioning_differences,fresh_cases=fresh.get('reuse',{}).get('fresh_cases'),reused_cases=fresh.get('reuse',{}).get('reused_cases'))
    write(out/'comparison_to_baseline.json',result);print(json.dumps(result))
    return result


def running_status(stage):
    out=args.results.resolve();out.mkdir(exist_ok=True)
    stamp=now()
    write(out/'status.json',dict(generated_at=stamp,complete=False,stage=stage,source_revision=protocol().get('revision',9),
                              total_planned_cases=len(tasks()),active_assays=protocol()['active_assays']))
    (out/'index.html').write_text(f'<html><body><h1>Normalized gradient reversal and finite countercue</h1><p>Generated {stamp}; producer running: {stage}. Final figures pending.</p></body></html>')


def finalize():
    verify();records=collect();out=args.results.resolve();out.mkdir(exist_ok=True)
    if (out/'records.json.gz').exists():raise RuntimeError('Results already exist; choose --output or preserve them first')
    with gzip.open(out/'records.json.gz','wt') as f:json.dump(records,f,default=convert,separators=(',',':'),allow_nan=False)
    if 'steady_state' in protocol():
        destination=CODE/'reference/chained_v11_20261001'
        destination.mkdir(exist_ok=False)
        shutil.copytree(HERE/'checkpoints',destination/'profiles')
        shutil.copytree(HERE/'preparations',destination/'preparations')
        shutil.copy2(HERE/'frozen_provenance.json',destination/'provenance.json')
        write(destination/'manifest.json',dict(generated_at=now(),source_revision=11,rule=protocol()['steady_state'],
              chain={'conditioning_input':'native resting state','release_input':'settled conditioning','gradient_input':'settled conditioning','countercue_input':'settled release'},
              profiles={p.name:sha(p) for p in (destination/'profiles').glob('*.npz')},
              histories=str((HERE/'histories').relative_to(ROOT)),records=str((out/'records.json.gz').relative_to(ROOT))))
    write(out/'run_provenance.json',records['provenance'])
    subprocess.run([sys.executable,str(CODE/'plot.py'),'--results',str(out)],check=True)
    comparison=compare(out)
    endpoints=list(csv.DictReader((out/'endpoints.csv').open()))
    write(out/'status.json',dict(generated_at=now(),complete=True,status_counts=dict(Counter(r['status'] for r in endpoints)),comparison_passed=comparison['passed'],source_revision=protocol().get('revision',9),total_endpoints=len(records['cases']),total_preparations=len(records['preparations']),active_assays=protocol()['active_assays']))
    state=json.loads((HERE/'launch.json').read_text());state.update(stage='finished',finished_at=now(),comparison_passed=comparison['passed'],results=str(out.relative_to(ROOT)));write(HERE/'launch.json',state)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['launch','prepare','respond','stage2','finalize','compare','status']);ap.add_argument('--work',type=Path,default=CODE/'work');ap.add_argument('--slot',type=int);ap.add_argument('--results',type=Path,default=CODE/'results')
    args=ap.parse_args();HERE=args.work.resolve()
    if args.action=='launch':launch()
    elif args.action=='prepare':run_preparation(args.slot)
    elif args.action=='respond':
        if sha(HERE/'response_manifest.json')!=json.loads((HERE/'response_manifest_sha256.json').read_text())['sha256']:raise RuntimeError('Response plan changed')
        run_response(args.slot)
    elif args.action=='stage2':stage2()
    elif args.action=='finalize':finalize()
    elif args.action=='compare':compare(args.results.resolve())
    else:print((HERE/'launch.json').read_text())
