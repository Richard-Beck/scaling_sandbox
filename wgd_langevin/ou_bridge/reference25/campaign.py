"""Fixed-length polarity representation experiment. All protocols enter every fit."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
from datetime import datetime, timezone

import numpy as np
from numba import njit
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import design
import reduced
from spatial import halfmass
from polarity1d.registry import make_model
from polarity1d.simulation import simulate
from polarity1d.metrics import spatial_contrast, profile_domain_count

OUT = HERE / 'results'
LENGTH = 25.
MODELS = design.MODELS
FAMILIES = design.FAMILIES


def suite(model):
    unit, condition, pulse, peak, width = design.clocks(model)
    result = []
    # Segment: start, stop, geometry (0 endpoint cosine, 1 linear gradient,
    # 2 oscillating gradient), left amplitude, right amplitude, period.
    def add(name, group, segments, end, focus=0.):
        result.append(dict(name=name, group=group, segments=segments, end=end,
                           focus=focus, unit=unit, width=width))
    def local(start, duration, left, right):
        return [start, start+duration, 0, peak*left, peak*right, 0]
    for label, strength, duration in [('weak', .5, 1), ('standard', 1, 1),
                                      ('strong', 2, 1), ('short', 1, .5), ('long', 1, 2)]:
        add('pulse_'+label, 'pulse', [local(0, pulse*duration, strength, 0)],
            pulse*duration+condition)
    add('pulse_right', 'pulse', [local(0, pulse, 0, 1)], pulse+condition)
    tonic = [0, condition, 1, .9*peak, .1*peak, 0]
    for label, delta in [('plus', .05), ('minus', -.05)]:
        add('weak_'+label, 'weak', [tonic,
            [condition, condition+unit, 1, peak*(.9-delta*.4), peak*(.1+delta*.4), 0],
            [condition+unit, condition+10*unit, 1, .9*peak, .1*peak, 0]],
            condition+10*unit, condition)
    for label, factor, strength in [('short', .5, 1), ('standard', 2, 1), ('long', 8, 1),
                                     ('strong', 2, 2)]:
        duration = pulse*factor
        add('counter_'+label, 'counter', [local(0, pulse, 1, 0),
            local(condition, duration, 0, strength)], condition+duration+condition,
            condition)
    add('counter_mirror', 'counter', [local(0, pulse, 0, 1),
        local(condition, 2*pulse, 1, 0)], 2*condition+2*pulse, condition)
    for factor in (1, 4, 16):
        period = factor*unit
        add(f'tracking_{factor}', 'tracking', [tonic,
            [condition, condition+4*period, 2, .9*peak, .1*peak, period]],
            condition+4*period, condition)
    for label, left, right in [('left', 1.01, .99), ('right', .99, 1.01)]:
        add('competition_'+label, 'competition', [local(0, pulse, left, right)],
            pulse+condition)
    return result


class Protocol:
    def __init__(self, case):
        self.case = case
        self.event_times = tuple(sorted({v for s in case['segments'] for v in s[:2]}))

    def __call__(self, t, x, length):
        width = self.case['width']
        for start, stop, kind, left, right, period in self.case['segments']:
            if start <= t < stop:
                if kind == 0:
                    z = x/length
                    a = np.where(z <= width, .5*(1+np.cos(np.pi*z/width)), 0)
                    b = np.where(1-z <= width, .5*(1+np.cos(np.pi*(1-z)/width)), 0)
                    return left*a+right*b
                if kind == 2:
                    mean, amplitude = (left+right)/2, (left-right)/2
                    left = mean+amplitude*np.cos(2*np.pi*(t-start)/period)
                    right = 2*mean-left
                return left+(right-left)*x/length
        return np.zeros_like(x)


@njit(cache=True)
def forcing(t, segments, width):
    for s in segments:
        start, stop, kind, left, right, period = s
        if start <= t < stop:
            if kind == 0:
                return width*(right-left)
            if kind == 2:
                return .5*(right-left)*np.cos(2*np.pi*(t-start)/period)
            return .5*(right-left)
    return 0.


@njit(cache=True)
def integrate(theta, times, segments, width, family, refinement=1.):
    out = np.zeros((len(times), 2))
    p = q = 0.
    t = 0.
    dt = min(.1, .2*theta[0], .2*theta[1], .2/(theta[6]+1e-6))/refinement
    for i in range(len(times)):
        while t < times[i]-1e-10:
            h = min(dt, times[i]-t)
            # Resolve the sharp self-activation transition near p=0. A fixed
            # step based only on relaxation times can shift competition onset.
            if family == 3 and theta[6] > 0:
                h = min(h, .1*(theta[7]+abs(p))/(theta[6]*refinement))
            for s in segments:
                for event in s[:2]:
                    if event > t+1e-10:
                        h = min(h, event-t)
            a,b = reduced.drift(p,q,theta,family,forcing(t+1e-11,segments,width))
            c,d = reduced.drift(p+h*a/2,q+h*b/2,theta,family,forcing(t+h/2,segments,width))
            e,f = reduced.drift(p+h*c/2,q+h*d/2,theta,family,forcing(t+h/2,segments,width))
            g,j = reduced.drift(p+h*e,q+h*f,theta,family,forcing(t+h-1e-11,segments,width))
            p += h*(a+2*c+2*e+g)/6
            q += h*(b+2*d+2*f+j)/6
            t += h
            if not np.isfinite(p+q) or abs(p)+abs(q)>1e5:
                out[i:] = 1e5
                return out
        out[i,0], out[i,1] = p,q
    return out


def now():
    return datetime.now(timezone.utc).isoformat()


def normalized_segments(case):
    s = np.array(case['segments'], dtype=float)
    s[:, [0,1,5]] /= case['unit']
    return s


def sample_times(case):
    edges = sorted({0., case['end'], *Protocol(case).event_times})
    parts = []
    for a,b in zip(edges[:-1], edges[1:]):
        if b>a:
            parts.append(np.linspace(a,b,121))
            parts.append(a+np.geomspace(min(.05, (b-a)/100), b-a, 100))
    if case['group']=='tracking':
        parts.append(np.linspace(case['focus'],case['end'],481))
    return np.unique(np.concatenate(parts))


SPATIAL_JOBS = [(m,c) for m in MODELS for c in suite(m)]
FIT_JOBS = [(m,f) for m in MODELS for f in FAMILIES]


def spatial_job(index):
    model, case = SPATIAL_JOBS[index]
    folder = OUT/'spatial'/model
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/(case['name']+'.npz')
    if path.exists():
        return
    options = dict(length=LENGTH, cells=126, times=sample_times(case))
    if model=='debelly':
        options['step'] = .5
    else:
        options.update(backend='bdf', max_step=2.5, rtol=2e-6)
    tick = time.monotonic()
    trajectory = simulate(make_model(model, design.DEBELLY if model=='debelly' else {}),
                          protocol=Protocol(case), **options)
    profiles = trajectory.observable
    if not np.all(np.isfinite(profiles)) or np.min(profiles)<-1e-6:
        raise ValueError('Invalid spatial concentrations')
    p = 2*halfmass(profiles)
    np.savez_compressed(path, time=trajectory.time, p=p, x=trajectory.x,
        profiles=profiles.astype(np.float32),
        contrast=np.array([spatial_contrast(row) for row in profiles]),
        fronts=np.array([profile_domain_count(row) for row in profiles]),
        zero_mass=np.abs(profiles).sum(axis=1)<1e-12,
        cue=np.array([forcing(t/case['unit'], normalized_segments(case),case['width'])
                      for t in trajectory.time]))
    design.save_json(path.with_suffix('.json'), dict(model=model, length=LENGTH,
        case=case, metadata=trajectory.metadata, generated=now(), seconds=time.monotonic()-tick))
    print(model,case['name'],time.monotonic()-tick,flush=True)


def load_cases(model):
    cases=[]
    for c in suite(model):
        c=c.copy()
        with np.load(OUT/'spatial'/model/(c['name']+'.npz')) as data:
            c.update({k:data[k] for k in ('time','p','cue')})
        c['t']=c['time']/c['unit']
        c['spec']=normalized_segments(c)
        cases.append(c)
    return cases


def prepare(cases):
    prepared=[]
    # Equal weight per stimulus family, then per trajectory, then per segment.
    # Each trajectory retains its absolute p scale. No fitted readout gain.
    for c in cases:
        edges=sorted({0.,c['end'],*Protocol(c).event_times})
        w=np.zeros(len(c['time']))
        group_count=sum(o['group']==c['group'] for o in cases)
        for a,b in zip(edges[:-1],edges[1:]):
            mask=(c['time']>=a)&(c['time']<b)
            if b==edges[-1]:
                mask|=c['time']==b
            w[mask]=1/np.sqrt(5*group_count*(len(edges)-1)*mask.sum())
        prepared.append((c,w))
    return prepared


def residual(z,family,prepared):
    theta=reduced.decode(z,family)
    return np.concatenate([(integrate(theta,c['t'],c['spec'],c['width'],
                            FAMILIES.index(family))[:,0]-c['p'])*w for c,w in prepared])


def fit_job(index):
    model,family=FIT_JOBS[index]
    path=OUT/'fits'/f'{model}_{family}.json'
    if path.exists():
        return
    cases=load_cases(model)
    prepared=prepare(cases)
    # Change the old half-range coordinate's amplitude bounds by exactly 2.
    lo,hi=reduced.bounds(family)
    for j,k in enumerate(reduced.active(family)):
        if k in (3,4,5,6,7):
            if reduced.LOG[k]:
                lo[j]+=np.log(2); hi[j]+=np.log(2)
            else:
                lo[j]*=2; hi[j]*=2
    rng=np.random.default_rng(260926+index)
    default=np.array([.5,5.,.8 if family=='memory' else 3.,.3,.3,.2,0.,.1])
    starts=[]
    tick=time.monotonic()
    for seed in range(4):
        theta=default.copy()
        if seed:
            theta[[0,1,3,4,5,7]]*=np.exp(rng.normal(0,1.2,6))
            theta[2]=rng.uniform(.1,.98) if family=='memory' else rng.uniform(.3,8)
            theta[6]=rng.uniform(0,.6) if family=='self_activation' else 0.
        result=least_squares(residual,np.clip(reduced.encode(theta,family),lo+1e-7,hi-1e-7),
            args=(family,prepared),bounds=(lo,hi),max_nfev=240,
            ftol=2e-6,xtol=2e-6,gtol=2e-6,diff_step=1e-4)
        starts.append(dict(seed=seed,theta=reduced.decode(result.x,family).tolist(),
            cost=float(result.cost),success=bool(result.success),message=result.message,
            nfev=result.nfev,active_bounds=result.active_mask.tolist()))
        print(model,family,seed,result.cost,result.success,flush=True)
        design.save_json(OUT/'progress'/f'{model}_{family}.json',dict(starts=starts,generated=now()))
    best=min(starts,key=lambda r:r['cost'])
    theta=np.array(best['theta'])
    metrics=[]
    pred={}
    for c,w in prepared:
        p=integrate(theta,c['t'],c['spec'],c['width'],FAMILIES.index(family))[:,0]
        pred[c['name']]=p
        metrics.append(dict(assay=c['name'],group=c['group'],
            rmse=float(np.sqrt(np.mean((p-c['p'])**2))),
            weighted_squared_error=float(np.sum(((p-c['p'])*w)**2)),
            observed_final=float(c['p'][-1]),predicted_final=float(p[-1]),
            max_abs_p=float(np.max(np.abs(p))),out_of_bounds=bool(np.any(np.abs(p)>1))))
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path.with_suffix('.npz'),**pred)
    design.save_json(path,dict(model=model,family=family,length=LENGTH,
        theta=theta.tolist(),parameter_names=reduced.NAMES,starts=starts,cost=best['cost'],
        weighted_rmse=float(np.sqrt(2*best['cost'])),metrics=metrics,
        generated=now(),seconds=time.monotonic()-tick))
    publish()


def publish():
    import subprocess
    subprocess.run([sys.executable,str(HERE/'report.py')],check=True)


def manifest():
    OUT.mkdir(exist_ok=True)
    paths=[Path(__file__),HERE.parent/'reduced.py',HERE.parent/'design.py',
           HERE.parent/'spatial.py']
    paths+=list((design.ROOT/'polarity1d/src/polarity1d').glob('*.py'))
    design.save_json(OUT/'manifest.json',dict(generated=now(),length=LENGTH,
        models=MODELS,families=FAMILIES,suites={m:suite(m) for m in MODELS},
        hashes={str(p.relative_to(design.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['manifest','spatial','fit'])
    parser.add_argument('index',type=int,nargs='?',default=0)
    args=parser.parse_args()
    if args.mode=='manifest': manifest()
    elif args.mode=='spatial': spatial_job(args.index)
    else: fit_job(args.index)
