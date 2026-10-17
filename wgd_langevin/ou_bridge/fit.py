"""Eight-start calibration, held-out protocols, honest optimizer diagnostics."""
import argparse
import time
from datetime import datetime, timezone
import numpy as np
from scipy.optimize import least_squares
from design import *
from reduced import *

def load_cases(model, L, training=False):
    cases = []
    for assay in CALIBRATION if training else ASSAYS:
        path = OUT/'spatial'/f'{model}_{L:g}'/f'{assay}.npz'
        with np.load(path) as data:
            c = {k: data[k].copy() for k in ('time', 'p', 'q', 'contrast', 'fronts')}
        protocol = Protocol(model, assay)
        if not np.isclose(c['time'][-1],protocol.end):
            raise ValueError(f'Stale spatial protocol: {path}')
        c.update(assay=assay, spec=spec_for(protocol), unit=protocol.unit)
        c['t'] = c['time']/protocol.unit
        cases.append(c)
    return cases

def prepare(cases):
    prepared = []
    # A fixed model/length scale from calibration only, not per held-out assay.
    scale = max(.02, max(np.max(np.abs(c['p'])) for c in cases))
    for c in cases:
        t = c['t']
        split = c['spec'][2]+5 if c['assay']=='pulse' else c['spec'][1]
        # Equal weighting of early and late phases prevents conditioning dominance.
        idx = np.unique(np.r_[0, np.linspace(0,len(t)-1,min(280,len(t))).astype(int),
                             np.flatnonzero((t >= split) & (t <= split+2))[::2], len(t)-1])
        weights = np.zeros(len(idx))
        for mask in (t[idx] <= split, t[idx] > split):
            if mask.any(): weights[mask] = 1/np.sqrt(2*len(cases)*mask.sum())/scale
        prepared.append((t[idx], c['p'][idx], c['spec'], weights))
    plus=next((c for c in cases if c['assay']=='weak_plus'),None)
    minus=next((c for c in cases if c['assay']=='weak_minus'),None)
    if plus is not None and minus is not None:
        idx=np.flatnonzero(plus['t']>=plus['spec'][1])[::2]
        idx=np.unique(np.r_[0,idx])
        delta=plus['p'][idx]-minus['p'][idx]
        delta_scale=max(1e-4,float(np.max(np.abs(delta))))
        # A fifth equally weighted observable isolates the weak response from
        # the much larger common conditioning baseline.
        prepared=[(t,y,s,w*np.sqrt(4/5))for t,y,s,w in prepared]
        prepared.append((plus['t'][idx],delta,np.stack([plus['spec'],minus['spec']]),
                         np.full(len(idx),1/(delta_scale*np.sqrt(5*len(idx))))))
    return prepared, scale

def residual(z, family, prepared):
    theta = decode(z, family)
    pieces=[]
    for t,y,spec,w in prepared:
        if spec.ndim==2:
            pred=(integrate(theta,t,spec[0],FAMILIES.index(family))[:,0]-
                  integrate(theta,t,spec[1],FAMILIES.index(family))[:,0])
        else: pred=integrate(theta,t,spec,FAMILIES.index(family))[:,0]
        pieces.append((pred-y)*w)
    return np.concatenate(pieces)

def evaluate(theta, family, cases, scale):
    rows, predictions = [], {}
    for c in cases:
        pred = integrate(theta, c['t'], c['spec'], FAMILIES.index(family))
        mask = c['t'] >= (c['spec'][2] if c['assay']=='pulse' else c['spec'][1])
        y, p = c['p'], pred[:,0]
        tail = c['t'] >= c['t'][-1]-1.
        final, actual = float(np.median(p[tail])), float(np.median(y[tail]))
        conditioned = c['t'] < c['spec'][1]
        pre = float(np.median(y[conditioned][-20:]))
        predicted_pre = float(np.median(p[conditioned][-20:]))
        def event(v, before):
            if abs(before)<.02: return 'not_conditioned', None
            good = (v*before < 0)&(np.abs(v)>=.02)
            for j in np.flatnonzero(mask & good):
                end = np.searchsorted(c['t'], c['t'][j]+1.)
                if end < len(v) and np.all(good[j:end+1]):
                    return 'reversed', float((c['t'][j]-c['spec'][1])*c['unit'])
            return 'right_censored', None
        status, latency = event(y,pre)
        pstatus, platency = event(p,predicted_pre)
        rows.append(dict(assay=c['assay'], rmse=float(np.sqrt(np.mean((p-y)**2))),
                         post_nrmse=float(np.sqrt(np.mean((p[mask]-y[mask])**2))/scale),
                         final_actual=actual, final_predicted=final,
                         retained_actual=abs(actual)>=.02, retained_predicted=abs(final)>=.02,
                         physical_range=bool(np.max(np.abs(p))<=.5),
                         reversal_actual=status if c['assay']=='reversal' else None,
                         reversal_predicted=pstatus if c['assay']=='reversal' else None,
                         latency_actual=latency if c['assay']=='reversal' else None,
                         latency_predicted=platency if c['assay']=='reversal' else None,
                         multifront_fraction=float(np.mean(c['fronts']>=2))))
        predictions[c['assay']] = pred
    return rows, predictions

def run(index):
    model,L,family = fit_jobs()[index]
    path = OUT/'fits'/f'{model}_{L:g}_{family}.json'
    if path.exists(): return
    tick = time.monotonic()
    cases = load_cases(model,L,True)
    prepared,scale = prepare(cases)
    rng = np.random.default_rng(20260925+index)
    lo, hi = bounds(family)
    candidates = []
    defaults = np.array([.5,5., .8 if family=='memory' else 3., .15,.15,.1,0.,.05])
    for start in range(8):
        theta = defaults.copy()
        if start:
            theta[[0,1,3,4,5,7]] *= np.exp(rng.normal(0,1.4,6))
            theta[2] = rng.uniform(.1,.98) if family=='memory' else rng.uniform(.3,8.)
            theta[6] = rng.uniform(0,.3) if family=='self_activation' else 0.
        z = np.clip(encode(theta,family),lo+1e-7,hi-1e-7)
        result = least_squares(residual,z,args=(family,prepared),bounds=(lo,hi),max_nfev=240,
                               ftol=2e-6,xtol=2e-6,gtol=2e-6,diff_step=1e-4)
        sv = np.linalg.svd(result.jac,compute_uv=False)
        candidates.append(dict(start=start,cost=float(result.cost),nfev=result.nfev,
                               success=bool(result.success),message=result.message,
                               theta=decode(result.x,family).tolist(),active_bounds=result.active_mask.tolist(),
                               jacobian_singular_values=sv.tolist()))
        print(model,L,family,start,result.cost,result.success,flush=True)
    best = min(candidates,key=lambda c:c['cost'])
    theta = np.array(best['theta'])
    allcases = load_cases(model,L)
    rows,pred = evaluate(theta,family,allcases,scale)
    refinement = max(float(np.max(np.abs(pred[c['assay']]-integrate(theta,c['t'],c['spec'],FAMILIES.index(family),2.))))
                     for c in allcases)
    save_json(path,dict(model=model,length=L,family=family,theta=theta.tolist(),parameter_names=NAMES,
               scale=scale,cost=best['cost'],starts=candidates,metrics=rows,
               max_step_halving_error=refinement,wall_seconds=time.monotonic()-tick,
               generated=datetime.now(timezone.utc).isoformat(),
               source_hashes=digest([HERE/'fit.py',HERE/'reduced.py',HERE/'design.py'])))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('index',type=int)
    run(parser.parse_args().index)
