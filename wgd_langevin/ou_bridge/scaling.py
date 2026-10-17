"""Joint calibration with at most two length exponents; untouched test lengths."""
import argparse
import numpy as np
from scipy.optimize import least_squares
from design import *
from reduced import *
from fit import load_cases, prepare, residual, evaluate

LAWS = {'constant': [], 'tp_power': [0], 'two_times_power': [0,1], 'tp_gain_power': [0,5]}

def parameters(record, L):
    theta = np.array(record['theta'])
    for index, exponent in zip(record['indices'],record['exponents']):
        theta[index] *= (L/record['reference_length'])**exponent
    return theta

def run(index):
    model, family = [(m,f) for m in MODELS for f in FAMILIES][index]
    target=OUT/'scaling'/f'{model}_{family}.json'
    if target.exists(): return
    Ls=lengths(model,True);ref=Ls[1]
    fits=[json.loads((OUT/'fits'/f'{model}_{L:g}_{family}.json').read_text()) for L in Ls]
    prepared=[prepare(load_cases(model,L,True))[0] for L in Ls]
    scale=max(f['scale'] for f in fits)
    lo,hi=bounds(family)
    records=[]
    for law,indices in LAWS.items():
        if family=='gaussian' and law=='two_times_power': continue
        def decode_law(z):
            theta=decode(z[:len(lo)],family)
            exponents=z[len(lo):]
            # Keep physical parameters within the pilot bounds at every tested
            # length. Without this, a legal reference tp and exponent could
            # imply an arbitrarily stiff, out-of-design tp at the shortest L.
            for parameter,exponent in zip(indices,exponents):
                factors=(np.array(lengths(model))/ref)**exponent
                lower=LOW[parameter]/factors.min();upper=HIGH[parameter]/factors.max()
                slot=active(family).index(parameter)
                fraction=(z[slot]-lo[slot])/(hi[slot]-lo[slot])
                theta[parameter]=np.exp(np.log(lower)+fraction*np.log(upper/lower))
            return dict(theta=theta.tolist(),indices=indices,
                        exponents=exponents.tolist(),reference_length=ref)
        def objective(z):
            record=decode_law(z)
            pieces=[]
            for L,data in zip(Ls,prepared):
                theta=parameters(record,L)
                pieces.append(residual(encode(theta,family),family,data)/np.sqrt(len(Ls)))
            return np.concatenate(pieces)
        candidates=[]
        for start in [1,0,2]:
            z=np.r_[np.clip(encode(fits[start]['theta'],family),lo+1e-7,hi-1e-7),np.zeros(len(indices))]
            res=least_squares(objective,z,bounds=(np.r_[lo,np.full(len(indices),-3.)],
                                                 np.r_[hi,np.full(len(indices),3.)]),
                              max_nfev=240,ftol=3e-6,xtol=3e-6,gtol=3e-6,diff_step=1e-4)
            candidates.append(dict(**decode_law(res.x),cost=float(res.cost),success=bool(res.success),
                                   nfev=res.nfev,active_bounds=res.active_mask.tolist()))
        best=min(candidates,key=lambda r:r['cost'])
        record=dict(**best,law=law,starts=candidates)
        # Only calibration cost selects the law. No test-length or held-out response tuning.
        record['validation']=[]
        for L in lengths(model):
            rows,_=evaluate(parameters(record,L),family,load_cases(model,L),scale)
            record['validation'].extend([dict(length=L,length_role='calibration' if L in Ls else 'held_out',**r) for r in rows])
        records.append(record)
        print(model,family,law,best['cost'],flush=True)
    selected=records[0]
    # Require >=10% calibration improvement for each extra exponent beyond current choice.
    for r in sorted(records[1:],key=lambda r:(len(r['indices']),r['cost'])):
        if r['cost'] < .9*selected['cost']: selected=r
    save_json(target,dict(model=model,family=family,scale=scale,laws=records,
                         selected_law=selected['law'],selection='>=10% improvement in calibration cost; no held-out selection',
                         source_hashes=digest([HERE/'scaling.py',HERE/'reduced.py',HERE/'design.py'])))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int);run(p.parse_args().index)
