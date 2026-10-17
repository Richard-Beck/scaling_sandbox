"""Independent fits at nine lengths; all seven protocols used, no length constraint."""
import sys,json,time,argparse
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
BASE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(BASE))
import design,fit,spatial
from reduced import *
OUT=Path(__file__).resolve().parent/'results';OLD=design.OUT
LENGTHS={m:([21.,23.,25.,27.,29.,31.,33.,36.,39.]if m=='debelly'else[5.,7.5,10.,12.5,15.,20.,25.,32.,39.])for m in design.MODELS}
JOBS=[(m,L)for m in design.MODELS for L in LENGTHS[m]]
FITS=[(m,L,f)for m,L in JOBS for f in design.FAMILIES]

def generate(index):
    m,L=JOBS[index];old=OLD/'spatial'/f'{m}_{L:g}';target=OUT/'spatial'/f'{m}_{L:g}'
    target.parent.mkdir(parents=True,exist_ok=True)
    if old.exists():
        if not target.exists():target.symlink_to(old,target_is_directory=True)
        return
    spatial.OUT=OUT;spatial.jobs=lambda:JOBS;spatial.run(index)

def run(index):
    m,L,f=FITS[index];path=OUT/'fits'/f'{m}_{L:g}_{f}.json'
    if path.exists():return
    fit.OUT=OUT;cases=fit.load_cases(m,L);prepared,scale=fit.prepare(cases)
    # prepare() assigns 80% to trajectories, 20% to weak contrast. Rebalance
    # to eight equal observables: seven trajectories and one weak contrast.
    prepared=[(t,y,s,w*np.sqrt((7/8)/.8 if s.ndim==1 else (1/8)/.2))for t,y,s,w in prepared]
    lo,hi=bounds(f);rng=np.random.default_rng(20260927+index)
    default=np.array([.5,5.,.8 if f=='memory'else 3.,.15,.15,.1,0.,.05]);seeds=[encode(default,f)]
    nearest=min(design.lengths(m,True),key=lambda x:abs(x-L));old=OLD/'fits'/f'{m}_{nearest:g}_{f}.json'
    previous=json.loads(old.read_text());seeds+=[encode(s['theta'],f)for s in sorted(previous['starts'],key=lambda s:s['cost'])[:2]]
    while len(seeds)<8:
        t=default.copy();t[[0,1,3,4,5,7]]*=np.exp(rng.normal(0,1.4,6));t[2]=rng.uniform(.1,.98)if f=='memory'else rng.uniform(.3,8.);t[6]=rng.uniform(0,.3)if f=='self_activation'else 0.;seeds.append(encode(t,f))
    tick=time.monotonic();starts=[]
    for i,z in enumerate(seeds):
        res=least_squares(fit.residual,np.clip(z,lo+1e-7,hi-1e-7),args=(f,prepared),bounds=(lo,hi),max_nfev=240,ftol=2e-6,xtol=2e-6,gtol=2e-6,diff_step=1e-4)
        starts.append(dict(start=i,theta=decode(res.x,f).tolist(),cost=float(res.cost),success=bool(res.success),nfev=res.nfev,active_bounds=res.active_mask.tolist(),singular_values=np.linalg.svd(res.jac,compute_uv=False).tolist()))
        print(m,L,f,i,res.cost,res.success,flush=True)
    best=min(starts,key=lambda s:s['cost']);metrics,pred=fit.evaluate(best['theta'],f,cases,scale)
    design.save_json(path,dict(model=m,length=L,family=f,theta=best['theta'],cost=best['cost'],scale=scale,starts=starts,metrics=metrics,assays=design.ASSAYS,parameter_names=NAMES,generated=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-tick))


def polish(index):
    m,L,f=FITS[index];path=OUT/'fits'/f'{m}_{L:g}_{f}.json';record=json.loads(path.read_text())
    best=min(record['starts'],key=lambda s:s['cost'])
    if best['success']:return
    fit.OUT=OUT;cases=fit.load_cases(m,L);prepared,scale=fit.prepare(cases)
    prepared=[(t,y,s,w*np.sqrt((7/8)/.8 if s.ndim==1 else (1/8)/.2))for t,y,s,w in prepared]
    lo,hi=bounds(f)
    res=least_squares(fit.residual,np.clip(encode(best['theta'],f),lo+1e-9,hi-1e-9),args=(f,prepared),bounds=(lo,hi),max_nfev=960,ftol=2e-6,xtol=2e-6,gtol=2e-6,diff_step=1e-4)
    record['starts'].append(dict(start=len(record['starts']),theta=decode(res.x,f).tolist(),cost=float(res.cost),success=bool(res.success),nfev=res.nfev,active_bounds=res.active_mask.tolist(),singular_values=np.linalg.svd(res.jac,compute_uv=False).tolist(),continuation=True))
    best=min(record['starts'],key=lambda s:s['cost']);record.update(theta=best['theta'],cost=best['cost'],generated=datetime.now(timezone.utc).isoformat())
    record['metrics']=fit.evaluate(best['theta'],f,cases,scale)[0]
    design.save_json(path,record);print('continued',m,L,f,res.cost,res.success,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['spatial','fit','polish']);p.add_argument('index',type=int);a=p.parse_args();{'spatial':generate,'fit':run,'polish':polish}[a.mode](a.index)
