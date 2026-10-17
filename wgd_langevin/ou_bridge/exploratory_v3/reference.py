"""Reference-only calibration followed by a fixed sqrt(2) length intervention."""
import sys,json,time,argparse
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent/'refit_v2'))
import campaign as c
OUT=HERE/'results';RATIO=np.sqrt(2.)
JOBS=[(f,s)for s in c.STUDIES for f in c.FAMILIES]

def curves(z,names,f,study,n=1024,seed=117):
    # A common regular sampling grid permits time-averaged MSD at many lags.
    cad=1/6;dur=24.;theta,eta,force=c.decode(z,names,0)
    lags=np.unique(np.r_[np.geomspace(1,144,36).round().astype(int),3,6,12,24,48,72,144])*cad
    if f in ('gaussian','memory'):
        rates,w=c.modes(theta,f);mean=theta[0]*force/(1-theta[2]if f=='memory'else 1.)
        return lags,2*c.displacement_variance(lags,rates,w)+(mean*lags)**2
    # Longest time scales can increase after transfer; burn-in scales with them.
    dt=min(.01,.1*theta[0],.1*theta[1]);burn=max(96.,16*max(theta[0],theta[1]))
    inc,_=c.ensemble(theta,eta,force,cad,dur,1.,n,seed,dt,burn)
    x=np.concatenate([np.zeros((n,1,2)),inc.cumsum(axis=1)],axis=1)
    vals=[]
    for lag in lags:
        k=int(round(lag/cad));vals.append(float(np.mean(np.sum((x[:,k:]-x[:,:-k])**2,axis=2))))
    return lags,np.array(vals)

def transfer(z,names,source):
    result=np.array(z).copy();multipliers={};inactive=[]
    for i,e in zip(source['indices'],source['exponents']):
        name={0:'tp',1:'tm',5:'force'}[i];factor=float(RATIO**e)
        if name in names:result[names.index(name)]+=np.log(factor);multipliers[name]=factor
        else:inactive.append(name)
    return result,multipliers,inactive

def run(index,continuation=False):
    f,study=JOBS[index];path=OUT/'reference'/f'{study}_{f}.json'
    old=json.loads(path.read_text())if path.exists()else None
    if old and not continuation:return
    names,lo,hi=c.specification(f,(),study);obs=[r for r in c.observations(study)if r['group']==0]
    linear=f in ('gaussian','memory');n=512 if linear else 96
    default=dict(tp=.6,tm=2.,r=.6 if f=='memory'else 2.,vs=.5,eta=.5,self_a=.1,self_K=.2,force=.1)
    base=np.log([default[k]for k in names]);rng=np.random.default_rng(20260928+index)
    seeds=[base]
    while len(seeds)<6:seeds.append(base+rng.normal(0,1.8,len(base)))
    records=[];tick=time.monotonic()
    if old:seeds=[];records=old['starts']
    for start,z in enumerate(seeds):
        fit=least_squares(lambda v:c.score(c.predictions(v,names,f,study,n=n),obs,study)[0],np.clip(z,lo+1e-6,hi-1e-6),bounds=(lo,hi),max_nfev=90,diff_step=.003,ftol=1e-5,xtol=1e-5,gtol=1e-5)
        records.append(dict(z=fit.x.tolist(),cost=float(fit.cost),success=bool(fit.success),start=start))
        print(study,f,start,fit.cost,flush=True)
    best=min(records,key=lambda r:r['cost']);polishn=2048 if linear else 384
    fit=least_squares(lambda v:c.score(c.predictions(v,names,f,study,n=polishn),obs,study)[0],np.log([old['parameters'][k]for k in names])if old else best['z'],bounds=(lo,hi),max_nfev=240 if old else 60,diff_step=.002,ftol=2e-5,xtol=2e-5,gtol=2e-5)
    z=fit.x;res,amp,rows,_=c.score(c.predictions(z,names,f,study,n=polishn),obs,study)
    check=c.score(c.predictions(z,names,f,study,n=4096 if linear else 1024,seed=937,dt=.01,burn=96.),obs,study,amp)
    # Alternative reference fits are descriptive multistart spread, not a CI.
    alternatives=[]
    for r in sorted(records,key=lambda r:r['cost']):
        if r['cost']>best['cost']+.5:continue
        rr,aa,_,_=c.score(c.predictions(r['z'],names,f,study,n=polishn),obs,study)
        if rr@rr>res@res+1.:continue
        alternatives.append((np.array(r['z']),aa))
        if len(alternatives)==3:break
    sources=[s for s in c.selections()if s['family']==f]
    lags,ref=curves(z,names,f,study);predictions=[]
    alternative_refs=[]
    for az,aa in alternatives:
        _,v=curves(az,names,f,study,n=512);alternative_refs.append((v*aa**2).tolist())
    for source in sources:
        wz,mult,inactive=transfer(z,names,source);_,wgd=curves(wz,names,f,study)
        alt=[]
        for az,aa in alternatives:
            w,_,_=transfer(az,names,source);_,v=curves(w,names,f,study,n=512);alt.append((v*aa**2).tolist())
        predictions.append(dict(model=source['model'],law=source['law'],exponents=source['exponents'],indices=source['indices'],multipliers=mult,inactive=inactive,spatial_calibration_improvement=source['improvement'],wgd_msd=(wgd*amp**2).tolist(),alternative_wgd_msd=alt))
    c.save_json(path,dict(study=study,family=f,length_ratio=RATIO,names=names,parameters=dict(zip(names,np.exp(z).tolist())),shared_amplitude=amp,reference_score=float(res@res),reference_rows=rows,reference_validation=check[2],starts=records,polish_success=bool(fit.success),lags=lags.tolist(),reference_msd=(ref*amp**2).tolist(),alternative_reference_msd=alternative_refs,predictions=predictions,generated=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-tick))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int);p.add_argument('--continue-fit',action='store_true');a=p.parse_args();run(a.index,a.continue_fit)
