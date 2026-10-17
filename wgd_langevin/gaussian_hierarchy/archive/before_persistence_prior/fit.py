"""Analytical optimization; frozen-parameter trajectory checks only afterward."""
import argparse,json,os,time
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import rice
from model import *

def analytic(z):
    pars,drifts,_=decode(z);out={}
    for i,s in enumerate(STUDIES):
        out[s]=[]
        cad,T=SCHEDULE[s]
        for g in (0,1):
            tau,sigma=pars[i,g]
            var=lambda t:float(displacement_variance(t,np.array([1/tau]),np.array([sigma*sigma*tau/2])))
            v=var(cad);vT=var(T)
            if s!='Hosny_baseline':
                speed=np.sqrt(np.pi*v/2)/cad
                p=dict(speed=speed,path_length=speed*T,endpoint=np.sqrt(np.pi*vT/2),msd=2*var(11.1111111111),early_msd=2*v)
                p['directness']=p['endpoint']/p['path_length']
            else:
                mu=drifts[g];step=np.sqrt(v);end=np.sqrt(vT)
                path=float(rice.mean(abs(mu)*cad/step,scale=step))*T/cad
                endpoint=rice.ppf([.25,.5,.75],abs(mu)*T/end,scale=end)
                parallel=(mu*T+ndtri([.25,.5,.75])*end)/path
                p=dict(path_length=path,speed=path/T,quantiles=dict(path_length=[path]*3,directness=(endpoint/path).tolist(),parallel=parallel.tolist(),endpoint=endpoint.tolist(),perpendicular=(ndtri([.25,.5,.75])*end/path).tolist()))
            out[s].append(p)
    return out

def residual(z,scenario):
    r,rows=data_residuals(analytic(z))
    # Hosny fits medians only; quartiles are held out for trajectory validation.
    mask=[row['study']!='Hosny_baseline' or row.get('quantile')==1 for row in rows]
    return np.r_[r[mask],pooling_residuals(z,SCENARIOS[scenario])]

def atomic(path,obj):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False));tmp.replace(path)

def run(index):
    scenario=list(SCENARIOS)[index//8];start=index%8;t=time.time()
    sol=least_squares(residual,initial(start),args=(scenario,),bounds=(LOW,HIGH),max_nfev=1200,ftol=1e-10,xtol=1e-10,gtol=1e-8)
    r,rows=data_residuals(analytic(sol.x));rows=[x for x in rows if x['study']!='Hosny_baseline' or x.get('quantile')==1]
    pool=pooling_residuals(sol.x,SCENARIOS[scenario])
    obj=dict(scenario=scenario,start=start,z=sol.x.tolist(),total_score=float(sol.fun@sol.fun),data_score=sum(x['residual']**2 for x in rows),pooling_penalty=float(pool@pool),rows=rows,hyperparameters=hyperparameters(sol.x,SCENARIOS[scenario]),success=bool(sol.success),message=sol.message,nfev=sol.nfev,optimality=float(sol.optimality),seconds=time.time()-t,generated=datetime.now(timezone.utc).isoformat(),provenance=provenance())
    atomic(OUT/f'fit_{index:02}.json',obj)
    print(json.dumps({k:obj[k] for k in ['scenario','start','total_score','success','seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int);a=p.parse_args();run(a.index)
