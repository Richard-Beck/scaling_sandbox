"""Bounded, conditional literature compatibility, with fixed polarity length laws.

Linear families use exact stationary Gaussian sampled-displacement covariance.
Nonlinear families use an isotropic p/q vector extension with Gaussian forcing
on p only. This is a declared closure, not noise inferred from spatial residuals.
The screen is descriptive; absence of a grid hit is not a global impossibility.
"""
import argparse
import sys
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.linalg import toeplitz
from numba import njit
from design import *
from scaling import parameters

OBS = pd.read_csv(ROOT/'wgd_langevin/data/observations.csv')
BANK_STEPS=2**23
BURN_STEPS=2**19
WINDOW_FRACTION=64
sys.path.insert(0,str(ROOT/'wgd_langevin/hosny'))
from data import DATA as HOSNY

def modes(theta,family):
    tp,tm,r=theta[:3]
    if family=='gaussian' or r<1e-8:
        return np.array([1/tp]),np.array([tp/2])
    A=np.array([[-1/tp,r/tp],[1/tm,-1/tm]])
    ev=np.sort(np.linalg.eigvals(A).real)
    slow,fast=-ev[1],-ev[0]
    ws=(fast-1/tp)/(fast-slow)
    # Fixed instantaneous p-noise coefficient eta=1 across lengths;
    # substitute closure B=eta^2*tp/2. Overall scale profiled in um/h.
    return np.array([slow,fast]),tp/2*np.array([ws,1-ws])

def displacement_variance(t,rates,weights):
    t=np.asarray(t)
    out=np.zeros_like(t,dtype=float)
    for rate,w in zip(rates,weights):
        z=rate*t
        shape=np.where(np.abs(z)<1e-3,z*z/2-z**3/6+z**4/24-z**5/120,z+np.expm1(-z))
        out+=2*w*shape/(rate*rate)
    return out

def summaries(increments):
    path=np.linalg.norm(increments,axis=2).sum(axis=1)
    net=increments.sum(axis=1);end=np.linalg.norm(net,axis=1)
    return dict(path_length=path,endpoint=end,directness=end/np.maximum(path,1e-20),
                parallel=net[:,1]/np.maximum(path,1e-20),perpendicular=net[:,0]/np.maximum(path,1e-20))

def linear_tracks(theta,family,k,cadence,duration,n=1024,seed=431):
    rates,w=modes(theta,family)
    steps=round(duration/cadence)
    lag=np.arange(steps+1)*cadence/k
    var=displacement_variance(lag,rates,w)*k*k
    first=np.r_[var[1],.5*(var[2:]+var[:-2]-2*var[1:-1])]
    cov=toeplitz(first)
    eigen,vectors=np.linalg.eigh(cov)
    if eigen.min() < -1e-7*max(var[1],1e-30): raise ValueError('invalid increment covariance')
    root=vectors*np.sqrt(np.maximum(eigen,0))[None,:]
    rng=np.random.default_rng(seed)
    return np.einsum('ij,njk->nik',root,rng.normal(size=(n,steps,2)))

@njit(cache=True)
def nonlinear_bank(theta,dt,noise,steps,seed,burn_steps=2**19):
    np.random.seed(seed)
    tp,tm,r,ps,vs,g,a,K=theta
    p=np.array([.05,0.]);q=np.array([.02,0.]);x=np.zeros(2)
    stride=16
    out=np.empty((steps//stride+1,2));pol=np.empty_like(out)
    out[0]=x;pol[0]=p
    for i in range(burn_steps+steps):
        pn=np.sqrt((p*p).sum());qn=np.sqrt((q*q).sum())
        dp=(-p+r*q/(1+qn/vs))/tp+a*p*pn/(K*K+pn*pn)
        dq=(p/(1+pn/ps)-q)/tm
        dW=noise*np.sqrt(dt)*np.random.normal(0.,1.,2)
        pp=p+dt*dp+dW;qq=q+dt*dq
        pnn=np.sqrt((pp*pp).sum());qnn=np.sqrt((qq*qq).sum())
        dp2=(-pp+r*qq/(1+qnn/vs))/tp+a*pp*pnn/(K*K+pnn*pnn)
        dq2=(pp/(1+pnn/ps)-qq)/tm
        newp=p+.5*dt*(dp+dp2)+dW
        q=q+.5*dt*(dq+dq2)
        if i>=burn_steps:
            x+=.5*dt*(p+newp)
            j=i-burn_steps+1
            if j%stride==0:
                out[j//stride]=x;pol[j//stride]=newp
        p=newp
    pol[0]=pol[1]
    return out,pol

def bank_tracks(bank,dt,k,cadence,duration,n=1024,seed=431):
    # Random overlapping windows: Monte Carlo sampling, not independent cells.
    step=cadence/k/dt;frames=round(duration/cadence)
    span=frames*step
    if step<1 or span>len(bank)/WINDOW_FRACTION: return None
    rng=np.random.default_rng(seed)
    start=rng.uniform(0,len(bank)-span-2,n)
    t=start[:,None]+np.arange(frames+1)[None,:]*step
    idx=np.floor(t).astype(int);fraction=(t-idx)[...,None]
    pos=(bank[idx]*(1-fraction)+bank[idx+1]*fraction)*k
    return np.diff(pos,axis=1)

def fit_amplitude(items,shared=True):
    # Each item: group, metric, prediction, observed, low, high.
    numerator={};denominator={}
    for group,metric,pred,obs,low,high in items:
        key='shared' if shared else group
        power=2 if metric=='msd' else 0 if metric=='directness' else 1
        width=np.log(high/low)/2
        numerator[key]=numerator.get(key,0)+power*np.log(obs/max(pred,1e-30))/width**2
        denominator[key]=denominator.get(key,0)+power**2/width**2
    amp={key:float(np.exp(numerator[key]/denominator[key])) for key in numerator}
    rows=[];score=0.
    for group,metric,pred,obs,low,high in items:
        power=2 if metric=='msd' else 0 if metric=='directness' else 1
        value=pred*amp['shared' if shared else group]**power
        score+=(np.log(max(value,1e-30)/obs)/(np.log(high/low)/2))**2
        rows.append(dict(group=group,metric=metric,predicted=float(value),observed=obs,low=low,high=high,
                         within=bool(low<=value<=high)))
    return float(score),amp,rows

def run(index):
    model,family=[(m,f) for m in MODELS for f in FAMILIES][index]
    source=OUT/'scaling'/f'{model}_{family}.json'
    record=json.loads(source.read_text());law=next(x for x in record['laws'] if x['law']==record['selected_law'])
    Ls=lengths(model);pairs=[(a,b)for a in Ls for b in Ls if b>a]
    noise_levels=[1.] if family in ('gaussian','memory') else [.05,.2,.6]
    banks={};audits=[]
    # Time conversion is shared between populations and freely scanned; seconds
    # are not transferred literally between unrelated cell systems.
    kgrid=np.geomspace(.0001,100.,31)
    if family not in ('gaussian','memory'):
        thetas=[parameters(law,L) for L in Ls]
        dt=min(min(.05*t[0],.05*t[1],.05/(t[6]+1e-6))/np.sqrt(max(1.,t[2])) for t in thetas)
        steps=BANK_STEPS
        parameter_banks={}
        for L,theta in zip(Ls,thetas):
            for noise in noise_levels:
                key=(tuple(theta),noise)
                if key not in parameter_banks:
                    parameter_banks[key]=nonlinear_bank(theta,dt,noise,steps,20260925,BURN_STEPS)
                bank,pol=parameter_banks[key]
                banks[L,noise]=bank
                h=len(pol)//2
                v1=float(np.mean((pol[:h]**2).sum(1)));v2=float(np.mean((pol[h:]**2).sum(1)))
                mean=float(np.linalg.norm(pol.mean(0))/np.sqrt(max((v1+v2)/2,1e-30)))
                audits.append(dict(length=L,noise=noise,half_variance_ratio=v2/max(v1,1e-30),
                                   normalized_mean=mean,bank_duration=steps*dt,
                                   passes=bool(.75<=v2/max(v1,1e-30)<=1.333 and mean<.2)))
        saved_dt=16*dt
        # Ensure a nonempty common camera-resolution range; 24h with 0.5h cadence.
        kgrid=np.geomspace(24/(steps*dt/WINDOW_FRACTION),.5/saved_dt,25)

    @lru_cache(maxsize=4096)
    def track(L,noise,k,cad,dur):
        if family in ('gaussian','memory'):
            return linear_tracks(parameters(law,L),family,k,cad,dur,n=512)
        return bank_tracks(banks[L,noise],saved_dt,k,cad,dur,n=512)

    results=[]
    for study,cad,dur in [('Xuan2018',1/6,12.),('Mallin2023',.5,24.),('Jemaa2017',.5,24.),('Jemaa2023',1.,24.)]:
        obs=OBS[(OBS.study==study)&(OBS.role=='fit')]
        best={};feasible={}
        for k in kgrid:
            for noise in noise_levels:
                cache={}
                for L in Ls:
                    inc=track(L,noise,float(k),cad,dur)
                    if inc is None: continue
                    sm=summaries(inc)
                    means={name:float(x.mean()) for name,x in sm.items()}
                    means['speed']=means['path_length']/dur
                    if study=='Xuan2018':
                        lag=float(obs[obs.metric=='msd'].lag_h.iloc[0])
                        if family in ('gaussian','memory'):
                            rates,w=modes(parameters(law,L),family)
                            means['msd']=float(2*k*k*displacement_variance(lag/k,rates,w))
                        else:
                            laginc=bank_tracks(banks[L,noise],saved_dt,float(k),lag,lag,n=1024)
                            if laginc is None: continue
                            means['msd']=float(np.mean(np.sum(laginc[:,0]**2,axis=1)))
                    if family in ('gaussian','memory'):
                        rates,w=modes(parameters(law,L),family)
                        means['speed']=float(np.sqrt(np.pi*k*k*displacement_variance(cad/k,rates,w)/2)/cad)
                        means['path_length']=means['speed']*dur
                        means['endpoint']=float(np.sqrt(np.pi*k*k*displacement_variance(dur/k,rates,w)/2))
                    cache[L]=means
                for ref,poly in pairs:
                    if ref not in cache or poly not in cache: continue
                    valid=(not audits) or all(next(a['passes'] for a in audits if a['length']==L and a['noise']==noise) for L in [ref,poly])
                    items=[]
                    for row in obs.itertuples():
                        L=ref if row.group=='reference' else poly
                        items.append((row.group,row.metric,cache[L][row.metric],row.value,row.low,row.high))
                    for shared in [True,False]:
                        hypothesis='length_only' if shared else 'length_plus_population_amplitude'
                        score,amp,pred=fit_amplitude(items,shared)
                        compatible=all(r['within'] for r in pred)
                        feasible[hypothesis]=feasible.get(hypothesis,0)+int(compatible and valid)
                        if hypothesis not in best or (not valid,score)<(not best[hypothesis]['simulation_valid'],best[hypothesis]['score']):
                            checks=[]
                            for row in OBS[(OBS.study==study)&(OBS.role!='fit')].itertuples():
                                L=ref if row.group=='reference' else poly
                                if row.metric not in cache[L] or row.group not in ['reference','polyploid']: continue
                                power=0 if row.metric=='directness' else 1
                                value=cache[L][row.metric]*amp['shared' if shared else row.group]**power
                                checks.append(dict(metric=row.metric,group=row.group,predicted=value,
                                                   low=row.low,high=row.high,within=bool(row.low<=value<=row.high)))
                            best[hypothesis]=dict(study=study,hypothesis=hypothesis,score=score,
                                reference_length=ref,polyploid_length=poly,k_hours_per_unit=float(k),noise=noise,
                                amplitude=amp,predictions=pred,checks=checks,grid_hit=compatible)
                            best[hypothesis]['simulation_valid']=valid
        for hypothesis,r in best.items():
            r['compatible_grid_points']=feasible[hypothesis]
            r['status']='simulation_unresolved' if not r['simulation_valid'] else 'insufficiently_constrained' if study.startswith('Jemaa') else (
                'compatible_on_grid' if feasible[hypothesis] else 'no_compatible_point_on_tested_grid')
            results.append(r)
        print(model,family,study,[(h,r['score'])for h,r in best.items()],flush=True)
        track.cache_clear()

    # Hosny: common baseline across three control sources; held-out KD conditions
    # test unchanged parameters. A common directed velocity offset is an explicit
    # nuisance hypothesis, not an inferred oxygen-to-polarity mechanism.
    best=None
    control_ids={'N','P','HN0','HP0','RN0','RP0'}
    for k in kgrid[::2]:
        for noise in noise_levels:
            for ref,poly in pairs:
                incs=[track(L,noise,float(k),.5,16.) for L in (ref,poly)]
                if any(x is None for x in incs): continue
                base_sd=np.sqrt(np.mean(incs[0]**2))/.5
                for bias in [0.,.1,.25,.5,1.]:
                    mm={}
                    for group,inc in zip(['non-PACC','PACC'],incs):
                        shifted=inc.copy();shifted[:,:,1]+=.5*bias*base_sd
                        mm[group]={name:np.quantile(x,[.25,.5,.75]) for name,x in summaries(shifted).items()}
                    num=den=0.
                    for row in HOSNY:
                        if row['id'] not in control_ids:continue
                        target=np.array(row['path']);unit=mm[row['group']]['path_length']
                        width=target[2]-target[0]
                        num+=float(unit@target)/width**2;den+=float(unit@unit)/width**2
                    amp=num/den;score=0.;preds=[]
                    for row in HOSNY:
                        if row['group']=='normoxia':continue
                        for metric in ['path','directness','parallel','endpoint','perpendicular']:
                            values=mm[row['group']]['path_length' if metric=='path' else metric]*(amp if metric in ['path','endpoint'] else 1.)
                            target=np.array(row[metric]);mask=np.isfinite(target)
                            width=max(float(target[-1]-target[0]) if np.all(np.isfinite(target[[0,-1]])) else .1,.1)
                            err=float(np.mean(((values[mask]-target[mask])/width)**2))
                            role='fit' if row['id'] in control_ids and metric in ['path','directness','parallel'] else 'held_out'
                            if role=='fit':score+=err/18
                            preds.append(dict(condition=row['id'],metric=metric,role=role,
                                              predicted=values.tolist(),observed=[float(v) if np.isfinite(v) else None for v in target],
                                              scaled_error=err))
                    valid=(not audits) or all(next(a['passes'] for a in audits if a['length']==L and a['noise']==noise) for L in [ref,poly])
                    if best is None or (not valid,score)<(not best['simulation_valid'],best['score']):
                        best=dict(study='Hosny2026',hypothesis='length_plus_shared_directional_offset',score=score,
                            reference_length=ref,polyploid_length=poly,k_hours_per_unit=float(k),noise=noise,
                            bias=bias,amplitude=amp,predictions=preds,simulation_valid=valid,
                            status='descriptive_only_no_sampling_uncertainty' if valid else 'simulation_unresolved')
    if best is not None:results.append(best)
    save_json(OUT/'literature'/f'{model}_{family}.json',dict(model=model,family=family,
              selected_law=law['law'],results=results,stationarity_audit=audits,
              nonlinear_closure='isotropic radial p/q, additive white noise on p only; fixed coefficient across lengths',
              linear_closure='Gaussian eta fixed across lengths; memory substitute covariance with B=tp/2',
              limitations=['Effective lengths are scenario coordinates, not inferred biological sizes.',
                           'Finite grids and conditional nuisance assumptions; no global exclusion or confidence intervals.',
                           'Hosny normoxia is unpaired and excluded; unchanged KD predictions are checks, not fitted mechanisms.',
                           'Nonlinear bank windows overlap; stationarity failures invalidate compatibility claims.'],
              source_hashes=digest([Path(__file__),source,ROOT/'wgd_langevin/data/observations.csv',ROOT/'wgd_langevin/hosny/data.py'])))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int);run(p.parse_args().index)
