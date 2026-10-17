"""Independent-seed and timestep sensitivity of selected literature scenarios."""
import argparse
from functools import lru_cache
import numpy as np
from design import *
from scaling import parameters
from literature import linear_tracks,nonlinear_bank,bank_tracks,summaries,modes,displacement_variance,BANK_STEPS,BURN_STEPS

def run(index):
    model,family=[(m,f) for m in MODELS for f in FAMILIES][index]
    result=json.loads((OUT/'literature'/f'{model}_{family}.json').read_text())
    scaling=json.loads((OUT/'scaling'/f'{model}_{family}.json').read_text())
    law=next(r for r in scaling['laws'] if r['law']==scaling['selected_law'])
    theta_all=[parameters(law,L)for L in lengths(model)]
    dt=min(min(.05*t[0],.05*t[1],.05/(t[6]+1e-6))/np.sqrt(max(1.,t[2]))for t in theta_all)
    @lru_cache(maxsize=30)
    def bank_by_theta(theta_tuple,noise,refinement):
        # Same seed is not a common-Brownian-path refinement: changes combine
        # time-step and finite-bank variability, explicitly reported as such.
        x,p=nonlinear_bank(np.array(theta_tuple),dt/refinement,noise,BANK_STEPS*refinement,91531,BURN_STEPS*refinement)
        h=len(p)//2;v1=np.mean((p[:h]**2).sum(1));v2=np.mean((p[h:]**2).sum(1))
        stationarity=bool(.75<=v2/max(v1,1e-30)<=1.333 and np.linalg.norm(p.mean(0))/np.sqrt(max((v1+v2)/2,1e-30))<.2)
        return x,16*dt/refinement,stationarity
    def bank(L,noise,refinement):
        return bank_by_theta(tuple(parameters(law,L)),noise,refinement)
    checks=[]
    for r in result['results']:
        if r['study']=='Hosny2026':
            k=r['k_hours_per_unit'];incs=[];stationary=True
            for L in [r['reference_length'],r['polyploid_length']]:
                if family in ['gaussian','memory']:
                    inc=linear_tracks(parameters(law,L),family,k,.5,16.,n=16384,seed=91531)
                else:
                    x,saved_dt,valid=bank(L,r['noise'],1);stationary &= valid
                    inc=bank_tracks(x,saved_dt,k,.5,16.,n=4096,seed=95131)
                incs.append(inc)
            if any(inc is None for inc in incs):continue
            sd=np.sqrt(np.mean(incs[0]**2))/.5
            for group,inc in zip(['non-PACC','PACC'],incs):
                inc[:,:,1]+=.5*r['bias']*sd
                quant={name:np.quantile(x,[.25,.5,.75])for name,x in summaries(inc).items()}
                for old in r['predictions']:
                    if (group=='PACC') != ('P' in old['condition']):continue
                    metric=old['metric'];v=quant['path_length' if metric=='path' else metric]
                    v=v*(r['amplitude'] if metric in ['path','endpoint']else 1.)
                    for q,value,original in zip([25,50,75],v,old['predicted']):
                        checks.append(dict(study='Hosny2026',hypothesis=r['hypothesis'],group=group,
                            metric=f"{old['condition']}:{metric}:q{q}",refinement=1,original=original,
                            independent=float(value),within=None,stationary=stationary,
                            relative_difference=float(abs(value-original)/max(abs(original),1e-12))))
            continue
        if r['study'] not in ['Xuan2018','Mallin2023']:continue
        cad,dur=(1/6,12.) if r['study']=='Xuan2018' else (.5,24.)
        k=r['k_hours_per_unit'];shared=r['hypothesis']=='length_only'
        for group,L in [('reference',r['reference_length']),('polyploid',r['polyploid_length'])]:
            amp=r['amplitude']['shared' if shared else group]
            for refinement in ([1] if family in ['gaussian','memory']else [1,2]):
                if family in ['gaussian','memory']:
                    inc=linear_tracks(parameters(law,L),family,k,cad,dur,n=16384,seed=91531)
                    stationary=True
                else:
                    x,saved_dt,stationary=bank(L,r['noise'],refinement)
                    inc=bank_tracks(x,saved_dt,k,cad,dur,n=4096,seed=95131)
                if inc is None:continue
                means={name:float(a.mean())for name,a in summaries(inc).items()}
                means['speed']=means['path_length']/dur
                if r['study']=='Xuan2018':
                    lag=11.1111111111
                    if family in ['gaussian','memory']:
                        rates,w=modes(parameters(law,L),family)
                        means['msd']=float(2*k*k*displacement_variance(lag/k,rates,w))
                    else:
                        laginc=bank_tracks(x,saved_dt,k,lag,lag,n=4096,seed=95131)
                        if laginc is None:continue
                        means['msd']=float(np.mean(np.sum(laginc[:,0]**2,axis=1)))
                for old in r['predictions']:
                    if old['group']!=group:continue
                    metric=old['metric'];power=2 if metric=='msd' else 0 if metric=='directness' else 1
                    value=means[metric]*amp**power
                    checks.append(dict(study=r['study'],hypothesis=r['hypothesis'],group=group,metric=metric,
                           refinement=refinement,original=old['predicted'],independent=value,
                           within=bool(old['low']<=value<=old['high']),stationary=stationary,
                           relative_difference=abs(value-old['predicted'])/max(abs(old['predicted']),1e-12)))
    save_json(OUT/'literature_validation'/f'{model}_{family}.json',dict(model=model,family=family,checks=checks,
              interpretation='Independent realization; nonlinear half-step differences also include finite-bank variation.',
              source_hashes=digest([Path(__file__),HERE/'literature.py'])))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int);run(p.parse_args().index)
