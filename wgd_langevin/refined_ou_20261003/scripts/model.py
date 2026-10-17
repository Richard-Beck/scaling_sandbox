"""Descriptive hierarchical Gaussian OU fit with explicit shrinkage scales."""
import json,sys,hashlib
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky,toeplitz
from scipy.special import ndtri
from scipy.stats import qmc
HERE=Path(__file__).resolve().parent;RUN=HERE.parent;OUT=RUN/'results';ROOT=RUN.parents[1]
def displacement_variance(t,rates,weights):
    """Integrated stationary OU covariance; variance of one displacement axis."""
    t=np.asarray(t);out=np.zeros_like(t,dtype=float)
    for rate,w in zip(rates,weights):
        z=rate*t
        shape=np.where(np.abs(z)<1e-3,z*z/2-z**3/6+z**4/24-z**5/120,z+np.expm1(-z))
        out+=2*w*shape/(rate*rate)
    return out

def save_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.json.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(path)

sys.path.insert(0,str(RUN/'data'))
from adapter import observations,curve_observations
import pandas as pd
CURVE=pd.read_csv(RUN/'data/xuan_msd_curve.csv')
CURVE_OBS=[dict(r,metric='msd_curve',objective_weight=r['weight'],scale_interpretation='legacy descriptive log tolerance; pixel bounds are diagnostic only') for r in curve_observations()]
STUDIES=['Xuan2018','Mallin2023','Jemaa2017','Jemaa2023','Hosny_baseline']
SCHEDULE={'Xuan2018':(1/6,12.),'Mallin2023':(.5,24.),'Jemaa2017':(.5,24.),'Jemaa2023':(1.,24.),'Hosny_baseline':(.5,16.)}
FIT_OBS={s:observations(s) for s in STUDIES};CHECK_OBS={s:observations(s,True) for s in STUDIES}
SCENARIOS={'primary':1.,'prior5':1.,'prior20':1.,'strong':.5,'weak':2.}
PRIOR_MINUTES={s:(5. if s=='prior5' else 20. if s=='prior20' else 10.) for s in SCENARIOS}
PRIOR_LOG_SD=np.log(10.)/1.959963984540054
SCENARIO_LABELS={s:f'{PRIOR_MINUTES[s]:g}-minute prior; '+('strong' if s=='strong' else 'weak' if s=='weak' else 'primary')+' pooling' for s in SCENARIOS}

def persistence_residuals(z,scenario):
    logs=np.asarray(z[:20]).reshape(5,2,2)[:,:,0]
    return ((logs-np.log(PRIOR_MINUTES[scenario]/60.))/PRIOR_LOG_SD).ravel()
BASE_SD=np.log(10.);EFFECT_SD=np.log(2.)
# z[study, group, parameter]: log(tau in h), log(sigma in um/h^(3/2)).
# Final two coordinates are signed Hosny mean velocities (um/h).
LOW=np.r_[np.tile([np.log(1e-6),np.log(1e-4)],10),[-200.,-200.]]
HIGH=np.r_[np.tile([np.log(12.-1e-8),np.log(1e8)],10),[200.,200.]]

def provenance():
    paths=[HERE/'model.py',HERE/'fit.py',RUN/'data/adapter.py',RUN/'data/observations.csv',RUN/'data/xuan_msd_curve.csv',RUN/'data/hosny_observations.csv']
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def decode(z):
    logpars=np.asarray(z[:20]).reshape(5,2,2)
    return np.exp(logpars),np.asarray(z[20:22]),logpars

def hyperparameters(z,scale=1.):
    pars,drifts,logs=decode(z);base=logs[:,0,:];effects=logs[:,1,:]-base
    output={}
    for j,name in enumerate(['tau','sigma']):
        output[name]=dict(population_log_baseline=float(base[:,j].mean()),population_baseline=float(np.exp(base[:,j].mean())),
            population_log_wgd_effect=float(effects[:,j].mean()),population_wgd_fold=float(np.exp(effects[:,j].mean())),
            fixed_log_baseline_sd=float(BASE_SD*scale),fixed_log_effect_sd=float(EFFECT_SD*scale),
            fitted_baseline_log_spread=float(base[:,j].std(ddof=1)),fitted_effect_log_spread=float(effects[:,j].std(ddof=1)))
    return output

def pooling_residuals(z,scale):
    _,_,logs=decode(z);base=logs[:,0,:];effects=logs[:,1,:]-base
    return np.r_[((base-base.mean(axis=0))/(BASE_SD*scale)).ravel(),((effects-effects.mean(axis=0))/(EFFECT_SD*scale)).ravel()]

@lru_cache(maxsize=12)
def normal_bank(frames,n,seed):
    sampler=qmc.Sobol(d=2*frames,scramble=True,seed=seed)
    uniform=sampler.random_base2(int(np.log2(n)))
    assert len(uniform)==n
    z=ndtri(np.clip(uniform,1e-12,1-1e-12)).reshape(n,frames,2)
    return z.transpose(1,0,2).reshape(frames,n*2)

@lru_cache(maxsize=96)
def unit_increments(tau,cad,duration,n,seed):
    frames=round(duration/cad)
    v=displacement_variance(np.arange(frames+1)*cad,np.array([1/tau]),np.array([tau/2]))
    first=np.r_[v[1],.5*(v[2:]+v[:-2]-2*v[1:-1])]
    cov=toeplitz(first);root=cholesky(cov+np.eye(frames)*max(first[0]*1e-12,1e-300),lower=True)
    return (root@normal_bank(frames,n,seed)).reshape(frames,n,2).transpose(1,0,2)

@lru_cache(maxsize=96)
def directness(tau,cad,duration,n,seed):
    inc=unit_increments(tau,cad,duration,n,seed)
    path=np.linalg.norm(inc,axis=2).sum(axis=1);end=np.linalg.norm(inc.sum(axis=1),axis=1)
    return float(np.mean(end/path))

def group_prediction(study,tau,sigma,drift,n,seed):
    cad,duration=SCHEDULE[study];rates=np.array([1/tau]);w=np.array([tau/2])
    stepvar=float(displacement_variance(cad,rates,w))*sigma**2
    p=dict(speed=float(np.sqrt(np.pi/2*stepvar)/cad))
    p['path_length']=p['speed']*duration
    p['endpoint']=float(np.sqrt(np.pi/2*displacement_variance(duration,rates,w))*sigma)
    p['msd']=float(2*displacement_variance(11.1111111111,rates,w)*sigma**2)
    p['early_msd']=float(2*stepvar)
    if study=='Xuan2018':
        lags=np.sort(CURVE.lag_h.unique());p['curve_lags_h']=lags.tolist()
        p['msd_curve']=(2*displacement_variance(lags,rates,w)*sigma**2).tolist()
    if study=='Mallin2023':p['directness']=directness(float(tau),cad,duration,n,seed)
    if study=='Hosny_baseline':
        inc=unit_increments(float(tau),cad,duration,n,seed)*sigma
        inc[:,:,1]+=drift*cad
        path=np.linalg.norm(inc,axis=2).sum(axis=1);net=inc.sum(axis=1);end=np.linalg.norm(net,axis=1)
        sm=dict(path_length=path,endpoint=end,directness=end/path,parallel=net[:,1]/path,perpendicular=net[:,0]/path)
        p={k:float(v.mean()) for k,v in sm.items()}
        p['speed']=p['path_length']/duration
        p['quantiles']={k:np.quantile(v,[.25,.5,.75]).tolist() for k,v in sm.items()}
    return p

def predict(z,n=1024,seed=19431):
    pars,drifts,_=decode(z)
    return {s:[group_prediction(s,*pars[i,g],drifts[g] if s=='Hosny_baseline' else 0.,n,seed+100*i) for g in (0,1)] for i,s in enumerate(STUDIES)}

def data_residuals(pred,checks=False):
    residuals=[];rows=[]
    for study in STUDIES:
        study_obs=(CHECK_OBS if checks else FIT_OBS)[study]
        if study=='Xuan2018' and not checks:study_obs=study_obs+CURVE_OBS
        for o in study_obs:
            g=o['group'];metric=o['metric'];p=pred[study][g]
            value=p['quantiles'][metric][o['quantile']] if 'quantile' in o else float(np.interp(o['lag_h'],p['curve_lags_h'],p['msd_curve'])) if metric=='msd_curve' else p['early_msd'] if checks and study=='Xuan2018' and metric=='msd' else p[metric]
            if study=='Hosny_baseline':res=(value-o['obs'])/o['width'];within=None
            else:res=np.log(max(value,1e-300)/o['obs'])/o.get('log_scale',np.log(o['high']/o['low'])/2);within=bool(o['low']<=value<=o['high'])
            raw_res=float(res);res*=np.sqrt(o.get('objective_weight',1.))
            residuals.append(res);rows.append(dict(study=study,**o,predicted=float(value),residual=float(res),unweighted_residual=raw_res,within=within))
    return np.array(residuals),rows

def objective(z,scale,n=1024,seed=19431):
    r,_=data_residuals(predict(z,n,seed))
    return np.r_[r,pooling_residuals(z,scale)]

def initial(index):
    # Generic physical seeds anchored only to observed speed/path, not former fits.
    times=np.array([[.12,.8],[.08,1.8],[.5,1.],[.5,1.],[.8,1.8]])
    rng=np.random.default_rng(202609290+index)
    if index:
        times*=np.exp(rng.normal(0,1.4,size=(5,2)))
        if index>=4:times=np.exp(rng.uniform(np.log(.01),np.log(6.),size=(5,2)))
    times=np.clip(times,1e-5,10.)
    logpars=np.zeros((5,2,2))
    speeds={'Xuan2018':[23.,14.75],'Mallin2023':[305/24,165/24],'Jemaa2017':[16.,27.],'Jemaa2023':[2.75,7.25],'Hosny_baseline':[105/16,175/16]}
    for i,study in enumerate(STUDIES):
        cad,_=SCHEDULE[study]
        for g in (0,1):
            tau=times[i,g];unit=np.sqrt(np.pi/2*displacement_variance(cad,np.array([1/tau]),np.array([tau/2])))/cad
            sigma=speeds[study][g]/unit
            logpars[i,g]=np.log([tau,sigma])
    z=np.r_[logpars.ravel(),[1.,2.5]]
    if index:z[20:]+=rng.normal(0,1,2)
    return np.clip(z,LOW+1e-8,HIGH-1e-8)

def record(z,scenario,n,seed=19431):
    scale=SCENARIOS[scenario];pred=predict(z,n,seed);res,rows=data_residuals(pred);pool=pooling_residuals(z,scale)
    pars,drifts,_=decode(z);parameters=[]
    for i,study in enumerate(STUDIES):
        parameters.append(dict(study=study,tau_reference=float(pars[i,0,0]),tau_wgd=float(pars[i,1,0]),tau_fold=float(pars[i,1,0]/pars[i,0,0]),
            sigma_reference=float(pars[i,0,1]),sigma_wgd=float(pars[i,1,1]),sigma_fold=float(pars[i,1,1]/pars[i,0,1]),
            drift_reference=float(drifts[0]) if study=='Hosny_baseline' else 0.,drift_wgd=float(drifts[1]) if study=='Hosny_baseline' else 0.,
            drift_change=float(drifts[1]-drifts[0]) if study=='Hosny_baseline' else 0.))
    return dict(z=np.asarray(z).tolist(),scenario=scenario,pooling_scale=scale,hyperparameters=hyperparameters(z,scale),parameters=parameters,
        data_score=float(res@res),pooling_penalty=float(pool@pool),total_score=float(res@res+pool@pool),rows=rows,n=n,seed=seed)
