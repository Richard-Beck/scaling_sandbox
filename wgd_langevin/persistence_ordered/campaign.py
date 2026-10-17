"""All four WGD variants, positive times below 12 h, WGD+ times strictly larger."""
import sys,os,json,time,hashlib,argparse
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.special import expit,logit
from scipy.optimize import least_squares
from scipy.linalg import cholesky,toeplitz
HERE=Path(__file__).resolve().parent;OUT=HERE/'results';ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'persistence_only'))
import importlib.util
_spec=importlib.util.spec_from_file_location('wgd_original_campaign',HERE.parent/'persistence_only/campaign.py')
core=importlib.util.module_from_spec(_spec);sys.modules[_spec.name]=core;_spec.loader.exec_module(core)
FLOOR=1e-6;CAP=12.;MAX_STEPS_PER_CELL=2_000_000
STARTS=8;SEARCH_EVALS=120;POLISH_EVALS=80

def jobs():return core.jobs()
def key(f,mask,study):return f+'_'+('ordered' if mask else 'shared')+'_'+study

def specification(f,mask,study):
    names,lo,hi=core.specification(f,(),study)
    for name in ('tp','tm'):
        if name in names:
            i=names.index(name);lo[i]=np.log(FLOOR);hi[i]=np.log(CAP-1e-8)
    names=list(names);lo=list(lo);hi=list(hi)
    if mask:
        for name in ('tp','tm'):
            if name in names:names.append(name+'_gap');lo.append(-12.);hi.append(12.)
    return names,np.array(lo),np.array(hi)

def physical(z,names):
    values={k:float(np.exp(v)) for k,v in zip(names,z) if not k.endswith('_gap')}
    for k,v in zip(names,z):
        if k.endswith('_gap'):
            name=k[:-4];a=values[name];values[name+'_poly']=float(a+(CAP-a)*expit(v))
    return values

def theta_values(z,names):
    d=physical(z,names);pn=list(d);pz=np.log(list(d.values()))
    return d,pn,pz,[core.decode(pz,pn,g) for g in (0,1)]

class NumericalUnresolved(RuntimeError):pass

def predict(z,names,f,study,n=96,seed=431,dt=.01,burn=48.,lag_override=None):
    d,pn,pz,decoded=theta_values(z,names)
    cad,dur,lag=(1/6,12.,11.1111111111) if study=='Xuan2018' else (.5,16.,1.) if study.startswith('Hosny') else (1. if study=='Jemaa2023' else .5,24.,1.)
    if lag_override is not None:lag=lag_override
    linear=f in ('gaussian','memory')
    if not linear:
        # Shared grid for both populations; refine dyadically for explicit stability.
        stability=min(min(.2*t[0],.2*t[1],.15/(t[6]/max(t[7],1e-12)+1e-12)) for t,eta,force in decoded)
        while dt>stability:dt/=2
        steps=int(np.ceil(cad/dt));actual=cad/steps
        if int(np.ceil(burn/actual))+round(dur/cad)*steps>MAX_STEPS_PER_CELL:
            raise NumericalUnresolved('Integration work cap: short-time candidate not evaluated; not biological infeasibility')
    output=[]
    for theta,eta,force in decoded:
        if linear:
            rates,w=core.modes(theta,f)
            mean=theta[0]*force/(1-theta[2] if f=='memory' else 1.)
            msd=float(2*core.displacement_variance(lag,rates,w)+(mean*lag)**2)
            if study in ('Xuan2018','Jemaa2017','Jemaa2023'):
                # Only analytic speed and MSD enter these objectives/checks.
                output.append(dict(speed=float(np.sqrt(np.pi/2*core.displacement_variance(cad,rates,w))/cad),msd=msd))
                continue
            v=core.displacement_variance(np.arange(round(dur/cad)+1)*cad,rates,w)
            first=np.r_[v[1],.5*(v[2:]+v[:-2]-2*v[1:-1])]
            cov=toeplitz(first)
            root=cholesky(cov+np.eye(len(first))*max(first[0]*1e-12,1e-300),lower=True)
            rng=np.random.default_rng(seed)
            inc=(root@rng.normal(size=(len(first),n*2))).reshape(len(first),n,2).transpose(1,0,2)
            inc[:,:,1]+=mean*cad
        else:
            inc,lag2=core.ensemble(theta,eta,force,cad,dur,lag,n,seed,dt,burn)
            if not np.isfinite(inc).all() or not np.isfinite(lag2).all():raise NumericalUnresolved('Non-finite nonlinear simulation')
            msd=float(lag2.mean())
        sm=core.summaries(inc);p={k:float(v.mean()) for k,v in sm.items()}
        p.update(speed=p['path_length']/dur,msd=msd)
        if linear and force==0:
            p['speed']=float(np.sqrt(np.pi/2*core.displacement_variance(cad,rates,w))/cad)
            p['path_length']=p['speed']*dur
            p['endpoint']=float(np.sqrt(np.pi/2*core.displacement_variance(dur,rates,w)))
        p['quantiles']={k:np.quantile(v,[.25,.5,.75]).tolist() for k,v in sm.items()}
        output.append(p)
    return output

def provenance():
    paths=[Path(__file__),HERE.parent/'persistence_only/campaign.py',ROOT/'wgd_langevin/data/observations.csv',ROOT/'wgd_langevin/hosny/data.py']
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def make_record(z,names,f,mask,study,obs,n,starts,stage,tick,rejections):
    pred=predict(z,names,f,study,n=n);res,amp,rows,feasible=core.score(pred,obs,study)
    d,pn,pz,decoded=theta_values(z,names);groups=[]
    for theta,eta,force in decoded:
        groups.append(dict(tp=float(theta[0]),tm=float(theta[1]),r=float(theta[2]),ps=float(theta[3]),vs=float(theta[4]),self_a=float(theta[6]),self_K=float(theta[7]),eta=float(eta),force=float(force),amplitude=amp))
    a,b=groups
    assert all(a[k]==b[k] for k in a if k not in ('tp','tm'))
    times=['tp'] if f=='gaussian' else ['tp','tm']
    assert all(0<a[k]<CAP and 0<b[k]<CAP for k in times)
    assert all(b[k]>a[k] for k in times) if mask else a==b
    return dict(key=key(f,mask,study),family=f,study=study,hypothesis='ordered_times' if mask else 'shared_null',stage=stage,complete=False,
        generated=datetime.now(timezone.utc).isoformat(),parameters=d,shared_amplitude=amp,ratios={k:b[k]/a[k] for k in times},score=float(res@res),rows=rows,
        all_reading_intervals=all(o['within'] for o in rows) if not study.startswith('Hosny') else None,
        amplitude_interval_feasible=feasible,coordinate_names=names,z=z.tolist(),starts=starts,starts_completed=len(starts),
        constraint_audit=dict(passes=True,group_parameters=groups),numerically_unresolved_evaluations=rejections,
        numerical_time_floor_h=FLOOR,time_upper_open_h=CAP,slurm_job=os.environ.get('SLURM_JOB_ID'),wall_seconds=time.monotonic()-tick,
        source_hashes=provenance(),sample_n=n)

def run(index):
    f,mask,study=jobs()[index];dest=OUT/(key(f,mask,study)+'.json');tick=time.monotonic()
    if dest.exists():
        old=json.loads(dest.read_text())
        if old['source_hashes']!=provenance():raise RuntimeError('Source changed; use a fresh output directory')
        if old['complete']:return
    names,lo,hi=specification(f,mask,study);obs=core.observations(study);linear=f in ('gaussian','memory')
    n=512 if linear else 96;polish_n=2048 if linear else 512
    default=dict(tp=.6,tm=2.,r=.6 if f=='memory' else 2.,vs=.5,eta=.5,self_a=.1,self_K=.2,force=.1)
    base=np.array([-4. if k.endswith('_gap') else np.log(default[k]) for k in names])
    rng=np.random.default_rng(202609280+index);seeds=[base]
    for i in range(STARTS-1):
        z=base+rng.normal(0,1.4,len(names))
        for j,k in enumerate(names):
            if k in ('tp','tm'):z[j]=rng.uniform(np.log(.025),np.log(8.))
            elif k.endswith('_gap'):z[j]=rng.uniform(-8,3)
        seeds.append(np.clip(z,lo+1e-6,hi-1e-6))
    if mask:
        nullpath=OUT/(key(f,(),study)+'.json')
        null=json.loads(nullpath.read_text())
        if not null['complete']:raise RuntimeError('Shared null has not completed')
        if null['source_hashes']!=provenance():raise RuntimeError('Null provenance differs')
        seeds[0]=np.array([-11.9 if k.endswith('_gap') else np.log(null['parameters'][k]) for k in names])
    rejections=0;starts=[];best=None;best_value=np.inf
    def objective(z,count=n):
        nonlocal rejections
        try:
            rr=core.score(predict(z,names,f,study,n=count),obs,study)[0]
            if not np.isfinite(rr).all():raise NumericalUnresolved('Nonfinite residual')
            return rr
        except (NumericalUnresolved,np.linalg.LinAlgError,FloatingPointError,ValueError):
            rejections+=1;return np.full(len(obs),1e6)
    # Save a valid initial checkpoint before long optimization calls.
    initial=make_record(base,names,f,mask,study,obs,n,starts,'initial',tick,rejections)
    core.save_json(dest,initial)
    best=base.copy();best_value=initial['score']
    for i,z in enumerate(seeds):
        z=np.clip(z,lo+1e-7,hi-1e-7)
        fit=least_squares(objective,z,bounds=(lo,hi),max_nfev=SEARCH_EVALS,diff_step=.003 if not linear else .0001,ftol=1e-5,xtol=1e-5,gtol=1e-5)
        value=float(2*fit.cost)
        starts.append(dict(index=i,z=fit.x.tolist(),score=value,success=bool(fit.success),nfev=fit.nfev,numerically_valid=value<1e11))
        if value<best_value:best=fit.x.copy();best_value=value
        core.save_json(dest,make_record(best,names,f,mask,study,obs,n,starts,'search',tick,rejections))
        print('search',i,value,fit.success,flush=True)
    candidates=[np.array(s['z']) for s in sorted(starts,key=lambda x:x['score'])[:2]]+[best]
    ranked=[]
    for i,z in enumerate(candidates):
        residual=objective(z,polish_n);ranked.append((float(residual@residual),i,z))
    ranked.sort(key=lambda x:x[0])
    best_value,_,best=ranked[0];best=best.copy();polishes=[]
    if best_value>=1e11:raise RuntimeError('No numerically evaluable polish candidate; search checkpoint retained')
    # Polish the two best candidates at the larger ensemble, checkpoint each result.
    for _,i,z in ranked[:2]:
        fit=least_squares(lambda zz:objective(zz,polish_n),np.clip(z,lo+1e-8,hi-1e-8),bounds=(lo,hi),max_nfev=POLISH_EVALS,diff_step=.003 if not linear else .0001,ftol=1e-5,xtol=1e-5,gtol=1e-5)
        value=float(2*fit.cost);polishes.append(dict(score=value,success=bool(fit.success),nfev=fit.nfev))
        if value<best_value:best=fit.x.copy();best_value=value
        core.save_json(dest,make_record(best,names,f,mask,study,obs,polish_n,starts,'polish',tick,rejections))
        print('polish',i,value,flush=True)
    record=make_record(best,names,f,mask,study,obs,polish_n,starts,'validation',tick,rejections);record['polishes']=polishes
    validation_n=8192 if linear else 2048;validations=[]
    for seed,dt,burn in ([(937,.01,48.)] if linear else [(937,.01,48.),(937,.005,48.),(937,.005,96.)]):
        try:
            pred=predict(best,names,f,study,n=validation_n,seed=seed,dt=dt,burn=burn)
            rr,_,rows,feasible=core.score(pred,obs,study,record['shared_amplitude'])
            validations.append(dict(seed=seed,max_dt=dt,burn=burn,n=validation_n,score=float(rr@rr),rows=rows,numerically_resolved=True))
        except NumericalUnresolved as exc:validations.append(dict(seed=seed,max_dt=dt,burn=burn,numerically_resolved=False,reason=str(exc)))
        record['validations']=validations;core.save_json(dest,record)
    checks=core.observations(study,True)
    if checks:
        lag=next((o['lag_h'] for o in checks if o.get('lag_h') is not None),None)
        try:
            pred=predict(best,names,f,study,n=validation_n,seed=1937,dt=.005,burn=96.,lag_override=lag)
            record['checks']=core.score(pred,checks,study,record['shared_amplitude'])[2]
        except NumericalUnresolved as exc:record['checks_unresolved']=str(exc)
    record.update(complete=True,stage='complete',generated=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-tick,
                  numerically_unresolved_evaluations=rejections)
    core.save_json(dest,record)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int,nargs='?');p.add_argument('--manifest',action='store_true');args=p.parse_args()
    if args.manifest:
        core.save_json(OUT/'manifest.json',dict(generated=datetime.now(timezone.utc).isoformat(),jobs=[dict(index=i,family=f,ordered=bool(mask),study=s,key=key(f,mask,s)) for i,(f,mask,s) in enumerate(jobs())],source_hashes=provenance(),
            constraints='0 < tp− < tp+ < 12 h; 0 < tm− < tm+ < 12 h when present; all other parameters shared. Equal-time nulls separate.',
            numerical_floor_h=FLOOR,starts=STARTS,search_evaluations=SEARCH_EVALS,polish_evaluations=POLISH_EVALS,reports=False))
    else:run(args.index)
