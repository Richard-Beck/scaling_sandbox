"""Direct WGD fits: group differences restricted to persistence times."""
import sys,json,argparse,time,hashlib,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.linalg import cholesky,toeplitz
from numba import njit
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'results'
FAMILIES=['gaussian','memory','coupled','self_activation']
HOSNY_PAIRS={
 'Hosny_baseline':('N','P'),
 'Hosny_HIF_control':('HN0','HP0'),
 'Hosny_HIF_KD':('HN1','HP1'),
 'Hosny_RhoA_control':('RN0','RP0'),
 'Hosny_RhoA_KD':('RN1','RP1'),
}
STUDIES=['Xuan2018','Mallin2023','Jemaa2017','Jemaa2023',*HOSNY_PAIRS]
OBS=pd.read_csv(ROOT/'wgd_langevin/data/observations.csv')
sys.path.insert(0,str(ROOT/'wgd_langevin/hosny'))
from data import DATA as HOSNY

def save_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.json.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(path)

def jobs():
    return [(f,mask,s) for varying in (False,True) for f in FAMILIES for s in STUDIES
            for mask in [(() if not varying else (0,) if f=='gaussian' else (0,1))]]

def jobkey(f,mask,study):return f+'_'+('persistence' if mask else 'shared')+'_'+study

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

def specification(f,mask,study):
    names=['tp'];lo=[.05];hi=[12.]
    if f!='gaussian':names+=['tm','r'];lo += [.05,.005];hi += [12.,.99 if f=='memory' else 8.]
    if f in ('coupled','self_activation'):names+=['vs','eta'];lo += [.03,.03];hi += [10.,5.]
    if f=='self_activation':names+=['self_a','self_K'];lo += [.0001,.02];hi += [2.,2.]
    if study.startswith('Hosny'):names+=['force'];lo += [.0001];hi += [10.]
    split={0:'tp',1:'tm'}
    for i in mask:
        name=split[i]
        if name in names:
            j=names.index(name);names.append(name+'_poly');lo.append(lo[j]);hi.append(hi[j])
    return names,np.log(lo),np.log(hi)

def decode(z,names,group):
    d=dict(zip(names,np.exp(z)))
    if group:
        for k in list(d):
            if k.endswith('_poly'):d[k[:-5]]=d[k]
    theta=np.array([d['tp'],d.get('tm',1.),d.get('r',0.),1.,d.get('vs',1.),1.,d.get('self_a',0.),d.get('self_K',.1)])
    return theta,d.get('eta',1.),d.get('force',0.)

@njit(cache=True)
def ensemble(theta,eta,force,cad,duration,lag,n,seed,dtmax=.01,burn=48.):
    # Independent cells, common random numbers across optimization and groups.
    np.random.seed(seed)
    tp,tm,r,ps,vs,g,a,K=theta
    steps=int(np.ceil(cad/dtmax));dt=cad/steps;nb=int(np.ceil(burn/dt));frames=int(round(duration/cad))
    inc=np.zeros((n,frames,2));lag2=np.zeros(n)
    for cell in range(n):
        px=0.;py=0.;qx=0.;qy=0.;xx=0.;xy=0.;lastx=0.;lasty=0.
        for j in range(nb+frames*steps):
            pn=np.sqrt(px*px+py*py);qn=np.sqrt(qx*qx+qy*qy)
            dx=(-px+r*qx/(1+qn/vs))/tp+a*px*pn/(K*K+pn*pn)
            dy=(-py+r*qy/(1+qn/vs))/tp+a*py*pn/(K*K+pn*pn)+force
            ex=(px/(1+pn)-qx)/tm;ey=(py/(1+pn)-qy)/tm
            nx=eta*np.sqrt(dt)*np.random.normal();ny=eta*np.sqrt(dt)*np.random.normal()
            ux=px+dt*dx+nx;uy=py+dt*dy+ny;vx=qx+dt*ex;vy=qy+dt*ey
            un=np.sqrt(ux*ux+uy*uy);vn=np.sqrt(vx*vx+vy*vy)
            dx2=(-ux+r*vx/(1+vn/vs))/tp+a*ux*un/(K*K+un*un)
            dy2=(-uy+r*vy/(1+vn/vs))/tp+a*uy*un/(K*K+un*un)+force
            newx=px+.5*dt*(dx+dx2)+nx;newy=py+.5*dt*(dy+dy2)+ny
            qx+=.5*dt*(ex+(ux/(1+un)-vx)/tm);qy+=.5*dt*(ey+(uy/(1+un)-vy)/tm)
            if j>=nb:
                k=j-nb;oldx=xx;oldy=xy
                xx+=.5*dt*(px+newx);xy+=.5*dt*(py+newy)
                if k*dt < lag <= (k+1)*dt:
                    frac=(lag-k*dt)/dt;lx=oldx+frac*(xx-oldx);ly=oldy+frac*(xy-oldy);lag2[cell]=lx*lx+ly*ly
                if (k+1)%steps==0:
                    frame=(k+1)//steps-1;inc[cell,frame,0]=xx-lastx;inc[cell,frame,1]=xy-lasty;lastx=xx;lasty=xy
            px=newx;py=newy
    return inc,lag2

def observations(study,checks=False):
    if not study.startswith('Hosny'):
        return [dict(group=int(r.group=='polyploid'),metric=r.metric,obs=float(r.value),
                low=float(r.low),high=float(r.high),lag_h=float(r.lag_h) if r.metric=='msd' else None,
                power=2 if r.metric=='msd' else 0 if r.metric=='directness' else 1)
                for r in OBS[(OBS.study==study)&(OBS.role==('check' if checks else 'fit'))].itertuples()
                if r.group in ('reference','polyploid') and r.metric in ('speed','msd','path_length','directness','endpoint')]
    rows=[]
    for r in HOSNY:
        if r['id'] not in HOSNY_PAIRS[study]:continue
        for metric in (('endpoint','perpendicular') if checks else ('path','directness','parallel')):
            q=np.array(r[metric]);width=max(float(q[2]-q[0]),.04) if np.isfinite(q[[0,2]]).all() else .1
            for j in range(3):
                if np.isfinite(q[j]):rows.append(dict(group=int(r['group']=='PACC'),metric='path_length' if metric=='path' else metric,
                    quantile=j,obs=float(q[j]),width=width,power=int(metric in ('path','endpoint')),block=r['id']))
    return rows

def predictions(z,names,f,study,n=96,seed=431,dt=.01,burn=48.,lag_override=None):
    cad,dur,lag=(1/6,12.,11.1111111111) if study=='Xuan2018' else (.5,16.,1.) if study.startswith('Hosny') else (1. if study=='Jemaa2023' else .5,24.,1.)
    if lag_override is not None:lag=lag_override
    out=[]
    for group in range(2):
        theta,eta,force=decode(z,names,group)
        if f in ('gaussian','memory'):
            rates,w=modes(theta,f);v=displacement_variance(np.arange(round(dur/cad)+1)*cad,rates,w)
            first=np.r_[v[1],.5*(v[2:]+v[:-2]-2*v[1:-1])]
            root=cholesky(toeplitz(first)+np.eye(len(first))*1e-12,lower=True)
            rng=np.random.default_rng(seed);inc=np.einsum('ij,njk->nik',root,rng.normal(size=(n,len(first),2)))
            mean=theta[0]*force/(1-theta[2] if f=='memory' else 1.)
            inc[:,:,1]+=mean*cad
            msd=float(2*displacement_variance(lag,rates,w)+(mean*lag)**2)
        else:inc,lag2=ensemble(theta,eta,force,cad,dur,lag,n,seed,dt,burn);msd=float(lag2.mean())
        s=summaries(inc)
        pred={k:float(v.mean()) for k,v in s.items()};pred['speed']=pred['path_length']/dur;pred['msd']=msd
        if f in ('gaussian','memory') and force==0:
            pred['speed']=float(np.sqrt(np.pi/2*displacement_variance(cad,rates,w))/cad)
            pred['path_length']=pred['speed']*dur
            pred['endpoint']=float(np.sqrt(np.pi/2*displacement_variance(dur,rates,w)))
        pred['quantiles']={k:np.quantile(v,[.25,.5,.75]).tolist() for k,v in s.items()}
        out.append(pred)
    return out

def score(pred,obs,study,amplitude=None):
    raw=np.array([pred[r['group']]['quantiles'][r['metric']][r['quantile']] if 'quantile' in r else pred[r['group']][r['metric']] for r in obs])
    power=np.array([r['power'] for r in obs]);target=np.array([r['obs'] for r in obs])
    if study.startswith('Hosny'):
        width=np.array([r['width'] for r in obs]);sel=power==1
        amp=max(1e-12,float(np.sum(raw[sel]*target[sel]/width[sel]**2)/np.sum((raw[sel]/width[sel])**2))) if amplitude is None else amplitude
        values=raw*amp**power;res=(values-target)/width
    else:
        width=np.array([np.log(r['high']/r['low'])/2 for r in obs])
        amp=float(np.exp(np.sum(power*np.log(target/np.maximum(raw,1e-30))/width**2)/np.sum(power**2/width**2))) if amplitude is None else amplitude
        values=raw*amp**power;res=np.log(np.maximum(values,1e-30)/target)/width
    rows=[dict(**r,predicted=float(v),within=bool(r.get('low',-np.inf)<=v<=r.get('high',np.inf)) if not study.startswith('Hosny') else None)for r,v in zip(obs,values)]
    feasible=None
    if not study.startswith('Hosny'):
        amin,amax=0.,float('inf');ok=True
        for r,p in zip(obs,raw):
            if r['power']:
                amin=max(amin,(r['low']/max(p,1e-30))**(1/r['power']));amax=min(amax,(r['high']/max(p,1e-30))**(1/r['power']))
            else:ok &= r['low']<=p<=r['high']
        feasible=bool(ok and amin<=amax)
    return res,amp,rows,feasible

def provenance():
    paths=[HERE/'campaign.py',ROOT/'wgd_langevin/data/observations.csv',ROOT/'wgd_langevin/hosny/data.py']
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def run(index):
    f,mask,study=jobs()[index];key=jobkey(f,mask,study);dest=OUT/(key+'.json')
    if dest.exists():
        record=json.loads(dest.read_text())
        if record['source_hashes']!=provenance():raise RuntimeError('Source changed: use a fresh result directory')
        return
    tick=time.monotonic();names,lo,hi=specification(f,mask,study);obs=observations(study)
    linear=f in ('gaussian','memory');n=512 if linear else 96
    rng=np.random.default_rng(20260926+index);starts=[]
    default=dict(tp=.6,tm=2.,r=.6 if f=='memory' else 2.,vs=.5,eta=.5,self_a=.1,self_K=.2,force=.1)
    base=np.log([default[k.replace('_poly','')] for k in names])
    seeds=[base]
    # Independent of all polarity calibrations, length laws and former WGD fits.
    for i in range(7):
        z=base+rng.normal(0,1.8,len(base))
        if i>=4:z=rng.uniform(lo,hi)
        seeds.append(np.clip(z,lo+1e-6,hi-1e-6))
    nullpath=OUT/(jobkey(f,(),study)+'.json')
    nullz=None
    if mask:
        if not nullpath.exists():raise RuntimeError('Fit shared null first')
        null=json.loads(nullpath.read_text())
        if null['source_hashes']!=provenance():raise RuntimeError('Incompatible null provenance')
        nullz=np.log([null['parameters'][k.replace('_poly','')] for k in names]);seeds[0]=nullz
    def objective(z,count=n):return score(predictions(z,names,f,study,n=count),obs,study)[0]
    for i,z in enumerate(seeds):
        z=np.clip(z,lo+1e-6,hi-1e-6)
        if i:
            for j,k in enumerate(names):
                if k.endswith('_poly'):z[j]=np.clip(z[j]+(-1 if i%2 else 1)*np.log(2),lo[j]+1e-6,hi[j]-1e-6)
        fit=least_squares(objective,z,bounds=(lo,hi),max_nfev=120,diff_step=.003 if not linear else .0001,ftol=1e-5,xtol=1e-5,gtol=1e-5)
        starts.append(dict(z=fit.x.tolist(),cost=float(fit.cost),success=bool(fit.success),nfev=fit.nfev,active_bounds=fit.active_mask.tolist()))
        print(key,i,fit.cost,fit.success,flush=True)
    polish_n=2048 if linear else 512
    candidates=sorted(starts,key=lambda r:r['cost'])[:2]
    if nullz is not None:candidates.append(dict(z=nullz.tolist(),cost=0))
    polished=[]
    for candidate in candidates:
        z=np.clip(candidate['z'],lo+1e-8,hi-1e-8)
        fit=least_squares(lambda zz:objective(zz,polish_n),z,bounds=(lo,hi),max_nfev=80,diff_step=.003 if not linear else .0001,ftol=1e-5,xtol=1e-5,gtol=1e-5)
        polished.append(fit)
    fit=min(polished,key=lambda x:x.cost);z=fit.x
    res,amp,rows,feasible=score(predictions(z,names,f,study,n=polish_n),obs,study)
    validations=[];validation_n=8192 if linear else 2048
    settings=[(937,.01,48.)] if linear else [(937,.01,48.),(937,.005,48.),(937,.005,96.)]
    for seed,dt,burn in settings:
        pred=predictions(z,names,f,study,n=validation_n,seed=seed,dt=dt,burn=burn)
        rr,_,vr,vf=score(pred,obs,study,amp)
        validations.append(dict(seed=seed,dt=dt,burn=burn,n=validation_n,score=float(rr@rr),rows=vr,interval_feasible_if_amplitude_reprofiled=vf))
    checks=observations(study,True);checkrows=[]
    if checks:
        lag=next((r['lag_h'] for r in checks if r.get('lag_h') is not None),None)
        pred=predictions(z,names,f,study,n=validation_n,seed=1937,dt=.005,burn=96.,lag_override=lag)
        checkrows=score(pred,checks,study,amp)[2]
    ratios={k[:-5]:float(np.exp(z[j]-z[names.index(k[:-5])]))for j,k in enumerate(names)if k.endswith('_poly')}
    # Local sensitivity after profiling amplitude. Rank is an identifiability diagnostic,
    # never a confidence interval, and can be affected by Monte Carlo derivatives.
    sv=np.linalg.svd(fit.jac,compute_uv=False);rank=int(np.sum(sv>max(sv[0]*1e-5,1e-7))) if len(sv) else 0
    nearbest=[r for r in starts if 2*r['cost']<=min(2*s['cost'] for s in starts)+1.]
    ranges={k[:-5]:[float(min(np.exp(np.array(r['z'])[j]-np.array(r['z'])[names.index(k[:-5])]) for r in nearbest)),
                      float(max(np.exp(np.array(r['z'])[j]-np.array(r['z'])[names.index(k[:-5])]) for r in nearbest))]
            for j,k in enumerate(names) if k.endswith('_poly')}
    near_bounds=[k for j,k in enumerate(names) if min(z[j]-lo[j],hi[j]-z[j])<.03*(hi[j]-lo[j])]
    diffs=[]
    for group in (0,1):
        theta,eta,force=decode(z,names,group);diffs.append(dict(tp=theta[0],tm=theta[1],r=theta[2],ps=theta[3],vs=theta[4],gain=theta[5],self_a=theta[6],self_K=theta[7],eta=eta,force=force,amplitude=amp))
    changed=[k for k in diffs[0] if diffs[0][k]!=diffs[1][k]]
    assert set(changed)<=({'tp','tm'} if mask else set())
    result=dict(key=key,family=f,indices=list(mask),study=study,names=names,z=z.tolist(),parameters=dict(zip(names,np.exp(z).tolist())),
        bounds=dict(zip(names,zip(np.exp(lo).tolist(),np.exp(hi).tolist()))),shared_amplitude=amp,ratios=ratios,score=float(res@res),rows=rows,
        all_reading_intervals=bool(all(r['within'] for r in rows)) if not study.startswith('Hosny') else None,
        amplitude_interval_feasible=feasible,starts=starts,polish_success=bool(fit.success),validations=validations,checks=checkrows,
        limited_information=study.startswith('Jemaa'),wall_seconds=time.monotonic()-tick,generated=datetime.now(timezone.utc).isoformat(),
        source_hashes=provenance(),slurm_job=os.environ.get('SLURM_JOB_ID'),
        identifiability=dict(profiled_jacobian_singular_values=sv.tolist(),local_rank=rank,shape_parameters=len(names),
            dimension_lower_bound=max(0,len(names)+1-len(obs)),nearbest_start_ratio_ranges=ranges,near_bounds=near_bounds,
            note='Local numerical rank and near-best starts are diagnostics, not confidence intervals.'),
        constraint_audit=dict(passes=True,changed_parameters=changed,group_parameters=diffs),
        simulation=dict(search_n=n,polish_n=polish_n,validation_n=validation_n))
    save_json(dest,result)
    from publish import publish
    publish()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int,nargs='?');p.add_argument('--manifest',action='store_true');a=p.parse_args()
    if a.manifest:
        import scipy,numba
        save_json(OUT/'manifest.json',dict(jobs=[dict(index=i,family=f,indices=list(mask),study=s,key=jobkey(f,mask,s)) for i,(f,mask,s) in enumerate(jobs())],
            generated=datetime.now(timezone.utc).isoformat(),source_hashes=provenance(),
            packages=dict(numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__),
            scope='All non-persistence parameters shared within each experiment; experiments fitted independently. No polarity fits used.'))
        print(len(jobs()))
    else:run(a.index)
