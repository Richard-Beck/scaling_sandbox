"""Joint WGD refits. Polarity transfers parameter eligibility, not constants/signs."""
import sys, json, argparse, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import cholesky, toeplitz
from numba import njit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from design import OUT as OLD, MODELS, FAMILIES, save_json, lengths
from literature import OBS, HOSNY, modes, displacement_variance, summaries
HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
STUDIES=['Xuan2018','Mallin2023','Jemaa2017','Jemaa2023','Hosny2026']

def selections():
    rows=[]
    for m in MODELS:
        for f in FAMILIES:
            s=json.loads((OLD/'scaling'/f'{m}_{f}.json').read_text())
            law=min((r for r in s['laws'] if r['indices']),key=lambda r:r['cost'])
            const=next(r for r in s['laws'] if not r['indices'])
            rows.append(dict(model=m,family=f,law=law['law'],indices=law['indices'],exponents=law['exponents'],
                improvement=1-law['cost']/max(const['cost'],1e-30),theta=law['theta'],
                strong_calibration_support=bool(law['cost']<.9*const['cost']),
                length_span=max(lengths(m))/min(lengths(m))))
    return rows

def jobs():
    keys=set()
    for r in selections():
        for study in STUDIES:
            mask=tuple(i for i in r['indices'] if study=='Hosny2026' or i!=5)
            keys.add((r['family'],mask,study));keys.add((r['family'],(),study))
    return sorted(keys)

def jobkey(f,mask,study):return f+'_'+('shared' if not mask else '-'.join(map(str,mask)))+'_'+study

# ps=1 fixes an exact latent-state gauge; readout scale remains free/shared.
# Parameter vector uses logs of positive physical parameters (hours and latent units).
def specification(f,mask,study):
    names=['tp'];lo=[.05];hi=[12.]
    if f!='gaussian':names+=['tm','r'];lo += [.05,.005];hi += [12.,.99 if f=='memory' else 8.]
    if f in ('coupled','self_activation'):names+=['vs','eta'];lo += [.03,.03];hi += [10.,5.]
    if f=='self_activation':names+=['self_a','self_K'];lo += [.0001,.02];hi += [2.,2.]
    if study=='Hosny2026':names+=['force'];lo += [.0001];hi += [10.]
    split={0:'tp',1:'tm',5:'force'}
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
def ensemble(theta,eta,force,cad,duration,lag,n,seed,dtmax=.02,burn=48.):
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
    if study!='Hosny2026':
        return [dict(group=int(r.group=='polyploid'),metric=r.metric,obs=float(r.value),low=float(r.low),high=float(r.high),
                power=2 if r.metric=='msd' else 0 if r.metric=='directness' else 1)
                for r in OBS[(OBS.study==study)&(OBS.role==('check' if checks else 'fit'))].itertuples()
                if r.metric in ('speed','msd','path_length','directness','endpoint') and (not checks or r.metric!='msd')]
    rows=[]
    for r in HOSNY:
        if r['id'] not in ('N','P','HN0','HP0','RN0','RP0'):continue
        for metric in ('path','directness','parallel','endpoint','perpendicular'):
            q=np.array(r[metric]);width=max(float(q[2]-q[0]),.04) if np.isfinite(q[[0,2]]).all() else .1
            for j in range(3):
                if np.isfinite(q[j]):rows.append(dict(group=int(r['group']=='PACC'),metric='path_length' if metric=='path' else metric,
                    quantile=j,obs=float(q[j]),width=width,power=int(metric in ('path','endpoint')),block=r['id']))
    return rows

def predictions(z,names,f,study,n=96,seed=431,dt=.02,burn=48.):
    cad,dur,lag=(1/6,12.,11.1111111111) if study=='Xuan2018' else (.5,16.,1.) if study=='Hosny2026' else (1. if study=='Jemaa2023' else .5,24.,1.)
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
        pred['quantiles']={k:np.quantile(v,[.25,.5,.75]).tolist() for k,v in s.items()}
        out.append(pred)
    return out

def score(pred,obs,study,amplitude=None):
    raw=np.array([pred[r['group']]['quantiles'][r['metric']][r['quantile']] if 'quantile' in r else pred[r['group']][r['metric']] for r in obs])
    power=np.array([r['power'] for r in obs]);target=np.array([r['obs'] for r in obs])
    if study=='Hosny2026':
        width=np.array([r['width'] for r in obs]);sel=power==1
        amp=max(1e-12,float(np.sum(raw[sel]*target[sel]/width[sel]**2)/np.sum((raw[sel]/width[sel])**2))) if amplitude is None else amplitude
        values=raw*amp**power;res=(values-target)/width
    else:
        width=np.array([np.log(r['high']/r['low'])/2 for r in obs])
        amp=float(np.exp(np.sum(power*np.log(target/np.maximum(raw,1e-30))/width**2)/np.sum(power**2/width**2))) if amplitude is None else amplitude
        values=raw*amp**power;res=np.log(np.maximum(values,1e-30)/target)/width
    rows=[dict(**r,predicted=float(v),within=bool(r.get('low',-np.inf)<=v<=r.get('high',np.inf)) if study!='Hosny2026' else None)for r,v in zip(obs,values)]
    feasible=None
    if study!='Hosny2026':
        amin,amax=0.,float('inf');ok=True
        for r,p in zip(obs,raw):
            if r['power']:
                amin=max(amin,(r['low']/max(p,1e-30))**(1/r['power']));amax=min(amax,(r['high']/max(p,1e-30))**(1/r['power']))
            else:ok &= r['low']<=p<=r['high']
        feasible=bool(ok and amin<=amax)
    return res,amp,rows,feasible

def run(index):
    f,mask,study=jobs()[index];key=jobkey(f,mask,study);dest=OUT/(key+'.json')
    if dest.exists():return
    tick=time.monotonic();names,lo,hi=specification(f,mask,study);obs=observations(study)
    linear=f in ('gaussian','memory');n=256 if linear else 64
    rng=np.random.default_rng(20260926+index);starts=[]
    default=dict(tp=.6,tm=2.,r=.6 if f=='memory' else 2.,vs=.5,eta=.5,self_a=.1,self_K=.2,force=.1)
    base=np.log([default[k.replace('_poly','')] for k in names])
    seeds=[base]
    # Include matched polarity-derived shapes after changing clock/state gauge.
    for source in selections():
        if source['family']!=f:continue
        if tuple(i for i in source['indices'] if study=='Hosny2026' or i!=5)!=mask:continue
        t=source['theta'];clock=.6/t[0];d=dict(default)
        d.update(tp=.6,tm=t[1]*clock,r=t[2],vs=t[4]/t[3],self_a=t[6]/clock,self_K=t[7]/t[3])
        seeds.append(np.log([max(d[k.replace('_poly','')],1e-10) for k in names]))
    while len(seeds)<4:seeds.append(base+rng.normal(0,1.3,len(base)))
    seeds=seeds[:4]
    nullpath=OUT/(jobkey(f,(),study)+'.json')
    if mask and nullpath.exists():
        null=json.loads(nullpath.read_text())['parameters']
        seeds[0]=np.log([null[k.replace('_poly','')] for k in names])
    def objective(z):return score(predictions(z,names,f,study,n=n),obs,study)[0]
    for i,z in enumerate(seeds):
        z=np.clip(z,lo+1e-6,hi-1e-6)
        for j,k in enumerate(names):
            if k.endswith('_poly') and not (i==0 and mask and nullpath.exists()):z[j]=np.clip(z[j]+(-1 if i%2 else 1)*np.log(2),lo[j]+1e-6,hi[j]-1e-6)
        fit=least_squares(objective,z,bounds=(lo,hi),max_nfev=65,diff_step=.003 if not linear else .0001,ftol=1e-4,xtol=1e-4,gtol=1e-4)
        starts.append(dict(z=fit.x.tolist(),cost=float(fit.cost),success=bool(fit.success),nfev=fit.nfev,active_bounds=fit.active_mask.tolist()))
        print(key,i,fit.cost,fit.success,flush=True)
    best=min(starts,key=lambda r:r['cost']);z=np.array(best['z'])
    # Final fit uses more cells than the search; independent validation is separate.
    polish_n=1024 if linear else 256
    fit=least_squares(lambda zz:score(predictions(zz,names,f,study,n=polish_n),obs,study)[0],z,bounds=(lo,hi),max_nfev=40,diff_step=.003,ftol=2e-4,xtol=2e-4)
    z=fit.x;res,amp,rows,feasible=score(predictions(z,names,f,study,n=polish_n),obs,study)
    validations=[]
    for seed,dt,burn in [(937,.02,48.),(137,.01,96.)]:
        pred=predictions(z,names,f,study,n=4096 if linear else 1024,seed=seed,dt=dt,burn=burn)
        rr,_,vr,vf=score(pred,obs,study,amp)
        validations.append(dict(seed=seed,dt=dt,burn=burn,score=float(rr@rr),rows=vr,interval_feasible_if_amplitude_reprofiled=vf))
    checks=observations(study,True) if study=='Mallin2023' else []
    checkrows=score(predictions(z,names,f,study,n=4096 if linear else 1024,seed=937),checks,study,amp)[2] if checks else []
    ratios={k[:-5]:float(np.exp(z[j]-z[names.index(k[:-5])]))for j,k in enumerate(names)if k.endswith('_poly')}
    save_json(dest,dict(key=key,family=f,indices=list(mask),study=study,names=names,z=z.tolist(),parameters=dict(zip(names,np.exp(z).tolist())),
        bounds=dict(zip(names,zip(np.exp(lo).tolist(),np.exp(hi).tolist()))),shared_amplitude=amp,ratios=ratios,score=float(res@res),rows=rows,
        all_reading_intervals=bool(all(r['within'] for r in rows)) if study!='Hosny2026' else None,
        amplitude_interval_feasible=feasible,starts=starts,polish_success=bool(fit.success),validations=validations,checks=checkrows,
        limited_information=study.startswith('Jemaa'),wall_seconds=time.monotonic()-tick,generated=datetime.now(timezone.utc).isoformat()))
    from publish import publish
    publish()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('index',type=int,nargs='?');p.add_argument('--manifest',action='store_true');a=p.parse_args()
    if a.manifest:
        save_json(OUT/'manifest.json',dict(jobs=[dict(index=i,family=f,indices=list(mask),study=s,key=jobkey(f,mask,s))for i,(f,mask,s) in enumerate(jobs())],selections=selections(),generated=datetime.now(timezone.utc).isoformat()))
        print(len(jobs()))
    else:run(a.index)
