"""Reanalyse the original uniform-medium studies with the SAME Hosny Poisson model."""
from pathlib import Path
import sys,json,csv,time,hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import differential_evolution,brentq
ROOT=Path(__file__).resolve().parent;WGD=ROOT.parent;MODEL=WGD/'hosny/poisson_velocity'
sys.path.insert(0,str(MODEL));from model import simulate,metrics,msd
from tests import run_tests
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
TAU=np.geomspace(.005,48,27);RATE=np.geomspace(.05,64,29)
STUDIES={'Xuan2018':(1/6,12),'Mallin2023':(.5,24),'Jemaa2017':(.5,24),'Jemaa2023':(1,24)}
GROUPS=('reference','polyploid');KEYS=('speed','path_length','directness','endpoint')
OBS=list(csv.DictReader((WGD/'data/observations.csv').open()))
for r in OBS:
    for k in ('value','low','high','lag_h'):r[k]=float(r[k]) if r[k] else None
    if r['metric']=='slow_fraction':r['notes']='Slow fraction: one average cell length in 1 h. Numerical threshold and precise implementation unavailable; group-specific thresholds are not established.'

def write(name,rows):
    if not rows:return
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)

def bank_task(args):
    assay,i=args;cad,dur=assay;qq=np.zeros((len(RATE),4))
    for j,rate in enumerate(RATE):
        m=metrics(simulate(TAU[i],rate,1,0,n=2048,cadence=cad,duration=dur,seed=92314))
        qq[j]=[m['path'].mean()/dur,m['path'].mean(),m['directness'].mean(),m['endpoint'].mean()]
    return assay,i,qq

class Surface:
    def __init__(self,study,q):
        self.study=study;self.q=q;self.interp=RegularGridInterpolator((np.log(TAU),np.log(RATE)),q,bounds_error=True)
        self.obs={g:[r for r in OBS if r['study']==study and r['group']==g and r['role']=='fit'] for g in GROUPS}
    def units(self,z):
        tau,rate=np.exp(z);q=self.interp([z])[0];return dict(zip(KEYS,q)),tau,rate
    def fit_eval(self,groups,shapes,ampids,return_pars=False):
        info=[];num={};den={}
        for g,z,aid in zip(groups,shapes,ampids):
            vals,tau,rate=self.units(z);obs=self.obs[g];items=[]
            for r in obs:
                power=2 if r['metric']=='msd' else (0 if r['metric']=='directness' else 1)
                unit=float(msd(r['lag_h'],tau,rate,1,0)) if power==2 else vals[r['metric']]
                width=max(np.log(r['high']/r['low'])/2,.001);target=np.log(r['value']/max(unit,1e-300));weight=1/width**2
                num[aid]=num.get(aid,0)+weight*power*target;den[aid]=den.get(aid,0)+weight*power*power
                items.append((r,power,unit,width))
            info.append((g,tau,rate,aid,items))
        amplitudes={i:np.exp(num[i]/den[i]) for i in num};total=0.;pars=[]
        for g,tau,rate,aid,items in info:
            A=amplitudes[aid];score=sum((np.log(unit*A**power/r['value'])/width)**2 for r,power,unit,width in items);total+=score
            pars.append(dict(study=self.study,group=g,tau=tau,gamma=1/tau,rate=rate,amplitude=A,kappa=0.,rms_velocity=A*np.sqrt(rate*tau),score=score))
        return (total,pars) if return_pars else total
    def fit(self,g,fixed_rate=None):
        bounds=[(np.log(TAU[0]),np.log(TAU[-1])),(np.log(RATE[0]),np.log(RATE[-1]))]
        if fixed_rate is not None:bounds[1]=(np.log(fixed_rate),np.log(fixed_rate))
        rr=differential_evolution(lambda z:self.fit_eval([g],[z],[0]),bounds,seed=32,popsize=12,maxiter=180,tol=1e-8,polish=True)
        return self.fit_eval([g],[rr.x],[0],True)[1][0]

def compare(surface,mask):
    bounds=[(np.log(TAU[0]),np.log(TAU[-1])),(np.log(RATE[0]),np.log(RATE[-1]))];mapping=[0,1]
    for j in (0,1):
        if mask&(1<<j):mapping[j]=len(bounds);bounds.append(bounds[j])
    def unpack(z):return [z[:2],z[mapping]],[0,1 if mask&4 else 0]
    def obj(z):sh,aa=unpack(z);return surface.fit_eval(GROUPS,sh,aa)
    best=None
    for seed in (17,43):
        rr=differential_evolution(obj,bounds,seed=seed,popsize=14,maxiter=220,tol=1e-8,polish=True)
        if best is None or rr.fun<best.fun:best=rr
    sh,aa=unpack(best.x);score,pars=surface.fit_eval(GROUPS,sh,aa,True)
    return dict(study=surface.study,mask=mask,changes='+'.join(k for i,k in enumerate(('tau','rate','amplitude')) if mask&(1<<i)) or 'none',extra_parameters=mask.bit_count(),score=score,parameters=pars)

def main():
    start=time.time();write('simulator_checks.csv',run_tests());write('observations.csv',OBS)
    assays=set(STUDIES.values());banks={};tasks=[]
    for assay in assays:
        path=OUT/f'bank_{assay[0]:.6f}_{assay[1]}.npz'
        if path.exists():banks[assay]=np.load(path)['q']
        else:banks[assay]=np.empty((len(TAU),len(RATE),4));tasks += [(assay,i) for i in range(len(TAU))]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for count,(assay,i,qq) in enumerate(pool.map(bank_task,tasks)):
            banks[assay][i]=qq
            if count%9==0:print('Bank slice',count+1,'/',len(tasks),flush=True)
    for assay,q in banks.items():np.savez_compressed(OUT/f'bank_{assay[0]:.6f}_{assay[1]}.npz',q=q,tau=TAU,rate=RATE,cells=2048)
    fits=[];surfaces=[];rate_profiles=[];comp=[];pred=[];tracks={};timepred=[];slow=[];slopes=[];validation=[];samples={}
    for study,assay in STUDIES.items():
        s=Surface(study,banks[assay]);cad,dur=assay
        if study.startswith('Jemaa'):
            # Speed alone: conditional bands rather than spurious optimal tau/rate.
            for g in GROUPS:
                r=s.obs[g][0]
                for i,tau in enumerate(TAU):
                    for j,rate in enumerate(RATE):
                        u=s.q[i,j,0];surfaces.append(dict(study=study,group=g,tau=tau,rate=rate,amplitude=r['value']/u,amplitude_low=r['low']/u,amplitude_high=r['high']/u,score=0.,compatible=True))
            continue
        for g in GROUPS:
            f=s.fit(g);fits.append(f);print('FIT',f,flush=True)
            for i,tau in enumerate(TAU):
                for j,rate in enumerate(RATE):
                    sc,pp=s.fit_eval([g],[[np.log(tau),np.log(rate)]],[0],True);p=pp[0];v,_,_=s.units(np.log([tau,rate]));lo=0.;hi=np.inf;ok=True
                    for r in s.obs[g]:
                        k=r['metric'];power=2 if k=='msd' else(0 if k=='directness' else 1)
                        unit=msd(r['lag_h'],tau,rate,1,0) if power==2 else v[k]
                        if power==0:ok &= r['low']<=unit<=r['high']
                        else:lo=max(lo,(r['low']/unit)**(1/power));hi=min(hi,(r['high']/unit)**(1/power))
                    surfaces.append(dict(study=study,group=g,tau=tau,rate=rate,amplitude=p['amplitude'],amplitude_low=lo,amplitude_high=hi,score=sc,compatible=bool(ok and lo<=hi)))
            for rate in np.geomspace(RATE[0],RATE[-1],19):rate_profiles.append(s.fit(g,float(rate)))
            p=simulate(f['tau'],f['rate'],f['amplitude'],0,n=65536,cadence=cad,duration=dur,seed=77131);m=metrics(p);m['speed']=m['path']/dur;m['path_length']=m['path'];samples[study,g]=(p,m);tracks[study+'_'+g]=p[:16]
            for r in OBS:
                if r['study']!=study or r['group']!=g or r['metric'] not in ('speed','path_length','directness','endpoint','msd'):continue
                k=r['metric'];val=float(msd(r['lag_h'],f['tau'],f['rate'],f['amplitude'],0)) if k=='msd' else float(m[k].mean())
                pred.append(dict(id=r['id'],study=study,group=g,metric=k,role=r['role'],observed=r['value'],low=r['low'],high=r['high'],predicted=val,within_interval=bool(r['low']<=val<=r['high'])))
            vals,_,_=s.units(np.log([f['tau'],f['rate']]))
            for k in KEYS:validation.append(dict(study=study,group=g,metric=k,bank=vals[k]*(1 if k=='directness' else f['amplitude']),independent=float(m[k].mean()),simulation_mean_se=float(m[k].std()/np.sqrt(len(p)))))
            for t in np.geomspace(.05,36,140):timepred.append(dict(study=study,group=g,time=t,msd=float(msd(t,f['tau'],f['rate'],f['amplitude'],0))))
            if study=='Xuan2018':
                # Stationary overlapping one-hour windows: per-cell fraction of time, averaged across cells.
                steps=round(1/cad);net=np.linalg.norm(p[:,steps:]-p[:,:-steps],axis=2);path=np.stack([m0 for m0 in [np.linalg.norm(np.diff(p,axis=1),axis=2)]])[0]
                accumulated=np.cumsum(np.pad(path,((0,0),(1,0))),axis=1);hour_path=accumulated[:,steps:]-accumulated[:,:-steps]
                for mode,window in [('net_displacement',net),('sampled_path',hour_path)]:
                    flat=np.sort(window.ravel())
                    for L in np.linspace(0,60,241):slow.append(dict(group=g,interpretation=mode,threshold_um=L,slow_fraction=float(np.searchsorted(flat,L,side='right')/len(flat))))
                for low,high in [(1/6,6),(1,6),(1/6,12)]:
                    tt=np.geomspace(low,high,80);yy=msd(tt,f['tau'],f['rate'],f['amplitude'],0);slopes.append(dict(group=g,lag_min=low,lag_max=high,slope=float(np.polyfit(np.log(tt),np.log(yy),1)[0])))
        for mask in range(8):comp.append(compare(s,mask))
    write('fits.csv',fits);write('parameter_surfaces.csv',surfaces);write('rate_profiles.csv',rate_profiles);write('predictions.csv',pred);write('numerical_validation.csv',validation);write('time_predictions.csv',timepred);write('slow_fraction_sensitivity.csv',slow);write('slope_windows.csv',slopes);np.savez_compressed(OUT/'tracks.npz',**tracks)
    (OUT/'model_comparison.json').write_text(json.dumps(comp,indent=2));write('model_comparison.csv',[{k:v for k,v in r.items() if k!='parameters'} for r in comp])
    # A common-tau model has proportional MSD curves and cannot cross, irrespective of kick rates/amplitudes.
    xf={f['group']:f for f in fits if f['study']=='Xuan2018'}
    def curve(g,t):f=xf[g];return msd(t,f['tau'],f['rate'],f['amplitude'],0)
    ts=np.geomspace(1e-4,36,500);cross=[]
    for a,b in zip(ts[:-1],ts[1:]):
        if (curve('polyploid',a)-curve('reference',a))*(curve('polyploid',b)-curve('reference',b))<0:cross.append(brentq(lambda t:curve('polyploid',t)-curve('reference',t),a,b))
    (OUT/'checks.json').write_text(json.dumps(dict(xuan_crossovers_h=cross,xuan_reported_crossover_interval=[2,3],common_tau_MSD='Curves proportional; crossover impossible for any stationary rates and amplitudes at shared tau.',slow_fraction='Common unknown length threshold; net-displacement and sampled-path interpretations checked separately. Threshold calibration is not an independent fit or a switching-state model.'),indent=2))
    sources=[ROOT/'run.py',MODEL/'model.py',WGD/'data/observations.csv']+[WGD.parent/'literature_review'/f'{s}.md' for s in ('xuan_2018','mallin_2023','jemaa_2017','jemaa_2023')]
    (OUT/'metadata.json').write_text(json.dumps(dict(model_path=str(MODEL/'model.py'),model='Identical stationary compound-Poisson OU velocity model as Hosny; kappa=0 in uniform medium; homogeneous populations.',loss='Original mean-summary squared log residual / half-log-reading-width squared. Reading intervals are not SEM/CI; metric definitions retained.',bounds=dict(tau=TAU[[0,-1]].tolist(),rate=RATE[[0,-1]].tolist()),fit_targets='Xuan: sampled mean speed and one late MSD. Mallin: mean sampled path and mean directness. Jemaa: speed-only viable families.',scope='Original uniform-medium contrasts only; no added Mallin gradient/acrylamide datasets, no new directional fits, no cell heterogeneity.',wall_seconds=time.time()-start,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2));print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
