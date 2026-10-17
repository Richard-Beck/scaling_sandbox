"""Stationary homogeneous Gaussian velocity OU fits to supplied literature summaries.
Noise coefficient eta: dv=-v/tau dt + eta dW. Simulator uses s=eta*sqrt(tau/2).
"""
from pathlib import Path
import sys,csv,json,time,hashlib
import numpy as np
from scipy.optimize import differential_evolution,brentq
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parent;WGD=ROOT.parent
sys.path.insert(0,str(WGD/'hosny/ou_velocity'))
from model import simulate,metrics,msd
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
TAU=np.geomspace(.0005,48,161);GROUPS=('reference','polyploid')
STUDIES={'Xuan2018':(1/6,12),'Mallin2023':(.5,24),'Jemaa2017':(.5,24),'Jemaa2023':(1,24)}
OBS=list(csv.DictReader((WGD/'data/observations.csv').open()))
for r in OBS:
 for k in ('value','low','high','lag_h'):r[k]=float(r[k]) if r[k] else None
 if r['metric']=='slow_fraction':r['notes']='One average cell length in one hour; numerical length and exact implementation unknown. Group-specific thresholds are not established.'
def write(name,rows):
 if not rows:return
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def sdev(tau,eta=1):return eta*np.sqrt(tau/2)
def moment(t,tau,eta=1):return msd(t,sdev(tau,eta),tau,0)
def raymean(t,tau):return np.sqrt(np.pi*moment(t,tau))/2

def checks():
 rows=[]
 for tau in (.0005,.1,3.):
  p,v=simulate(sdev(tau,4),tau,0,n=32768,duration=4,seed=456,return_velocity=True)
  for name,x,expected in [('MSD',(p[:,-1]**2).sum(1),float(moment(4,tau,4))),('velocity variance',(v[:,0]**2).sum(1),16*tau),('mean sampled speed',metrics(p)['path']/4,4*raymean(.5,tau)/.5)]:
   se=x.std()/np.sqrt(len(x));rows.append(dict(test=name,tau=tau,actual=float(x.mean()),expected=expected,passed=bool(abs(x.mean()-expected)<5*se)))
 write('simulator_checks.csv',rows)
 assert all(r['passed'] for r in rows)

def main():
 start=time.time();checks();write('observations.csv',OBS)
 # Only mean per-track directness needs simulation during fitting. Other means are exact.
 bank=OUT/'directness_bank.npz'
 if bank.exists():
  b=np.load(bank);assert np.array_equal(b['tau'],TAU);direct=b['directness']
 else:
  direct=[]
  for i,tau in enumerate(TAU):
   p=simulate(1,tau,0,n=16384,cadence=.5,duration=24,seed=92314)
   direct.append(metrics(p)['directness'].mean())
   if i%40==0:print('Directness bank',i+1,'/',len(TAU),flush=True)
  direct=np.array(direct);np.savez_compressed(bank,tau=TAU,directness=direct,cells=16384)
 interp=PchipInterpolator(np.log(TAU),direct)
 def unit(study,k,tau,lag=None):
  cad,dur=STUDIES[study]
  return {'speed':lambda:raymean(cad,tau)/cad,'path_length':lambda:dur*raymean(cad,tau)/cad,'endpoint':lambda:raymean(dur,tau),'directness':lambda:float(interp(np.log(tau))),'msd':lambda:float(moment(lag,tau))}[k]()
 def evaluate(study,groups,taus,ampids):
  items=[];num={};den={}
  for g,tau,aid in zip(groups,taus,ampids):
   for r in OBS:
    if r['study']!=study or r['group']!=g or r['role']!='fit':continue
    power=2 if r['metric']=='msd' else 0 if r['metric']=='directness' else 1
    u=unit(study,r['metric'],tau,r['lag_h']);width=np.log(r['high']/r['low'])/2
    num[aid]=num.get(aid,0)+power*np.log(r['value']/u)/width**2;den[aid]=den.get(aid,0)+power**2/width**2
    items.append((g,aid,r,power,u,width))
  amps={a:np.exp(num[a]/den[a]) for a in num};pars=[]
  for g,tau,aid in zip(groups,taus,ampids):
   eta=amps[aid];score=sum((np.log(u*eta**power/r['value'])/width)**2 for gg,a,r,power,u,width in items if gg==g)
   pars.append(dict(study=study,group=g,tau=float(tau),gamma=float(1/tau),noise_coefficient=float(eta),stationary_axis_sd=float(sdev(tau,eta)),rms_velocity=float(eta*np.sqrt(tau)),score=float(score),tau_at_bound=bool(np.isclose(tau,TAU[0]) or np.isclose(tau,TAU[-1]))))
  return sum(p['score'] for p in pars),pars
 def fit(study,groups,mask=3):
  nt=2 if len(groups)==2 and mask&1 else 1;aa=[0,1] if len(groups)==2 and mask&2 else [0]*len(groups)
  def ev(z):return evaluate(study,groups,np.exp(z) if nt==2 else [np.exp(z[0])]*len(groups),aa)
  r=differential_evolution(lambda z:ev(z)[0],[(np.log(TAU[0]),np.log(TAU[-1]))]*nt,seed=42,popsize=15,maxiter=300,tol=1e-10,polish=True)
  return ev(r.x)
 fits=[];comp=[];pred=[];surfaces=[];times=[];tracks={};slow=[];slopes=[];validation=[]
 for study,(cad,dur) in STUDIES.items():
  for g in GROUPS:
   if study.startswith('Jemaa'):
    r=next(r for r in OBS if r['study']==study and r['group']==g and r['role']=='fit')
    for tau in TAU:
     u=unit(study,'speed',tau);surfaces.append(dict(study=study,group=g,tau=tau,rms_low=r['low']/u*np.sqrt(tau),rms_high=r['high']/u*np.sqrt(tau),rms=r['value']/u*np.sqrt(tau),compatible=True))
    continue
   _,ps=fit(study,[g]);f=ps[0];fits.append(f);print('FIT',f,flush=True)
   tau=f['tau'];eta=f['noise_coefficient'];sd=f['stationary_axis_sd']
   p=simulate(sd,tau,0,n=65536,cadence=cad,duration=dur,seed=77131);m=metrics(p);m['speed']=m['path']/dur;m['path_length']=m['path'];tracks[study+'_'+g]=p[:16]
   for r in OBS:
    if r['study']!=study or r['group']!=g or r['metric'] not in ('speed','path_length','endpoint','directness','msd'):continue
    k=r['metric'];power=2 if k=='msd' else 0 if k=='directness' else 1
    val=float(m[k].mean()) if k=='directness' else float(unit(study,k,tau,r['lag_h'])*eta**power)
    pred.append(dict(id=r['id'],study=study,group=g,metric=k,role=r['role'],observed=r['value'],low=r['low'],high=r['high'],predicted=val,within_interval=bool(r['low']<=val<=r['high'])))
   for k in ('speed','path_length','directness','endpoint'):
    if k=='directness' and study!='Mallin2023':continue
    validation.append(dict(study=study,group=g,metric=k,fitting_prediction=float(unit(study,k,tau)*eta**(k!='directness')),independent=float(m[k].mean()),simulation_mean_se=float(m[k].std()/256)))
   for t in np.geomspace(.05,36,140):times.append(dict(study=study,group=g,time=t,msd=float(moment(t,tau,eta))))
   for tt in np.geomspace(TAU[0],TAU[-1],1601):
    lo=0.;hi=np.inf;ok=True
    for r in [r for r in OBS if r['study']==study and r['group']==g and r['role']=='fit']:
     power=2 if r['metric']=='msd' else 0 if r['metric']=='directness' else 1;u=unit(study,r['metric'],tt,r['lag_h'])
     if power==0:ok &= r['low']<=u<=r['high']
     else:lo=max(lo,(r['low']/u)**(1/power));hi=min(hi,(r['high']/u)**(1/power))
    surfaces.append(dict(study=study,group=g,tau=tt,rms_low=lo*np.sqrt(tt),rms_high=hi*np.sqrt(tt),rms=np.sqrt(lo*hi*tt),compatible=bool(ok and lo<=hi)))
   if study=='Xuan2018':
    steps=round(1/cad);net=np.linalg.norm(p[:,steps:]-p[:,:-steps],axis=2);acc=np.cumsum(np.pad(np.linalg.norm(np.diff(p,axis=1),axis=2),((0,0),(1,0))),axis=1)
    for mode,win in [('net_displacement',net),('sampled_path',acc[:,steps:]-acc[:,:-steps])]:
     flat=np.sort(win.ravel())
     for L in np.linspace(0,60,241):slow.append(dict(group=g,interpretation=mode,threshold_um=L,slow_fraction=float(np.searchsorted(flat,L,side='right')/len(flat))))
    for low,high in [(1/6,6),(1,6),(1/6,12)]:
     tt=np.geomspace(low,high,80);slopes.append(dict(group=g,lag_min=low,lag_max=high,slope=float(np.polyfit(np.log(tt),np.log(moment(tt,tau,eta)),1)[0])))
  if not study.startswith('Jemaa'):
   for mask in range(4):
    score,pars=fit(study,GROUPS,mask);comp.append(dict(study=study,mask=mask,changes=['none','tau','noise coefficient','both'][mask],score=score,parameters=pars))
 for name,rows in [('fits',fits),('predictions',pred),('parameter_surfaces',surfaces),('time_predictions',times),('slow_fraction_sensitivity',slow),('slope_windows',slopes),('numerical_validation',validation)]:write(name+'.csv',rows)
 write('model_comparison.csv',[{k:v for k,v in r.items() if k!='parameters'} for r in comp]);(OUT/'model_comparison.json').write_text(json.dumps(comp,indent=2));np.savez_compressed(OUT/'tracks.npz',**tracks)
 xf={f['group']:f for f in fits if f['study']=='Xuan2018'}
 def difference(t):
  return float(moment(t,xf['polyploid']['tau'],xf['polyploid']['noise_coefficient'])-moment(t,xf['reference']['tau'],xf['reference']['noise_coefficient']))
 ts=np.geomspace(1e-4,36,500);cross=[brentq(difference,a,b) for a,b in zip(ts[:-1],ts[1:]) if difference(a)*difference(b)<0]
 (OUT/'checks.json').write_text(json.dumps(dict(xuan_crossovers_h=cross,xuan_reported_crossover_interval=[2,3]),indent=2))
 sources=[Path(__file__),WGD/'hosny/ou_velocity/model.py',WGD/'data/observations.csv']
 (OUT/'metadata.json').write_text(json.dumps(dict(model='dv=-v/tau dt + eta dW; dx=v dt; stationary homogeneous isotropic Gaussian OU',parameterization='eta is Wiener coefficient; simulator sigma=eta*sqrt(tau/2); RMS speed=eta*sqrt(tau)',bounds_tau=TAU[[0,-1]].tolist(),fit='Original Xuan speed + late MSD; Mallin path + mean per-track directness. Reading intervals are not confidence intervals.',observation='Exact Gaussian moments for mean sampled speed/path, endpoint and MSD; simulated mean per-track directness. No ratio-of-means approximation.',hosny='Existing Gaussian simulator imported read-only; Hosny analyses and outputs unchanged.',wall_seconds=time.time()-start,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2));print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
