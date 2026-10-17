"""Figure-only assays: conserved pools, finite opposing cues, spring-model competition."""
from pathlib import Path
import sys,json,argparse,time,hashlib
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[2]; HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'polarity1d/src'));sys.path.insert(0,str(ROOT/'reports/trial_report_assets/zmurchok2020'))
from polarity1d.registry import make_model
from polarity1d.contract import Grid1D
from polarity1d.initial import noisy_initial_state
from polarity1d.simulation import simulate
from polarity1d.campaign import _make_protocol
from polarity1d.metrics import phenotype_summary
from polarity1d.protocols import NoStimulus,LinearGradient
MODELS=['wave_pinning','otsuji','goryachev']; LENGTHS=[9,12,15,18,21,27,33,39]
SOURCE=json.loads((ROOT/'polarity1d/campaigns/wp1_figure1_common_horizon_v1.json').read_text())

def config(model,assay,L,n,seed=1):
 t=next(t for t in SOURCE['tasks'] if t['model']==model and t['assay']==assay and t['simulation']['length']==15)
 sim=dict(t['simulation'],length=L,cells=n);protocol=dict(t['protocol'])
 if assay.startswith('gradient'):protocol['slope']=.1/L
 initial=noisy_initial_state(make_model(model),Grid1D(L,n),amplitude=.01 if assay.startswith('competing') else 0,seed=seed)
 return sim,protocol,initial,t['metrics']

def sim(model,initial,L,n,protocol,end,interval=5,max_step=10):
 return simulate(make_model(model),initial_state=initial,length=L,cells=n,protocol=protocol,end_time=end,output_interval=interval,max_step=max_step,backend='bdf')

def scalar_metrics(tr):
 u=tr.fields['active'];v=tr.fields['inactive'];L=tr.metadata['length_um'] if 'length_um' in tr.metadata else (tr.x[-1]+tr.x[0]);m=(u+v).mean(axis=1)*L
 return {'mass_error':float(np.max(abs(m/m[0]-1))),'initial_mass':float(m[0])}

def pool_task(model,L,seed=1,n=201):
 results=[]
 for assay,anchor in [('competing_fronts_extended_4x',200),('gradient_reversal_long_conditioning',1997.5)]:
  cfg,pc,initial,metrics=config(model,assay,L,n,seed)
  _,rpc,ri,_=config(model,assay,15,n,seed)
  ref=sim(model,ri,15,n,_make_protocol(rpc),anchor,2.5)
  uref=ref.fields['active'][-1];vref=ref.fields['inactive'][-1];target=float(vref.mean());refmass=float((ri.sum(axis=0)).mean()*15)
  # Chiou Figure 6: supply the extended cytoplasmic background, retaining peak excess.
  # Include the low active branch for models with nonzero basal active concentration.
  bulk=float(uref.min()+vref.max()); predicted=(refmass+bulk*(L-15))/L
  mean0=float(initial.sum(axis=0).mean())
  cache={}
  def probe(density):
   if density not in cache:
    tr=sim(model,initial*(density/mean0),L,n,_make_protocol(pc),anchor,2.5)
    cache[density]=float(tr.fields['inactive'][-1].mean())-target
   return cache[density]
  # Track the nearest root to the bulk-fill prediction; record ambiguity and mismatch.
  candidates=np.unique(np.r_[np.geomspace(max(.01,predicted*.3),max(mean0*1.5,predicted*2),12),mean0,predicted])
  brackets=[]
  for a,b in zip(candidates[:-1],candidates[1:]):
   fa,fb=probe(float(a)),probe(float(b))
   if fa*fb<=0:brackets.append((float(a),float(b)))
  root=None
  if L==15:root=mean0
  elif brackets:
   bracket=min(brackets,key=lambda b:abs(np.mean(b)-predicted));root=brentq(probe,*bracket,xtol=1e-6)
  arms=[('concentration',mean0),('bulk_fill',predicted)]+([] if root is None else [('matched_pool',root)])
  for arm,density in arms:
   tr=simulate(make_model(model),initial_state=initial*(density/mean0),protocol=_make_protocol(pc),**cfg)
   r=phenotype_summary(tr,**metrics);r.update(model=model,L=L,n=n,seed=seed,assay=assay,arm=arm,target_inactive_mean=target,anchor_s=anchor,root_count=len(brackets),initial_density=density,bulk_background_density=bulk,**scalar_metrics(tr))
   idx=int(np.argmin(abs(tr.time-anchor)));u=tr.fields['active'][idx];v=tr.fields['inactive'][idx]
   r.update(anchor_inactive_mean=float(v.mean()),inactive_mean_relative_error=float(abs(v.mean()/target-1)),anchor_active_mass=float(u.mean()*L),reference_active_mass=float(uref.mean()*15),anchor_peak_height=float(u.max()),anchor_peak_excess_mass=float((u-u.min()).mean()*L),reference_peak_excess_mass=float((uref-uref.min()).mean()*15))
   tr.save(HERE/f'data/pool_{model}_L{L}_s{seed}_n{n}_{assay}_{arm}.npz');results.append(r)
  if root is None:results.append(dict(model=model,L=L,n=n,seed=seed,assay=assay,arm='matched_pool',status='no_matching_root',target_inactive_mean=target,probes=cache))
 dest=HERE/f'data/pool_{model}_L{L}_s{seed}_n{n}.json';dest.write_text(json.dumps(results,indent=2))

class OpposingCue:
 def __init__(self,L,d):self.L=L;self.d=d
 @property
 def event_times(self):return (0.,float(self.d))
 def __call__(self,t,x,L):return .1*x/L if t<self.d else np.zeros_like(x)

def orientation(tr):
 u=tr.fields['active'];x=tr.x;return np.sum(u*(2*x/(x[-1]+x[0])-1),axis=1)/u.sum(axis=1)

def finite_task(model,L,n=201):
 m=make_model(model);initial=m.initial_state(Grid1D(L,n))
 # Establish polarity, remove the conditioning stimulus, and confirm its persistence.
 prep=sim(model,initial,L,n,LinearGradient(.1/L,'left'),2000,5)
 state=np.array([prep.fields[k][-1] for k in m.state_names])
 free=sim(model,state,L,n,NoStimulus(),2000,5);state=np.array([free.fields[k][-1] for k in m.state_names]);p=orientation(free)
 eligible=bool(p[-1]<-.05 and np.ptp(p[-101:])<.005)
 free.save(HERE/f'data/finite_{model}_L{L}_n{n}_conditioning.npz')
 rows=[];durations=[10,20,40,80,160,320,640,1280,2560,5120]
 for duration in durations if eligible else []:
  pulse=sim(model,state,L,n,OpposingCue(L,duration),duration,min(2.5,duration/10))
  after=np.array([pulse.fields[k][-1] for k in m.state_names]);recovery=sim(model,after,L,n,NoStimulus(),5000,5)
  pr=orientation(recovery);pp=orientation(pulse);settled=np.ptp(pr[-101:])<.005
  success=bool(pr[-1]>.05 and settled);mass=np.r_[(pulse.fields['active']+pulse.fields['inactive']).mean(axis=1),(recovery.fields['active']+recovery.fields['inactive']).mean(axis=1)]
  r=dict(model=model,L=L,n=n,duration=duration,eligible=eligible,initial_polarity=float(p[-1]),final_polarity=float(pr[-1]),settled=bool(settled),success=success,transient_max=float(pp.max()),mass_error=float(np.max(abs(mass/mass[0]-1))));rows.append(r)
  np.savez_compressed(HERE/f'data/finite_{model}_L{L}_n{n}_d{duration:g}.npz',time=np.r_[free.time-2000,pulse.time[1:],duration+recovery.time[1:]],polarity=np.r_[p,pp[1:],pr[1:]],x=free.x,active=np.vstack([free.fields['active'],pulse.fields['active'][1:],recovery.fields['active'][1:]]))
  if len(rows)>=10:
   yes=[r['duration'] for r in rows if r['success']]
   if yes:
    hi=min(yes);no=[r['duration'] for r in rows if not r['success'] and r['duration']<hi]
    if no and hi-max(no)>max(1,.01*hi):durations.append((hi+max(no))/2)
 (HERE/f'data/finite_{model}_L{L}_n{n}.json').write_text(json.dumps(dict(model=model,L=L,n=n,eligible=eligible,initial_polarity=float(p[-1]),conditioning_range=float(np.ptp(p[-101:])),assays=rows),indent=2))

def spring_task(L,n=201):
 import run as spring
 from scipy.integrate import solve_ivp
 def integrate(y,end,cue=False):
  def rhs(t,y):
   out=spring.rhs(t,y,L0=L)
   if cue and t<40:
    x=np.linspace(0,1,n);activation=4*(1.01*np.exp(-(x/.2)**2)+.99*np.exp(-((x-1)/.2)**2));u=y[:-2:2];v=y[1:-2:2];length=y[-1]-y[-2];delta=3+80*(length-L)
    transfer=spring.auth['f'](u,v,4+activation,5,6,2,delta)-spring.auth['f'](u,v,4,5,6,2,delta)
    out[:-2:2]+=transfer;out[1:-2:2]-=transfer
   return out
  times=np.arange(0,end+1,1.)
  s=solve_ivp(rhs,(0,end),y,method='LSODA',t_eval=times,rtol=2e-7,atol=2e-9,max_step=2)
  if not s.success:raise RuntimeError(s.message)
  return s.t,s.y.T
 # Start with two physical peaks, then use finite dual cues; a controlled imbalance selects a winner.
 y=spring.initial(n,L,0);x=np.linspace(0,1,n);u=.05645+.5*(1.01*np.exp(-(x/.1)**2)+.99*np.exp(-((x-1)/.1)**2));y[:-2:2]=u;y[1:-2:2]=2-u
 t1,y1=integrate(y,40,True);t2,y2=integrate(y1[-1],3000)
 t=np.r_[t1,t2[1:]+40];y=np.vstack([y1,y2[1:]]);p,m,length=spring.metrics(y);u=y[:,:-2:2];left=u[:,:n//3].max(axis=1);right=u[:,-n//3:].max(axis=1);base=u.min(axis=1);loser=np.minimum(left,right)-base;winner=np.maximum(left,right)-base
 two=(loser>.1*winner)&(winner>.1);est=bool(np.any(two[t<=40]));resolved=(loser<.1*winner)&(winner>.1)&(t>=40)
 persistent=np.logical_and.accumulate(resolved[::-1])[::-1];ix=np.flatnonzero(persistent&(t<=t[-1]-100));latency=None if not est or not len(ix) else float(t[ix[0]]-40)
 np.savez_compressed(HERE/f'data/spring_L{L:g}_n{n}.npz',time=t,y=y,polarity=p,length=length,stimulus=np.where(t[:,None]<40,4*(1.01*np.exp(-(x/.2)**2)+.99*np.exp(-((x-1)/.2)**2)),0),protrusion_left=spring.auth['F'](u[:,0]),protrusion_right=spring.auth['F'](u[:,-1]),retraction=.01*(length-L),tension=80*(length-L),loser_ratio=loser/np.maximum(winner,1e-12))
 (HERE/f'data/spring_L{L:g}_n{n}.json').write_text(json.dumps(dict(L=L,n=n,established=est,resolution_latency=latency,censored=bool(est and latency is None),final_polarity=float(p[-1]),final_contrast=float(np.ptp(u[-1])),mass_error=float(np.max(abs(m/m[0]-1))),time_unit='authors model time'),indent=2))

def tasks():
 out=[]
 for model in MODELS:
  for L in LENGTHS:
   for seed in [1,2,3]:out.append(('pool',model,L,seed,201))
   out.append(('finite',model,L,201))
 for L in [.75,1,1.25,1.5,2,2.5,3]:
  for n in [101,201]:out.append(('spring',L,n))
 # Matched spatial refinement at both endpoints for representative competition/reversal/finite assays.
 for model in MODELS:
  for L in [15,39]:out.extend([('pool',model,L,1,401),('finite',model,L,401)])
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--index',type=int);p.add_argument('--manifest',action='store_true');a=p.parse_args()
 if a.manifest:(HERE/'manifest.json').write_text(json.dumps(tasks(),indent=2));print(len(tasks()))
 else:
  task=tasks()[a.index];start=time.time();print(task,flush=True)
  {'pool':pool_task,'finite':finite_task,'spring':spring_task}[task[0]](*task[1:]);print('done',time.time()-start,flush=True)
