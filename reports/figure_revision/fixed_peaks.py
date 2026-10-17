"""Chiou-inspired physical peak preservation at the assay restart, distinct from a mean-pool match."""
from pathlib import Path
import sys,json,argparse
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));from run_assays import sim,MODELS,NoStimulus,make_model,Grid1D,phenotype_summary,scalar_metrics
from polarity1d.protocols import ReversingGradient
LENGTHS=[15,18,21,27,33,39]
def tasks():
 return [(m,L,s,201) for m in MODELS for L in LENGTHS for s in [1,2,3]]+[(m,L,1,401) for m in MODELS for L in [15,39]]
def state_at(path,anchor):
 z=np.load(path);i=np.argmin(abs(z['time']-anchor));return z['x'],np.array([z['field__active'][i],z['field__inactive'][i]])
def extend(x,state,L,assay,n):
 # Keep spatial resolution at the reference value as separation increases.
 actual=int(round(L*n/15));grid=Grid1D(L,actual);xx=grid.x;u,v=state
 if assay.startswith('competing'):
  distance=np.minimum(xx,L-xx);left=xx<=L/2
  out=np.empty((2,actual));mid=len(x)//2;bg=state[:,mid]
  for j in range(2):
   out[j,left]=np.interp(distance[left],x[:mid+1],state[j,:mid+1],right=bg[j])
   out[j,~left]=np.interp(distance[~left],(15-x[mid:])[::-1],state[j,mid:][::-1],right=bg[j])
 else:
  out=np.array([np.interp(xx,x,a,right=a[-1]) for a in state]);bg=state[:,-1]
 # Cell-integrated active peak excess should match to numerical precision.
 target=(u-bg[0]).sum()*15/n;actual_excess=(out[0]-bg[0]).sum()*L/actual
 if actual_excess>0:out[0]=bg[0]+(out[0]-bg[0])*target/actual_excess
 return out,actual,bg

def run(model,L,seed,n):
 results=[]
 for assay,anchor in [('competing_fronts_extended_4x',200),('gradient_reversal_long_conditioning',2000)]:
  # Source trajectories are the new, complete native campaign with matched solver settings.
  refpath=HERE/f'data/pool_{model}_L15_s{seed}_n{n}_{assay}_concentration.npz'
  nativepath=HERE/f'data/pool_{model}_L{L}_s{seed}_n{n}_{assay}_concentration.npz'
  xr,sr=state_at(refpath,anchor);xn,sn=state_at(nativepath,anchor);extended,actual,bg=extend(xr,sr,L,assay,n)
  for arm,initial,cells in [('native_restart',sn,n),('fixed_peak',extended,actual)]:
   if assay.startswith('competing'):
    protocol=NoStimulus();metrics=dict(analysis_start=0,orientation_threshold=.05,polarization_threshold=.1,sustained_duration=50);end=20000-anchor
   else:
    protocol=ReversingGradient(.1/L,5,20);metrics=dict(reversal_time=5,analysis_start=0,orientation_threshold=.05,polarization_threshold=.1,sustained_duration=20);end=20005
   tr=sim(model,initial,L,cells,protocol,end,2.5 if assay.startswith('gradient') else 5)
   r=phenotype_summary(tr,**metrics);r.update(model=model,L=L,n=n,actual_cells=cells,seed=seed,arm=arm,assay=assay,anchor_s=anchor,**scalar_metrics(tr))
   r.update(initial_active_peak_mass=float((initial[0]-bg[0]).sum()*L/cells),reference_active_peak_mass=float((sr[0]-bg[0]).sum()*15/n),initial_inactive_mean=float(initial[1].mean()),reference_inactive_mean=float(sr[1].mean()),initial_peak_height=float(initial[0].max()),reference_peak_height=float(sr[0].max()),background_active=float(bg[0]),background_inactive=float(bg[1]))
   if assay.startswith('competing') and r.get('competition_resolution_latency_s') is not None:r['competition_latency_from_original_withdrawal_s']=r['competition_resolution_latency_s']+anchor-40
   tr.save(HERE/f'data/fixed_{model}_L{L}_s{seed}_n{n}_{assay}_{arm}.npz');results.append(r)
 (HERE/f'data/fixed_{model}_L{L}_s{seed}_n{n}.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--index',type=int);a=p.parse_args();run(*tasks()[a.index]);print('complete',tasks()[a.index],flush=True)
