"""Round6 Figure3 deterministic response experiments. Frozen native round5 source.
No kinetics, noise law or physical protocol is fitted or retuned.
"""
from pathlib import Path
import os,sys,json,time,subprocess,argparse,shlex
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parents[1]
ROOT=next(p for p in OUT.parents if (p/'polarity1d/manuscript').is_dir())
CODE=ROOT/'polarity1d/manuscript/code/mechanisms';sys.path.insert(0,str(CODE))
import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.optimize._numdiff import approx_derivative
import figure3_round5_edge_20261023 as old
RUN=old.RUN;PYTHON=old.source.PYTHON
NATIVE=OUT/'supporting/native';META=OUT/'supporting/metadata'
MODELS=('wave_pinning','holmes_model3','spring_ring','goryachev')

def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [clean(x) for x in v]
 if isinstance(v,np.ndarray):return clean(v.tolist())
 if isinstance(v,(np.bool_,)):return bool(v)
 if isinstance(v,(np.integer,)):return int(v)
 if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
 return v

def write(p,v):old.write(p,clean(v))
def read(p):return old.read(p)
def wrap(d,C):return (np.asarray(d)+C/2)%C-C/2

def prepare():
 m=read(RUN/'manifest.json');jobs=[]
 for model in MODELS:
  basals=[c for c in m['jobs'] if c['model']==model and c['fraction']==.001]
  for size in ([1.,1.5,2.,2.5,3.,3.15,3.97] if model=='spring_ring' else [20.,30.,37.8,40.,47.6,50.,60.]):
   exact=next((c for c in basals if c['size']==size),None)
   nearest=exact or min(basals,key=lambda c:abs(c['size']-size))
   jobs.append(dict(slot=len(jobs),model=model,size=size,cells=round(size*(100 if model=='spring_ring' else 20)),platform=nearest['platform'],mean=nearest['mean'],source_slot=nearest['slot'],source_size=nearest['size'],input=nearest['input'],exact=exact is not None,edge_and_center_only=exact is None,source_input_sha256=old.sha(RUN/nearest['input'])))
 manifest=dict(generated_utc=old.now(),jobs=jobs,source_run=str(RUN),source_manifest_sha256=old.sha(RUN/'manifest.json'),source_code_sha256={str(p.relative_to(ROOT)):old.sha(p) for p in (Path(__file__),CODE/'figure3_round5_edge_20261023.py')},pulse_sd={'chemical':1.5,'spring':.15},pulse_duration={'chemical':50,'spring':5},fraction=.001,solver={'method':'BDF','rtol':2e-10,'atol':2e-12},preparation='Exact sizes reuse round5 settled f=.001 state; matched sizes periodically resample nearest original basal at unchanged native cell density then settle 16000 chemical/1600 Spring. Same residual, contrast and tail gates.',amplitude_stop='First >10% terminal shared-shift gain departure from nominal, topology change, or B peak; amplitude doubles from nominal/2. B peak itself included if next doubling would cross it.',duration='Fixed nominal dose; tau=12.5,25,50,100,200 chemical, Spring /10; observation ends at 20 nominal tau.',storage_estimate_bytes=28*3000000)
 write(META/'native_manifest.json',manifest)
 return manifest

def worker(slot):
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Native integrations require SLURM')
 m=read(META/'native_manifest.json');c=m['jobs'][slot];folder=NATIVE/f'{slot:03d}';folder.mkdir(parents=True,exist_ok=True)
 if (folder/'status.json').exists() and read(folder/'status.json').get('complete'):return
 start=time.monotonic();status=dict(condition=c,generated_utc=old.now(),complete=False,slurm_job_id=os.environ['SLURM_JOB_ID']);write(folder/'status.json',status)
 try:
  native,_=old.source.native_module(RUN,c['platform']);spring=c['platform']=='spring';n=c['cells'];size=c['size'];mean=c['mean'];ell=.15 if spring else 1.5;tau=5. if spring else 50.;horizon=20*tau;a0=old.REF[c['model']]*.0001
  if spring:
   closure=native.Closure(size,n);rhs=lambda y,a:closure.rhs(y,a);jac=lambda y:sparse.csc_matrix(closure.jacobian(y,mean));sparsity=None
  else:
   grid=native.ring.RingGrid(size,n);model=native.ring.native_model(c['model']);shape=(len(model.state_names),n);rhs=lambda y,a:native.ring.native_rhs(model,grid,y.reshape(shape),np.broadcast_to(a,(n,))).ravel();sparsity=native.ring.periodic_sparsity(shape[0],n);jac=lambda y:approx_derivative(lambda v:rhs(v,mean),y,method='3-point',sparsity=sparsity).tocsc()
  def solve(initial,a,end,times,maxstep):
   opts=dict(method='BDF',rtol=2e-10,atol=2e-12,max_step=maxstep,t_eval=times)
   if spring:opts['jac']=lambda t,y:sparse.csc_matrix(closure.jacobian(y,a))
   else:opts['jac_sparsity']=sparsity
   sol=solve_ivp(lambda t,y:rhs(y,a),(0.,end),initial,**opts)
   if not sol.success:raise RuntimeError(sol.message)
   return sol.y.T
  source_folder=RUN/'workers'/f"{c['source_slot']:03d}"
  if c['exact']:
   y=np.load(source_folder/'mean_state.npz')['state'];warm_tail=read(source_folder/'status.json')['mean_tail_scaled_change'];basal_C=read(RUN/'manifest.json')['jobs'][c['source_slot']]['basal_circumference']
  else:
   z=np.load(RUN/c['input'])['state'].ravel();old_n=round(c['source_size']*(100 if spring else 20));species=(len(z)-1)//old_n if spring else len(z)//old_n
   u=z[:species*old_n].reshape(species,old_n);old_x=(np.arange(old_n)+.5)/old_n;new_x=(np.arange(n)+.5)/n
   fields=np.array([np.interp(new_x,old_x,v,period=1) for v in u]);y=fields.ravel()
   if spring:y=np.r_[y,z[-1]*size/c['source_size']]
   basal_C=float(y[-1]) if spring else size
   end=1600 if spring else 16000;w=solve(y,mean,end,np.array([0,.9*end,end]),50 if spring else 500);y=w[-1];warm_tail=float(np.max(abs(w[-1]-w[-2]))/max(1,np.max(abs(y))))
  C=float(y[-1]) if spring else size;e0=old.edges(y[:n],C);residual=float(np.max(abs(rhs(y,mean))));contrast=float(np.ptp(y[:n]));status.update(mean_circumference=C,basal_circumference=basal_C,mean_stationary_rhs=residual,mean_tail_scaled_change=warm_tail,mean_contrast=contrast)
  np.savez_compressed(folder/'mean_state.npz',state=y)
  if e0 is None or contrast<.05 or residual>1e-6 or warm_tail>1e-4:
   status.update(complete=True,defined=False,reason='Fails unchanged unique-region, contrast, RHS or tail stationarity gate');return
  if c['exact']:
   z=np.load(source_folder/'pulse_responses.npz');g=z['phase_sensitivity'];ci=z['input_covariance_integral'];linear=read(source_folder/'status.json')['linear_response'];scale=np.mean([v['measured_edge_gain']/v['predicted_phase_gain'] for v in linear['edge_gain_validation']])
  else:
   g,linear=old.phase_kernel(rhs,jac,y,n,C,spring);ci,_=old.covariance_integral(n,basal_C if spring else C,ell,tau,read(RUN/'manifest.json')['native_protocols'][c['platform']]['noise']['knot']);scale=1.
  width=(e0[1]-e0[0])%C;center=(e0[0]+width/2)%C;x=(np.arange(n)+.5)*C/n
  standard=np.exp(-.5*(wrap(x[:,None]-x[None,:],C)/ell)**2);pulse_map=g@standard
  covg=np.fft.irfft(np.fft.rfft(ci)*np.fft.rfft(g*scale),n=n);contribution=(g*scale)*covg*mean**2;Ds=float(contribution.sum())
  np.savez_compressed(folder/'map.npz',x=x,center=center,edges=e0,activity=y[:n],g=g,pulse_map=pulse_map,input_covariance_integral=ci,contribution=contribution,Ds=Ds,scale=scale)
  nominal_times=np.unique(np.r_[np.linspace(0,tau,51),np.linspace(tau,horizon,381)])
  control_cache={};experiments=[];curves=[]
  def control(times):
   key=tuple(times)
   if key not in control_cache:control_cache[key]=solve(y,mean,horizon,times,tau/2)
   return control_cache[key]
  def measure(kind,label,position,amplitude,duration,multiplier=1):
   times=np.unique(np.r_[np.linspace(0,duration,51),np.linspace(duration,horizon,381)])
   base=control(times);d=wrap(x-position,C);profile=np.exp(-.5*(d/ell)**2)
   on_times=times[times<=duration];off_times=times[times>duration]-duration
   on=solve(y,mean+amplitude*profile,duration,on_times,min(tau,duration)/20)
   off=solve(on[-1],mean,horizon-duration,np.r_[0.,off_times],tau/2)[1:];states=np.vstack([on,off]);all_edges=[];centroids=[];widths=[];counts=[]
   for state in states:
    cc=float(state[-1]) if spring else C;u=state[:n];ee=old.edges(u,cc);counts.append(int(np.sum((u>=(u.min()+u.max())/2)&~np.roll(u>=(u.min()+u.max())/2,1))))
    all_edges.append(ee/cc*C if ee is not None else np.array([np.nan,np.nan]));moment=np.mean(u*np.exp(2j*np.pi*(np.arange(n)+.5)/n));centroids.append(np.angle(moment)*C/(2*np.pi))
    widths.append(((ee[1]-ee[0])%cc)/cc*C if ee is not None else np.nan)
   control_edges=np.array([old.edges(v[:n],float(v[-1]) if spring else C)/(float(v[-1]) if spring else C)*C for v in base]);delta=wrap(np.array(all_edges)-control_edges,C)
   control_cent=np.array([np.angle(np.mean(v[:n]*np.exp(2j*np.pi*(np.arange(n)+.5)/n)))*C/(2*np.pi) for v in base]);dc=wrap(np.array(centroids)-control_cent,C);dw=np.array(widths)-((control_edges[:,1]-control_edges[:,0])%C)
   gain=delta[-1]/(amplitude*duration);shared=float(np.mean(gain));pred=float(g@profile);topology=any(k!=1 for k in counts)
   r=dict(experiment=len(experiments),kind=kind,label=label,position=float(position%C),centered_position=float(wrap(position-center,C)),amplitude=amplitude,multiplier=multiplier,duration=duration,dose=amplitude*duration,left_gain=gain[0],right_gain=gain[1],shared_gain=shared,centroid_gain=dc[-1]/(amplitude*duration),terminal_width_change=dw[-1],predicted_gain=pred,topology_change=topology,maximum_region_count=max(counts),relative_map_error=abs(shared-pred)/max(abs(pred),np.max(abs(pulse_map))*1e-6),relaxation_tail_relative=float(np.nanmax(abs(delta[-1]-delta[-11]))/max(np.nanmax(abs(delta[-1])),1e-20)),integrated_shared_gain=np.trapezoid(np.mean(delta,axis=1),times)/(amplitude*duration))
   experiments.append(r)
   curves.extend(dict(experiment=r['experiment'],time=float(t),left_displacement=dd[0],right_displacement=dd[1],centroid_displacement=cc,width_change=ww,region_count=int(kk)) for t,dd,cc,ww,kk in zip(times,delta,dc,dw,counts))
   return r
  # Every condition has edge+centre; the two endpoint sizes get the full map test.
  positions=[('centre',center),('right_edge',e0[1])]
  endpoint=size in ([1.,3.] if spring else [20.,60.])
  if endpoint:
   positions += [('halfway_centre_edge',center+width/4),('edge_minus_3sd',e0[1]-3*ell),('edge_minus_1sd',e0[1]-ell),('edge_plus_1sd',e0[1]+ell),('edge_plus_3sd',e0[1]+3*ell),('basal_quarter',e0[1]+(C-width)/4),('antipode',center+C/2)]
  for label,position in positions:measure('map_validation',label,position,a0,tau)
  nominal=next(r for r in experiments if r['label']=='right_edge')['shared_gain']
  if c['exact']:
   amplitudes=[.5*a0]+[a0*2**k for k in range(0,15) if a0*2**k<old.REF[c['model']]]+[old.REF[c['model']]]
   previous=None
   for amp in amplitudes:
    r=measure('amplitude','right_edge',e0[1],amp,tau,amp/a0);departure=abs(r['shared_gain']/nominal-1) if np.isfinite(r['shared_gain']) else np.inf;r['gain_departure_from_nominal']=departure
    stop='topology_change' if r['topology_change'] else 'gain_departure_gt_10pct' if departure>.1 else 'continuous_cue_peak' if amp==old.REF[c['model']] else ''
    if stop:
     status['amplitude_limit']=dict(last_linear_amplitude=previous,first_stopping_amplitude=amp,stop_reason=stop);break
    previous=amp
   for duration in np.array([12.5,25,50,100,200])/(10 if spring else 1):measure('duration','right_edge',e0[1],a0*tau/duration,float(duration),tau/duration)
  import pandas as pd
  pd.DataFrame(clean(experiments)).to_csv(folder/'experiments.csv',index=False);pd.DataFrame(clean(curves)).to_csv(folder/'curves.csv.gz',index=False)
  status.update(complete=True,defined=True,n_experiments=len(experiments),Ds=Ds,native_gain_scale=scale,linear_response=linear)
 except Exception as exc:
  import traceback;status.update(error=str(exc),traceback=traceback.format_exc());raise
 finally:status.update(generated_utc=old.now(),wall_seconds=time.monotonic()-start);write(folder/'status.json',status)

def submit():
 manifest=read(META/'native_manifest.json');worker=OUT/'scripts/worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+shlex.join([PYTHON,'-B',str(Path(__file__)),'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
 cmd=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time=00:30:00','--job-name=f3round6','--array=0-'+str(len(manifest['jobs'])-1),'--output='+str(OUT/'logs/%A_%a.out'),str(worker)];p=subprocess.run(cmd,text=True,capture_output=True,check=True);write(META/'submission.json',dict(generated_utc=old.now(),command=cmd,stdout=p.stdout,stderr=p.stderr));print(p.stdout)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','submit','worker']);p.add_argument('--slot',type=int);a=p.parse_args();prepare() if a.action=='prepare' else submit() if a.action=='submit' else worker(a.slot)
