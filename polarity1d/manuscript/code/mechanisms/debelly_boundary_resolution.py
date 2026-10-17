"""Resolve centering as a boundary-condition change; audit independent-end variant.

The scientific solver is unchanged. The corrected exploratory adapter disables
in-loop projection and retains each end's own published Rac/tension command.
An inactive end is explicitly constrained; this is not a free migrating cell.
"""
from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace
import json,sys,numpy as np
from debelly_controls import n,r,base
OUT=n.DATA/'runs/debelly_boundary_resolution_20261008'
LENGTHS=[23.,27.,31.,35.,39.]
VARIANTS=['legacy_centered','independent_endpoints']
def model(variant):
 return n.ds.DeBelly(replace(n.model_for('debelly').parameters,center_boundary_motion=variant=='legacy_centered'))
def run(variant,L,state,cue,end,interval=5):
 tr=n.ds.simulate_debelly(model(variant),protocol=cue,length=L,cells=201,end_time=end,output_interval=interval,step=.5,initial_checkpoint=state)
 fields=tr.fields;final={k:fields[k][-1].copy() for k in n.DB_NAMES}
 return SimpleNamespace(time=tr.time,x=tr.x,fields=fields,state=final)
def polarity(y,x,L):return np.sum(y*(2*x/L-1),axis=-1)/np.sum(y,axis=-1)
def save(path,tr):np.savez_compressed(path,time=tr.time,x=tr.x,**{'field_'+k:v for k,v in tr.fields.items()},**{'state_'+k:v for k,v in tr.state.items()})
def snapshot(variant,tr,i,L):
 p=model(variant).parameters;f=tr.fields;dx=L/201;b=f['mca_bound'][i];xi=p.linker_friction*b;R,H=n.ds._homogeneous_rac_rho(p,.5);rac=f['rac'][i-1];ten=f['tension'][i-1]
 raw=np.array([p.protrusion_speed,-p.protrusion_speed])*n.ds.smooth_threshold(rac[[0,-1]],p.protrusion_switch_scale,1.3*R)/(1+ten[[0,-1]]**2/p.protrusion_tension_scale)
 bc=raw-raw.mean() if variant=='legacy_centered' else raw
 vm=(f['displacement'][i]-f['displacement'][i-1])/.5;force=p.active_stress*n.ds._gradient(f['rho'][i-1],dx);solve=n.ds._dirichlet_elliptic
 parts={'endpoint':solve(xi,p.cortical_viscosity,np.zeros(201),dx,*bc),'direct_rho':solve(xi,p.cortical_viscosity,force,dx,0,0),'membrane_drag':solve(xi,p.cortical_viscosity,xi*vm,dx,0,0)}
 err=float(np.max(abs(sum(parts.values())-f['velocity'][i])));assert err<1e-11
 P=polarity(f['mca_bound'],tr.x,L);rate=-float((P[i+1]-P[i])/.5)
 return dict(time=float(tr.time[i]),mca_polarity=float(P[i]),centroid_speed=L/2*rate,reorientation_rate=rate,raw_left_velocity=float(raw[0]),raw_right_velocity=float(raw[1]),used_left_velocity=float(bc[0]),used_right_velocity=float(bc[1]),left_tension=float(ten[0]),right_tension=float(ten[-1]),direct_rho_rms=float(np.sqrt(np.mean(parts['direct_rho']**2))),endpoint_rms=float(np.sqrt(np.mean(parts['endpoint']**2))),membrane_drag_rms=float(np.sqrt(np.mean(parts['membrane_drag']**2))),reconstruction_error=err)

def setup():
 OUT.mkdir(parents=True,exist_ok=True);jobs=[dict(length=L,variant=v) for L in LENGTHS for v in VARIANTS]
 r.write(OUT/'jobs.json',jobs)
 r.write(OUT/'provenance.json',dict(generated_at=n.now(),sources={str(p):n.sha(p) for p in [Path(__file__),Path(n.ds.__file__),Path(base.__code__.co_filename)]},preparation='Independent homogeneous basal rest plus identical -1e-5 cosine bound-MCA seed; 60000 basal units. Not substituted into the retained manuscript assay.',challenge='Continuous gradient 0.24*x for 1500 s from each variant own prepared state; no old checkpoint is silently reinterpreted.',correction='center_boundary_motion=False; independent local prescribed end commands, inactive end constrained. No substrate or free-cell claim.'))
 print(len(jobs))
def worker(slot):
 for path,digest in r.read(OUT/'provenance.json')['sources'].items():assert n.sha(Path(path))==digest
 job=r.read(OUT/'jobs.json')[slot];L=job['length'];v=job['variant'];state=base(L);state['mca_bound']-=1e-5*np.cos(np.pi*n.rd.Grid1D(L,201).x/L)
 prep=run(v,L,state,n.ConstantInput(0),60000.,interval=20.);save(OUT/f'{v}_L{L:g}_preparation.npz',prep)
 # Last 10000 is an observational holdout: no convergence time is inferred.
 mask=prep.time>=50000.;tailchange={k:float(np.max(np.ptp(prep.fields[k][mask],axis=0))) for k in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension']}
 tr=run(v,L,prep.state,lambda_cue(),1500.,interval=.5);save(OUT/f'{v}_L{L:g}_cue.npz',tr);P=polarity(tr.fields['mca_bound'],tr.x,L);cross=np.flatnonzero(P<=0)
 indices={'early':100,'mid':int(np.argmin(abs(P-.05))),'zero':int(cross[0]) if len(cross) else 2998};rows=[]
 for stage,i in indices.items():
  i=min(max(1,i),len(tr.time)-2);row=snapshot(v,tr,i,L);row['stage']=stage;rows.append(row)
 finalP=float(polarity(prep.fields['rac'][-1],prep.x,L));savepath=OUT/f'{v}_L{L:g}_summary.json'
 r.write(savepath,dict(generated_at=n.now(),variant=v,length=L,preparation_left_polarity=finalP,preparation_mca_polarity=float(P[0]),preparation_tail_change=tailchange,prepared_polarized=bool(finalP<-.05 and np.ptp(prep.fields['rac'][-1])>=.05),crossing_time=float(tr.time[cross[0]]) if len(cross) else None,rows=rows,preparation_sha256=n.sha(OUT/f'{v}_L{L:g}_preparation.npz'),cue_sha256=n.sha(OUT/f'{v}_L{L:g}_cue.npz')))
 print(v,L,finalP,rows,flush=True)
class lambda_cue:
 event_times=()
 def __call__(self,time,x,length):return .24*np.asarray(x)
if __name__=='__main__':
 if sys.argv[1]=='setup':setup()
 else:worker(int(sys.argv[2]))
