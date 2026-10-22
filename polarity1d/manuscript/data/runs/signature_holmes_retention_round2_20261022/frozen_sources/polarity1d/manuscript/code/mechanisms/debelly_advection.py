"""Native cue trajectories and exact instantaneous mechanical/MCA decomposition."""
from pathlib import Path
import sys,json,numpy as np
from debelly_controls import n,r,g
OUT=n.DATA/'runs/debelly_advection_20261008'
LENGTHS=list(np.arange(23.,40.,2.))
STAGES=['early','mid','zero']
def moment(b,L,x):return b@(2*x/L-1)/b.sum(axis=-1)
def boundary(rac,tension,p):
 R,H=n.ds._homogeneous_rac_rho(p,.5);threshold=1.3*R
 a=p.rac_to_rho_coupling*p.protrusion_speed*float(n.ds.smooth_threshold(rac[0],p.protrusion_switch_scale,threshold))/(1+tension[0]**2/p.protrusion_tension_scale)
 b=-p.rac_to_rho_coupling*p.protrusion_speed*float(n.ds.smooth_threshold(rac[-1],p.protrusion_switch_scale,threshold))/(1+tension[-1]**2/p.protrusion_tension_scale)
 return (a-b)/2,(b-a)/2

def adv(y,v,totalv,bc,dx):
 # Freeze native upwind donor selection on the full velocity: component fluxes then add exactly.
 faces=.5*(v[:-1]+v[1:]);totalfaces=.5*(totalv[:-1]+totalv[1:]);flux=faces*np.where(totalfaces>=0,y[:-1],y[1:]);res=np.zeros(len(y));res[:-1]-=flux/dx;res[1:]+=flux/dx
 res[0]+=bc[0]*y[0]/dx;res[-1]-=bc[1]*y[-1]/dx
 return res

def decompose(tr,i,L,cue):
 p=n.model_for('debelly').parameters;N=len(tr.x);dx=L/N;f=tr.fields;x=tr.x;b=f['mca_bound'][i];xi=p.linker_friction*b
 vm=(f['displacement'][i]-f['displacement'][i-1])/.5
 bc=boundary(f['rac'][i-1],f['tension'][i-1],p)
 force=p.active_stress*p.rho_to_rac_coupling*n.ds._gradient(f['rho'][i-1],dx)
 solve=n.ds._dirichlet_elliptic
 velocities={
  'boundary':solve(xi,p.cortical_viscosity,np.zeros(N),dx,*bc),
  'rho':solve(xi,p.cortical_viscosity,force,dx,0.,0.),
  'membrane':solve(xi,p.cortical_viscosity,xi*vm,dx,0.,0.)}
 v=f['velocity'][i];res=float(np.max(abs(sum(velocities.values())-v)));assert res<1e-11,(L,i,res)
 state={k:f[k][i].copy() for k in n.DB_NAMES};calls=[];original=n.ds._transport
 def traced(values,velocity,diffusion,step,dx,**kwargs):
  y=original(values,velocity,diffusion,step,dx,**kwargs)
  calls.append(dict(old=values.copy(),new=y.copy(),v=velocity.copy(),D=diffusion,dt=step,dx=dx,**kwargs));return y
 n.ds._transport=traced
 try:step=n.run('debelly',L,state,cue,.5,interval=.5)
 finally:n.ds._transport=original
 c=calls[0];new=step.state['mca_bound'];pre=c['new'];P=float(moment(b,L,x));q=2*x/L-1;w=(q-P)/new.sum();terms={}
 for name,vpart in velocities.items():terms[name]=-float(adv(pre,vpart,c['v'],c['boundary_velocity'] if name=='boundary' else (0.,0.),dx)@w)
 dif=n.ds._neumann_transport_matrix(N,dx,1.,c['D'],np.zeros(N));rhs=-(dif[1]-1)*pre;rhs[1:]-=dif[2,:-1]*pre[:-1];rhs[:-1]-=dif[0,1:]*pre[1:]
 terms['diffusion']=-float(rhs@w)
 exchange=-c['decay']*pre+c['source'];terms['exchange']=-float(exchange@w)
 closure=(new-pre)/.5;terms['replenishment']=-float(closure@w)
 actual=-float((moment(new,L,x)-P)/.5);bres=float(actual-sum(terms.values()));assert abs(bres)<1e-11,(L,i,bres)
 row=dict(length=L,time=float(tr.time[i]),mca_polarity=P,centroid_speed=L/2*actual,reorientation_rate=actual,mechanical_residual=res,budget_residual=bres,
  boundary_speed=float(abs(bc[1])),left_tension=float(f['tension'][i-1,0]),right_tension=float(f['tension'][i-1,-1]),left_rac=float(f['rac'][i-1,0]),right_rac=float(f['rac'][i-1,-1]),
  rho_difference=float(f['rho'][i-1,0]-f['rho'][i-1,-1]),stress_gradient_rms=float(np.sqrt(np.mean(force**2))),
  screening_length=float(np.sqrt(p.cortical_viscosity/(p.linker_friction*b.mean()))),
  rho_velocity_rms=float(np.sqrt(np.mean(velocities['rho']**2))),boundary_velocity_rms=float(np.sqrt(np.mean(velocities['boundary']**2))),
  **{'rate_'+k:value for k,value in terms.items()})
 profiles={'x_fraction':x/L,'rho':f['rho'][i-1],'mca':b,'stress_gradient':force,'total':v,**velocities}
 return row,profiles

def setup():
 OUT.mkdir(parents=True,exist_ok=True)
 manifest=r.read(n.DATA/'runs/three_test_reversal_20261006/manifest.json')
 cases=[c for c in manifest['cases'] if c['model']=='debelly' and c['condition']=='fixed_slope' and c['length'] in LENGTHS]
 assert len(cases)==len(LENGTHS)
 r.write(OUT/'cases.json',cases)
 r.write(OUT/'provenance.json',dict(generated_at=n.now(),sources={str(p):n.sha(p) for p in [Path(__file__),Path(n.ds.__file__),n.CODE/'native.py']},definition='Exact current relative20 basal checkpoints; continuous native fixed-slope gradient0.24, step0.5, 201 cells. Mechanical superposition at saved native states, not independently evolved interventions. Bound-MCA native-step budget uses frozen full-flow upwind donors and exact changing normalization.'))
 print(len(cases))
def worker(slot):
 for path,digest in r.read(OUT/'provenance.json')['sources'].items():assert n.sha(Path(path))==digest
 case=r.read(OUT/'cases.json')[slot];L=case['length'];cp=n.ROOT/case['checkpoint'];assert n.sha(cp)==case['checkpoint_sha256'];state,_=n.read_checkpoint(cp,'debelly');cue=g.response_cue(case)
 tr=n.run('debelly',L,state,cue,1200.,interval=.5);P=moment(tr.fields['mca_bound'],L,tr.x);cross=int(np.flatnonzero(P<=0)[0]);assert tr.time[cross]<1200
 indices={'early':int(np.searchsorted(tr.time,50.)),'mid':int(np.argmin(abs(P-.05))),'zero':cross}
 rows=[];profiles={}
 for stage,i in indices.items():
  row,prof=decompose(tr,i,L,cue);row['stage']=stage;rows.append(row);profiles.update({stage+'_'+k:v for k,v in prof.items()})
 np.savez_compressed(OUT/f'L{L:g}_profiles.npz',**profiles)
 n.save_checkpoint(OUT/f'L{L:g}_native.npz','debelly',tr.state,tr)
 r.write(OUT/f'L{L:g}_summary.json',dict(generated_at=n.now(),length=L,checkpoint_sha256=n.sha(cp),crossing_time=float(tr.time[cross]),crossing_bracket=[float(tr.time[cross-1]),float(tr.time[cross])],rows=rows,native_sha256=n.sha(OUT/f'L{L:g}_native.npz'),profiles_sha256=n.sha(OUT/f'L{L:g}_profiles.npz')))
 print(L,rows,flush=True)
if __name__=='__main__':
 if sys.argv[1]=='setup':setup()
 else:worker(int(sys.argv[2]))
