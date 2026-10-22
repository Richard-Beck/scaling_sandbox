"""Independent array/force-balance audit, without using decomposition routines."""
from pathlib import Path
import numpy as np
from debelly_corrected_advection import OUT,n,r

def main():
 cases=r.read(OUT/'cases.json');prov=r.read(OUT/'provenance.json')
 for path,digest in prov['sources'].items():assert n.sha(Path(path))==digest
 p=n.model_for('debelly').parameters;maxres=0.;maxdrift=0.;maxrate=0.;switchmin=1.;archives=[]
 R,H=n.ds._homogeneous_rac_rho(p,.5)
 for case in cases:
  L=case['length'];cp=n.ROOT/case['checkpoint'];assert n.sha(cp)==case['checkpoint_sha256'];state,_=n.read_checkpoint(cp,'debelly')
  summary=r.read(OUT/f'L{L:g}_summary.json');file=OUT/f'L{L:g}_native.npz';assert n.sha(file)==summary['native_sha256'];z=np.load(file);prof=np.load(OUT/f'L{L:g}_profiles.npz');x=z['x'];N=len(x);dx=L/N;t=z['time'];q=2*x/L-1
  assert np.all(np.diff(t)==.5) and t[-1]==1200 and N==201
  for k in n.DB_NAMES:
   np.testing.assert_array_equal(z['field_'+k][0],state[k]);assert np.all(np.isfinite(z['field_'+k]))
  for k in ['rac','rho','mca_bound','mca_unbound']:assert np.min(z['field_'+k])>=0
  b=z['field_mca_bound'];u=z['field_mca_unbound'];mass=(b+u).sum(1);maxdrift=max(maxdrift,float(np.max(abs(mass/mass[0]-1))));assert maxdrift<1e-10
  P=np.sum(b*q,axis=1)/np.sum(b,axis=1);cross=np.flatnonzero(P<=0)[0];assert t[cross]==summary['crossing_time']
  for a in summary['rows']:
   stage=a['stage'];i=int(np.searchsorted(t,a['time']));xi=p.linker_friction*b[i];vm=(z['field_displacement'][i]-z['field_displacement'][i-1])/.5
   assert abs(P[i]-a['mca_polarity'])<1e-12
   rate=-(P[i+1]-P[i])/.5;maxrate=max(maxrate,abs(rate-a['reorientation_rate']));assert maxrate<1e-11
   assert abs(sum(a['rate_'+k] for k in ['boundary','rho','membrane','diffusion','exchange','replenishment'])-rate)<1e-11
   total=sum(prof[stage+'_'+k] for k in ['boundary','rho','membrane']);np.testing.assert_allclose(total,z['field_velocity'][i],atol=1e-12,rtol=1e-10)
   rac=z['field_rac'][i-1];ten=z['field_tension'][i-1];switch=n.ds.smooth_threshold(rac[[0,-1]],p.protrusion_switch_scale,1.3*R);switchmin=min(switchmin,float(switch.min()))
   raw=np.array([p.protrusion_speed,-p.protrusion_speed])*switch/(1+ten[[0,-1]]**2/p.protrusion_tension_scale);bc=raw
   rho=z['field_rho'][i-1];grad=np.empty(N);grad[1:-1]=(rho[2:]-rho[:-2])/(2*dx);grad[0]=(rho[1]-rho[0])/dx;grad[-1]=(rho[-1]-rho[-2])/dx;force=p.active_stress*p.rho_to_rac_coupling*grad
   for component in ['boundary','rho','membrane']:
    v=prof[stage+'_'+component];ratio=p.cortical_viscosity/dx**2;lhs=(xi+2*ratio)*v;lhs[[0,-1]]+=ratio*v[[0,-1]];lhs[1:]-=ratio*v[:-1];lhs[:-1]-=ratio*v[1:]
    rhs=np.zeros(N)
    if component=='rho':rhs=force
    elif component=='membrane':rhs=xi*vm
    else:rhs[[0,-1]]=2*ratio*bc
    residual=float(np.linalg.norm(lhs-rhs)/(1+np.linalg.norm(rhs)));maxres=max(maxres,residual);assert residual<1e-9,(L,stage,component,residual)
  archives.append(dict(length=L,sha256=n.sha(file),checkpoint_sha256=n.sha(cp)))
 audit=dict(generated_at=n.now(),passed=True,completed_lengths=len(cases),completed_states=3*len(cases),exact_native_initial_states=True,positive_finite_fields=True,max_relative_total_mca_drift=maxdrift,max_independent_force_balance_residual=maxres,max_archived_next_step_rate_error=maxrate,minimum_endpoint_rac_switch=switchmin,scientific_sources_unchanged=True,archives=archives)
 r.write(OUT/'independent_audit.json',audit);print(audit)
if __name__=='__main__':main()
