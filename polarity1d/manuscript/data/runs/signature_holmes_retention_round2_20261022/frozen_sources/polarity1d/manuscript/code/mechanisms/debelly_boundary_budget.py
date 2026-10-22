"""Correctly labeled instantaneous transport attribution for both boundary variants."""
import numpy as np
from debelly_boundary_resolution import OUT,n,r,model,run,lambda_cue,polarity
from debelly_advection import adv

def main():
 rows=[]
 for job in r.read(OUT/'jobs.json'):
  L=job['length'];variant=job['variant'];p=model(variant).parameters;z=np.load(OUT/f'{variant}_L{L:g}_cue.npz');summary=r.read(OUT/f'{variant}_L{L:g}_summary.json');x=z['x'];q=2*x/L-1;dx=L/len(x)
  for record in summary['rows']:
   stage=record['stage'];i=int(np.searchsorted(z['time'],record['time']));b=z['field_mca_bound'][i];xi=p.linker_friction*b;vm=(z['field_displacement'][i]-z['field_displacement'][i-1])/.5;force=p.active_stress*n.ds._gradient(z['field_rho'][i-1],dx);bc=(record['used_left_velocity'],record['used_right_velocity']);solve=n.ds._dirichlet_elliptic
   velocities={'prescribed_endpoints':solve(xi,p.cortical_viscosity,np.zeros(len(x)),dx,*bc),'direct_rho_stress':solve(xi,p.cortical_viscosity,force,dx,0,0),'membrane_drag':solve(xi,p.cortical_viscosity,xi*vm,dx,0,0)}
   actualv=z['field_velocity'][i];assert np.max(abs(sum(velocities.values())-actualv))<1e-11
   state={k:z['field_'+k][i].copy() for k in n.DB_NAMES};calls=[];original=n.ds._transport
   def traced(values,velocity,diffusion,step,dx,**kwargs):
    y=original(values,velocity,diffusion,step,dx,**kwargs);calls.append(dict(new=y.copy(),v=velocity.copy(),D=diffusion,dt=step,dx=dx,**kwargs));return y
   n.ds._transport=traced
   try:tr=run(variant,L,state,lambda_cue(),.5,interval=.5)
   finally:n.ds._transport=original
   c=calls[0];new=tr.state['mca_bound'];pre=c['new'];P=polarity(b,x,L);w=(q-P)/new.sum();terms={}
   for name,v in velocities.items():terms[name]=-float(adv(pre,v,c['v'],c['boundary_velocity'] if name=='prescribed_endpoints' else (0,0),dx)@w)
   A=n.ds._neumann_transport_matrix(len(x),dx,1,c['D'],np.zeros(len(x)));d=-(A[1]-1)*pre;d[1:]-=A[2,:-1]*pre[:-1];d[:-1]-=A[0,1:]*pre[1:];terms['bound_diffusion']=-float(d@w);terms['exchange']=-float((-c['decay']*pre+c['source'])@w);terms['uniform_replenishment']=-float(((new-pre)/.5)@w)
   actual=-float((polarity(new,x,L)-P)/.5);assert abs(sum(terms.values())-actual)<1e-11;np.testing.assert_allclose(new,z['field_mca_bound'][i+1],atol=1e-12,rtol=1e-10)
   rows.append(dict(variant=variant,length=L,stage=stage,time=record['time'],**terms,total=actual,residual=actual-sum(terms.values())))
 n.csvwrite(OUT/'instantaneous_budget.csv',rows);r.write(OUT/'budget_audit.json',dict(generated_at=n.now(),passed=True,states=len(rows),max_residual=max(abs(a['residual']) for a in rows),agrees_with_stored_next_native_step=True,labels='Direct source decomposition with frozen state; not causal removal or independently advected components. Total-flow donors retained.'))
 for a in rows:
  if a['variant']=='independent_endpoints' and a['length'] in [23.,39.]:print(a)
if __name__=='__main__':main()
