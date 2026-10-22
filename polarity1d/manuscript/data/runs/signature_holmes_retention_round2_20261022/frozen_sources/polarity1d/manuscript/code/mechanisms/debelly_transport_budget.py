"""Exact native half-step MCA moment budget at common pulse times."""
from pathlib import Path
import sys,numpy as np
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE));import gradient_reversal as g
r,n,e=g.r,g.n,g.e;OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly'
rows=[]
def op(y,v,D,bound,dx):
 A=n.ds._neumann_transport_matrix(len(y),dx,1.,D,v,bound);ay=(A[1]-1)*y;ay[1:]+=A[2,:-1]*y[:-1];ay[:-1]+=A[0,1:]*y[1:];return -ay
for L in [18,23,39]:
 job=next(j for j in r.read(OUT/'jobs.json') if j['case']['length']==L and j['label']=='pass');z=np.load(OUT/f'L{L}_pass.npz');x=z['x'];q=2*x/L-1;dx=L/201
 for tt in [25.,50.,100.,150.]:
  i=int(np.argmin(abs(z['pulse_time']-tt)));state={k:z['pulse_field_'+k][i].copy() for k in n.DB_NAMES};calls=[];original=n.ds._transport
  def traced(values,velocity,diffusion,step,dx,*,decay=0.,source=0.,boundary_velocity=(0.,0.)):
   result=original(values,velocity,diffusion,step,dx,decay=decay,source=source,boundary_velocity=boundary_velocity)
   calls.append(dict(old=values.copy(),new=result.copy(),velocity=velocity.copy(),diffusion=diffusion,dt=step,decay=decay,source=np.asarray(source),boundary=boundary_velocity));return result
  n.ds._transport=traced
  tr=n.run('debelly',L,state,g.response_cue(job['case']),.5,interval=.5)
  n.ds._transport=original;total=sum(state[k].sum() for k in ['mca_bound','mca_unbound'])
  adv=sum(op(c['new'],c['velocity'],0.,c['boundary'],dx)@q/total for c in calls)
  dif=sum(op(c['new'],np.zeros(201),c['diffusion'],(0.,0.),dx)@q/total for c in calls)
  exc=sum((-c['decay']*c['new']+c['source'])@q/total for c in calls)
  actual=sum((tr.state[k]-state[k])@q for k in ['mca_bound','mca_unbound'])/(total*.5)
  row=dict(length=L,time=float(z['pulse_time'][i]),advection_moment_rate=float(adv),diffusion_moment_rate=float(dif),staggered_exchange_moment_rate=float(exc),net_native_moment_rate=float(actual),budget_residual=float(actual-adv-dif-exc))
  assert abs(row['budget_residual'])<1e-11;rows.append(row)
r.write(OUT/'transport_budget.json',dict(generated_at=n.now(),definition='Exact native 0.5-unit split update from saved pulse states. Total-MCA normalized nominal first moment; backward-Euler spatial terms evaluated on each newly solved field. Binding/unbinding source terms are lagged in the manuscript solver, so they leave a transient splitting contribution. The subsequent spatially uniform mass closure has zero first moment on the symmetric grid. No equations or native parameters changed.',rows=rows));print(rows)
