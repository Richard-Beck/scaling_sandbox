"""Diagnostic necessity/sufficiency probe of initial MCA imbalance at L39."""
from pathlib import Path
import sys,numpy as np
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE));import gradient_reversal as g
r,n,e=g.r,g.n,g.e;OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly'
jobs=r.read(OUT/'jobs.json');case=jobs[4]['case'];state,prep=n.read_checkpoint(n.ROOT/case['checkpoint'],'debelly');assert n.sha(n.ROOT/case['checkpoint'])==case['checkpoint_sha256']
q=2*prep.x/39-1;M0=sum(state[k]@q for k in ['mca_bound','mca_unbound'])/sum(state[k].sum() for k in ['mca_bound','mca_unbound'])
with np.load(OUT/'L23_pass.npz') as z:
 b=z['pulse_field_mca_bound'][0]+z['pulse_field_mca_unbound'][0];q23=2*z['x']/23-1;target=b@q23/b.sum()
alpha=target/M0
for k in ['mca_bound','mca_unbound']:
 even=.5*(state[k]+state[k][::-1]);odd=.5*(state[k]-state[k][::-1]);state[k]=even+alpha*odd
assert all(np.min(state[k])>=0 for k in ['mca_bound','mca_unbound'])
d=float(sys.argv[1]) if len(sys.argv)>1 else jobs[2]['duration']
pulse=n.run('debelly',39,state,g.response_cue(case),d,interval=2.5);recovery=n.run('debelly',39,pulse.state,e.withdrawal(case),10000.,interval=2.5)
arrays={}
for phase,tr in [('pulse',pulse),('recovery',recovery)]:
 for key in ['time','orientation','contrast','x']:arrays[phase+'_'+key]=getattr(tr,key)
 for key,value in tr.fields.items():arrays[phase+'_field_'+key]=value
path=OUT/('L39_reduced_mca_memory'+(f'_pulse_{d:g}' if len(sys.argv)>1 else '')+'.npz');np.savez_compressed(path,**arrays)
r.write(path.with_suffix('.json'),dict(generated_at=n.now(),length=39,duration=d,diagnostic='Explicit initial-state perturbation only: reduce odd reflection component of bound and unbound MCA to match L23 initial total-MCA first moment, preserve both even components, total mass and every other original full-state field; then apply original native fixed slope and basal recovery. This is not a manuscript preparation or endpoint.',original_checkpoint_sha256=case['checkpoint_sha256'],original_mca_moment=float(M0),target_mca_moment=float(target),odd_component_multiplier=float(alpha),on_cue_orientation=float(pulse.orientation[-1]),tail_min_orientation=float(recovery.orientation[-201:].min()),tail_max_orientation=float(recovery.orientation[-201:].max()),sha256=n.sha(path)))
print(float(alpha),float(recovery.orientation[-1]))
