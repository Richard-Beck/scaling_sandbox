"""Causal diagnostic: hold MCA fixed during cue only, native basal recovery."""
from pathlib import Path
import sys,numpy as np
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE));import gradient_reversal as g
r,n,e=g.r,g.n,g.e
OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly'
slot=int(sys.argv[1]);job=r.read(OUT/'jobs.json')[2*slot];case=job['case'];L=case['length'];d=job['duration']
state,prep=n.read_checkpoint(n.ROOT/case['checkpoint'],'debelly')
assert n.sha(n.ROOT/case['checkpoint'])==case['checkpoint_sha256']
original=n.ds._transport
# Suppress MCA exchange, diffusion, advection during pulse; this isolates the
# chemical/transient mechanical response without rewriting adhesion memory.
n.ds._transport=lambda values,*args,**kwargs:values.copy()
pulse=n.run('debelly',L,state,g.response_cue(case),d,interval=2.5)
n.ds._transport=original
recovery=n.run('debelly',L,pulse.state,e.withdrawal(case),10000.,interval=2.5)
arrays={}
for phase,tr in [('pulse',pulse),('recovery',recovery)]:
 for key in ['time','orientation','contrast','x']:arrays[phase+'_'+key]=getattr(tr,key)
 for key,value in tr.fields.items():arrays[phase+'_field_'+key]=value
path=OUT/f'L{L:g}_mca_clamp.npz';np.savez_compressed(path,**arrays)
r.write(path.with_suffix('.json'),dict(generated_at=n.now(),length=L,duration=d,diagnostic='Explicit causal perturbation: all MCA redistribution/exchange suppressed during pulse only; original native dynamics restored for basal recovery. The original exact checkpoint, cue, dose, Rac/Rho/mechanics and recovery horizon are retained. This is not part of the response comparison.',checkpoint_sha256=case['checkpoint_sha256'],sha256=n.sha(path),on_cue_orientation=float(pulse.orientation[-1]),final_orientation=float(recovery.orientation[-1]),tail_min_orientation=float(recovery.orientation[-201:].min()),tail_max_orientation=float(recovery.orientation[-201:].max())))
print(path,float(pulse.orientation[-1]),float(recovery.orientation[-1]))
