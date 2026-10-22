from pathlib import Path
import csv,json,hashlib,sys
import numpy as np
import pandas as pd
root=Path('/home/4473331/projects/scaling_sandbox/polarity1d/manuscript')
out=root/'figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/supporting/data'
linear_path=root/'data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv'
function_path=root/'figures/length_mechanisms_20261008/dev/linear_then_ring_checkpoint_20261017/supporting/data/conditional_functional_readouts.csv'
run=root/'data/runs/chemical_ring_functional_v2_20261017'
a=pd.read_csv(linear_path);b=pd.read_csv(run/'results/basal.csv');m=pd.read_csv(run/'results/matrix.csv');f=pd.read_csv(function_path)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
comparisons=[]
for mod in ['wave_pinning','goryachev','holmes_model3']:
 for _,line in a[a.model==mod].iterrows():
  matches=b[(b.model==mod)&(b.epsilon==.005)&(b.circumference==2*line.domain_length)]
  if not len(matches):continue
  ring=matches.iloc[0];status=json.loads(Path(ring.status_path).read_text());lp=root/line.checkpoint;rp=Path(status['checkpoint'])
  assert sha(lp)==line.checkpoint_sha256;assert sha(rp)==ring.checkpoint_sha256
  with np.load(lp) as z:lstate=z['state'];ln=lstate.shape[-1]
  with np.load(rp) as z:state=z['state'];n=state.shape[-1]
  center=(status['core']['peak_cell']+n//2)%n
  norm=7.5 if mod=='holmes_model3' else 1.
  pointgain=float(state[1,center])/norm
  patchgain=m[(m.model==mod)&(m.circumference==2*line.domain_length)&(m.fraction==.5)].iloc[0].initial_cue_weighted_gain
  linrepair=f[(f.model==mod)&(f.length_native==line.domain_length)]
  rec=dict(model=mod,linear_L=line.domain_length,ring_C=ring.circumference,linear_cells=ln,linear_dx=line.domain_length/ln,ring_cells=n,ring_dx=ring.dx,
   linear_peak=line.old_pole_peak_active,ring_peak=ring.basal_peak,peak_ring_over_linear=ring.basal_peak/line.old_pole_peak_active,
   linear_core_width=line.old_pole_width,ring_full_core_width=ring.basal_width,ring_half_width_over_linear=.5*ring.basal_width/line.old_pole_width,
   linear_excess_core_amount=line.old_pole_excess_active,ring_full_excess_core_amount=ring.basal_core_excess_amount,ring_half_excess_over_linear=.5*ring.basal_core_excess_amount/line.old_pole_excess_active,
   linear_absolute_core_amount=line.old_pole_total_active,ring_full_absolute_core_amount=ring.basal_core_amount,ring_half_absolute_over_linear=.5*ring.basal_core_amount/line.old_pole_total_active,
   linear_opposite_point_gain=line.opposite_inactive/norm,ring_antipodal_point_gain=pointgain,ring_antipodal_patch_gain=patchgain,
   point_gain_ring_over_linear=pointgain/(line.opposite_inactive/norm),patch_gain_over_ring_point=patchgain/pointgain,
   linear_conditional_repair=None if not len(linrepair) else linrepair.iloc[0].immediate_net_core_repair_per_damage_native,
   ring_basal_conditional_repair=ring.net_repair_per_damage,
   repair_ring_over_linear=None if not len(linrepair) else ring.net_repair_per_damage/linrepair.iloc[0].immediate_net_core_repair_per_damage_native,
   linear_checkpoint=str(lp),linear_checkpoint_sha256=sha(lp),ring_checkpoint=str(rp),ring_checkpoint_sha256=sha(rp))
  comparisons.append(rec)
p=pd.DataFrame(comparisons);p.to_csv(out/'ring_functional_ring_linear_mirror_comparison.csv',index=False)
cases=[('wave_pinning',40,.4,'common-C40-d16'),('goryachev',40,.4,'common-C40-d16'),('holmes_model3',40,.4,'common-C40-d16'),('holmes_model3',40,.3,'same-state-Holmes-distance-control'),('holmes_model3',30,.4,'Holmes-cue-site-single-peak-size-control'),('goryachev',40,.5,'Goryachev-antipodal-smaller'),('goryachev',80,.5,'Goryachev-antipodal-larger-return'),('holmes_model3',40,.5,'antipodal-two-site-vector-cancellation-diagnostic')]
records=[]
for mod,C,frac,role in cases:
 x=m[(m.model==mod)&(m.circumference==C)&(m.fraction==frac)].iloc[0];assert x.status=='complete';status=json.loads(Path(x.status_path).read_text())
 for phase,col in [('cue','cue_path'),('recovery','recovery_path')]:
  path=Path(x[col]);
  with np.load(path) as z:
   species=z['states'].shape[1] if 'states' in z.files else z['snapshot_states'].shape[1]
   records.append(dict(model=mod,circumference=C,fraction=frac,short_route=x.short_route,long_route=x.long_route,role=role,
    phase=phase,phase_duration=float(z['time'][-1]),phase_onset_from_challenge=0 if phase=='cue' else 1000,
    samples=len(z['time']),cells=int(x.cells),native_fields=species,schema='full native states' if 'states' in z.files else 'primary history plus native snapshots',
    stored_array_shape=str(z['states'].shape if 'states' in z.files else z['primary'].shape),path=str(path),logical_bytes=path.stat().st_size,allocated_bytes=path.stat().st_blocks*512,sha256=sha(path),
    checkpoint_sha256=x.checkpoint_sha256,incumbent_peak_cell=status['core']['peak_cell'],cue_center_cell=int(x.cue_center_cell),
    cue_peak={'wave_pinning':.1,'goryachev':6.4,'holmes_model3':1}[mod],cue_full_support=3,
    final_location_association=x.final_location_association))
e=pd.DataFrame(records);e.to_csv(out/'ring_functional_vector_example_set.csv',index=False)
print('13 native mirror-matched basal pairs; 8 retained response conditions,16 phase files,bytes',e.logical_bytes.sum(),'allocated',e.allocated_bytes.sum())
print('Max relative deviations percent:')
for mod,x in p.groupby('model'):
 print(mod,{k:float(np.nanmax(abs(x[k]-1))*100) for k in ['peak_ring_over_linear','ring_half_width_over_linear','ring_half_excess_over_linear','point_gain_ring_over_linear','patch_gain_over_ring_point']})
 print('repair',None if x.repair_ring_over_linear.isna().all() else float(np.nanmax(abs(x.repair_ring_over_linear-1))*100))
