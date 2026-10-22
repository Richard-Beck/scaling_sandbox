"""Independent array audit and concise comparison for the boundary correction."""
from pathlib import Path
import json,numpy as np
from debelly_boundary_resolution import OUT,n,r,base

def main():
 jobs=r.read(OUT/'jobs.json');maxdrift=0.;records=[]
 for path,digest in r.read(OUT/'provenance.json')['sources'].items():assert n.sha(Path(path))==digest
 for job in jobs:
  L=job['length'];v=job['variant'];s=r.read(OUT/f'{v}_L{L:g}_summary.json');re=r.read(OUT/f'{v}_L{L:g}_recovery.json')
  for label,digest in [('preparation',s['preparation_sha256']),('cue',s['cue_sha256']),('recovery',re['sha256'])]:assert n.sha(OUT/f'{v}_L{L:g}_{label}.npz')==digest
  prep=np.load(OUT/f'{v}_L{L:g}_preparation.npz');cue=np.load(OUT/f'{v}_L{L:g}_cue.npz');rec=np.load(OUT/f'{v}_L{L:g}_recovery.npz');x=cue['x'];q=2*x/L-1
  initial=base(L);initial['mca_bound']-=1e-5*np.cos(np.pi*x/L)
  i=int(np.searchsorted(cue['time'],s['crossing_time']));assert cue['time'][i]==s['crossing_time']
  for k in n.DB_NAMES:
   np.testing.assert_array_equal(prep['field_'+k][0],initial[k]);np.testing.assert_array_equal(cue['field_'+k][0],prep['state_'+k]);np.testing.assert_array_equal(rec['field_'+k][0],cue['field_'+k][i])
  for z in [prep,cue,rec]:
   for k in ['rac','rho','mca_bound','mca_unbound']:assert np.all(np.isfinite(z['field_'+k])) and np.min(z['field_'+k])>=0
   mass=(z['field_mca_bound']+z['field_mca_unbound']).sum(1);maxdrift=max(maxdrift,float(np.max(abs(mass/mass[0]-1))))
  assert maxdrift<1e-10
  hold=prep['time']>=50000
  holdmax=max(float(np.max(np.ptp(prep['field_'+k][hold],axis=0))) for k in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension']);assert holdmax<1e-8
  P=np.sum(cue['field_mca_bound']*q,axis=1)/cue['field_mca_bound'].sum(1);assert P[i]<=0 and P[i-1]>0
  tail=rec['time']>=9500;rac=rec['field_rac'];PR=(rac*q).sum(1)/rac.sum(1);passed=bool(PR[tail].min()>.05 and np.ptp(PR[tail])<.005 and np.ptp(rac[tail],axis=1).min()>=.05);assert passed==re['retained_right_orientation']
  for a in s['rows']:
   j=int(np.searchsorted(cue['time'],a['time']));actual=-(P[j+1]-P[j])/.5;assert abs(actual-a['reorientation_rate'])<1e-11
   if v=='independent_endpoints':assert a['used_left_velocity']==a['raw_left_velocity'] and a['used_right_velocity']==a['raw_right_velocity']
   # Displacement-mapped centroid, plus midpoint-relative movement: distinct from nominal first moment.
   b=cue['field_mca_bound'][j:j+2];u=cue['field_displacement'][j:j+2];mapped=((x+u)*b).sum(1)/b.sum(1);mid=.5*(L+cue['field_left_boundary_displacement'][j:j+2,0]+cue['field_right_boundary_displacement'][j:j+2,0]);lab=-float(np.diff(mapped)[0]/.5);relative=-float(np.diff(mapped-mid)[0]/.5)
   records.append(dict(variant=v,length=L,stage=a['stage'],cue_time=a['time'],crossing_time=s['crossing_time'],prepared_left_rac=s['preparation_left_polarity'],nominal_centroid_speed=a['centroid_speed'],displacement_mapped_lab_centroid_speed=lab,displacement_mapped_midpoint_relative_speed=relative,direct_rho_velocity_rms=a['direct_rho_rms'],endpoint_velocity_rms=a['endpoint_rms'],membrane_drag_rms=a['membrane_drag_rms'],retained_right_after_zero=passed,final_right_rac=float(PR[-1]),mechanical_reconstruction_error=a['reconstruction_error']))
 n.csvwrite(OUT/'comparison.csv',records)
 audit=dict(generated_at=n.now(),passed=True,producer='finished',variants=2,lengths=5,preparations=10,cue_trials=10,recoveries=10,identical_initial_seed_verified=True,each_response_restarts_own_variant_state=True,all_preparation_10000_holdout_excursions_below_1e8=True,positive_finite_chemical_fields=True,maximum_relative_total_mca_drift=maxdrift,independent_crossing_and_recovery_gates=True,corrected_endpoint_commands_unmodified=True,source_hashes_unchanged=True,all_ten_withdrawals_at_zero_retain_right_orientation=all(a['retained_right_after_zero'] for a in records))
 r.write(OUT/'independent_audit.json',audit);print(json.dumps(audit,indent=2))
 for L in [23.,39.]:
  for v in ['legacy_centered','independent_endpoints']:
   a=next(a for a in records if a['length']==L and a['variant']==v and a['stage']=='zero');print(a)
if __name__=='__main__':main()
