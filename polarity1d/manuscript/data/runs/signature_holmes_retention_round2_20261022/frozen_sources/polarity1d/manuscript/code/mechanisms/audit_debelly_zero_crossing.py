"""Independent saved-array and checkpoint audit of threshold experiment."""
from pathlib import Path
import sys,json,numpy as np
from debelly_zero_crossing import OUT,n,r,base

def main():
 prov=r.read(OUT/'provenance.json');assert prov['model_sha256']==n.sha(Path(n.ds.__file__))
 jobs=r.read(n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly/jobs.json');summaries=r.read(OUT/'withdrawal_summary.json');checked=[]
 for row in summaries:
  L=row['length'];case=next(j['case'] for j in jobs if j['case']['length']==L);cp=n.ROOT/case['checkpoint'];assert n.sha(cp)==case['checkpoint_sha256'];state,_=n.read_checkpoint(cp,'debelly')
  file=OUT/f'withdrawal_L{L:g}_offset{row["offset"]:g}.npz'
  with np.load(file) as z:
   x=z['x'];D=float(z['duration']);assert z['pulse_time'][-1]==D and z['recovery_time'][-1]==10000
   for k in n.DB_NAMES:
    np.testing.assert_array_equal(z['pulse_field_'+k][0],state[k]);np.testing.assert_array_equal(z['recovery_field_'+k][0],z['pulse_field_'+k][-1]);assert np.all(np.isfinite(z['recovery_field_'+k]))
   for k in ['rac','mca_bound']:
    u=z['recovery_field_'+k][-1];v=float(np.sum(u*(2*x/L-1))/np.sum(u));assert abs(v-row['final_'+('rac' if k=='rac' else 'mca')])<1e-12
   b=z['pulse_field_mca_bound'][-1];v=float(np.sum(b*(2*x/L-1))/np.sum(b));assert (v<=0) if row['offset']==0 else (v>0)
  checked.append(dict(length=L,offset=row['offset'],sha256=n.sha(file)))
 # A relaxed state must remain fixed before a perturbation is added.
 maxerr=0.
 for L in [12.,19.5,23.,39.]:
  b=base(L);tr=n.run('debelly',L,b,n.ConstantInput(0),.5,interval=.5)
  for k in n.DB_NAMES:maxerr=max(maxerr,float(np.max(abs(tr.state[k]-b[k]))))
 assert maxerr<1e-12,maxerr
 audit=dict(generated_at=n.now(),passed=True,exact_preparation_checkpoints=True,pulse_recovery_continuity=True,finite_fields=True,independent_final_moments=True,withdrawal_signs_correct=True,unperturbed_rest_single_step_max_error=maxerr,native_source_unchanged=True,withdrawal_archives=checked)
 r.write(OUT/'independent_audit.json',audit);print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
