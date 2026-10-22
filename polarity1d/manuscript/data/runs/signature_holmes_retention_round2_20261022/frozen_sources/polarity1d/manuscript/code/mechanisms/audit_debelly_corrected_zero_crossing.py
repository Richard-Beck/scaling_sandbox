"""Independent saved-array and checkpoint audit of threshold experiment."""
from pathlib import Path
import sys,json,numpy as np
from debelly_corrected_zero_crossing import OUT,n,r,base,LENGTHS
def z_lengths():return LENGTHS
import debelly_corrected_reversal as corrected

def main():
 wp=r.read(OUT/'withdrawal_provenance.json');assert n.sha(Path(wp['source']))==wp['sha256']==n.sha(OUT/'withdrawal_producer.py')
 prov=r.read(OUT/'provenance.json');assert prov['model_sha256']==n.sha(Path(n.ds.__file__))
 cases=r.read(corrected.DEFAULT/'manifest.json')['cases'];summaries=[];checked=[]
 for L in [18.,23.,39.]:
  for delta in [-100.,-50.,0.]:
   with np.load(OUT/f'withdrawal_L{L:g}_offset{delta:g}.npz') as zz:
    x=zz['x'];q=2*x/L-1;rr=zz['recovery_field_rac'][-1];bb=zz['recovery_field_mca_bound'][-1]
    summaries.append(dict(length=L,offset=delta,final_rac=float(rr@q/rr.sum()),final_mca=float(bb@q/bb.sum()),qualified_reversal=('preparation_ineligible' not in zz.files)))
 for row in summaries:
  L=row['length'];case=next(c for c in cases if c['model']=='debelly' and c['condition']=='fixed_slope' and c['length']==L);short=not case['eligible'];cp=n.ROOT/case['checkpoint'] if not short else None
  if not short:
   assert n.sha(cp)==case['checkpoint_sha256'];state,_=n.read_checkpoint(cp,'debelly')
  file=OUT/f'withdrawal_L{L:g}_offset{row["offset"]:g}.npz'
  with np.load(file) as z:
   x=z['x'];D=float(z['duration']);assert z['pulse_time'][-1]==D and z['recovery_time'][-1]==10000
   for k in n.DB_NAMES:
    
    if not short:np.testing.assert_array_equal(z['pulse_field_'+k][0],state[k])
    np.testing.assert_array_equal(z['recovery_field_'+k][0],z['pulse_field_'+k][-1]);assert np.all(np.isfinite(z['recovery_field_'+k]))
   for k in ['rac','mca_bound']:
    u=z['recovery_field_'+k][-1];v=float(np.sum(u*(2*x/L-1))/np.sum(u));assert abs(v-row['final_'+('rac' if k=='rac' else 'mca')])<1e-12
   b=z['pulse_field_mca_bound'][-1];v=float(np.sum(b*(2*x/L-1))/np.sum(b));
   if not short:assert (v<=0) if row['offset']==0 else (v>0)
   else:assert bool(z['preparation_ineligible']) and not row['qualified_reversal']
  checked.append(dict(length=L,offset=row['offset'],sha256=n.sha(file)))
 # A relaxed state must remain fixed before a perturbation is added.
 maxerr=0.
 for L in [12.,19.5,23.,39.]:
  b=base(L);tr=n.run('debelly',L,b,n.ConstantInput(0),.5,interval=.5)
  for k in n.DB_NAMES:maxerr=max(maxerr,float(np.max(abs(tr.state[k]-b[k]))))
 assert maxerr<1e-12,maxerr
 for L in z_lengths():
  for eps in [1e-5,1e-6]:
   with np.load(OUT/f'seed_L{L:g}_e{eps:g}.npz') as zz:
    expected=base(L);expected['mca_bound']+=eps*np.cos(np.pi*zz['x']/L)
    for k in n.DB_NAMES:
     np.testing.assert_allclose(zz['field_'+k][0],expected[k],atol=1e-13)
     assert np.all(np.isfinite(zz['field_'+k]))
    for k in ['rac','rho','mca_bound','mca_unbound']:assert np.min(zz['field_'+k])>=0
    mass=(zz['field_mca_bound']+zz['field_mca_unbound']).sum(1);assert np.max(abs(mass/mass[0]-1))<1e-10
 audit=dict(generated_at=n.now(),passed=True,exact_preparation_checkpoints=True,pulse_recovery_continuity=True,finite_fields=True,independent_final_moments=True,withdrawal_signs_correct=True,unperturbed_rest_single_step_max_error=maxerr,native_source_unchanged=True,withdrawal_archives=checked)
 r.write(OUT/'independent_audit.json',audit);print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
