"""Independently check completed width companion and publish completion metadata."""
from pathlib import Path
import gzip,json
import numpy as np
import countercue_full_reorientation as f
import countercue_pilot as p
import runner as r

def main():
 m=f.verify();out=f.RESULTS
 records=json.load(gzip.open(out/'records.json.gz','rt'));validation=f.read(out/'validation.json');checks=[];comparison=[]
 for c in records['cases']:
  tests=c['tests'];q=c['result'];passing=[]
  for t in tests:
   error=max(abs(t['minimum_normalized_polarity_last_500']-1),abs(t['maximum_normalized_polarity_last_500']-1))
   np.testing.assert_allclose(error,t['maximum_relative_error_last_500'],rtol=1e-12,atol=1e-12)
   assert t['full_success']==bool(t['success'] and error<=.01)
   if t['full_success']:passing.append(t)
  assert q['endpoint']==(min(t['duration'] for t in passing) if passing else None)
  assert q['status']==('success' if passing else 'censored')
  relaxed=[t['duration'] for t in tests if t['success']]
  relaxed_endpoint=min(relaxed) if relaxed else None
  if passing:assert q['endpoint']>=relaxed_endpoint
  comparison.append(dict(q,original_gate_endpoint_on_combined_grid=relaxed_endpoint,full_to_original_gate_ratio_on_combined_grid=q['endpoint']/relaxed_endpoint if passing else None))
  if passing:
   assert q['relative_bracket_width']<=.01000000001
   assert not q['lower_grid_limited']
  selected=next(t for t in tests if t['duration']==q['selected_duration'])
  with np.load(r.ROOT/q['selected_trace']) as z:
   prefix=f"trial{q['selected_trace_row']}_"
   initial=float(z[prefix+'pulse_orientation'][0]);tail=z[prefix+'recovery_tail_orientation']
   checked=f.score({k:v for k,v in selected.items() if k not in ['full_success','initial_orientation','maximum_relative_error_last_500','minimum_normalized_polarity_last_500','maximum_normalized_polarity_last_500']},tail,initial)
   assert checked['full_success']==selected['full_success']
   np.testing.assert_allclose(initial,c['initial_orientation'],rtol=1e-8,atol=1e-10)
   assert all(np.all(np.isfinite(z[name])) for name in z.files if name.startswith(prefix))
  checks.append(dict(case_index=c['case_index'],selected_trial_verified=True,all_selected_native_arrays_finite=True))
 assert len(checks)==126 and validation['source_data_unchanged']
 r.csvwrite(out/'endpoints.csv',comparison)
 validation.update(original_and_full_gates_compared_on_same_grid=True,independent_endpoint_audit_passed=True,all_selected_traces_verified=True,all_selected_native_arrays_finite=True)
 f.write(out/'independent_audit.json',dict(generated_at=r.now(),passed=True,source_sha256=r.sha(Path(__file__).resolve()),cases=checks))
 f.write(out/'validation.json',validation)
 # The original finalizer duplicated generated_at in its status constructor after
 # writing and validating all scientific artifacts. Publish status after this audit.
 f.write(out/'status.json',dict(validation,generated_at=r.now(),complete=True,producer='finished',stage='complete'))
 f.write(out/'completion_metadata_provenance.json',dict(generated_at=r.now(),reason='Publish completion after independent audit; repair finalizer duplicate-key error in status metadata only.',analysis_source_preserved=True,analysis_source_sha256=m['analysis_source_sha256'],audit_source_sha256=r.sha(Path(__file__).resolve())))
 from countercue_pilot_plot_revision import dashboard
 dashboard(r.now())
 f.write(out/'artifact_hashes.json',{x.name:r.sha(x) for x in sorted(out.iterdir()) if x.suffix in ['.png','.pdf','.svg','.csv']})
 print(json.dumps(validation),flush=True)

if __name__=='__main__':main()
