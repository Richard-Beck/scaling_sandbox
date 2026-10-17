"""Complete the validated campaign; excluded rows count as false in status aggregation."""
from countercue_full_lengths import *

def finalize(m,groups,round):
 cases,rows=endpoints(m,groups);assert len(cases)==462 and len(rows)==924
 selected_checks=[]
 for c in cases:
  if not c['eligible'] or c['restricted']:continue
  ts=c['tests'];assert all(t['success'] for t in ts if t['full_success'])
  for name,q in c['result'].items():
   if q['status']=='success':assert q['relative_bracket_width']<=.01000000001 and not q['lower_grid_limited']
   if c['canonical_index']!=c['case_index']:continue
   selected=next(t for t in ts if t['duration']==q['selected_duration']);prefix=f"trial{q['selected_trace_row']}_"
   with np.load(r.ROOT/q['selected_trace']) as z:
    arrays={k[len(prefix):]:z[k] for k in z.files if k.startswith(prefix)}
    assert all(np.all(np.isfinite(a)) for a in arrays.values())
    expected=score({k:v for k,v in selected.items() if k not in ['full_success','initial_orientation','maximum_relative_error_last_500']},arrays,c['initial_orientation'])
    assert expected['full_success']==selected['full_success']
   selected_checks.append(dict(case_index=c['case_index'],definition=name,verified=True))
  if c['result']['full_polarity']['status']=='success':assert c['result']['full_polarity']['endpoint']>=c['result']['orientation']['endpoint']
 # All newly required original duration grids must be present. Refinements add to these.
 for c in m['cases']:
  if c['eligible'] and not c['restricted']:
   assert set(r.PROTOCOL['duration_grids'][c['model']]['refined_durations']).issubset({t['duration'] for t in groups[c['canonical_index']]})
  assert r.sha(r.ROOT/c['checkpoint'])==c['checkpoint_sha256']
 with gzip.open(RESULTS/'records.json.gz','wt') as stream:json.dump(dict(generated_at=r.now(),cases=cases,manifest=m),stream,default=r.convert,allow_nan=False)
 r.csvwrite(RESULTS/'endpoints.csv',rows)
 validation=dict(generated_at=r.now(),passed=True,conditions=len(cases),readouts=len(rows),distinct_simulated_conditions=len(groups),distinct_trials=sum(len(v) for v in groups.values()),new_base_trials=sum(len(j['durations']) for j in m['base_jobs']),new_refinement_trials=sum(len(j['durations']) for k in range(round) for j in read(OUT/f'jobs_{k}.json')),coverage_complete=True,exact_checkpoint_hashes_verified=True,source_hashes_verified=True,cue_checks_passed=True,full_gate_implies_original_gate=True,selected_trials_verified=len(selected_checks),selected_native_arrays_all_finite=True,status_counts_by_definition={name:dict(Counter(x['status'] for x in rows if x['definition']==name)) for name in ['orientation','full_polarity']},nonmonotonic_readouts=sum(bool(q.get('nonmonotonic_success_pattern',False)) for q in rows),successful_brackets_all_within_one_percent=True,all_gate_samples_checked=True)
 write(RESULTS/'validation.json',validation);write(RESULTS/'selected_trace_checks.json',selected_checks)
 status=dict(validation,complete=True,producer='finished',stage='complete')
 from countercue_full_lengths_plot import plot
 plot(cases,status)
 write(RESULTS/'status.json',status)
 write(RESULTS/'artifact_hashes.json',{f.name:r.sha(f) for f in sorted(RESULTS.iterdir()) if f.suffix in ['.png','.pdf','.svg','.csv']})
 print(json.dumps(validation),flush=True)


if __name__=="__main__":
 m,groups,completed,missing=collect(1)
 finalize(m,groups,1)
 write(RESULTS/"completion_provenance.json",dict(generated_at=r.now(),source_sha256=r.sha(Path(__file__).resolve()),frozen_analysis_source_preserved=True,change="Count absent nonmonotonic flags on ineligible/restricted CSV rows as false; scientific definitions and simulation data unchanged."))
 launch=read(OUT/"launch.json");launch.update(stage="complete",producer="finished");write(OUT/"launch.json",launch)
