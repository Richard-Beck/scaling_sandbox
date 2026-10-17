"""Supplement remaining >5% threshold brackets without changing frozen pilot sources."""
from pathlib import Path
import gzip,json,math
import numpy as np
import countercue_pilot as p
import runner as r

p.verify()
with gzip.open(p.RESULTS/'records.json.gz','rt') as f:records=json.load(f)
jobs=p.read(p.OUT/'refinement_jobs.json');previous_sha=r.sha(p.OUT/'refinement_jobs.json');start=len(jobs);targets=[]
for c in records['cases']:
 if c['case_index']!=c['canonical_index']:continue
 q=c['result']
 if q['status']!='success':continue
 lo=q['last_failing_duration_below_endpoint'];hi=q['endpoint']
 if lo<=0 or (hi-lo)/hi<=.05:continue
 steps=math.ceil((hi-lo)/(.01*lo));durations=np.linspace(lo,hi,steps+1)[1:-1].tolist();existing={t['duration'] for t in c['tests']};durations=[d for d in durations if d not in existing]
 for i in range(0,len(durations),10):jobs.append(dict(slot=len(jobs),case_index=c['case_index'],durations=durations[i:i+10]))
 targets.append(dict(case_index=c['case_index'],model=c['model'],length=c['length'],condition=c['condition'],previous_bracket=[lo,hi],additional_trials=len(durations)))
assert targets
p.write(p.OUT/'refinement_jobs.json',jobs)
p.write(p.OUT/'supplemental_refinement_provenance.json',dict(generated_at=r.now(),source_sha256=r.sha(Path(__file__).resolve()),previous_jobs_sha256=previous_sha,updated_jobs_sha256=r.sha(p.OUT/'refinement_jobs.json'),reason='Resolve >5% remaining brackets; uniform local grid, no assumption of globally monotone success',targets=targets))
p.status('final_bracket_refinement_running',base_complete=True,additional_jobs=len(jobs)-start)
array=p.submit('refine_worker',array=f'{start}-{len(jobs)-1}');final=p.submit('finalize',dependency=array,time='00:30:00')
launch=p.read(p.OUT/'launch.json');launch.update(stage='final_bracket_refinement_running',supplemental_refinement_array=array,supplemental_finalize_job=final);p.write(p.OUT/'launch.json',launch)
print(json.dumps(dict(array=array,final=final,targets=targets)),flush=True)
