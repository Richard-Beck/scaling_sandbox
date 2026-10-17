"""Verify completed assay coverage, conserved pools, controls and protected report hashes."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE));import run_assays;import fixed_peaks
for task in run_assays.tasks():
 if task[0]=='pool':_,model,L,seed,n=task;p=HERE/f'data/pool_{model}_L{L}_s{seed}_n{n}.json'
 elif task[0]=='finite':_,model,L,n=task;p=HERE/f'data/finite_{model}_L{L}_n{n}.json'
 else:_,L,n=task;p=HERE/f'data/spring_L{L:g}_n{n}.json'
 assert p.exists(),p
for model,L,seed,n in fixed_peaks.tasks():assert (HERE/f'data/fixed_{model}_L{L}_s{seed}_n{n}.json').exists()
pool=[r for p in (HERE/'data').glob('pool_*.json') for r in json.loads(p.read_text())]
fixed=[r for p in (HERE/'data').glob('fixed_*.json') for r in json.loads(p.read_text())]
finite=[json.loads(p.read_text()) for p in (HERE/'data').glob('finite_*.json')]
spring=[json.loads(p.read_text()) for p in (HERE/'data').glob('spring_*.json')]
assert max(r.get('mass_error',0) for r in pool+fixed+spring+[a for r in finite for a in r['assays']])<1e-5
control=[r for r in fixed if r['arm']=='fixed_peak']
assert max(abs(r['initial_active_peak_mass']/r['reference_active_peak_mass']-1) for r in control)<1e-10
assert max(abs(r['initial_peak_height']/r['reference_peak_height']-1) for r in control)<1e-3
assert max(r['inactive_mean_relative_error'] for r in pool if r['arm']=='matched_pool')<1e-5
refinements=[]
for r in fixed:
 if r['n']!=401:continue
 c=next(x for x in fixed if x['n']==201 and all(x[k]==r[k] for k in ['model','L','seed','arm','assay']))
 key='competition_resolution_latency_s' if r['assay'].startswith('competing') else 'reversal_latency_s'
 assert (c.get(key) is None)==(r.get(key) is None)
 if r['assay'].startswith('gradient'):assert c['reversal_conditioned']==r['reversal_conditioned']
 err=None if r.get(key) is None else abs(c[key]-r[key])/max(r[key],1)
 assert err is None or err<.06
 refinements.append(dict(model=r['model'],L=r['L'],arm=r['arm'],assay=r['assay'],coarse=c.get(key),fine=r.get(key),relative_error=err))
for r in finite:
 ds=r['assays'];successes=[a['duration'] for a in ds if a['success']]
 if successes:
  upper=min(successes)
  assert all(a['success'] for a in ds if a['duration']>=upper)
for r in spring:
 c=next(x for x in spring if x['L']==r['L'] and x['n']==101)
 assert c['resolution_latency']==r['resolution_latency'] and c['established']==r['established']
protected={}
for line in (HERE/'report_hashes_before.txt').read_text().splitlines():
 digest,name=line.split(maxsplit=1);current=hashlib.sha256((ROOT/name).read_bytes()).hexdigest();assert digest==current,name;protected[name]=current
for stem in ['figure1','figure2','supplement1','supplement2','supplement3_tracking']:
 for ext in ['png','pdf']:assert (HERE/f'figures/{stem}.{ext}').stat().st_size>1000
assert (HERE/'figures/figure3.png').stat().st_size>1000
summary=dict(assay_tasks=122,fixed_peak_tasks=60,pool_summaries=len(pool),finite_pulse_assays=sum(len(r['assays']) for r in finite),spring_competition_assays=len(spring),maximum_relative_mass_error=max(r.get('mass_error',0) for r in pool+fixed+spring+[a for r in finite for a in r['assays']]),maximum_fixed_peak_mass_error=max(abs(r['initial_active_peak_mass']/r['reference_active_peak_mass']-1) for r in control),maximum_fixed_peak_height_error=max(abs(r['initial_peak_height']/r['reference_peak_height']-1) for r in control),maximum_matched_mean_error=max(r['inactive_mean_relative_error'] for r in pool if r['arm']=='matched_pool'),refinement=refinements,protected_report_hashes=protected,passed=True)
(HERE/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['refinement','protected_report_hashes']},indent=2))
