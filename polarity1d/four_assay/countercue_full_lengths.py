"""Three patch cues, both polarity gates, all canonical lengths and exact v11 release states."""
from pathlib import Path
import argparse,gzip,json,math,os,shlex,subprocess,sys,shutil
from collections import Counter
import numpy as np
import runner as r
import countercue_pilot as p

CODE=Path(__file__).resolve().parent
OUT=CODE/'pilots/width_full_lengths_20261002';WORK=OUT/'work';RESULTS=OUT/'results'
REUSE=CODE/'pilots/full_reorientation_width_20261002/results/records.json.gz'
PLOT=CODE/'countercue_full_lengths_plot.py'
EXPERIMENTS=[('fractional','peak_matched'),('physical','peak_matched'),('physical','integral_matched')]

class Cue:
 def __init__(self,case,duration):self.case=case;self.stop=float(duration);self.event_times=(self.stop,)
 def __call__(self,t,x,L):
  x=np.asarray(x,float);c=self.case
  if t>=self.stop:return np.zeros_like(x)
  w=.2*(L if c['width_mode']=='fractional' else c['reference_domain_length'])
  if c.get('physical_patch_policy')=='restrict' and w>L/2+1e-8:raise ValueError('Physical patch exceeds half the instantaneous domain')
  peak=c['A']
  if c['dose_mode']=='integral_matched':
   if L>=w:peak*=L/c['reference_domain_length']
   else:
    # Integral of the physical cosine patch on the actual cell domain.
    area=.5*(L+w/np.pi*np.sin(np.pi*L/w))
    peak*=.1*L/area
  distance=(L-x)/w
  return np.where((distance>=0)&(distance<=1),peak*.5*(1+np.cos(np.pi*np.clip(distance,0,1))),0.)

def write(path,v):p.write(path,v)
def read(path):return p.read(path)
def hashes():
 h=p.source_hashes()
 for path in [Path(__file__).resolve(),PLOT,REUSE]:h[str(path.relative_to(r.ROOT))]=r.sha(path)
 return h

def verify():
 p.verify();prov=read(OUT/'provenance.json');assert hashes()==prov['sources'],'Frozen source/input changed'
 assert r.sha(OUT/'manifest.json')==prov['manifest_sha256']
 return read(OUT/'manifest.json')

def submit(action,array=None,dependency=None,round=0):
 command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__).resolve()),action,'--round',str(round)])
 if array is not None:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
 opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=02:00:00','--job-name=width-all-'+action,'--chdir='+str(r.ROOT),'--output='+str(WORK/'logs'/f'{action}_r{round}_%A_%a.out'),'--error='+str(WORK/'logs'/f'{action}_r{round}_%A_%a.err')]
 if array is not None:opts+=['--array='+array]
 if dependency:opts+=['--dependency=afterany:'+dependency]
 return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]

def score(t,z,initial):
 tail=np.asarray(z['recovery_tail_orientation']);assert len(tail)==201 and np.all(np.isfinite(tail))
 np.testing.assert_allclose(float(z['pulse_orientation'][0]),initial,rtol=1e-8,atol=1e-10)
 error=float(np.max(abs(tail/abs(initial)-1)))
 assert initial<-.05
 return dict(**t,full_success=bool(t['success'] and error<=.01),initial_orientation=initial,maximum_relative_error_last_500=error)

def launch(policy):
 assert not OUT.exists()
 for path in [RESULTS,WORK/'chunks',WORK/'logs',OUT/'frozen_sources']:path.mkdir(parents=True,exist_ok=True)
 p.verify();original=p.records();preps={(c['model'],c['length']):c for c in original['preparations'] if c['assay']=='finite_countercue'}
 scales=read(CODE/'reference/continuation_20261001/fixed_scales.json')
 reuse=json.load(gzip.open(REUSE,'rt'));lookup={(c['model'],c['length'],c['width_mode'],c['dose_mode']):c for c in reuse['cases'] if c['fraction']==.2}
 cases=[];reused={};jobs=[];checks=[]
 for model in r.MODELS:
  ref=p.LENGTHS[model][1];_,tr=r.read_checkpoint(r.ROOT/preps[model,ref]['checkpoint'],model)
  domain=float(tr.fields['length'][-1]) if model=='spring' else ref
  for task in [t for t in r.tasks() if t['model']==model and t['assay']=='finite_countercue']:
   L=task['length'];prep=preps[model,L];_,tr=r.read_checkpoint(r.ROOT/prep['checkpoint'],model)
   initial=float(tr.orientation[-1]);actual=float(tr.fields['length'][-1]) if model=='spring' else L
   fixed=next(s['fixed_scales'] for s in scales if s['model']==model and s['length']==L and s['event']=='release')
   for mode,dose in EXPERIMENTS:
    eligible=bool(prep['eligible']);restricted=bool(policy=='restrict' and mode=='physical' and .2*domain>(L if model=='spring' else actual)/2+1e-8)
    c=dict(case_index=len(cases),model=model,length=L,cells=task['cells'],A=r.PROTOCOL['normalization']['amplitude_by_model'][model],family='width',fraction=.2,width_mode=mode,dose_mode=dose,condition=mode+'_'+dose,reference_length=ref,reference_domain_length=domain,checkpoint=prep['checkpoint'],checkpoint_sha256=prep['checkpoint_sha256'],fixed_scales=fixed,physical_patch_policy=policy,initial_orientation=initial,initial_domain_length=actual,eligible=eligible,restricted=restricted,patch_exceeds_half_initial_domain=bool(mode=='physical' and .2*domain>actual/2),patch_exceeds_initial_domain=bool(mode=='physical' and .2*domain>actual))
    c['canonical_index']=c['case_index']
    if model!='spring' and L==ref:
     c['canonical_index']=next((d['case_index'] for d in cases if d['model']==model and d['length']==L and d['width_mode']=='fractional'),c['case_index'])
    cases.append(c)
    if not eligible or restricted or c['canonical_index']!=c['case_index']:continue
    if (model,L,mode,dose) in lookup:
     old=lookup[model,L,mode,dose];np.testing.assert_allclose(old['initial_orientation'],initial,rtol=1e-8,atol=1e-10)
     assert old['checkpoint_sha256']==c['checkpoint_sha256']
     reused[str(c['case_index'])]=old['tests']
    else:
     ds=r.PROTOCOL['duration_grids'][model]['refined_durations']
     for start in range(0,len(ds),10):jobs.append(dict(slot=len(jobs),case_index=c['case_index'],durations=ds[start:start+10]))
    # Check cues on actual and bounded Spring domains, including truncation.
    for length in ([actual,L+.1] if model=='spring' else [L]):
     x=np.linspace(0,length,20001);cue=Cue(c,10)(0,x,length)
     w=.2*(length if mode=='fractional' else domain);span=min(length,w)
     area=.5*(span+w/np.pi*np.sin(np.pi*span/w))
     if c.get('physical_patch_policy')=='restrict' and w>L/2+1e-8:raise ValueError('Physical patch exceeds half the instantaneous domain')
     peak=c['A'] if dose=='peak_matched' else c['A']*.1*length/area
     np.testing.assert_allclose(np.trapezoid(cue,x),peak*area,rtol=2e-7,atol=1e-12)
     if dose=='integral_matched':np.testing.assert_allclose(np.trapezoid(cue,x),c['A']*.1*length,rtol=2e-7,atol=1e-12)
     assert np.all(cue>=0) and np.all(Cue(c,10)(10,x,length)==0)
     if (model,L,mode,dose) in lookup:np.testing.assert_allclose(cue,p.Cue(c,10)(0,x,length),rtol=1e-14,atol=1e-14)
     checks.append(dict(case_index=c['case_index'],domain_length=length,integral=float(np.trapezoid(cue,x)),expected_integral=peak*area))
 manifest=dict(created_at=r.now(),cases=cases,base_jobs=jobs,recovery_horizon=10000,gate_window=500,full_gate='All 201 original tail orientation samples within 0.99–1.01 times abs(initial orientation), plus original success gate',original_gate='Original orientation >0.05 throughout tail, contrast >=0.05, orientation temporal range <0.005',bracket_relative_resolution=.01,physical_patch_policy=policy,physical_patch_rule='Original cosine patch restricted to actual cell domain. Integral-matched peak uses actual truncated cosine integral if domain shorter than support width.',reference_lengths={model:p.LENGTHS[model][1] for model in r.MODELS},alias_rule='Static reference-length cues are identical; Spring retains distinct dynamic-domain cues',reuse_records_sha256=r.sha(REUSE),initial_states='Exact canonical revision-11 release checkpoints, no reconditioning',cue_checks=checks)
 write(OUT/'manifest.json',manifest)
 with gzip.open(OUT/'reused_tests.json.gz','wt') as stream:json.dump(reused,stream,default=r.convert,allow_nan=False)
 source_hashes=hashes()
 for rel in source_hashes:
  src=r.ROOT/rel;dest=OUT/'frozen_sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
 write(OUT/'provenance.json',dict(generated_at=r.now(),sources=source_hashes,manifest_sha256=r.sha(OUT/'manifest.json'),reused_tests_sha256=r.sha(OUT/'reused_tests.json.gz')))
 progress('base_trials_running',manifest)
 array=submit('worker',array=f'0-{len(jobs)-1}',round=-1);controller=submit('plan',dependency=array,round=0)
 write(OUT/'launch.json',dict(generated_at=r.now(),base_array=array,controller=controller,base_jobs=len(jobs),base_trials=sum(len(j['durations']) for j in jobs),reused_trials=sum(len(t) for t in reused.values()),stage='base_trials_running'))
 print(json.dumps(read(OUT/'launch.json')),flush=True)

def collect(round,partial=False):
 m=verify();prov=read(OUT/'provenance.json');assert r.sha(OUT/'reused_tests.json.gz')==prov['reused_tests_sha256']
 groups={int(k):v for k,v in json.load(gzip.open(OUT/'reused_tests.json.gz','rt')).items()};missing=[];expected=len(m['base_jobs'])
 batches=[('base',m['base_jobs'])]+[(f'r{k}',read(OUT/f'jobs_{k}.json')) for k in range(round)]
 completed=0
 for prefix,jobs in batches:
  for job in jobs:
   path=WORK/'chunks'/f'{prefix}_{job["slot"]:05d}.json'
   if not path.exists():missing.append(str(path));continue
   row=read(path);assert row['case_index']==job['case_index'];assert [t['duration'] for t in row['tests']]==job['durations']
   groups.setdefault(job['case_index'],[]).extend(row['tests']);completed+=1
 if missing and not partial:raise RuntimeError(f'Missing {len(missing)} execution chunks: '+missing[0])
 for tests in groups.values():
  tests.sort(key=lambda t:t['duration']);ds=[t['duration'] for t in tests];assert len(ds)==len(set(ds))
 return m,groups,completed,missing

def worker(slot,round):
 m=verify();base=round==-1;jobs=m['base_jobs'] if base else read(OUT/f'jobs_{round}.json');job=jobs[slot];case=m['cases'][job['case_index']]
 state,_=p.state_meta(case);p.Cue=Cue # Explicit cue adapter; native trial/recovery/diagnostics unchanged.
 prefix='base' if base else f'r{round}';path=WORK/'chunks'/f'{prefix}_{slot:05d}';arrays={};rows=[]
 for i,D in enumerate(job['durations']):
  t,z=p.trial(case,D);t=score(t,z,case['initial_orientation']);t.update(trace_row=i,trace=str(path.with_suffix('.npz').relative_to(r.ROOT)));rows.append(t)
  arrays.update({f'trial{i}_{k}':v for k,v in z.items()})
 np.savez_compressed(path.with_suffix('.npz'),**arrays);verify()
 write(path.with_suffix('.json'),dict(case_index=case['case_index'],tests=rows,generated_at=r.now()))


def endpoints(m,groups,partial=False):
 cases=[];rows=[]
 for c in m['cases']:
  tests=groups.get(c['canonical_index'],[])
  result={}
  for name,gate in [('orientation','success'),('full_polarity','full_success')]:
   if not c['eligible']:q=dict(status='ineligible',endpoint=None)
   elif c['restricted']:q=dict(status='restricted',endpoint=None)
   elif not tests:q=dict(status='pending',endpoint=None)
   else:
    passing=[t for t in tests if t[gate]];selected=passing[0] if passing else tests[-1];endpoint=selected['duration'] if passing else None
    lo=max((t['duration'] for t in tests if endpoint is not None and t['duration']<endpoint),default=0.)
    q=dict(status='success' if passing else ('pending' if partial else 'censored'),endpoint=endpoint,last_failing_duration_below_endpoint=lo if passing else None,relative_bracket_width=(endpoint-lo)/endpoint if passing else None,lower_grid_limited=bool(passing and lo==0),maximum_tested_duration=max(t['duration'] for t in tests),tested_count=len(tests),selected_duration=selected['duration'],selected_trace=selected['trace'],selected_trace_row=selected['trace_row'],selected_morphology=selected['morphology'],selected_tail_field_pass_0_1_percent=selected['tail_field_pass_0_1_percent'],selected_maximum_relative_error_last_500=selected['maximum_relative_error_last_500'],nonmonotonic_success_pattern=any(tests[k][gate] and not tests[k+1][gate] for k in range(len(tests)-1)))
   result[name]=q
   rows.append(dict(case_index=c['case_index'],canonical_index=c['canonical_index'],model=c['model'],length=c['length'],condition=c['condition'],width_mode=c['width_mode'],dose_mode=c['dose_mode'],definition=name,**q))
  cases.append(dict(**c,result=result,tests=tests))
 keys=list(dict.fromkeys(k for row in rows for k in row))
 rows=[{key:row.get(key) for key in keys} for row in rows]
 return cases,rows

def progress(stage,m=None,groups=None,completed=0,missing=None,round=0):
 if m is None:m,groups,completed,missing=collect(round,partial=True)
 if groups is None:groups={int(k):v for k,v in json.load(gzip.open(OUT/'reused_tests.json.gz','rt')).items()}
 cases,rows=endpoints(m,groups,partial=True)
 with gzip.open(RESULTS/'records.json.gz','wt') as stream:json.dump(dict(generated_at=r.now(),cases=cases,manifest=m),stream,default=r.convert,allow_nan=False)
 r.csvwrite(RESULTS/'endpoints.csv',rows)
 status=dict(generated_at=r.now(),complete=False,producer='running',stage=stage,round=round,conditions=len(cases),readouts=len(rows),completed_new_chunks=completed,pending_chunks=len(missing) if missing is not None else len(m['base_jobs']),tested_distinct_trials=sum(len(v) for v in groups.values()))
 write(RESULTS/'status.json',status)
 from countercue_full_lengths_plot import plot
 plot(cases,status)
 return status

def plan(round):
 m,groups,completed,missing=collect(round);jobs=[];targets=[]
 for c in m['cases']:
  if not c['eligible'] or c['restricted'] or c['case_index']!=c['canonical_index']:continue
  tests=groups[c['case_index']];existing={t['duration'] for t in tests};extra=set()
  for gate in ['success','full_success']:
   passing=[t for t in tests if t[gate]]
   if not passing:continue
   hi=passing[0]['duration'];lo=max((d for d in existing if d<hi),default=0.)
   if lo>0 and (hi-lo)/hi<=.01:continue
   if lo==0:ds=[hi/2**k for k in range(1,11)]
   else:
    steps=math.ceil((hi-lo)/(.01*lo));ds=np.linspace(lo,hi,steps+1)[1:-1].tolist()
   extra.update(float(d) for d in ds)
   targets.append(dict(case_index=c['case_index'],definition=gate,bracket=[lo,hi]))
  durations=sorted(extra-existing)
  for i in range(0,len(durations),10):jobs.append(dict(slot=len(jobs),case_index=c['case_index'],durations=durations[i:i+10]))
 if jobs:
  if round>=4:raise RuntimeError('Threshold refinement incomplete after four rounds; inspect unresolved targets')
  write(OUT/f'jobs_{round}.json',jobs);write(OUT/f'targets_{round}.json',targets)
  progress('threshold_refinements_running',m,groups,completed,missing,round)
  array=submit('worker',array=f'0-{len(jobs)-1}',round=round);controller=submit('plan',dependency=array,round=round+1)
  launch=read(OUT/'launch.json');launch.update(stage='threshold_refinements_running',round=round,array=array,controller=controller);write(OUT/'launch.json',launch)
  print(json.dumps(launch),flush=True)
 else:finalize(m,groups,round)

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
 validation=dict(generated_at=r.now(),passed=True,conditions=len(cases),readouts=len(rows),distinct_simulated_conditions=len(groups),distinct_trials=sum(len(v) for v in groups.values()),new_base_trials=sum(len(j['durations']) for j in m['base_jobs']),new_refinement_trials=sum(len(j['durations']) for k in range(round) for j in read(OUT/f'jobs_{k}.json')),coverage_complete=True,exact_checkpoint_hashes_verified=True,source_hashes_verified=True,cue_checks_passed=True,full_gate_implies_original_gate=True,selected_trials_verified=len(selected_checks),selected_native_arrays_all_finite=True,status_counts_by_definition={name:dict(Counter(x['status'] for x in rows if x['definition']==name)) for name in ['orientation','full_polarity']},nonmonotonic_readouts=sum(q.get('nonmonotonic_success_pattern',False) for q in rows),successful_brackets_all_within_one_percent=True,all_gate_samples_checked=True)
 write(RESULTS/'validation.json',validation);write(RESULTS/'selected_trace_checks.json',selected_checks)
 status=dict(validation,complete=True,producer='finished',stage='complete')
 from countercue_full_lengths_plot import plot
 plot(cases,status)
 write(RESULTS/'status.json',status)
 write(RESULTS/'artifact_hashes.json',{f.name:r.sha(f) for f in sorted(RESULTS.iterdir()) if f.suffix in ['.png','.pdf','.svg','.csv']})
 print(json.dumps(validation),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['launch','worker','plan','progress']);ap.add_argument('--slot',type=int);ap.add_argument('--round',type=int,default=-1);ap.add_argument('--policy',choices=['allow','restrict']);args=ap.parse_args()
 if args.action=='launch':
  assert args.policy is not None;launch(args.policy)
 elif args.action=='worker':worker(args.slot,args.round)
 elif args.action=='plan':plan(args.round)
 else:print(json.dumps(progress('base_trials_running' if args.round==0 else 'threshold_refinements_running',round=args.round)))
