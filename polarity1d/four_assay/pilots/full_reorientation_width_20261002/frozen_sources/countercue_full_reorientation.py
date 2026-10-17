"""Companion width analysis: sustained polarity within 1% of reflected initial strength."""
from pathlib import Path
import argparse,gzip,json,math,os,shlex,subprocess,sys
from collections import Counter
import numpy as np
import countercue_pilot as p
import runner as r

CODE=Path(__file__).resolve().parent
OUT=CODE/'pilots/full_reorientation_width_20261002'
RESULTS=OUT/'results';WORK=OUT/'work'
SOURCE=p.RESULTS/'records.json.gz'

def write(path,value):p.write(path,value)
def read(path):return p.read(path)
def verify():
 p.verify();m=read(OUT/'manifest.json')
 assert r.sha(SOURCE)==m['source_records_sha256']
 assert r.sha(Path(__file__).resolve())==m['analysis_source_sha256']
 for c in m['cases']:assert r.sha(r.ROOT/c['checkpoint'])==c['checkpoint_sha256']
 return m

def submit(action,array=None,dependency=None,round=0):
 command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__).resolve()),action,'--round',str(round)])
 if array is not None:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
 opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=02:00:00','--job-name=full-reorient-'+action,'--chdir='+str(r.ROOT),'--output='+str(WORK/'logs'/f'{action}_r{round}_%A_%a.out'),'--error='+str(WORK/'logs'/f'{action}_r{round}_%A_%a.err')]
 if array is not None:opts+=['--array='+array]
 if dependency:opts+=['--dependency=afterany:'+dependency]
 return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]

def status(stage,**kw):
 stamp=r.now();write(RESULTS/'status.json',dict(generated_at=stamp,complete=False,producer='running',stage=stage,**kw))
 (RESULTS/'index.html').write_text(f'<html><body><h1>Width companion: full polarity reorientation</h1><p>Generated {stamp}; producer running: {stage}.</p><p>Final 500-unit recovered polarity must remain within 1% of reversed initial strength; original stability and contrast gates retained. Recovery horizon 10000.</p></body></html>')

def score(t,tail,initial):
 tail=np.asarray(tail);assert len(tail)==201 and np.all(np.isfinite(tail));assert initial<-.05
 relative=tail/abs(initial)
 error=float(np.max(abs(relative-1)))
 stable=bool(np.ptp(tail)<.005 and t['minimum_contrast_last_500']>=.05)
 assert stable==t['settled']
 return dict(**t,full_success=bool(t['success'] and error<=.01),initial_orientation=initial,maximum_relative_error_last_500=error,minimum_normalized_polarity_last_500=float(np.min(relative)),maximum_normalized_polarity_last_500=float(np.max(relative)))

def launch():
 assert not OUT.exists()
 for d in [RESULTS,WORK/'logs',WORK/'scans',WORK/'chunks',OUT/'frozen_sources']:d.mkdir(parents=True,exist_ok=True)
 p.verify();source=json.load(gzip.open(SOURCE,'rt'))
 cases=[{k:v for k,v in c.items() if k not in ['tests','result']} for c in source['cases'] if c['family']=='width']
 for c in cases:
  _,tr=p.state_meta(c);c['initial_orientation']=float(tr.orientation[-1])
  c['old_result']=next(x['result'] for x in source['cases'] if x['case_index']==c['case_index'])
 manifest=dict(created_at=r.now(),cases=cases,source_records_sha256=r.sha(SOURCE),analysis_source_sha256=r.sha(Path(__file__).resolve()),definition='Every original 2.5-unit sample in final 500 recovery units must be within [0.99,1.01] of abs(initial settled orientation), with original success gate retained.',recovery_horizon=10000,bracket_relative_resolution=.01,maximum_duration='Original pilot duration maximum, no horizon extension',nonmonotonic='All available tested durations rescored; shortest passing tested pulse, no global monotonicity assumption')
 write(OUT/'manifest.json',manifest)
 import shutil
 shutil.copy2(Path(__file__).resolve(),OUT/'frozen_sources'/Path(__file__).name)
 status('rescoring_saved_width_trials')
 array=submit('scan',array='0-6');controller=submit('plan',dependency=array)
 write(OUT/'launch.json',dict(scan_array=array,controller=controller,stage='rescoring_saved_width_trials'))
 print(json.dumps(read(OUT/'launch.json')),flush=True)

def scan(slot):
 m=verify();model=r.MODELS[slot];source=json.load(gzip.open(SOURCE,'rt'));lookup={c['case_index']:c for c in source['cases']}
 cs=[c for c in m['cases'] if c['model']==model and c['case_index']==c['canonical_index']]
 groups={}
 for c in cs:
  original=lookup[c['case_index']];bytrace={}
  for t in original['tests']:bytrace.setdefault(t['trace'],[]).append(t)
  tests=[]
  for trace,rows in bytrace.items():
   with np.load(r.ROOT/trace) as z:
    for t in rows:
     initial=float(z[f"trial{t['trace_row']}_pulse_orientation"][0])
     np.testing.assert_allclose(initial,c['initial_orientation'],rtol=1e-8,atol=1e-10)
     tests.append(score(t,z[f"trial{t['trace_row']}_recovery_tail_orientation"],initial))
  assert len(tests)==len(original['tests'])
  groups[str(c['case_index'])]=sorted(tests,key=lambda t:t['duration'])
 write(WORK/'scans'/f'{model}.json',dict(model=model,groups=groups,complete=True))
 print(json.dumps(dict(model=model,cases=len(groups),trials=sum(len(t) for t in groups.values()))),flush=True)

def collected(round):
 groups={}
 for model in r.MODELS:
  part=read(WORK/'scans'/f'{model}.json');assert part['complete'];groups.update({int(k):v for k,v in part['groups'].items()})
 for previous in range(round):
  jobs=read(OUT/f'jobs_{previous}.json')
  for job in jobs:
   path=WORK/'chunks'/f'r{previous}_{job["slot"]:05d}.json'
   row=read(path);assert row['case_index']==job['case_index']
   assert [t['duration'] for t in row['tests']]==job['durations']
   groups[job['case_index']].extend(row['tests'])
 for tests in groups.values():
  tests.sort(key=lambda t:t['duration']);ds=[t['duration'] for t in tests];assert len(ds)==len(set(ds))
 return groups

def plan(round):
 m=verify();groups=collected(round);jobs=[];targets=[]
 for c in m['cases']:
  if c['case_index']!=c['canonical_index']:continue
  tests=groups[c['case_index']];passing=[t for t in tests if t['full_success']]
  if not passing:continue
  hi=passing[0]['duration'];lower=[t['duration'] for t in tests if t['duration']<hi];lo=max(lower,default=0.)
  if lo>0 and (hi-lo)/hi<=.01:continue
  if lo==0:ds=[hi/2**k for k in range(1,11)]
  else:
   steps=math.ceil((hi-lo)/(.01*lo));ds=np.linspace(lo,hi,steps+1)[1:-1].tolist()
  existing={t['duration'] for t in tests};ds=sorted(set(float(d) for d in ds)-existing)
  for i in range(0,len(ds),10):jobs.append(dict(slot=len(jobs),case_index=c['case_index'],durations=ds[i:i+10]))
  targets.append(dict(case_index=c['case_index'],model=c['model'],length=c['length'],bracket=[lo,hi],additional_trials=len(ds)))
 if jobs and round<4:
  write(OUT/f'jobs_{round}.json',jobs);write(OUT/f'targets_{round}.json',targets)
  status('refining_full_reorientation_thresholds',round=round,additional_jobs=len(jobs),additional_trials=sum(len(j['durations']) for j in jobs))
  array=submit('worker',array=f'0-{len(jobs)-1}',round=round);controller=submit('plan',dependency=array,round=round+1)
  launch=read(OUT/'launch.json');launch.update(stage='refining_full_reorientation_thresholds',round=round,array=array,controller=controller);write(OUT/'launch.json',launch)
  print(json.dumps(launch),flush=True)
 else:finalize(m,groups,round)

def worker(slot,round):
 m=verify();job=read(OUT/f'jobs_{round}.json')[slot];c=next(c for c in m['cases'] if c['case_index']==job['case_index']);rows=[];arrays={}
 path=WORK/'chunks'/f'r{round}_{slot:05d}'
 for i,D in enumerate(job['durations']):
  t,z=p.trial(c,D);t=score(t,z['recovery_tail_orientation'],c['initial_orientation']);t.update(trace_row=i,trace=str(path.with_suffix('.npz').relative_to(r.ROOT)))
  rows.append(t);arrays.update({f'trial{i}_{k}':v for k,v in z.items()})
 np.savez_compressed(path.with_suffix('.npz'),**arrays)
 verify();write(path.with_suffix('.json'),dict(case_index=c['case_index'],tests=rows,generated_at=r.now()))


def finalize(m,groups,round):
 rows=[];cases=[]
 for c in m['cases']:
  tests=groups[c['canonical_index']];passing=[t for t in tests if t['full_success']];selected=passing[0] if passing else tests[-1]
  endpoint=selected['duration'] if passing else None;lo=max((t['duration'] for t in tests if endpoint is not None and t['duration']<endpoint),default=0.)
  row=dict(case_index=c['case_index'],canonical_index=c['canonical_index'],model=c['model'],length=c['length'],condition=c['condition'],fraction=c['fraction'],width_mode=c['width_mode'],dose_mode=c['dose_mode'],status='success' if passing else 'censored',endpoint=endpoint,last_failing_duration_below_endpoint=lo if passing else None,relative_bracket_width=(endpoint-lo)/endpoint if passing else None,lower_grid_limited=bool(passing and lo==0),maximum_tested_duration=max(t['duration'] for t in tests),tested_count=len(tests),original_gate_endpoint=c['old_result']['endpoint'],original_gate_status=c['old_result']['status'],nonmonotonic_success_pattern=any(tests[i]['full_success'] and not tests[i+1]['full_success'] for i in range(len(tests)-1)),selected_duration=selected['duration'],selected_trace=selected['trace'],selected_trace_row=selected['trace_row'],selected_morphology=selected['morphology'],selected_tail_field_pass_0_1_percent=selected['tail_field_pass_0_1_percent'],selected_maximum_relative_error_last_500=selected['maximum_relative_error_last_500'],initial_orientation=c['initial_orientation'])
  rows.append(row);cases.append(dict(**c,result=row,tests=tests))
 with gzip.open(RESULTS/'records.json.gz','wt') as f:json.dump(dict(generated_at=r.now(),cases=cases,manifest=m),f,default=r.convert,allow_nan=False)
 r.csvwrite(RESULTS/'endpoints.csv',rows)
 assert all(t['success'] for tests in groups.values() for t in tests if t['full_success'])
 additional=sum(len(j['durations']) for k in range(round) for j in read(OUT/f'jobs_{k}.json'))
 validation=dict(generated_at=r.now(),passed=True,conditions=len(rows),distinct_conditions=len(groups),rescored_saved_trials=sum(len(t) for t in groups.values())-additional,additional_trials=additional,status_counts=dict(Counter(x['status'] for x in rows)),lower_grid_limited=sum(x['lower_grid_limited'] for x in rows),brackets_wider_than_one_percent=sum(x['status']=='success' and x['relative_bracket_width']>.01000000001 for x in rows),nonmonotonic_conditions=sum(x['nonmonotonic_success_pattern'] for x in rows),source_data_unchanged=r.sha(SOURCE)==m['source_records_sha256'],all_201_tail_samples_checked=True,exact_release_checkpoints_verified=True,original_gate_preserved=True)
 write(RESULTS/'validation.json',validation)
 plot(cases,validation)
 write(RESULTS/'status.json',dict(generated_at=r.now(),complete=True,producer='finished',**validation))
 print(json.dumps(validation),flush=True)


def plot(cases,validation):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.lines import Line2D
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 labels=['Wave-pinning','Otsuji','Goryachev','LEGI','deBelly','Holmes 3','Spring']
 fig,axes=plt.subplots(7,2,figsize=(14,20))
 for i,(model,label) in enumerate(zip(r.MODELS,labels)):
  for j,dose in enumerate(['peak_matched','integral_matched']):
   ax=axes[i,j]
   for fraction,color in [(.1,'#176d9c'),(.2,'#db7b24')]:
    for mode,ls in [('fractional','--'),('physical','-')]:
     chosen=sorted([c for c in cases if c['model']==model and c['fraction']==fraction and c['width_mode']==mode and c['dose_mode']==('peak_matched' if mode=='fractional' else dose)],key=lambda c:c['length'])
     x=[c['length'] for c in chosen];y=[c['result']['endpoint'] if c['result']['status']=='success' else np.nan for c in chosen]
     ax.plot(x,y,color=color,ls=ls,lw=1.7,label=f'{mode}, {100*fraction:g}%')
     for c in chosen:
      q=c['result']
      if q['status']=='censored':ax.plot(c['length'],q['maximum_tested_duration'],'^',mfc='white',mec=color,ms=7,mew=1.2)
      else:
       marker='s' if q['selected_morphology']=='nonmonotone' else 'o'
       ax.plot(c['length'],q['endpoint'],marker,color=color,mfc=color if q['selected_tail_field_pass_0_1_percent'] else 'white',ms=5)
   ax.set_yscale('log');ax.set_xticks(p.LENGTHS[model]);ax.set_ylabel('Shortest qualifying tested pulse\n(native time units)');ax.set_xlabel('Resting L0' if model=='spring' else 'Nominal cell length')
   ax.set_title(label+(' — fixed peak A' if j==0 else ' — equal spatial integral'))
   ax.grid(axis='y',color='#eeeeee',lw=.6)
   if i==0:ax.legend(ncol=2,fontsize=8,loc='upper left')
 fig.suptitle('Finite countercue: full polarity reorientation across lengths',fontsize=15,y=.985)
 fig.text(.5,.958,'Recovered polarity / |initial settled polarity| must remain within 0.99–1.01 throughout the final 500 units.\nOriginal orientation stability and contrast gates retained; recovery lasts 10000 units.',ha='center',fontsize=10)
 fig.text(.5,.027,'Hollow upward triangle: no qualifying tested pulse. Square: retained nonmonotone active profile.\nHollow circle/square: native fields still vary >0.1% fixed scale in the tail. Missing endpoints break lines.\nIntegral matching is between physical and fractional patches at each length. Polarity recovery does not imply complete pattern replacement.\nGenerated '+validation['generated_at']+'; producer finished. Original width thresholds preserved.',ha='center',fontsize=9)
 fig.subplots_adjust(left=.11,right=.98,top=.915,bottom=.105,hspace=.60,wspace=.28)
 for ext in ['png','pdf','svg']:fig.savefig(RESULTS/f'width_full_reorientation.{ext}',dpi=180)
 plt.close(fig)
 (RESULTS/'index.html').write_text(f'<html><head><meta charset="utf-8"><title>Full polarity reorientation</title><style>body{{font-family:sans-serif;max-width:1400px;margin:30px auto}}img{{width:100%}}</style></head><body><h1>Width companion: full polarity reorientation</h1><p>Generated {validation["generated_at"]}; producer finished. {validation["status_counts"]}</p><p>All 201 original samples in the final 500 recovery units must have polarity within 1% of the reflected initial magnitude. Original gates and 10000 recovery horizon retained. This is polarity recovery, not complete spatial pattern replacement.</p><p><a href="endpoints.csv">Endpoints and original threshold comparison</a> · <a href="validation.json">Validation</a> · <a href="../manifest.json">Definition and provenance</a></p><a href="width_full_reorientation.pdf">PDF</a> · <a href="width_full_reorientation.svg">SVG</a><img src="width_full_reorientation.png"></body></html>')
 write(RESULTS/'artifact_hashes.json',{f.name:r.sha(f) for f in sorted(RESULTS.iterdir()) if f.suffix in ['.png','.pdf','.svg','.csv']})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['launch','scan','plan','worker']);ap.add_argument('--slot',type=int);ap.add_argument('--round',type=int,default=0);args=ap.parse_args()
 if args.action=='launch':launch()
 elif args.action=='scan':scan(args.slot)
 elif args.action=='worker':worker(args.slot,args.round)
 else:plan(args.round)
