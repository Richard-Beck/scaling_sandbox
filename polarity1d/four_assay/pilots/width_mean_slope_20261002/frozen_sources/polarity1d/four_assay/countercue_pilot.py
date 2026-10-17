"""Finite-countercue width and independent mean/slope pilot, exact v11 release states."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,os,shlex,subprocess,sys
from collections import Counter
import numpy as np
import runner as r
import settling
CODE=Path(__file__).resolve().parent
OUT=CODE/'pilots/width_mean_slope_20261002'
WORK=OUT/'work'; RESULTS=OUT/'results'
LENGTHS={'wave_pinning':[9.,21.,39.], 'otsuji':[9.,21.,39.], 'goryachev':[9.,21.,39.], 'legi':[9.,17.,31.], 'debelly':[19.,27.,39.], 'holmes_model3':[14.,30.,50.], 'spring':[.8,1.6,2.6]}
BASE=CODE/'results/records.json.gz'
PLOT=CODE/'countercue_pilot_plot.py'

def write(p,v):p.write_text(json.dumps(v,indent=2,default=r.convert,allow_nan=False)+'\n')
def read(p):return json.loads(p.read_text())
def records():return json.load(gzip.open(BASE,'rt'))
def source_hashes():
 h=r.sources()
 for p in [Path(__file__).resolve(),PLOT,BASE]:h[str(p.relative_to(r.ROOT))]=r.sha(p)
 return h

def verify():
 frozen=read(OUT/'provenance.json')
 assert source_hashes()==frozen['sources'],'Campaign source changed'
 assert r.sha(OUT/'manifest.json')==frozen['manifest_sha256']
 return frozen

class Cue:
 def __init__(self,case,duration):
  self.case=case;self.stop=float(duration);self.event_times=(self.stop,)
 def __call__(self,t,x,L):
  x=np.asarray(x,float);c=self.case
  if t>=self.stop:return np.zeros_like(x)
  if c['family']=='mean_slope':
   if L>c['length_bound']+1e-8:raise ValueError('Physical length exceeds linear-cue positivity bound')
   y=c['A']*(c['mean_factor']+c['slope_factor']*(x-L/2)/c['length_bound'])
  else:
   w=c['fraction']*(L if c['width_mode']=='fractional' else c['reference_domain_length'])
   B=c['A']*(L/c['reference_domain_length'] if c['dose_mode']=='integral_matched' else 1.)
   if w>L/2+1e-8:raise ValueError('Patch exceeds half the cell domain')
   z=(L-x)/w
   y=B*.5*(1+np.cos(np.pi*np.minimum(np.maximum(z,0),1)))
   y=np.where((z>=0)&(z<=1),y,0.)
  if np.min(y)<-1e-12:raise ValueError('Negative cue; clipping is prohibited')
  return y

def submit(action,array=None,dependency=None,time='02:00:00'):
 cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__).resolve()),action])
 if array:cmd+=' --slot "$SLURM_ARRAY_TASK_ID"'
 opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time='+time,'--job-name=cue-pilot-'+action,'--chdir='+str(r.ROOT),'--output='+str(WORK/'logs'/f'{action}_%A_%a.out'),'--error='+str(WORK/'logs'/f'{action}_%A_%a.err')]
 if array:opts.append('--array='+array)
 if dependency:opts.append('--dependency=afterany:'+dependency)
 return subprocess.check_output(opts+['--wrap',cmd],text=True).strip().split(';')[0]

def state_meta(case):
 p=r.ROOT/case['checkpoint'];assert r.sha(p)==case['checkpoint_sha256']
 state,tr=r.read_checkpoint(p,case['model'])
 return state,tr

def status(stage,**extra):
 stamp=r.now();write(RESULTS/'status.json',dict(generated_at=stamp,complete=False,producer='running',stage=stage,**extra))
 (RESULTS/'index.html').write_text(f'<html><body><h1>Finite-countercue pilot</h1><p>Generated {stamp}; producer running: {stage}. Exact original release states. Final plots pending.</p></body></html>')

def launch():
 OUT.mkdir(exist_ok=False)
 for p in [WORK/'chunks',WORK/'logs',RESULTS,OUT/'frozen_sources']:p.mkdir(parents=True)
 old=records();preps={(p['model'],p['length']):p for p in old['preparations'] if p['assay']=='finite_countercue'}
 tasks={(t['model'],t['length']):t for t in r.tasks() if t['assay']=='finite_countercue'}
 scale_rows=read(CODE/'reference/continuation_20261001/fixed_scales.json')
 cases=[]
 for model,lengths in LENGTHS.items():
  A=r.PROTOCOL['normalization']['amplitude_by_model'][model];Lref=lengths[1]
  ref=preps[model,Lref];_,tr=r.read_checkpoint(r.ROOT/ref['checkpoint'],model)
  ref_domain=float(tr.fields['length'][-1]) if model=='spring' else Lref
  # For Spring, d(length)/dt=-.02*(length-L0)+Fleft+Fright, each F<=.001.
  # Hence a release extension <=.1 remains <=.1 for all cues.
  Lbound=lengths[-1]+.1 if model=='spring' else lengths[-1]
  for L in lengths:
   prep=preps[model,L];assert prep['eligible']
   assert r.sha(r.ROOT/prep['checkpoint'])==prep['checkpoint_sha256']
   fixed=next(z['fixed_scales'] for z in scale_rows if z['model']==model and z['length']==L and z['event']=='release')
   base=dict(model=model,length=L,cells=tasks[model,L]['cells'],A=A,reference_length=Lref,reference_domain_length=ref_domain,length_bound=Lbound,checkpoint=prep['checkpoint'],checkpoint_sha256=prep['checkpoint_sha256'],fixed_scales=fixed)
   for f in [.1,.2]:
    for mode,dose in [('fractional','peak_matched'),('physical','peak_matched'),('physical','integral_matched')]:
     cases.append(dict(**base,case_index=len(cases),family='width',fraction=f,width_mode=mode,dose_mode=dose,condition=f'{mode}_{dose}_f{f:g}'))
   for mean in [.5,1.,2.]:
    for slope in [.25,.5,1.]:
     cases.append(dict(**base,case_index=len(cases),family='mean_slope',mean_factor=mean,slope_factor=slope,condition=f'mean{mean:g}_slope{slope:g}'))
 # Deduplicate cues that are exactly identical at the static reference length.
 for c in cases:
  canonical=c['case_index']
  if c['family']=='width' and c['model']!='spring' and c['length']==c['reference_length']:
   canonical=next(d['case_index'] for d in cases if d['model']==c['model'] and d['length']==c['length'] and d['family']=='width' and d['fraction']==c['fraction'] and d['width_mode']=='fractional')
  c['canonical_index']=canonical
 jobs=[]
 for c in cases:
  if c['case_index']!=c['canonical_index']:continue
  D=r.PROTOCOL['duration_grids'][c['model']]['refined_durations']
  for i in range(0,len(D),10):jobs.append(dict(slot=len(jobs),case_index=c['case_index'],durations=D[i:i+10]))
 manifest=dict(created_at=r.now(),lengths=LENGTHS,cases=cases,base_jobs=jobs,recovery=10000.,gate_window=500.,threshold_definition='Shortest passing tested pulse; original gate, with explicit lower-grid refinement',refinement_rule=dict(floor_halvings=10,bracket_points=9),spring_length_bound='L0+.1 follows bounded native protrusion force; runtime checked',cue_coordinates='Native cue coordinates; Spring instantaneous physical domain, deBelly original reference domain',histories='All native fields sampled at 21 pulse points, logarithmic early-recovery points and 21 tail points; dense original solver samples used for tail diagnostics')
 # Check analytic cues, positivity, integral matching and aliases before submission.
 checks=0
 for c in cases:
  L=c['length'] if c['model']!='spring' else c['length']+.05
  x=np.linspace(0,L,20001);q=Cue(c,10)(0,x,L);assert np.min(q)>=-1e-12
  if c['family']=='mean_slope':
   np.testing.assert_allclose(np.trapezoid(q,x)/L,c['A']*c['mean_factor'],rtol=1e-12)
  else:
   w=c['fraction']*(L if c['width_mode']=='fractional' else c['reference_domain_length'])
   peak=c['A']*(L/c['reference_domain_length'] if c['dose_mode']=='integral_matched' else 1.)
   np.testing.assert_allclose(np.trapezoid(q,x),peak*w/2,rtol=1e-7)
  assert np.all(Cue(c,10)(10,x,L)==0)
  if c['canonical_index']!=c['case_index']:
   np.testing.assert_array_equal(q,Cue(cases[c['canonical_index']],10)(0,x,L))
  checks+=1
 write(OUT/'manifest.json',manifest)
 hashes=source_hashes()
 for rel in hashes:
  src=r.ROOT/rel;dest=OUT/'frozen_sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(src.read_bytes())
 write(OUT/'provenance.json',dict(generated_at=r.now(),sources=hashes,manifest_sha256=r.sha(OUT/'manifest.json'),cue_checks=checks,initial_preparations='Exact original v11 cue-free release checkpoints',original_results_sha256=r.sha(BASE)))
 status('base_trials_running',conditions=len(cases),distinct_conditions=sum(c['case_index']==c['canonical_index'] for c in cases),base_jobs=len(jobs))
 array=submit('worker',array=f'0-{len(jobs)-1}')
 controller=submit('refine_plan',dependency=array)
 write(OUT/'launch.json',dict(generated_at=r.now(),stage='base_trials_running',base_array=array,refinement_controller=controller,base_jobs=len(jobs),conditions=len(cases)))
 print(json.dumps(read(OUT/'launch.json')),flush=True)

def spatial_integral(y,L,model):
 return float(np.mean(y)*L) if model!='spring' else float(np.trapezoid(y,np.linspace(0,L,len(y))))

def samples(tr,requested):
 idx=np.unique(np.clip(np.searchsorted(tr.time,np.asarray(requested)),0,len(tr.time)-1))
 return idx

def first_time(tr,good):
 idx=np.flatnonzero(good);return float(tr.time[idx[0]]) if len(idx) else None

def trial(case,D):
 model,L=case['model'],case['length'];state,initial=state_meta(case)
 cue=Cue(case,D)
 pulse=r.run(model,L,state,cue,D,interval=min(2.5,D/20),n=case['cells'])
 recovery=r.run(model,L,pulse.state,r.Cue('zero',model),10000.,interval=2.5,n=case['cells'])
 keep=recovery.time>=9500.;o=recovery.orientation[keep];contrast=recovery.contrast[keep]
 stable=bool(np.ptp(o)<.005 and np.min(contrast)>=.05);success=bool(stable and np.min(o)>.05)
 monitored=settling.monitored(recovery,model,L);tail_metrics={}
 for name,y in monitored.items():
  scale=case['fixed_scales'][name];z=y[keep]
  variation=float(np.max(np.ptp(z,axis=0))/scale)
  rate=float(np.max(abs(z[-1]-z[-2]))/(recovery.time[-1]-recovery.time[-2])*500/scale)
  tail_metrics[name]=dict(range_scaled=variation,projected_drift_500=rate)
 fieldtail=max(max(v.values()) for v in tail_metrics.values())
 scale=case['fixed_scales'][{'legi':'response','holmes_model3':'rac_active','debelly':'rac'}.get(model,'active')]
 y=recovery.active[-1];aligned=y if o[-1]>0 else y[::-1]
 drawdown=float(np.max(np.maximum.accumulate(aligned)-aligned)/scale)
 start_right=float(initial.right[-1]);start_left=float(initial.left[-1]);field_contrast=float(initial.contrast[-1])
 # Cue-driven rates use the original reaction law and the sampled instantaneous fields.
 pidx=samples(pulse,np.linspace(0,D,21));times=pulse.time[pidx]
 activation=[];cue_integrals=[];domain_lengths=[]
 for i in pidx:
  physical=float(pulse.fields['length'][i]) if model=='spring' else L
  x=pulse.x if model!='spring' else np.linspace(0,physical,case['cells'])
  if model=='debelly' and r.model_for(model).parameters.input_sampling_mode=='julia_element_right_edge':x=np.minimum(x+.5*L/case['cells'],L)
  # Sample the pulse's left limit at withdrawal for temporal quadrature.
  profile=cue(min(float(pulse.time[i]),np.nextafter(D,0)),x,physical)
  if model=='debelly':
   par=r.model_for(model).parameters;rate=par.rac_deactivation*profile/(1+pulse.fields['rho'][i]**2) if par.opto_target=='rac' else np.zeros_like(profile)
  elif model=='spring':rate=profile*pulse.fields['inactive'][i]
  else:
   m=r.model_for(model);s=np.vstack([pulse.fields[name][i] for name in m.state_names]);rate=m.reaction(pulse.time[i],x,s,profile)[0]-m.reaction(pulse.time[i],x,s,np.zeros_like(profile))[0]
  activation.append(spatial_integral(rate,physical,model));cue_integrals.append(spatial_integral(profile,physical,model));domain_lengths.append(physical)
 row=dict(duration=D,success=success,settled=stable,final_orientation=float(o[-1]),orientation_range_last_500=float(np.ptp(o)),minimum_contrast_last_500=float(np.min(contrast)),pulse_end_orientation=float(pulse.orientation[-1]),pulse_end_contrast=float(pulse.contrast[-1]),pulse_right_orientation_time=first_time(pulse,pulse.orientation>.05),new_end_rise_time=first_time(pulse,pulse.right>=start_right+.1*field_contrast),old_end_decay_time=first_time(pulse,pulse.left<=start_left-.1*field_contrast),sampled_cue_space_time_integral=float(np.trapezoid(cue_integrals,times)),sampled_cue_driven_activation_integral=float(np.trapezoid(activation,times)),minimum_pulse_domain_length=min(domain_lengths),maximum_pulse_domain_length=max(domain_lengths),tail_field_max_scaled=fieldtail,tail_field_metrics=tail_metrics,tail_field_pass_0_1_percent=fieldtail<=.001,morphology='approximately_monotone_end' if drawdown<=.01 else 'nonmonotone',morphology_backward_excursion_scaled=drawdown)
 ridx=samples(recovery,np.r_[0.,np.geomspace(2.5,9500,25),np.linspace(9500,10000,21)])
 arrays=dict(pulse_time=pulse.time[pidx],recovery_time=recovery.time[ridx],pulse_orientation=pulse.orientation[pidx],recovery_orientation=recovery.orientation[ridx],recovery_tail_orientation=o,pulse_contrast=pulse.contrast[pidx],recovery_contrast=recovery.contrast[ridx],pulse_active=pulse.active[pidx],recovery_active=recovery.active[ridx],x=pulse.x,cue_integral_times=times,cue_spatial_integral=np.asarray(cue_integrals),cue_driven_activation_rate=np.asarray(activation))
 for name,values in pulse.fields.items():arrays['pulse_field_'+name]=np.asarray(values)[pidx]
 for name,values in recovery.fields.items():arrays['recovery_field_'+name]=np.asarray(values)[ridx]
 for phase,tr in [('pulse',pulse),('recovery',recovery)]:
  if isinstance(tr.state,dict):
   for name,value in tr.state.items():arrays[phase+'_state_'+name]=value
  else:arrays[phase+'_state']=tr.state
 return row,arrays

def worker(slot,refined=False):
 verify();manifest=read(OUT/'manifest.json');jobs=read(OUT/'refinement_jobs.json') if refined else manifest['base_jobs'];job=jobs[slot];case=manifest['cases'][job['case_index']]
 prefix='refine' if refined else 'base';path=WORK/'chunks'/f'{prefix}_{slot:05d}'
 rows=[];arrays={}
 for i,D in enumerate(job['durations']):
  row,a=trial(case,float(D));row['trace_row']=i;rows.append(row)
  arrays.update({f'trial{i}_{k}':v for k,v in a.items()})
 np.savez_compressed(path.with_suffix('.npz'),**arrays)
 verify();write(path.with_suffix('.json'),dict(case_index=case['case_index'],slot=slot,phase=prefix,tests=rows,trace=str(path.with_suffix('.npz').relative_to(r.ROOT)),generated_at=r.now()))
 print(json.dumps(dict(case_index=case['case_index'],slot=slot,phase=prefix,successes=sum(t['success'] for t in rows))),flush=True)

def collected(manifest,refined=False):
 groups={}
 jobs=read(OUT/'refinement_jobs.json') if refined else manifest['base_jobs'];prefix='refine' if refined else 'base'
 for j in jobs:
  p=WORK/'chunks'/f'{prefix}_{j["slot"]:05d}.json'
  if not p.exists():raise RuntimeError(f'Missing execution output: {p}')
  row=read(p);assert row['case_index']==j['case_index'];assert [t['duration'] for t in row['tests']]==[float(d) for d in j['durations']]
  for t in row['tests']:
   t=dict(t,trace=row['trace']);groups.setdefault(j['case_index'],[]).append(t)
 return groups

def refinement_plan():
 verify();m=read(OUT/'manifest.json');groups=collected(m);jobs=[]
 for idx,tests in groups.items():
  ordered=sorted(tests,key=lambda t:t['duration']);passing=[t for t in ordered if t['success']]
  if not passing:continue
  first=passing[0]['duration'];prior=[t['duration'] for t in ordered if t['duration']<first]
  if not prior:
   durations=[first/(2**k) for k in range(1,11)]
  else:
   lower=max(prior);durations=np.linspace(lower,first,11)[1:-1].tolist()
  durations=sorted(set(float(d) for d in durations)-{t['duration'] for t in ordered})
  if durations:jobs.append(dict(slot=len(jobs),case_index=idx,durations=durations))
 write(OUT/'refinement_jobs.json',jobs);status('threshold_refinements_running',base_complete=True,refinement_jobs=len(jobs))
 launch=read(OUT/'launch.json')
 if jobs:
  array=submit('refine_worker',array=f'0-{len(jobs)-1}');final=submit('finalize',dependency=array,time='00:30:00');launch.update(stage='threshold_refinements_running',refinement_array=array,finalize_job=final,refinement_jobs=len(jobs));write(OUT/'launch.json',launch)
 else:finalize()
 print(json.dumps(launch),flush=True)

def finalize():
 verify();m=read(OUT/'manifest.json');groups=collected(m)
 if (OUT/'refinement_jobs.json').exists():
  for idx,tests in collected(m,refined=True).items():groups[idx].extend(tests)
 old=records();baseline={(c['model'],c['length']):c for c in old['cases'] if c['assay']=='finite_countercue'}
 cases=[];rows=[]
 for c in m['cases']:
  tests=sorted(groups[c['canonical_index']],key=lambda t:t['duration']);ds=[t['duration'] for t in tests];assert len(ds)==len(set(ds))
  passed=[t for t in tests if t['success']];endpoint=passed[0]['duration'] if passed else None
  bracket=max((t['duration'] for t in tests if endpoint is not None and t['duration']<endpoint),default=0.)
  selected=passed[0] if passed else tests[-1]
  row=dict(case_index=c['case_index'],canonical_index=c['canonical_index'],model=c['model'],length=c['length'],family=c['family'],condition=c['condition'],fraction=c.get('fraction'),width_mode=c.get('width_mode'),dose_mode=c.get('dose_mode'),mean_factor=c.get('mean_factor'),slope_factor=c.get('slope_factor'),status='success' if passed else 'censored',endpoint=endpoint,last_failing_duration_below_endpoint=bracket if passed else None,lower_grid_limited=bool(passed and endpoint==ds[0]),maximum_tested_duration=max(ds),tested_count=len(tests),nonmonotonic_success_pattern=any(tests[i]['success'] and not tests[i+1]['success'] for i in range(len(tests)-1)),selected_duration=selected['duration'],selected_trace=selected['trace'],selected_trace_row=selected['trace_row'],selected_tail_field_max_scaled=selected['tail_field_max_scaled'],selected_tail_field_pass_0_1_percent=selected['tail_field_pass_0_1_percent'],selected_morphology=selected['morphology'],selected_final_orientation=selected['final_orientation'],selected_activation_integral=selected['sampled_cue_driven_activation_integral'],original_endpoint=baseline[c['model'],c['length']].get('endpoint'),original_status=baseline[c['model'],c['length']]['status'])
  rows.append(row);cases.append(dict(**c,result=row,tests=tests))
 with gzip.open(RESULTS/'records.json.gz','wt') as f:json.dump(dict(generated_at=r.now(),cases=cases,manifest=m,provenance=read(OUT/'provenance.json')),f,default=r.convert,allow_nan=False)
 r.csvwrite(RESULTS/'endpoints.csv',rows)
 # Exact controls at the longest static length: mean A/2 and slope A/Lmax.
 controls=[]
 for c in cases:
  if c['family']=='mean_slope' and c['mean_factor']==.5 and c['slope_factor']==1 and c['length']==c['length_bound'] and c['model']!='spring':
   original={t['duration']:t for t in baseline[c['model'],c['length']]['tests']}
   mismatches=[t['duration'] for t in c['tests'] if t['duration'] in original and t['success']!=original[t['duration']]['success']]
   assert not mismatches,(c['model'],mismatches)
   controls.append(dict(model=c['model'],length=c['length'],passing_pattern_matches=True,original_grid_trials_checked=len(original)))
 validation=dict(generated_at=r.now(),passed=True,conditions=len(cases),distinct_conditions=len(groups),base_trial_count=sum(len(j['durations']) for j in m['base_jobs']),refinement_trial_count=sum(len(j['durations']) for j in read(OUT/'refinement_jobs.json')),complete_coverage=True,exact_checkpoint_hashes_verified=True,source_hashes_verified=True,original_cue_controls=controls,aliases_reused_exactly=True,full_final_states_and_sampled_native_field_histories_saved=True,status_counts=dict(Counter(t['status'] for t in rows)),lower_grid_limited=sum(t['lower_grid_limited'] for t in rows),gate_success_with_field_tail_above_0_1_percent=sum(t['status']=='success' and not t['selected_tail_field_pass_0_1_percent'] for t in rows),nonmonotonic_conditions=sum(t['nonmonotonic_success_pattern'] for t in rows))
 write(RESULTS/'validation.json',validation)
 subprocess.run([sys.executable,str(PLOT)],check=True)
 write(RESULTS/'status.json',dict(generated_at=r.now(),complete=True,producer='finished',validation_passed=True,**{k:validation[k] for k in ['conditions','distinct_conditions','base_trial_count','refinement_trial_count','status_counts','lower_grid_limited']}))
 launch=read(OUT/'launch.json');launch.update(stage='finished',finished_at=r.now());write(OUT/'launch.json',launch)
 print(json.dumps(validation),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['launch','worker','refine_plan','refine_worker','finalize','status']);ap.add_argument('--slot',type=int);args=ap.parse_args()
 if args.action=='launch':launch()
 elif args.action=='worker':worker(args.slot)
 elif args.action=='refine_worker':worker(args.slot,refined=True)
 elif args.action=='refine_plan':refinement_plan()
 elif args.action=='finalize':finalize()
 else:print((OUT/'launch.json').read_text())
