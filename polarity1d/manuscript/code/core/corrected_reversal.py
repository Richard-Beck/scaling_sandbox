"""Matched five-model countercue from exact corrected basal-establishment endpoints.

No legacy preparation or trial reuse. Sparse field audit statistics keep trial
storage compact; all scalar pulse/recovery samples remain available.
"""
from __future__ import annotations
import argparse,gzip,json,math,os,shlex,shutil,subprocess,sys
from pathlib import Path
from collections import Counter
import numpy as np
import native as n
import polarity_establishment as e
import settling

MODELS=['wave_pinning','goryachev','debelly','holmes_model3','spring']
DEFAULT=n.DATA/'runs/corrected_reversal_20261003'
ESTABLISHMENT=n.DATA/'establishment_basal'
GRID=n.DATA/'inputs/reversal_duration_grids.json'
DEFINITIONS=['orientation','full_polarity']
RECOVERY=10000.;STEP=2.5;TAIL=500.;MAX_ROUNDS=12


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    n.write(path,value)


def read(path):return json.loads(path.read_text())


def source_hashes():
    sources=n.sources()
    for path in [GRID,ESTABLISHMENT/'manifest.json',ESTABLISHMENT/'results/endpoints.csv',n.ROOT/'code/audit_reversal.py']:
        sources[str(path.relative_to(n.ROOT))]=n.sha(path)
    return sources


def verify(out):
    provenance=read(out/'provenance.json')
    if source_hashes()!=provenance['sources']:
        raise RuntimeError('Frozen source/input hashes changed')
    if n.sha(out/'manifest.json')!=provenance['manifest_sha256']:
        raise RuntimeError('Frozen experiment manifest changed')
    return read(out/'manifest.json')


def initialize(out):
    n.protect_retained_output(out)
    if out.exists():raise ValueError('Use a fresh output directory')
    import csv
    original=read(ESTABLISHMENT/'manifest.json')
    summaries={int(c['case_index']):c for c in csv.DictReader((ESTABLISHMENT/'results/endpoints.csv').open())}
    grids=read(GRID)['durations_by_model'];cases=[]
    for original_case in original['cases']:
        if original_case['model'] not in MODELS or original_case['amplitude_factor']!=1.:continue
        c=dict(original_case);s=summaries[c['case_index']]
        c.update(establishment_case_index=c['case_index'],case_index=len(cases),preparation_status=s['status'],
                 checkpoint=None,checkpoint_sha256=None,eligible=False,initial_orientation=None,initial_contrast=None)
        if s['status']=='polarized':
            checkpoint=ESTABLISHMENT/s['selected_checkpoint']
            if n.sha(checkpoint)!=s['selected_checkpoint_sha256']:raise RuntimeError('Selected input checkpoint hash mismatch')
            _,tr=n.read_checkpoint(checkpoint,c['model'])
            ori=float(tr.orientation[-1]);contrast=float(tr.contrast[-1])
            c.update(checkpoint=str(checkpoint.relative_to(n.ROOT)),checkpoint_sha256=s['selected_checkpoint_sha256'],
                     eligible=bool(ori<-.05 and contrast>=.05 and np.isfinite(ori)),initial_orientation=ori,initial_contrast=contrast)
            if not c['eligible']:c['preparation_status']='checkpoint_failed_polarity_gate'
        c['durations']=grids[c['model']]
        cases.append(c)
    if len(cases)!=330:raise RuntimeError('Expected five-model 330-condition coverage')
    out.mkdir(parents=True)
    for folder in ['work/chunks','work/jobs','work/logs','results','frozen_sources']:(out/folder).mkdir(parents=True,exist_ok=True)
    manifest=dict(created_at=n.now(),models=MODELS,excluded_models=['otsuji','legi'],cases=cases,
                  amplitude_factor=1.,recovery_horizon=RECOVERY,sample_interval=STEP,gate_window=TAIL,
                  original_gate='All tail orientation >0.05, contrast >=0.05, orientation temporal range <0.005.',
                  full_polarity_gate='Original gate plus all tail orientation within 1% of abs(initial orientation).',
                  field_tail_gate='Per-field pointwise tail temporal range and last-sample drift projected over 500, divided by frozen establishment scales, both <=0.001; separate flag, not a scalar response gate.',
                  no_cue_gate='Duration-zero recovery must retain its initial full native profile within 1% of frozen establishment scales, and its final500 samples must stay left-oriented below -0.05, contrast>=0.05 and orientation range<0.005. Exact original input checkpoint remains the start.',
                  cue='Mirror of the establishment cosine patch: same extra peak, width, width mode and background; withdraw extra patch after each independent pulse.',
                  conditions=original['conditions'],reference_lengths=original['reference_lengths'],
                  refinement=dict(relative_bracket_tolerance=.01,maximum_rounds=MAX_ROUNDS,
                                  rule='Refine every adjacent sampled fail-to-pass transition separately for each gate; no global monotonicity assumption; zero-bound transitions extend lower using ten halvings.'),
                  data_policy='All scalar time samples; final native states and pointwise field-tail extrema/last two samples. Full native spatiotemporal histories only for selected post-final reconstructions.')
    write(out/'manifest.json',manifest)
    hashes=source_hashes()
    for rel in hashes:
        src=n.ROOT/rel;dest=out/'frozen_sources'/rel
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    write(out/'provenance.json',dict(generated_at=n.now(),sources=hashes,manifest_sha256=n.sha(out/'manifest.json')))
    # Stage zero audits every eligible preparation with an independent no-cue control.
    jobs=[dict(slot=i,case_index=c['case_index'],durations=[0.]) for i,c in enumerate(c for c in cases if c['eligible'])]
    save_jobs(out,0,jobs)
    write(out/'launch.json',dict(generated_at=n.now(),stage='initialized',round=0,arrays=[]))
    publish(out,'initialized',complete=False)
    return manifest



def refresh_initialized(out):
    """Refresh only an unlaunched initialization after development checks."""
    state=read(out/'launch.json')
    if state['stage']!='initialized' or state.get('arrays') or any((out/'work/chunks').iterdir()):
        raise RuntimeError('Only an unlaunched, trial-free initialization can be refreshed')
    manifest=read(out/'manifest.json');grids=read(GRID)['durations_by_model']
    for c in manifest['cases']:c['durations']=grids[c['model']]
    write(out/'manifest.json',manifest)
    hashes=source_hashes()
    for rel in hashes:
        src=n.ROOT/rel;dest=out/'frozen_sources'/rel
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    write(out/'provenance.json',dict(generated_at=n.now(),sources=hashes,manifest_sha256=n.sha(out/'manifest.json')))
    publish(out,'initialized')


def save_jobs(out,round,jobs):
    path=out/'work/jobs'/f'r{round}.json';write(path,jobs)
    write(path.with_suffix('.sha256.json'),dict(sha256=n.sha(path)))


def jobs(out,round):
    path=out/'work/jobs'/f'r{round}.json'
    if n.sha(path)!=read(path.with_suffix('.sha256.json'))['sha256']:raise RuntimeError('Job coverage changed')
    return read(path)


def trial(case,duration):
    checkpoint=n.ROOT/case['checkpoint']
    if n.sha(checkpoint)!=case['checkpoint_sha256']:raise RuntimeError('Checkpoint changed')
    state,prep=n.read_checkpoint(checkpoint,case['model'])
    pulse=None
    if duration>0:
        # Integrate through the left limit, then restart under basal withdrawal.
        cue=e.Patch(case['amplitude'],case['width'],side='right',width_mode=case['width_mode'],
                    fraction=case['width_fraction_at_reference'],background=case.get('background_stimulus',0.))
        pulse=n.run(case['model'],case['length'],state,cue,duration,interval=min(STEP,duration/20),n=case['cells'])
        state=pulse.state
    recovery=n.run(case['model'],case['length'],state,e.withdrawal(case),RECOVERY,interval=STEP,n=case['cells'])
    tail=recovery.time>=RECOVERY-TAIL;o=recovery.orientation[tail];contrast=recovery.contrast[tail]
    if len(o)!=201 or not all(np.isfinite(v).all() for v in [o,contrast]):raise RuntimeError('Invalid tail samples')
    orientation=bool(np.min(o)>.05 and np.min(contrast)>=.05 and np.ptp(o)<.005)
    error=float(np.max(abs(o/abs(case['initial_orientation'])-1)))
    full=bool(orientation and error<=.01)
    fields=settling.monitored(recovery,case['model'],case['length']);initial_fields=settling.monitored(prep,case['model'],case['length'])
    arrays=dict(x=recovery.x,pulse_time=pulse.time if pulse else np.array([0.]),
                pulse_orientation=pulse.orientation if pulse else np.array([case['initial_orientation']]),
                pulse_contrast=pulse.contrast if pulse else np.array([case['initial_contrast']]),
                recovery_time=recovery.time,recovery_orientation=recovery.orientation,recovery_contrast=recovery.contrast)
    maximum=0.;no_cue_excursion=0.;metrics={}
    for name,y in fields.items():
        if not np.isfinite(y).all():raise RuntimeError('Nonfinite full native field')
        scale=case['fixed_scales'][name];a=y[tail]
        minimum=a.min(axis=0);maximum_field=a.max(axis=0)
        change=float(np.max(maximum_field-minimum))/scale
        drift=float(np.max(abs(a[-1]-a[-2])))/STEP*TAIL/scale
        maximum=max(maximum,change,drift);metrics[name]=dict(range_scaled=change,projected_drift_500=drift)
        for suffix,v in [('min',minimum),('max',maximum_field),('penultimate',a[-2]),('last',a[-1])]:arrays['fieldtail_'+name+'_'+suffix]=v
        if duration==0:
            no_cue_excursion=max(no_cue_excursion,float(np.max(abs(y-initial_fields[name][-1])))/scale)
            arrays['no_cue_'+name+'_min']=y.min(axis=0)
            arrays['no_cue_'+name+'_max']=y.max(axis=0)
            arrays['no_cue_'+name+'_initial']=initial_fields[name][-1]
    if isinstance(recovery.state,dict):arrays.update({'finalstate_'+k:v for k,v in recovery.state.items()})
    else:arrays['finalstate']=recovery.state
    row=dict(duration=float(duration),orientation_success=orientation,full_polarity_success=full,
             final_orientation=float(o[-1]),minimum_contrast_last_500=float(np.min(contrast)),orientation_range_last_500=float(np.ptp(o)),
             maximum_relative_polarity_error_last_500=error,field_tail_pass=bool(maximum<=.001),
             field_tail_max_scaled=maximum,field_tail_metrics=metrics,
             no_cue_maximum_scaled_excursion=no_cue_excursion if duration==0 else None,
             no_cue_left_tail_valid=bool(np.max(o)<-.05 and np.ptp(o)<.005 and np.min(contrast)>=.05) if duration==0 else None,
             no_cue_preparation_valid=bool(no_cue_excursion<=.01 and not orientation and np.max(o)<-.05 and np.ptp(o)<.005 and np.min(contrast)>=.05) if duration==0 else None)
    return row,arrays


def worker(out,round,slot):
    manifest=verify(out);job=jobs(out,round)[slot];case=manifest['cases'][job['case_index']]
    rows=[];arrays={}
    for i,d in enumerate(job['durations']):
        row,z=trial(case,d);prefix=f'trial{i}_';arrays.update({prefix+k:v for k,v in z.items()})
        row.update(trace=f'work/chunks/r{round}_{slot:05d}.npz',trace_prefix=prefix)
        rows.append(row)
    path=out/'work/chunks'/f'r{round}_{slot:05d}'
    tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**arrays);tmp.replace(path.with_suffix('.npz'))
    write(path.with_suffix('.json'),dict(generated_at=n.now(),round=round,slot=slot,case_index=case['case_index'],
          checkpoint_sha256=case['checkpoint_sha256'],trials=rows,trace_sha256=n.sha(path.with_suffix('.npz'))))


def collect(out,partial=False):
    manifest=verify(out);groups={c['case_index']:[] for c in manifest['cases']};pending=[];errors=[]
    for path in sorted((out/'work/jobs').glob('r*.json')):
        if path.name.endswith('.sha256.json'):continue
        round=int(path.stem[1:])
        for job in jobs(out,round):
            result=out/'work/chunks'/f'r{round}_{job["slot"]:05d}.json'
            error=result.with_suffix('.error.json')
            if error.exists():errors.append(str(error.relative_to(out)));continue
            if not result.exists():pending.append(dict(round=round,slot=job['slot']));continue
            data=read(result)
            if data['checkpoint_sha256']!=manifest['cases'][job['case_index']]['checkpoint_sha256']:raise RuntimeError('Worker used a different checkpoint')
            if data['case_index']!=job['case_index'] or [t['duration'] for t in data['trials']]!=job['durations']:raise RuntimeError('Trial coverage differs from frozen job')
            if not result.with_suffix('.npz').exists():raise RuntimeError('Missing trial trace')
            groups[job['case_index']].extend(data['trials'])
    for ts in groups.values():
        ts.sort(key=lambda t:t['duration'])
        if len(ts)!=len({t['duration'] for t in ts}):raise RuntimeError('Duplicate duration per case')
    if not partial and (pending or errors):raise RuntimeError(f'Incomplete execution: {len(pending)} missing, {len(errors)} errors')
    return manifest,groups,pending,errors


def preparation_valid(case,trials):
    if not case['eligible']:return False
    zero=next((t for t in trials if t['duration']==0),None)
    return bool(zero and zero['no_cue_preparation_valid'])


def rows_and_records(manifest,groups,complete):
    rows=[];cases=[]
    for c in manifest['cases']:
        ts=groups[c['case_index']];valid=preparation_valid(c,ts);results={}
        for gate in DEFINITIONS:
            passing=[t for t in ts if t['duration']>0 and t[gate+'_success']]
            best=passing[0] if valid and passing else None
            status='ineligible' if not c['eligible'] else ('preparation_invalid' if ts and not valid else ('success' if best else 'censored' if complete else 'pending'))
            lower=max((t['duration'] for t in ts if best and t['duration']<best['duration'] and not t[gate+'_success']),default=0.) if best else None
            bracket=(best['duration']-lower)/best['duration'] if best else None
            q=dict(status=status,endpoint=best['duration'] if best else None,maximum_tested_duration=max((t['duration'] for t in ts),default=None),
                   lower_failing_duration=lower,relative_bracket_width=bracket,
                   selected_trace=best['trace'] if best else None,selected_trace_prefix=best['trace_prefix'] if best else None,
                   selected_tail_field_pass_0_1_percent=best['field_tail_pass'] if best else None,
                   selected_field_tail_pass=best['field_tail_pass'] if best else None,
                   nonmonotonic_success_pattern=any(a[gate+'_success'] and not b[gate+'_success'] for a,b in zip(ts,ts[1:])),
                   duration_boundary_unresolved=bool(best and bracket>.01000000001),
                   unresolved_transition_count=sum((not a[gate+'_success']) and b[gate+'_success'] and (b['duration']-a['duration'])/b['duration']>.01000000001 for a,b in zip(ts,ts[1:])))
            results[gate]=q
            rows.append(dict(**{k:c[k] for k in ['case_index','establishment_case_index','model','length','condition','width_mode','width','amplitude','initial_orientation']},definition=gate,**q))
        cases.append(dict(**c,preparation_valid=valid,result=results,tests=ts))
    return rows,cases


def publish(out,stage,complete=False):
    import fcntl
    with (out/'results/.publish.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        existing=out/'results/status.json'
        if existing.exists() and read(existing).get('complete') and not complete:
            return read(existing)
        launch=read(out/'launch.json')
        if launch['stage']=='complete':stage='complete';complete=True
        elif launch['stage'] in ['execution_incomplete','refinement_unresolved']:stage=launch['stage']
        return _publish_locked(out,stage,complete)


def _publish_locked(out,stage,complete=False):
    manifest,groups,pending,errors=collect(out,partial=True)
    rows,cases=rows_and_records(manifest,groups,complete)
    stamp=n.now();results=out/'results'
    n.csvwrite(results/'endpoints.csv',rows)
    with gzip.open(results/'records.json.gz','wt') as f:json.dump(dict(generated_at=stamp,cases=cases,manifest=manifest),f,default=n.convert,allow_nan=False)
    status=dict(generated_at=stamp,complete=bool(complete and not pending and not errors),producer='finished' if complete else 'stopped with error' if stage in ['execution_incomplete','refinement_unresolved'] else 'not started' if stage=='initialized' else 'running',
                stage=stage,cases=len(cases),readouts=len(rows),pending_chunks=len(pending),execution_errors=len(errors),
                completed_trials=sum(len(t) for t in groups.values()),status_counts=dict(Counter(r['status'] for r in rows)),
                figure_state='Final figures only after complete audit; interim data remain explicit.')
    write(results/'status.json',status)
    (results/'index.html').write_text(f'<!doctype html><html><head><meta charset="utf-8"><title>Corrected matched reversal</title></head><body><h1>Five-model matched basal countercue</h1><p>Data generated {stamp}; producer {status["producer"]}; stage {stage}.</p><p>{len(cases)} conditions; {status["completed_trials"]} completed independent trials; {len(pending)} pending chunks; {len(errors)} execution errors.</p><p>Same exact corrected basal endpoints, nominal peak, three current widths. Scalar response and native-field stability are separate. No universal scaling law is assumed.</p><a href="endpoints.csv">Current endpoints</a> · <a href="status.json">Status</a> · <a href="../manifest.json">Frozen experiment</a></body></html>')
    return status


def submit(out,action,round,offset=0,array=None,dependency=None):
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__).resolve()),action,'--out',str(out),'--round',str(round),'--offset',str(offset)])
    if array:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=18:00:00' if action=='monitor' else '--time=12:00:00' if action=='worker' else '--time=01:00:00',
          '--job-name=corrected-reversal-'+action,'--chdir='+str(n.ROOT),
          '--output='+str(out/'work/logs'/f'{action}_r{round}_%A_%a.out'),'--error='+str(out/'work/logs'/f'{action}_r{round}_%A_%a.err')]
    if array:opts.append('--array='+array)
    if dependency:opts.append('--dependency=afterany:'+dependency)
    return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]


def launch_round(out,round):
    verify(out);planned=jobs(out,round)
    if not planned:return advance(out,round)
    text=subprocess.check_output(['scontrol','show','config'],text=True)
    import re
    limit=int(re.search(r'MaxArraySize\s*=\s*(\d+)',text).group(1));arrays=[]
    try:
        for offset in range(0,len(planned),limit):
            count=min(limit,len(planned)-offset)
            arrays.append(submit(out,'worker',round,offset,array=f'0-{count-1}'))
            write(out/'launch.json',dict(generated_at=n.now(),stage='submitting',round=round,arrays=arrays))
    except Exception:
        if arrays:subprocess.run(['scancel',*arrays],check=True)
        write(out/'launch.json',dict(generated_at=n.now(),stage='submission_failed_workers_cancelled',round=round,arrays=arrays));raise
    controller=submit(out,'advance',round,dependency=':'.join(arrays))
    write(out/'launch.json',dict(generated_at=n.now(),stage='no_cue_controls_running' if round==0 else 'duration_grid_running' if round==1 else 'duration_refinement_running',round=round,arrays=arrays,controller=controller))
    publish(out,'no_cue_controls_running' if round==0 else 'duration_grid_running' if round==1 else 'duration_refinement_running')


def refinement(trials):
    existing={t['duration'] for t in trials};extra=set()
    for gate in DEFINITIONS:
        for a,b in zip(trials,trials[1:]):
            if a[gate+'_success'] or not b[gate+'_success']:continue
            lo,hi=a['duration'],b['duration']
            if (hi-lo)/hi<=.01:continue
            candidates=[hi/2**k for k in range(1,11)] if lo==0 else [(lo+hi)/2]
            extra.update(float(d) for d in candidates if d not in existing and d>0)
    return sorted(extra)



def explanation(out,cases):
    case=next(c for c in cases if c['model']=='wave_pinning' and c['length']==15 and c['condition']=='fixed_20pct')
    q=case['result']['full_polarity']
    if q['status']!='success':return
    passing=float(q['endpoint'])
    lower=[t for t in case['tests'] if t['duration']<passing and not t['full_polarity_success']]
    failed=max(lower,key=lambda t:t['duration'])['duration'] if lower else 0.
    state,prep=n.read_checkpoint(n.ROOT/case['checkpoint'],case['model'])
    arrays=dict(passing_duration=passing,failed_duration=failed,establishment_case_index=case['establishment_case_index'])
    for label,d in [('passing',passing),('failed',failed)]:
        cue=e.Patch(case['amplitude'],case['width'],side='right',width_mode=case['width_mode'],fraction=case['width_fraction_at_reference'],background=case.get('background_stimulus',0.))
        pulse=n.run(case['model'],case['length'],state,cue,d,interval=min(STEP,d/20),n=case['cells']) if d>0 else None
        recovery=n.run(case['model'],case['length'],pulse.state if pulse else state,e.withdrawal(case),RECOVERY,interval=STEP,n=case['cells'])
        full=n.concatenate(pulse,recovery,d) if pulse else recovery
        prefix='' if label=='passing' else 'failed_'
        arrays.update({prefix+k:getattr(full,k) for k in ['time','x','active','orientation','contrast']})
    np.savez_compressed(out/'results/explanation_wave_pinning.npz',**arrays)


def advance(out,round):
    manifest,groups,pending,errors=collect(out,partial=True)
    if pending or errors:
        write(out/'launch.json',dict(generated_at=n.now(),stage='execution_incomplete',round=round,arrays=[]))
        publish(out,'execution_incomplete');raise RuntimeError(f'Missing {len(pending)} chunks; {len(errors)} execution errors; no biological censoring assigned.')
    new=[]
    if round==0:
        for c in manifest['cases']:
            if preparation_valid(c,groups[c['case_index']]):
                durations=[d for d in c['durations'] if d>0]
                for start in range(0,len(durations),10):new.append(dict(slot=len(new),case_index=c['case_index'],durations=durations[start:start+10]))
    else:
        for c in manifest['cases']:
            if not preparation_valid(c,groups[c['case_index']]):continue
            durations=refinement(groups[c['case_index']])
            for start in range(0,len(durations),10):new.append(dict(slot=len(new),case_index=c['case_index'],durations=durations[start:start+10]))
    if new and round<=MAX_ROUNDS:
        save_jobs(out,round+1,new);return launch_round(out,round+1)
    rows,cases=rows_and_records(manifest,groups,True)
    unresolved=sum(q['duration_boundary_unresolved'] for c in cases for q in c['result'].values())
    write(out/'results/validation.json',dict(generated_at=n.now(),passed=unresolved==0,
          complete_execution=True,source_hashes_verified=True,checkpoint_hashes_verified=True,
          scalar_gates_from_all_original_tail_samples=True,full_gate_implies_orientation=True,
          shortest_tested_durations_selected=True,unresolved_readouts=unresolved,
          refinement_rounds=round,all_sampled_boundaries_resolved=not new))
    if unresolved or new:
        write(out/'launch.json',dict(generated_at=n.now(),stage='refinement_unresolved',round=round,arrays=[]))
        publish(out,'refinement_unresolved');raise RuntimeError('Bounded refinement exhausted; explicit unresolved readouts, no final figure.')
    audit_path=out/'results/independent_audit.json'
    subprocess.run([sys.executable,str(n.ROOT/'code/audit_reversal.py'),'--campaign',str(out),'--before-publication','--output',str(audit_path)],check=True)
    independent=read(audit_path)
    if not independent.get('passed') or not independent.get('complete'):raise RuntimeError('Independent complete trial audit failed')
    candidate={(r['case_index'],r['definition']):r for r in rows}
    for audited in independent['selected_readouts']:
        expected=candidate[audited['case_index'],audited['definition']]
        for key in ['status','endpoint','selected_tail_field_pass_0_1_percent']:
            if expected[key]!=audited[key]:raise RuntimeError('Independent selected readout differs: '+key)
    validation=read(out/'results/validation.json');validation['independent_audit_passed']=True
    validation['independent_audit_sha256']=n.sha(audit_path);write(out/'results/validation.json',validation)
    explanation(out,cases)
    write(out/'launch.json',dict(generated_at=n.now(),stage='complete',round=round,arrays=[]))
    publish(out,'complete',complete=True)
    subprocess.run([sys.executable,str(n.ROOT/'code/figure_builder.py'),'--campaign',str(out)],check=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['init','launch','worker','advance','status','smoke','monitor','refresh-initialized'])
    parser.add_argument('--out',type=Path,default=DEFAULT)
    parser.add_argument('--round',type=int,default=0);parser.add_argument('--slot',type=int);parser.add_argument('--offset',type=int,default=0)
    args=parser.parse_args();out=args.out.resolve();n.protect_retained_output(out)
    if args.action=='init':initialize(out)
    elif args.action=='refresh-initialized':refresh_initialized(out)
    elif args.action=='launch':
        state=read(out/'launch.json')
        if state['stage']!='initialized':raise ValueError('Campaign already launched; inspect its scheduler state before resubmission')
        launch_round(out,0)
        submit(out,'monitor',0)
    elif args.action=='worker':
        slot=args.offset+args.slot
        try:worker(out,args.round,slot)
        except Exception:
            import traceback
            write(out/'work/chunks'/f'r{args.round}_{slot:05d}.error.json',dict(generated_at=n.now(),traceback=traceback.format_exc()));raise
    elif args.action=='advance':
        try:advance(out,args.round)
        except Exception:
            import traceback
            state=read(out/'launch.json')
            if state['stage'] not in ['complete','refinement_unresolved']:
                write(out/'launch.json',dict(generated_at=n.now(),stage='execution_incomplete',round=args.round,arrays=[],error=traceback.format_exc()))
                publish(out,'execution_incomplete')
            raise
    elif args.action=='monitor':
        import time
        for _ in range(2000):
            state=read(out/'launch.json')
            if state['stage'] in ['complete','execution_incomplete','refinement_unresolved']:return
            publish(out,state['stage'])
            time.sleep(30)
        submit(out,'monitor',args.round)
    elif args.action=='status':print(json.dumps(publish(out,read(out/'launch.json')['stage'],complete=read(out/'launch.json')['stage']=='complete')))
    else:
        manifest=verify(out);case=next(c for c in manifest['cases'] if c['model']=='wave_pinning' and c['length']==15 and c['condition']=='fixed_20pct')
        rows=[]
        for duration in [0.,1.]:
            row,_=trial(case,duration);rows.append(row)
        if not rows[0]['no_cue_preparation_valid']:raise RuntimeError('Representative no-cue control failed')
        print(json.dumps(dict(passed=True,case_index=case['case_index'],trials=rows),indent=2))

if __name__=='__main__':main()
