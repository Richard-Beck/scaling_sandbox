"""Continual fixed-slope input: first sampled recovery of initial scalar strength.

Reuse exact manuscript checkpoints, native solvers and the unchanged g*x cue.
There is no withdrawal, post-cue recovery, or persistence gate in this assay.
"""
from pathlib import Path
import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import traceback
import numpy as np

CODE=Path(__file__).resolve().parent
sys.path.insert(0,str(CODE))
import gradient_reversal as g
n=g.n
DEFAULT=n.DATA/'runs/continuous_fixed_slope_20261006'
SOURCE=n.DATA/'runs/fixed_slope_reversal_20261006'
HORIZON=20000.


def read(path):return json.loads(Path(path).read_text())


def sources():
    result=n.sources()
    paths=[Path(__file__),CODE/'continuous_gradient_audit.py',CODE/'continuous_gradient_plot.py',
           CODE/'gradient_reversal.py',SOURCE/'manifest.json',SOURCE/'results/endpoints.csv']
    for path in paths:result[str(path.relative_to(n.ROOT))]=n.sha(path)
    return result


def verify(out):
    p=read(out/'provenance.json')
    if sources()!=p['sources']:raise RuntimeError('Frozen scientific source/input closure changed')
    if n.sha(out/'manifest.json')!=p['manifest_sha256']:raise RuntimeError('Frozen manifest changed')
    return read(out/'manifest.json')


def initialize(out):
    n.protect_retained_output(out)
    if out.exists():raise ValueError('Fresh destination required')
    source=read(SOURCE/'manifest.json')
    cases=[]
    for original in source['cases']:
        if original['condition']!='relative_20pct':continue
        c={k:v for k,v in original.items() if k not in ['durations','base_durations']}
        c.update(source_case_index=c['case_index'],case_index=len(cases),
                 target_orientation=abs(c['initial_orientation']) if c['eligible'] else None,
                 first_segment_end=10*min(d for d in original['durations'] if d>0))
        if c['eligible']:
            p=n.ROOT/c['checkpoint']
            if n.sha(p)!=c['checkpoint_sha256']:raise RuntimeError('Starting checkpoint changed')
        cases.append(c)
    if len(cases)!=110 or sum(c['eligible'] for c in cases)!=97:raise RuntimeError('Unexpected coverage')
    for folder in ['work/cases','work/logs','results','frozen_sources']:(out/folder).mkdir(parents=True,exist_ok=True)
    manifest=dict(created_at=n.now(),experiment='continuous_fixed_physical_slope',models=source['models'],
                  slope_by_model=source['slope_by_model'],preparation='relative_20pct',cases=cases,
                  maximum_application_duration=HORIZON,
                  input_rule='Same added g*x throughout; native basal input remains. No crossfade, withdrawal, reconditioning or slope retuning. Spring follows instantaneous physical position; deBelly keeps native nominal cue coordinates.',
                  outcome='First recorded time with rightward orientation P >= abs(initial left P). No 1% band, contrast/stability condition, confirmation window or basal recovery is imposed; contrast is retained as a diagnostic.',
                  sampling='Geometric segments: first end is 10 times the inherited minimum positive pulse coordinate; multiply end by 10 until common cap 20000. Sample interval=end/2000, with deBelly intervals rounded up to a multiple of its native 0.5 step. Native numerical parameters are unchanged.',
                  reported_time='Earliest passing native sample (upper crossing bracket). Also retain last prior sample and a linearly interpolated estimate; no interpolation changes the success decision.',
                  data_policy='All native fields and orientation/contrast samples through the segment containing first success, or through common horizon for censored cases.',
                  source_campaign_manifest_sha256=n.sha(SOURCE/'manifest.json'))
    n.write(out/'manifest.json',manifest)
    hashes=sources()
    for rel in hashes:
        dest=out/'frozen_sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/rel,dest)
    n.write(out/'provenance.json',dict(generated_at=n.now(),sources=hashes,manifest_sha256=n.sha(out/'manifest.json')))
    collect(out)
    print(json.dumps(dict(cases=len(cases),eligible=sum(c['eligible'] for c in cases),out=str(out))))


def measurement(time,orientation,target):
    hits=np.flatnonzero(orientation>=target)
    if not len(hits):return dict(status='not_reached_within_horizon',duration=None,lower_time=None,upper_time=None,interpolated_time=None,crossing_index=None)
    i=int(hits[0])
    if i==0:raise RuntimeError('Left preparation already passes rightward target')
    lo=float(time[i-1]);hi=float(time[i]);a=float(orientation[i-1]);b=float(orientation[i])
    estimate=lo+(target-a)/(b-a)*(hi-lo)
    return dict(status='success',duration=hi,lower_time=lo,upper_time=hi,interpolated_time=float(estimate),crossing_index=i)


def worker(out,slot):
    manifest=verify(out);c=manifest['cases'][slot]
    if not c['eligible']:return
    dest=out/'work/cases'/f'case_{slot:03d}'
    dest.mkdir(exist_ok=True)
    if (dest/'result.json').exists():raise RuntimeError('Completed case must not be overwritten')
    try:
        p=n.ROOT/c['checkpoint']
        if n.sha(p)!=c['checkpoint_sha256']:raise RuntimeError('Checkpoint hash mismatch')
        state,prep=n.read_checkpoint(p,c['model'])
        start=0.;end=c['first_segment_end'];tr=None;segments=[]
        while True:
            end=min(end,HORIZON)
            interval=end/2000.
            if c['model']=='debelly':interval=.5*max(1,int(np.ceil(interval/.5)))
            segment=n.run(c['model'],c['length'],state,g.response_cue(c),end-start,
                          interval=interval,n=c['cells'],max_step=10.)
            tr=n.concatenate(tr,segment,start)
            segments.append(dict(start=start,end=end,sample_interval=interval))
            if not np.isfinite(tr.orientation).all():raise RuntimeError('Nonfinite orientation')
            m=measurement(tr.time,tr.orientation,c['target_orientation'])
            if m['status']=='success' or end==HORIZON:break
            state=tr.state;start=end;end*=10.
        archive=dest/'trajectory.npz'
        n.save_checkpoint(archive,c['model'],tr.state,tr)
        idx=m['crossing_index']
        row=dict(case_index=slot,source_case_index=c['source_case_index'],model=c['model'],length=c['length'],
                 condition=c['condition'],checkpoint_sha256=c['checkpoint_sha256'],
                 gradient_slope=c['gradient_slope'],initial_orientation=c['initial_orientation'],
                 target_orientation=c['target_orientation'],**m,
                 orientation_at_crossing=float(tr.orientation[idx]) if idx is not None else None,
                 contrast_at_crossing=float(tr.contrast[idx]) if idx is not None else None,
                 final_orientation=float(tr.orientation[-1]),maximum_observed_orientation=float(tr.orientation.max()),
                 simulated_until=float(tr.time[-1]),segments=segments,archive_sha256=n.sha(archive),generated_at=n.now())
        n.write(dest/'result.json',row)
        print(json.dumps({k:row[k] for k in ['case_index','model','length','status','duration','target_orientation','final_orientation']}))
    except Exception:
        n.write(dest/'error.json',dict(generated_at=n.now(),error=traceback.format_exc()))
        raise


def collect(out):
    manifest=read(out/'manifest.json');rows=[];pending=[];errors=[]
    for c in manifest['cases']:
        dest=out/'work/cases'/f'case_{c["case_index"]:03d}'
        if not c['eligible']:
            rows.append(dict(case_index=c['case_index'],model=c['model'],length=c['length'],status='preparation_ineligible',duration=None))
        elif (dest/'result.json').exists():rows.append(read(dest/'result.json'))
        else:
            pending.append(c['case_index'])
            if (dest/'error.json').exists():errors.append(dict(case_index=c['case_index'],**read(dest/'error.json')))
    status=dict(generated_at=n.now(),producer='finished' if not pending else 'running',
                cases=len(manifest['cases']),eligible=sum(c['eligible'] for c in manifest['cases']),
                completed=sum(r['status']!='preparation_ineligible' for r in rows),pending=pending,errors=errors,
                success=sum(r['status']=='success' for r in rows),not_reached=sum(r['status']=='not_reached_within_horizon' for r in rows),
                preparation_ineligible=sum(r['status']=='preparation_ineligible' for r in rows))
    n.write(out/'results/status.json',status)
    n.write(out/'results/records.json',dict(generated_at=n.now(),cases=rows))
    columns=['case_index','model','length','status','duration','lower_time','upper_time','interpolated_time',
             'target_orientation','orientation_at_crossing','contrast_at_crossing','final_orientation',
             'maximum_observed_orientation','simulated_until','gradient_slope']
    if rows:n.csvwrite(out/'results/endpoints.csv',[{k:r.get(k) for k in columns} for r in rows])
    (out/'results/index.html').write_text(f'<!doctype html><html><head><meta charset="utf-8"><title>Continuous gradient status</title></head><body><h1>Continuous fixed-slope experiment</h1><p>Data generated {status["generated_at"]}; producer {status["producer"]}.</p><p>{status["completed"]}/{status["eligible"]} integrations complete; {status["success"]} reach original scalar strength; {status["not_reached"]} do not reach within 20,000. {len(errors)} errors.</p><a href="endpoints.csv">Readouts</a> · <a href="status.json">Status</a></body></html>')
    return status


def submit(out,action,array=None,dependency=None):
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__).resolve()),action,'--out',str(out)])
    if array:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G',
          '--time=02:00:00','--job-name=continual-gradient-'+action,'--chdir='+str(n.ROOT),
          '--output='+str(out/'work/logs'/f'{action}_%A_%a.out'),'--error='+str(out/'work/logs'/f'{action}_%A_%a.err')]
    if array:opts.append('--array='+array)
    if dependency:opts.append('--dependency=afterany:'+dependency)
    return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]


def launch(out):
    m=verify(out)
    if (out/'launch.json').exists():raise RuntimeError('Already submitted')
    slots=[str(c['case_index']) for c in m['cases'] if c['eligible']]
    array=submit(out,'worker',array=','.join(slots))
    finish=submit(out,'finalize',dependency=array)
    n.write(out/'launch.json',dict(generated_at=n.now(),worker_array=array,finalizer=finish))
    print(json.dumps(read(out/'launch.json')))


def finalize(out):
    verify(out);status=collect(out)
    if status['pending'] or status['errors']:raise RuntimeError('Incomplete execution: '+json.dumps(status))
    import continuous_gradient_audit as audit
    audit.audit(out)
    import continuous_gradient_plot as plot
    plot.build(out)
    print(json.dumps(status))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','launch','worker','status','finalize']);p.add_argument('--out',type=Path,default=DEFAULT);p.add_argument('--slot',type=int)
    a=p.parse_args();out=a.out.resolve()
    if a.action=='init':initialize(out)
    elif a.action=='launch':launch(out)
    elif a.action=='worker':worker(out,a.slot)
    elif a.action=='finalize':finalize(out)
    else:print(json.dumps(collect(out),indent=2))
