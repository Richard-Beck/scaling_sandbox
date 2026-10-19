"""Bounded 3x native-input stress screen; separate numerical companions."""
from pathlib import Path
from types import SimpleNamespace
from contextlib import contextmanager
import argparse
import copy
import resource
import shutil
import time
import traceback
import numpy as np
import revision_five_model_compensation_20261019 as p

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'data/runs/revision_five_model_high_strength_screen_20261019'


def initialize(out):
    p.initialize(out)
    manifest=p.read(out/'manifest.json');provenance=p.read(out/'provenance.json');checks=[]
    for c in manifest['cases']:
        c['baseline_response_amplitude']=c['response_amplitude'];c['response_amplitude']*=3;c['native_peak_multiplier']=3.
        state,prep=p.n.read_checkpoint(out/c['frozen_checkpoint'],c['model']);z=p.regions(c,state,prep)
        path=out/c['regional_input'];np.savez_compressed(path,**z)
        c['regional_input_sha256']=p.n.sha(path);provenance['frozen_inputs'][c['regional_input']]=p.n.sha(path)
        checks.append(dict(slot=c['slot'],**p.channel_check(c,state,z)))
    manifest.update(experiment='five_model_high_strength_native_and_numerical_screen',
        primary_peak_multiplier=3.,numerical_screen=dict(
            models_with_tolerance_companion=['wave_pinning','goryachev','holmes_model3','spring'],
            tolerance_factor=.1,max_step_unchanged=True,output_times_unchanged=True,
            prefix_times='all five planned withdrawals',off_reference_times='first and last positive pulse',
            comparison_tolerance_fixed_scales=1e-4,
            debelly='Original half-unit staggered update retained. Step refinement changes per-update boundary relaxation in the accepted model, so no temporal-convergence claim is made; high-dose full-state/channel/finite/restart checks still apply.'))
    rel=str(Path(__file__).relative_to(ROOT));provenance['frozen_sources'][rel]=p.n.sha(Path(__file__))
    shutil.copy2(Path(__file__),out/'frozen_sources'/rel)
    p.write(out/'manifest.json',manifest);provenance['manifest_sha256']=p.n.sha(out/'manifest.json');p.write(out/'provenance.json',provenance)
    p.csvwrite(out/'results/native_channel_preflight.csv',checks)
    worker=out/'worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(out/'frozen_sources/code/mechanisms'/Path(__file__).name)+' worker --out '+str(out)+' --slot "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
    cmd='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:30:00 --array=0-9 --job-name=F2-high-strength-screen --output='+str(out/'logs/%A_%a.out')+' --error='+str(out/'logs/%A_%a.err')+' '+str(worker)
    (out/'submission_command.txt').write_text(cmd+'\n')
    print(manifest['experiment'],p.n.sha(out/'manifest.json'))


def load_phase(path,key):
    with np.load(path) as z:
        fields={name[6:]:z[name].copy() for name in z.files if name.startswith('field_')}
        tr=SimpleNamespace(time=z['time'].copy(),x=z['x'].copy(),orientation=z['orientation'].copy(),
             contrast=z['contrast'].copy(),fields=fields,active=fields[p.active_key(key)])
        if 'raw_boundaries' in z.files:tr.raw_boundaries=z['raw_boundaries'].copy()
    tr.state=p.state_at(tr,key,-1);return tr


@contextmanager
def tighter_solver(calls):
    import scipy.integrate
    import polarity1d.solvers as solvers
    original_scipy=scipy.integrate.solve_ivp;original_rd=solvers.solve_ivp
    def invoke(fun,span,y0,*args,**kwargs):
        before=dict(method=kwargs['method'],rtol=kwargs['rtol'],atol=kwargs['atol'],max_step=kwargs['max_step'])
        kwargs['rtol']*=.1;kwargs['atol']*=.1
        calls.append(dict(native=before,companion=dict(method=kwargs['method'],rtol=kwargs['rtol'],atol=kwargs['atol'],max_step=kwargs['max_step'])))
        return original_scipy(fun,span,y0,*args,**kwargs)
    scipy.integrate.solve_ivp=invoke;solvers.solve_ivp=invoke
    try:yield
    finally:scipy.integrate.solve_ivp=original_scipy;solvers.solve_ivp=original_rd


def worker(out,slot):
    started=time.monotonic();p.worker(out,slot)
    m=p.verify(out);c=m['cases'][slot];key=c['model'];folder=out/'workers'/f'{slot:03d}'
    result=p.read(folder/'result.json');state,prep=p.n.read_checkpoint(out/c['frozen_checkpoint'],key)
    if key=='debelly':
        summary=dict(generated_at=p.n.now(),model=key,length=c['length'],native_update_step=.5,
            status='native_high_dose_restart_checks_pass_no_step_convergence_claim',
            reason='The accepted native per-step boundary relaxation is part of the model definition; halving dt would change it.',
            primary_native_result_sha256=p.n.sha(folder/'result.json'),
            full_state_prefix_checks_pass=result['all_prefix_comparisons_pass'])
    else:
        p.write(folder/'status.json',dict(generated_at=p.n.now(),state='numerical_companion',slot=slot))
        target=folder/'numerical_companion';target.mkdir(exist_ok=True);calls=[];checks=[]
        nominal=load_phase(folder/'cue.npz',key)
        with tighter_solver(calls):
            tight=p.sampled_run(c,state,p.cue(c),p.sample_times(key,max(c['pulse_durations'])))
            p.save_phase(target/'cue.npz',tight)
            for duration in c['pulse_durations']:
                index=int(np.flatnonzero(tight.time==duration)[0]);base_index=int(np.flatnonzero(nominal.time==duration)[0])
                error=p.state_error(p.state_at(tight,key,index),p.state_at(nominal,key,base_index),key,c,remove_spring_origin=key=='spring')
                checks.append(dict(phase='cue_on',duration=duration,maximum_fixed_scaled_difference=max(error.values()),fieldwise_fixed_scaled_difference=error))
            for duration in [c['pulse_durations'][0],c['pulse_durations'][-1]]:
                index=int(np.flatnonzero(tight.time==duration)[0]);restart=p.state_at(tight,key,index)
                off=p.sampled_run(c,restart,p.e.withdrawal(c),p.sample_times(key,p.OFF,off=True))
                j=next(i for i,w in enumerate(result['withdrawals']) if w['duration']==duration)
                base=load_phase(out/result['withdrawals'][j]['archive'],key)
                error=p.state_error(off.state,base.state,key,c,remove_spring_origin=key=='spring')
                diag=p.gates(off,c,prep);original=result['withdrawals'][j]
                agree=all(diag[g]==original[g] for g in ['orientation_success','full_polarity_success','field_tail_pass'])
                p.save_phase(target/f'off_{j:02d}.npz',off)
                checks.append(dict(phase='off',duration=duration,maximum_fixed_scaled_difference=max(error.values()),fieldwise_fixed_scaled_difference=error,
                    all_scalar_and_field_tail_flags_agree=agree,companion_diagnostics=diag,
                    native_flags={g:original[g] for g in ['orientation_success','full_polarity_success','field_tail_pass']}))
                del off,base
        passed=all(q['maximum_fixed_scaled_difference']<=1e-4 and q.get('all_scalar_and_field_tail_flags_agree',True) for q in checks)
        summary=dict(generated_at=p.n.now(),model=key,length=c['length'],status='passes' if passed else 'requires_diagnosis',
            numerical_companion_tolerance_factor=.1,native_max_step_and_output_times_unchanged=True,
            calls=calls,checks=checks,primary_native_result_sha256=p.n.sha(folder/'result.json'),
            numerical_archives={str(q.relative_to(out)):p.n.sha(q) for q in target.glob('*.npz')})
    summary.update(total_wall_seconds=time.monotonic()-started,peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p.write(folder/'numerical_screen.json',summary)
    p.write(folder/'status.json',dict(generated_at=p.n.now(),state='screen_complete',slot=slot,numerical_status=summary['status']))
    assert summary['status']!='requires_diagnosis',summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','preflight','worker']);parser.add_argument('--out',type=Path,default=DEFAULT);parser.add_argument('--slot',type=int)
    a=parser.parse_args();out=a.out.resolve()
    if a.action=='init':initialize(out)
    elif a.action=='preflight':p.verify(out);print('passed',p.n.sha(out/'manifest.json'))
    else:
        try:worker(out,a.slot)
        except Exception:p.write(out/'workers'/f'{a.slot:03d}'/'error.json',dict(generated_at=p.n.now(),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
