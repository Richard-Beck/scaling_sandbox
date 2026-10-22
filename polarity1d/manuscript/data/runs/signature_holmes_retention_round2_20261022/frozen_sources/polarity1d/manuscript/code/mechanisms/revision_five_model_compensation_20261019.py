"""Five-model native compensation pilot and immutable length-series inventory.

No launch action: root submits after independent preflight. Shared constant-ON
prefix states are compared with separately stopped native runs, not declared
bit-identical. The corrected deBelly kernel and its half-unit update are used.
"""
from pathlib import Path
from types import SimpleNamespace
from dataclasses import replace
import argparse
import copy
import csv
import json
import os
import shutil
import resource
import sys
import time
import traceback
import numpy as np

PACKAGE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PACKAGE/'code/core'))
import native as n
import polarity_establishment as e
import settling

SOURCE=PACKAGE/'data/runs/debelly_corrected_three_test_20261008'
DEFAULT=PACKAGE/'data/runs/revision_five_model_compensation_pilot_20261019'
MODELS=['wave_pinning','goryachev','holmes_model3','spring','debelly']
PULSES={'wave_pinning':[1.,10.,100.,1000.,10000.],
        'goryachev':[1.,10.,100.,1000.,5000.],
        'holmes_model3':[1.,10.,100.,1000.,10000.],
        'spring':[.01,.1,1.,10.,50.],
        'debelly':[10.,50.,200.,800.,2000.]}
OFF=10000.
PREFIX_TOLERANCE=1e-4
_model_for=n.model_for


def model_for(key):
    model=_model_for(key)
    return n.ds.DeBelly(replace(model.parameters,center_boundary_motion=False)) if key=='debelly' else model


n.model_for=model_for


def read(p):return json.loads(p.read_text())


def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True);n.write(p,d)


def csvwrite(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp.csv')
    n.csvwrite(tmp,rows);tmp.replace(p)


def cue(c):
    return e.Patch(c['response_amplitude'],c['width'],side='right',width_mode='physical',
                   background=c.get('background_stimulus',0.))


def sample_times(key,end,off=False):
    if key=='debelly':
        # Every retained time lies on the original half-unit native step lattice.
        fine=np.arange(0,min(100.,end)+.25,.5)
        coarse=np.arange(0,end+1.25,2.5)
    elif key=='spring':
        fine=np.r_[np.linspace(0,min(.1,end),101),np.arange(0,min(1.,end)+.005,.01),
                   np.arange(0,min(10.,end)+.05,.1)]
        coarse=np.arange(0,end+1.25,2.5) if off else np.arange(0,end+.25,.5)
    else:
        fine=np.r_[np.arange(0,min(1.,end)+.0005,.001),np.arange(0,min(10.,end)+.005,.01),
                   np.arange(0,min(100.,end)+.05,.1)]
        coarse=np.arange(0,end+1.25,2.5)
    times=np.r_[fine,coarse,0.,end]
    if not off:times=np.r_[times,PULSES[key]]
    times=np.unique(np.round(times[(times>=0)&(times<=end)],12))
    assert times[0]==0 and times[-1]==end and np.all(np.diff(times)>0)
    return times


def active_key(key):return 'rac' if key=='debelly' else 'rac_active' if key=='holmes_model3' else 'active'


def state_at(tr,key,index):
    if key=='debelly':return {name:tr.fields[name][index].copy() for name in n.DB_NAMES}
    if key=='spring':
        y=np.empty(2*tr.active.shape[1]+2)
        y[:-2:2]=tr.fields['active'][index];y[1:-2:2]=tr.fields['inactive'][index]
        y[-2:]=tr.raw_boundaries[index]
        return y
    return np.vstack([tr.fields[name][index] for name in n.model_for(key).state_names])


def sampled_run(c,state,stimulus,times):
    key,L,cells=c['model'],c['length'],c['cells']
    if key=='debelly':
        tr=n.ds.simulate_debelly(n.model_for(key),protocol=stimulus,length=L,cells=cells,
             end_time=float(times[-1]),times=times,step=.5,initial_checkpoint=state)
        fields=tr.fields;active=fields['rac'];ori=n.ds._orientation(active)
        result=SimpleNamespace(time=tr.time,x=tr.x,fields=fields,active=active,
             orientation=ori,contrast=np.ptp(np.maximum(active,0),axis=1))
    elif key=='spring':
        # Capture native LSODA raw state, changing t_eval only. Preserve native
        # output_interval=2.5's max_step=2.5 and native tolerances/RHS closure.
        import scipy.integrate
        original=scipy.integrate.solve_ivp;captures=[]
        def output_capture(fun,span,y0,*args,**kwargs):
            assert kwargs['method']=='LSODA' and kwargs['max_step']==2.5
            assert kwargs['rtol']==2e-7 and kwargs['atol']==2e-9
            keep=times[(times>=span[0])&(times<=span[1])]
            kwargs['t_eval']=np.unique(np.r_[span[0],keep,span[1]])
            sol=original(fun,span,y0,*args,**kwargs);captures.append(sol);return sol
        scipy.integrate.solve_ivp=output_capture
        try:
            tr=n.sp.run_spring(length=L,stimulus=stimulus,end_time=float(times[-1]),
                  output_interval=2.5,cells=cells,initial_state=state)
        finally:scipy.integrate.solve_ivp=original
        assert len(captures)==1
        raw=captures[0].y.T
        fields={name:tr[name] for name in ['active','inactive','length']}
        result=SimpleNamespace(time=tr['time'],x=tr['x_fraction']*L,fields=fields,
             active=tr['active'],orientation=tr['orientation'],contrast=np.ptp(tr['active'],axis=1),
             raw_boundaries=raw[:,-2:].copy())
        np.testing.assert_allclose(result.time,times,rtol=0,atol=1e-12)
    else:
        model=n.model_for(key);state=np.asarray(state)
        assert state.min()>=-1e-8
        tr=n.simulate(model,initial_state=np.maximum(state,0.),length=L,cells=cells,
             protocol=stimulus,times=times,end_time=float(times[-1]),max_step=10.,backend='bdf')
        fields={name:tr.fields[name] for name in model.state_names};active=fields[active_key(key)]
        result=SimpleNamespace(time=tr.time,x=tr.x,fields=fields,active=active,
             orientation=n.ho._orientation(tr) if key=='holmes_model3' else n.rd._orientation(tr),
             contrast=np.ptp(np.maximum(active,0),axis=1))
    result.state=state_at(result,key,-1)
    return result


def save_phase(path,tr):
    data=dict(time=tr.time,x=tr.x,orientation=tr.orientation,contrast=tr.contrast)
    data.update({'field_'+name:value for name,value in tr.fields.items()})
    if hasattr(tr,'raw_boundaries'):data['raw_boundaries']=tr.raw_boundaries
    tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**data);tmp.replace(path)


def save_state(path,key,state):
    data={'state_'+name:value for name,value in state.items()} if key=='debelly' else dict(state=state)
    np.savez_compressed(path,**data)


def state_error(a,b,key,c,remove_spring_origin=False):
    if key=='debelly':
        return {name:float(np.max(abs(a[name]-b[name]))/c['fixed_scales'].get(name,c['length'])) for name in n.DB_NAMES}
    if key=='spring':
        aa=np.array(a,copy=True);bb=np.array(b,copy=True)
        if remove_spring_origin:aa[-2:]-=aa[-2];bb[-2:]-=bb[-2]
        return dict(active=float(np.max(abs(aa[:-2:2]-bb[:-2:2]))/c['fixed_scales']['active']),
                    inactive=float(np.max(abs(aa[1:-2:2]-bb[1:-2:2]))/c['fixed_scales']['inactive']),
                    boundaries=float(np.max(abs(aa[-2:]-bb[-2:]))/c['length']))
    return {name:float(np.max(abs(a[i]-b[i]))/c['fixed_scales'][name]) for i,name in enumerate(n.model_for(key).state_names)}


def regions(c,state,prep):
    key=c['model'];u=state['rac'] if key=='debelly' else state[:-2:2] if key=='spring' else state[0]
    half=float(u.min()+.5*np.ptp(u));above=u>half;peak=int(np.argmax(u));left=right=peak
    while left>0 and above[left-1]:left-=1
    while right<len(u)-1 and above[right+1]:right+=1
    core=np.zeros(len(u),bool);core[left:right+1]=True
    physical_length=float(state[-1]-state[-2]) if key=='spring' else c['length']
    q=np.linspace(0,1,len(u)) if key=='spring' else (np.arange(len(u))+.5)/len(u)
    xx=q*physical_length
    weights=np.full(len(u),physical_length/(len(u)-1)) if key=='spring' else np.full(len(u),c['length']/len(u))
    if key=='spring':weights[[0,-1]]*=.5
    footprint=(physical_length-xx<=c['width'])
    signal=cue(c)(0.,xx,physical_length)
    return dict(x_native=prep.x,material_fraction=q,initial_active=u.copy(),initial_core_mask=core,
        initial_weights=weights,initial_cue_footprint=footprint,initial_signal=signal,
        initial_halfcontrast_threshold=np.array(half),initial_physical_length=np.array(physical_length))


def channel_check(c,state,z):
    key=c['model'];S=z['initial_signal'];x=z['x_native']
    if key=='debelly':
        model=n.model_for(key);p=model.parameters
        assert p.center_boundary_motion is False and p.opto_target=='rac'
        al,be=n.ds._derived_rates(p,state['mca_bound'],state['tension'],S)
        a0,b0=n.ds._derived_rates(p,state['mca_bound'],state['tension'],np.zeros_like(S))
        np.testing.assert_allclose(al-a0,S,rtol=2e-14,atol=2e-14);np.testing.assert_array_equal(be,b0)
        channel='rac_deactivation * S / (1 + rho^2); no explicit inactive Rac'
        extra=p.rac_deactivation*S/(1+state['rho']**2)
    elif key=='spring':
        u,v=state[:-2:2],state[1:-2:2];length=state[-1]-state[-2];delta=3+80*(length-c['length'])
        extra=n.sp.auth['f'](u,v,4+S,5,6,2,delta)-n.sp.auth['f'](u,v,4,5,6,2,delta)
        np.testing.assert_allclose(extra,S*v,rtol=1e-13,atol=1e-13);channel='S * inactive; native stretch/dilution/feedback retained'
    else:
        model=n.model_for(key);change=model.reaction(0.,x,state,S)-model.reaction(0.,x,state,np.zeros_like(S))
        expected=S*state[1]/model.rac_total if key=='holmes_model3' else S*state[1]
        np.testing.assert_allclose(change[0],expected,rtol=1e-13,atol=1e-13)
        np.testing.assert_allclose(change[0]+change[1],0.,atol=1e-14)
        extra=change[0];channel='S * inactive_Rac / Rac_total' if key=='holmes_model3' else 'S * inactive'
    w=z['initial_weights'];core=z['initial_core_mask'];fp=z['initial_cue_footprint']
    return dict(model=key,length=c['length'],checkpoint_sha256=c['checkpoint_sha256'],native_channel=channel,
        sampled_peak=float(S.max()),sampled_spatial_integral=float(np.dot(S,w)),
        initial_extra_native_activity_source_integral=float(np.dot(extra,w)),
        initial_core_footprint_overlap_fraction=float(np.dot(w,core&fp)/np.dot(w,core)),
        initial_core_cells=int(core.sum()),initial_domain_length=float(z['initial_physical_length']),
        support_exceeds_initial_domain=bool(c['width']>float(z['initial_physical_length'])),
        separated_challenger_readout_available=bool(not np.any(core&fp)),
        cue_coordinate_system='instantaneous_physical' if key=='spring' else 'fixed_nominal_solver' if key=='debelly' else 'fixed_physical')


def gates(tr,c,prep,zero=False):
    target=np.arange(0,OFF+1.25,2.5);ii=np.searchsorted(tr.time,target)
    np.testing.assert_allclose(tr.time[ii],target,rtol=0,atol=1e-12)
    tail=target>=OFF-500.;o=tr.orientation[ii][tail];contrast=tr.contrast[ii][tail]
    orientation=bool(o.min()>.05 and contrast.min()>=.05 and np.ptp(o)<.005)
    error=float(np.max(abs(o/abs(c['initial_orientation'])-1)))
    fields=settling.monitored(tr,c['model'],c['length']);initial=settling.monitored(prep,c['model'],c['length'])
    maximum=0.;excursion=0.;detail={}
    for name,y in fields.items():
        assert np.isfinite(y).all();a=y[ii][tail];scale=c['fixed_scales'][name]
        variation=float(np.max(np.ptp(a,axis=0))/scale)
        drift=float(np.max(abs(a[-1]-a[-2]))/2.5*500/scale)
        maximum=max(maximum,variation,drift);detail[name]=dict(range_scaled=variation,projected_drift_500=drift)
        if zero:excursion=max(excursion,float(np.max(abs(y-initial[name][-1]))/scale))
    return dict(orientation_success=orientation,full_polarity_success=bool(orientation and error<=.01),
        final_orientation=float(o[-1]),maximum_relative_polarity_error_last_500=error,
        field_tail_pass=bool(maximum<=.001),field_tail_max_scaled=maximum,field_tail_metrics=detail,
        no_cue_maximum_scaled_excursion=excursion if zero else None,
        no_cue_preparation_valid=bool(excursion<=.01 and not orientation and o.max()<-.05 and np.ptp(o)<.005 and contrast.min()>=.05) if zero else None)


def verify(out):
    p=read(out/'provenance.json');assert n.sha(out/'manifest.json')==p['manifest_sha256']
    for path,digest in p['frozen_sources'].items():assert n.sha(out/'frozen_sources'/path)==digest,path
    for path,digest in p['frozen_inputs'].items():assert n.sha(out/path)==digest,path
    return read(out/'manifest.json')


def initialize(out):
    assert not out.exists(),'Fresh campaign required.'
    original=read(SOURCE/'manifest.json');assert read(SOURCE/'results/status.json')['complete']
    inventory=[];pilot=[]
    for model in MODELS:
        group=sorted([c for c in original['cases'] if c['model']==model and c['condition']=='fixed_physical' and c['eligible']],key=lambda c:c['length'])
        for c in group:
            inventory.append(dict(model=model,length=c['length'],checkpoint=c['checkpoint'],checkpoint_sha256=c['checkpoint_sha256'],
                native_baseline_peak=c['response_amplitude'],widths=[.1,.2,.4] if model=='spring' else [1.5,3.,5.],
                pulse_grid=PULSES[model],boundary_method=c.get('boundary_method','not_applicable')))
        for old in [group[0],group[-1]]:
            c={k:copy.deepcopy(v) for k,v in old.items() if k not in ['durations','base_durations','reused_positive_trials']}
            c.update(slot=len(pilot),source_case_index=old['case_index'],pulse_durations=PULSES[model],
                     width=.2 if model=='spring' else 3.,width_mode='physical',native_peak_multiplier=1.)
            if model=='debelly':assert c['center_boundary_motion'] is False and c['boundary_method']=='independent_endpoints'
            pilot.append(c)
    assert len(inventory)==97 and len(pilot)==10
    for path in ['inputs/basals','inputs/regions','frozen_sources','workers','logs','results']:(out/path).mkdir(parents=True,exist_ok=True)
    inputs={};checks=[]
    for c in pilot:
        src=PACKAGE/c['checkpoint'];assert n.sha(src)==c['checkpoint_sha256']
        dst=out/'inputs/basals'/f'{c["model"]}_L{c["length"]:g}.npz';shutil.copy2(src,dst)
        c['frozen_checkpoint']=str(dst.relative_to(out));inputs[c['frozen_checkpoint']]=n.sha(dst)
        state,prep=n.read_checkpoint(dst,c['model']);z=regions(c,state,prep)
        path=out/'inputs/regions'/f'{c["slot"]:03d}.npz';np.savez_compressed(path,**z)
        c['regional_input']=str(path.relative_to(out));c['regional_input_sha256']=n.sha(path);inputs[c['regional_input']]=n.sha(path)
        checks.append(dict(slot=c['slot'],**channel_check(c,state,z)))
    hashes=n.sources()
    for path in [Path(__file__),PACKAGE/'code/debelly_corrected_reversal.py',SOURCE/'manifest.json',SOURCE/'results/endpoints.csv']:
        hashes[str(path.relative_to(PACKAGE))]=n.sha(path)
    for path in hashes:
        dst=out/'frozen_sources'/path;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(PACKAGE/path,dst)
    manifest=dict(created_at=n.now(),experiment='five_model_compensation_platform_pilot',cases=pilot,
        length_series_inventory=inventory,models=MODELS,qualified_lengths=97,pilot_protocol_rows=60,
        full_stages=dict(baseline_width_duration=1455,middlewidth_strength_duration=970,zero_controls=97,total=2522),
        off_duration=OFF,prefix_comparison_tolerance_on_inherited_fixed_scales=PREFIX_TOLERANCE,
        restart_roundoff_tolerance='32 eps * max(1, initial magnitude), with actual fieldwise errors retained',
        direct_off_reference='First/last positive pulse per worker, plus middle pulses with ON difference >1e-5 fixed scales or near a retained/scalar/tail gate; same final-state tolerance1e-4 and exact gate agreement.',
        design='Exact basal → constant cue ON prefix → five exact-state withdrawal branches + exact-basal no-cue control.',
        numerical='Native BDF/LSODA tolerances and max-step semantics; deBelly original half-unit staggered update, independent local endpoints.',
        scan_limits='No tested success up to maximum is not proof of a minimum beyond it; fail/fail gaps may hide success islands.',
        source_manifest_sha256=n.sha(SOURCE/'manifest.json'))
    write(out/'manifest.json',manifest)
    write(out/'provenance.json',dict(created_at=n.now(),manifest_sha256=n.sha(out/'manifest.json'),frozen_sources=hashes,frozen_inputs=inputs))
    csvwrite(out/'results/native_channel_preflight.csv',checks)
    write(out/'results/latest_status.json',dict(generated_at=n.now(),producer='frozen_not_submitted',complete=False,total_workers=10,completed_workers=0))
    worker=out/'worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(out/'frozen_sources/code/mechanisms'/Path(__file__).name)+' worker --out '+str(out)+' --slot "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
    command='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=02:00:00 --array=0-9 --job-name=F2-five-model-pilot --output='+str(out/'logs/%A_%a.out')+' --error='+str(out/'logs/%A_%a.err')+' '+str(worker)
    (out/'submission_command.txt').write_text(command+'\n')
    print(json.dumps(dict(output=str(out),manifest_sha256=n.sha(out/'manifest.json'),workers=10,protocol_rows=60),indent=2))


def worker(out,slot):
    started=time.monotonic();m=verify(out);c=m['cases'][slot];key=c['model'];folder=out/'workers'/f'{slot:03d}';folder.mkdir(exist_ok=True)
    assert not (folder/'result.json').exists(),'Refusing completed overwrite.'
    state,prep=n.read_checkpoint(out/c['frozen_checkpoint'],key)
    write(folder/'status.json',dict(generated_at=n.now(),state='cue_on',slot=slot))
    phase_started=time.monotonic()
    on=sampled_run(c,state,cue(c),sample_times(key,max(c['pulse_durations'])))
    on_seconds=time.monotonic()-phase_started
    save_phase(folder/'cue.npz',on);records=[];comparisons=[]
    for j,duration in enumerate([0.]+c['pulse_durations']):
        if duration==0:restart=copy.deepcopy(state)
        else:
            ids=np.flatnonzero(on.time==duration);assert len(ids)==1
            restart=state_at(on,key,int(ids[0]))
            direct=n.run(key,c['length'],copy.deepcopy(state),cue(c),duration,
                         interval=min(2.5,duration/20),n=c['cells'])
            reference=direct.state
            error=state_error(restart,reference,key,c,remove_spring_origin=(key=='spring'))
            save_state(folder/f'direct_prefix_{j:02d}.npz',key,reference)
            comparisons.append(dict(duration=duration,maximum_scaled_state_difference=max(error.values()),fieldwise_scaled_difference=error,
                direct_native_sample_interval=min(2.5,duration/20),spring_origin_removed_for_comparison=key=='spring'))
            assert max(error.values())<=PREFIX_TOLERANCE,(key,duration,error)
        save_state(folder/f'restart_{j:02d}.npz',key,restart)
        write(folder/'status.json',dict(generated_at=n.now(),state='withdrawal',slot=slot,duration=duration))
        phase_started=time.monotonic()
        off=sampled_run(c,restart,e.withdrawal(c),sample_times(key,OFF,off=True))
        off_seconds=time.monotonic()-phase_started
        first=state_at(off,key,0);roundoff=state_error(first,restart,key,c)
        # Arraywise scale-free roundoff test, separate from fixed scientific scales.
        aa=list(first.values()) if key=='debelly' else [first];bb=list(restart.values()) if key=='debelly' else [restart]
        absolute=max(float(np.max(abs(a-b))) for a,b in zip(aa,bb))
        bound=32*np.finfo(float).eps*max(1.,max(float(np.max(abs(b))) for b in bb))
        assert absolute<=bound,(absolute,bound)
        path=folder/f'off_{j:02d}.npz';save_phase(path,off)
        diag=gates(off,c,prep,zero=duration==0)
        reference_record=None
        near_gate=(abs(diag['final_orientation']-.05)<.01 or
                   abs(diag['maximum_relative_polarity_error_last_500']-.01)<.005 or
                   abs(diag['field_tail_max_scaled']-.001)<.0005)
        if duration>0 and (j in (1,len(c['pulse_durations'])) or max(error.values())>1e-5 or near_gate):
            phase_started=time.monotonic()
            ref_off=n.run(key,c['length'],copy.deepcopy(reference),e.withdrawal(c),OFF,interval=2.5,n=c['cells'])
            ref_diag=gates(ref_off,c,prep)
            final_error=state_error(off.state,ref_off.state,key,c,remove_spring_origin=(key=='spring'))
            for gate in ['orientation_success','full_polarity_success','field_tail_pass']:assert diag[gate]==ref_diag[gate],(duration,gate)
            assert max(final_error.values())<=PREFIX_TOLERANCE,(key,duration,final_error)
            reference_path=folder/f'direct_off_final_{j:02d}.npz';save_state(reference_path,key,ref_off.state)
            reference_record=dict(reason='first_or_last_pulse' if j in (1,len(c['pulse_durations'])) else 'appreciable_ON_difference_or_near_gate',
                wall_seconds=time.monotonic()-phase_started,
                final_state=str(reference_path.relative_to(out)),final_state_sha256=n.sha(reference_path),
                final_fieldwise_scaled_difference=final_error,maximum_final_scaled_difference=max(final_error.values()),
                diagnostics=ref_diag,all_scalar_and_field_tail_flags_agree=True,
                final_orientation_absolute_difference=abs(diag['final_orientation']-ref_diag['final_orientation']))
            del ref_off
        records.append(dict(duration=duration,archive=str(path.relative_to(out)),archive_sha256=n.sha(path),
            restart=str((folder/f'restart_{j:02d}.npz').relative_to(out)),restart_sha256=n.sha(folder/f'restart_{j:02d}.npz'),
            saved_initial_max_absolute_roundoff=absolute,saved_initial_roundoff_bound=bound,
            saved_initial_fieldwise_scaled_roundoff=roundoff,native_off_wall_seconds=off_seconds,
            off_raw_field_bytes=sum(v.nbytes for v in off.fields.values()),direct_native_OFF_reference=reference_record,**diag))
        del off
    allocated=sum(p.stat().st_blocks*512 for p in folder.rglob('*') if p.is_file())
    raw_field_bytes=sum(v.nbytes for v in on.fields.values())
    result=dict(generated_at=n.now(),case=c,manifest_sha256=n.sha(out/'manifest.json'),cue_archive_sha256=n.sha(folder/'cue.npz'),
        prefixes=comparisons,withdrawals=records,all_prefix_comparisons_pass=True,
        wall_seconds=time.monotonic()-started,native_on_wall_seconds=on_seconds,
        peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        allocated_worker_bytes=allocated,on_raw_field_bytes=raw_field_bytes)
    write(folder/'result.json',result);write(folder/'status.json',dict(generated_at=n.now(),state='complete',slot=slot))


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','preflight','worker']);p.add_argument('--out',type=Path,default=DEFAULT);p.add_argument('--slot',type=int)
    a=p.parse_args();out=a.out.resolve()
    if a.action=='init':initialize(out)
    elif a.action=='preflight':
        m=verify(out);print(json.dumps(dict(passed=True,workers=len(m['cases']),manifest_sha256=n.sha(out/'manifest.json')),indent=2))
    else:
        assert a.slot is not None and 0<=a.slot<10
        try:worker(out,a.slot)
        except Exception:
            write(out/'workers'/f'{a.slot:03d}'/'error.json',dict(generated_at=n.now(),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
