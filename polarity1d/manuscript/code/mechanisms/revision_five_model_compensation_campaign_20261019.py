"""Expanded all-qualified-length native cue campaign; freeze/preflight/worker only."""
from pathlib import Path
import argparse
import copy
import json
import resource
import shutil
import time
import traceback
import numpy as np
import revision_five_model_compensation_20261019 as p

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
PILOT=ROOT/'data/runs/revision_five_model_compensation_pilot_v2_20261019'
DEFAULT=ROOT/'data/runs/revision_five_model_compensation_baseline_v2_20261019'


def save_npz(path,**arrays):
    tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**arrays);tmp.replace(path)


def event(time,good):
    hits=np.flatnonzero(good)
    if not len(hits):return dict(status='no_sampled_entry_through_finite_horizon',lower=float(time[-1]),upper=None)
    i=int(hits[0]);bad=np.flatnonzero(~good);last=int(bad[-1]+1) if len(bad) else 0
    return dict(status='already_present' if i==0 else 'sampled_entry',lower=float(time[max(0,i-1)]),upper=float(time[i]),
        retained_status='not_retained_through_last_sample' if last==len(time) else 'present_all_samples' if last==0 else 'sampled_retained_entry',
        retained_lower=None if last==len(time) else float(time[max(0,last-1)]),
        retained_upper=None if last==len(time) else float(time[last]))


def regional(tr,c,z,phase):
    u=tr.active;count=len(tr.time);core=z['initial_core_mask'];q=z['material_fraction']
    if c['model']=='spring':
        length=np.asarray(tr.fields['length']);base=np.ones(len(q))/(len(q)-1);base[[0,-1]]*=.5
        weights=length[:,None]*base;fp=(1-q[None,:])*length[:,None]<=c['width']
    else:
        length=np.full(count,c['length']);weights=np.broadcast_to(z['initial_weights'],u.shape)
        fp=np.broadcast_to(z['initial_cue_footprint'],u.shape)
        if c['model']=='debelly':length=c['length']+np.asarray(tr.fields['boundary_separation'])[:,0]
    initial=float(np.dot(z['initial_active']*core,z['initial_weights']))
    old=np.sum(u*weights*core,axis=1);overlap=np.sum(weights*core*fp,axis=1)/np.sum(weights*core,axis=1)
    peak=np.max(np.where(fp,u,-np.inf),axis=1)
    assert np.isfinite(peak).all()
    arrays=dict(time=tr.time,original_core_active_amount=old,original_core_active_fraction=old/initial,
        actual_cue_region_active_peak=peak,actual_cue_region_active_amount=np.sum(u*weights*fp,axis=1),
        fixed_distal_material_region_active_amount=np.sum(u*weights*(q>=.5),axis=1),
        actual_cue_original_core_overlap_fraction=overlap,actual_domain_length=length,
        original_core_physical_or_nominal_extent=np.sum(weights*core,axis=1),
        global_orientation=tr.orientation,active_contrast=tr.contrast)
    inactive='rac_inactive' if c['model']=='holmes_model3' else 'inactive'
    if inactive in tr.fields:
        arrays['actual_cue_region_inactive_mean']=np.sum(tr.fields[inactive]*weights*fp,axis=1)/np.sum(weights*fp,axis=1)
    events=[]
    for loss in [.1,.25,.5]:events.append(dict(phase=phase,observable='fixed_original_core_active_fraction',threshold=1-loss,**event(tr.time,old/initial<=1-loss)))
    if np.any(z['initial_core_mask']&z['initial_cue_footprint']):
        events.append(dict(phase=phase,observable='separated_challenger_peak_above_initial_halfcontrast',threshold=float(z['initial_halfcontrast_threshold']),status='unavailable_initial_cue_core_overlap',lower=None,upper=None))
    else:
        events.append(dict(phase=phase,observable='separated_challenger_peak_above_initial_halfcontrast',threshold=float(z['initial_halfcontrast_threshold']),**event(tr.time,(peak>=z['initial_halfcontrast_threshold'])&(overlap==0))))
    return arrays,events


def scientific_key(c):return (c['model'],c['length'],c['checkpoint_sha256'],c['width'],c['response_amplitude'])


def geometry(c,state,prep,z):
    u=z['initial_active'];above=u>u.min()+.1*np.ptp(u);peak=int(np.argmax(u));left=right=peak
    while left>0 and above[left-1]:left-=1
    while right<len(u)-1 and above[right+1]:right+=1
    core10=np.zeros(len(u),bool);core10[left:right+1]=True
    xp=z['material_fraction']*float(z['initial_physical_length'])
    if c['model']=='debelly':
        xp=prep.x+state['displacement']
        physical_length=c['length']+float(state['right_boundary_displacement'][0]-state['left_boundary_displacement'][0])
    else:physical_length=float(z['initial_physical_length'])
    z.update(initial_core_10pc_mask=core10,initial_physical_positions=xp,
             initial_actual_physical_domain_length=np.array(physical_length))
    return dict(initial_solver_domain_length=c['length'],initial_actual_physical_domain_length=physical_length,
        initial_core_10pc_footprint_overlap_fraction=float(np.dot(z['initial_weights'],core10&z['initial_cue_footprint'])/np.dot(z['initial_weights'],core10)),
        initial_physical_mapping_min_neighbor_spacing=float(np.min(np.diff(xp))),
        initial_requested_support_native_left=max(0.,float(z['initial_physical_length'])-c['width']),
        initial_requested_support_native_right=float(z['initial_physical_length']),
        initial_sampled_support_physical_min=float(xp[z['initial_cue_footprint']].min()),
        initial_sampled_support_physical_max=float(xp[z['initial_cue_footprint']].max()))


def exact_pilot_reuses():
    m=p.verify(PILOT);lookup={};zeros={}
    current=p.n.sources();closure=p.read(PILOT/'provenance.json')['frozen_sources']
    assert all(closure[k]==h for k,h in current.items())
    assert p.n.sha(ROOT/'code/mechanisms/revision_five_model_compensation_20261019.py')==closure['code/mechanisms/revision_five_model_compensation_20261019.py']
    for c in m['cases']:
        folder=PILOT/'workers'/f'{c["slot"]:03d}';result=folder/'result.json';r=p.read(result)
        assert r['case']==c and r['manifest_sha256']==p.n.sha(PILOT/'manifest.json')
        on=folder/'cue.npz';assert p.n.sha(on)==r['cue_archive_sha256']
        reuse=dict(source_manifest=str((PILOT/'manifest.json').relative_to(ROOT)),source_manifest_sha256=p.n.sha(PILOT/'manifest.json'),
            source_result=str(result.relative_to(ROOT)),source_result_sha256=p.n.sha(result),
            cue_archive=str(on.relative_to(ROOT)),cue_archive_sha256=p.n.sha(on),withdrawals=[])
        for w in r['withdrawals']:
            item=copy.deepcopy(w)
            for name,hashname in [('archive','archive_sha256'),('restart','restart_sha256')]:
                path=PILOT/item[name];assert p.n.sha(path)==item[hashname]
                item[name]=str(path.relative_to(ROOT))
            if w['duration']==0:
                zeros[(c['model'],c['length'],c['checkpoint_sha256'])]=dict(source_result=reuse['source_result'],source_result_sha256=reuse['source_result_sha256'],withdrawals=[item])
            else:reuse['withdrawals'].append(item)
        assert [w['duration'] for w in reuse['withdrawals']]==p.PULSES[c['model']]
        lookup[scientific_key(c)]=reuse
    return lookup,zeros


def initialize(out,stage):
    assert not out.exists(),'Fresh campaign required.'
    original=p.read(SOURCE/'manifest.json');assert p.read(SOURCE/'results/status.json')['complete']
    bases=[c for c in original['cases'] if c['eligible'] and c['condition']=='fixed_physical' and c['model'] in p.MODELS]
    assert len(bases)==97
    reuse_on,reuse_zero=exact_pilot_reuses()
    for name in ['inputs/basals','inputs/regions','workers','logs','results','frozen_sources']:(out/name).mkdir(parents=True,exist_ok=True)
    cases=[];jobs=[];inputs={};checks=[]
    for old in sorted(bases,key=lambda c:(p.MODELS.index(c['model']),c['length'])):
        src=ROOT/old['checkpoint'];assert p.n.sha(src)==old['checkpoint_sha256']
        dst=out/'inputs/basals'/f'{old["model"]}_L{old["length"]:g}.npz';shutil.copy2(src,dst)
        inputs[str(dst.relative_to(out))]=p.n.sha(dst)
        state,prep=p.n.read_checkpoint(dst,old['model'])
        widths=[.1,.2,.4] if old['model']=='spring' else [1.5,3.,5.]
        specs=[(width,1.,'constant_on') for width in widths]+[(widths[1],0.,'zero_control')] if stage=='baseline' else [(widths[1],factor,'constant_on') for factor in [.03,.1,.3,3.]]
        for width,factor,kind in specs:
            c={k:copy.deepcopy(v) for k,v in old.items() if k not in ['durations','base_durations','reused_positive_trials']}
            c.update(slot=len(cases),source_case_index=old['case_index'],case_kind=kind,
                width=width,width_mode='physical',native_peak_multiplier=factor,
                baseline_response_amplitude=old['response_amplitude'],response_amplitude=old['response_amplitude']*factor,
                pulse_durations=p.PULSES[old['model']] if kind=='constant_on' else [],
                frozen_checkpoint=str(dst.relative_to(out)),condition='compensation_fixed_native_width',
                assay_coordinate_system='instantaneous_physical' if old['model']=='spring' else 'fixed_nominal_solver' if old['model']=='debelly' else 'fixed_physical')
            if c['model']=='debelly':assert c['center_boundary_motion'] is False
            z=p.regions(c,state,prep);extra_geometry=geometry(c,state,prep,z)
            rp=out/'inputs/regions'/f'{c["slot"]:04d}.npz';save_npz(rp,**z)
            c['regional_input']=str(rp.relative_to(out));c['regional_input_sha256']=p.n.sha(rp);inputs[c['regional_input']]=c['regional_input_sha256']
            check=dict(**p.channel_check(c,state,z),**extra_geometry);checks.append(dict(slot=c['slot'],case_kind=kind,**check))
            c['initial_geometry']=check
            c['exact_reuse']=reuse_zero.get((c['model'],c['length'],c['checkpoint_sha256'])) if kind=='zero_control' else reuse_on.get(scientific_key(c))
            cases.append(c)
            if not c['exact_reuse']:jobs.append(dict(task=len(jobs),case_slot=c['slot']))
    expected=388 if stage=='baseline' else 388
    assert len(cases)==expected
    logical=sum(1 if c['case_kind']=='zero_control' else 5 for c in cases)
    assert logical==(1552 if stage=='baseline' else 1940)
    hashes=p.n.sources()
    for path in [Path(__file__),ROOT/'code/mechanisms/revision_five_model_compensation_20261019.py',ROOT/'code/debelly_corrected_reversal.py',SOURCE/'manifest.json']:
        hashes[str(path.relative_to(ROOT))]=p.n.sha(path)
    for rel in hashes:
        dest=out/'frozen_sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
    manifest=dict(created_at=p.n.now(),experiment='five_model_compensation_'+stage,stage=stage,
        cases=cases,jobs=jobs,qualified_lengths=97,models=p.MODELS,logical_protocol_rows=logical,
        constant_on_groups=sum(c['case_kind']=='constant_on' for c in cases),zero_controls=sum(c['case_kind']=='zero_control' for c in cases),
        exact_reuse_case_groups=sum(bool(c['exact_reuse']) for c in cases),new_native_workers=len(jobs),
        off_duration=p.OFF,primary_full_field_archives=True,native_engine='Accepted platform pilot v2 sampled_run/state_at/save_phase/gates verbatim.',
        primary_readouts='Original-core active amount/fraction, actual cue-region concentration and amount, fixed distal material region amount, actual cue/core overlap, orientation and native coupled fields.',
        scan_limits='No tested success up to finite maximum; no global duration monotonicity or absent success-island claim. Refine all observed transitions.',
        operating_point='Requested peak factors are native model-specific inputs, not matched optical exposure. Width also changes integrated input.',
        numerical='Original model parameters and solver defaults, exact deBelly half-unit update lattice/independent endpoints; Spring raw u/v/boundary coordinates and stretch retained.',
        resource_forecast=dict(full3492_allocated_archive_bytes=73807613952,full3492_raw_field_payload_bytes=145356894456,
            full3492_single_CPU_elapsed_hours=9.0036,reserve_multiplier=[2,3],nonbinding=True,
            basis='Measured baseline shortest/longest platform pilot; excludes future refinements/high-dose numerical companions and dose-dependent stiffness/compression.'),
        scientific_source_manifest_sha256=p.n.sha(SOURCE/'manifest.json'),platform_pilot_manifest_sha256=p.n.sha(PILOT/'manifest.json'))
    p.write(out/'manifest.json',manifest);p.write(out/'provenance.json',dict(created_at=p.n.now(),manifest_sha256=p.n.sha(out/'manifest.json'),frozen_sources=hashes,frozen_inputs=inputs))
    p.csvwrite(out/'results/native_channel_geometry_preflight.csv',checks)
    p.write(out/'results/latest_status.json',dict(generated_at=p.n.now(),producer='frozen_not_submitted',complete=False,total_case_groups=len(cases),total_native_workers=len(jobs),logical_protocol_rows=logical))
    worker=out/'worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(out/'frozen_sources/code/mechanisms'/Path(__file__).name)+' worker --out '+str(out)+' --task "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
    cmd='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 --array=0-'+str(len(jobs)-1)+' --job-name=F2-compensation-'+stage+' --output='+str(out/'logs/%A_%a.out')+' --error='+str(out/'logs/%A_%a.err')+' '+str(worker)
    (out/'submission_command.txt').write_text(cmd+'\n')
    print(json.dumps(dict(out=str(out),manifest_sha256=p.n.sha(out/'manifest.json'),case_groups=len(cases),new_workers=len(jobs),logical_protocols=logical),indent=2))


def worker(out,task):
    started=time.monotonic();m=p.verify(out);c=m['cases'][m['jobs'][task]['case_slot']];assert c['exact_reuse'] is None
    folder=out/'workers'/f'{c["slot"]:04d}';folder.mkdir(exist_ok=True);assert not (folder/'result.json').exists()
    state,prep=p.n.read_checkpoint(out/c['frozen_checkpoint'],c['model'])
    with np.load(out/c['regional_input']) as zz:z={k:zz[k].copy() for k in zz.files}
    on=None;events=[];records=[];on_record=None
    if c['case_kind']=='constant_on':
        p.write(folder/'status.json',dict(generated_at=p.n.now(),state='cue_on',task=task))
        on=p.sampled_run(c,state,p.cue(c),p.sample_times(c['model'],max(c['pulse_durations'])))
        p.save_phase(folder/'cue.npz',on);arrays,ee=regional(on,c,z,'cue_on');events.extend(ee)
        save_npz(folder/'cue_regional.npz',**arrays)
        on_record=dict(archive=str((folder/'cue.npz').relative_to(out)),archive_sha256=p.n.sha(folder/'cue.npz'),regional=str((folder/'cue_regional.npz').relative_to(out)),regional_sha256=p.n.sha(folder/'cue_regional.npz'))
    for j,duration in enumerate(c['pulse_durations'] if on is not None else [0.]):
        restart=copy.deepcopy(state) if duration==0 else p.state_at(on,c['model'],int(np.flatnonzero(on.time==duration)[0]))
        rp=folder/f'restart_{j:02d}.npz';p.save_state(rp,c['model'],restart)
        p.write(folder/'status.json',dict(generated_at=p.n.now(),state='withdrawal',task=task,duration=duration))
        off=p.sampled_run(c,restart,p.e.withdrawal(c),p.sample_times(c['model'],p.OFF,off=True))
        first=p.state_at(off,c['model'],0);aa=list(first.values()) if c['model']=='debelly' else [first];bb=list(restart.values()) if c['model']=='debelly' else [restart]
        absolute=max(float(np.max(abs(a-b))) for a,b in zip(aa,bb));bound=32*np.finfo(float).eps*max(1.,max(float(np.max(abs(b))) for b in bb));assert absolute<=bound
        path=folder/f'off_{j:02d}.npz';p.save_phase(path,off)
        arrays,ee=regional(off,c,z,'off');events.extend([dict(duration=duration,**row) for row in ee])
        reg=folder/f'off_regional_{j:02d}.npz';save_npz(reg,**arrays)
        endpoint={name:float(values[-1]) for name,values in arrays.items() if name!='time'}
        records.append(dict(duration=duration,archive=str(path.relative_to(out)),archive_sha256=p.n.sha(path),
            restart=str(rp.relative_to(out)),restart_sha256=p.n.sha(rp),regional=str(reg.relative_to(out)),regional_sha256=p.n.sha(reg),
            saved_initial_max_absolute_roundoff=absolute,saved_initial_roundoff_bound=bound,regional_endpoint=endpoint,
            **p.gates(off,c,prep,zero=duration==0)))
        del off
    p.write(folder/'result.json',dict(generated_at=p.n.now(),case=c,manifest_sha256=p.n.sha(out/'manifest.json'),
        on=on_record,withdrawals=records,events=events,wall_seconds=time.monotonic()-started,
        peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    p.write(folder/'status.json',dict(generated_at=p.n.now(),state='complete',task=task))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','preflight','worker']);parser.add_argument('--stage',choices=['baseline','strength'],default='baseline');parser.add_argument('--out',type=Path,default=DEFAULT);parser.add_argument('--task',type=int)
    a=parser.parse_args();out=a.out.resolve()
    if a.action=='init':initialize(out,a.stage)
    elif a.action=='preflight':
        m=p.verify(out);print(json.dumps(dict(passed=True,manifest_sha256=p.n.sha(out/'manifest.json'),workers=len(m['jobs']),protocol_rows=m['logical_protocol_rows']),indent=2))
    else:
        try:worker(out,a.task)
        except Exception:
            m=p.read(out/'manifest.json');slot=m['jobs'][a.task]['case_slot'];p.write(out/'workers'/f'{slot:04d}'/'error.json',dict(generated_at=p.n.now(),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
