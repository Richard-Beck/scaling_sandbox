"""Derived fixed-initial-core activity/contrast: no native integrations.

Compare retained fractional half-height organization ROI with the native assay's
whole-cell ROI, retaining raw amounts and three explicitly named floor rules.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import os
import resource
import shutil
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
NATIVE=ROOT/'data/runs/revision_five_model_compensation_baseline_v2_20261019'
DEFAULT=ROOT/'data/runs/revision_compensation_excess_baseline_20261019'
ORGANIZATION=ROOT/'data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv'
MODELS=['wave_pinning','goryachev','holmes_model3','spring','debelly']
FLOORS=['current_minimum','current_sample_quantile_0p05','frozen_initial_minimum']


def now():return datetime.now(timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(path):return json.loads(path.read_text())
def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def csvwrite(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def save(path,**a):
    tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**a);tmp.replace(path)


def initial_core(z,c):
    u=z['initial_active'];q=z['material_fraction'];threshold=u.min()+.5*np.ptp(u)
    assert u[0]>=threshold and np.ptp(u)>0
    below=np.flatnonzero(u<threshold)
    i=int(below[0]) if len(below) else len(u)
    edge=1. if i==len(u) else float(q[i-1]+(threshold-u[i-1])/(u[i]-u[i-1])*(q[i]-q[i-1]))
    if c['model']=='spring':
        dq=np.diff(q);portion=np.clip(edge-q[:-1],0,dq)
        weights=np.zeros(len(q));weights[:-1]+=portion-portion**2/(2*dq);weights[1:]+=portion**2/(2*dq)
        full=np.zeros(len(q));full[:-1]+=dq/2;full[1:]+=dq/2
    else:
        dx=c['length']/len(u);weights=np.clip(edge*c['length']-np.arange(len(u))*dx,0,dx)
        full=np.full(len(u),dx)
    # Spring weights here are material fractions, multiplied by actual L(t).
    initial_length=float(z['initial_physical_length']) if c['model']=='spring' else 1.
    physical=weights*initial_length
    return edge,weights,full,dict(
        basal_background=float(u.min()),initial_profile_peak=float(u.max()),initial_profile_contrast=float(np.ptp(u)),
        organization_core_width=float(physical.sum()),organization_core_raw_amount=float(u@physical),
        organization_core_excess_amount=float(np.maximum(u-u.min(),0)@physical),
        worker_core_raw_amount=float((u*z['initial_core_mask'])@z['initial_weights']),
        worker_core_excess_amount=float((np.maximum(u-u.min(),0)*z['initial_core_mask'])@z['initial_weights']))


def measures(path,c,z):
    key='rac' if c['model']=='debelly' else 'rac_active' if c['model']=='holmes_model3' else 'active'
    with np.load(path) as a:
        t=a['time'].copy();u=a['field_'+key]
        lengths=a['field_length'].copy() if c['model']=='spring' else np.ones(len(t))
    edge,weights,full,initial=initial_core(z,c);core=z['initial_core_mask']
    maskweights=full*core;initial_length=float(z['initial_physical_length']) if c['model']=='spring' else 1.
    basalu=z['initial_active'];basalfloor=initial['basal_background']
    out=dict(time=t,current_global_minimum=u.min(axis=1),current_global_sample_quantile_0p05=np.quantile(u,.05,axis=1),
        global_active_contrast=np.ptp(u,axis=1),fixed_initial_core_peak=np.max(u[:,core],axis=1),
        organization_core_raw_amount=(u@weights)*lengths,worker_core_raw_amount=(u@maskweights)*lengths,
        organization_core_extent=weights.sum()*lengths,worker_core_extent=maskweights.sum()*lengths)
    floors=[out['current_global_minimum'],out['current_global_sample_quantile_0p05'],np.full(len(t),basalfloor)]
    for name,floor in zip(FLOORS,floors):
        excess=np.maximum(u-floor[:,None],0)
        # Normalize every convention by its own basal value at the same floor rule.
        b0=float(np.quantile(basalu,.05)) if name=='current_sample_quantile_0p05' else basalfloor
        for roi,w in [('organization_core',weights),('worker_core',maskweights)]:
            amount=(excess@w)*lengths;denom=float(np.maximum(basalu-b0,0)@w)*initial_length
            assert denom>0
            out[roi+'_'+name+'_excess_amount']=amount
            out[roi+'_'+name+'_excess_fraction']=amount/denom
        peak=np.maximum(out['fixed_initial_core_peak']-floor,0)
        out['worker_core_'+name+'_peak_contrast']=peak
        out['worker_core_'+name+'_peak_contrast_fraction']=peak/(basalu[core].max()-b0)
    out['organization_core_raw_fraction']=out['organization_core_raw_amount']/initial['organization_core_raw_amount']
    out['worker_core_raw_fraction']=out['worker_core_raw_amount']/initial['worker_core_raw_amount']
    return out,initial


def event(t,good):
    first=np.flatnonzero(good)
    if not len(first):return dict(status='no_sampled_entry_through_finite_horizon',lower=float(t[-1]),upper=None)
    i=int(first[0]);bad=np.flatnonzero(~good);last=int(bad[-1]+1) if len(bad) else 0
    return dict(status='already_present' if i==0 else 'sampled_entry',lower=float(t[max(0,i-1)]),upper=float(t[i]),
        retained_status='not_retained_through_last_sample' if last==len(t) else 'present_all_samples' if last==0 else 'sampled_retained_entry',
        retained_lower=None if last==len(t) else float(t[max(0,last-1)]),retained_upper=None if last==len(t) else float(t[last]))


def initialize(out,native):
    assert not out.exists();m=read(native/'manifest.json');prov=read(native/'provenance.json')
    assert sha(native/'manifest.json')==prov['manifest_sha256']
    metrics=list(csv.DictReader(ORGANIZATION.open()));assert len(metrics)==97
    refs=[];sources={str(native/'manifest.json'):sha(native/'manifest.json'),str(native/'provenance.json'):sha(native/'provenance.json'),str(ORGANIZATION):sha(ORGANIZATION)}
    for c in m['cases']:
        reuse=c.get('exact_reuse');resultpath=ROOT/reuse['source_result'] if reuse else native/'workers'/f'{c["slot"]:04d}'/'result.json'
        r=read(resultpath);sources[str(resultpath)]=sha(resultpath)
        if reuse:assert sources[str(resultpath)]==reuse['source_result_sha256']
        else:assert r['case']==c and r['manifest_sha256']==sha(native/'manifest.json')
        record=dict(case=c,result_path=str(resultpath),result_sha256=sha(resultpath),regional_input=str(native/c['regional_input']),
                    regional_input_sha256=c['regional_input_sha256'],phases=[])
        base=ROOT if reuse else native
        on=None if c['case_kind']=='zero_control' else dict(archive=reuse['cue_archive'],archive_sha256=reuse['cue_archive_sha256']) if reuse else r['on']
        if on:record['phases'].append(dict(phase='cue_on',pulse_duration=None,archive=str(base/on['archive']),archive_sha256=on['archive_sha256']))
        for w in reuse['withdrawals'] if reuse else r['withdrawals']:
            record['phases'].append(dict(phase='after_withdrawal',pulse_duration=w['duration'],archive=str(base/w['archive']),archive_sha256=w['archive_sha256'],
                scalar_orientation_success=w['orientation_success'],scalar_full_polarity_success=w['full_polarity_success'],native_field_tail_pass=w['field_tail_pass']))
        refs.append(record)
    (out/'workers').mkdir(parents=True);(out/'logs').mkdir();(out/'frozen_sources').mkdir()
    frozen=out/'frozen_sources'/Path(__file__).name;shutil.copy2(Path(__file__),frozen)
    sources[str(Path(__file__))]=sha(Path(__file__))
    for rel in ['code/mechanisms/initial_polarity_report.py','code/mechanisms/revision_five_model_compensation_campaign_20261019.py']:
        sources[str(ROOT/rel)]=sha(ROOT/rel)
    manifest=dict(created_at=now(),native_run=str(native),native_manifest_sha256=sha(native/'manifest.json'),
        cases=refs,models=MODELS,organization_metrics=str(ORGANIZATION),organization_metrics_sha256=sha(ORGANIZATION),
        conventions=dict(core='Retained organization: fixed basal interpolated half-height boundary with fractional FV/trapezoid integration; worker: original frozen above-halfheight whole-cell set. Spring fixed material core with evolving physical weights; deBelly nominal ROI.',
            floors=FLOORS,excess='Positive part of active concentration minus named background; own basal denominator for each rule. Current floors are descriptive profile references, not separately observed molecular pools.',
            quantile='Unweighted native sample 5th percentile, explicitly a spatial-floor sensitivity.',limits='Regional organization loss/redistribution is distinct from molecular damage or retained replacement.'),
        source_hashes=sources,no_native_integrations=True)
    write(out/'manifest.json',manifest)
    shell=out/'worker.sh';shell.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(frozen)+' worker --out '+str(out)+' --task "$SLURM_ARRAY_TASK_ID"\n');shell.chmod(0o750)
    cmd='sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:30:00 --array=0-4 --job-name=F2-excess-derived --output='+str(out/'logs/%A_%a.out')+' --error='+str(out/'logs/%A_%a.err')+' '+str(shell)
    (out/'submission_command.txt').write_text(cmd+'\n')
    print(json.dumps(dict(run=str(out),manifest_sha256=sha(out/'manifest.json'),native_case_groups=len(refs),derived_workers=5,command=cmd),indent=2))


def worker(out,task):
    started=time.monotonic();m=read(out/'manifest.json');assert m['no_native_integrations']
    for rel,h in m['source_hashes'].items():assert sha(Path(rel))==h,rel
    model=m['models'][task];folder=out/'workers'/model;folder.mkdir(exist_ok=True);assert not (folder/'result.json').exists()
    metrics=list(csv.DictReader(Path(m['organization_metrics']).open()));initials=[];endpoints=[];events=[];sources={};artifacts={};maximum_org_error=0.
    cases=[r for r in m['cases'] if r['case']['model']==model]
    for ref in cases:
        c=ref['case'];slot=c['slot'];rp=Path(ref['regional_input']);assert sha(rp)==ref['regional_input_sha256']
        assert sha(Path(ref['result_path']))==ref['result_sha256']
        with np.load(rp) as zz:z={k:zz[k].copy() for k in zz.files}
        edge,w,full,initial=initial_core(z,c)
        original=next(r for r in metrics if r['model']==model and float(r['domain_length'])==c['length'])
        assert original['checkpoint_sha256']==c['checkpoint_sha256']
        for key,other in [('organization_core_width','old_pole_width'),('organization_core_excess_amount','old_pole_excess_active'),('initial_profile_peak','old_pole_peak_active')]:
            discrepancy=abs(initial[key]-float(original[other]))/max(1.,abs(float(original[other])))
            maximum_org_error=max(maximum_org_error,discrepancy);assert discrepancy<=2e-12,(model,c['length'],key,discrepancy)
        common=dict(case_slot=slot,model=model,length=c['length'],width=c['width'],native_peak_multiplier=c['native_peak_multiplier'],
            requested_native_peak=c['response_amplitude'],coordinate_system=c['assay_coordinate_system'],initial_core_cells=int(z['initial_core_mask'].sum()),
            initial_cue_core_overlap_fraction=c['initial_geometry']['initial_core_footprint_overlap_fraction'])
        initials.append(dict(**common,organization_fixed_material_edge_fraction=edge,**initial))
        for j,phase in enumerate(ref['phases']):
            path=Path(phase['archive']);assert sha(path)==phase['archive_sha256'],str(path);sources[str(path)]=phase['archive_sha256']
            a,_=measures(path,c,z);target=folder/f'{slot:04d}_{j:02d}_{phase["phase"]}.npz';save(target,**a);artifacts[str(target.relative_to(out))]=sha(target)
            # Independent direct reconstruction of the existing whole-cell amount.
            if phase['phase']=='cue_on':windows=[(t,{k:v[a['time']<=t] for k,v in a.items()}) for t in c['pulse_durations']]
            else:windows=[(phase['pulse_duration'],a)]
            for duration,window in windows:
                t=window['time'];assert t[-1]==duration if phase['phase']=='cue_on' else t[-1]==10000.
                base=dict(**common,phase=phase['phase'],pulse_duration=duration,derived_trace=str(target.relative_to(out)),native_archive=str(path),native_archive_sha256=phase['archive_sha256'])
                end=dict(**base,**{k+'_first':float(v[0]) for k,v in window.items() if k!='time'},**{k+'_last':float(v[-1]) for k,v in window.items() if k!='time'})
                for name in ['scalar_orientation_success','scalar_full_polarity_success','native_field_tail_pass']:
                    end[name]=phase.get(name)
                endpoints.append(end)
                for roi in ['organization_core','worker_core']:
                    for floor in FLOORS:
                        fraction=window[roi+'_'+floor+'_excess_fraction']
                        for loss in [.1,.25,.5,.75,.9]:
                            events.append(dict(**base,roi=roi,background_rule=floor,observable='excess_fraction',loss_fraction=loss,
                                observation_end=float(t[-1]),min_sample_interval=float(np.min(np.diff(t))),max_sample_interval=float(np.max(np.diff(t))),
                                **event(t,fraction<=1-loss)))
                    for loss in [.1,.25,.5]:events.append(dict(**base,roi=roi,background_rule='none_raw_activity',observable='raw_fraction',loss_fraction=loss,
                        observation_end=float(t[-1]),**event(t,window[roi+'_raw_fraction']<=1-loss)))
        write(folder/'status.json',dict(generated_at=now(),model=model,completed_cases=len(initials),total_cases=len(cases),state='working'))
    for name,rows in [('initial_core_definitions.csv',initials),('pulse_endpoints.csv',endpoints),('sampled_loss_events.csv',events)]:
        target=folder/name;csvwrite(target,rows);artifacts[str(target.relative_to(out))]=sha(target)
    write(folder/'result.json',dict(generated_at=now(),model=model,completed_case_groups=len(cases),manifest_sha256=sha(out/'manifest.json'),
        native_source_hashes=sources,derived_artifact_hashes=artifacts,maximum_relative_organization_replay_error=maximum_org_error,
        endpoint_rows=len(endpoints),event_rows=len(events),wall_seconds=time.monotonic()-started,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,no_native_integrations=True))
    write(folder/'status.json',dict(generated_at=now(),model=model,completed_cases=len(cases),total_cases=len(cases),state='complete'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','worker']);parser.add_argument('--out',type=Path,default=DEFAULT);parser.add_argument('--native',type=Path,default=NATIVE);parser.add_argument('--task',type=int)
    a=parser.parse_args()
    if a.action=='init':initialize(a.out.resolve(),a.native.resolve())
    else:worker(a.out.resolve(),a.task)
