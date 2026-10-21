"""Bounded numerical companions; exact main 201-cell starts remain canonical.

Stages are separately frozen and never submitted here: ten remapped/relaxed
401-cell basals, twenty width protocols, then ten actual-schedule train protocols.
The six deBelly protocols also receive a 201-cell, rate-preserving quarter-step
companion. A stronger temporal arm separately qualifies quarter-step basals,
reports their shift from the native state, and repeats those six protocols.
Exact-start temporal sensitivities and their failed controls remain separate.
All fate, field qualification and event-bracket differences remain visible.
This module imports no solver until a native worker is invoked.
"""
from __future__ import annotations

import argparse
import copy
import csv
from contextlib import contextmanager, nullcontext
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import sys
import time
import traceback

import numpy as np

sys.dont_write_bytecode = True
ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'polarity1d/manuscript').is_dir())
PACKAGE = ROOT / 'polarity1d/manuscript'
PYTHON = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'
MODELS = ('wave_pinning', 'goryachev', 'debelly', 'holmes_model3', 'spring')
LENGTHS = (19, 39)
DB_NAMES = ('rac', 'rho', 'mca_bound', 'mca_unbound', 'displacement', 'velocity', 'tension',
            'raw_left_displacement', 'raw_right_displacement', 'left_boundary_displacement',
            'right_boundary_displacement')
DB_CONTROLLERS = DB_NAMES[-4:]
RULES = dict(basal_chunk=1000., basal_consecutive_windows=2, basal_tolerance=.0001,
             basal_holdout=10000., basal_holdout_tolerance=.001, basal_cap=60000.,
             zero_excursion_tolerance=.01, tail_window=500., tail_tolerance=.001,
             absolute_D_difference=.05, relative_crossing_midpoint_difference=.1,
             maximum_endpoint_scaled_field_difference=.05,
             observation_cadences=(5., 12.5))


def now(): return datetime.now(timezone.utc).isoformat()
def read(path): return json.loads(Path(path).read_text())
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def write(path, value):
    path=Path(path);path.parent.mkdir(parents=True, exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n');temporary.replace(path)


def modules():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    return tuple(importlib.import_module(name) for name in ('signature_native_20261020',
        'signature_inputs_20261020', 'signature_readouts_20261020', 'signature_campaign_20261020'))


def fv_remap(values, new_cells):
    """Exact overlap transfer of cell averages on the same normalized interval."""
    values=np.asarray(values, dtype=float)
    if values.ndim < 1 or values.shape[-1] < 2 or new_cells < 2 or not np.isfinite(values).all():
        raise ValueError('Finite cell averages and at least two cells are required')
    old_edges=np.linspace(0., 1., values.shape[-1]+1)
    new_edges=np.linspace(0., 1., int(new_cells)+1)
    overlaps=np.maximum(0., np.minimum(new_edges[1:,None], old_edges[None,1:])-
                              np.maximum(new_edges[:-1,None], old_edges[None,:-1]))
    return values @ (overlaps/np.diff(new_edges)[:,None]).T


def remap_state(model, state, new_cells=401):
    """Transfer every state component, preserving positive densities and pools."""
    if model=='debelly':
        if set(state)!=set(DB_NAMES): raise ValueError('All eleven deBelly state fields are required')
        result={}
        for name in DB_NAMES:
            value=np.asarray(state[name], dtype=float)
            if name in DB_CONTROLLERS:
                if not np.all(value==value[0]): raise ValueError('Endpoint-controller arrays must be constant')
                result[name]=np.full(new_cells, value[0])
            else: result[name]=fv_remap(value, new_cells)
        old_pool=float(np.mean(state['mca_bound']+state['mca_unbound']))
        new_pool=float(np.mean(result['mca_bound']+result['mca_unbound']))
        for name in ('rac','rho','mca_bound','mca_unbound'):
            if np.min(state[name])<0 or np.min(result[name])<0: raise ValueError('Negative remapped density')
        pools={'bound_plus_unbound_MCA_mean':dict(before=old_pool, after=new_pool)}
        description='Conservative finite-volume overlap transfer; all four raw/local endpoint arrays copied exactly.'
    elif model=='spring':
        state=np.asarray(state, dtype=float);old_cells=(len(state)-2)//2
        if len(state)!=2*old_cells+2 or old_cells<2: raise ValueError('Invalid Spring state')
        result=np.empty(2*new_cells+2);old_q=np.linspace(0.,1.,old_cells);new_q=np.linspace(0.,1.,new_cells)
        length=float(state[-1]-state[-2])
        if length<=0: raise ValueError('Spring endpoints must remain ordered')
        pools={}
        for name,component in (('active',0),('inactive',1)):
            old=state[component:-2:2]
            if old.min()<0: raise ValueError('Negative Spring density')
            new=np.interp(new_q, old_q, old)
            old_mass=float(np.trapezoid(old,old_q)*length)
            new_mass=float(np.trapezoid(new,new_q)*length)
            if old_mass>0:
                if new_mass<=0: raise ValueError('Positive Spring pool was lost')
                new*=old_mass/new_mass
            result[component:-2:2]=new
            pools[name+'_physical_amount']=dict(before=old_mass,after=float(np.trapezoid(new,new_q)*length))
        result[-2:]=state[-2:]
        description='Nodal linear interpolation; separate positive-density physical amounts corrected; exact raw endpoints retained.'
    else:
        state=np.asarray(state,dtype=float)
        expected=4 if model=='holmes_model3' else 2
        if model not in ('wave_pinning','goryachev','holmes_model3') or state.ndim!=2 or len(state)!=expected:
            raise ValueError('Invalid reaction-diffusion state')
        if state.min()<0: raise ValueError('Negative source density')
        result=fv_remap(state,new_cells)
        groups=(('Rac',0,1),('Rho',2,3)) if model=='holmes_model3' else (('active_plus_inactive',0,1),)
        pools={name:dict(before=float(np.mean(state[a]+state[b])),after=float(np.mean(result[a]+result[b])))
               for name,a,b in groups}
        description='Conservative finite-volume overlap transfer of each species; no inactive offset or kinetic retuning.'
    for item in pools.values():
        if not np.isclose(item['before'],item['after'],rtol=2e-13,atol=2e-14):
            raise ValueError('Remapping changed a declared pool')
    return result,dict(method=description,pools=pools,new_cells=int(new_cells),
                       exact_main_checkpoint_replaced=False,companion_only=True)


def rate_preserving_relaxation(reference=.98, step=.25, reference_step=.5):
    if not 0<reference<=1 or not np.isfinite(step) or step<=0 or reference_step<=0:
        raise ValueError('Positive steps and a relaxation factor in (0,1] are required')
    return reference**(step/reference_step)


@contextmanager
def quarter_step_controller(native):
    original=native.model_for
    baseline=original('debelly')
    configured=native.n.ds.DeBelly(replace(baseline.parameters,
        boundary_relaxation=rate_preserving_relaxation(baseline.parameters.boundary_relaxation)))
    native.model_for=lambda key:configured if key=='debelly' else original(key)
    try: yield dict(native_step=.25,controller_relaxation=configured.parameters.boundary_relaxation,
        reference_relaxation=baseline.parameters.boundary_relaxation,reference_step=.5,
        conversion='rho(dt)=rho(0.5)**(dt/0.5); physical controller decay rate preserved',
        native_baseline_changed=False)
    finally: native.model_for=original


def observation_indices(time, events, cadence):
    """Coarsen each phase clock while keeping every boundary and the tail end."""
    time=np.asarray(time,dtype=float)
    if time.ndim!=1 or np.any(np.diff(time)<=0) or cadence<=0: raise ValueError('Invalid observation clock')
    wanted=[time[0],time[-1]]
    for event in events:
        wanted.extend((event['start'],event['stop']))
        wanted.extend(event['start']+np.arange(0.,event['stop']-event['start'],cadence))
    positions=np.searchsorted(time,np.asarray(wanted));positions=np.clip(positions,0,len(time)-1)
    previous=np.maximum(positions-1,0)
    choose=np.abs(time[previous]-wanted)<np.abs(time[positions]-wanted)
    positions=np.where(choose,previous,positions)
    if not np.allclose(time[positions],wanted,rtol=0,atol=1e-8):
        raise ValueError('Requested coarsening is absent from the native saved clock')
    return np.unique(positions)


def event_comparison(reference, companion, tolerance=.1):
    if reference['observed'] != companion['observed']:
        return dict(agrees=False,classification='crossing_censoring_changed')
    if not reference['observed']:
        return dict(agrees=True,classification='both_no_sampled_crossing',relative_midpoint_difference=None)
    if reference.get('lower_time') is None or companion.get('lower_time') is None:
        agrees=reference.get('left_censored')==companion.get('left_censored')
        return dict(agrees=agrees,classification='left_censored_comparison',relative_midpoint_difference=None)
    a,b=reference['lower_time'],reference['upper_time'];c,d=companion['lower_time'],companion['upper_time']
    overlap=max(a,c)<=min(b,d)
    midpoint=.5*(a+b);relative=abs(.5*(c+d)-midpoint)/max(abs(midpoint),1e-12)
    return dict(agrees=bool(overlap or relative<=tolerance),brackets_overlap=bool(overlap),
                relative_midpoint_difference=float(relative),classification='observed_brackets_compared')


def comparison_classification(reference, companion, event):
    if reference['readout_eligible']!=companion['readout_eligible']: return 'readout_eligibility_changed'
    if (reference['final_d'] is not None and reference['final_d']>=.5)!=(companion['final_d'] is not None and companion['final_d']>=.5): return 'replacement_fate_changed'
    if reference['field_tail_pass']!=companion['field_tail_pass']: return 'finite_field_tail_status_changed'
    if not event['agrees']: return 'crossing_timing_or_censoring_changed'
    if reference['final_d'] is not None and abs(reference['final_d']-companion['final_d'])>RULES['absolute_D_difference']: return 'final_D_difference_exceeds_declared_tolerance'
    if companion.get('maximum_endpoint_scaled_field_difference',0.)>RULES['maximum_endpoint_scaled_field_difference']: return 'endpoint_native_field_difference_exceeds_declared_tolerance'
    return 'agreement_within_declared_comparison_tolerances'


def phase_endpoint_comparison(reference, companion, events):
    """Keep all phase differences, including those hidden by equal final fate."""
    rows=[]
    for event in events:
        endpoints=[]
        for measurement in (reference,companion):
            series=measurement['series'];time=np.asarray(series['time'])
            selected=np.flatnonzero(abs(time-event['stop'])<=1e-8)
            if len(selected)!=1:raise ValueError('Phase endpoint is absent or ambiguous in saved native readouts')
            i=int(selected[0])
            endpoints.append({key:series[key][i] for key in ('D_half','D_half_raw','D_half_eligible',
                'share_half','contrast_fraction_of_basal','amount_fraction_of_basal','ineligibility_reason')})
        coarse,fine=endpoints;difference=(fine['D_half']-coarse['D_half']
            if coarse['D_half'] is not None and fine['D_half'] is not None else None)
        eligibility_agrees=coarse['D_half_eligible']==fine['D_half_eligible']
        sign_agrees=(coarse['D_half']>=0)==(fine['D_half']>=0) if difference is not None else None
        rows.append(dict(name=event['name'],time=event['stop'],main201=coarse,companion=fine,
            D_difference=difference,eligibility_agrees=eligibility_agrees,challenger_bias_sign_agrees=sign_agrees,
            within_D_tolerance=difference is None or abs(difference)<=RULES['absolute_D_difference']))
    values=[abs(row['D_difference']) for row in rows if row['D_difference'] is not None]
    return dict(endpoints=rows,maximum_absolute_eligible_D_difference=max(values,default=None),
        all_endpoint_eligibility_agrees=all(row['eligibility_agrees'] for row in rows),
        all_eligible_endpoint_D_within_tolerance=all(row['within_D_tolerance'] for row in rows),
        all_available_challenger_bias_signs_agree=all(row['challenger_bias_sign_agrees'] is not False for row in rows),
        absolute_D_tolerance=RULES['absolute_D_difference'])


def _source_basals(paths):
    result={}
    for path in paths:
        controls_path=Path(path).parent/'control_flags.json'
        controls=read(controls_path)['controls'] if controls_path.exists() else []
        for row in read(path)['basals']:
            if row['requested_reference_length'] not in LENGTHS: continue
            row=dict(row)
            control=row.get('no_cue_control')
            if control is None:
                matches=[item for item in controls if item['basal_id']==row['basal_id']]
                if len(matches)!=1:raise ValueError('Missing or ambiguous exact-source basal control')
                control=matches[0]
                row.update(no_cue_control=control,source_control_path=str(controls_path.resolve()),
                           source_control_sha256=sha(controls_path))
            if not row.get('qualified') or not control.get('passed'):
                raise ValueError('Source basal lacks a passed native full-field/D control')
            checkpoint=Path(row['checkpoint'])
            if sha(checkpoint)!=row['checkpoint_sha256']: raise ValueError('Source basal hash changed')
            key=row['basal_id']
            if key in result and result[key]['checkpoint_sha256']!=row['checkpoint_sha256']:
                raise ValueError('Conflicting exact 201-cell basals')
            result[key]=row
    expected={f'{model}:L{length}' for model in MODELS for length in LENGTHS}
    if set(result)!=expected: raise ValueError('All ten exact model/length starts are required')
    return result


def _main_cases(paths, stage, models=MODELS, unavailable=()):
    _,_,_,campaign=modules();result={}
    for path in paths:
        run=Path(path).resolve();manifest=campaign.verify(run)
        for case in manifest['cases']:
            if case['model'] not in models:continue
            if case['requested_reference_length'] not in LENGTHS: continue
            if stage=='width':
                if case['stage']!='width' or case['integrated_area_factor']!=1. or case['requested_width_label'] not in (.95,7.6):continue
                assay='width_'+str(case['requested_width_label'])
            else:
                if case['stage'] not in ('train','train_resource_pilot') or (case['pulse_count'],case['gap_multiple_of_pulse'],case['total_on_factor_of_calibration'])!=(4,4.,1.5):continue
                assay='train_4_gap4_on1p5'
            worker=run/'workers'/f"{case['task']:04d}";record=read(worker/'result.json')
            campaign.verify_worker_artifacts(run,case,record)
            key=(case['model'],case['requested_reference_length'],assay)
            if key in result: raise ValueError('Ambiguous main source protocol')
            result[key]=dict(case=case,trial=record['trial'],worker=str(worker),
                checkpoint=record['basal']['checkpoint'],checkpoint_sha256=record['basal']['checkpoint_sha256'],
                main_manifest_sha256=sha(run/'manifest.json'),main_worker_result_sha256=sha(worker/'result.json'))
    assays=('width_0.95','width_7.6') if stage=='width' else ('train_4_gap4_on1p5',)
    expected={(model,length,assay) for model in models for length in LENGTHS for assay in assays}
    absent={(row['model'],row['requested_reference_length'],'train_4_gap4_on1p5') for row in unavailable}
    if not absent.issubset(expected) or set(result)!=(expected-absent):
        raise ValueError('Selected main protocols do not match the complete declared reference inventory')
    return result


def _unavailable_train_sources(path):
    if path is None:return []
    collection=read(path)
    if collection.get('native_missing_or_failed') is not False:
        raise ValueError('Missing or failed execution is not a calibration-unavailable reference')
    selected=[row for row in collection['conditions'] if row['model'] in MODELS
        and row['requested_reference_length'] in LENGTHS
        and (row['pulse_count'],row['gap_multiple_of_pulse'],row['total_on_factor_of_calibration'])==(4,4.,1.5)]
    keys=[(row['model'],row['requested_reference_length']) for row in selected]
    if len(keys)!=len(set(keys)):raise ValueError('Duplicate calibration-unavailable references')
    for row in selected:
        source=Path(row['calibration_result_path']);result=read(source)
        if sha(source)!=row['calibration_result_sha256'] or not row['calibration_artifact_hashes_verified']:
            raise ValueError('Unavailable reference lacks verified calibration identity')
        if result['basal']['checkpoint_sha256']!=row['exact_basal_sha256'] or result['manifest_sha256']!=row['calibration_manifest_sha256']:
            raise ValueError('Unavailable reference changed its exact basal or source campaign')
        if row['train_executed'] or row['calibration_duration'] is not None or result['calibration']['train_ready']:
            raise ValueError('A defined train cannot be labelled calibration-unavailable')
        if result['calibration']['status']!='no_sampled_commitment_through_finite_horizon':
            raise ValueError('Unsupported calibration-unavailable reason')
        receipt=source.parents[2]/'results/artifact_verification.json'
        matches=[item for item in read(receipt)['workers'] if item['worker_result_sha256']==sha(source)]
        if len(matches)!=1 or not matches[0]['passed'] or not matches[0]['all_scientific_artifacts_hashed']:
            raise ValueError('Unavailable calibration artifact audit did not pass')
    return selected


def _freeze(out, stage, cases, copies, unavailable=()):
    out=Path(out).resolve()
    if out.exists() or not out.is_relative_to(PACKAGE/'data/runs'):
        raise ValueError('Choose a fresh numerical campaign directory beneath data/runs')
    for folder in ('inputs','frozen_sources','workers','logs','results'):(out/folder).mkdir(parents=True)
    inputs={}
    for relative,source in copies.items():
        source=Path(source);before=sha(source);destination=out/relative
        destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
        if sha(source)!=before or sha(destination)!=before:raise ValueError('Input changed while freezing')
        inputs[relative]=before
    sources={}
    for source in sorted((PACKAGE/'code').rglob('*')):
        if not source.is_file() or '__pycache__' in source.parts or source.suffix not in ('.py','.json','.txt','.md'):continue
        relative=str(source.relative_to(ROOT));destination=out/'frozen_sources'/relative
        destination.parent.mkdir(parents=True,exist_ok=True);before=sha(source);shutil.copy2(source,destination)
        if sha(source)!=before or sha(destination)!=before:raise ValueError('Source changed while freezing')
        sources[relative]=before
    for source in (PACKAGE/'data/inputs/fixed_scales.json',
                   PACKAGE/'data/runs/debelly_corrected_three_test_20261008/manifest.json'):
        relative=str(source.relative_to(ROOT));destination=out/'frozen_sources'/relative
        destination.parent.mkdir(parents=True,exist_ok=True);before=sha(source);shutil.copy2(source,destination)
        if sha(source)!=before or sha(destination)!=before:raise ValueError('Runtime scientific input changed while freezing')
        sources[relative]=before
    manifest=dict(created_utc=now(),stage=stage,cases=cases,workers=len(cases),source_sha256=sources,
        input_sha256=inputs,rules=RULES,native_main_states_replaced=False,
        scope='Numerical companions only; all differences and failed qualifications remain explicit.',
        unavailable_reference_cases=list(unavailable),requested_reference_slots=len(cases)+len(unavailable),
        temporal_conversion='deBelly quarter-step uses sqrt(0.98) to preserve controller physical decay time.',
        additional_DB_quarter_step_protocols=sum(c['model']=='debelly' for c in cases) if stage in ('width','train') else 0)
    write(out/'manifest.json',manifest);write(out/'freeze.json',dict(manifest_sha256=sha(out/'manifest.json')))
    script=out/'frozen_sources'/Path(__file__).relative_to(ROOT)
    worker=out/'worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+shlex.quote(PYTHON)+' -B '+shlex.quote(str(script))+' worker --out '+shlex.quote(str(out))+' --task "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
    command=['sbatch','--parsable','--partition=red','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=4G','--time=01:00:00',f'--array=0-{len(cases)-1}','--job-name=F2-signature-numerics-'+stage,'--output='+str(out/'logs/%A_%a.out'),'--error='+str(out/'logs/%A_%a.err'),str(worker)]
    (out/'submission_command.txt').write_text(shlex.join(command)+'\n');verify(out)
    return dict(output=str(out),stage=stage,workers=len(cases),submitted=False,manifest_sha256=sha(out/'manifest.json'))


def prepare_basals(out, collections, temporal=False):
    basals=_source_basals(collections);cases=[];copies={}
    for model in MODELS:
        if temporal and model!='debelly':continue
        for length in LENGTHS:
            source=basals[f'{model}:L{length}'];relative=f'inputs/{model}_L{length}_201.npz'
            copies[relative]=source['checkpoint']
            cases.append(dict(task=len(cases),model=model,requested_reference_length=length,
                native_length=source['native_length'],fixed_scales=source['fixed_scales'],
                source_checkpoint=relative,source_checkpoint_sha256=source['checkpoint_sha256'],cells=201 if temporal else 401))
    for i,path in enumerate(collections):
        copies[f'inputs/source_basals_{i}.json']=path
        controls=Path(path).parent/'control_flags.json'
        if controls.exists():copies[f'inputs/source_controls_{i}.json']=controls
    return _freeze(out,'temporal_basals' if temporal else 'basals',cases,copies)


def prepare_protocols(out, stage, numerical_basals, main_runs, temporal=False, unavailable_trains=None):
    basals={row['basal_id']:row for row in read(numerical_basals)['basals']}
    unavailable=_unavailable_train_sources(unavailable_trains)
    if unavailable and (stage!='train' or temporal):raise ValueError('Unavailable-reference records apply to primary mesh train comparisons only')
    selected=_main_cases(main_runs,stage,('debelly',) if temporal else MODELS,unavailable)
    copies={'inputs/numerical_basals.json':numerical_basals};cases=[]
    if unavailable_trains is not None:copies['inputs/unavailable_trains.json']=unavailable_trains
    for i,row in enumerate(unavailable):
        if basals[f"{row['model']}:L{row['requested_reference_length']}"]['source_checkpoint_sha256']!=row['exact_basal_sha256']:
            raise ValueError('Unavailable reference does not share its numerical companion source basal')
        copies[f'inputs/unavailable_{i}_calibration_result.json']=row['calibration_result_path']
        copies[f'inputs/unavailable_{i}_artifact_verification.json']=Path(row['calibration_result_path']).parents[2]/'results/artifact_verification.json'
    for (model,length,assay),source in sorted(selected.items()):
        basal=basals[f'{model}:L{length}']
        if not basal['qualified'] or basal['source_checkpoint_sha256']!=source['checkpoint_sha256']:
            raise ValueError('Missing qualified same-branch companion basal')
        key=f'{model}_L{length}_{assay}';base=f'inputs/{key}'
        copies[base+'/main201.npz']=Path(source['worker'])/'trial/trace.npz'
        copies[base+'/main_readouts.json']=Path(source['worker'])/'trial/readouts.json'
        copies[base+'/main_basal201.npz']=source['checkpoint']
        copies[base+('/temporal_basal201.npz' if temporal else '/basal401.npz')]=basal['checkpoint']
        cases.append(dict(task=len(cases),model=model,requested_reference_length=length,
            native_length=source['case']['native_length'],assay=assay,
            main_case=source['case'],main_trial=source['trial'],input_base=base,
            source_checkpoint_sha256=source['checkpoint_sha256'],source_main_manifest_sha256=source['main_manifest_sha256'],
            source_main_result_sha256=source['main_worker_result_sha256'],cells=201 if temporal else 401,
            companion_basal_shift=basal.get('source_basal_field_shift'),
            companion_basal_no_cue_control=basal.get('no_cue_control')))
    return _freeze(out,('temporal_' if temporal else '')+stage,cases,copies,unavailable)


def verify(out):
    out=Path(out);manifest=read(out/'manifest.json')
    if sha(out/'manifest.json')!=read(out/'freeze.json')['manifest_sha256']:raise ValueError('Frozen manifest changed')
    for relative,digest in manifest['source_sha256'].items():
        if sha(out/'frozen_sources'/relative)!=digest:raise ValueError('Frozen source changed')
    for relative,digest in manifest['input_sha256'].items():
        if sha(out/relative)!=digest:raise ValueError('Frozen input changed')
    if [case['task'] for case in manifest['cases']]!=list(range(manifest['workers'])):raise ValueError('Case inventory changed')
    return manifest


def _relax(native, model, length, state, scales, folder, cells=401, debelly_step=.5):
    zero=native.n.ConstantInput(0.);elapsed=0.;consecutive=0;checks=[];last=None
    while elapsed+RULES['basal_chunk']<=RULES['basal_cap']:
        trace=native.run_native(model,length,state,zero,
            native.sample_times(model,RULES['basal_chunk'],early=False,debelly_step=debelly_step),cells=cells,debelly_step=debelly_step)
        native.save_phase(folder/f'chunk_{len(checks):03d}.npz',trace)
        gate=native.native_field_tail_gate(trace,model,length,scales,window=RULES['basal_chunk'],tolerance=RULES['basal_tolerance'])
        elapsed+=RULES['basal_chunk'];checks.append(dict(elapsed=elapsed,kind='convergence',gate=gate))
        state=trace.state;last=trace
        consecutive=consecutive+1 if gate['field_tail_pass'] else 0
        if consecutive<RULES['basal_consecutive_windows']:continue
        if elapsed+RULES['basal_holdout']>RULES['basal_cap']:break
        hold=native.run_native(model,length,state,zero,
            native.sample_times(model,RULES['basal_holdout'],early=False,debelly_step=debelly_step),cells=cells,debelly_step=debelly_step)
        native.save_phase(folder/f'holdout_{len(checks):03d}.npz',hold)
        initial=native.settling.monitored(trace,model,length);fields=native.settling.monitored(hold,model,length)
        excursion=max(float(np.max(abs(value-initial[name][-1]))/scales[name]) for name,value in fields.items())
        final=native.native_field_tail_gate(hold,model,length,scales,window=RULES['basal_chunk'],tolerance=RULES['basal_tolerance'])
        elapsed+=RULES['basal_holdout'];checks.append(dict(elapsed=elapsed,kind='holdout',maximum_scaled_excursion=excursion,gate=final))
        last,state=hold,hold.state
        if excursion<=RULES['basal_holdout_tolerance'] and final['field_tail_pass']:
            return last,dict(qualified=True,elapsed=elapsed,checks=checks)
        consecutive=0
    return last,dict(qualified=False,elapsed=elapsed,checks=checks,reason='native_relaxation_did_not_qualify_within_cap')


def _basal_shift(native, source, companion, case):
    old=native.settling.monitored(source,case['model'],case['native_length'])
    new=native.settling.monitored(companion,case['model'],case['native_length']);values={}
    for name,reference in old.items():
        reference=np.asarray(reference[-1]);actual=np.asarray(new[name][-1])
        if reference.shape!=actual.shape:
            actual=(np.interp(np.linspace(0.,1.,len(reference)),np.linspace(0.,1.,len(actual)),actual)
                    if case['model']=='spring' else fv_remap(actual,len(reference)))
        values[name]=float(np.max(abs(actual-reference))/case['fixed_scales'][name])
    return dict(maximum_scaled_difference=max(values.values()),per_field_scaled_difference=values,
        comparison='Native monitored fields; no fitted spatial shift or amplitude normalization.')


def _debelly_pool_gate(trace, reference_mean):
    pools=np.mean(trace.fields['mca_bound']+trace.fields['mca_unbound'],axis=1)
    error=float(np.max(abs(pools-reference_mean))/max(abs(reference_mean),1e-12))
    minimum=min(float(np.min(trace.fields[name])) for name in ('rac','rho','mca_bound','mca_unbound'))
    return dict(passed=bool(error<=1e-10 and minimum>=-1e-10),maximum_relative_MCA_pool_error=error,
        minimum_chemical_density=minimum,reference_MCA_mean=float(reference_mean),
        pool='Native nominal-coordinate bound-plus-unbound MCA; mechanical boundary displacement is separate.',
        mass_tolerance=1e-10,positivity_tolerance=1e-10)


def _readouts(native,campaign,readouts,trace,basal,case):
    return campaign._analyze(trace,basal,case,readouts)


def _crossing(readouts,measurement,case,events):
    times=np.asarray(measurement['series']['time']);values=np.asarray(measurement['series']['D_half'],dtype=float)
    if case['assay'].startswith('width_'):
        keep=times<=events[0]['stop']+1e-10;times,values=times[keep],values[keep]
    return readouts.first_crossing(times,values,.5)


def _measure_variant(native,campaign,readouts,trace,basal,case,events,output):
    main=case['main_case'];measurement=_readouts(native,campaign,readouts,trace,basal,main)
    tail=native.native_field_tail_gate(trace,main['model'],main['native_length'],main['fixed_scales'])
    event=_crossing(readouts,measurement,case,events);sampling=[]
    for cadence in RULES['observation_cadences']:
        index=observation_indices(trace.time,events,cadence)
        sample=dict(series={key:[value[i] for i in index] for key,value in measurement['series'].items()})
        sampled_event=_crossing(readouts,sample,case,events)
        sampled_trace=copy.copy(trace)
        for name in ('time','active','orientation','contrast','left','right','physical_length'):
            setattr(sampled_trace,name,np.asarray(getattr(trace,name))[index])
        sampled_trace.fields={name:value[index] for name,value in trace.fields.items()}
        if hasattr(trace,'raw_boundaries'):sampled_trace.raw_boundaries=trace.raw_boundaries[index]
        sampled_tail=native.native_field_tail_gate(sampled_trace,main['model'],main['native_length'],main['fixed_scales'])
        sampling.append(dict(cadence=cadence,frames=len(index),all_phase_boundaries_retained=True,
            event=sampled_event,event_comparison=event_comparison(event,sampled_event),final_D=measurement['final_d'],
            sampled_field_tail=sampled_tail,field_tail_flag_agrees=sampled_tail['field_tail_pass']==tail['field_tail_pass']))
    write(output/'readouts.json',measurement)
    return dict(final_d=measurement['final_d'],readout_eligible=measurement['readout_eligible'],
        ineligibility_reason=measurement['ineligibility_reason'],field_tail_pass=tail['field_tail_pass'],
        field_tail=tail,first_sampled_D_event=event,observation_sensitivity=sampling)


def _endpoint_fields(native,trace,case,source_prefix):
    """Native-coordinate endpoint differences, without fitted alignment/shift."""
    model,length=case['model'],case['native_length']
    with np.load(Path(source_prefix)/'main201.npz',allow_pickle=False) as archive:
        coarse={name[6:]:archive[name][-1] for name in archive.files if name.startswith('field_')}
    if model=='debelly':
        coarse={name:coarse[name] for name in ('rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','boundary_separation')}
        coarse['displacement']=coarse['displacement']-coarse['displacement'].mean()
        coarse['physical_length']=length+coarse.pop('boundary_separation')
    fine={name:value[-1] for name,value in native.settling.monitored(trace,model,length).items()}
    scales=case['main_case']['fixed_scales'];detail={}
    for name,reference in coarse.items():
        value=np.asarray(fine[name]);reference=np.asarray(reference)
        if value.ndim==1 and value.shape!=reference.shape:
            value=(np.interp(np.linspace(0.,1.,len(reference)),np.linspace(0.,1.,len(value)),value)
                   if model=='spring' else fv_remap(value,len(reference)))
        if value.shape!=reference.shape:raise ValueError('Native endpoint field dimensions differ')
        detail[name]=float(np.max(abs(value-reference))/scales[name])
    return dict(maximum_scaled_difference=max(detail.values()),per_field_scaled_difference=detail,
        comparison='Same reference fraction; FV overlap projection or Spring nodal interpolation; no fitted shift or amplitude normalization.')


def _phases(inputs,case,step=.5):
    main=case['main_case'];model=main['model'];kind='node' if model=='spring' else 'cell'
    cue=inputs.AreaCue(main['cue_area'],main['cue_width'],grid_kind=kind);zero=inputs.AreaCue(0.,main['cue_width'],grid_kind=kind)
    phases=[]
    for event in case['main_trial']['events']:
        duration=event['stop']-event['start']
        if duration<=0:continue
        if model=='debelly' and not np.isclose(duration/step,round(duration/step),rtol=0,atol=1e-8):raise ValueError('Source phase is off the requested native lattice')
        name=event['name'];on=name=='cue_on' or name.endswith('_on')
        phases.append(dict(name=name,duration=duration,cue=cue if on else zero))
    return phases


def worker(out,task):
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    out=Path(out);manifest=verify(out);case=manifest['cases'][task];folder=out/'workers'/f'{task:04d}'
    if (folder/'result.json').exists():raise FileExistsError('Completed numerical workers are immutable')
    folder.mkdir(parents=True,exist_ok=True);write(folder/'status.json',dict(state='running',generated_utc=now()))
    started=time.monotonic();native,inputs,readouts,campaign=modules();model=case['model'];length=case['native_length']
    if manifest['stage'] in ('basals','temporal_basals'):
        original,source=native.read_checkpoint(out/case['source_checkpoint'],model)
        temporal=manifest['stage']=='temporal_basals';cells=case['cells'];step=.25 if temporal else .5
        remapped,remapping=remap_state(model,original,cells)
        with quarter_step_controller(native) if temporal else nullcontext(None) as conversion:
            relaxed,qualification=_relax(native,model,length,remapped,case['fixed_scales'],folder/'relaxation',cells,step)
        result=dict(case=case,remapping=remapping,qualification=qualification,qualified=False)
        if temporal:result['controller_conversion']=conversion
        if qualification['qualified']:
            checkpoint=folder/('qualified_basal_quarter_step201.npz' if temporal else 'qualified_basal401.npz')
            native.save_checkpoint(checkpoint,model,relaxed.state,relaxed)
            state,basal=native.read_checkpoint(checkpoint,model)
            with quarter_step_controller(native) if temporal else nullcontext(None):
                zero=native.run_native(model,length,state,native.n.ConstantInput(0.),
                    native.sample_times(model,4000. if model=='debelly' else 10000.,debelly_step=step),cells=cells,debelly_step=step)
            native.save_phase(folder/'zero_control.npz',zero)
            control=native.no_cue_control_gate(zero,basal,case['fixed_scales'])
            gate=native.native_field_tail_gate(zero,model,length,case['fixed_scales'])
            readout_case=dict(model=model,native_length=length,cue_width=.1 if model=='spring' else 3.8)
            zero_readout=campaign._analyze(zero,basal,readout_case,readouts)
            write(folder/'zero_readouts.json',zero_readout)
            D_values=[value for value in zero_readout['series']['D_half'] if value is not None]
            D_flipped=any(value>=.5 for value in D_values)
            polarized=bool(relaxed.orientation[-1]<-.05 and relaxed.contrast[-1]>=.05)
            qualified=bool(control['no_cue_control_pass'] and gate['field_tail_pass'] and polarized
                           and not D_flipped and zero_readout['readout_eligible'])
            if temporal:
                source_pool=float(np.mean(original['mca_bound']+original['mca_unbound']))
                pool=_debelly_pool_gate(zero,source_pool)
                for phase in (folder/'relaxation').glob('*.npz'):
                    _,history=native.read_checkpoint(phase,model)
                    checked=_debelly_pool_gate(history,source_pool)
                    pool['maximum_relative_MCA_pool_error']=max(pool['maximum_relative_MCA_pool_error'],checked['maximum_relative_MCA_pool_error'])
                    pool['minimum_chemical_density']=min(pool['minimum_chemical_density'],checked['minimum_chemical_density'])
                    pool['passed']=pool['passed'] and checked['passed']
                qualified=qualified and pool['passed'];result['MCA_conservation_and_positivity']=pool
            result.update(qualified=qualified,checkpoint=str(checkpoint.resolve()),checkpoint_sha256=sha(checkpoint),
                no_cue_control=control,field_tail=gate,qualified_left_polarity=polarized,
                source_basal_field_shift=_basal_shift(native,source,relaxed,case),
                spontaneous_D_replacement=D_flipped,
                basal_id=f"{model}:L{case['requested_reference_length']}",source_checkpoint_sha256=case['source_checkpoint_sha256'],cells=cells)
    elif manifest['stage'].startswith('temporal_'):
        prefix=out/case['input_base'];state,basal=native.read_checkpoint(prefix/'temporal_basal201.npz',model)
        with quarter_step_controller(native) as conversion:
            trace,events=native.run_phases(model,length,state,_phases(inputs,case,.25),folder/'qualified_quarter_step201',cells=201,debelly_step=.25)
        measured=_measure_variant(native,campaign,readouts,trace,basal,case,events,folder/'qualified_quarter_step201')
        endpoint=_endpoint_fields(native,trace,case,prefix)
        measured.update(endpoint_native_fields=endpoint,maximum_endpoint_scaled_field_difference=endpoint['maximum_scaled_difference'])
        reference=read(prefix/'main_readouts.json');source_event=_crossing(readouts,reference,case,case['main_trial']['events'])
        coarse=dict(final_d=reference['final_d'],readout_eligible=reference['readout_eligible'],field_tail_pass=case['main_trial']['field_tail_pass'])
        event=event_comparison(source_event,measured['first_sampled_D_event'])
        result=dict(case=case,qualified_quarter_step201=measured,controller_conversion=conversion,
            main201=coarse,main201_event=source_event,event_comparison=event,
            classification=comparison_classification(coarse,measured,event),
            source_basal_field_shift=case['companion_basal_shift'],
            no_cue_control=case['companion_basal_no_cue_control'],
            MCA_conservation_and_positivity=_debelly_pool_gate(trace,float(np.mean(state['mca_bound']+state['mca_unbound']))))
    else:
        prefix=out/case['input_base'];state,basal=native.read_checkpoint(prefix/'basal401.npz',model)
        trace,events=native.run_phases(model,length,state,_phases(inputs,case),folder/'mesh401',cells=401)
        measured=_measure_variant(native,campaign,readouts,trace,basal,case,events,folder/'mesh401')
        endpoint=_endpoint_fields(native,trace,case,prefix)
        measured.update(endpoint_native_fields=endpoint,maximum_endpoint_scaled_field_difference=endpoint['maximum_scaled_difference'])
        reference=read(prefix/'main_readouts.json');source_event=_crossing(readouts,reference,case,case['main_trial']['events'])
        coarse=dict(final_d=reference['final_d'],readout_eligible=reference['readout_eligible'],field_tail_pass=case['main_trial']['field_tail_pass'])
        event=event_comparison(source_event,measured['first_sampled_D_event'])
        result=dict(case=case,mesh401=measured,main201=coarse,main201_event=source_event,
                    event_comparison=event,classification=comparison_classification(coarse,measured,event))
        if model=='debelly':
            state201,basal201=native.read_checkpoint(prefix/'main_basal201.npz',model)
            with quarter_step_controller(native) as conversion:
                zero=native.run_native(model,length,state201,native.n.ConstantInput(0.),native.sample_times(model,4000.,debelly_step=.25),cells=201,debelly_step=.25)
                native.save_phase(folder/'quarter_step_zero.npz',zero)
                zero_control=native.no_cue_control_gate(zero,basal201,case['main_case']['fixed_scales'])
                temporal,temporal_events=native.run_phases(model,length,state201,_phases(inputs,case,.25),folder/'quarter_step201',cells=201,debelly_step=.25)
            temporal_measurement=_measure_variant(native,campaign,readouts,temporal,basal201,case,temporal_events,folder/'quarter_step201')
            endpoint=_endpoint_fields(native,temporal,case,prefix)
            temporal_measurement.update(endpoint_native_fields=endpoint,maximum_endpoint_scaled_field_difference=endpoint['maximum_scaled_difference'])
            temporal_event=event_comparison(source_event,temporal_measurement['first_sampled_D_event'])
            result['quarter_step201']=dict(conversion=conversion,no_cue_control=zero_control,
                measurement=temporal_measurement,event_comparison=temporal_event,
                classification=comparison_classification(coarse,temporal_measurement,temporal_event))
    result.update(manifest_sha256=sha(out/'manifest.json'),generated_utc=now(),wall_seconds=time.monotonic()-started,
        peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result['artifact_sha256']={str(path.relative_to(out)):sha(path) for path in folder.rglob('*')
        if path.is_file() and path.name not in ('status.json','result.json','failed.json')}
    write(folder/'result.json',result);write(folder/'status.json',dict(state='complete',generated_utc=now()))
    return dict(task=task,complete=True,qualified=result.get('qualified'),classification=result.get('classification'))


def collect(out):
    out=Path(out);manifest=verify(out);records=[];missing=[];failures=[]
    for case in manifest['cases']:
        path=out/'workers'/f"{case['task']:04d}"/'result.json'
        if not path.exists():
            missing.append(case['task'])
            failed=path.with_name('failed.json')
            if failed.exists():failures.append(dict(task=case['task'],failure=read(failed)))
            continue
        result=read(path)
        if result['manifest_sha256']!=sha(out/'manifest.json') or result['case']!=case:raise ValueError('Worker provenance changed')
        for relative,digest in result['artifact_sha256'].items():
            if sha(out/relative)!=digest:raise ValueError('Numerical artifact changed')
        records.append(result)
    status=dict(generated_utc=now(),stage=manifest['stage'],complete=not missing,workers=len(manifest['cases']),completed=len(records),missing_tasks=missing,
        recorded_failed_tasks=failures,
        unavailable_reference_cases=manifest.get('unavailable_reference_cases',[]),
        native_main_states_replaced=False)
    if not missing and manifest['stage'] in ('basals','temporal_basals'):
        rows=[dict(basal_id=r.get('basal_id',f"{r['case']['model']}:L{r['case']['requested_reference_length']}"),
            qualified=r['qualified'],checkpoint=r.get('checkpoint'),checkpoint_sha256=r.get('checkpoint_sha256'),
            source_checkpoint_sha256=r['case']['source_checkpoint_sha256'],model=r['case']['model'],
            native_length=r['case']['native_length'],fixed_scales=r['case']['fixed_scales'],cells=r['case']['cells'],
            source_basal_field_shift=r.get('source_basal_field_shift'),no_cue_control=r.get('no_cue_control'),
            MCA_conservation_and_positivity=r.get('MCA_conservation_and_positivity')) for r in records]
        name='basals_quarter_step201.json' if manifest['stage']=='temporal_basals' else 'basals401.json'
        write(out/'results'/name,dict(basals=rows,all_qualified=all(r['qualified'] for r in rows),manifest_sha256=sha(out/'manifest.json')))
        status['all_companion_basals_qualified']=all(r['qualified'] for r in rows)
    if not missing:
        write(out/'results/comparisons.json',dict(records=records,manifest_sha256=sha(out/'manifest.json')))
        if manifest['stage'] not in ('basals','temporal_basals'):
            rows=[]
            for record in records:
                common=dict(task=record['case']['task'],model=record['case']['model'],
                    requested_reference_length=record['case']['requested_reference_length'],assay=record['case']['assay'],
                    native_length=record['case']['native_length'],main201=record['main201'],main201_event=record['main201_event'])
                if 'mesh401' in record:
                    rows.append(dict(common,variant='401_cells',comparison=record['mesh401'],
                        event_comparison=record['event_comparison'],classification=record['classification'],
                        separately_qualified_basal=True,companion_control_pass=True,
                        interpretation_status='qualified_spatial_companion'))
                if 'qualified_quarter_step201' in record:
                    rows.append(dict(common,variant='qualified_201_cells_quarter_step_rate_preserving_controller',
                        comparison=record['qualified_quarter_step201'],event_comparison=record['event_comparison'],
                        classification=record['classification'],controller_conversion=record['controller_conversion'],
                        source_basal_field_shift=record['source_basal_field_shift'],
                        no_cue_control=record['no_cue_control'],companion_control_pass=record['no_cue_control']['no_cue_control_pass'],
                        MCA_conservation_and_positivity=record['MCA_conservation_and_positivity'],
                        interpretation_status='separately_qualified_temporal_companion'))
                if 'quarter_step201' in record:
                    temporal=record['quarter_step201']
                    rows.append(dict(common,variant='201_cells_quarter_step_rate_preserving_controller',
                        comparison=temporal['measurement'],event_comparison=temporal['event_comparison'],
                        classification=temporal['classification'],controller_conversion=temporal['conversion'],
                        no_cue_control=temporal['no_cue_control'],
                        companion_control_pass=temporal['no_cue_control']['no_cue_control_pass'],
                        interpretation_status=('exact_main_start_temporal_companion_control_pass'
                            if temporal['no_cue_control']['no_cue_control_pass'] else
                            'unqualified_temporal_diagnostic_no_cue_control_failed')))
            counts={}
            for row in rows:counts[row['classification']]=counts.get(row['classification'],0)+1
            summary=dict(generated_utc=now(),manifest_sha256=sha(out/'manifest.json'),
                actual_companion_protocols=len(rows),classification_counts=counts,rows=rows,
                unavailable_reference_cases=manifest.get('unavailable_reference_cases',[]),
                qualified_companion_protocols=sum(row['companion_control_pass'] for row in rows),
                failed_companion_control_protocols=sum(not row['companion_control_pass'] for row in rows),
                all_differences_retained=True,native_main_states_replaced=False)
            write(out/'results/comparison_summary.json',summary)
    write(out/'results/status.json',status);return status


def audit_basals(out):
    """Verify remapping, native restart continuity and gates from saved fields."""
    from signature_analysis_20261020 import _same_numeric_tree
    out=Path(out);manifest=verify(out)
    if manifest['stage'] not in ('basals','temporal_basals'):raise ValueError('Expected a numerical-basal campaign')
    native,inputs,readouts,campaign=modules();records=[]
    for case in manifest['cases']:
        folder=out/'workers'/f"{case['task']:04d}";result=read(folder/'result.json')
        for relative,digest in result['artifact_sha256'].items():
            if sha(out/relative)!=digest:raise ValueError('Numerical basal artifact changed')
        source,source_trace=native.read_checkpoint(out/case['source_checkpoint'],case['model'])
        previous,receipt=remap_state(case['model'],source,case['cells'])
        _same_numeric_tree(receipt,result['remapping'],'conservative_remap')
        previous_trace=None;gates=[];pool_checks=[]
        for i,check in enumerate(result['qualification']['checks']):
            name='chunk' if check['kind']=='convergence' else 'holdout'
            _,trace=native.read_checkpoint(folder/f'relaxation/{name}_{i:03d}.npz',case['model'])
            start=native.restart_error(previous,trace)
            if not start['passed']:raise ValueError('401-cell relaxation lost an exact native restart')
            gate=native.native_field_tail_gate(trace,case['model'],case['native_length'],case['fixed_scales'],
                window=RULES['basal_chunk'],tolerance=RULES['basal_tolerance'])
            _same_numeric_tree(gate,check['gate'],'relaxation_gate')
            if check['kind']=='holdout':
                target=native.settling.monitored(previous_trace,case['model'],case['native_length'])
                actual=native.settling.monitored(trace,case['model'],case['native_length'])
                excursion=max(float(np.max(abs(value-target[key][-1]))/case['fixed_scales'][key])
                              for key,value in actual.items())
                _same_numeric_tree(excursion,check['maximum_scaled_excursion'],'holdout_excursion')
            previous=native.state_at(trace);previous_trace=trace;gates.append(dict(kind=name,passed=True))
            if manifest['stage']=='temporal_basals':
                pool_checks.append(_debelly_pool_gate(trace,float(np.mean(source['mca_bound']+source['mca_unbound']))))
        state,basal=native.read_checkpoint(result['checkpoint'],case['model'])
        if not native.restart_error(previous,basal)['passed']:raise ValueError('Qualified checkpoint differs from final relaxation state')
        _,zero=native.read_checkpoint(folder/'zero_control.npz',case['model'])
        if not native.restart_error(state,zero)['passed']:raise ValueError('401-cell no-cue control changed its native start')
        control=native.no_cue_control_gate(zero,basal,case['fixed_scales'])
        _same_numeric_tree(control,result['no_cue_control'],'no_cue_control')
        field=native.native_field_tail_gate(zero,case['model'],case['native_length'],case['fixed_scales'])
        _same_numeric_tree(field,result['field_tail'],'no_cue_tail')
        readout_case=dict(model=case['model'],native_length=case['native_length'],cue_width=.1 if case['model']=='spring' else 3.8)
        measured=campaign._analyze(zero,basal,readout_case,readouts)
        _same_numeric_tree(measured,read(folder/'zero_readouts.json'),'no_cue_readouts')
        if 'source_basal_field_shift' in result:
            _same_numeric_tree(_basal_shift(native,source_trace,basal,case),result['source_basal_field_shift'],'basal_field_shift')
        if manifest['stage']=='temporal_basals':
            pool=_debelly_pool_gate(zero,float(np.mean(source['mca_bound']+source['mca_unbound'])))
            for checked in pool_checks:
                pool['maximum_relative_MCA_pool_error']=max(pool['maximum_relative_MCA_pool_error'],checked['maximum_relative_MCA_pool_error'])
                pool['minimum_chemical_density']=min(pool['minimum_chemical_density'],checked['minimum_chemical_density'])
                pool['passed']=pool['passed'] and checked['passed']
            _same_numeric_tree(pool,result['MCA_conservation_and_positivity'],'MCA_conservation_and_positivity')
            if not pool['passed']:raise ValueError('Quarter-step companion lost conserved MCA or positivity')
        if not result['qualified'] or not control['no_cue_control_pass'] or not field['field_tail_pass']:
            raise ValueError('401-cell companion was not independently qualified')
        records.append(dict(task=case['task'],model=case['model'],requested_reference_length=case['requested_reference_length'],
            passed=True,relaxation_checks=gates,native_restart_continuity=True,
            maximum_no_cue_scaled_excursion=control['no_cue_maximum_scaled_excursion'],
            native_main_checkpoint_replaced=False,source_basal_field_shift=result.get('source_basal_field_shift')))
    summary=dict(generated_utc=now(),passed=True,manifest_sha256=sha(out/'manifest.json'),
        audit_source_sha256=sha(Path(__file__)),qualified_basals=len(records),records=records)
    write(out/'results/independent_basal_audit.json',summary);return summary


def audit_protocols(out):
    """Recompute companion measurements from every archived native field trace."""
    from signature_analysis_20261020 import _same_numeric_tree, audit_native_phase_join, declared_state_fields
    out=Path(out).resolve();manifest=verify(out)
    if manifest['stage'] not in ('width','train','temporal_width','temporal_train'):raise ValueError('Expected numerical companion protocols')
    native,inputs,readouts,campaign=modules();records=[]
    for case in manifest['cases']:
        folder=out/'workers'/f"{case['task']:04d}";result=read(folder/'result.json');prefix=out/case['input_base']
        if result['case']!=case or result['manifest_sha256']!=sha(out/'manifest.json'):
            raise ValueError('Numerical worker identity differs from its frozen case')
        for relative,digest in result['artifact_sha256'].items():
            if sha(out/relative)!=digest:raise ValueError('Numerical companion artifact changed')
        variants=([('qualified_quarter_step201','qualified_quarter_step201','temporal_basal201.npz',201,.25)]
                  if manifest['stage'].startswith('temporal_') else [('mesh401','mesh401','basal401.npz',401,.5)])
        if case['model']=='debelly' and not manifest['stage'].startswith('temporal_'):
            variants.append(('quarter_step201','quarter_step201','main_basal201.npz',201,.25))
        for variant,directory,basal_name,cells,step in variants:
            _,basal=native.read_checkpoint(prefix/basal_name,case['model'])
            state=native.state_at(basal)
            _,trace=native.read_checkpoint(folder/directory/'trace.npz',case['model'])
            if trace.cells!=cells or not native.restart_error(state,trace)['passed']:
                raise ValueError('Numerical protocol lost its qualified exact native start')
            measured=campaign._analyze(trace,basal,case['main_case'],readouts)
            _same_numeric_tree(measured,read(folder/directory/'readouts.json'),'companion_readouts')
            stored=result[variant]['measurement'] if variant=='quarter_step201' else result[variant]
            tail=native.native_field_tail_gate(trace,case['model'],case['native_length'],case['main_case']['fixed_scales'])
            _same_numeric_tree(tail,stored['field_tail'],'companion_native_tail')
            phases=[];elapsed=0.
            for i,phase in enumerate(_phases(inputs,case,step)):
                path=folder/directory/f'phase_{i:03d}.npz'
                phases.append(dict(index=i,name=phase['name'],start=elapsed,stop=elapsed+phase['duration'],
                                   archive=dict(path=str(path),sha256=sha(path))))
                elapsed+=phase['duration']
            state_fields,declaration=declared_state_fields(case['model'],out,manifest)
            joined=audit_native_phase_join(folder/directory/'trace.npz',phases,dict(time=trace.time,active=trace.active),state_fields=state_fields)
            event=_crossing(readouts,measured,case,phases)
            _same_numeric_tree(event,stored['first_sampled_D_event'],'companion_sampled_crossing')
            endpoint=_endpoint_fields(native,trace,case,prefix)
            _same_numeric_tree(endpoint,stored['endpoint_native_fields'],'companion_endpoint_fields')
            records.append(dict(task=case['task'],variant=variant,model=case['model'],
                requested_reference_length=case['requested_reference_length'],assay=case['assay'],
                passed=True,full_native_phase_join=joined,exact_native_basal_restart=True,
                native_readout_and_tail_reconstruction=True,native_state_declaration=declaration))
    summary=dict(generated_utc=now(),passed=True,manifest_sha256=sha(out/'manifest.json'),
        analysis_source_sha256=sha(Path(__file__)),audited_actual_protocols=len(records),records=records)
    write(out/'results/independent_protocol_audit.json',summary);return summary


def summarize(out, runs):
    """Expose longitudinal sensitivity separately from final replacement fate."""
    out=Path(out).resolve()
    if out.exists():raise FileExistsError('Choose a fresh numerical-summary folder')
    out.mkdir(parents=True);records=[];phase_rows=[];profile_rows=[];sources=[];unavailable=[]
    source_hash=sha(Path(__file__))
    for run in runs:
        run=Path(run).resolve();manifest=verify(run)
        status=read(run/'results/status.json');audit=read(run/'results/independent_protocol_audit.json')
        if not status['complete'] or not audit['passed']:raise ValueError('Numerical producer and independent native audit must complete first')
        sources.append(dict(run=str(run),manifest_sha256=sha(run/'manifest.json'),
            native_audit_sha256=sha(run/'results/independent_protocol_audit.json')))
        unavailable.extend(manifest.get('unavailable_reference_cases',[]))
        for case in manifest['cases']:
            folder=run/'workers'/f"{case['task']:04d}";result=read(folder/'result.json');prefix=run/case['input_base']
            reference=read(prefix/'main_readouts.json');events=case['main_trial']['events']
            variants=([('qualified_quarter_step201','qualified_quarter_step201',result['qualified_quarter_step201'],True)]
                      if 'qualified_quarter_step201' in result else [('mesh401','mesh401',result['mesh401'],True)])
            if 'quarter_step201' in result:
                temporal=result['quarter_step201'];variants.append(('quarter_step201','quarter_step201',temporal['measurement'],temporal['no_cue_control']['no_cue_control_pass']))
            for variant,directory,measurement,control_pass in variants:
                companion=read(folder/directory/'readouts.json');phase=phase_endpoint_comparison(reference,companion,events)
                issues=[]
                event_check=(result['quarter_step201']['event_comparison'] if variant=='quarter_step201' else result['event_comparison'])
                if not event_check['agrees']:issues.append('replacement_crossing_timing_or_censoring_changed')
                if ((result['main201']['final_d'] is not None and result['main201']['final_d']>=.5)!=
                    (measurement['final_d'] is not None and measurement['final_d']>=.5)):
                    issues.append('final_replacement_fate_changed')
                if not control_pass:issues.append('companion_no_cue_control_failed')
                if result['main201']['field_tail_pass']!=measurement['field_tail_pass']:issues.append('finite_field_tail_status_changed')
                if measurement['maximum_endpoint_scaled_field_difference']>RULES['maximum_endpoint_scaled_field_difference']:
                    issues.append('native_endpoint_field_difference_above_tolerance')
                if not phase['all_endpoint_eligibility_agrees']:issues.append('phase_endpoint_eligibility_changed')
                if not phase['all_eligible_endpoint_D_within_tolerance']:issues.append('phase_endpoint_D_difference_above_tolerance')
                if not phase['all_available_challenger_bias_signs_agree']:issues.append('phase_endpoint_challenger_bias_changed')
                early=[]
                for i,endpoint in enumerate(phase['endpoints']):
                    common=dict(run=str(run),task=case['task'],model=case['model'],
                        requested_reference_length=case['requested_reference_length'],assay=case['assay'],variant=variant)
                    flat=dict(common,name=endpoint['name'],time=endpoint['time'],D_difference=endpoint['D_difference'],
                              eligibility_agrees=endpoint['eligibility_agrees'])
                    for label in ('main201','companion'):
                        for key,value in endpoint[label].items():flat[label+'_'+key]=value
                    phase_rows.append(flat)
                    if 'gap' in endpoint['name'] and len(early)<2:
                        loss=[]
                        for label in ('main201','companion'):
                            pulse=phase['endpoints'][i-1][label];gap=endpoint[label]
                            loss.append(None if not pulse['D_half_eligible'] or pulse['D_half']<0 else
                                        not gap['D_half_eligible'] or gap['D_half']<0)
                        early.append(dict(name=endpoint['name'],main_challenger_bias_lost=loss[0],
                            companion_challenger_bias_lost=loss[1],classification_agrees=loss[0]==loss[1]))
                if case['model']=='debelly':
                    for label,path in (('main201',prefix/'main201.npz'),(variant,folder/directory/'trace.npz')):
                        with np.load(path,allow_pickle=False) as archive:
                            values={key:archive['field_'+key][-1].copy() for key in DB_NAMES}
                            x=archive['x'].copy()
                        for j,position in enumerate(x):
                            profile_rows.append(dict(run=str(run),task=case['task'],requested_reference_length=case['requested_reference_length'],
                                assay=case['assay'],variant=variant,profile=label,cell=j,nominal_x=float(position),
                                **{key:float(value[j]) for key,value in values.items()}))
                records.append(dict(run=str(run),task=case['task'],model=case['model'],
                    requested_reference_length=case['requested_reference_length'],assay=case['assay'],variant=variant,
                    companion_no_cue_control_pass=control_pass,phase_comparison=phase,first_two_gap_bias_loss=early,
                    first_replacement_event_comparison=event_check,
                    main201_final=result['main201'],companion_final_D=measurement['final_d'],
                    companion_field_tail_pass=measurement['field_tail_pass'],
                    endpoint_native_fields=measurement['endpoint_native_fields'],
                    source_basal_field_shift=result.get('source_basal_field_shift',case.get('companion_basal_shift')),
                    identical_realized_phase_schedule=True,realized_events=events,
                    issues=issues,all_declared_readout_comparisons_agree=not issues,
                    native_worker_result_sha256=sha(folder/'result.json')))
    for filename,rows in (('phase_endpoints.csv',phase_rows),('debelly_native_endpoint_profiles.csv',profile_rows)):
        keys=sorted({key for row in rows for key in row})
        with (out/filename).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader();writer.writerows(rows)
    if sha(Path(__file__))!=source_hash:raise ValueError('Numerical summary source changed during execution')
    summary=dict(generated_utc=now(),analysis_source_sha256=source_hash,sources=sources,
        independent_native_identity_audits_passed=True,actual_protocol_comparisons=len(records),
        numerical_agreement_for_all_requested_readouts=all(row['all_declared_readout_comparisons_agree'] for row in records),
        unavailable_reference_cases=unavailable,records=records,
        outputs={name:sha(out/name) for name in ('phase_endpoints.csv','debelly_native_endpoint_profiles.csv')},
        interpretation='Identity audits are distinct from numerical agreement. Final fate agreement does not establish agreement in pulse/gap timing or coupled native fields.')
    write(out/'summary.json',summary);return dict(output=str(out),actual_protocol_comparisons=len(records),
        numerical_agreement_for_all_requested_readouts=summary['numerical_agreement_for_all_requested_readouts'],summary_sha256=sha(out/'summary.json'))


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for command in ('prepare-basals','prepare-temporal-basals'):
        prep=sub.add_parser(command);prep.add_argument('--out',type=Path,required=True);prep.add_argument('--basals',nargs='+',type=Path,required=True)
    for stage in ('width','train','temporal-width','temporal-train'):
        p=sub.add_parser('prepare-'+stage);p.add_argument('--out',type=Path,required=True);p.add_argument('--numerical-basals',type=Path,required=True);p.add_argument('--main-runs',nargs='+',type=Path,required=True)
        p.add_argument('--unavailable-trains',type=Path)
    for command in ('verify','collect','audit-basals','audit-protocols'):
        sub.add_parser(command).add_argument('--out',type=Path,required=True)
    work=sub.add_parser('worker');work.add_argument('--out',type=Path,required=True);work.add_argument('--task',type=int,required=True)
    report=sub.add_parser('summarize');report.add_argument('--out',type=Path,required=True);report.add_argument('--runs',type=Path,nargs='+',required=True)
    args=parser.parse_args()
    if args.command in ('prepare-basals','prepare-temporal-basals'):
        result=prepare_basals(args.out,args.basals,temporal=args.command=='prepare-temporal-basals')
    elif args.command.startswith('prepare-'):
        stage=args.command[8:];temporal=stage.startswith('temporal-')
        result=prepare_protocols(args.out,stage.removeprefix('temporal-'),args.numerical_basals,args.main_runs,
                                 temporal=temporal,unavailable_trains=args.unavailable_trains)
    elif args.command=='verify':m=verify(args.out);result=dict(verified=True,workers=m['workers'],stage=m['stage'])
    elif args.command=='collect':result=collect(args.out)
    elif args.command=='audit-basals':result=audit_basals(args.out)
    elif args.command=='audit-protocols':result=audit_protocols(args.out)
    elif args.command=='summarize':result=summarize(args.out,args.runs)
    else:
        try:result=worker(args.out,args.task)
        except Exception:
            write(args.out/'workers'/f'{args.task:04d}'/'failed.json',dict(generated_utc=now(),traceback=traceback.format_exc()));raise
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
