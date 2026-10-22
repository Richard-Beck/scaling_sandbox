"""Bind the three predeclared S11 contexts before selection, then run exports.

This controller performs read-only native identity/observation audits. It never
integrates a model, chooses a main example, or treats an observation companion
as an independent stochastic realization. Source inventories must already pass
the complete native/stage gates used by the Figure 4 canonical normalizer.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import fcntl
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'polarity1d/manuscript').is_dir())
RUNS=REPO/'polarity1d/manuscript/data/runs'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
DEFAULT=RUNS/'signature_figure4_native_atlas_final_20261020'
REQUIRED=('common','E1_extra','E4_long')
LOGICAL=REPO/'docs/tasks/briefs/20261020_figure4_s11_final_coverage_decisions.json'
IDENTITY=REPO/'docs/tasks/briefs/20261020_figure4_s11_identity_decisions.json'
BASE_RULES=REPO/'docs/tasks/briefs/20261020_figure4_s11_decisions.json'
MODULES=['signature_figure4_atlas'+suffix+'_20261020.py' for suffix in
         ['', '_audit','_pdf','_csv','_package','_text_audit','_finalize']]
FRONTEND='signature_figure4_atlas_frontend_audit_20261020.js'
LAW_KEYS=('generator','gaussian_kernel_sha256','physical_covariance_length',
          'covariance_definition','retained_spatial_zero_mode','knot_interval',
          'tau','rho','coupling_knot','interpolation')
SPRING_LAW_KEYS=('basal_C','basal_physical_ell','tau','knot','coupling_knot',
                 'kernel_sha256','bridge_base_interval','spatial_zero_mode_retained',
                 'geometry','interpolation')


def stamp():return datetime.now(timezone.utc).isoformat()
def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
    return h.hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise FileExistsError('Frozen identity/completion receipts are immutable: '+str(path))
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def verified(path,expected):
    path=Path(path).resolve()
    if not expected or sha(path)!=expected:raise ValueError('Source identity failed: '+str(path))
    return path


def verified_native_artifact(case,name):
    folder=Path(case['scalar_path']).parent;status=read(case['source_status_path'])
    expected=status.get('output_sha256',{}).get(name)
    if expected is None:expected=read(folder/'result.json')['artifact_sha256'][name]
    return verified(folder/name,expected)


def fullstate_identity(path,model,size):
    """Canonical numeric restart identity, independent of NPZ zip metadata.

    Every numerical checkpoint array is retained; in particular inactive
    species and the complete Spring state (including circumference) are used.
    Model/size, array names/shapes and canonical float64 numeric bytes bind the
    geometry. Container hashes and nonnumeric metadata remain provenance.
    """
    with np.load(path,allow_pickle=False) as z:
        if 'state' not in z:raise ValueError('A complete native state is required')
        state=z['state']
        if model=='spring_ring':
            if 'material_positions' not in z:raise ValueError('Spring material geometry is required')
            n=len(z['material_positions'])
            if state.shape!=(2*n+1,):raise ValueError('Spring restart must include both fields and circumference')
            required=('state','material_positions')
        else:
            species={'wave_pinning':2,'goryachev':2,'holmes_model3':4}.get(model)
            if species is None or 'x' not in z or state.shape!=(species,len(z['x'])):
                raise ValueError('Complete model-native species and ring grid are required')
            required=('state','x')
        arrays={}
        for key in sorted(z.files):
            value=z[key]
            if value.dtype.kind not in 'biuf':
                if value.dtype.kind in 'cO':raise ValueError('Undeclared complex/object restart variable')
                continue
            value=np.asarray(value,dtype='<f8',order='C')
            if not np.isfinite(value).all():raise ValueError('Finite complete initial state required')
            arrays[key]=dict(shape=list(value.shape),numeric_float64_sha256=hashlib.sha256(value.tobytes(order='C')).hexdigest())
        if not set(required)<=set(arrays):raise ValueError('Incomplete restart identity')
    payload=dict(schema='all_numeric_native_restart_arrays_v1',model=model,C_label=float(size),arrays=arrays)
    return dict(fullstate_and_geometry_sha256=digest(payload),**payload)


def forcing_identity(case,status,law):
    job=status.get('job',{})
    kind=case['kind'];noisy=kind=='noise'
    mean=job.get('mean',law.get('theoretical_input_mean',law.get('mean')))
    if mean is None or not np.isfinite(mean):raise ValueError('Declared native input mean required')
    if kind=='zero':mean=0.
    if noisy:
        keys=SPRING_LAW_KEYS if case.get('model')=='spring_ring' else LAW_KEYS
        if any(key not in law for key in keys):
            raise ValueError('Native Gaussian-law declaration incomplete: '+str([k for k in keys if k not in law]))
        if not law.get('gaussian_knots_sha256'):raise ValueError('An actual native Gaussian-knot hash is required')
        if int(law.get('seed',case['seed']))!=int(case['seed']):raise ValueError('Effective forcing seed disagrees with native receipt')
        if float(law.get('sigma',job.get('sigma',1.)))!=float(case['sigma']):raise ValueError('Forcing sigma disagrees with its native law/receipt')
    return dict(kind=kind,native_mean=float(mean),sigma=float(case['sigma']) if noisy else 0.,
        effective_stochastic_seed=int(case['seed']) if noisy else None,
        law={key:law[key] for key in keys} if noisy else dict(law='uniform_constant_native_input'))


def run_identities(case,initial,forcing,law):
    payload=dict(model=case['model'],C_label=float(case['C_label']),
        initial_fullstate_and_geometry_sha256=initial['fullstate_and_geometry_sha256'],
        forcing=forcing,native_horizon=float(case['native_horizon']),
        actual_gaussian_knots_sha256=law['gaussian_knots_sha256'] if case['kind']=='noise' else None)
    physical='physical_'+digest(payload);cohort=None
    if case['kind']=='noise':
        cohort='gaussian_'+digest(dict(model=case['model'],C_label=float(case['C_label']),f=float(case['f']),
            initial_fullstate_and_geometry_sha256=initial['fullstate_and_geometry_sha256'],
            effective_seed=int(case['seed']),native_horizon=float(case['native_horizon']),
            law=forcing['law'],actual_gaussian_knots_sha256=law['gaussian_knots_sha256']))
    return physical,cohort,payload


def observation_equivalence(original,fresh,*,rtol=0.,atol=0.):
    """Check all original numerical fields at shared native saved times.

    Scalars use exact equality. A profile/final-state replay may use only the
    explicit tolerance in its preexisting qualified native-equivalence receipt.
    Static grids are always exact. Nothing is interpolated.
    """
    records=[]
    with np.load(original,allow_pickle=False) as a,np.load(fresh,allow_pickle=False) as b:
        required={k for k in a.files if a[k].dtype.kind in 'biufc'}
        if not required<=set(b.files):raise ValueError('Companion dropped original numerical variables')
        indices=None
        if 'time' in required:
            if not np.array_equal(a['time'],b['time']):
                indices=np.searchsorted(b['time'],a['time'])
                if np.any(indices>=len(b['time'])) or not np.array_equal(b['time'][indices],a['time']):
                    raise ValueError('Companion does not retain every original native saved time')
        for key in sorted(required):
            x,y=a[key],b[key]
            if indices is not None and key not in ('x','material_positions') and x.ndim and x.shape[0]==len(a['time']):y=y[indices]
            exact=key in ('x','material_positions') or not (rtol or atol)
            if exact:np.testing.assert_array_equal(x,y)
            else:np.testing.assert_allclose(x,y,rtol=rtol,atol=atol,equal_nan=True)
            finite=np.isfinite(x)&np.isfinite(y)
            difference=float(np.max(abs(x[finite]-y[finite]))) if np.any(finite) and x.dtype.kind!='b' else 0.
            records.append(dict(field=key,shape=list(x.shape),bitwise_equal=x.dtype==y.dtype and x.tobytes()==y.tobytes(),
                matching_nonfinite_values=int((~finite).sum()),maximum_absolute_difference=difference,rtol=0. if exact else rtol,atol=0. if exact else atol))
    return records


def old_observation_identity(case,declaration,current_state):
    """Independently bind a declared old observation without feature loading."""
    old=declaration.get('original_case',declaration)
    if not isinstance(old,dict) or not all(old.get(k) for k in ('trace','profile')):
        raise ValueError('Observation/prefix declaration lacks original native archives')
    old_status_path=Path(declaration.get('original_status',Path(old['trace']).parent/'status.json'))
    old_status=read(old_status_path)
    if old_status.get('complete') is not True:raise ValueError('Original observation is not complete')
    if declaration.get('original_status_sha256'):verified(old_status_path,declaration['original_status_sha256'])
    job=old_status['job'];location=job.get('input',job.get('path'))
    if not location:raise ValueError('Native restart job has no declared checkpoint location')
    checkpoint=Path(location)
    if not checkpoint.is_absolute():checkpoint=old_status_path.parents[2]/checkpoint
    expected=job.get('input_sha256',job.get('sha256'))
    verified(checkpoint,expected)
    state=fullstate_identity(checkpoint,case['model'],case['C_label'])
    if state['fullstate_and_geometry_sha256']!=current_state['fullstate_and_geometry_sha256']:
        raise ValueError('Observation alias changed a complete basal state or geometry')
    source=[]
    for key,other in [('trace','scalar'),('profile','profile')]:
        verified(old[key],old[key+'_sha256'])
        source.append(dict(role=key,path=old[key],sha256=old[key+'_sha256']))
    law_path=Path(old['trace']).parent/'noise_audit.json';law=read(law_path)
    old_case=dict(case,sigma=float(job.get('sigma',1.)) if case['kind']=='noise' else 0.,seed=job['seed'],
        native_horizon=float(old.get('horizon',old_status['native_time'])))
    physical,cohort,payload=run_identities(old_case,state,forcing_identity(old_case,old_status,law),law)
    return dict(physical_run_id=physical,payload=payload,source_archives=source,
        original_status=dict(path=str(old_status_path),sha256=sha(old_status_path)),
        basal_checkpoint_container_sha256=expected,complete_initial_state_equal=True),old


def bind_case(case,state_cache):
    row=dict(case)
    for role in ('scalar','profile','source_status','basal_checkpoint','basal_profile'):
        verified(row[role+'_path'],row[role+'_sha256'])
    state_key=(row['basal_checkpoint_sha256'],row['model'],float(row['C_label']))
    if state_key not in state_cache:state_cache[state_key]=fullstate_identity(row['basal_checkpoint_path'],row['model'],row['C_label'])
    initial=state_cache[state_key];status=read(row['source_status_path'])
    if status.get('complete') is not True:raise ValueError('Complete native observation required')
    law_path=Path(row.get('noise_audit_path',Path(row['scalar_path']).parent/'noise_audit.json'))
    law=read(law_path)
    expected=row.get('noise_audit_sha256',status.get('output_sha256',{}).get('noise_audit.json'))
    if not expected:
        result_path=Path(row['scalar_path']).parent/'result.json'
        expected=read(result_path)['artifact_sha256']['noise_audit.json']
    verified(law_path,expected)
    forcing=forcing_identity(row,status,law)
    physical,cohort,payload=run_identities(row,initial,forcing,law)
    relation=dict(new_gaussian_path=row.get('new_gaussian_path'),independent_noise_seed=row.get('independent_noise_seed'),
        retained_prefix=row.get('retained_prefix'),paired_sigma1_source=row.get('paired_sigma1_source'))
    alias_proofs=[]
    companion=row.get('observation_companion_of')
    if companion:
        old_identity,old=old_observation_identity(row,companion,initial)
        # Extended control observations are explicitly different physical runs.
        extended=float(old_identity['payload']['native_horizon'])!=float(row['native_horizon'])
        proof=dict(original_identity=old_identity,extended_horizon=extended)
        if not extended:
            if old_identity['physical_run_id']!=physical:raise ValueError('Declared same-horizon companion changed physical forcing/state identity')
            proof['scalar_checks']=observation_equivalence(old['trace'],row['scalar_path'])
            replay_path=Path(row['scalar_path']).parent/'sigma1_evolution_equivalence.json'
            rtol=atol=0.
            if replay_path.exists():
                verified_native_artifact(row,'sigma1_evolution_equivalence.json')
                replay=read(replay_path)
                if replay.get('passed') is not True:raise ValueError('Native replay equivalence failed')
                rtol,atol=float(replay['rtol']),float(replay['atol'])
                proof['native_replay_receipt']=dict(path=str(replay_path),sha256=sha(replay_path),rtol=rtol,atol=atol)
            proof['profile_checks']=observation_equivalence(old['profile'],row['profile_path'],rtol=rtol,atol=atol)
        alias_proofs.append(proof)
    prefix=row.get('retained_prefix')
    if prefix:
        old_identity,old=old_observation_identity(row,prefix,initial)
        receipt_path=Path(row['scalar_path']).parent/'retained_native_prefix_equivalence.json'
        verified_native_artifact(row,'retained_native_prefix_equivalence.json')
        receipt=read(receipt_path)
        if receipt.get('passed') is not True or status.get('retained_native_prefix_equivalence_passed') is not True:
            raise ValueError('Long observation lacks its qualified full-state native-prefix proof')
        declared=row.get('native_metadata',{}).get('retained_native_prefix_equivalence')
        if declared:verified(receipt_path,declared['sha256'])
        relation.update(original_physical_run_id=old_identity['physical_run_id'],
            native_prefix_equivalence_receipt=dict(path=str(receipt_path),sha256=sha(receipt_path),
                rule=receipt['rule'],all_fields_bitwise_equal=all(c['bitwise_equal'] for c in receipt['checks'])),
            dependence='Same declared Gaussian prefix; longer physical run is not another independent seed lineage')
    if row.get('origin')=='E4' and row['kind']=='noise':
        if type(row.get('new_gaussian_path')) is not bool or (not row['new_gaussian_path'] and not prefix):
            raise ValueError('E4 requires explicit original-prefix versus independent top-up identity')
    row.update(source_observation_id=row['id'],id=physical,physical_run_id=physical,
        independent_realization_id=physical,latent_Gaussian_cohort_id=cohort,
        independent_stochastic_seed_count=int(row['kind']=='noise'),deterministic_control=row['kind']!='noise',
        prefix_and_topup_relation=relation)
    proof=dict(source_observation_id=row['source_observation_id'],physical_run_id=physical,
        physical_payload=payload,latent_Gaussian_cohort_id=cohort,
        basal_checkpoint_container_sha256=row['basal_checkpoint_sha256'],full_initial_state=initial,
        native_Gaussian_audit=dict(path=str(law_path),sha256=expected,actual_knots_sha256=law.get('gaussian_knots_sha256')),
        companion_identity_and_equivalence=alias_proofs,prefix_and_topup_relation=relation,
        selection_features_loaded=False)
    return row,proof


def validate_contexts(contexts):
    if set(contexts)!=set(REQUIRED):raise ValueError('Exactly the three frozen final contexts are required')
    counts={}
    for context,rows in contexts.items():
        groups=defaultdict(list)
        for r in rows:
            groups[(r['model'],float(r['C_label']),float(r['f']),r['kind'],float(r['sigma']),float(r['native_horizon']))].append(r)
        if any(len(v)!=1 for k,v in groups.items() if k[3]!='noise'):raise ValueError('Each deterministic condition has one control')
        if any(len({r['seed'] for r in v})!=len(v) or len(v)<2 for k,v in groups.items() if k[3]=='noise'):
            raise ValueError('Every noisy condition requires all declared distinct seeds and two selectable observations')
        if context in ('E1_extra','E4_long') and any(len(v)!=12 for k,v in groups.items() if k[3]=='noise'):
            raise ValueError('Every E1-extra/E4-long condition requires twelve stochastic realizations')
        counts[context]=dict(physical_runs=len(rows),groups=len(groups),group_counts=dict(Counter(k[3] for k in groups)),
            expected_pages=len({(r['model'],float(r['C_label']),float(r['f']),r['kind'],float(r['sigma']),float(r['native_horizon']),float(r['profile_cadence'])) for r in rows}),
            noisy_realizations=sum(r['kind']=='noise' for r in rows),deterministic_controls=sum(r['kind']!='noise' for r in rows))
        if context=='common' and any(r['sigma']!=1. for r in rows if r['kind']=='noise'):raise ValueError('Common uses sigma1 only')
        if context=='E1_extra':
            if set(r['kind'] for r in rows)!={'noise'} or set(r['sigma'] for r in rows)!={.25,.5,2.} or len(groups)!=72:
                raise ValueError('E1_extra requires exactly the72 noisy extra-sigma groups')
        if context=='E4_long' and counts[context]['group_counts']!={'noise':52,'mean':52,'zero':12}:
            raise ValueError('E4_long requires all52 noisy,52 mean and12 zero groups')
    return counts


def bind(run,stage_roots):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Complete native identity binding uses SLURM')
    from signature_figure4_regions_stage_20261020 import stage_rows,choose_arm
    if (run/'physical_identity_binding.json').exists():raise FileExistsError('Final physical map is already frozen')
    logical,identity=read(LOGICAL),read(IDENTITY)
    if tuple(c['atlas_context'] for c in logical['required_contexts'])!=REQUIRED:raise ValueError('Logical context contract changed')
    for path,expected in [(LOGICAL,read(run/'coverage_master.json')['logical_rules_sha256']),
                          (IDENTITY,read(run/'coverage_master.json')['identity_rules_sha256'])]:verified(path,expected)
    by_stage={};authorities=[]
    for root in stage_roots:
        stage,rows,authority=stage_rows(root)
        if stage in by_stage:raise ValueError('Exactly one canonical inventory per accepted stage')
        by_stage[stage]=rows;authorities.append(authority)
    if set(by_stage)!={'E1','E2','E3','E4'}:raise ValueError('All four complete normalized stages are required before final selection')
    common,replacements=choose_arm({k:by_stage[k] for k in ('E2','E3')},'common')
    contexts=dict(common=common,E1_extra=[r for r in by_stage['E1'] if r['kind']=='noise' and r['sigma'] in (.25,.5,2.)],E4_long=by_stage['E4'])
    counts=validate_contexts(contexts)
    mapped={};proofs=[];cache={};preferred={};duplicates=[]
    for context,rows in contexts.items():
        mapped[context]=[]
        for r in rows:
            row,proof=bind_case(r,cache);row['atlas_context']=context
            key=(context,row['physical_run_id'])
            if key in preferred:raise ValueError('Canonical input still duplicates a physical run inside a context')
            preferred[key]=row;mapped[context].append(row);proofs.append(dict(atlas_context=context,**proof))
    physical={r['physical_run_id'] for rows in mapped.values() for r in rows}
    if len(physical)!=sum(len(rows) for rows in mapped.values()):raise ValueError('Final contexts unexpectedly duplicate physical runs')
    # Every E1 sigma branch shares the verified actual Gaussian cohort with
    # its other sigma branches; common sigma1 supplies references, not pages.
    common_by_key={(r['model'],r['C_label'],r['f'],r['seed'],r['native_horizon']):r for r in mapped['common'] if r['kind']=='noise'}
    common_by_physical={r['physical_run_id']:r for r in mapped['common']}
    common_means={(r['model'],r['C_label'],r['f'],r['native_horizon']):r for r in mapped['common'] if r['kind']=='mean'}
    common_zeros={(r['model'],r['C_label'],r['native_horizon']):r for r in mapped['common'] if r['kind']=='zero'}
    references=[]
    for row in mapped['E1_extra']:
        key=(row['model'],row['C_label'],row['f'],row['seed'],row['native_horizon'])
        canonical=common_by_key.get(key)
        if canonical is None or canonical['latent_Gaussian_cohort_id']!=row['latent_Gaussian_cohort_id']:
            raise ValueError('Every E1 extra-sigma observation requires its exact common sigma1 latent-path reference')
        mean=common_means.get((row['model'],row['C_label'],row['f'],row['native_horizon']))
        zero=common_zeros.get((row['model'],row['C_label'],row['native_horizon']))
        if mean is None or zero is None:raise ValueError('Every E1 condition requires canonical common sigma0 mean and zero controls')
        references.append(dict(source_physical_run_id=row['physical_run_id'],common_sigma1_physical_run_id=canonical['physical_run_id'],
            common_sigma0_mean_physical_run_id=mean['physical_run_id'],common_zero_physical_run_id=zero['physical_run_id'],
            common_condition_key=dict(model=row['model'],C_label=row['C_label'],f=row['f'],native_horizon=row['native_horizon']),
            reference_scope='The common condition page displays its declared medoid/random; every paired seed remains in its candidate/native source inventory.',
            latent_Gaussian_cohort_id=row['latent_Gaussian_cohort_id']))
    for rows in mapped.values():
        for row in rows:
            lineage=row['latent_Gaussian_cohort_id']
            prefix_relation=row['prefix_and_topup_relation']
            if row.get('retained_prefix'):
                original=common_by_physical.get(prefix_relation['original_physical_run_id'])
                if original is None:raise ValueError('An E4 retained prefix has no canonical common physical-run reference')
                lineage=original['latent_Gaussian_cohort_id']
                prefix_relation['original_latent_Gaussian_cohort_id']=lineage
            row['latent_seed_lineage_id']=lineage
    proof_by_source={(p['atlas_context'],p['source_observation_id']):p for p in proofs}
    for context,rows in mapped.items():
        for row in rows:
            proof_by_source[(context,row['source_observation_id'])]['latent_seed_lineage_id']=row['latent_seed_lineage_id']
    for context,rows in mapped.items():
        write(run/'canonical_inputs'/context/'inventory.json',dict(atlas_context=context,cases=rows,
            source_stage_authorities=authorities,selection_features_loaded=False,source_case_count=len(rows)))
        write(run/'canonical_inputs'/context/'native_identity_verification.json',dict(passed=True,complete=True,cases=len(rows),
            inventory_sha256=sha(run/'canonical_inputs'/context/'inventory.json'),scope='Complete native/stage source gates, full numeric basal-state identity, forcing-law/path identities and declared observation equivalence before selection.'))
    result=dict(passed=True,complete=True,frozen_utc=stamp(),required_contexts=list(REQUIRED),source_stages=authorities,
        source_case_inventory_hashes={context:sha(run/'canonical_inputs'/context/'inventory.json') for context in REQUIRED},
        logical_rules_sha256=sha(LOGICAL),identity_rules_sha256=sha(IDENTITY),controller_source_sha256=sha(__file__),
        selection_features_loaded=False,native_integrations=0,main_example_selection=False,contexts=counts,
        exact_physical_identity_records=proofs,E1_common_reference_crosswalk=references,
        canonical_stage_observation_replacements=replacements,
        independence_and_reuse_census=dict(physical_runs=len(physical),native_observations=len(proofs),
            noisy_physical_runs=sum(r['kind']=='noise' for rows in mapped.values() for r in rows),
            deterministic_control_runs=sum(r['kind']!='noise' for rows in mapped.values() for r in rows),
            unique_latent_Gaussian_cohorts=len({r['latent_Gaussian_cohort_id'] for rows in mapped.values() for r in rows if r['kind']=='noise'}),
            unique_latent_seed_lineages=len({r['latent_seed_lineage_id'] for rows in mapped.values() for r in rows if r['kind']=='noise'}),
            E4_retained_prefix_runs=sum(bool(r.get('retained_prefix')) for r in mapped['E4_long']),
            E4_independent_topup_runs=sum(r.get('new_gaussian_path') is True for r in mapped['E4_long'] if r['kind']=='noise'),
            stochastic_seed_count_for_deterministic_controls=0))
    write(run/'physical_identity_binding.json',result)
    shutil.copy2(__file__,run/'frozen_identity_controller.py')
    refresh(run)
    return dict(passed=True,physical_runs=len(physical),contexts=counts,selection_features_loaded=False)


def prepare(run,context,metrics,memory,qos):
    from signature_figure4_atlas_20261020 import freeze
    binding=read(run/'physical_identity_binding.json')
    if binding.get('passed') is not True or binding['required_contexts']!=list(REQUIRED):raise ValueError('Complete frozen three-context identity binding is required')
    inputs=run/'canonical_inputs'/context
    verified(inputs/'inventory.json',binding['source_case_inventory_hashes'][context])
    target=run/'contexts'/context
    result=freeze(target,inputs/'inventory.json',inputs/'native_identity_verification.json',metrics,BASE_RULES,context,memory,qos)
    for source,name in [(run/'physical_identity_binding.json','physical_identity_binding.json'),
                        (LOGICAL,'logical_coverage_freeze.json'),(IDENTITY,'identity_binding_freeze.json')]:shutil.copy2(source,target/'inputs'/name)
    pipeline=target/'inputs/pipeline';pipeline.mkdir()
    for name in MODULES+[FRONTEND]:shutil.copy2(HERE/name,pipeline/name)
    plan=read(target/'manifest.json');by_model=defaultdict(list)
    for job in plan['jobs']:by_model[job['condition']['model']].append(job)
    if len(by_model)!=4:raise ValueError('Four representative model/resource pilots are required')
    quad=lambda j:2*int(round(j['condition']['native_horizon']/j['condition']['profile_cadence']))*max(c['native_cells'] for c in j['cases'])*min(2,len(j['cases']))
    pilots=sorted(max(jobs,key=lambda j:(quad(j),-j['slot']))['slot'] for jobs in by_model.values())
    write(target/'inputs/pipeline_freeze.json',dict(context=context,physical_binding_sha256=sha(run/'physical_identity_binding.json'),
        sources_sha256={name:sha(pipeline/name) for name in MODULES+[FRONTEND]},atlas_manifest_sha256=sha(target/'manifest.json'),
        pilot_slots=pilots,pilot_rule='Largest native quadrilateral count per model; all four must complete before the remaining production array.',
        identity_frozen_before_candidate_features=True,condition_metrics_path=str(metrics),condition_metrics_sha256=sha(metrics)))
    return dict(result,pilot_slots=pilots,context=context)


def sbatch(run,context,role,command,*,array=None,dependency=None,memory=None,time='04:00:00'):
    target=run/'contexts'/context;plan=read(target/'manifest.json');qos=plan['resources']['qos']
    args=['sbatch','--parsable','--partition=red','--qos='+qos,'--ntasks=1','--cpus-per-task=1',
        '--mem='+str(memory or plan['resources']['memory_gib'])+'G','--time='+time,
        '--job-name=F4-S11-'+context+'-'+role,'--output='+str(target/'logs'/(role+'_%A_%a.out')),
        '--error='+str(target/'logs'/(role+'_%A_%a.err'))]
    if array:args.append('--array='+','.join(str(n) for n in array))
    if dependency:args.append('--dependency=afterok:'+':'.join(str(n) for n in dependency))
    args.extend(['--wrap','exec env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '+command])
    receipt=target/'inputs'/('submission_'+role+'.json')
    if receipt.exists():
        previous=read(receipt)
        if previous['args']!=args:raise ValueError('An existing producer has a different frozen payload')
        return previous['job_id']
    identifier=subprocess.check_output(args,text=True).strip().split(';')[0]
    write(receipt,dict(submitted_utc=stamp(),job_id=identifier,args=args,role=role,context=context))
    return identifier


def production_resource_choice(pilot_memory,pilot_records,requested=None):
    """Execution-only sizing from the largest measured native render pilots."""
    peaks=[r['resource']['peak_resident_memory_KiB'] for r in pilot_records]
    if not peaks or any(not np.isfinite(p) or p<=0 for p in peaks):
        raise ValueError('Finite positive native pilot memory measurements are required')
    peak=max(peaks)
    minimum=max(2,math.ceil(1.5*peak/1024**2+1.))
    selected=int(pilot_memory if requested is None else requested)
    if selected<minimum:
        raise ValueError('Production memory must retain the measured pilot margin: '+str(minimum)+' GiB required')
    return dict(execution_only=True,scientific_inputs_and_native_values_unchanged=True,
        pilot_allocation_GiB=int(pilot_memory),maximum_measured_RSS_KiB=int(peak),
        margin_rule='At least 1.5 times maximum measured pilot RSS plus 1 GiB, rounded up; minimum 2 GiB.',
        minimum_with_margin_GiB=minimum,production_allocation_GiB=selected,
        explicit_post_pilot_override=requested is not None,
        worker_receipts=[dict(slot=r['slot'],worker_sha256=r['worker_sha256'],resource=r['resource']) for r in pilot_records])


def launch(run,context,pilots_only=False,production_memory=None):
    target=run/'contexts'/context;plan=read(target/'manifest.json');pipeline=target/'inputs/pipeline';freeze=read(target/'inputs/pipeline_freeze.json')
    verified(target/'manifest.json',freeze['atlas_manifest_sha256'])
    for name,value in freeze['sources_sha256'].items():verified(pipeline/name,value)
    q=shlex.quote;py=q(PYTHON)+' -B '
    atlas=py+q(str(target/'inputs/frozen_atlas.py'))
    if pilots_only:
        if production_memory is not None:raise ValueError('Production resources are chosen after the native pilots complete')
        return dict(pilots=sbatch(run,context,'pilots',atlas+' worker --run '+q(str(target))+' --slot "$SLURM_ARRAY_TASK_ID"',array=freeze['pilot_slots']))
    pilot_records=[]
    for slot in freeze['pilot_slots']:
        path=target/'workers'/f'{slot:04d}.json';record=read(path)
        if record.get('complete') is not True or record['condition']!=plan['jobs'][slot]['condition']:raise ValueError('Every representative resource pilot must pass first')
        if not record.get('resource') or record['resource']['peak_resident_memory_KiB']>=plan['resources']['memory_gib']*1024**2:
            raise ValueError('Representative pilot must report successful execution within the declared memory allocation')
        pilot_records.append(dict(slot=slot,worker_sha256=sha(path),resource=record.get('resource'),passed=True))
    pilot_gate=dict(passed=True,records=pilot_records,declared_memory_gib=plan['resources']['memory_gib'])
    pilot_path=target/'inputs/resource_pilot_gate.json'
    if pilot_path.exists():
        if read(pilot_path)!=pilot_gate:raise ValueError('Frozen resource pilot evidence changed')
    else:write(pilot_path,pilot_gate)
    resource_path=target/'inputs/production_resource_receipt.json'
    if resource_path.exists() and production_memory is None:
        previous=read(resource_path)
        if previous.get('explicit_post_pilot_override'):production_memory=previous['production_allocation_GiB']
    resources=production_resource_choice(plan['resources']['memory_gib'],pilot_records,production_memory)
    if resource_path.exists():
        if read(resource_path)!=resources:raise ValueError('Frozen post-pilot production resources changed')
    else:write(resource_path,resources)
    jobs={};remaining=[j['slot'] for j in plan['jobs'] if j['slot'] not in freeze['pilot_slots']]
    if remaining:jobs['render']=sbatch(run,context,'render',atlas+' worker --run '+q(str(target))+' --slot "$SLURM_ARRAY_TASK_ID"',array=remaining,memory=resources['production_allocation_GiB'])
    jobs['collect']=sbatch(run,context,'collect',atlas+' collect --run '+q(str(target)),dependency=[jobs['render']] if remaining else None,memory=4)
    for role,suffix in [('native_audit','audit'),('pdf','pdf'),('text_audit','text_audit')]:
        jobs[role]=sbatch(run,context,role,py+q(str(pipeline/('signature_figure4_atlas_'+suffix+'_20261020.py')))+' --run '+q(str(target)),dependency=[jobs['collect']],memory=2 if role=='text_audit' else 4)
    csv=py+q(str(pipeline/'signature_figure4_atlas_csv_20261020.py'))
    jobs['csv']=sbatch(run,context,'csv',csv+' worker --run '+q(str(target))+' --slot "$SLURM_ARRAY_TASK_ID"',array=list(range(len(plan['jobs']))),dependency=[jobs['collect']],memory=4)
    jobs['csv_collect']=sbatch(run,context,'csv_collect',csv+' collect --run '+q(str(target)),dependency=[jobs['csv'],jobs['pdf']],memory=4)
    controller=py+q(str(run/'frozen_identity_controller.py'))
    jobs['package']=sbatch(run,context,'package',controller+' package --run '+q(str(run))+' --context '+context,
        dependency=[jobs['native_audit'],jobs['text_audit'],jobs['csv_collect']],memory=4)
    return jobs


def start(run,stage_roots,metrics_plan,memory,qos):
    """One accepted-input launch; all subsequent steps follow SLURM gates."""
    specification=read(metrics_plan)
    if set(specification['contexts'])!=set(REQUIRED):raise ValueError('All three accepted forced-condition sources must be explicit')
    for context,source in specification['contexts'].items():
        verified(source['path'],source['sha256']);verified(source['status_path'],source['status_sha256'])
        gate=read(source['status_path'])
        if gate.get('passed') is not True or gate.get('complete') is not True:raise ValueError('Accepted complete condition-source gate required')
    if not (run/'physical_identity_binding.json').exists():bind(run,stage_roots)
    binding=read(run/'physical_identity_binding.json')
    if binding['required_contexts']!=list(REQUIRED) or binding.get('passed') is not True:raise ValueError('Final native identity binding failed')
    plan_path=run/'accepted_condition_sources.json'
    if plan_path.exists():
        if read(plan_path)!=specification:raise ValueError('Accepted condition sources changed after selection preparation')
    else:write(plan_path,specification)
    jobs={};q=shlex.quote
    for context in REQUIRED:
        metrics=Path(specification['contexts'][context]['path']).resolve();target=run/'contexts'/context
        if not target.exists():prepare(run,context,metrics,memory,qos)
        else:
            frozen=read(target/'inputs/pipeline_freeze.json')
            if frozen['physical_binding_sha256']!=sha(run/'physical_identity_binding.json') or frozen['condition_metrics_sha256']!=sha(metrics):
                raise ValueError('Existing prepared context has different native or forcing-classification sources')
        pilot=launch(run,context,True)['pilots']
        dispatch=q(PYTHON)+' -B '+q(str(run/'frozen_identity_controller.py'))+' launch --run '+q(str(run))+' --context '+context
        dispatcher=sbatch(run,context,'dispatch',dispatch,dependency=[pilot],memory=2,time='00:10:00')
        jobs[context]=dict(resource_pilots=pilot,dependent_production_dispatch=dispatcher)
    return dict(submitted=True,contexts=jobs,physical_identity_binding_sha256=sha(run/'physical_identity_binding.json'),
        native_integrations=0,main_example_selection=False,completion_master=str(run/'coverage_master.json'))


def package(run,context):
    target=run/'contexts'/context;pipeline=target/'inputs/pipeline';portable=target/'results/portable_export'
    for suffix in ('package',):subprocess.run([PYTHON,'-B',str(pipeline/f'signature_figure4_atlas_{suffix}_20261020.py'),'--run',str(target)],check=True)
    subprocess.run(['node',str(pipeline/FRONTEND),str(portable)],check=True)
    subprocess.run([PYTHON,'-B',str(pipeline/'signature_figure4_atlas_finalize_20261020.py'),'--run',str(target),'--package',str(portable)],check=True)
    return refresh(run)


def refresh(run):
    # Context package jobs can finish concurrently. Serialize the source-bound
    # master snapshot and immutable final gate, including the final SHA binding.
    with (run/'.coverage_master.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _refresh(run)


def _refresh(run):
    binding=read(run/'physical_identity_binding.json');indices=[];incomplete=[]
    for context in REQUIRED:
        package=run/'contexts'/context/'results/portable_export';gate_path=package/'completion_gate.json'
        if not gate_path.exists():incomplete.append(context);continue
        gate=read(gate_path)
        if gate.get('passed') is not True or gate['analysis_arm']!=context:raise ValueError('A required context completion gate failed')
        expected=binding['contexts'][context]
        if gate['groups']!=expected['expected_pages'] or gate['independent_candidates']!=expected['physical_runs']:
            raise ValueError('A required final context omits source pages or candidate physical runs')
        inventory_path=run/'canonical_inputs'/context/'inventory.json'
        verified(inventory_path,binding['source_case_inventory_hashes'][context])
        required_cases=read(inventory_path)['cases']
        plan=read(package/'metadata/source_manifest.json')
        if sorted(c['physical_run_id'] for j in plan['jobs'] for c in j['cases'])!=sorted(c['physical_run_id'] for c in required_cases):
            raise ValueError('Rendered candidate census differs from its preselection canonical binding')
        for name,value in gate['final_portable_assets_sha256'].items():verified(package/name,value)
        index_path=package/'page_data/page_data_index.json';verified(index_path,gate['native_CSV_page_index_sha256']);index=read(index_path)
        # Bind every page, PDF volume, HPC CSV shard and portable companion.
        for page in index['pages']:
            for part in page['native_field_CSV_shards']:verified(part['path'],part['sha256'])
            for item in [page['selection_CSV'],*page['companions'].values()]:verified(package/'page_data'/Path(item['path']).name,item['sha256'])
        volumes=read(package/'pdf/lossless_pdf_audit.json')['volumes']
        for volume in volumes:verified(package/'pdf'/volume['file'],volume['sha256'])
        indices.append(dict(atlas_context=context,path=str(index_path),sha256=sha(index_path),pages=index['summary']['pages'],
            completion_gate=dict(path=str(gate_path),sha256=sha(gate_path)),volumes=volumes,
            required_case_inventory_sha256=binding['source_case_inventory_hashes'][context],
            page_and_shard_map=index['pages']))
    long_inventory=read(run/'canonical_inputs/E4_long/inventory.json')['cases']
    long_conditions={(r['model'],r['C_label'],r['f']) for r in long_inventory if r['kind']=='noise'}
    common_inventory=read(run/'canonical_inputs/common/inventory.json')['cases']
    selected_noisy=[r for r in common_inventory if r['kind']=='noise' and (r['model'],r['C_label'],r['f']) not in long_conditions]+[r for r in long_inventory if r['kind']=='noise']
    zero_needed={(r['model'],r['C_label'],r['native_horizon']) for r in selected_noisy}
    longest_pages=[]
    for entry in indices:
        if entry['atlas_context'] not in ('common','E4_long'):continue
        for page in entry['page_and_shard_map']:
            c=page['condition']
            included=(c['model'],c['C_label'],c['native_horizon']) in zero_needed if c['kind']=='zero' else entry['atlas_context']=='E4_long' or (c['model'],c['C_label'],c['f']) not in long_conditions
            if included:longest_pages.append(dict(atlas_context=entry['atlas_context'],pdf_volume=page['pdf_volume'],pdf_page=page['pdf_page'],condition=c,selected_cases=page['selected_cases']))
    result=dict(generated_utc=stamp(),passed=not incomplete,all_required_contexts_complete=not incomplete,
        status='complete' if not incomplete else 'required_contexts_running_or_pending',required_contexts=list(REQUIRED),
        selection_features_loaded=any((run/'contexts'/c/'workers').exists() and any((run/'contexts'/c/'workers').glob('*.json')) for c in REQUIRED),
        incomplete_contexts=incomplete,page_data_indices=indices,source_case_inventory_hashes=binding['source_case_inventory_hashes'],
        physical_identity_binding_sha256=sha(run/'physical_identity_binding.json'),logical_rules_sha256=sha(LOGICAL),identity_rules_sha256=sha(IDENTITY),
        independence_and_reuse_census=binding['independence_and_reuse_census'],
        longest_view='References common and E4_long physical pages; never an additional independent census.',
        longest_page_selector=longest_pages,
        E1_common_reference_crosswalk=binding['E1_common_reference_crosswalk'],prototype_counts_toward_required_coverage=False,
        native_integrations=0,main_example_selection=False)
    # Progress master is replaceable until complete; a final immutable gate is
    # written once, so concurrent context completions cannot hide a finished arm.
    if (run/'complete.json').exists():return read(run/'complete.json')
    temp=run/f'coverage_master.{os.getpid()}.tmp';temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temp.replace(run/'coverage_master.json')
    if not incomplete:write(run/'complete.json',dict(result,coverage_master_sha256=sha(run/'coverage_master.json')))
    return dict(passed=result['passed'],required_contexts_complete=len(indices),incomplete_contexts=incomplete)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['bind','prepare','pilots','launch','package','refresh','start'])
    p.add_argument('--run',type=Path,default=DEFAULT);p.add_argument('--stage-analyses',type=Path,nargs='+')
    p.add_argument('--context',choices=REQUIRED);p.add_argument('--condition-metrics',type=Path)
    p.add_argument('--metrics-plan',type=Path)
    p.add_argument('--production-memory-gib',type=int,help='Execution-only render allocation after native resource pilots and their review.')
    p.add_argument('--memory-gib',type=int,default=16);p.add_argument('--qos',default='small');args=p.parse_args();run=args.run.resolve()
    if not run.is_relative_to(RUNS.resolve()):raise ValueError('Use manuscript native run storage')
    if args.action=='start':
        if not args.metrics_plan:p.error('start requires the accepted three-context condition-source receipt')
        result=start(run,args.stage_analyses or [],args.metrics_plan.resolve(),args.memory_gib,args.qos)
    elif args.action=='bind':result=bind(run,args.stage_analyses or [])
    elif args.action=='prepare':
        if not args.context or not args.condition_metrics:p.error('prepare requires context and accepted condition metrics')
        result=prepare(run,args.context,args.condition_metrics.resolve(),args.memory_gib,args.qos)
    elif args.action in ('pilots','launch','package'):
        if not args.context:p.error('This action requires a frozen context')
        result=package(run,args.context) if args.action=='package' else launch(run,args.context,args.action=='pilots',args.production_memory_gib)
    else:result=refresh(run)
    print(json.dumps(result))


if __name__=='__main__':main()
