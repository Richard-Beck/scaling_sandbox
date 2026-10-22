"""Future E1--E4 A1--A4/A7 adapters, frozen estimands and independent audits.

Preparation and verification are gated by delivered corrected E0 artifacts.
No model integration is implemented here. E1 sigma arms share whole-seed
bootstrap weights; long E4 observations are never pooled into E0 windows.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import numpy as np
import signature_figure4_analysis_20261020 as A
import signature_figure4_audit_20261020 as Q
import signature_figure4_e0_20261020 as E
import signature_figure4_inventory_20261020 as I


def case_key(row):
    return (row['model'],float(row['C_label']),float(row['f']),row['kind'],
            row['seed'] if row['kind']=='noise' else None,float(row['sigma']),float(row['native_horizon']))


def scalar_identity(first,second):
    with np.load(first['scalar_path'],allow_pickle=False) as a, np.load(second['scalar_path'],allow_pickle=False) as b:
        required={'time','P_real','P_imag','contrast','circumference','current_half_contrast_regions'}
        if not required<=set(a.files)&set(b.files):raise ValueError('An observation companion lacks required native scalar fields')
        for key in set(a.files)&set(b.files):
            numeric=a[key].dtype.kind in 'biufc' and b[key].dtype.kind in 'biufc'
            if not np.array_equal(a[key],b[key],equal_nan=numeric):
                raise ValueError('Declared observation companions have different native scalar samples: '+key)


def deduplicate(rows):
    bykey={};replacements=[]
    for row in rows:
        key=case_key(row)
        if key not in bykey:
            bykey[key]=row;continue
        prior=bykey[key]
        samefiles=(prior['scalar_sha256'],prior['profile_sha256'])==(row['scalar_sha256'],row['profile_sha256'])
        if samefiles:continue
        if row['kind']=='noise' and not (row.get('observation_companion_of') or prior.get('observation_companion_of')):
            raise ValueError('Duplicate independent noisy identity without an explicit observation-companion declaration')
        scalar_identity(prior,row)
        chosen=row if row['profile_cadence']<prior['profile_cadence'] else prior
        bykey[key]=chosen
        replacements.append(dict(native_key=key,chosen_id=chosen['id'],reason='Exact scalar identity; prefer finer saved profiles, never count another independent seed'))
    return sorted(bykey.values(),key=lambda r:case_key(r)[:3]+(r['kind'],str(r['seed']),r['sigma'],r['native_horizon'])),replacements


def select_stage(rows,stage):
    if stage=='E1':
        keys={(r['model'],r['C_label'],r['f'],r['native_horizon']) for r in rows if r.get('origin')=='E1' and r['kind']=='noise'}
        selected=[r for r in rows if (r['kind']=='noise' and (r['model'],r['C_label'],r['f'],r['native_horizon']) in keys)
                  or (r['kind']=='mean' and (r['model'],r['C_label'],r['f'],r['native_horizon']) in keys)
                  or (r['kind']=='zero' and any((r['model'],r['C_label'],r['native_horizon'])==(m,C,T) for m,C,f,T in keys))]
        for key in keys:
            seeds=[]
            latent=[]
            for sigma in (.25,.5,1.,2.):
                members=[r for r in selected if r['kind']=='noise' and (r['model'],r['C_label'],r['f'],r['native_horizon'])==key and r['sigma']==sigma]
                members=sorted(members,key=lambda r:r['seed'])
                seeds.append([r['seed'] for r in members])
                latent.append([r.get('latent_gaussian_sha256') for r in members])
            if not all(len(s)==12 and s==seeds[0] for s in seeds):
                raise ValueError('E1 requires twelve exact common latent-seed blocks at sigma .25/.5/1/2, including E2 companions where needed')
            if any(not h for h in latent[0]) or not all(h==latent[0] for h in latent):
                raise ValueError('Matching integer seeds are insufficient: exact latent Gaussian-knot hashes must agree across E1 sigma')
    elif stage in ('E2','E3'):
        selected=[r for r in rows if r.get('origin') in ('E0','core','refinement',stage)]
    elif stage=='E4':
        selected=[r for r in rows if r.get('origin')=='E4']
    else:raise ValueError('Only E1--E4 are supported')
    if not any(r['kind']=='noise' and r.get('origin')==stage for r in selected):
        raise ValueError('The requested native stage has no completed noisy cases')
    return selected


def condition_groups(rows,stage,sigma):
    noisy=[r for r in rows if r['kind']=='noise' and r['sigma']==sigma]
    keys=sorted({(r['model'],r['C_label'],r['f'],r['native_horizon']) for r in noisy})
    if len({k[:3] for k in keys})!=len(keys):raise ValueError('Do not pool different observation horizons in one arm')
    groups=[];used={}
    for model,C,f,T in keys:
        members=sorted([r for r in noisy if (r['model'],r['C_label'],r['f'],r['native_horizon'])==(model,C,f,T)],key=lambda r:r['seed'])
        means=[r for r in rows if (r['model'],r['C_label'],r['f'],r['native_horizon'],r['kind'])==(model,C,f,T,'mean')]
        zero=[r for r in rows if (r['model'],r['C_label'],r['native_horizon'],r['kind'])==(model,C,T,'zero')]
        if len(means)!=1 or len(zero)!=1 or len(members) not in (2,12):
            raise ValueError('A condition lacks exact horizon-matched controls or independent seeds')
        if len({r['scalar_cadence'] for r in members})!=1 or len({r['profile_cadence'] for r in members})!=1:
            raise ValueError('One condition must have one explicit scalar/profile cadence; do not mix observation grids')
        key=f'{stage}:{model}:{C}:{f}' if stage=='E1' else f'{model}:{C}:{f}'
        groups.append(dict(task=len(groups),model=model,C_label=C,f=f,sigma=sigma,stage=stage,
            noise=members,mean=means[0],zero=zero[0],bootstrap_pairing_key=key))
        for r in members+means+zero:used[r['id']]=r
    return groups,list(used.values())


def abnormal_persistence(metrics):
    """Root-approved observed-control E1 rule, never a sigma-to-zero limit."""
    required=(0.,.25,.5,1.)
    if not all(s in metrics for s in required):return dict(status='unavailable_missing_sigma_arm',forced_E1=None)
    label=metrics[1.]['regime']
    persists=label in ('multiple_fronts','low_contrast') and all(metrics[s]['regime']==label for s in (.25,.5))
    flag=persists and metrics[0.]['regime']==label
    return dict(status='evaluated_abnormal_majority_with_observed_mean',forced_E1=flag,
                low_sigma_abnormal_persistence=persists,mean_agrees_with_sigma1_abnormal=label in ('multiple_fronts','low_contrast') and metrics[0.]['regime']==label,
                abnormal_regime=label if persists else None,
                rule='Same abnormal majority at sigma1,.5,.25 and matched sigma0 mean; healthy single_front never triggers')


def release_check(path):
    gate=E.read_json(path)
    if gate.get('e0_interim_delivered') is not True:
        raise ValueError('Deliver the corrected E0 interim before E1/E2 analysis')
    artifacts=gate.get('interim_artifacts',gate.get('artifacts',{}))
    if not artifacts:raise ValueError('Concrete delivered E0 artifacts and hashes are required')
    for artifact,digest in artifacts.items():E.require_hash(artifact,digest)
    return gate


def native_rows(run,out):
    manifest=E.read_json(run/'manifest.json');status=E.read_json(run/'results/status.json')
    census=E.read_json(run/'results/case_census.json')
    if not status.get('complete') or status.get('missing_slots') or status.get('failed_slots'):
        raise ValueError('Missing native work is not censoring')
    if census['native_run_manifest_sha256']!=E.digest(run/'manifest.json'):
        raise ValueError('Native census/manifest mismatch')
    total=len(manifest['jobs'])
    if len(census['cases'])!=total or status.get('completed')!=total or status.get('total')!=total:
        raise ValueError('Every declared native job must appear in the completed census')
    basals={(b['model'],float(b['size'])):b for b in manifest['basals']};rows=[]
    jobs={int(j['slot']):j for j in manifest['jobs']};seen_slots=set()
    for case in census['cases']:
        slot=int(Path(case['trace']).parent.name)
        if slot not in jobs or slot in seen_slots:raise ValueError('Native case slots must identify every frozen job exactly once')
        seen_slots.add(slot);job=jobs[slot]
        if any(case[k]!=job[k] for k in ('model','size','fraction','kind','seed','sigma')):
            raise ValueError('Native census labels disagree with their exact frozen job')
        for a,b in (('horizon','horizon'),('scalar_interval','scalar_interval'),('profile_interval','profile_interval')):
            if case[a]!=job['protocol'][b]:raise ValueError('Native observation schedule disagrees with its frozen job')
        E.require_hash(Path(case['trace']).parent/'result.json',case['native_result_sha256'])
        model=case['model'];C=float(case['size']);b=basals[(model,C)]
        if case['basal_input_sha256']!=b['input_sha256']:raise ValueError('Native case uses another basal input')
        checkpoint=(run/b['input']).resolve();E.require_hash(checkpoint,b['input_sha256'])
        active,q,actualC=I.basal_profile(checkpoint,model,C)
        folder=out/'inputs/basals';folder.mkdir(parents=True,exist_ok=True)
        compact=folder/f'{model}_C{C:g}_{b["input_sha256"][:16]}.npz'
        if not compact.exists():np.savez_compressed(compact,primary=active,material_positions=q,circumference=actualC,source_checkpoint_sha256=np.array(b['input_sha256']))
        native_status=Path(case['trace']).parent/'status.json'
        if not E.read_json(native_status).get('complete'):raise ValueError('Native case is not terminal-complete')
        rows.append(dict(id=run.name+'__'+case['case_id'],model=model,C_label=C,f=float(case['fraction']),
            kind='zero' if case['kind']=='basal' else case['kind'],seed=case['seed'],sigma=float(case['sigma']),origin=manifest['stage'],
            scalar_path=case['trace'],scalar_sha256=case['trace_sha256'],profile_path=case['profile'],profile_sha256=case['profile_sha256'],
            source_status_path=str(native_status),source_status_sha256=E.digest(native_status),source_manifest_sha256=E.digest(run/'manifest.json'),
            scalar_cadence=float(case['scalar_interval']),profile_cadence=float(case['profile_interval']),native_horizon=float(case['horizon']),
            basal_checkpoint_path=str(checkpoint),basal_checkpoint_sha256=b['input_sha256'],basal_profile_path=str(compact),basal_profile_sha256=E.digest(compact),
            basal_contrast=float(np.ptp(active)),actualC0=actualC,native_cells=len(active),
            observation_companion_of=case.get('observation_companion_of'),independent_noise_seed=case.get('independent_noise_seed'),
            new_gaussian_path=case.get('new_gaussian_path'),retained_prefix=case.get('retained_prefix'),
            paired_sigma1_source=case.get('paired_sigma1_source'),
            native_metadata={k:v for k,v in case.items() if k not in ('trace','profile')}))
    return rows


def require_complete_native_gate(run):
    """Engineering identity only; no tau, regime or hypothesis is inspected."""
    run=Path(run).resolve();path=run/'results/complete_native_independent_audit_20261020.json'
    manifest=E.read_json(run/'manifest.json');gate=E.read_json(path)
    if not gate.get('passed') or not gate.get('full_native_run_covered') or gate.get('missing_tasks') or gate.get('failed_tasks') or gate.get('native_manifest_sha256')!=E.digest(run/'manifest.json'):
        raise ValueError('A complete independent native engineering gate is required before response analysis')
    rows=gate.get('receipts',[]);slots=[int(row['slot']) for row in rows]
    if len(slots)!=len(manifest['jobs']) or set(slots)!=set(range(len(manifest['jobs']))):raise ValueError('The independent native gate does not cover every frozen worker exactly once')
    for receipt in rows:
        p=Path(receipt['path']);E.require_hash(p,receipt['sha256']);row=E.read_json(p);slot=int(receipt['slot']);folder=run/f'workers/{slot:04d}'
        if receipt.get('passed') is not True or row.get('passed') is not True or row.get('native_manifest_sha256')!=E.digest(run/'manifest.json') or row.get('slot')!=slot:raise ValueError('The independent native worker audit belongs to different inputs')
        for name,key in [('result.json','native_result_sha256'),('status.json','native_status_sha256')]:E.require_hash(folder/name,row[key])
        for filename,snapshot in row.get('hashed_native_file_snapshots',{}).items():
            info=Path(filename).stat();actual=dict(inode=info.st_ino,size=info.st_size,mtime_ns=info.st_mtime_ns,ctime_ns=info.st_ctime_ns)
            if actual!=snapshot:raise ValueError('An independently hashed native file changed after engineering qualification')
    return dict(path=str(path),sha256=E.digest(path),native_manifest_sha256=E.digest(run/'manifest.json'),workers=len(slots))


def declared_coverage(rows,stage,base,selection=None,triggers=None):
    observed={(r['model'],r['C_label'],r['f'],r['native_horizon']) for r in rows if r['kind']=='noise' and r.get('origin')==stage}
    basals={(r['model'],r['C_label'],r['native_horizon']) for r in base if r['kind']=='zero'}
    if stage=='E1':
        expected={(model,C,f,T) for model,C,T in basals for f in selection['e1_strengths'][model]}
    elif stage=='E2':
        expected={(model,C,f,T) for model,C,T in basals for f in (.005,.02,.05)}
    elif stage=='E3':
        horizons={model:T for model,C,T in basals}
        expected={(model,C,f,T) for model,T in horizons.items() for C in ((1.5,2.5) if model=='spring_ring' else (30.,50.)) for f in (.001,.003,.01,.03,.1)}
    elif stage=='E4':
        expected={(r['model'],float(r['C_label']),float(r['f']),4*float(r['original_horizon'])) for r in triggers if str(r['trigger'])=='True'}
    else:raise ValueError('Unknown future stage')
    if observed!=expected:raise ValueError(f'{stage} full declared condition coverage is incomplete: missing {sorted(expected-observed)}, unexpected {sorted(observed-expected)}')


def prepare(out,stage,e0,runs,gate,forced_rule=None):
    release=release_check(gate)
    if not out.is_relative_to((I.PACKAGE/'data/runs').resolve()):raise ValueError('Prepare future analysis under manuscript/data/runs')
    if out.exists():raise FileExistsError('Use a fresh future-stage analysis directory')
    if release.get('corrected_basal_inventory_sha256'):
        E.require_hash(e0/'inventory.json',release['corrected_basal_inventory_sha256'])
    if (e0/'SUPERSEDED.json').exists():raise ValueError('A superseded E0 inventory is not eligible')
    corrected=E.read_json(e0/'inventory.json')
    if len(corrected.get('basal_native_initial_checks',[]))!=12 or not all(r['exact_native_initial_primary'] for r in corrected['basal_native_initial_checks']):
        raise ValueError('The corrected E0 inventory requires all twelve exact native basal-profile checks')
    if stage=='E1' and not (forced_rule and forced_rule.get('root_approved') is True and forced_rule.get('rule')=='same_abnormal_majority_at_sigma1_both_lower_and_mean'):
        raise ValueError('Freeze the root-approved observed-mean E1 forced-state rule before analysis')
    native_qualification={str(run):require_complete_native_gate(run) for run in runs}
    out.mkdir(parents=True);(out/'results').mkdir();(out/'inputs').mkdir()
    base=E.read_json(e0/'inventory.json')
    rows=[dict(r,origin='E0') for r in base['cases']]
    for run in runs:rows.extend(native_rows(run,out))
    authority={};selection=None;triggers=None
    if stage=='E1':
        path=I.ROOT/'docs/tasks/briefs/20261020_figure4_e1_reviewed_selection.json';selection=E.read_json(path);authority[str(path)]=E.digest(path)
    if stage=='E4':
        path=e0/'results/interim_handoff/e4_trigger_inventory.csv';triggers=list(csv.DictReader(path.open()));authority[str(path)]=E.digest(path)
    declared_coverage(rows,stage,base['cases'],selection,triggers)
    for row in rows:
        if row['kind']!='noise':continue
        path=Path(row['scalar_path']).parent/'noise_audit.json'
        if path.is_file():
            noise=E.read_json(path)
            row.update(noise_audit_path=str(path),noise_audit_sha256=E.digest(path),latent_gaussian_sha256=noise.get('gaussian_knots_sha256'))
    rows,replacements=deduplicate(rows)
    selected=select_stage(rows,stage)
    newgroups=Counter((r['model'],r['C_label'],r['f'],r['sigma'],r['native_horizon']) for r in selected if r['kind']=='noise' and r.get('origin')==stage)
    if any(n!=12 for n in newgroups.values()):raise ValueError('Every new future-stage condition requires twelve native seed outcomes')
    inventory=dict(created_utc=E.stamp(),stage=stage,cases=selected,counts=dict(Counter(r['kind'] for r in selected)),
        no_native_integrations=True,observation_replacements=replacements,source_censuses={str(run/'results/case_census.json'):E.digest(run/'results/case_census.json') for run in runs},
        corrected_E0_inventory_sha256=E.digest(e0/'inventory.json'))
    E.write_json(out/'inventory.json',inventory)
    frozen=out/'frozen_sources';frozen.mkdir()
    for path in (Path(__file__),Path(A.__file__),Path(E.__file__),Path(I.__file__),Path(Q.__file__)):
        shutil.copy2(path,frozen/path.name)
    E.write_json(out/'freeze.json',dict(stage=stage,inventory_sha256=E.digest(out/'inventory.json'),
        source_sha256={str(p):E.digest(p) for p in frozen.glob('*.py')},release_path=str(gate),release_sha256=E.digest(gate),
        estimands_path=str(e0/'inputs/estimands.json'),estimands_sha256=E.digest(e0/'inputs/estimands.json'),
        authority_sha256=authority,native_qualification=native_qualification,E1_forced_rule=forced_rule,forced_rule_status='root_approved_observed_control_definition' if forced_rule else 'not_activated'))
    return dict(prepared=True,stage=stage,counts=inventory['counts'],next='SLURM verify; then worker arrays and independent audits')


def prepare_original_seed_sensitivity(out,parent):
    """A separate long-observation arm for the exact retained old seed cohort."""
    parent=parent.resolve();out=out.resolve()
    if not out.is_relative_to((I.PACKAGE/'data/runs').resolve()):raise ValueError('Prepare read-only analyses under manuscript/data/runs')
    if out.exists():raise FileExistsError('Use a fresh original-seed sensitivity directory')
    frozen=E.read_json(parent/'freeze.json');status=E.read_json(parent/'results/analysis/status.json')
    if frozen['stage']!='E4' or not status['passed'] or not status['complete'] or status['stage']!='E4':raise ValueError('Complete E4 analysis is required before its exact-original-seed sensitivity')
    for name,digest in status['outputs'].items():E.require_hash(parent/'results/analysis'/name,digest)
    source=E.read_json(parent/'inventory.json');E.require_hash(parent/'inventory.json',frozen['inventory_sha256'])
    noises=[r for r in source['cases'] if r['kind']=='noise']
    if any(type(r.get('new_gaussian_path')) is not bool for r in noises):raise ValueError('E4 must explicitly distinguish original Gaussian paths from new top-up seeds')
    retained=[dict(r,observation_arm='E4_original_seed_sensitivity') for r in noises if r['new_gaussian_path'] is False]
    if any(not r.get('retained_prefix') for r in retained):raise ValueError('Every retained old seed needs its exact original-prefix identity metadata')
    keys={(r['model'],r['C_label'],r['f'],r['native_horizon']) for r in retained}
    expected={(r['model'],r['C_label'],r['f'],r['native_horizon']) for r in noises}
    if keys!=expected:raise ValueError('Every E4 condition must retain its exact original-seed cohort')
    controls=[dict(r,observation_arm='E4_original_seed_sensitivity') for r in source['cases'] if r['kind']!='noise']
    rows=retained+controls
    for key in keys:
        if sum((r['model'],r['C_label'],r['f'],r['native_horizon'])==key for r in retained) not in (2,12):raise ValueError('Retained old cohorts must match the original two- or twelve-seed design')
    stage='E4_original_seed_sensitivity';out.mkdir();(out/'results').mkdir();(out/'inputs').mkdir()
    E.write_json(out/'inventory.json',dict(created_utc=E.stamp(),stage=stage,cases=rows,counts=dict(Counter(r['kind'] for r in rows)),
        parent_inventory_path=str(parent/'inventory.json'),parent_inventory_sha256=E.digest(parent/'inventory.json'),
        parent_analysis_status_path=str(parent/'results/analysis/status.json'),parent_analysis_status_sha256=E.digest(parent/'results/analysis/status.json'),
        observation_arm=stage,scientific_replication='Subset of the same12 long runs; not additional independent realizations',no_native_integrations=True))
    closure=out/'frozen_sources';closure.mkdir()
    for path in (Path(__file__),Path(A.__file__),Path(E.__file__),Path(I.__file__),Path(Q.__file__)):shutil.copy2(path,closure/path.name)
    E.write_json(out/'freeze.json',dict(stage=stage,inventory_sha256=E.digest(out/'inventory.json'),source_sha256={str(p):E.digest(p) for p in closure.glob('*.py')},
        release_path=frozen['release_path'],release_sha256=frozen['release_sha256'],estimands_path=frozen['estimands_path'],estimands_sha256=frozen['estimands_sha256'],
        authority_sha256={str(parent/'inventory.json'):E.digest(parent/'inventory.json'),str(parent/'results/analysis/status.json'):E.digest(parent/'results/analysis/status.json')},
        forced_rule_status='Not an additional sigma arm; uses the same frozen E1 forced-state qualification in final integration',observation_arm=stage))
    return dict(prepared=True,stage=stage,counts=dict(Counter(r['kind'] for r in rows)),condition_groups=len(keys),next='SLURM verify; then worker arrays and independent audits')


def verify(out):
    frozen=E.read_json(out/'freeze.json');E.require_hash(out/'inventory.json',frozen['inventory_sha256'])
    for p,d in frozen['source_sha256'].items():E.require_hash(p,d)
    for p,d in frozen.get('authority_sha256',{}).items():E.require_hash(p,d)
    for run,expected in frozen.get('native_qualification',{}).items():
        actual=require_complete_native_gate(run)
        if actual!=expected:raise ValueError('The accepted native engineering qualification changed after intake')
    E.require_hash(frozen['release_path'],frozen['release_sha256']);release_check(frozen['release_path'])
    inventory=E.read_json(out/'inventory.json');checked={}
    for row in inventory['cases']:
        for role in ('scalar','profile','source_status','basal_checkpoint','basal_profile'):
            p=row[role+'_path'];expected=row[role+'_sha256']
            if p not in checked:checked[p]=E.digest(p)
            if checked[p]!=expected:raise ValueError('Native scientific input changed: '+p)
        if row.get('noise_audit_path'):
            E.require_hash(row['noise_audit_path'],row['noise_audit_sha256'])
        with np.load(row['basal_profile_path'],allow_pickle=False) as z:basal=z['primary']
        with np.load(row['profile_path'],allow_pickle=False) as z:primary=z['primary'][0]
        if not np.array_equal(basal,primary):raise ValueError('Every native initial profile must match its exact compact basal')
    E.write_json(out/'results/native_identity_verification.json',dict(passed=True,complete=True,inventory_sha256=E.digest(out/'inventory.json'),
        native_files=len(checked),cases=len(inventory['cases']),all_initial_native_profiles_match_compact_basals=True))
    stage=frozen['stage'];arms=[]
    for sigma in sorted({r['sigma'] for r in inventory['cases'] if r['kind']=='noise'}):
        child=out/'arms'/('sigma_'+str(sigma).replace('.','p'));child.mkdir(parents=True)
        (child/'results').mkdir();(child/'inputs').mkdir();(child/'logs').mkdir();(child/'analysis_workers').mkdir()
        groups,cases=condition_groups(inventory['cases'],stage,sigma)
        E.write_json(child/'inventory.json',dict(cases=cases,counts=dict(Counter(r['kind'] for r in cases)),parent_inventory_sha256=E.digest(out/'inventory.json')))
        E.write_json(child/'results/native_identity_verification.json',dict(passed=True,complete=True,inventory_sha256=E.digest(child/'inventory.json'),parent_identity_sha256=E.digest(out/'results/native_identity_verification.json')))
        cfg=E.read_json(frozen['estimands_path']);cfg.update(source_inventory=str(child/'inventory.json'),stage=stage,sigma=sigma,
            bootstrap_pairing='E1 same twelve exact latent seed blocks across sigma arms' if stage=='E1' else 'Independent conditions; same condition weights across views',required_fhalf_baseline=.0001)
        E.write_json(child/'inputs/estimands.json',cfg)
        m=dict(created_utc=E.stamp(),stage=stage,sigma=sigma,scope=f'{stage} A1-A4/A7 descriptive stage-local analysis; integrated hypothesis decisions remain separate',
            inventory_sha256=E.digest(child/'inventory.json'),native_identity_sha256=E.digest(child/'results/native_identity_verification.json'),
            estimands_sha256=E.digest(child/'inputs/estimands.json'),source_sha256=frozen['source_sha256'],groups=groups,views=E.VIEWS,
            census_counts=dict(Counter(r['kind'] for r in cases)),required_fhalf_baseline=.0001,bootstrap_resamples=2000)
        E.write_json(child/'analysis_manifest.json',m)
        E.write_json(child/'collector_freeze_v3.json',dict(source_sha256=frozen['source_sha256']))
        for role,module,command,memory in [('analysis','signature_figure4_e0_20261020.py','worker',4),('audit','signature_figure4_audit_20261020.py','worker',4),('collect','signature_figure4_e0_20261020.py','collect',8)]:
            lines=['#!/bin/bash','#SBATCH --partition=red','#SBATCH --qos=small','#SBATCH --cpus-per-task=1',f'#SBATCH --mem={memory}G','#SBATCH --time=01:00:00',f'#SBATCH --output={child}/logs/{role}_%A_%a.log']
            if role!='collect':lines.append(f'#SBATCH --array=0-{len(groups)-1}')
            lines+=['set -euo pipefail','export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1',
                f'/home/4473331/.conda/envs/scaling_sandbox/bin/python -B {out}/frozen_sources/{module} {command} --out {child}'+(' --task "$SLURM_ARRAY_TASK_ID"' if role!='collect' else '')]
            (child/f'submit_{role}.slurm').write_text('\n'.join(lines)+'\n')
        arms.append(dict(sigma=sigma,path=str(child),condition_groups=len(groups),analysis_manifest_sha256=E.digest(child/'analysis_manifest.json')))
    E.write_json(out/'results/arms.json',dict(stage=stage,arms=arms,verified=True,no_native_integrations=True))
    return dict(verified=True,stage=stage,arms=arms)


def collect(out):
    freeze=E.read_json(out/'freeze.json');arms=E.read_json(out/'results/arms.json');tables={};receipts=[];metrics={};paired_arrays={}
    for path,digest in freeze['source_sha256'].items():E.require_hash(path,digest)
    if str(Path(__file__).resolve()) not in freeze['source_sha256']:raise ValueError('Use the actual frozen stage collector')
    for arm in arms['arms']:
        child=Path(arm['path']);gate=E.read_json(child/'results/analysis/status.json');audit=E.read_json(child/'results/analysis/independent_estimand_audit.json')
        am=E.read_json(child/'analysis_manifest.json')
        hmap={(g['model'],g['C_label'],g['f']):g['noise'][0]['native_horizon'] for g in am['groups']}
        nmap={(g['model'],g['C_label'],g['f']):len(g['noise']) for g in am['groups']}
        hm={model:{v for (m,C,f),v in hmap.items() if m==model} for model,_,_ in hmap}
        if not(gate['passed'] and gate['complete'] and audit['passed'] and audit['condition_groups']==arm['condition_groups']):raise ValueError('Complete independently audited sigma arms are required')
        for name,digest in gate['outputs'].items():
            E.require_hash(child/'results/analysis'/name,digest)
            rows=list(csv.DictReader((child/'results/analysis'/name).open()))
            annotated=[]
            for row in rows:
                item=dict(row,stage=arms['stage'],experiment=arms['stage'],sigma=arm['sigma'],analysis_sigma=arm['sigma'],observation_arm=f"{arms['stage']}_sigma{arm['sigma']}")
                if row.get('model') in hm and len(hm[row['model']])==1:item['native_horizon']=next(iter(hm[row['model']]))
                if row.get('C_label') and row.get('f'):
                    key=(row['model'],float(row['C_label']),float(row['f']))
                    if key in nmap:item.setdefault('n_seeds',nmap[key])
                item.setdefault('heading_endpoint_R_threshold',.1)
                item.setdefault('relative_contrast_threshold',.1)
                item.setdefault('raw_region_contour',.5)
                item.setdefault('replication_unit','whole independent noisy realization; exact paired latent-seed blocks across E1 sigma')
                annotated.append(item)
            tables.setdefault(name,[]).extend(annotated)
        receipts.append(dict(sigma=arm['sigma'],status_sha256=E.digest(child/'results/analysis/status.json'),audit_sha256=E.digest(child/'results/analysis/independent_estimand_audit.json')))
        for p in (child/'analysis_workers').glob('*/result.json'):
            result=E.read_json(p)
            with np.load(p.parent/'arrays.npz',allow_pickle=False) as z:
                boot={name:z[name] for name in z.files if name.endswith(('__tau_boot_lower','__tau_boot_upper','__Dtheta_boot','__regime_boot'))}
            for view,m in result['metrics'].items():
                key=(m['model'],m['C_label'],m['f'],view)
                values=metrics.setdefault(key,{})
                values[arm['sigma']]=m
                paired_arrays[(key,arm['sigma'])]=boot
                means=[x for x in result['case_metrics'] if x['kind']=='mean' and x['view']==view]
                if len(means)!=1:raise ValueError('One exact matched mean metric per condition/view required')
                if 0. in values and values[0.]!=means[0]:raise ValueError('Matched mean differs across paired sigma arms')
                values[0.]=means[0]
    flags=[]
    if arms['stage']=='E1':
        rule=freeze.get('E1_forced_rule')
        approved=rule and rule.get('root_approved') is True and rule.get('rule')=='same_abnormal_majority_at_sigma1_both_lower_and_mean'
        for (model,C,f,view),values in metrics.items():
            flag=abnormal_persistence(values) if approved else dict(status='unavailable_rule_not_approved',forced_E1=None)
            flags.append(dict(model=model,C_label=C,f=f,view=view,**flag))
        tables['E1_forced_flags.csv']=flags
        comparisons=[]
        for key,values in metrics.items():
            reference=values[1.];ref=paired_arrays[(key,1.)];prefix=key[3]+'__'
            for sigma in (.25,.5,2.):
                current=values[sigma];data=paired_arrays[(key,sigma)]
                lo,hi=A.interval_ratio(data[prefix+'tau_boot_lower'],data[prefix+'tau_boot_upper'],ref[prefix+'tau_boot_lower'],ref[prefix+'tau_boot_upper'])
                ci=A.interval_quantiles(lo,hi)
                dn=data[prefix+'Dtheta_boot'];dd=ref[prefix+'Dtheta_boot']
                mobility=A.ratio(dn,dd);mobility[(dn<=0)|(dd<=0)]=np.nan
                dci=A.finite_quantiles(mobility)
                row=dict(model=key[0],C_label=key[1],f=key[2],view=key[3],sigma=sigma,reference_sigma=1.,
                    replication='paired twelve retained latent Gaussian seed blocks',paired_seed_blocks=12,bootstrap_resamples=2000,
                    tau_ratio=current['tau']/reference['tau'] if current['tau'] is not None and reference['tau'] is not None else None,
                    tau_ratio_ci_lower=ci[0],tau_ratio_ci_upper=ci[1],tau_ratio_ci_upper_unbounded=bool(np.isinf(ci[1])),
                    Dtheta_ratio=current['Dtheta']/reference['Dtheta'] if current['Dtheta'] is not None and reference['Dtheta'] is not None and min(current['Dtheta'],reference['Dtheta'])>0 else None,
                    Dtheta_ratio_ci_lower=float(dci[0]),Dtheta_ratio_ci_upper=float(dci[1]),mobility_valid_bootstrap_fraction=float(np.isfinite(mobility).mean()))
                difference=data[prefix+'regime_boot']-ref[prefix+'regime_boot'];fci=A.finite_quantiles(difference)
                for j,name in enumerate(E.REGIMES):row.update({name+'_difference':current[name+'_fraction']-reference[name+'_fraction'],name+'_difference_ci_lower':float(fci[0,j]),name+'_difference_ci_upper':float(fci[1,j])})
                comparisons.append(row)
        tables['sigma_comparisons.csv']=comparisons
        flag_lookup={(r['model'],r['C_label'],r['f'],r['view']):r for r in flags}
        for row in tables['condition_metrics.csv']:
            flag=flag_lookup[(row['model'],float(row['C_label']),float(row['f']),row['view'])]['forced_E1']
            mean=row['forced']=='True'
            row.update(forced_mean=mean,forced_E1=flag,forced_combined=True if mean or flag is True else None if flag is None else False)
    controls={}
    for row in tables['case_metrics.csv']:
        if row['kind'] not in ('mean','zero'):continue
        key=(row['id'],row['view'])
        if key in controls:continue
        record=dict(row,sigma=0.,analysis_sigma=0.,deterministic_control=True,bootstrap_resamples=0,
            n_seeds=0,n_deterministic_controls=1,replication_unit='one deterministic matched control; no bootstrap')
        for suffix in ('status','','lower','upper','upper_unbounded','coverage'):
            source='per_case_tau'+('_'+suffix if suffix else '')
            target='tau'+('_'+suffix if suffix else '')
            record[target]=record[source]
        controls[key]=record
    tables['control_metrics.csv']=list(controls.values())
    target=out/'results/analysis';target.mkdir(exist_ok=False)
    for name,rows in tables.items():E.csv_write(target/name,rows)
    E.write_json(target/'status.json',dict(generated_utc=E.stamp(),passed=True,complete=True,stage=arms['stage'],
        arms=receipts,outputs={name:E.digest(target/name) for name in tables},scope='Stage-local A1-A4/A7. Integrated H1/H2/H3 and A6 decisions remain separate.',
        E1_additional_forced_rule_activated=bool(flags and flags[0]['status']!='unavailable_rule_not_approved'),no_native_integrations=True))
    return dict(passed=True,complete=True,stage=arms['stage'],output=str(target))


if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True);a.add_argument('--stage',choices=('E1','E2','E3','E4'),required=True);a.add_argument('--e0',type=Path,required=True);a.add_argument('--native-runs',type=Path,nargs='+',required=True);a.add_argument('--gate',type=Path,required=True);a.add_argument('--forced-rule',type=Path)
    a=sub.add_parser('prepare-original-seeds');a.add_argument('--out',type=Path,required=True);a.add_argument('--parent',type=Path,required=True)
    for name in ('verify','collect'):sub.add_parser(name).add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    result=prepare(args.out.resolve(),args.stage,args.e0.resolve(),[r.resolve() for r in args.native_runs],args.gate.resolve(),E.read_json(args.forced_rule) if args.forced_rule else None) if args.command=='prepare' else prepare_original_seed_sensitivity(args.out,args.parent) if args.command=='prepare-original-seeds' else verify(args.out.resolve()) if args.command=='verify' else collect(args.out.resolve())
    print(json.dumps(result,indent=2,allow_nan=False))
