"""Frozen final Figure 4 criteria, finite-window qualification and seed effects.

This module does not integrate models. All decisions are explicit pass/fail/
inconclusive records. Censoring and unavailable observations remain distinct.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import gzip
import io
from datetime import datetime,timezone
import numpy as np
import signature_figure4_analysis_20261020 as A

VIEWS=('main','analysis_first_half','analysis_second_half','late_quarter','cadence_x2','cadence_x4','cadence_x8')
MODES=('multiple_fronts','low_contrast')
MODELS=('wave_pinning','goryachev','holmes_model3','spring_ring')
PRIMARY_A6=('multiple_front_fraction','region_count_mean')
SOURCE_FREEZE='source_freeze_v3.json'
FROZEN_SOURCES='frozen_sources_v3'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4<<20),b''):h.update(block)
    return h.hexdigest()


def number(value):
    if value is None or value=='':return None
    x=float(value)
    return None if math.isnan(x) else x


def boolean(value):
    if value is True or value=='True':return True
    if value is False or value=='False':return False
    return None


def interval(row,stem):
    lower=number(row.get(stem+'_lower'));upper=number(row.get(stem+'_upper'))
    if boolean(row.get(stem+'_upper_unbounded')) is True:upper=math.inf
    return (lower,upper) if lower is not None and upper is not None and lower<=upper else None


def decision(name,status,reason,**values):
    if status not in ('pass','fail','inconclusive'):raise ValueError('Explicit three-way criterion status required')
    return dict(criterion=name,status=status,reason=reason,**values)


def conjunction(name,parts):
    statuses=[p['status'] for p in parts]
    status='fail' if 'fail' in statuses else 'inconclusive' if not parts or 'inconclusive' in statuses else 'pass'
    return decision(name,status,'All required criteria are reported separately',components=parts)


def at_least(name,parts,required):
    passed=sum(p['status']=='pass' for p in parts);unknown=sum(p['status']=='inconclusive' for p in parts)
    status='pass' if passed>=required else 'fail' if passed+unknown<required else 'inconclusive'
    return decision(name,status,'Predeclared model count',passes=passed,inconclusive=unknown,required=required,total=len(parts))


def forced(row):
    if 'forced_combined' in row:return boolean(row['forced_combined'])
    mean=boolean(row.get('forced_mean',row.get('forced')))
    additional=boolean(row.get('forced_E1'))
    required=boolean(row.get('E1_required')) is True
    if mean is True or additional is True:return True
    if mean is None or (required and additional is None):return None
    return False


def persistence_eligible(row):
    coverage=number(row.get('tau_coverage'))
    if coverage is None or coverage<.5 or interval(row,'tau') is None:return False
    if row.get('tau_status') in ('resolved','right_censored'):return True
    if row.get('tau_status')!='unresolved_below_two_cadences' or number(row.get('tau')) is not None:return False
    cadence=number(row.get('cadence'));observed=interval(row,'tau');bracket=row.get('tau_crossing_bracket')
    if isinstance(bracket,str):
        try:bracket=json.loads(bracket)
        except (ValueError,TypeError):return False
    if cadence is None or not math.isfinite(cadence) or cadence<=0 or not isinstance(bracket,(list,tuple)) or len(bracket)!=2:return False
    try:lo,hi=map(float,bracket)
    except (ValueError,TypeError):return False
    return (observed[0]==0 and math.isfinite(observed[1]) and math.isclose(observed[1],2*cadence,rel_tol=1e-12,abs_tol=1e-12)
        and math.isfinite(lo) and math.isfinite(hi) and 0<=lo<hi<=2*cadence+1e-12)


def comparison_eligible(row):
    if boolean(row.get('evaluable')) is not True:return False
    return min(number(row.get('coverage_small')) or 0.,number(row.get('coverage_large')) or 0.)>=.5


def include(row,inclusive):
    return True if inclusive else forced(row) is False


def h1_model(ratios,small_tolerance,large_tolerance,*,inclusive=False):
    selected=[r for r in ratios if comparison_eligible(r) and include(r,inclusive)]
    unknown_forced=any(comparison_eligible(r) and forced(r) is None for r in ratios) and not inclusive
    intervals=[interval(r,'tau_ratio_ci') for r in selected]
    reverse=[r for r,b in zip(selected,intervals) if b and b[1]<1]
    gain=bool(selected) and not unknown_forced and all(b is not None and b[0]>1 for b in intervals)
    gain_status='pass' if gain else 'fail' if reverse else 'inconclusive'
    gains=decision('H1_all_evaluable_gains',gain_status,'Every evaluable ratio CI must exclude1 above; undefined or missing qualification cannot pass',
        evaluable_inputs=len(selected),gain_inputs=sum(b is not None and b[0]>1 for b in intervals),reversal_inputs=[r.get('f') for r in reverse],unknown_forced_qualification=unknown_forced)
    no_reverse=decision('H1_no_supported_reversal','fail' if reverse else 'inconclusive' if unknown_forced or not selected else 'pass',
        'All designed comparisons remain in the coverage table',reversal_inputs=[r.get('f') for r in reverse])
    a=interval(small_tolerance,'fhalf_ci') if small_tolerance else None
    b=interval(large_tolerance,'fhalf_ci') if large_tolerance else None
    tolerance=decision('H1_noise_tolerance','pass' if a and b and b[0]>a[1] else 'fail' if a and b and b[1]<=a[0] else 'inconclusive',
        'Largest/smallest weakest-input tolerance CIs must be non-overlapping in the increasing direction',small_ci=a,large_ci=b)
    return conjunction('H1_model',[gains,no_reverse,tolerance])


def h2_model(ratios,*,inclusive=False):
    ordered=sorted(ratios,key=lambda r:float(r['f']))
    selected=[r for r in ordered if comparison_eligible(r) and include(r,inclusive)]
    reversals=[r for r in selected if interval(r,'tau_ratio_ci') and interval(r,'tau_ratio_ci')[1]<1]
    ranges=[]
    for later in selected:
        ci=interval(later,'tau_ratio_ci')
        if ci is None or not ci[0]<=1<=ci[1]:continue
        for earlier in selected:
            if float(earlier['f'])>=float(later['f']):continue
            span=[r for r in ordered if float(earlier['f'])<=float(r['f'])<=float(later['f'])]
            if not all(comparison_eligible(r) and include(r,inclusive) for r in span):continue
            if not all(number(r.get('large_multiple_front_fraction')) is not None and float(r['large_multiple_front_fraction'])>.5 for r in span):continue
            peak=number(earlier.get('tau_ratio'));end=number(later.get('tau_ratio'))
            observed=peak is not None and end is not None and end>0 and peak>=2*end
            eb=interval(earlier,'tau_ratio');lb=interval(later,'tau_ratio')
            certified=eb is not None and lb is not None and math.isfinite(lb[1]) and lb[1]>0 and eb[0]>=2*lb[1]
            if observed or certified:ranges.append(dict(earlier_f=earlier['f'],later_f=later['f'],fold_drop=peak/end if observed else None,observation_bound_certified=certified))
    return decision('H2_model','pass' if reversals or ranges else 'inconclusive',
        'Explicit existence rule over all designed input pairs; no forced or unavailable gap is silently bridged',
        reversal_inputs=[r['f'] for r in reversals],peak_drop_ranges=ranges,evaluable_inputs=len(selected))


def a6_gate(features,*,inclusive=False):
    records=[]
    for name in PRIMARY_A6:
        members=[r for r in features if r.get('statistic',r.get('feature'))==name and r.get('outcome')=='direct_tau'
                 and boolean(r.get('forced_inclusive')) is inclusive
                 and (inclusive or boolean(r.get('primary_view')) is True)
                 and r.get('imaging_view')=='primary' and boolean(r.get('primary_H2_feature')) is True]
        if len(members)!=1:
            records.append(decision(name,'inconclusive','Exactly one audited primary direct-tau predictor row is required'));continue
        r=members[0];rho=number(r.get('pooled_spearman',r.get('pooled_rho')))
        success=int(r.get('successful_folds',0));unavailable=int(r.get('unavailable_folds',0));planned=int(r.get('planned_folds',4))
        if planned!=4:raise ValueError('A6 denominator is the four predeclared models')
        status='pass' if rho is not None and rho<=-.6 and success>=3 else 'fail' if rho is not None and (rho>-.6 or success+unavailable<3) else 'inconclusive'
        records.append(decision(name,status,'Only the two frozen direct-tau negative-association predictors can support H2',pooled_spearman=rho,successful_folds=success,unavailable_folds=unavailable,planned_folds=4))
    status='pass' if any(r['status']=='pass' for r in records) else 'fail' if all(r['status']=='fail' for r in records) else 'inconclusive'
    return decision('H2_A6_predictor_gate',status,'Missing folds are not successes; incomplete coverage is distinguished from evaluated disagreement',features=records)


def strong_input(strong,weak,*,inclusive=False,small_strong=None,small_weak=None):
    if number(weak.get('f'))!=.0001:return decision('H3_strong_input','inconclusive','The actual weakest noisy baseline is required')
    if number(strong.get('f')) is None or float(strong['f'])<=float(weak['f']):return decision('H3_strong_input','fail','The tested input must exceed the weakest input')
    eligibility=persistence_eligible(strong) and persistence_eligible(weak)
    if small_strong is not None:eligibility=eligibility and persistence_eligible(small_strong)
    if small_weak is not None:eligibility=eligibility and persistence_eligible(small_weak)
    relevant=[strong,weak]+([small_strong] if small_strong is not None else [])+([small_weak] if small_weak is not None else [])
    if not eligibility or not all(include(r,inclusive) for r in relevant):
        return decision('H3_strong_input','inconclusive','Coverage, empirical persistence or forced qualification is unavailable/excluded; contrast loss is retained separately')
    a,b=interval(strong,'tau'),interval(weak,'tau')
    passed=a[1]<.5*b[0]
    contradicted=math.isfinite(b[1]) and a[0]>=.5*b[1]
    ci_a,ci_b=interval(strong,'tau_ci'),interval(weak,'tau_ci')
    return decision('H3_strong_input','pass' if passed else 'fail' if contradicted else 'inconclusive',
        'Conservative empirical observation-bound test; a censored weak baseline is a lower bound',
        strong_bounds=a,weak_bounds=b,confidence_bound_qualifier=bool(ci_a and ci_b and ci_a[1]<.5*ci_b[0]),
        strong_regime=strong.get('regime'),strong_tau_status=strong.get('tau_status'),weak_tau_status=weak.get('tau_status'),
        strong_point_tau=number(strong.get('tau')),weak_point_tau=number(weak.get('tau')))


def e1_mode_corroboration(rows,mode,*,inclusive=False):
    if len({float(r['sigma']) for r in rows})!=len(rows):raise ValueError('Duplicate sigma observations require prior declared assembly')
    identities={(r.get('model'),number(r.get('C_label')),number(r.get('f'))) for r in rows}
    if len(identities)>1:raise ValueError('E1 corroboration must hold at exactly the same model/f/C')
    bysigma={float(r['sigma']):r for r in rows}
    baseline=bysigma.get(1.);reports=[]
    if baseline is None or mode not in MODES or baseline.get('regime')!=mode or not include(baseline,inclusive):
        return decision('H3_E1_corroboration','inconclusive','A qualified sigma-one abnormal majority is required',sigma_checks=[])
    for sigma in (.25,.5,2.):
        r=bysigma.get(sigma)
        reports.append(dict(sigma=sigma,present=r is not None,qualified=r is not None and include(r,inclusive),regime=r.get('regime') if r else None,
            same_abnormal_majority=r is not None and include(r,inclusive) and r.get('regime')==mode))
    matched=any(r['same_abnormal_majority'] for r in reports)
    all_present=all(r['present'] for r in reports)
    complete=all(r['present'] and r['qualified'] for r in reports)
    return decision('H3_E1_corroboration','pass' if matched and all_present else 'fail' if complete else 'inconclusive',
        'Preset existence rule over all three other positive sigmas; every sigma is shown',sigma_checks=reports)


def time_weighted_seed_summary(time,count,contrast,basal_contrast):
    time=np.asarray(time,float);count=np.asarray(count,float);contrast=np.asarray(contrast,float)
    if time.ndim!=1 or len(time)<2 or count.shape!=time.shape or contrast.shape!=time.shape or np.any(np.diff(time)<=0):raise ValueError('Complete ordered frame histories required')
    if basal_contrast<=0 or not all(np.isfinite(x).all() for x in (time,count,contrast)) or np.any(count<0) or np.any(contrast<0):raise ValueError('Finite observed counts and contrasts required')
    dt=np.diff(time);weight=dt/dt.sum()
    return np.array([weight@count[:-1],weight@(contrast[:-1]/basal_contrast)])


def fingerprint_effects(small,large,*,draws=2000,seed=2026102027):
    small=np.asarray(small,float);large=np.asarray(large,float)
    if small.ndim!=2 or large.ndim!=2 or small.shape[1]!=2 or large.shape[1]!=2 or min(len(small),len(large))<1:raise ValueError('One two-component summary per whole seed required')
    if not np.isfinite(small).all() or not np.isfinite(large).all() or np.any(small<0) or np.any(large<0):raise ValueError('Nonnegative finite whole-seed summaries required')
    rng=np.random.default_rng(seed);ws=rng.multinomial(len(small),np.ones(len(small))/len(small),size=draws);wl=rng.multinomial(len(large),np.ones(len(large))/len(large),size=draws)
    a=ws@small/len(small);b=wl@large/len(large);results={}
    for j,name in enumerate(('region_count','normalized_contrast')):
        denominator=small[:,j].mean();numerator=large[:,j].mean();valid=a[:,j]>0
        ratio=np.divide(b[:,j],a[:,j],out=np.zeros(draws),where=valid)
        lower=ratio.copy();upper=ratio.copy();lower[~valid]=0.;upper[~valid]=math.inf
        results[name]=dict(estimate=float(numerator/denominator) if denominator>0 else None,
            ci_lower=float(np.quantile(lower,.025,method='inverted_cdf')),ci_upper=float(np.quantile(upper,.975,method='inverted_cdf')),
            ci_upper_unbounded=bool(np.isinf(np.quantile(upper,.975,method='inverted_cdf'))),
            valid_bootstrap_fraction=float(valid.mean()),independent_small_seeds=len(small),independent_large_seeds=len(large),
            uncertainty_unit='whole seed; independent sizes; joint two-component weights',available=bool(denominator>0))
    return results


def fingerprints_separate(first,second):
    first=first or {};second=second or {}
    checks=[]
    for name in ('region_count','normalized_contrast'):
        a,b=first.get(name),second.get(name)
        ia=interval(a,'ci') if a else None;ib=interval(b,'ci') if b else None
        available=a is not None and b is not None and a.get('available') is True and b.get('available') is True and ia is not None and ib is not None
        disjoint=available and (ia[1]<ib[0] or ib[1]<ia[0])
        checks.append(decision(name,'pass' if disjoint else 'inconclusive','Only dimensionless within-model size effects are compared across models',first_ci=ia,second_ci=ib))
    return conjunction('H3_fingerprint_separation',checks)


def h3_overall(candidates):
    pairs=[]
    for i,a in enumerate(candidates):
        for b in candidates[i+1:]:
            if a['model']==b['model'] or a['mode']==b['mode'] or a['mode'] not in MODES or b['mode'] not in MODES:continue
            fingerprint=fingerprints_separate(a['fingerprints'],b['fingerprints'])
            combined=conjunction('H3_distinct_model_pair',[a['strong'],b['strong'],a['E1'],b['E1'],fingerprint])
            pairs.append(dict(first_model=a['model'],first_f=a['f'],first_mode=a['mode'],second_model=b['model'],second_f=b['f'],second_mode=b['mode'],**combined))
    status='pass' if any(p['status']=='pass' for p in pairs) else 'inconclusive'
    return decision('H3_overall',status,'At least two different models must show distinct corroborated, non-forced failure modes with both seed fingerprints separated',pairs=pairs)


def interval_and_label_compatibility(first,second,*,label):
    if first is None or second is None:
        return decision(label,'inconclusive','A required empirical persistence/coverage view is unavailable')
    if not first.get('regime') or not second.get('regime'):
        return decision(label,'inconclusive','A required A4 regime label is unavailable')
    changed=first.get('regime')!=second.get('regime')
    if not persistence_eligible(first) or not persistence_eligible(second):
        return decision(label,'fail' if changed else 'inconclusive','Observed majority changes remain visible even if the persistence comparison is unavailable',
            regime_changed=changed,first_regime=first.get('regime'),second_regime=second.get('regime'),persistence_comparison_unavailable=True)
    a,b=interval(first,'tau_ci'),interval(second,'tau_ci')
    if a is None or b is None:return decision(label,'fail' if changed else 'inconclusive','Observed majority changes remain visible when required whole-seed persistence intervals are unavailable',
        regime_changed=changed,first_regime=first.get('regime'),second_regime=second.get('regime'),persistence_comparison_unavailable=True)
    disjoint=a[1]<b[0] or b[1]<a[0]
    da,db=interval(first,'Dtheta_ci'),interval(second,'Dtheta_ci')
    mobility_disjoint=bool(da and db and (da[1]<db[0] or db[1]<da[0]))
    return decision(label,'fail' if disjoint or changed else 'pass','Interval and A4-label sensitivity; mobility is reported separately',
        first_tau_ci=a,second_tau_ci=b,persistence_intervals_disjoint=disjoint,regime_changed=changed,
        first_regime=first.get('regime'),second_regime=second.get('regime'),mobility_intervals_disjoint=mobility_disjoint)


def window_qualification(views):
    parts=[interval_and_label_compatibility(views.get(a),views.get(b),label=a+'_versus_'+b) for a,b in
        [('analysis_first_half','analysis_second_half'),('main','late_quarter')]]
    for view in ('cadence_x2','cadence_x4','cadence_x8'):
        parts.append(interval_and_label_compatibility(views.get('main'),views.get(view),label='main_versus_'+view))
    return conjunction('A7_window_qualification',parts)


def support_label(common,longest,window,*,common_only_censored=False,horizon=None):
    incompatible=horizon is not None and horizon['status']=='fail'
    if longest['status']=='pass' and (window['status']=='fail' or incompatible):return 'sensitive_support'
    if common['status']=='pass' and longest['status']=='pass' and window['status']=='pass' and (horizon is None or horizon['status']=='pass'):
        return 'robust_common_and_long'
    if longest['status']=='pass' and common['status']=='inconclusive' and common_only_censored and window['status']=='pass' and (horizon is None or horizon['status']=='pass'):
        return 'qualified_long_horizon_only'
    if longest['status']=='pass' and common['status']=='inconclusive' and window['status']=='pass' and (horizon is None or horizon['status']=='pass'):
        return 'qualified_long_horizon_common_inconclusive'
    if longest['status']=='pass' and common['status']=='fail':return 'long_support_with_contrary_common_evidence'
    if longest['status']=='pass':return 'long_support_with_unresolved_qualification'
    if common['status']=='pass':return 'common_only_support_not_confirmed_long'
    return 'not_established'


def assemble_metrics(e2,e3,e4,triggers):
    """Predeclared observation choice, independent of response values."""
    def key(r):return (r['model'],float(r['C_label']),float(r['f']),r['view'])
    original={model:{float(r['C_label']) for r in e2 if r['model']==model} for model in MODELS}
    common={key(r):dict(r,selected_observation_arm='common_E2_augmented_E0') for r in e2}
    if len(common)!=len(e2):raise ValueError('Duplicate common observations require prior declared deduplication')
    reference_checks=[]
    for row in e3:
        k=key(row)
        if float(row['C_label']) in original.get(row['model'],set()):
            reference_checks.append(dict(key=k,matched_in_common=k in common,source='retained_E3_reference'))
        else:
            if k in common:raise ValueError('A new E3 size is duplicated')
            common[k]=dict(row,selected_observation_arm='common_new_E3_size')
    selected={(r['model'],float(r['C_label']),float(r['f'])) for r in triggers if boolean(r.get('trigger')) is True}
    longmap={key(r):r for r in e4};longest={}
    if len(longmap)!=len(e4):raise ValueError('Duplicate long observations cannot select an outcome')
    if selected-{k[:3] for k in common}:raise ValueError('A selected E4 condition is missing its common observation')
    for k,row in common.items():
        if k[:3] in selected:
            if k not in longmap:raise ValueError('A preselected long observation is missing; common fallback is prohibited')
            if boolean(longmap[k].get('long_arm_qualified')) is not True:raise ValueError('A preselected long observation failed/misses native, control or prefix qualification')
            longest[k]=dict(longmap[k],selected_observation_arm='preselected_E4',common_observation_arm=row['selected_observation_arm'])
        else:longest[k]=dict(row,selected_observation_arm='sole_common_arm_not_selected_for_E4',common_observation_arm=row['selected_observation_arm'])
    unexpected={k[:3] for k in longmap}-selected
    if unexpected:raise ValueError('Unplanned E4 observations cannot replace selected arms')
    return list(common.values()),list(longest.values()),reference_checks


def metric_key(row):
    return (row['model'],float(row['C_label']),float(row['f']),row['view'])


def condition_key(row):
    return (row['model'],float(row['C_label']),float(row['f']))


def native_key(case):
    """Exact observation identity, including its processed basal authority."""
    return tuple(case[name+'_sha256'] for name in ('scalar','profile','basal_profile'))


def checked_json(path,expected=None):
    path=Path(path)
    if expected is not None and sha(path)!=expected:raise ValueError('Input identity changed: '+str(path))
    return json.loads(path.read_text())


def checked_table(directory,status,name):
    path=Path(directory)/name
    if name not in status['outputs'] or sha(path)!=status['outputs'][name]:raise ValueError('Audited table identity changed: '+str(path))
    with path.open(newline='') as stream:return list(csv.DictReader(stream))


def load_analysis_child(root):
    """Read only completed, independently audited native-estimand exports."""
    root=Path(root);directory=root/'results/analysis'
    status=checked_json(directory/'status.json');audit=checked_json(directory/'independent_estimand_audit.json')
    manifest=checked_json(root/'analysis_manifest.json',status['analysis_manifest_sha256'])
    if not all(boolean(r.get('passed')) is True and boolean(r.get('complete')) is True for r in (status,audit)):
        raise ValueError('A complete, independently audited analysis arm is required')
    if audit['analysis_manifest_sha256']!=sha(root/'analysis_manifest.json') or audit['inventory_sha256']!=sha(root/'inventory.json'):
        raise ValueError('Independent estimand audit belongs to different inputs')
    if len(manifest['groups'])!=status['condition_groups'] or audit['condition_groups']!=status['condition_groups']:
        raise ValueError('Missing condition work is an error, not censoring')
    identity=checked_json(root/'results/native_identity_verification.json',manifest['native_identity_sha256'])
    if not identity['passed'] or identity['inventory_sha256']!=sha(root/'inventory.json'):
        raise ValueError('Native identity gate is not accepted')
    for p,digest in manifest['source_sha256'].items():
        if sha(p)!=digest:raise ValueError('Frozen estimand producer changed')
    sigma=float(manifest.get('sigma',1.));arrays={};groups={};receipts=[]
    audited={int(r['task']):r for r in audit['receipts']}
    for group in manifest['groups']:
        task=int(group['task']);worker=root/f'analysis_workers/{task:04d}'
        result=checked_json(worker/'result.json');ar=checked_json(root/f'estimand_audit_workers/{task:04d}.json',audited[task]['receipt_sha256'])
        if not result['passed'] or result['analysis_manifest_sha256']!=sha(root/'analysis_manifest.json') or ar['producer_result_sha256']!=sha(worker/'result.json'):
            raise ValueError('A worker changed after its independent audit')
        if sha(worker/'arrays.npz')!=result['arrays_sha256']:raise ValueError('Whole-seed bootstrap archive changed')
        key=group['model'],float(group['C_label']),float(group['f']),sigma
        if key in arrays:raise ValueError('Duplicate condition observations require declared deduplication')
        with np.load(worker/'arrays.npz',allow_pickle=False) as archive:
            arrays[key]={name:archive[name] for name in archive.files if name.endswith(('__tau_boot_lower','__tau_boot_upper','__Dtheta_boot','__Dx_boot','__Cmean_boot'))}
        groups[key]=group
        receipts.append(dict(path=str(worker/'result.json'),sha256=sha(worker/'result.json'),arrays_path=str(worker/'arrays.npz'),arrays_sha256=result['arrays_sha256']))
    rows=checked_table(directory,status,'condition_metrics.csv')
    if len(rows)!=len(groups)*len(VIEWS):raise ValueError('Every required view must be exported')
    if len({metric_key(r) for r in rows})!=len(rows):raise ValueError('Duplicate metric rows')
    return dict(metrics=rows,arrays=arrays,groups=groups,receipts=receipts,
        gates=[dict(path=str(directory/'status.json'),sha256=sha(directory/'status.json')),dict(path=str(directory/'independent_estimand_audit.json'),sha256=sha(directory/'independent_estimand_audit.json'))])


def load_stage(root,expected_stage,expected_status_sha256):
    root=Path(root);directory=root/'results/analysis'
    status=checked_json(directory/'status.json',expected_status_sha256)
    if not status['passed'] or not status['complete'] or status['stage']!=expected_stage:raise ValueError('The declared completed stage is required')
    frozen=checked_json(root/'freeze.json');arms=checked_json(root/'results/arms.json')
    if sha(root/'inventory.json')!=frozen['inventory_sha256']:raise ValueError('Stage inventory changed')
    for p,digest in frozen['source_sha256'].items():
        if sha(p)!=digest:raise ValueError('Frozen stage source changed')
    for p,digest in frozen.get('authority_sha256',{}).items():
        if sha(p)!=digest:raise ValueError('Frozen stage design changed')
    if not arms['verified'] or arms['stage']!=expected_stage:raise ValueError('Stage intake verification is required')
    datasets=[]
    for arm in arms['arms']:
        child=Path(arm['path'])
        if sha(child/'analysis_manifest.json')!=arm['analysis_manifest_sha256']:raise ValueError('Arm manifest changed')
        datasets.append(load_analysis_child(child))
    arrays={};groups={};receipts=[];gates=[]
    for dataset in datasets:
        if arrays.keys()&dataset['arrays'].keys():raise ValueError('Duplicate sigma/condition arms')
        arrays.update(dataset['arrays']);groups.update(dataset['groups']);receipts.extend(dataset['receipts']);gates.extend(dataset['gates'])
    metrics=checked_table(directory,status,'condition_metrics.csv')
    if len(metrics)!=len(groups)*len(VIEWS):raise ValueError('Stage-wide condition coverage is incomplete')
    flags=checked_table(directory,status,'E1_forced_flags.csv') if expected_stage=='E1' else []
    if expected_stage=='E1' and (len({condition_key(r) for r in metrics})!=24 or len(groups)!=96):
        raise ValueError('Full E1 requires all24 conditions and all four positive sigma arms')
    gates.append(dict(path=str(directory/'status.json'),sha256=sha(directory/'status.json')))
    return dict(metrics=metrics,arrays=arrays,groups=groups,receipts=receipts,gates=gates,flags=flags)


def attach_forced_flags(metrics,flags,required_conditions):
    lookup={metric_key(r):r for r in flags}
    if len(lookup)!=len(flags):raise ValueError('Duplicate E1 forced-state flags')
    required=set(required_conditions);out=[]
    for row in metrics:
        r=dict(row);key=metric_key(r);needed=key[:3] in required
        r['regime_descriptor']='Native raw half-contrast region descriptor with exact native basal contrast qualification; separate from A5 optics'
        r['persistence_bound_eligible']=persistence_eligible(r)
        r['forced_mean']=boolean(r.get('forced_mean',r.get('forced')))
        r['E1_required']=needed;r['forced_E1']=boolean(lookup.get(key,{}).get('forced_E1')) if needed else False
        if needed and key not in lookup:raise ValueError('A required E1 qualification is missing work')
        r.pop('forced_combined',None);r['forced_combined']=forced(r);out.append(r)
    return out


def bootstrap_ratio_row(small,large,sa,la):
    """Compare independent whole-seed ensembles, retaining interval draws."""
    view=small['view'];prefix=view+'__';sb=interval(small,'tau');lb=interval(large,'tau')
    lower,upper=A.interval_ratio(*(lb or (0.,math.inf)),*(sb or (0.,math.inf)))
    blo,bhi=A.interval_ratio(la[prefix+'tau_boot_lower'],la[prefix+'tau_boot_upper'],sa[prefix+'tau_boot_lower'],sa[prefix+'tau_boot_upper'])
    ci=A.interval_quantiles(blo,bhi);sp=number(small.get('tau'));lp=number(large.get('tau'))
    force=any(forced(r) is True for r in (small,large))
    force=None if not force and any(forced(r) is None for r in (small,large)) else force
    sd=number(small.get('Dtheta'));ld=number(large.get('Dtheta'))
    sd_boot=sa[prefix+'Dtheta_boot'];ld_boot=la[prefix+'Dtheta_boot']
    sC=number(small.get('actual_C_mean'));lC=number(large.get('actual_C_mean'))
    if sC is None or lC is None or min(sC,lC)<=0:raise ValueError('Actual circumference is required for geometry and mobility')
    geometry=(lC/sC)**2;geometry_boot=(la[prefix+'Cmean_boot']/sa[prefix+'Cmean_boot'])**2
    mobility=A.ratio(sd_boot,ld_boot);mobility[(sd_boot<=0)|(ld_boot<=0)]=np.nan
    mci=A.finite_quantiles(mobility);excess_ci=A.interval_quantiles(blo/geometry_boot,bhi/geometry_boot)
    dx=A.ratio(ld_boot*la[prefix+'Cmean_boot']**2,sd_boot*sa[prefix+'Cmean_boot']**2)
    dx[(sd_boot<=0)|(ld_boot<=0)]=np.nan;dxci=A.finite_quantiles(dx)
    return dict(model=small['model'],f=float(small['f']),view=view,small_C_label=float(small['C_label']),large_C_label=float(large['C_label']),
        tau_ratio=lp/sp if lp is not None and sp is not None and sp>0 else None,
        tau_ratio_lower=float(lower),tau_ratio_upper=float(upper),tau_ratio_upper_unbounded=bool(math.isinf(float(upper))),
        tau_ratio_ci_lower=float(ci[0]),tau_ratio_ci_upper=float(ci[1]),tau_ratio_ci_upper_unbounded=bool(math.isinf(float(ci[1]))),
        evaluable=persistence_eligible(small) and persistence_eligible(large),coverage_small=number(small.get('tau_coverage')),coverage_large=number(large.get('tau_coverage')),
        tau_status_small=small['tau_status'],tau_status_large=large['tau_status'],forced_combined=force,
        finite_resolution_bound_small=small['tau_status']=='unresolved_below_two_cadences' and persistence_eligible(small),
        finite_resolution_bound_large=large['tau_status']=='unresolved_below_two_cadences' and persistence_eligible(large),
        large_multiple_front_fraction=number(large['multiple_fronts_fraction']),
        Dtheta_mobility_proxy=sd/ld if sd is not None and ld is not None and min(sd,ld)>0 else None,
        Dtheta_mobility_proxy_ci_lower=float(mci[0]),Dtheta_mobility_proxy_ci_upper=float(mci[1]),mobility_proxy_valid_bootstrap_fraction=float(np.isfinite(mobility).mean()),
        reference_geometry_ratio=9.,actual_geometry_ratio=geometry,
        excess_persistence=lp/sp/geometry if lp is not None and sp is not None and sp>0 else None,
        excess_persistence_ci_lower=float(excess_ci[0]),excess_persistence_ci_upper=float(excess_ci[1]),excess_persistence_ci_upper_unbounded=bool(math.isinf(excess_ci[1])),
        Dx_large_over_small=ld/sd*geometry if sd is not None and ld is not None and min(sd,ld)>0 else None,Dx_ratio_ci_lower=float(dxci[0]),Dx_ratio_ci_upper=float(dxci[1]),
        small_observation_arm=small.get('selected_observation_arm'),large_observation_arm=large.get('selected_observation_arm'),
        native_horizon_small=number(small.get('native_horizon')),native_horizon_large=number(large.get('native_horizon')),
        n_seeds_small=int(small['n_seeds']),n_seeds_large=int(large['n_seeds']),replication_unit='independent whole seed ensembles, independent resampling across sizes')


def noise_tolerance(metrics,arrays,model,C,view,*,inclusive=False):
    rows=sorted([r for r in metrics if (r['model'],float(r['C_label']),r['view'])==(model,C,view)],key=lambda r:float(r['f']))
    if not rows:raise ValueError('A declared model/size/view has no conditions')
    selected=[r for r in rows if include(r,inclusive)]
    missing_qualification=not inclusive and any(forced(r) is None for r in rows)
    available=len(selected)>=2 and float(rows[0]['f'])==.0001 and float(selected[0]['f'])==.0001 and not missing_qualification
    if not available:
        reason='unavailable_required_forced_qualification' if missing_qualification else 'unavailable_missing_declared_weak_input' if float(rows[0]['f'])!=.0001 else 'unavailable_after_forced_exclusion'
        return dict(model=model,C_label=C,view=view,forced_inclusive=inclusive,fhalf_status=reason,fhalf=None,fhalf_lower=.0001,fhalf_upper=math.inf,
            fhalf_ci_lower=.0001,fhalf_ci_upper=math.inf,fhalf_ci_upper_unbounded=True,weakest_f=.0001,unavailable_is_not_censoring=True)
    f=np.array([float(r['f']) for r in selected]);bounds=[interval(r,'tau') if persistence_eligible(r) else (0.,math.inf) for r in selected]
    point=A.conservative_fhalf(f,[b[0] for b in bounds],[b[1] for b in bounds])
    los=np.array([arrays[condition_key(r)][view+'__tau_boot_lower'] if persistence_eligible(r) else np.zeros(2000) for r in selected])
    his=np.array([arrays[condition_key(r)][view+'__tau_boot_upper'] if persistence_eligible(r) else np.full(2000,math.inf) for r in selected])
    if los.shape!=his.shape or los.shape[1]!=2000:raise ValueError('Exactly2000 whole-seed interval bootstrap draws required')
    draws=[A.conservative_fhalf(f,los[:,b],his[:,b]) for b in range(2000)];ci=A.interval_quantiles([r['lower'] for r in draws],[r['upper'] for r in draws])
    return dict(model=model,C_label=C,view=view,forced_inclusive=inclusive,fhalf_status=point['status'],fhalf=point['estimate'],
        fhalf_lower=point['lower'],fhalf_upper=point['upper'],fhalf_upper_unbounded=math.isinf(point['upper']),fhalf_ci_lower=ci[0],fhalf_ci_upper=ci[1],
        fhalf_ci_upper_unbounded=math.isinf(ci[1]),weakest_f=.0001,weakest_tau_status=selected[0]['tau_status'],
        excluded_forced_strengths=[float(r['f']) for r in rows if not include(r,inclusive)],bootstrap_resamples=2000)


def physical_seed_summary(frames,basal_contrast,start,stop):
    """Restrict observed constant-frame intervals without interpolating fields."""
    time=np.array([r['time'] for r in frames],float)
    count=np.array([r['region_count'] for r in frames],float);contrast=np.array([r['contrast'] for r in frames],float)
    if len(time)<2 or np.any(np.diff(time)<=0) or start<time[0]-1e-8 or stop>time[-1]+1e-8 or start>=stop:
        raise ValueError('A fingerprint requires complete observed coverage of its physical window')
    if basal_contrast<=0 or not all(np.isfinite(v).all() for v in (time,count,contrast)) or np.any(count<0) or np.any(contrast<0):raise ValueError('Invalid observed fingerprint frame')
    overlap=np.maximum(0.,np.minimum(time[1:],stop)-np.maximum(time[:-1],start))
    if not np.isclose(overlap.sum(),stop-start,rtol=1e-10,atol=1e-9):raise ValueError('The physical fingerprint window is incomplete')
    return np.array([overlap@count[:-1]/overlap.sum(),overlap@(contrast[:-1]/basal_contrast)/overlap.sum()])


def validate_A5_audits(audits,inventory_sha256,canonical_ids,accepted_receipts):
    """Bind passed audit masters or their hashed shards to every exact case."""
    if not audits:raise ValueError('Independent native A5 frame audits are required')
    expected=set(canonical_ids)
    if len(expected)!=len(canonical_ids) or set(accepted_receipts)!=expected:raise ValueError('Canonical A5 case coverage is duplicated or differs from accepted measurement receipts')
    verified=[]
    for spec in audits:
        master=checked_json(spec['path'],spec['sha256'])
        if master.get('passed') is not True:raise ValueError('A5 independent audit failed')
        if master.get('inventory_sha256') is not None and master['inventory_sha256']!=inventory_sha256:raise ValueError('A5 independent audit belongs to a different canonical inventory')
        count=master.get('cases',master.get('native_cases'))
        if count is not None and int(count)!=len(expected):raise ValueError('A5 independent audit canonical case count is incomplete')
        shards=master.get('source_shards',master.get('shards',[]))
        if isinstance(shards,list) and shards:
            payloads=[checked_json(row['path'],row['sha256']) for row in shards]
        else:payloads=[master]
        records={}
        for payload in payloads:
            if payload.get('passed') is not True or payload.get('inventory_sha256')!=inventory_sha256:raise ValueError('A5 independent audit evidence is not bound to the selected inventory')
            raw=payload.get('records')
            if isinstance(raw,dict):raw=[dict(row,case_id=key) for key,row in raw.items()]
            if not isinstance(raw,list):raise ValueError('A5 independent audit lacks exact canonical case records')
            if payload.get('cases') is not None and int(payload['cases'])!=len(raw):raise ValueError('A5 audit shard case count disagrees with its records')
            for row in raw:
                case_id=row.get('case_id',row.get('id'))
                if case_id in records or case_id not in expected or row.get('passed') is not True:raise ValueError('A5 independent audit has duplicate, unrelated or failed canonical case records')
                measurement=row.get('measurement_sha256',row.get('measured_result_sha256'))
                if measurement is not None and measurement!=accepted_receipts[case_id]['sha256']:raise ValueError('A5 measured result changed after its independent audit')
                records[case_id]=row
        if set(records)!=expected:raise ValueError('A5 independent audit does not cover every canonical noisy case and control')
        if isinstance(shards,list) and shards and master.get('records') is not None:
            direct=master['records']
            if isinstance(direct,dict):direct=[dict(row,case_id=key) for key,row in direct.items()]
            if not isinstance(direct,list) or len(direct)!=len(expected):raise ValueError('A5 master records disagree with their verified shard coverage')
            master_ids=[]
            for row in direct:
                case_id=row.get('case_id',row.get('id'));master_ids.append(case_id)
                measurement=row.get('measurement_sha256',row.get('measured_result_sha256'))
                if case_id not in expected or row.get('passed') is not True or (measurement is not None and measurement!=accepted_receipts[case_id]['sha256']):raise ValueError('A5 master records disagree with their independently hashed measurement evidence')
                if row.get('primary_sha256') is not None and records[case_id].get('primary_sha256') not in (None,row['primary_sha256']):raise ValueError('A5 master primary hashes disagree with their verified shards')
            if set(master_ids)!=expected or len(set(master_ids))!=len(master_ids):raise ValueError('A5 master has duplicate or missing canonical case records')
        verified.append(records)
    return verified


def load_fingerprint_cases(inventory,worker_root,*,expected_inventory_sha256,complete,audits):
    """Only accepted A5 observations are eligible for final H3 fingerprints."""
    inventory=Path(inventory);worker_root=Path(worker_root)
    inv=checked_json(inventory,expected_inventory_sha256)
    completion=checked_json(complete['path'],complete['sha256'])
    if not completion.get('passed',completion.get('complete',False)) or completion.get('complete') is False:raise ValueError('Accepted A5 completion is required')
    if completion.get('inventory_sha256')!=sha(inventory) or completion.get('cases')!=len(inv['cases']):raise ValueError('A5 completion belongs to a different canonical inventory')
    receipt_spec=completion.get('seed_receipts',{})
    if not receipt_spec.get('path') or not receipt_spec.get('sha256'):raise ValueError('A5 completed native receipt index is required')
    accepted_receipts=checked_json(receipt_spec['path'],receipt_spec['sha256'])
    accepted={r['id']:r for r in accepted_receipts}
    if len(accepted)!=len(inv['cases']):raise ValueError('A5 native receipt index is incomplete or duplicated')
    audit_records=validate_A5_audits(audits,sha(inventory),[case['id'] for case in inv['cases']],accepted)
    rows={};receipts=[]
    for case in inv['cases']:
        if case['kind']!='noise':continue
        path=worker_root/case['id']/'result.json'
        if case['id'] not in accepted or Path(accepted[case['id']]['path']).resolve()!=path.resolve():raise ValueError('A5 accepted case path differs from selected measurement')
        result=checked_json(path,accepted[case['id']]['sha256'])
        if not result['passed'] or result['inventory_sha256']!=sha(inventory) or native_key(case)!=native_key(result['case']):raise ValueError('A5 measurement does not match the exact canonical case')
        primary=result['primary'];observed=Path(primary['path'])
        if sha(observed)!=primary['sha256']:raise ValueError('Observed A5 frames changed')
        for audited in audit_records:
            if audited[case['id']].get('primary_sha256') not in (None,primary['sha256']):raise ValueError('A5 primary observation changed after its independent audit')
        if not primary['smoothed'] or float(primary['threshold'])!=.5:raise ValueError('H3 uses the frozen primary optical half-contour view')
        cadence=8. if case['model']=='spring_ring' else 80.
        if float(result['common_imaging_cadence'])!=cadence:raise ValueError('Fingerprint cadence is not the common physical A5 cadence')
        with gzip.open(observed,'rt') as stream:variant=json.load(stream)
        record=dict(case=case,frames=variant['result']['frames'],basal_contrast=float(primary['observed_basal_contrast']),
            result_path=str(path),result_sha256=sha(path),primary_path=str(observed),primary_sha256=sha(observed),imaging_cadence=cadence)
        key=native_key(case)
        if key in rows:raise ValueError('Independent A5 observations are duplicated before selection')
        rows[key]=record;receipts.append({k:record[k] for k in ('result_path','result_sha256','primary_path','primary_sha256')})
    return rows,receipts


def fingerprint_summaries(groups,observations,metrics):
    fractions={'main':(.5,1.),'analysis_first_half':(.5,.75),'analysis_second_half':(.75,1.),'late_quarter':(.875,1.),'cadence_x2':(.5,1.),'cadence_x4':(.5,1.),'cadence_x8':(.5,1.)}
    summaries={};rows=[]
    for metric in metrics:
        key=condition_key(metric);group=groups[key];view=metric['view'];values=[];seeds=set()
        for case in group['noise']:
            if case['seed'] in seeds:raise ValueError('Repeated seed cannot be an independent fingerprint replicate')
            seeds.add(case['seed']);record=observations.get(native_key(case))
            if record is None:raise ValueError('The selected native observation lacks its accepted A5 fingerprint')
            horizon=float(case['native_horizon']);a,b=fractions[view]
            value=physical_seed_summary(record['frames'],record['basal_contrast'],a*horizon,b*horizon);values.append(value)
            rows.append(dict(model=key[0],C_label=key[1],f=key[2],view=view,id=case['id'],seed=case['seed'],native_horizon=horizon,
                observed_region_count_mean=float(value[0]),basal_normalized_contrast_mean=float(value[1]),
                optical_basal_contrast=record['basal_contrast'],imaging_cadence=record['imaging_cadence'],
                fingerprint_descriptor='A5 primary optical half-contour count and optical contrast / identically processed basal contrast',
                heading_cadence_view_does_not_change_A5_primary_imaging=view.startswith('cadence_'),
                physical_window_start=a*horizon,physical_window_stop=b*horizon,uncertainty_unit='one whole independent seed',
                A5_result_sha256=record['result_sha256'],A5_primary_sha256=record['primary_sha256']))
        summaries[metric_key(metric)]=np.array(values)
    return summaries,rows


def comparison_tables(metrics,arrays,*,arm):
    lookup={metric_key(r):r for r in metrics};ratios=[];tolerances=[]
    sizes={model:sorted({float(r['C_label']) for r in metrics if r['model']==model}) for model in MODELS}
    for view in VIEWS:
        for model in MODELS:
            if len(sizes[model])<2:raise ValueError('At least two declared sizes required for every model')
            sm,lg=sizes[model][0],sizes[model][-1]
            small_f={float(r['f']) for r in metrics if (r['model'],float(r['C_label']),r['view'])==(model,sm,view)}
            large_f={float(r['f']) for r in metrics if (r['model'],float(r['C_label']),r['view'])==(model,lg,view)}
            for f in sorted(small_f|large_f):
                if f not in small_f or f not in large_f:
                    present=lookup[(model,sm if f in small_f else lg,f,view)]
                    ratios.append(dict(model=model,f=f,view=view,small_C_label=sm,large_C_label=lg,arm=arm,analysis_arm=arm,sigma=1.,evaluable=False,
                        native_horizon_small=number(present.get('native_horizon')) if f in small_f else None,
                        native_horizon_large=number(present.get('native_horizon')) if f in large_f else None,
                        comparison_status='unavailable_unmatched_designed_input',missing_size=sm if f not in small_f else lg,forced_combined=None));continue
                small=lookup[(model,sm,f,view)];large=lookup[(model,lg,f,view)]
                ratios.append(dict(bootstrap_ratio_row(small,large,arrays[(model,sm,f)],arrays[(model,lg,f)]),arm=arm,analysis_arm=arm,sigma=1.))
            for C in sizes[model]:
                for inclusive in (False,True):tolerances.append(dict(noise_tolerance(metrics,arrays,model,C,view,inclusive=inclusive),arm=arm))
    return ratios,tolerances


def evaluate_arm(metrics,arrays,e1_metrics,seed_summaries,a6_features,*,arm):
    """All models, windows, inputs and forced views; no outcome-based selection."""
    if set(r['model'] for r in metrics)!=set(MODELS):raise ValueError('Final decisions require all four predeclared models')
    lookup={metric_key(r):r for r in metrics}
    if len(lookup)!=len(metrics):raise ValueError('Duplicate observations cannot select a favorable outcome')
    conditions={condition_key(r) for r in metrics}
    if any((m,C,f,v) not in lookup for m,C,f in conditions for v in VIEWS):raise ValueError('A required window is missing work')
    sizes={model:sorted({float(r['C_label']) for r in metrics if r['model']==model}) for model in MODELS}
    ratios,tolerances=comparison_tables(metrics,arrays,arm=arm)
    criteria=[];candidates=[];fingerprint_rows=[];all_decisions=[]
    for view in VIEWS:
        for model in MODELS:
            sm,lg=sizes[model][0],sizes[model][-1]
            if sm==lg:raise ValueError('At least two declared sizes required')
            fs=sorted({float(r['f']) for r in metrics if (r['model'],float(r['C_label']),r['view'])==(model,sm,view)} &
                      {float(r['f']) for r in metrics if (r['model'],float(r['C_label']),r['view'])==(model,lg,view)})
            for f in fs:
                a=seed_summaries.get((model,sm,f,view));b=seed_summaries.get((model,lg,f,view))
                effects=fingerprint_effects(a,b) if a is not None and b is not None else {}
                for statistic,values in effects.items():fingerprint_rows.append(dict(arm=arm,model=model,f=f,view=view,statistic=statistic,**values))
        for inclusive in (False,True):
            per_h1=[];per_h2=[];per_candidates=[]
            for model in MODELS:
                sm,lg=sizes[model][0],sizes[model][-1]
                model_ratios=[r for r in ratios if r['model']==model and r['view']==view]
                def tolerance(C):return next(r for r in tolerances if (r['model'],r['C_label'],r['view'],r['forced_inclusive'])==(model,C,view,inclusive))
                h1=h1_model(model_ratios,tolerance(sm),tolerance(lg),inclusive=inclusive);h2=h2_model(model_ratios,inclusive=inclusive)
                per_h1.append(h1);per_h2.append(h2)
                for result in (h1,h2):
                    item=dict(arm=arm,view=view,forced_inclusive=inclusive,model=model,**result);criteria.append(item);all_decisions.append(item)
                weak=lookup.get((model,lg,.0001,view))
                strong_rows=sorted([r for r in metrics if (r['model'],float(r['C_label']),r['view'])==(model,lg,view) and float(r['f'])>.0001],key=lambda r:float(r['f']))
                for row in strong_rows:
                    f=float(row['f']);mode=row['regime'];small=lookup.get((model,sm,f,view))
                    strong=strong_input(row,weak,inclusive=inclusive) if weak else decision('H3_strong_input','inconclusive','Declared weakest input missing at largest size')
                    matches=[r for r in e1_metrics if (r['model'],float(r['C_label']),float(r['f']),r['view'])==(model,lg,f,view)]
                    corroboration=e1_mode_corroboration(matches,mode,inclusive=inclusive)
                    a=seed_summaries.get((model,sm,f,view));b=seed_summaries.get((model,lg,f,view))
                    effects=fingerprint_effects(a,b) if a is not None and b is not None and small is not None and include(small,inclusive) and include(row,inclusive) else {}
                    candidate=dict(model=model,f=f,mode=mode,mode_descriptor='A4/E1 native raw descriptor; fingerprints separately use A5 optical observations',strong=strong,E1=corroboration,fingerprints=effects)
                    per_candidates.append(candidate)
                    candidates.append(dict(arm=arm,view=view,forced_inclusive=inclusive,**candidate))
            A6=a6_gate(a6_features,inclusive=inclusive)
            results=[at_least('H1_overall',per_h1,3),conjunction('H2_overall',[at_least('H2_model_count',per_h2,2),A6]),h3_overall(per_candidates)]
            for result in results:
                item=dict(arm=arm,view=view,forced_inclusive=inclusive,model='all',A6_scope='primary imaging and direct tau, main view; not reestimated as an additional window gate',**result)
                criteria.append(item);all_decisions.append(item)
    windows=[]
    for key in sorted(conditions):
        result=window_qualification({v:lookup[key+(v,)] for v in VIEWS})
        windows.append(dict(arm=arm,model=key[0],C_label=key[1],f=key[2],**result))
    return dict(criteria=criteria,strong_input_candidates=candidates,size_ratios=ratios,noise_tolerance=tolerances,
        fingerprint_effects=fingerprint_rows,window_qualification=windows,decisions=all_decisions)


def qualify_long_rows(rows,gate,triggers):
    selected={condition_key(r) for r in triggers if boolean(r.get('trigger')) is True}
    qualified={condition_key(r) for r in gate.get('qualified_conditions',[])}
    required=('passed','complete','all_preselected_arms_complete','all_required_prefix_checks_passed','all_required_controls_passed')
    if not all(gate.get(name) is True for name in required) or qualified!=selected:
        raise ValueError('All preselected E4 conditions require accepted native, control and original-prefix gates')
    if {condition_key(r) for r in rows}!=selected:raise ValueError('The long observations must exactly cover the preselected condition grid')
    return [dict(r,long_arm_qualified=True) for r in rows]


def select_sigma_one(dataset):
    metrics=[r for r in dataset['metrics'] if float(r.get('sigma',1.))==1.]
    arrays={k[:3]:v for k,v in dataset['arrays'].items() if k[3]==1.}
    groups={k[:3]:v for k,v in dataset['groups'].items() if k[3]==1.}
    if {condition_key(r) for r in metrics}!=set(groups):raise ValueError('The sigma-one condition map is incomplete')
    return dict(metrics=metrics,arrays=arrays,groups=groups)


def select_observation_inputs(metrics,e2,e3,e4,*,longest):
    arrays={};groups={}
    for row in metrics:
        key=condition_key(row)
        origin=row.get('common_observation_arm',row['selected_observation_arm'])
        source=e4 if row['selected_observation_arm']=='preselected_E4' else e3 if origin=='common_new_E3_size' else e2
        if not longest and row['selected_observation_arm']=='preselected_E4':raise ValueError('Common horizon cannot silently use a long observation')
        if key not in source['arrays'] or key not in source['groups']:raise ValueError('The selected estimator lacks its exact whole-seed ensemble')
        arrays[key]=source['arrays'][key];groups[key]=source['groups'][key]
    return arrays,groups


def horizon_qualifications(common,longest):
    lookup={metric_key(r):r for r in longest};out=[]
    for first in common:
        key=metric_key(first);last=lookup.get(key)
        if last is None:raise ValueError('A declared common observation has no selected longest arm')
        result=interval_and_label_compatibility(first,last,label='common_versus_longest')
        out.append(dict(model=key[0],C_label=key[1],f=key[2],view=key[3],
            common_native_horizon=number(first.get('native_horizon')),long_native_horizon=number(last.get('native_horizon')),
            selected_long_observation_arm=last.get('selected_observation_arm'),**result))
    return out


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not math.isfinite(value):return None
    return value


def write_json(path,value):
    Path(path).write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n')


def csv_bytes(rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
    for row in rows:writer.writerow({k:json.dumps(clean(v),separators=(',',':'),allow_nan=False) if isinstance(v,(dict,list,tuple)) else v for k,v in row.items()})
    return stream.getvalue().encode('utf-8')


def write_csv(path,rows):
    Path(path).write_bytes(csv_bytes(rows))


def load_a6(spec,observation_status):
    root=Path(spec['root']);complete=checked_json(root/'complete.json',spec['complete_sha256'])
    audit=checked_json(root/'independent_A6_audit.json',spec['independent_audit_sha256'])
    if not complete.get('passed',complete.get('complete',False)) or audit.get('passed') is not True:
        raise ValueError('Independent completed A6 gate required')
    source=complete.get('sources',{});observation_status=Path(observation_status)
    if source.get('ratio_status_sha256')!=sha(observation_status) or source.get('size_ratios_sha256')!=sha(observation_status.parent/'size_ratios.csv'):
        raise ValueError('A6 must use this exact accepted common or longest persistence arm; E0 cannot be relabelled as final')
    path=root/'H2_primary_features.csv'
    expected=complete.get('outputs',{}).get(path.name)
    if isinstance(expected,dict):expected=expected.get('sha256')
    if expected is None:expected=spec.get('features_sha256')
    if not expected or sha(path)!=expected:raise ValueError('Audited A6 predictor table identity is required')
    with path.open(newline='') as stream:features=list(csv.DictReader(stream))
    for inclusive in (False,True):
        selected=[r for r in features if r.get('statistic') in PRIMARY_A6 and r.get('imaging_view')=='primary' and r.get('outcome')=='direct_tau' and boolean(r.get('forced_inclusive')) is inclusive and boolean(r.get('primary_H2_feature')) is True]
        if {r['statistic'] for r in selected}!=set(PRIMARY_A6) or len(selected)!=2:
            raise ValueError('Both frozen A6 predictors and forced views must remain visible')
    return features


def assemble_input_data(spec):
    if set(spec.get('stages',{}))!={'E1','E2','E3','E4'}:raise ValueError('Full final decisions require all four future stages')
    datasets={name:load_stage(row['root'],name,row['status_sha256']) for name,row in spec['stages'].items()}
    stages={name:select_sigma_one(dataset) for name,dataset in datasets.items() if name!='E1'}
    trigger_spec=spec['E4_triggers'];triggerpath=Path(trigger_spec['path'])
    if sha(triggerpath)!=trigger_spec['sha256']:raise ValueError('The prespecified E4 trigger grid changed')
    with triggerpath.open(newline='') as stream:triggers=list(csv.DictReader(stream))
    if sum(boolean(r.get('trigger')) is True for r in triggers)!=52:raise ValueError('The frozen E4 design contains exactly52 selected conditions')
    gate_spec=spec['E4_native_gate'];long_gate=checked_json(gate_spec['path'],gate_spec['sha256'])
    stages['E4']['metrics']=qualify_long_rows(stages['E4']['metrics'],long_gate,triggers)
    common,longest,references=assemble_metrics(stages['E2']['metrics'],stages['E3']['metrics'],stages['E4']['metrics'],triggers)
    flags=datasets['E1']['flags'];required={condition_key(r) for r in flags}
    if len(required)!=24:raise ValueError('The full reviewed E1 forced-qualification grid must be present')
    common=attach_forced_flags(common,flags,required);longest=attach_forced_flags(longest,flags,required)
    e1=attach_forced_flags(datasets['E1']['metrics'],flags,required)
    return dict(datasets=datasets,stages=stages,common=common,longest=longest,references=references,e1=e1)


def prepare_observations(bundle,spec_path):
    """Export final arm-specific inputs before A6 fits, avoiding circular joins."""
    bundle=Path(bundle);verify_source_freeze(bundle);spec=checked_json(spec_path);data=assemble_input_data(spec)
    parent=bundle/'results/assembled_observations';parent.mkdir(exist_ok=False);receipts=[]
    for arm in ('common','longest'):
        metrics=data[arm];stages=data['stages']
        arrays,groups=select_observation_inputs(metrics,stages['E2'],stages['E3'],stages['E4'],longest=arm=='longest')
        ratios,tolerances=comparison_tables(metrics,arrays,arm=arm);out=parent/arm;out.mkdir()
        for name,rows in dict(condition_metrics=[dict(r,analysis_arm=arm) for r in metrics],size_ratios=ratios,noise_tolerance=tolerances).items():write_csv(out/(name+'.csv'),rows)
        write_json(out/'status.json',dict(passed=True,complete=True,analysis_arm=arm,sigma=1.,forcing_sigma=1.,condition_groups=len(groups),
            source_freeze_sha256=sha(bundle/SOURCE_FREEZE),rule_freeze_sha256=sha(bundle/'rule_freeze.json'),
            input_spec_sha256=sha(spec_path),source_stage_gates={name:dataset['gates'] for name,dataset in data['datasets'].items()},
            outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()},
            scope='Frozen observation choice and whole-seed persistence comparisons; A6 consumes this exact status and ratio table. Final hypothesis decisions remain pending accepted A5/A6 inputs.',no_native_integrations=True))
        receipts.append(dict(analysis_arm=arm,path=str(out),status_sha256=sha(out/'status.json')))
    write_json(parent/'complete.json',dict(passed=True,complete=True,arms=receipts,E3_reference_checks=data['references'],
        scope='Independent estimand audits verified upstream; final observation-arm selection is frozen and response-independent.'))
    return dict(passed=True,complete=True,arms=receipts)


def integrated_evaluation(bundle,spec_path):
    """Final release requires complete stage exports; no native work runs here."""
    bundle=Path(bundle);verify_source_freeze(bundle);spec=checked_json(spec_path);data=assemble_input_data(spec)
    datasets=data['datasets'];stages=data['stages'];common=data['common'];longest=data['longest'];references=data['references'];e1=data['e1']
    original_spec=spec['E4_original_seed_sensitivity']
    original_dataset=load_stage(original_spec['root'],'E4_original_seed_sensitivity',original_spec['status_sha256'])
    original=select_sigma_one(original_dataset)
    if set(original['groups'])!=set(stages['E4']['groups']):raise ValueError('Original-seed sensitivity must retain every preselected E4 condition')
    for key,group in original['groups'].items():
        if any(case.get('new_gaussian_path') is not False or not case.get('retained_prefix') for case in group['noise']):raise ValueError('An original-seed sensitivity contains a new Gaussian path or lacks its original identity')
        if not {native_key(r) for r in group['noise']}<={native_key(r) for r in stages['E4']['groups'][key]['noise']}:raise ValueError('Original-seed sensitivity is not an exact subset of the same long observations')
    original_metrics=attach_forced_flags(original['metrics'],datasets['E1']['flags'],{condition_key(r) for r in datasets['E1']['flags']})
    observations={};A5_receipts=[]
    for row in spec['A5_observation_sets']:
        loaded,receipts=load_fingerprint_cases(**row)
        for key,value in loaded.items():
            if key in observations and observations[key]['primary_sha256']!=value['primary_sha256']:
                raise ValueError('Competing A5 measurements require explicit prior observation selection')
            observations[key]=value
        A5_receipts.extend(receipts)
    tables={};results={};inputs={};seed_rows=[]
    for arm,metrics in (('common',common),('longest',longest)):
        arrays,groups=select_observation_inputs(metrics,stages['E2'],stages['E3'],stages['E4'],longest=arm=='longest')
        summaries,rows=fingerprint_summaries(groups,observations,metrics);seed_rows.extend(dict(r,arm=arm) for r in rows)
        observation_directory=bundle/'results/assembled_observations'/arm
        observation_status=checked_json(observation_directory/'status.json')
        if not observation_status['passed'] or not observation_status['complete']:raise ValueError('Prepare complete frozen observation arms before final A6 evaluation')
        prepared_metrics=checked_table(observation_directory,observation_status,'condition_metrics.csv')
        if {metric_key(r) for r in prepared_metrics}!={metric_key(r) for r in metrics}:raise ValueError('Final observations differ from the accepted A6 input arm')
        expected=hashlib.sha256(csv_bytes([dict(r,analysis_arm=arm) for r in metrics])).hexdigest()
        if observation_status['outputs']['condition_metrics.csv']!=expected:raise ValueError('Final condition values changed after frozen A6 inputs were prepared')
        features=load_a6(spec['A6'][arm],observation_directory/'status.json');result=evaluate_arm(metrics,arrays,e1,summaries,features,arm=arm);results[arm]=result
        for table_name in ('size_ratios','noise_tolerance'):
            expected=hashlib.sha256(csv_bytes(result[table_name])).hexdigest()
            if observation_status['outputs'][table_name+'.csv']!=expected:raise ValueError('Final persistence comparisons changed after accepted A6 input preparation')
        for name,rows in result.items():
            if name!='decisions':tables.setdefault(name,[]).extend(rows)
        inputs[arm]=dict(metrics=metrics,condition_groups=len(groups))
    tables['condition_metrics']=[dict(r,arm=arm) for arm,metrics in (('common',common),('longest',longest)) for r in metrics]
    tables['fingerprint_seed_summaries']=seed_rows;tables['E3_reference_checks']=references
    tables['E1_condition_metrics']=e1
    tables['original_seed_condition_metrics']=original_metrics
    tables['horizon_qualification']=horizon_qualifications(common,longest)
    longmap={metric_key(r):r for r in longest}
    tables['original_seed_sensitivity']=[dict(model=r['model'],C_label=r['C_label'],f=r['f'],view=r['view'],
        original_n_seeds=int(r['n_seeds']),all_long_n_seeds=int(longmap[metric_key(r)]['n_seeds']),
        cohort='Exact retained original Gaussian paths; same long horizon; not additional independent realizations',
        **interval_and_label_compatibility(r,longmap[metric_key(r)],label='original_seed_cohort_versus_all12_long')) for r in original_metrics]
    out=bundle/'results/final_analysis';out.mkdir(exist_ok=False)
    for name,rows in tables.items():write_csv(out/(name+'.csv'),rows)
    write_json(out/'decision_details.json',dict(results=results,inputs=inputs))
    source=checked_json(bundle/SOURCE_FREEZE)
    write_json(out/'complete.json',dict(passed=True,complete=True,generated_utc=datetime.now(timezone.utc).isoformat(),
        rule_freeze_sha256=sha(bundle/'rule_freeze.json'),source_freeze_sha256=sha(bundle/SOURCE_FREEZE),source_sha256=source['source_sha256'],
        input_spec_path=str(Path(spec_path).resolve()),input_spec_sha256=sha(spec_path),
        source_stage_gates={name:dataset['gates'] for name,dataset in datasets.items()},
        original_seed_sensitivity_gates=original_dataset['gates'],
        A5_receipts=A5_receipts,outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        decision_rows=len(tables['criteria']),whole_seed_bootstrap_resamples=2000,
        scope='Frozen H1/H2/H3 criterion tables for common and longest-qualified arms; all windows and forced views retained. A7 sensitivities remain explicit qualifications, not hidden criterion gates.',
        pending_independent_final_decision_audit=True,no_native_integrations=True))
    return dict(passed=True,complete=True,output=str(out),complete_sha256=sha(out/'complete.json'),pending_independent_final_decision_audit=True)


def verify_rule_freeze(bundle):
    receipt=json.loads((bundle/'rule_freeze.json').read_text())
    if not receipt.get('rules_frozen') or receipt.get('future_response_outcomes_or_correlations_loaded') is not False:raise ValueError('A pre-outcome rule freeze is required')
    for row in receipt['authorities']:
        if sha(row['frozen'])!=row['sha256']:raise ValueError('Frozen scientific decision authority changed')
    return receipt


def freeze_source(bundle):
    verify_rule_freeze(bundle)
    target=bundle/FROZEN_SOURCES;target.mkdir(exist_ok=False)
    sources={}
    audit_source=Path(__file__).with_name('signature_figure4_final_audit_20261020.py')
    for source in (Path(__file__),Path(A.__file__),audit_source):
        destination=target/source.name;shutil.copy2(source,destination);sources[str(destination.resolve())]=sha(destination)
    record=dict(generated_utc=datetime.now(timezone.utc).isoformat(),source=str((target/Path(__file__).name).resolve()),source_sha256=sources,
        source_version=3,supersedes_source_freeze_sha256=sha(bundle/'source_freeze_v2.json'),
        revision='Reviewed correction of resolved-or-bounded eligibility: use actual finite sub-two-cadence crossing bounds at adequate coverage, without inventing a point tau. A5 identity gates from v2 retained; original94d authority retained.',
        bound_eligibility_clarification=dict(path=str((bundle/'inputs/20261020_figure4_bound_eligibility_clarification.json').resolve()),sha256=sha(bundle/'inputs/20261020_figure4_bound_eligibility_clarification.json')),
        rule_freeze_sha256=sha(bundle/'rule_freeze.json'),future_scientific_evaluation_performed=False,
        accepted_E0_A6_was_loaded_after_v1_rule_and_source_freeze=True)
    (bundle/SOURCE_FREEZE).write_text(json.dumps(record,indent=2)+'\n')
    return record


def verify_source_freeze(bundle):
    verify_rule_freeze(bundle);record=checked_json(bundle/SOURCE_FREEZE)
    if record['rule_freeze_sha256']!=sha(bundle/'rule_freeze.json'):raise ValueError('Rule freeze changed after executable freeze')
    for path,digest in record['source_sha256'].items():
        if sha(path)!=digest:raise ValueError('The frozen final-decision executable changed')
    checked_json(record['bound_eligibility_clarification']['path'],record['bound_eligibility_clarification']['sha256'])
    actual=str(Path(__file__).resolve())
    if actual not in record['source_sha256'] or sha(actual)!=record['source_sha256'][actual]:
        raise ValueError('Run the actual frozen final-statistics source closure')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=('freeze-source','verify-rules','prepare-observations','evaluate'));p.add_argument('--bundle',type=Path,required=True);p.add_argument('--spec',type=Path);args=p.parse_args()
    result=freeze_source(args.bundle.resolve()) if args.command=='freeze-source' else integrated_evaluation(args.bundle.resolve(),args.spec) if args.command=='evaluate' else prepare_observations(args.bundle.resolve(),args.spec) if args.command=='prepare-observations' else verify_rule_freeze(args.bundle.resolve())
    print(json.dumps(result,indent=2,allow_nan=False))
