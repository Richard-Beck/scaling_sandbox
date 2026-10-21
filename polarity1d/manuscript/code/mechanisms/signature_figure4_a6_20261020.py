"""Predeclared A6 condition associations and held-out-model sign transfer.

Native times/lengths/rates are compared only within a model. Cross-model
predictors are dimensionless large/small quotients; direct and proxy outcomes
remain separate, with all four planned folds including unavailable folds.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import os
import shutil

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/"polarity1d/manuscript"
RUNS=ROOT/"data/runs"
DECISIONS=REPO/"docs/tasks/briefs/20261020_figure4_a6_decisions.json"
A5=RUNS/"signature_figure4_regions_strict_20261020/results/analysis"
RATIOS=RUNS/"signature_figure4_e0_corrected_basals_20261020/results/analysis"
DESTINATION=RUNS/"signature_figure4_a6_strict_20261020/results/analysis"
MODELS=['wave_pinning','goryachev','holmes_model3','spring']
PRIMARY=['multiple_front_fraction','region_count_mean']
RAW_TOTAL_COUNTS={'observed_births','observed_takeovers','observed_pairwise_takeovers','unscorable_takeover_tracks'}
EPSILON=1e-12


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def clean_json(value):
    if isinstance(value,dict):return {k:clean_json(v) for k,v in value.items()}
    if isinstance(value,list):return [clean_json(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    return value.item() if isinstance(value,np.generic) else value


def boolean(value):
    if isinstance(value,(bool,np.bool_)):return bool(value)
    if str(value).lower() in ('true','1'):return True
    if str(value).lower() in ('false','0'):return False
    raise ValueError('A canonical explicit boolean is required: '+repr(value))


def optional_boolean(value):
    return None if value is None or pd.isna(value) else boolean(value)


def numeric(value):
    return np.nan if value is None or pd.isna(value) else float(value)


def canonical(model):return str(model).removesuffix('_ring')


def rank_association(x,y,f):
    x,y,f=np.asarray(x,float),np.asarray(y,float),np.asarray(f,float)
    keep=np.isfinite(x)&np.isfinite(y)&np.isfinite(f)
    x,y,f=x[keep],y[keep],f[keep]
    result=dict(n=int(len(x)),distinct_f=int(len(np.unique(f))),spearman=None)
    if len(x)<3 or len(np.unique(f))<3:
        return dict(result,status='unavailable',reason='fewer_than_three_matched_f')
    for vector,name in [(x,'predictor'),(y,'outcome')]:
        if np.ptp(vector)<=EPSILON*max(1.,float(np.max(abs(vector)))):
            return dict(result,status='unavailable',reason='constant_'+name)
    rho=float(spearmanr(x,y).statistic)
    if not np.isfinite(rho):return dict(result,status='unavailable',reason='nonfinite_rank_correlation')
    return dict(result,status='observed',reason=None,spearman=0. if abs(rho)<=EPSILON else rho)


def quotient(large,small):
    if not np.isfinite(small) or not np.isfinite(large):return np.nan,'unavailable_feature_amount'
    if small<=0:return np.nan,'zero_or_nonpositive_small_denominator'
    return float(large/small),'observed'


def outcomes(row):
    direct=np.nan
    statuses=(row.get('tau_status_small'),row.get('tau_status_large'))
    coverage=optional_boolean(row.get('evaluable',False)) is True
    ratio=numeric(row.get('tau_ratio'));proxy=numeric(row.get('Dtheta_mobility_proxy'))
    if statuses==('resolved','resolved') and coverage and np.isfinite(ratio) and ratio>0:
        direct=float(np.log(ratio));direct_status='observed_resolved_direct_tau'
    else:direct_status='unavailable_censored_unresolved_or_low_coverage_direct_tau'
    proxy_valid=np.isfinite(proxy) and proxy>0 and numeric(row.get('mobility_proxy_valid_bootstrap_fraction'))>=.5
    return {'direct_tau':(direct,direct_status),
        'Dtheta_proxy':(float(np.log(proxy)) if proxy_valid else np.nan,
                        'observed_Dtheta_proxy_sensitivity' if proxy_valid else 'unavailable_Dtheta_proxy')}


def leave_one_model_out(points,primary=False):
    """Four sign-transfer folds; an unavailable model stays in the denominator."""
    rows=[]
    for model in MODELS:
        train=points[points.model!=model];test=points[points.model==model]
        tr=rank_association(train.feature_ratio,train.log_persistence_ratio,train.f)
        te=rank_association(test.feature_ratio,test.log_persistence_ratio,test.f)
        if tr['status']!='observed' or te['status']!='observed':
            status='unavailable';success=False
        else:
            success=tr['spearman']<0 and te['spearman']<0 if primary else tr['spearman']*te['spearman']>0
            status='success' if success else 'failure'
        rows.append(dict(held_out_model=model,status=status,sign_correct=success,
            expected_sign='negative' if primary else 'training_sign_exploratory',
            training_spearman=tr['spearman'],training_n=tr['n'],training_distinct_f=tr['distinct_f'],training_reason=tr['reason'],
            held_out_spearman=te['spearman'],held_out_n=te['n'],held_out_distinct_f=te['distinct_f'],held_out_reason=te['reason']))
    return rows


def primary_criterion(rho,folds):
    successes=sum(r['sign_correct'] for r in folds)
    threshold=rho is not None and rho<=-.6+EPSILON
    return dict(pooled_negative_threshold_pass=threshold,successful_folds=successes,planned_folds=4,
        unavailable_folds=sum(r['status']=='unavailable' for r in folds),
        A6_primary_predictor_support=bool(threshold and successes>=3))


def source_inputs(a5,ratios):
    a5gate=read(a5/'complete.json');ratiogate=read(ratios/'status.json')
    assert a5gate['passed'] and ratiogate['passed'] and ratiogate['complete']
    assert sha(a5/'metrics.csv')==a5gate['outputs']['metrics']['sha256']
    assert sha(ratios/'size_ratios.csv')==ratiogate['outputs']['size_ratios.csv']
    features=pd.read_csv(a5/'metrics.csv',float_precision='round_trip')
    responses=pd.read_csv(ratios/'size_ratios.csv',float_precision='round_trip')
    return features,responses,dict(A5_complete_sha256=sha(a5/'complete.json'),A5_metrics_sha256=sha(a5/'metrics.csv'),
        ratio_status_sha256=sha(ratios/'status.json'),size_ratios_sha256=sha(ratios/'size_ratios.csv'))


def matched_pair_context(features,response,experiment):
    """Match each size to its declared arm/horizon before seeing associations."""
    fields=['native_horizon_small','native_horizon_large']
    intentional=response.get('comparison_status')=='unavailable_unmatched_designed_input'
    missing_size=numeric(response.get('missing_size'))
    if intentional and (optional_boolean(response.get('evaluable')) is not False or not np.isfinite(missing_size)):
        raise ValueError('Intentional unmatched input needs an explicit missing size and evaluable=False')
    declared=all(k in response and pd.notna(response[k]) for k in fields)
    if experiment!='E0' and (response.get('analysis_arm')!=experiment or (not declared and not intentional)):
        raise ValueError('Expanded A6 needs explicit matching analysis_arm and both native horizons')
    if 'analysis_arm' in response and pd.notna(response['analysis_arm']) and response['analysis_arm']!=experiment:
        raise ValueError('Ratio analysis arm disagrees with requested A6 arm')
    sizes=[response['small_C_label'],response['large_C_label']]
    part=features[(features.model==response['model'])&np.isclose(features.f,response['f'],rtol=0,atol=1e-14)&
        (np.isclose(features.C_label,sizes[0],rtol=0,atol=1e-10)|np.isclose(features.C_label,sizes[1],rtol=0,atol=1e-10))]
    if intentional and not any(np.isclose(missing_size,size,rtol=0,atol=1e-10) for size in sizes):
        raise ValueError('Intentional missing size must be one member of the declared size pair')
    horizons=[None if pd.isna(response.get(k)) or response.get(k) is None else float(response[k]) for k in fields]
    if declared or experiment!='E0':
        if 'native_horizon' not in part:raise ValueError('A5 native horizons are required for expanded A6')
        for size,horizon in zip(sizes,horizons):
            observed=part[np.isclose(part.C_label,size,rtol=0,atol=1e-10)].native_horizon.unique()
            if intentional and np.isclose(size,missing_size,rtol=0,atol=1e-10):
                if len(observed) or horizon is not None:raise ValueError('Intentionally missing input must remain unavailable in A5 and its ratio horizon')
                continue
            if horizon is None or not np.isfinite(horizon) or horizon<=0:raise ValueError('Declared native horizons must be positive and finite')
            if len(observed)!=1 or not np.isfinite(observed[0]) or not np.isclose(observed[0],horizon,rtol=0,atol=1e-8):
                raise ValueError('A5 and ratio native horizon mismatch or mixed horizon rows: '+str((response['model'],response['f'],size)))
        if experiment=='common' and declared and not np.isclose(*horizons,rtol=0,atol=1e-8):
            raise ValueError('Common A6 requires equal observation horizons within each model size pair')
    return part,dict(analysis_arm=experiment,native_horizon_small=horizons[0],native_horizon_large=horizons[1],
        matched_horizons_across_sizes=bool(np.isclose(*horizons,rtol=0,atol=1e-8)) if declared else None,
        comparison_status=response.get('comparison_status','matched'),intentionally_missing_size=missing_size if intentional else None)


def condition_pairs(features,responses,*,experiment='E0',sigma=1.,hold_takeover=False):
    features=features.copy();responses=responses.copy()
    unknown=set(features.kind.unique())-{'noise','mean','zero','basal'}
    if unknown:
        raise ValueError('A5 condition-kind identity is invalid: '+repr(sorted(unknown)))
    features['model']=features.model.map(canonical);responses['model']=responses.model.map(canonical)
    features=features[(features.kind=='noise')&(features.experiment==experiment)&np.isclose(features.sigma,sigma)]
    responses=responses[responses.view=='main']
    if responses.duplicated(['model','f']).any():raise ValueError('One independent response condition per model/f and arm is required')
    if not len(features):raise ValueError('No canonical matching A5 forcing family')
    for r in responses.to_dict('records'):
        if 'sigma' in r and pd.notna(r['sigma']):
            if not np.isclose(float(r['sigma']),sigma,rtol=0,atol=1e-12):raise ValueError('A5 and response forcing sigma mismatch')
        elif experiment!='E0':raise ValueError('Expanded ratio tables must declare their forcing sigma')
    features['smoothed']=features.smoothed.map(boolean)
    features['primary']=features.primary.map(boolean)
    features['cadence_multiple']=features.observation_cadence/np.where(features.model=='spring',8.,80.)
    # Fail before fitting any association if a complete canonical noise view
    # is absent or mislabeled. Never infer a missing kind from other columns.
    for r in responses.to_dict('records'):
        pair,_=matched_pair_context(features,r,experiment);candidates=pair[pair.primary]
        for size in [r['small_C_label'],r['large_C_label']]:
            if r.get('comparison_status')=='unavailable_unmatched_designed_input' and np.isclose(size,numeric(r.get('missing_size')),rtol=0,atol=1e-10):continue
            for statistic in PRIMARY:
                required=candidates[np.isclose(candidates.C_label,size,rtol=0,atol=1e-10)&(candidates.statistic==statistic)]
                if len(required)!=1:
                    raise ValueError('Canonical primary A5 noise view missing or duplicated: '+str((r['model'],r['f'],size,statistic)))
    # Common physical imaging cadence is80 chemical /8 Spring; its quotient
    # groups equivalent observational views without pooling native time units.
    rows=[]
    for r in responses.to_dict('records'):
        part,pair_context=matched_pair_context(features,r,experiment)
        if not len(part):raise ValueError('A declared response condition lacks A5 measurements: '+str((r['model'],r['f'])))
        for key,group in part.groupby(['threshold','smoothed','cadence_multiple','statistic'],dropna=False):
            threshold,smoothed,cadence,statistic=key
            small=group[np.isclose(group.C_label,r['small_C_label'],rtol=0,atol=1e-10)]
            large=group[np.isclose(group.C_label,r['large_C_label'],rtol=0,atol=1e-10)]
            if len(small)>1 or len(large)>1:raise ValueError('Duplicate small or large A5 observational view')
            s,l=small.iloc[0] if len(small) else pd.Series(dtype=object),large.iloc[0] if len(large) else pd.Series(dtype=object)
            sv,lv=numeric(s.get('estimate')),numeric(l.get('estimate'))
            matched=len(small)==1 and len(large)==1
            feature_status='observed_condition_mean' if matched else 'unavailable_unmatched_observation_view'
            view_pair_status='matched' if matched else 'unavailable_small_view' if not len(small) else 'unavailable_large_view'
            if matched and hold_takeover and 'takeover' in statistic:
                sv,lv=np.nan,np.nan;feature_status='pending_strict_global_takeover_replay'
            elif matched and statistic=='observed_takeovers_per_birth':
                unscorable=part[(part.threshold==threshold)&(part.smoothed==smoothed)&(part.cadence_multiple==cadence)&
                    (part.statistic=='ever_unscorable_takeover_fraction')&part.C_label.isin([r['small_C_label'],r['large_C_label']])]
                if len(unscorable)!=2 or unscorable.estimate.isna().any() or (unscorable.estimate>0).any():
                    sv,lv=np.nan,np.nan;feature_status='unavailable_takeover_probability_with_unscorable_identity_or_amount'
            ratio,ratio_status=quotient(lv,sv)
            if not matched:ratio=np.nan;ratio_status='unavailable_unmatched_observation_view'
            unequal_exposure=(statistic in RAW_TOTAL_COUNTS and pair_context['matched_horizons_across_sizes'] is False)
            if unequal_exposure:
                ratio=np.nan;ratio_status='unavailable_unequal_horizon_exposure_raw_total_count'
            forced=optional_boolean(r.get('forced_combined') if 'forced_combined' in r else r.get('forced'))
            for outcome,(value,status) in outcomes(r).items():
                rows.append(dict(model=r['model'],f=r['f'],experiment=experiment,sigma=sigma,**pair_context,
                    threshold=threshold,smoothed=bool(smoothed),cadence_multiple=cadence,
                    imaging_view='primary' if s.get('primary',False) or l.get('primary',False) else f'threshold{threshold:g}_smooth{int(smoothed)}_cadence{cadence:g}',
                    observational_pair_status=view_pair_status,
                    statistic=statistic,small_C_label=r['small_C_label'],large_C_label=r['large_C_label'],
                    raw_total_count=statistic in RAW_TOTAL_COUNTS,raw_count_exposure_guard=unequal_exposure,
                    small_estimate=sv,large_estimate=lv,feature_ratio=ratio,feature_status=feature_status,feature_ratio_status=ratio_status,
                    small_ci_lower=s.get('ci_lower'),small_ci_upper=s.get('ci_upper'),large_ci_lower=l.get('ci_lower'),large_ci_upper=l.get('ci_upper'),
                    seeds_small=s.get('independent_seeds'),seeds_large=l.get('independent_seeds'),
                    native_profile_cadence_small=s.get('native_profile_cadence'),native_profile_cadence_large=l.get('native_profile_cadence'),
                    seeds_with_estimand_small=s.get('seeds_with_estimand'),seeds_with_estimand_large=l.get('seeds_with_estimand'),
                    outcome=outcome,log_persistence_ratio=value,outcome_status=status,forced=forced,
                    forcing_status='unavailable_unknown' if forced is None else 'known_forced' if forced else 'known_nonforced',
                    tau_status_small=r.get('tau_status_small'),tau_status_large=r.get('tau_status_large'),
                    source_tau_ratio_lower=r.get('tau_ratio_lower'),source_tau_ratio_upper=r.get('tau_ratio_upper'),
                    source_tau_ratio_ci_lower=r.get('tau_ratio_ci_lower'),source_tau_ratio_ci_upper=r.get('tau_ratio_ci_upper'),
                    source_Dtheta_proxy=r.get('Dtheta_mobility_proxy'),source_Dtheta_proxy_ci_lower=r.get('Dtheta_mobility_proxy_ci_lower'),
                    source_Dtheta_proxy_ci_upper=r.get('Dtheta_mobility_proxy_ci_upper')))
    return pd.DataFrame(rows)


def association_tables(points):
    if points.analysis_arm.nunique()!=1 or points.sigma.nunique()!=1:
        raise ValueError('Separate A6 runs are required for each analysis arm and forcing sigma')
    within,cross,folds,primary=[],[],[],[]
    for key,allpoints in points.groupby(['statistic','imaging_view','outcome'],dropna=False):
        statistic,view,outcome=key
        for inclusive in [False,True]:
            selected=allpoints if inclusive else allpoints[allpoints.forcing_status=='known_nonforced']
            horizons={m:sorted({('not_declared_or_unmatched' if pd.isna(r.native_horizon_small) or pd.isna(r.native_horizon_large) else f'{r.native_horizon_small:g}/{r.native_horizon_large:g}')
                for r in p.itertuples()}) for m,p in allpoints.groupby('model')}
            context=dict(statistic=statistic,imaging_view=view,outcome=outcome,forced_inclusive=inclusive,
                excluded_forced=int((allpoints.forcing_status=='known_forced').sum()) if not inclusive else 0,
                excluded_unknown_forcing=int((allpoints.forcing_status=='unavailable_unknown').sum()) if not inclusive else 0,
                primary_H2_feature=statistic in PRIMARY,analysis_arm=points.analysis_arm.iloc[0],sigma=float(points.sigma.iloc[0]),
                native_horizon_pairs_by_model=json.dumps(horizons,sort_keys=True),
                heterogeneous_size_horizons=bool(any(pd.notna(r.native_horizon_small) and pd.notna(r.native_horizon_large) and not r.matched_horizons_across_sizes for r in allpoints.itertuples())))
            for model in MODELS:
                part=selected[(selected.model==model)&(selected.observational_pair_status=='matched')]
                for transform,column in [('small','small_estimate'),('large','large_estimate'),('large_over_small','feature_ratio')]:
                    rank=rank_association(part[column],part.log_persistence_ratio,part.f)
                    within.append(dict(context,model=model,feature_transform=transform,**rank))
            valid=selected[(selected.observational_pair_status=='matched')&np.isfinite(selected.feature_ratio)&np.isfinite(selected.log_persistence_ratio)]
            pooled=rank_association(valid.feature_ratio,valid.log_persistence_ratio,valid.f)
            these_folds=leave_one_model_out(valid,statistic in PRIMARY)
            folds.extend(dict(context,**fold) for fold in these_folds)
            criterion=primary_criterion(pooled['spearman'],these_folds) if statistic in PRIMARY else dict(
                successful_folds=sum(r['sign_correct'] for r in these_folds),planned_folds=4,
                unavailable_folds=sum(r['status']=='unavailable' for r in these_folds),
                A6_primary_predictor_support=None,pooled_negative_threshold_pass=None)
            cross.append(dict(context,**pooled,**criterion,predictor='within_model_large_over_small_only',
                models_with_usable_points=int(valid.model.nunique()),usable_points_by_model=json.dumps(valid.groupby('model').size().to_dict(),sort_keys=True)))
            if statistic in PRIMARY:
                primary.append(dict(context,pooled_spearman=pooled['spearman'],pooled_status=pooled['status'],
                    pooled_reason=pooled['reason'],**criterion,
                    primary_view=bool(view=='primary' and outcome=='direct_tau' and not inclusive)))
    return dict(within_model=pd.DataFrame(within),cross_model=pd.DataFrame(cross),LOMO_folds=pd.DataFrame(folds),H2_primary_features=pd.DataFrame(primary))


def analyze(a5,ratios,output,*,experiment='E0',sigma=1.,rules_path=DECISIONS,hold_takeover=False):
    if output.exists():raise FileExistsError('A6 exports use a fresh immutable output directory')
    rules=read(rules_path);assert rules['status']=='predeclared_before_loading_A6_correlations'
    features,responses,sources=source_inputs(a5,ratios)
    points=condition_pairs(features,responses,experiment=experiment,sigma=sigma,hold_takeover=hold_takeover)
    tables=association_tables(points)
    output.mkdir(parents=True)
    shutil.copy2(rules_path,output/'predeclared_rules.json');shutil.copy2(__file__,output/'frozen_analysis.py')
    points.to_csv(output/'condition_pairs.csv',index=False,float_format='%.17g')
    for name,frame in tables.items():frame.to_csv(output/(name+'.csv'),index=False,float_format='%.17g')
    primary=tables['H2_primary_features'];main=primary[primary.primary_view]
    assert len(main)==2
    summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),complete=True,passed=True,native_integrations=0,
        analysis_arm=experiment,sigma=sigma,
        native_horizon_pairs=points[['model','f','native_horizon_small','native_horizon_large','matched_horizons_across_sizes']].drop_duplicates().to_dict('records'),
        source_sha256=sha(Path(__file__)),predeclared_rules_sha256=sha(rules_path),sources=sources,
        primary_A6_support=bool(main.A6_primary_predictor_support.any()),
        primary_features=main.where(pd.notnull(main),None).to_dict('records'),
        planned_models=MODELS,rank_uncertainty_scope=rules['rank_uncertainty_scope'],
        strict_global_takeover_features_pending=hold_takeover,
        scientific_scope='Associations of observed-region condition means and held-out-model sign transfer. E0 only unless an explicit matched forcing-family input is supplied; no causal or cell-velocity inference. Direct tau and censored Dtheta-proxy sensitivity remain separate.',
        outputs={p.name:sha(p) for p in output.glob('*.csv')})
    # DataFrames may preserve float NaNs when converting mixed numeric rows.
    summary=clean_json(summary)
    write(output/'complete.json',summary)
    return dict(output=str(output),primary_A6_support=summary['primary_A6_support'],primary_features=summary['primary_features'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--a5',type=Path,default=A5);parser.add_argument('--ratios',type=Path,default=RATIOS)
    parser.add_argument('--output',type=Path,default=DESTINATION);parser.add_argument('--experiment',default='E0');parser.add_argument('--sigma',type=float,default=1.)
    parser.add_argument('--rules',type=Path,default=DECISIONS)
    parser.add_argument('--hold-takeover',action='store_true',help='Retain unavailable takeover rows while corrected strict-global tracking is pending; never use pairwise-only values as primary global takeovers.')
    args=parser.parse_args();print(json.dumps(analyze(args.a5.resolve(),args.ratios.resolve(),args.output.resolve(),experiment=args.experiment,sigma=args.sigma,rules_path=args.rules.resolve(),hold_takeover=args.hold_takeover)))


if __name__=='__main__':main()
