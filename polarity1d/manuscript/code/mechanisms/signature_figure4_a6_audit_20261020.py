"""Independent A6 source joins and rank/fold audit; no estimator imports.

Average ranks and Pearson correlation are reconstructed directly instead of
calling the production scipy Spearman function. This checks derived A6 tables,
not native region-tracking identities, which retain their own native audits.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import math
import numpy as np
import pandas as pd

MODELS=('wave_pinning','goryachev','holmes_model3','spring')
PRIMARY=('multiple_front_fraction','region_count_mean')
RAW_COUNTS={'observed_births','observed_takeovers','observed_pairwise_takeovers','unscorable_takeover_tracks'}
EPS=1e-12


def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
    return h.hexdigest()
def number(value):return None if value is None or pd.isna(value) else float(value)
def boolean(value):
    if value is None or pd.isna(value):return None
    if str(value).lower() in ('true','1','1.0'):return True
    if str(value).lower() in ('false','0','0.0'):return False
    raise ValueError('Explicit canonical boolean required')
def same(actual,expected):
    a,b=number(actual),number(expected)
    if a is None or b is None:
        if a is not b:raise ValueError('Unavailable value replaced by a point')
    elif not math.isclose(a,b,rel_tol=4e-12,abs_tol=2e-13):raise ValueError('Independent numeric reconstruction differs: '+repr((a,b)))


def exact_source_value(actual,expected):
    a,b=number(actual),number(expected)
    if a!=b:raise ValueError('Canonical source value did not roundtrip exactly: '+repr((a,b)))


def average_ranks(values):
    v=np.asarray(values,float);order=np.argsort(v,kind='mergesort');ranks=np.empty(len(v),float);start=0
    while start<len(v):
        stop=start+1
        while stop<len(v) and v[order[stop]]==v[order[start]]:stop+=1
        ranks[order[start:stop]]=(start+stop+1)/2.
        start=stop
    return ranks


def independent_rank(x,y,f):
    a,b,inputs=[np.asarray(v,float) for v in (x,y,f)]
    keep=np.isfinite(a)&np.isfinite(b)&np.isfinite(inputs);a,b,inputs=a[keep],b[keep],inputs[keep]
    result=dict(n=len(a),distinct_f=len(set(inputs.tolist())),spearman=None)
    if result['n']<3 or result['distinct_f']<3:return dict(result,status='unavailable',reason='fewer_than_three_matched_f')
    for value,name in [(a,'predictor'),(b,'outcome')]:
        if float(value.max()-value.min())<=EPS*max(1.,float(abs(value).max())):
            return dict(result,status='unavailable',reason='constant_'+name)
    ar,br=average_ranks(a),average_ranks(b);ar-=ar.mean();br-=br.mean()
    rho=float(ar@br/math.sqrt(float(ar@ar)*float(br@br)))
    return dict(result,status='observed',reason=None,spearman=0. if abs(rho)<=EPS else rho)


def audit_source_points(points,features,responses,arm,sigma):
    features=features.copy();responses=responses.copy()
    for frame in (features,responses):frame['model']=frame.model.map(lambda value:str(value).removesuffix('_ring'))
    features=features[(features.kind=='noise')&(features.experiment==arm)&np.isclose(features.sigma,sigma)]
    responses=responses[responses.view=='main']
    key=lambda r:(r['model'],float(r['C_label']),float(r['f']),float(r['threshold']),boolean(r['smoothed']),
                  float(r['observation_cadence']),r['statistic'])
    lookup={}
    for row in features.to_dict('records'):
        k=key(row)
        if k in lookup:raise ValueError('Duplicate source A5 measurement')
        lookup[k]=row
    response_lookup={(r['model'],float(r['f'])):r for r in responses.to_dict('records')}
    if len(response_lookup)!=len(responses):raise ValueError('Duplicate response condition')
    seen=set();unmatched=unequal=unknown_forcing=0
    for p in points.to_dict('records'):
        condition=(p['model'],float(p['f']))
        if condition not in response_lookup:raise ValueError('A6 condition has no canonical response source')
        r=response_lookup[condition];dt=float(p['cadence_multiple'])*(8. if p['model']=='spring' else 80.)
        sides=[]
        for side in ('small','large'):
            k=(p['model'],float(p[side+'_C_label']),float(p['f']),float(p['threshold']),boolean(p['smoothed']),dt,p['statistic'])
            source=lookup.get(k);sides.append(source)
            if source:
                if arm!='E0':same(source['native_horizon'],r.get('native_horizon_'+side))
                same(source.get('native_profile_cadence'),p.get('native_profile_cadence_'+side))
            seen.add((k,p['outcome']))
        s,l=sides;matched=s is not None and l is not None
        primary=any(boolean(source.get('primary')) is True for source in sides if source is not None)
        if (p['imaging_view']=='primary')!=primary:raise ValueError('A sensitivity observation was mislabeled primary')
        if boolean(p.get('raw_count_exposure_guard')):
            if p['statistic'] not in RAW_COUNTS or number(p.get('native_horizon_small'))==number(p.get('native_horizon_large')):
                raise ValueError('Invalid raw-count exposure guard')
            unequal+=1
        sv=number(s.get('estimate')) if s else None;lv=number(l.get('estimate')) if l else None
        if matched and p['statistic']=='observed_takeovers_per_birth':
            ambiguity=[]
            for source in sides:
                akey=(p['model'],float(source['C_label']),float(p['f']),float(p['threshold']),boolean(p['smoothed']),dt,'ever_unscorable_takeover_fraction')
                value=number(lookup.get(akey,{}).get('estimate'));ambiguity.append(value)
            if any(v is None or v>0 for v in ambiguity):sv=lv=None
        if p['feature_status']=='pending_strict_global_takeover_replay':raise ValueError('Final independent A6 audit rejects pending takeover features')
        exact_source_value(p['small_estimate'],sv);exact_source_value(p['large_estimate'],lv)
        ratio=lv/sv if matched and sv is not None and lv is not None and sv>0 else None
        if p['statistic'] in RAW_COUNTS and number(p.get('native_horizon_small'))!=number(p.get('native_horizon_large')):ratio=None
        same(p['feature_ratio'],ratio)
        forced=boolean(r.get('forced_combined',r.get('forced')))
        if boolean(p['forced']) is not forced:raise ValueError('A6 forcing identity differs from canonical source')
        forcing_status='unavailable_unknown' if forced is None else 'known_forced' if forced else 'known_nonforced'
        if p['forcing_status']!=forcing_status:raise ValueError('Unknown forcing was converted into nonforced')
        unknown_forcing+=int(forced is None)
        tau=number(r.get('tau_ratio'));proxy=number(r.get('Dtheta_mobility_proxy'))
        if p['outcome']=='direct_tau':
            available=r.get('tau_status_small')=='resolved' and r.get('tau_status_large')=='resolved' and boolean(r.get('evaluable')) is True and tau is not None and tau>0
            expected=math.log(tau) if available else None
        else:
            fraction=number(r.get('mobility_proxy_valid_bootstrap_fraction'))
            expected=math.log(proxy) if proxy is not None and proxy>0 and fraction is not None and fraction>=.5 else None
        same(p['log_persistence_ratio'],expected)
        if not matched:
            unmatched+=1
            if p['observational_pair_status']=='matched':raise ValueError('Unmatched source view labeled matched')
    required={(key(r),outcome) for r in features.to_dict('records') if
              (r['model'],float(r['f'])) in response_lookup and any(math.isclose(float(r['C_label']),float(response_lookup[(r['model'],float(r['f']))][side+'_C_label']),abs_tol=1e-10,rel_tol=0) for side in ('small','large'))
              for outcome in ('direct_tau','Dtheta_proxy')}
    if not required<=seen:raise ValueError('A6 silently omitted source observational views')
    return dict(source_point_rows=len(points),source_native_views=len(lookup),all_source_observation_views_retained=True,
        unmatched_view_rows=unmatched,unequal_raw_count_rows=unequal,unknown_forcing_rows=unknown_forcing)


def verify_rank(row,expected):
    same(row['spearman'],expected['spearman'])
    if row['status']!=expected['status'] or int(row['n'])!=expected['n'] or int(row['distinct_f'])!=expected['distinct_f']:
        raise ValueError('Independent rank coverage/status differs')
    reason=None if pd.isna(row.get('reason')) else row['reason']
    if reason!=expected['reason']:raise ValueError('Independent unavailable-rank reason differs')


def audit_ranks(points,cross,within,folds,primary):
    ckeys=['statistic','imaging_view','outcome','forced_inclusive']
    cross_map={tuple(r[k] for k in ckeys):r for r in cross.to_dict('records')}
    within_map={tuple(r[k] for k in ckeys+['model','feature_transform']):r for r in within.to_dict('records')}
    fold_map={tuple(r[k] for k in ckeys+['held_out_model']):r for r in folds.to_dict('records')}
    primary_map={tuple(r[k] for k in ckeys):r for r in primary.to_dict('records')}
    if len(cross_map)!=len(cross) or len(within_map)!=len(within) or len(fold_map)!=len(folds) or len(primary_map)!=len(primary):
        raise ValueError('Duplicate association/fold rows')
    counts=CounterLike();primary_results=[]
    for group,raw in points.groupby(['statistic','imaging_view','outcome'],dropna=False):
        statistic,view,outcome=group
        for inclusive in (False,True):
            key=(*group,inclusive)
            selected=raw if inclusive else raw[raw.forcing_status=='known_nonforced']
            valid=selected[(selected.observational_pair_status=='matched')&np.isfinite(selected.feature_ratio)&np.isfinite(selected.log_persistence_ratio)]
            expected=independent_rank(valid.feature_ratio,valid.log_persistence_ratio,valid.f)
            verify_rank(cross_map[key],expected);counts['cross']+=1
            successes=unavailable=0
            for model in MODELS:
                training=valid[valid.model!=model];testing=valid[valid.model==model]
                train=independent_rank(training.feature_ratio,training.log_persistence_ratio,training.f)
                test=independent_rank(testing.feature_ratio,testing.log_persistence_ratio,testing.f)
                available=train['status']==test['status']=='observed'
                correct=(train['spearman']<0 and test['spearman']<0 if statistic in PRIMARY else train['spearman']*test['spearman']>0) if available else False
                status='success' if correct else 'failure' if available else 'unavailable'
                actual=fold_map[(*key,model)]
                if actual['status']!=status or boolean(actual['sign_correct']) is not correct:raise ValueError('Independent fold sign/availability differs')
                same(actual['training_spearman'],train['spearman']);same(actual['held_out_spearman'],test['spearman'])
                if int(actual['training_n'])!=train['n'] or int(actual['held_out_n'])!=test['n']:raise ValueError('Independent fold sample coverage differs')
                successes+=int(correct);unavailable+=int(not available);counts['folds']+=1
                part=selected[(selected.model==model)&(selected.observational_pair_status=='matched')]
                for transform,column in [('small','small_estimate'),('large','large_estimate'),('large_over_small','feature_ratio')]:
                    rank=independent_rank(part[column],part.log_persistence_ratio,part.f)
                    verify_rank(within_map[(*key,model,transform)],rank);counts['within']+=1
            pooled=cross_map[key]
            if int(pooled['successful_folds'])!=successes or int(pooled['unavailable_folds'])!=unavailable or int(pooled['planned_folds'])!=4:
                raise ValueError('All four folds must remain in the denominator')
            if statistic in PRIMARY:
                threshold=expected['spearman'] is not None and expected['spearman']<=-.6+EPS
                supported=bool(threshold and successes>=3);p=primary_map[key]
                same(p['pooled_spearman'],expected['spearman'])
                if boolean(p['pooled_negative_threshold_pass'])!=threshold or boolean(p['A6_primary_predictor_support'])!=supported:raise ValueError('Primary H2 criterion differs')
                if int(p['successful_folds'])!=successes or int(p['unavailable_folds'])!=unavailable or int(p['planned_folds'])!=4:raise ValueError('Primary fold denominator differs')
                main=view=='primary' and outcome=='direct_tau' and not inclusive
                if boolean(p['primary_view'])!=main:raise ValueError('A sensitivity view was relabeled primary')
                if main:primary_results.append(dict(statistic=statistic,pooled_spearman=expected['spearman'],successful_folds=successes,unavailable_folds=unavailable,A6_primary_predictor_support=supported))
    if counts['cross']!=len(cross) or counts['within']!=len(within) or counts['folds']!=len(folds):raise ValueError('Unexpected/omitted derived association rows')
    if len(primary_results)!=2:raise ValueError('Exactly two predeclared primary H2 features required')
    return dict(independent_associations=counts,independent_primary_rank_checks=primary_results,all_planned_folds_retained=True)


class CounterLike(dict):
    def __missing__(self,key):return 0


def audit(analysis,a5,ratios,output):
    target=output or analysis/'independent_A6_audit.json'
    if target.exists():raise FileExistsError('Independent A6 receipts are immutable')
    gate=read(analysis/'complete.json');source=gate['sources'];a5gate=read(a5/'complete.json');ratio_gate=read(ratios/'status.json')
    if not gate.get('passed') or not gate.get('complete') or gate.get('strict_global_takeover_features_pending'):raise ValueError('Complete strict A6 required')
    if a5gate.get('passed') is not True or ratio_gate.get('passed') is not True or ratio_gate.get('complete') is not True:raise ValueError('Complete accepted source gates required')
    if a5gate['outputs']['metrics']['sha256']!=source['A5_metrics_sha256'] or ratio_gate['outputs']['size_ratios.csv']!=source['size_ratios_sha256']:
        raise ValueError('A6 source bindings differ from their upstream complete gates')
    if sha(analysis/'predeclared_rules.json')!=gate['predeclared_rules_sha256']:raise ValueError('Frozen A6 rule copy changed')
    for path,expected in [(a5/'complete.json',source['A5_complete_sha256']),(a5/'metrics.csv',source['A5_metrics_sha256']),
                          (ratios/'status.json',source['ratio_status_sha256']),(ratios/'size_ratios.csv',source['size_ratios_sha256'])]:
        if sha(path)!=expected:raise ValueError('Canonical A6 source changed')
    tables={}
    for name,expected in gate['outputs'].items():
        if sha(analysis/name)!=expected:raise ValueError('A6 output hash mismatch')
        tables[Path(name).stem]=pd.read_csv(analysis/name,float_precision='round_trip')
    sources=audit_source_points(tables['condition_pairs'],pd.read_csv(a5/'metrics.csv',float_precision='round_trip'),pd.read_csv(ratios/'size_ratios.csv',float_precision='round_trip'),gate.get('analysis_arm','E0'),gate.get('sigma',1.))
    ranks=audit_ranks(tables['condition_pairs'],tables['cross_model'],tables['within_model'],tables['LOMO_folds'],tables['H2_primary_features'])
    result=dict(passed=True,complete=True,generated_utc=datetime.now(timezone.utc).isoformat(),analysis_arm=gate.get('analysis_arm','E0'),sigma=gate.get('sigma',1.),
        strict_global_A5_input=True,no_pending_takeover_replay=True,source_complete_sha256=sha(analysis/'complete.json'),sources=source,
        output_hashes_verified=len(tables),**sources,**ranks,audit_source_sha256=sha(__file__),
        rank_algorithm='Stable sorting, exact tied-value average ranks, centred Pearson correlation of those ranks; production scipy Spearman is not imported.',
        scope='Independent canonical-source joins, explicit forcing/horizon/cadence/unavailable guards and all associations/folds. Native field qualification and latent tracking identities remain scoped to their separate upstream audits.')
    target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--a5',type=Path,required=True);p.add_argument('--ratios',type=Path,required=True);p.add_argument('--output',type=Path)
    a=p.parse_args();print(json.dumps(audit(a.analysis.resolve(),a.a5.resolve(),a.ratios.resolve(),a.output.resolve() if a.output else None)))


if __name__=='__main__':main()
