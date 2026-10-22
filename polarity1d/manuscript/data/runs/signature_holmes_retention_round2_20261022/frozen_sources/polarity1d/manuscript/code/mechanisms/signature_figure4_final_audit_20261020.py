"""Independent final Figure 4 criterion and whole-seed fingerprint checks.

No producer helper is imported. Native scalar estimands and A6 fits already
have separate independent audits; this checks their final joins and decisions.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timezone
import numpy as np

MODELS=('wave_pinning','goryachev','holmes_model3','spring_ring')
VIEWS=('main','analysis_first_half','analysis_second_half','late_quarter','cadence_x2','cadence_x4','cadence_x8')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4<<20),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def table(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))


def numeric(value):
    if value is None or value=='':return None
    result=float(value)
    return None if math.isnan(result) else result


def flag(value):
    return True if value is True or value=='True' else False if value is False or value=='False' else None


def bounds(row,stem):
    lo=numeric(row.get(stem+'_lower'));hi=numeric(row.get(stem+'_upper'))
    if flag(row.get(stem+'_upper_unbounded')):hi=math.inf
    return (lo,hi) if lo is not None and hi is not None and lo<=hi else None


def combine(values):
    return 'fail' if 'fail' in values else 'inconclusive' if not values or 'inconclusive' in values else 'pass'


def count_status(values,required):
    passed=values.count('pass');unknown=values.count('inconclusive')
    return 'pass' if passed>=required else 'fail' if passed+unknown<required else 'inconclusive'


def allowed(row,inclusive):return inclusive or flag(row.get('forced_combined')) is False
def eligible(row):
    observed=bounds(row,'tau')
    if (numeric(row.get('tau_coverage')) or 0)<.5 or observed is None:return False
    status=row.get('tau_status')
    if status in ('resolved','right_censored'):return True
    if status!='unresolved_below_two_cadences' or numeric(row.get('tau')) is not None:return False
    cadence=numeric(row.get('cadence'));bracket=row.get('tau_crossing_bracket')
    if isinstance(bracket,str):
        try:bracket=json.loads(bracket)
        except (ValueError,TypeError):return False
    if cadence is None or not math.isfinite(cadence) or cadence<=0 or not isinstance(bracket,(list,tuple)) or len(bracket)!=2:return False
    try:lo,hi=map(float,bracket)
    except (ValueError,TypeError):return False
    return observed[0]==0 and math.isfinite(observed[1]) and math.isclose(observed[1],2*cadence,rel_tol=1e-12,abs_tol=1e-12) and math.isfinite(lo) and math.isfinite(hi) and 0<=lo<hi<=2*cadence+1e-12
def comparison(row):return flag(row.get('evaluable')) is True and min(numeric(row.get('coverage_small')) or 0,numeric(row.get('coverage_large')) or 0)>=.5


def h1(rows,small,large,inclusive):
    selected=[r for r in rows if comparison(r) and allowed(r,inclusive)]
    unknown=not inclusive and any(comparison(r) and flag(r.get('forced_combined')) is None for r in rows)
    cis=[bounds(r,'tau_ratio_ci') for r in selected];reverse=any(ci and ci[1]<1 for ci in cis)
    gain='pass' if selected and not unknown and all(ci and ci[0]>1 for ci in cis) else 'fail' if reverse else 'inconclusive'
    noreverse='fail' if reverse else 'inconclusive' if unknown or not selected else 'pass'
    a=bounds(small,'fhalf_ci');b=bounds(large,'fhalf_ci')
    tolerance='pass' if a and b and b[0]>a[1] else 'fail' if a and b and b[1]<=a[0] else 'inconclusive'
    return combine([gain,noreverse,tolerance])


def h2(rows,inclusive):
    rows=sorted(rows,key=lambda r:float(r['f']));selected=[r for r in rows if comparison(r) and allowed(r,inclusive)]
    if any(bounds(r,'tau_ratio_ci') and bounds(r,'tau_ratio_ci')[1]<1 for r in selected):return 'pass'
    for later in selected:
        ci=bounds(later,'tau_ratio_ci')
        if not ci or not ci[0]<=1<=ci[1]:continue
        for earlier in selected:
            if float(earlier['f'])>=float(later['f']):continue
            span=[r for r in rows if float(earlier['f'])<=float(r['f'])<=float(later['f'])]
            if not all(comparison(r) and allowed(r,inclusive) and (numeric(r.get('large_multiple_front_fraction')) or 0)>.5 for r in span):continue
            p=numeric(earlier.get('tau_ratio'));q=numeric(later.get('tau_ratio'));a=bounds(earlier,'tau_ratio');b=bounds(later,'tau_ratio')
            if p is not None and q is not None and q>0 and p>=2*q:return 'pass'
            if a and b and math.isfinite(b[1]) and b[1]>0 and a[0]>=2*b[1]:return 'pass'
    return 'inconclusive'


def predictor(rows,inclusive):
    statuses=[]
    for feature in ('multiple_front_fraction','region_count_mean'):
        selected=[r for r in rows if r['statistic']==feature and r['outcome']=='direct_tau' and r['imaging_view']=='primary' and flag(r.get('forced_inclusive')) is inclusive and flag(r.get('primary_H2_feature')) is True and (inclusive or flag(r.get('primary_view')) is True)]
        if len(selected)!=1:statuses.append('inconclusive');continue
        r=selected[0];rho=numeric(r.get('pooled_spearman'));n=int(r['successful_folds']);missing=int(r['unavailable_folds'])
        if int(r['planned_folds'])!=4:raise AssertionError('A6 fold denominator changed')
        statuses.append('pass' if rho is not None and rho<=-.6 and n>=3 else 'fail' if rho is not None and (rho>-.6 or n+missing<3) else 'inconclusive')
    return 'pass' if 'pass' in statuses else 'fail' if all(s=='fail' for s in statuses) else 'inconclusive'


def strong(strong_row,weak,inclusive):
    if weak is None or not eligible(strong_row) or not eligible(weak) or not allowed(strong_row,inclusive) or not allowed(weak,inclusive):return 'inconclusive'
    a=bounds(strong_row,'tau');b=bounds(weak,'tau')
    return 'pass' if a[1]<b[0]/2 else 'fail' if math.isfinite(b[1]) and a[0]>=b[1]/2 else 'inconclusive'


def corroborated(rows,mode,inclusive):
    observations={float(r['sigma']):r for r in rows}
    if len(rows)!=len(observations):raise AssertionError('Duplicate sigma observations')
    baseline=observations.get(1.)
    if mode not in ('multiple_fronts','low_contrast') or baseline is None or baseline['regime']!=mode or not allowed(baseline,inclusive):return 'inconclusive'
    all_present=all(s in observations for s in (.25,.5,2.))
    matched=any(s in observations and observations[s]['regime']==mode and allowed(observations[s],inclusive) for s in (.25,.5,2.))
    qualified=all(s in observations and allowed(observations[s],inclusive) for s in (.25,.5,2.))
    return 'pass' if all_present and matched else 'fail' if qualified else 'inconclusive'


def independent_effect(small,large):
    rng=np.random.default_rng(2026102027)
    left=rng.multinomial(len(small),np.ones(len(small))/len(small),size=2000)@small/len(small)
    right=rng.multinomial(len(large),np.ones(len(large))/len(large),size=2000)@large/len(large)
    result={}
    for j,name in enumerate(('region_count','normalized_contrast')):
        low=np.zeros(2000);high=np.full(2000,math.inf);valid=left[:,j]>0
        low[valid]=high[valid]=right[valid,j]/left[valid,j]
        result[name]=(float(np.quantile(low,.025,method='inverted_cdf')),float(np.quantile(high,.975,method='inverted_cdf')))
    return result


def window(first,last):
    if not first or not last or not first.get('regime') or not last.get('regime'):return 'inconclusive'
    change=first['regime']!=last['regime']
    if not eligible(first) or not eligible(last):return 'fail' if change else 'inconclusive'
    a=bounds(first,'tau_ci');b=bounds(last,'tau_ci')
    if not a or not b:return 'fail' if change else 'inconclusive'
    return 'fail' if change or a[1]<b[0] or b[1]<a[0] else 'pass'


def close(actual,expected,label):
    if not np.allclose(actual,expected,rtol=2e-12,atol=2e-12,equal_nan=True):raise AssertionError(label)


def fhalf(f,lo,hi):
    """Direct first-crossing bounds, including a censored weak baseline."""
    if hi[0]!=lo[0] or not math.isfinite(hi[0]):
        below=[i for i in range(1,len(f)) if hi[i]<lo[0]/2]
        return float(f[0]),float(f[below[0]]) if below else math.inf
    target=lo[0]/2
    for i in range(1,len(f)):
        if hi[i]<=target:
            if hi[i]!=lo[i] or hi[i-1]!=lo[i-1]:return float(f[i-1]),float(f[i])
            weight=(lo[i-1]-target)/(lo[i-1]-lo[i]);point=math.exp(math.log(f[i-1])+weight*math.log(f[i]/f[i-1]))
            return point,point
        if lo[i]<=target<hi[i]:return float(f[i-1]),math.inf
    return float(f[-1]),math.inf


def bootstrap_archive_index(spec):
    index={}
    for stage in ('E2','E3','E4'):
        root=Path(spec['stages'][stage]['root'])
        for arm in read(root/'results/arms.json')['arms']:
            if float(arm['sigma'])!=1.:continue
            child=Path(arm['path']);manifest=read(child/'analysis_manifest.json')
            for group in manifest['groups']:
                worker=child/f"analysis_workers/{int(group['task']):04d}";result=read(worker/'result.json');path=worker/'arrays.npz'
                if digest(path)!=result['arrays_sha256']:raise AssertionError('Audited seed-bootstrap archive changed')
                index[(stage,group['model'],float(group['C_label']),float(group['f']))]=path
    return index


def native_archive(row,archive_index,cache):
    arm=row.get('selected_observation_arm');origin=row.get('common_observation_arm',arm)
    stage='E4' if arm=='preselected_E4' else 'E3' if origin=='common_new_E3_size' else 'E2'
    key=(stage,row['model'],float(row['C_label']),float(row['f']));view=row['view'];cachekey=key+(view,)
    if cachekey not in cache:
        with np.load(archive_index[key],allow_pickle=False) as z:
            cache[cachekey]=(z[view+'__tau_boot_lower'],z[view+'__tau_boot_upper'])
    return cache[cachekey]


def audit(bundle):
    bundle=Path(bundle);out=bundle/'results/final_analysis';complete=read(out/'complete.json')
    if not complete['passed'] or not complete['complete']:raise AssertionError('Final producer is incomplete')
    for name,expected in complete['outputs'].items():
        if digest(out/name)!=expected:raise AssertionError('Final output changed: '+name)
    for path,expected in complete['source_sha256'].items():
        if digest(path)!=expected:raise AssertionError('Frozen producer changed')
    if complete['source_sha256'].get(str(Path(__file__).resolve()))!=digest(__file__):raise AssertionError('Use the actual frozen independent audit executable')
    spec=read(complete['input_spec_path'])
    if digest(complete['input_spec_path'])!=complete['input_spec_sha256']:raise AssertionError('Final input specification changed')
    rows={name:table(out/(name+'.csv')) for name in ('criteria','condition_metrics','size_ratios','noise_tolerance','strong_input_candidates','fingerprint_seed_summaries','fingerprint_effects','E1_condition_metrics','window_qualification','horizon_qualification','original_seed_sensitivity','original_seed_condition_metrics')}
    index={(r['arm'],r['view'],flag(r['forced_inclusive']),r['model'],r['criterion']):r for r in rows['criteria']}
    if len(index)!=308:raise AssertionError('All four models, seven views, both arms and both forced views must remain visible')
    metrics={(r['arm'],r['model'],float(r['C_label']),float(r['f']),r['view']):r for r in rows['condition_metrics']}
    checked=0;features={arm:table(Path(spec['A6'][arm]['root'])/'H2_primary_features.csv') for arm in ('common','longest')}
    archives=bootstrap_archive_index(spec);archive_cache={}
    for ratio in rows['size_ratios']:
        if ratio.get('comparison_status')=='unavailable_unmatched_designed_input':
            if flag(ratio['evaluable']) is not False:raise AssertionError('An unmatched designed input was made evaluable')
            continue
        prefix=(ratio['arm'],ratio['model']);suffix=(float(ratio['f']),ratio['view'])
        small=metrics[prefix+(float(ratio['small_C_label']),)+suffix];large=metrics[prefix+(float(ratio['large_C_label']),)+suffix]
        if flag(ratio['evaluable']) is not (eligible(small) and eligible(large)):raise AssertionError('Size-comparison coverage qualification disagrees')
        sl,su=native_archive(small,archives,archive_cache);ll,lu=native_archive(large,archives,archive_cache)
        with np.errstate(divide='ignore',invalid='ignore'):
            lo=np.where(np.isinf(su),0.,ll/su);hi=np.where((sl<=0)|np.isinf(lu),math.inf,lu/sl)
        expected=(np.quantile(lo,.025,method='inverted_cdf'),np.quantile(hi,.975,method='inverted_cdf'))
        close(bounds(ratio,'tau_ratio_ci'),expected,'Independent size-ratio seed interval arithmetic');checked+=1
    for tolerance in rows['noise_tolerance']:
        if tolerance['fhalf_status'].startswith('unavailable_'):continue
        key=(tolerance['arm'],tolerance['model'],float(tolerance['C_label']));view=tolerance['view'];inclusive=flag(tolerance['forced_inclusive'])
        selected=sorted([r for k,r in metrics.items() if k[:3]==key and k[4]==view and allowed(r,inclusive)],key=lambda r:float(r['f']))
        if not selected or float(selected[0]['f'])!=.0001:raise AssertionError('Tolerance changed its actual weakest-input baseline')
        strengths=[float(r['f']) for r in selected];point=[bounds(r,'tau') if eligible(r) else (0.,math.inf) for r in selected]
        point_bounds=fhalf(strengths,[p[0] for p in point],[p[1] for p in point]);close(bounds(tolerance,'fhalf'),point_bounds,'Independent tolerance observation bounds')
        draws=[native_archive(r,archives,archive_cache) if eligible(r) else (np.zeros(2000),np.full(2000,math.inf)) for r in selected]
        results=[fhalf(strengths,[x[0][b] for x in draws],[x[1][b] for x in draws]) for b in range(2000)]
        expected=(np.quantile([r[0] for r in results],.025,method='inverted_cdf'),np.quantile([r[1] for r in results],.975,method='inverted_cdf'))
        close(bounds(tolerance,'fhalf_ci'),expected,'Independent tolerance whole-seed bounds');checked+=1
    for arm in ('common','longest'):
        for view in VIEWS:
            for inclusive in (False,True):
                h1s=[];h2s=[]
                for model in MODELS:
                    ratios=[r for r in rows['size_ratios'] if (r['arm'],r['model'],r['view'])==(arm,model,view)]
                    tolerances=[r for r in rows['noise_tolerance'] if (r['arm'],r['model'],r['view'],flag(r['forced_inclusive']))==(arm,model,view,inclusive)]
                    tolerances.sort(key=lambda r:float(r['C_label']))
                    one=h1(ratios,tolerances[0],tolerances[-1],inclusive);two=h2(ratios,inclusive)
                    if index[(arm,view,inclusive,model,'H1_model')]['status']!=one or index[(arm,view,inclusive,model,'H2_model')]['status']!=two:raise AssertionError('Independent H1/H2 model criterion disagrees')
                    h1s.append(one);h2s.append(two);checked+=2
                if index[(arm,view,inclusive,'all','H1_overall')]['status']!=count_status(h1s,3):raise AssertionError('H1 model denominator/threshold disagrees')
                if index[(arm,view,inclusive,'all','H2_overall')]['status']!=combine([count_status(h2s,2),predictor(features[arm],inclusive)]):raise AssertionError('H2 predictor/model gate disagrees')
                checked+=2
    # Re-evaluate all candidate strong-input and preset-sigma tests directly.
    candidate_groups={}
    for candidate in rows['strong_input_candidates']:
        arm=candidate['arm'];view=candidate['view'];model=candidate['model'];f=float(candidate['f']);inclusive=flag(candidate['forced_inclusive'])
        C=max(k[2] for k in metrics if k[0]==arm and k[1]==model);row=metrics[(arm,model,C,f,view)];weak=metrics.get((arm,model,C,.0001,view))
        if candidate['mode']!=row['regime']:raise AssertionError('Candidate failure mode differs from its measured A4 majority')
        if json.loads(candidate['strong'])['status']!=strong(row,weak,inclusive):raise AssertionError('Conservative strong-input bound test disagrees')
        e1=[r for r in rows['E1_condition_metrics'] if (r['model'],float(r['C_label']),float(r['f']),r['view'])==(model,C,f,view)]
        if json.loads(candidate['E1'])['status']!=corroborated(e1,candidate['mode'],inclusive):raise AssertionError('Preset E1 existence rule disagrees')
        candidate_groups.setdefault((arm,view,inclusive),[]).append(candidate);checked+=2
    # Each fingerprint seed summary is checked against its native A5 frames.
    receiptmap={r['primary_sha256']:r['primary_path'] for r in complete['A5_receipts']};framecache={};summaries={}
    for row in rows['fingerprint_seed_summaries']:
        digest_value=row['A5_primary_sha256'];path=receiptmap[digest_value]
        if digest_value not in framecache:
            if digest(path)!=digest_value:raise AssertionError('A5 observed frames changed')
            with gzip.open(path,'rt') as stream:framecache[digest_value]=json.load(stream)['result']['frames']
        frames=framecache[digest_value];start=float(row['physical_window_start']);stop=float(row['physical_window_stop']);basal=float(row['optical_basal_contrast'])
        counts=contrasts=duration=0.
        for left,right in zip(frames[:-1],frames[1:]):
            weight=max(0.,min(float(right['time']),stop)-max(float(left['time']),start))
            duration+=weight;counts+=weight*left['region_count'];contrasts+=weight*left['contrast']/basal
        close(duration,stop-start,'Fingerprint window coverage')
        expected=np.array([counts/duration,contrasts/duration]);close(expected,[float(row['observed_region_count_mean']),float(row['basal_normalized_contrast_mean'])],'Independent whole-seed fingerprint summary')
        key=(row['arm'],row['model'],float(row['C_label']),float(row['f']),row['view'])
        group=summaries.setdefault(key,{})
        if row['seed'] in group:raise AssertionError('A frame or repeated seed was counted as an independent fingerprint realization')
        group[row['seed']]=expected;checked+=1
    effects={}
    for row in rows['fingerprint_effects']:
        key=(row['arm'],row['model'],float(row['f']),row['view']);Cs=sorted(k[2] for k in summaries if (k[0],k[1],k[3],k[4])==key)
        if key not in effects:effects[key]=independent_effect(np.array(list(summaries[(key[0],key[1],Cs[0],key[2],key[3])].values())),np.array(list(summaries[(key[0],key[1],Cs[-1],key[2],key[3])].values())))
        close(bounds(row,'ci'),effects[key][row['statistic']],'Independent seed-bootstrap fingerprint interval');checked+=1
    for (arm,view,inclusive),candidates in candidate_groups.items():
        support=False
        for a in candidates:
            fingerprints=json.loads(a['fingerprints'])
            for name,measurement in fingerprints.items():
                expected=effects[(arm,a['model'],float(a['f']),view)][name]
                close(bounds(measurement,'ci'),expected,'Candidate fingerprint differs from independently verified seed effect')
            for b in candidates:
                if a['model']==b['model'] or a['mode']==b['mode'] or a['mode'] not in ('multiple_fronts','low_contrast') or b['mode'] not in ('multiple_fronts','low_contrast'):continue
                if any(json.loads(r[name])['status']!='pass' for r in (a,b) for name in ('strong','E1')):continue
                af=json.loads(a['fingerprints']);bf=json.loads(b['fingerprints'])
                separated=True
                for name in ('region_count','normalized_contrast'):
                    if name not in af or name not in bf or not af[name]['available'] or not bf[name]['available']:separated=False;continue
                    ai=bounds(af[name],'ci');bi=bounds(bf[name],'ci')
                    separated=separated and bool(ai and bi and (ai[1]<bi[0] or bi[1]<ai[0]))
                support=support or separated
        if index[(arm,view,inclusive,'all','H3_overall')]['status']!=('pass' if support else 'inconclusive'):raise AssertionError('Distinct-model/mode H3 fingerprint gate disagrees')
        checked+=1
    for row in rows['window_qualification']:
        key=(row['arm'],row['model'],float(row['C_label']),float(row['f']))
        pairs=[('analysis_first_half','analysis_second_half'),('main','late_quarter')]+[('main',v) for v in ('cadence_x2','cadence_x4','cadence_x8')]
        if row['status']!=combine([window(metrics.get(key+(a,)),metrics.get(key+(b,))) for a,b in pairs]):raise AssertionError('Independent A7 condition qualification disagrees')
        checked+=1
    for row in rows['horizon_qualification']:
        key=(row['model'],float(row['C_label']),float(row['f']),row['view'])
        if row['status']!=window(metrics.get(('common',)+key),metrics.get(('longest',)+key)):raise AssertionError('Independent horizon qualification disagrees')
        checked+=1
    originals={(r['model'],float(r['C_label']),float(r['f']),r['view']):r for r in rows['original_seed_condition_metrics']}
    for row in rows['original_seed_sensitivity']:
        key=(row['model'],float(row['C_label']),float(row['f']),row['view'])
        if row['status']!=window(originals[key],metrics[('longest',)+key]):raise AssertionError('Independent original-cohort qualification disagrees')
        if int(row['original_n_seeds']) not in (2,12) or int(row['all_long_n_seeds'])!=12:raise AssertionError('Original versus all12 cohort labels disagree')
        checked+=1
    receipt=dict(passed=True,complete=True,generated_utc=datetime.now(timezone.utc).isoformat(),checks=checked,
        audit_source_sha256=digest(__file__),producer_complete_sha256=digest(out/'complete.json'),
        all_four_models_and_seven_views_and_forced_views=True,independent_criterion_arithmetic=True,
        all_seed_summaries_replayed_against_native_A5_frames=True,whole_seed_fingerprint_bootstrap_replayed=True,
        final_outputs_sha256=complete['outputs'],scope='Independent final criterion and fingerprint joins; upstream native estimands and A6 regressions retain their separate declared independent gates.')
    path=out/'independent_final_decision_audit.json'
    if path.exists():raise FileExistsError('Preserve the existing independent final audit')
    path.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    return dict(passed=True,checks=checked,path=str(path),sha256=digest(path))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--bundle',type=Path,required=True);args=parser.parse_args()
    print(json.dumps(audit(args.bundle.resolve()),indent=2))
