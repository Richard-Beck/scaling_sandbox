"""Portable scientific analysis and simulated point nulls for Figure 4 round 3.

Round-2 fallback is explicit; it never pretends to be the new 96-seed campaign.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path
import resource
import time
for _thread_variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_thread_variable]='1'
import numpy as np
import figure4_round3_estimators_20261022 as E

ROOT=Path(__file__).resolve().parents[2]
DEFAULT_OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4'
ROUND2=ROOT/'figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4'
MODELS=('wave_pinning','holmes_model3','spring_ring','goryachev')
NULLS=('drift','same_rate','perimeter')


def read_csv(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalized_P(archive):
    """Explicit round-3 definition; historical stored P used another denominator."""
    moment=archive['M1_real'].astype(float)+1j*archive['M1_imag'].astype(float)
    M0=archive['M0'].astype(float)
    return np.divide(moment,M0,out=np.full_like(moment,np.nan+1j*np.nan),where=M0>0)


def native_regime(regions,contrast_ratio):
    regions=np.asarray(regions);low=np.asarray(contrast_ratio)<.1
    fractions=dict(low_contrast=float(np.mean(low)),multiple_fronts=float(np.mean((~low)&(regions>=2))),single_front=float(np.mean((~low)&(regions==1))))
    return next((name for name,fraction in fractions.items() if fraction>.5),'mixed')


def control_flag(row,noisy_regime,lower_rows=()):
    """Logical round-2 union; unavailable low-sigma evidence stays explicit."""
    with np.load(row['compact'],allow_pickle=False) as z:
        mean_regime=native_regime(z['region_count'],z['contrast_ratio'])
        mean_fraction=float(np.mean((z['R']<.1)|(z['contrast_ratio']<.1)))
    mean_forced=mean_fraction>.5
    if mean_forced:return True,dict(forced_mean=True,forced_E1='unknown',mean_regime=mean_regime,noise_regime=noisy_regime,forcing_reason='matched_mean_union_fraction_above_half')
    if noisy_regime not in ['multiple_fronts','low_contrast'] or noisy_regime!=mean_regime:
        return False,dict(forced_mean=False,forced_E1=False,mean_regime=mean_regime,noise_regime=noisy_regime,forcing_reason='abnormal_same_majority_gate_false_before_lower_sigma')
    lower_labels={}
    for sigma in [.25,.5]:
        members=[r for r in lower_rows if float(r['sigma'])==sigma]
        if len(members)!=12 or {int(r['replicate']) for r in members}!=set(range(12)):break
        regions=[];contrast=[]
        for member in members:
            if member.get('compact_sha256') and digest(member['compact'])!=member['compact_sha256']:raise ValueError('Lower-sigma compact hash changed')
            with np.load(member['compact'],allow_pickle=False) as z:regions.extend(z['region_count']);contrast.extend(z['contrast_ratio'])
        lower_labels[sigma]=native_regime(regions,contrast)
    if len(lower_labels)==2:
        flag=all(label==noisy_regime for label in lower_labels.values())
        return flag,dict(forced_mean=False,forced_E1=flag,mean_regime=mean_regime,noise_regime=noisy_regime,forcing_reason='paired12_seed_sigma0.25_0.5_observed_rule',lower_sigma_regimes=lower_labels)
    return 'unknown',dict(forced_mean=False,forced_E1='unknown',mean_regime=mean_regime,noise_regime=noisy_regime,forcing_reason='same_abnormal_mean_and_sigma1_requires_sigma0.25_and_0.5')


def forcing_candidates(campaign,output):
    campaign=Path(campaign);control_path=campaign/'matched_mean_controls/results/compact_control_index.csv'
    controls=read_csv(control_path) if control_path.exists() else [];index=read_csv(campaign/'results/compact_index.csv');per={}
    for row in index:per.setdefault(row['condition_id'],[]).append(row)
    candidates=[]
    for row in controls:
        members=per.get(row['condition_id'],[])
        if len(members)<96:raise ValueError('Forcing diagnostics require the completed sigma1 cohort')
        regions=[];contrast=[]
        for member in members:
            with np.load(member['compact'],allow_pickle=False) as z:regions.extend(z['region_count']);contrast.extend(z['contrast_ratio'])
        label=native_regime(regions,contrast);flag,details=control_flag(row,label)
        if flag=='unknown':candidates.append(dict(condition_id=row['condition_id'],model=row['model'],size=row['size'],fraction=row['fraction'],**details))
    if candidates:E.write_csv(output,candidates)
    else:
        Path(output).parent.mkdir(parents=True,exist_ok=True);Path(output).write_text('condition_id,model,size,fraction,forcing_reason\n')
    return candidates


def tracked_uncertainty(seeds,pooled,summary):
    """Whole seeds, fixed observed fit cutoff, paired with the kobs bootstrap."""
    n=len(seeds);weights=np.random.default_rng(summary['bootstrap_seed']).multinomial(n,np.full(n,1/n),size=E.DRAW_COUNT)
    nums=[];counts=[];times=[];njumps=[];jump_loss=[];unscored=[]
    maxlag=pooled.get('MSD_max_lag',0)
    for seed in seeds:
        measure,msd=E.drift_jump(seed['time'],seed['C'],seed['ids'],seed['headings'],P=seed['P'],max_lag=maxlag)
        if not msd:
            m=int(round(maxlag/seed['cadence']));nums.append(np.zeros(m));counts.append(np.zeros(m))
        else:
            counts.append(np.array([r['eligible_pairs'] for r in msd],float));nums.append(np.nan_to_num([r['MSD']*r['eligible_pairs'] for r in msd]))
        times.append(measure['valid_time']);njumps.append(measure['n_switches']);angles=measure.get('jump_angles',[])
        jump_loss.append(float(np.sum(1-np.cos(angles))));unscored.append(measure.get('n_scored_switches',0)!=measure['n_switches'])
    lag=np.arange(1,len(nums[0])+1)*seeds[0]['cadence']
    with np.errstate(invalid='ignore',divide='ignore'):
        msdb=(weights@np.asarray(nums))/(weights@np.asarray(counts))
        finite=np.isfinite(msdb);Ds=np.nansum(msdb*lag,axis=1)/(2*np.sum(finite*lag**2,axis=1))
        Ds[finite.sum(1)<2]=np.nan
        exposure=weights@np.asarray(times);lamb=weights@np.asarray(njumps)/exposure
        eff=weights@np.asarray(jump_loss)/exposure
    eff[(weights@np.asarray(unscored,int))>0]=np.nan
    kp=(2*np.pi/summary['circumference_boot'])**2*Ds+eff;j=E.analysis_lag(summary)
    ratio=kp/summary['kboot'][:,j] if j is not None else np.full(E.DRAW_COUNT,np.nan)
    out=dict(tracked_bootstrap_resamples=E.DRAW_COUNT,tracked_bootstrap_fit_range='fixed_original_pooled_median_tenure_cutoff',
             tracked_bootstrap_valid_fraction=float(np.isfinite(ratio).mean()))
    for name,value in [('Ds',Ds),('lambda_s',lamb),('lambda_eff',eff),('kpred',kp),('kpred_over_kobs',ratio)]:
        lo,hi=E.percentile(value,True);out[name+'_ci_lower']=lo;out[name+'_ci_upper']=hi
    return out


def export_analysis(groups, contexts, data, prefix=''):
    """groups indexed (model,size,f); every missing designed cell gets a row."""
    lagrows=[];alphas=[];analysis=[];endpoints=[]
    fs_by_model={m:sorted({k[2] for k in groups if k[0]==m}) for m in set(k[0] for k in groups)}
    for key,s in groups.items():
        m,c,f=key;context=contexts.get(key,{})
        base=dict(model=m,C_label=c,f=f,C_mean=s['circumference'],n_seeds=s['n_seeds'],**context)
        for j,lag in enumerate(s['lags']):
            lo,hi=E.percentile(s['Cboot'][:,j]);klo,khi=E.percentile(s['kboot'][:,j],True)
            lagrows.append(dict(base,lag=lag,target_lag=s['target_lags'][j],C=s['C'][j],C_ci_lower=lo,C_ci_upper=hi,
                                k=s['k'][j],k_ci_lower=klo,k_ci_upper=khi,coverage=s['coverage'][j],
                                eligible_pairs=s['eligible_pairs'][j],possible_pairs=s['possible_pairs'][j],
                                bootstrap_valid_fraction=float(np.isfinite(s['kboot'][:,j]).mean())))
        j=E.analysis_lag(s)
        if j is None:analysis.append(dict(base,status='unavailable_no_analysis_lag',analysis_lag=math.nan,k=math.nan))
        else:
            lo,hi=E.percentile(s['kboot'][:,j],True)
            with np.errstate(invalid='ignore',divide='ignore'):SE=float(np.nanstd(np.log(s['kboot'][:,j]),ddof=1))
            analysis.append(dict(base,status='resolved',analysis_lag=s['lags'][j],k=s['k'][j],k_ci_lower=lo,k_ci_upper=hi,
                                 SE_log_k=SE,bootstrap_valid_fraction=float(np.isfinite(s['kboot'][:,j]).mean())))
    for m,fs in fs_by_model.items():
        sizes=E.SPRING_SIZES if m=='spring_ring' else E.CHEMICAL_SIZES
        for f in fs:
            neighborhood={c:groups[(m,c,f)] for c in sizes if (m,c,f) in groups}
            for c in sizes:
                ctx=dict(contexts.get((m,c,f),{}));flags=[contexts.get((m,n,f),{}).get('forced_combined') for n in E.fixed_neighborhood(c,m)]
                ctx['forced_combined']=any(x is True or x=='True' for x in flags) if all(x in [True,False,'True','False'] for x in flags) else True if any(x is True or x=='True' for x in flags) else 'unknown'
                alphas.append(dict(model=m,f=f,**ctx,**E.local_alpha(neighborhood,c,m)))
            if sizes[0] not in neighborhood or sizes[-1] not in neighborhood:
                endpoints.append(dict(model=m,f=f,status='unavailable_missing_endpoint'));continue
            small,large=neighborhood[sizes[0]],neighborhood[sizes[-1]]
            candidates=[]
            for lag in set(small['lags'])&set(large['lags']):
                a=int(np.flatnonzero(small['lags']==lag)[0]);b=int(np.flatnonzero(large['lags']==lag)[0])
                cs=np.array([small['C'][a],large['C'][b]])
                if np.all((cs>.03)&(cs<.995)) and min(small['coverage'][a],large['coverage'][b])>=.5:
                    score=np.max(abs(np.log(-np.log(cs))-math.log(math.log(2))));candidates.append((score,lag,a,b))
            if not candidates:endpoints.append(dict(model=m,f=f,status='unavailable_no_shared_endpoint_lag'));continue
            _,lag,a,b=min(candidates)
            with np.errstate(invalid='ignore',divide='ignore'):
                alpha=np.log(small['k'][a]/large['k'][b])/np.log(large['circumference']/small['circumference'])
                boot=np.log(small['kboot'][:,a]/large['kboot'][:,b])/np.log(large['circumference_boot']/small['circumference_boot'])
            lo,hi=E.percentile(boot,True)
            endpoints.append(dict(model=m,f=f,lag=lag,alpha=alpha,alpha_ci_lower=lo,alpha_ci_upper=hi,
                                  small_C_label=sizes[0],large_C_label=sizes[-1],bootstrap_valid_fraction=float(np.isfinite(boot).mean()),status='resolved'))
    tols=[];folds=[]
    for m,fs in fs_by_model.items():
        sizes=E.SPRING_SIZES if m=='spring_ring' else E.CHEMICAL_SIZES
        for ref in ((100.,200.,400.) if m=='spring_ring' else (1000.,2000.,4000.)):
            per={}
            for c in sizes:
                cs=[];boots=[];means=[]
                for f in fs:
                    s=groups.get((m,c,f))
                    if s is None or not np.any(np.isclose(s['lags'],ref)):
                        cs.append(math.nan);boots.append(np.full(E.DRAW_COUNT,np.nan));continue
                    j=int(np.flatnonzero(np.isclose(s['lags'],ref))[0])
                    cs.append(s['C'][j]);boots.append(s['Cboot'][:,j]);means.append(s['circumference'])
                value=E.tolerance(fs,cs,boots);per[c]=value
                row={k:v for k,v in value.items() if not k.startswith('bootstrap_') or k=='bootstrap_valid_fraction'}
                tols.append(dict(model=m,C_label=c,C_mean=float(np.mean(means)) if means else c,reference_lag=ref,
                                 forced_inclusive=True,bootstrap_resamples=E.DRAW_COUNT,designed_inputs=fs,**row))
            for c in sizes:
                fold=dict(fhalf_fold=1.,fhalf_fold_ci_lower=1.,fhalf_fold_ci_upper=1.,status='resolved',bootstrap_valid_fraction=1.) if c==sizes[0] else E.tolerance_fold(per[sizes[0]],per[c])
                if c==sizes[0] and not np.isfinite(per[c]['fhalf']):fold.update(fhalf_fold=math.nan,status='unavailable_reference')
                folds.append(dict(model=m,C_label=c,C_equivalent=c*20 if m=='spring_ring' else c,reference_lag=ref,**fold))
    validrows=[dict(model=k[0],C_label=k[1],f=k[2],n_seeds=s['n_seeds'],valid_fraction=contexts.get(k,{}).get('valid_fraction',1. if prefix.startswith('null') else math.nan),
                    source=contexts.get(k,{}).get('source',''),coverage_min=float(np.nanmin(s['coverage'])),coverage_max=float(np.nanmax(s['coverage']))) for k,s in groups.items()]
    for name,rows in [('condition_lag_metrics',lagrows),('local_alpha',alphas),('analysis_lag_metrics',analysis),('end_to_end_alpha',endpoints),('noise_tolerance',tols),('noise_tolerance_folds',folds),('valid_fraction',validrows)]:E.write_csv(data/(prefix+name+'.csv'),rows)
    return dict(alpha_rows=len(alphas),alpha_resolved=sum(r['status']=='resolved' for r in alphas),alpha_failures=[r for r in alphas if r['status']!='resolved'],
                conditions=len(groups),condition_lag_rows=len(lagrows),bootstrap_resamples=E.DRAW_COUNT,
                analysis_precision_unmet=[r for r in analysis if r.get('SE_log_k',math.inf)>.07])


def compact(out,campaign,tracked=True):
    campaign=Path(campaign);manifest=json.loads((campaign/'manifest.json').read_text());manifest_sha=digest(campaign/'manifest.json')
    status_path=campaign/'status.json';status=json.loads(status_path.read_text()) if status_path.exists() else {}
    if not status.get('complete') or not status.get('precision_qualified') or status.get('status')!='COMPLETE':
        raise ValueError('Compact campaign final collection and precision qualification must be COMPLETE before replacing provisional figures')
    allocation=json.loads((campaign/'results/allocation.json').read_text());targets=allocation['targets']
    precision_path=Path(allocation.get('precision_source',''))
    if not precision_path.is_file() or digest(precision_path)!=allocation.get('precision_sha256'):raise ValueError('Final whole-seed precision receipt is missing or changed')
    precision_rows=json.loads(precision_path.read_text()).get('conditions',[]);precision={r['condition_id']:r for r in precision_rows}
    if len(precision)!=len(manifest['conditions']):raise ValueError('Final precision receipt omits or duplicates designed conditions')
    for condition in manifest['conditions']:
        cid=condition['condition_id'];assessment=precision.get(cid,{})
        if assessment.get('n_seeds')!=targets[cid]:raise ValueError('Final precision cohort differs from the allocated independent seeds: '+cid)
        unresolved=(assessment.get('se_log_k') is not None and assessment['se_log_k']>.07) or assessment.get('reason')=='bootstrap_valid_fraction_below_0.975'
        if unresolved and targets[cid]<384:raise ValueError('Final cohort remains below the requested precision before384seed cap: '+cid)
    records=read_csv(campaign/'results/compact_index.csv');per={}
    for row in records:
        if row['status'].lower()!='complete':continue
        if row['manifest_sha256']!=manifest_sha:raise ValueError('Compact index belongs to a different campaign manifest')
        per.setdefault(row['condition_id'],[]).append(row)
    missing=[c['condition_id'] for c in manifest['conditions'] if len(per.get(c['condition_id'],[]))<96]
    if missing:raise ValueError(f'Complete compact campaign required before replacing provisional data: {len(missing)} conditions below96 seeds')
    for condition in manifest['conditions']:
        cid=condition['condition_id'];members=per[cid];target=targets[cid]
        if len(members)!=target or {int(r['replicate']) for r in members}!=set(range(target)):
            raise ValueError('Compact index does not match the final allocated independent seed cohort: '+cid)
    old=[r for r in read_csv(ROUND2/'supporting/data/old_tau_condition_metrics.csv') if r['arm']=='common' and r['view']=='main']
    flags={(r['model'],float(r['C_label']),float(r['f'])):r['forced_combined'] for r in old}
    control_path=campaign/'matched_mean_controls/results/compact_control_index.csv'
    controls={r['condition_id']:r for r in read_csv(control_path)} if control_path.exists() else {}
    sigma_path=campaign/'lower_sigma_controls/results/compact_sigma_index.csv';lower_rows=read_csv(sigma_path) if sigma_path.exists() else []
    groups={};first12={};contexts={};trackgroups={};extended=[];forcing_rows=[]
    for condition in manifest['conditions']:
        key=(condition['model'],float(condition['size']),float(condition['fraction']));rows=sorted(per[condition['condition_id']],key=lambda r:int(r['replicate']))
        moments=[];means=[];valid_n=0;frames_n=0;multiple_n=0;low_n=0;single_n=0;multiple_seed=[];frame_seed=[]
        for row in rows:
            path=Path(row['compact'])
            if row.get('compact_sha256') and digest(path)!=row['compact_sha256']:raise ValueError('Compact content hash changed: '+str(path))
            with np.load(path,allow_pickle=False) as z:
                t=z['time'].astype(float);P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float);C=z['circumference'].astype(float)
                cadence=float(np.median(np.diff(t)));target,lags=E.lag_grid(key[0],t[-1]-t[0],cadence)
                moments.append(E.seed_moments(P,abs(P)>=.1,cadence,lags));means.append(float(C.mean()));valid_n+=int(np.sum(abs(P)>=.1));frames_n+=len(P)
                mult=int(np.sum((z['region_count']>=2)&(z['contrast_ratio']>=.1)));multiple_n+=mult;multiple_seed.append(mult);frame_seed.append(len(P))
                low_n+=int(np.sum(z['contrast_ratio']<.1));single_n+=int(np.sum((z['region_count']==1)&(z['contrast_ratio']>=.1)))
        groups[key]=E.summarize(moments,means,('compact',key),target,lags)
        first12[key]=E.summarize(moments[:12],means[:12],('compact_first12',key),target,lags)
        weights=np.random.default_rng(groups[key]['bootstrap_seed']).multinomial(len(rows),np.full(len(rows),1/len(rows)),size=E.DRAW_COUNT)
        mlo,mhi=E.percentile((weights@np.asarray(multiple_seed))/(weights@np.asarray(frame_seed)))
        noisy_regime='low_contrast' if low_n/frames_n>.5 else 'multiple_fronts' if multiple_n/frames_n>.5 else 'single_front' if single_n/frames_n>.5 else 'mixed'
        flag=flags.get(key,'unknown');details={}
        if key not in flags and condition['condition_id'] in controls:flag,details=control_flag(controls[condition['condition_id']],noisy_regime,[r for r in lower_rows if r['condition_id']==condition['condition_id']])
        forcing_rows.append(dict(condition_id=condition['condition_id'],model=key[0],C_label=key[1],f=key[2],forced_combined=flag,**details))
        contexts[key]=dict(arm='common',view='main',source='completed_compact_campaign',forced_combined=flag,valid_fraction=valid_n/frames_n,multiple_front_fraction=multiple_n/frames_n,
                           multiple_front_fraction_ci_lower=mlo,multiple_front_fraction_ci_upper=mhi)
        if tracked:trackgroups[key]=rows
        max_horizon=max(float(r['horizon']) for r in rows);native=3200. if key[0]=='spring_ring' else 32000.
        if len(rows)>96 or max_horizon>native:
            extended.append(dict(model=key[0],C_label=key[1],f=key[2],n_seeds=len(rows),max_horizon=max_horizon,native_horizon=native,extended=max_horizon>native,over96=len(rows)>96))
        if len(groups)%20==0:print('compact conditions',len(groups),flush=True)
    data=out/'supporting/data';result=export_analysis(groups,contexts,data)
    export_analysis(first12,{k:dict(v,source='completed_compact_first12') for k,v in contexts.items()},data,prefix='round3_first12_')
    export_round2_comparison(data,first12,groups)
    E.write_csv(data/'extended_conditions.csv',extended if extended else [dict(status='none',n_seeds=96)])
    E.write_csv(data/'forcing_conditions.csv',forcing_rows)
    if tracked:result['tracked']=compact_tracked(out,groups,contexts,trackgroups)
    result.update(generated_utc=datetime.now(timezone.utc).isoformat(),source='completed_compact_campaign',completed_campaign=True,
                  campaign_manifest=str(campaign/'manifest.json'),campaign_manifest_sha256=manifest_sha,
                  compact_index_sha256=digest(campaign/'results/compact_index.csv'),seed_counts=sorted({s['n_seeds'] for s in groups.values()}),
                  unknown_forcing_conditions=[dict(model=k[0],C_label=k[1],f=k[2]) for k,v in contexts.items() if v['forced_combined']=='unknown'],
                  extended_or_over96=[dict(model=k[0],C_label=k[1],f=k[2],n_seeds=s['n_seeds']) for k,s in groups.items() if s['n_seeds']>96],
                  limitations=['Endpoint eligibility is abs(P)>=.1; tracked detector eligibility is distinct.','Forcing flags join unchanged round2 exact-condition control decisions; absent flags remain unknown.'])
    meta=out/'supporting/metadata';meta.mkdir(parents=True,exist_ok=True);(meta/'analysis_receipt.json').write_text(json.dumps(result,indent=2,allow_nan=True)+'\n')
    return result


def export_round2_comparison(data,first12,full):
    matched_path=data/'round2_matched_P_condition_lag_metrics.csv'
    matched=read_csv(matched_path) if matched_path.exists() else []
    original=read_csv(ROUND2/'supporting/data/condition_lag_metrics.csv')
    old_index={(r['model'],float(r['C_label']),float(r['f']),float(r['lag'])):r for r in original if r['arm']=='common' and r['view']=='main'}
    matched_index={(r['model'],float(r['C_label']),float(r['f']),float(r['lag'])):r for r in matched}
    rows=[]
    for key,s in first12.items():
        ref=200. if key[0]=='spring_ring' else 2000.;old=old_index.get((*key,ref));newP=matched_index.get((*key,ref))
        if old is None:continue
        j=int(np.flatnonzero(np.isclose(s['lags'],ref))[0]);lo,hi=E.percentile(s['kboot'][:,j],True)
        complete=full[key];k=int(np.flatnonzero(np.isclose(complete['lags'],ref))[0]);flo,fhi=E.percentile(complete['kboot'][:,k],True)
        comparison_values=np.array([float(old['k']),float(newP['k']) if newP else math.nan,complete['k'][k]])
        rows.append(dict(measure='k_at_fixed_reference',model=key[0],C_label=key[1],f=key[2],reference_lag=ref,n_old_seeds=int(old['n_seeds']),old_12seed_subset=int(old['n_seeds'])==12,n_new_seeds=complete['n_seeds'],n_new_first12_seeds=s['n_seeds'],
                         original_round2_value=old['k'],original_round2_ci_lower=old['k_ci_lower'],original_round2_ci_upper=old['k_ci_upper'],
                         matched_P_round2_value=newP['k'] if newP else math.nan,matched_P_round2_ci_lower=newP['k_ci_lower'] if newP else math.nan,matched_P_round2_ci_upper=newP['k_ci_upper'] if newP else math.nan,
                         round3_first12_value=s['k'][j],round3_first12_ci_lower=lo,round3_first12_ci_upper=hi,
                         full_round3_value=complete['k'][k],full_round3_ci_lower=flo,full_round3_ci_upper=fhi,
                         status='resolved' if np.all(np.isfinite(comparison_values)&(comparison_values>0)) else 'matched_P_baseline_unavailable' if newP is None else 'unavailable_nonpositive_or_undefined_rate'))
    original_tol=read_csv(ROUND2/'supporting/data/noise_tolerance.csv');matched_tol=read_csv(data/'round2_matched_P_noise_tolerance.csv') if (data/'round2_matched_P_noise_tolerance.csv').exists() else []
    orig_index={(r['model'],float(r['C_label']),float(r['reference_lag'])):r for r in original_tol if r['arm']=='common' and r['view']=='main' and r['forced_inclusive']=='True'}
    match_index={(r['model'],float(r['C_label']),float(r['reference_lag'])):r for r in matched_tol}
    full_tol={(r['model'],float(r['C_label']),float(r['reference_lag'])):r for r in read_csv(data/'noise_tolerance.csv')}
    for r in read_csv(data/'round3_first12_noise_tolerance.csv'):
        key=r['model'],float(r['C_label']),float(r['reference_lag']);old=orig_index.get(key);newP=match_index.get(key)
        if old is None:continue
        old_counts={str(f):int(record['n_seeds']) for (m,c,f,lag),record in old_index.items() if (m,c,lag)==key}
        new_counts={str(f):s['n_seeds'] for (m,c,f),s in full.items() if (m,c)==key[:2]};complete=full_tol.get(key,{})
        nold=next(iter(set(old_counts.values()))) if len(set(old_counts.values()))==1 else None
        nnew=next(iter(set(new_counts.values()))) if len(set(new_counts.values()))==1 else None
        rows.append(dict(measure='fhalf_first_downcrossing',model=key[0],C_label=key[1],reference_lag=key[2],n_old_seeds=nold,n_new_seeds=nnew,n_new_first12_seeds=12,
                         old_seed_counts_by_f=old_counts,new_seed_counts_by_f=new_counts,old_12seed_subset=all(n==12 for n in old_counts.values()),
                         original_round2_value=old['fhalf'],original_round2_ci_lower=old['fhalf_ci_lower'],original_round2_ci_upper=old['fhalf_ci_upper'],
                         matched_P_round2_value=newP['fhalf'] if newP else math.nan,matched_P_round2_ci_lower=newP['fhalf_ci_lower'] if newP else math.nan,matched_P_round2_ci_upper=newP['fhalf_ci_upper'] if newP else math.nan,
                         round3_first12_value=r['fhalf'],round3_first12_ci_lower=r['fhalf_ci_lower'],round3_first12_ci_upper=r['fhalf_ci_upper'],
                         full_round3_value=complete.get('fhalf',math.nan),full_round3_ci_lower=complete.get('fhalf_ci_lower',math.nan),full_round3_ci_upper=complete.get('fhalf_ci_upper',math.nan),
                         original_round2_status=old['status'],matched_P_round2_status=newP['status'] if newP else 'unavailable',full_round3_status=complete.get('status','unavailable'),
                         status='resolved' if newP and all(z.get('status')=='resolved' for z in [old,newP,complete]) else 'bounded_or_unavailable_fhalf'))
    E.write_csv(data/'round2_comparison.csv',rows)


def _compact_tracked_worker_limit(requested=None):
    """Stay within this SLURM task; ordinary interactive imports use one CPU."""
    allocated=1
    if os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOBID'):
        try:allocated=max(1,int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
        except ValueError:allocated=1
    return allocated if requested is None else max(1,min(allocated,int(requested)))


def _compact_tracked_worker_init():
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'


def _compact_seed0_switch_indices(seed):
    """Export adjacent accepted-dominant switches, with genuine ID gaps intact."""
    filtered=E.smooth_flicker(seed['ids'],5)
    return np.flatnonzero((filtered[1:]!=filtered[:-1])&(filtered[1:]>=0)&(filtered[:-1]>=0))+1


def _compact_tracked_condition(task):
    """Unchanged scientific operations for one condition, with local histories."""
    key,seed_rows,s,context=task
    rows=[];msdrows=[];mles=[];events=[]
    # Each worker loads one condition; full campaign histories stay off heap.
    seeds=[]
    for row in seed_rows:
        with np.load(row['compact'],allow_pickle=False) as z:
            t=z['time'].astype(float);C=z['circumference'].astype(float);P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float)
            offsets=z['region_offset'];region_ids=z['region_id'];head=z['region_heading']
            maps=[{int(region_ids[j]):float(head[j]) for j in range(int(a),int(b)) if np.isfinite(head[j])} for a,b in zip(offsets[:-1],offsets[1:])]
            seeds.append(dict(time=t,P=P,C=C,ids=z['dominant_track_id'].astype(int),headings=maps,cadence=float(np.median(np.diff(t)))))
    dt=seeds[0]['cadence'];ids=[];maps=[];C=[];P=[];seed_labels=[]
    for seed_index,seed in enumerate(seeds):
        if not np.isclose(seed['cadence'],dt):raise ValueError('Mixed compact cadence within a condition')
        ids.extend([-1]+list(seed['ids']));maps.extend([{}]+seed['headings']);C.extend([seed['C'][0]]+list(seed['C']));P.extend([complex(np.nan)]+list(seed['P']))
        seed_labels.extend([-1]+[seed_index]*len(seed['P']))
    P=np.asarray(P);C=np.asarray(C);t=np.arange(len(ids))*dt
    base=dict(model=key[0],C_label=key[1],f=key[2],C_mean=s['circumference'],n_seeds=len(seeds),**context)
    for hold,fraction in [(5,.25),(2,.25),(10,.25),(5,.125),(5,.5)]:
        measure,msd=E.drift_jump(t,C,ids,maps,hold=hold,range_fraction=fraction,P=P);measure.pop('jump_angles',None)
        j=E.analysis_lag(s);kobs=s['k'][j] if j is not None else math.nan;kp=(2*np.pi/base['C_mean'])**2*measure['Ds']+measure['lambda_eff']
        ci=tracked_uncertainty(seeds,measure,s) if hold==5 and fraction==.25 else {}
        rows.append(dict(base,**measure,**ci,kobs=kobs,kpred=kp,kpred_over_kobs=kp/kobs if kobs>0 else math.nan,analysis_lag=s['lags'][j] if j is not None else math.nan))
        if hold==5 and fraction==.25:msdrows.extend(dict(base,**m) for m in msd)
    for step in [1,2]:mles.append(dict(base,**E.wrapped_mixture_mle(P,base['C_mean'],dt,step,seed_labels)))
    seed=seeds[0]
    for index in _compact_seed0_switch_indices(seed):
        events.append(dict(model=key[0],C_label=key[1],f=key[2],replicate=0,time=seed['time'][index]))
    execution=dict(worker_pid=os.getpid(),peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rows,msdrows,mles,events,execution


def compact_tracked(out,groups,contexts,trackgroups,workers=None):
    rows=[];msdrows=[];mles=[];events=[];execution=[]
    worker_count=min(_compact_tracked_worker_limit(workers),max(1,len(trackgroups)))
    _compact_tracked_worker_init()
    tasks=((key,seed_rows,groups[key],contexts[key]) for key,seed_rows in trackgroups.items())
    def collect(results):
        for condition_rows,condition_msd,condition_mles,condition_events,record in results:
            rows.extend(condition_rows);msdrows.extend(condition_msd);mles.extend(condition_mles);events.extend(condition_events);execution.append(record)
    if worker_count==1:collect(map(_compact_tracked_condition,tasks))
    else:
        # Spawn workers inherit one BLAS thread before NumPy imports. map keeps
        # manifest condition order, preserving CSV bytes and all RNG identities.
        with ProcessPoolExecutor(max_workers=worker_count,mp_context=get_context('spawn'),initializer=_compact_tracked_worker_init) as pool:
            collect(pool.map(_compact_tracked_condition,tasks,chunksize=1))
    data=out/'supporting/data';E.write_csv(data/'drift_decomposition.csv',rows);E.write_csv(data/'within_tenure_msd.csv',msdrows);E.write_csv(data/'detector_free_mle.csv',mles);E.write_csv(data/'seed0_switch_events.csv',events)
    return dict(conditions=len(trackgroups),workers=worker_count,allocated_cpus=_compact_tracked_worker_limit(),BLAS_threads=1,
                worker_pids=sorted({r['worker_pid'] for r in execution}),maximum_worker_peak_resident_memory_KiB=max((r['peak_resident_memory_KiB'] for r in execution),default=0),
                jump_heading_source='whole_P',drift_coordinate='tracked-region physical centroid, local unwrap within tenure')


def fallback(out, tracked=True):
    import signature_figure4_e0_20261020 as OLD
    metrics=[r for r in read_csv(ROUND2/'supporting/data/old_tau_condition_metrics.csv') if r['arm']=='common' and r['view']=='main']
    receipt=json.loads((ROUND2/'supporting/metadata/analysis_receipt.json').read_text())
    if digest(ROUND2/'supporting/data/old_tau_condition_metrics.csv')!=receipt['outputs']['old_tau_condition_metrics.csv']:raise ValueError('Frozen round2 scalar-selection table changed')
    manifests={}
    for row in receipt['manifests']:
        if digest(row['path'])!=row['sha256']:raise ValueError('Frozen round2 analysis manifest changed')
        manifests[str(Path(row['path']).resolve())]=json.loads(Path(row['path']).read_text())
    groups={};contexts={};case_groups={};sources=[]
    for row in metrics:
        key=(row['model'],float(row['C_label']),float(row['f']))
        # Round-2 receipt gives selected observation provenance, not discovery order.
        stage='E3' if (row.get('common_observation_arm') or row['selected_observation_arm'])=='common_new_E3_size' else 'E2'
        # Explicitly load the immutable sigma1 stage manifest; avoids horizon arm mixing.
        manifest_path=ROOT/f'data/runs/signature_figure4_{stage}_analysis_20261020/arms/sigma_1p0/analysis_manifest.json'
        selected=manifests[str(manifest_path.resolve())]
        g=next(g for g in selected['groups'] if (g['model'],float(g['C_label']),float(g['f']))==key)
        moments=[];means=[];valid_n=0;frames_n=0
        for case in g['noise']:
            p=Path(case['scalar_path'])
            if digest(p)!=case['scalar_sha256']:raise ValueError('Round2 scalar hash changed: '+str(p))
            with np.load(p,allow_pickle=False) as z:
                t=z['time'];ix=OLD.select(t,'main',case['native_horizon']);t=t[ix]
                P=normalized_P(z)[ix];C=np.broadcast_to(z['circumference'],z['time'].shape)[ix]
                cadence=float(np.median(np.diff(t)));target,lags=E.lag_grid(key[0],t[-1]-t[0]+cadence,cadence)
                moments.append(E.seed_moments(P,abs(P)>=.1,cadence,lags));means.append(float(C.mean()))
                valid_n+=int(np.sum(abs(P)>=.1));frames_n+=len(P)
            sources.append(dict(path=str(p),sha256=case['scalar_sha256']))
        groups[key]=E.summarize(moments,means,key,target,lags)
        contexts[key]=dict(arm='common',view='main',source='provisional_round2',forced_combined=row['forced_combined'],multiple_front_fraction=float(row['multiple_fronts_fraction']),
                           multiple_front_fraction_ci_lower=float(row['multiple_fronts_ci_lower']),multiple_front_fraction_ci_upper=float(row['multiple_fronts_ci_upper']),valid_fraction=valid_n/frames_n)
        case_groups[key]=g['noise']
        if len(groups)%20==0:print('round2 conditions',len(groups),flush=True)
    data=out/'supporting/data';result=export_analysis(groups,contexts,data)
    for name in ['condition_lag_metrics','analysis_lag_metrics','noise_tolerance']:
        (data/('round2_matched_P_'+name+'.csv')).write_text((data/(name+'.csv')).read_text())
    result.update(generated_utc=datetime.now(timezone.utc).isoformat(),source='provisional_round2',completed_campaign=False,
                  seed_counts=sorted({g['n_seeds'] for g in groups.values()}),scalar_sources=sources,
                  limitations=['Round2 common main observation only; new compact campaign is pending.','Historical stored P was M1/M1_initial; explicit round3 P=M1/M0 reconstructed from native moments. This changes endpoint eligibility, while retaining the same heading.','Whole-seed pointwise bootstrap; no hypothesis gates.','sqrt2 target lags rounded to nearest saved endpoint; ties down, target/actual both retained.'])
    if tracked:result['tracked']=fallback_tracked(out,groups,contexts,case_groups)
    meta=out/'supporting/metadata';meta.mkdir(parents=True,exist_ok=True)
    (meta/'analysis_receipt.json').write_text(json.dumps(result,indent=2,allow_nan=True)+'\n')
    return result


def fallback_tracked(out,groups,contexts,case_groups):
    roots=[ROOT/'data/runs/signature_figure4_regions_common_20261020',ROOT/'data/runs/signature_figure4_regions_E2_20261020']
    rows=[];msdrows=[];mles=[];switches=[];missing=[];raw_receipts=[]
    for key,cases in case_groups.items():
        allids=[];allhead=[];allC=[];trackP=[];dt=None;ncase=0;wholeP=[];whole_labels=[];Pcad=None;retained_seeds=[]
        for case in cases:
            path=next((r/'workers'/case['id']/'result.json' for r in roots if (r/'workers'/case['id']/'result.json').exists()),None)
            if path is None:continue
            record=json.loads(path.read_text());raw=next((s for s in record['sensitivities'] if s['threshold']==.5 and s['downsample_factor']==1 and not s['smoothed']),None)
            if raw is None:continue
            if record['case']['scalar_sha256']!=case['scalar_sha256']:raise ValueError('Retained region track belongs to a different scalar history')
            if digest(raw['path'])!=raw['sha256']:raise ValueError('Retained raw region frame archive changed')
            raw_receipts.append(dict(path=raw['path'],sha256=raw['sha256'],case_id=case['id']))
            with gzip.open(raw['path'],'rt') as stream:frames=json.load(stream)['result']['frames']
            if len(frames)<3:continue
            cadence=frames[1]['time']-frames[0]['time']
            if dt is not None and not np.isclose(cadence,dt):continue
            dt=cadence;ncase+=1
            allids.extend([-1]+[f['dominant_track_id'] if f['dominant_track_id'] is not None else -1 for f in frames])
            allhead.extend([{}]+[{int(r['track_id']):r['heading'] for r in f['regions'] if r['heading'] is not None} for f in frames]);allC.extend([frames[0]['circumference']]+[f['circumference'] for f in frames])
            # Seed labels differ between archives; detector-free data remains exact scalar main view.
            with np.load(case['scalar_path'],allow_pickle=False) as z:
                indices=np.searchsorted(z['time'],[f['time'] for f in frames]);nativeP=normalized_P(z)[indices]
                trackP.extend([complex(np.nan)]+list(nativeP))
                retained_seeds.append(dict(time=np.array([f['time'] for f in frames]),C=np.array([f['circumference'] for f in frames]),P=nativeP,
                                           ids=np.array([f['dominant_track_id'] if f['dominant_track_id'] is not None else -1 for f in frames]),
                                           headings=[{int(r['track_id']):r['heading'] for r in f['regions'] if r['heading'] is not None} for f in frames],cadence=cadence))
                keep=z['time']>=case['native_horizon']*.5;P=normalized_P(z)[keep]
                pcad=float(np.median(np.diff(z['time'][keep])))
                if Pcad is None:Pcad=pcad
                if np.isclose(Pcad,pcad):
                    wholeP.extend([complex(np.nan)]+list(P));whole_labels.extend([-1]+[ncase]*len(P))
        base=dict(model=key[0],C_label=key[1],f=key[2],C_mean=groups[key]['circumference'],source='provisional_round2_raw_tracks',n_seeds=ncase,**{k:v for k,v in contexts[key].items() if k in ['forced_combined','multiple_front_fraction','multiple_front_fraction_ci_lower','multiple_front_fraction_ci_upper']})
        if not ncase:
            rows.append(dict(base,status='unavailable_no_retained_raw_tracks'));missing.append(key);continue
        t=np.arange(len(allids))*dt
        for hold,fraction in [(5,.25),(2,.25),(10,.25),(5,.125),(5,.5)]:
            measure,msd=E.drift_jump(t,allC,allids,allhead,hold=hold,range_fraction=fraction,P=np.asarray(trackP))
            # No observed jumps imply empirical jump contribution zero; uncertainty is not zero.
            if measure['n_switches']==0 and measure['valid_time']>0:measure['lambda_eff']=0.;measure['lambda_s']=0.
            s=groups[key];j=E.analysis_lag(s);kobs=s['k'][j] if j is not None else math.nan
            kp=(2*np.pi/base['C_mean'])**2*measure['Ds']+measure['lambda_eff']
            angles=measure.pop('jump_angles',[])
            ci=tracked_uncertainty(retained_seeds,measure,s) if hold==5 and fraction==.25 and len(retained_seeds)==s['n_seeds'] else {}
            rows.append(dict(base,**measure,**ci,kobs=kobs,kpred=kp,kpred_over_kobs=kp/kobs if kobs>0 else math.nan,analysis_lag=s['lags'][j] if j is not None else math.nan,
                             jump_no_event_upper95=3/measure['valid_time'] if measure['n_switches']==0 and measure['valid_time']>0 else math.nan))
            if hold==5 and fraction==.25:msdrows.extend(dict(base,**m) for m in msd)
        # Whole-condition mixture likelihood; no seed boundaries enter increment pairs.
        for step in [1,2]:mles.append(dict(base,**E.wrapped_mixture_mle(wholeP,base['C_mean'],Pcad,step,whole_labels)))
        if len(rows)%50==0:print('tracked conditions',len(rows)//5,flush=True)
    data=out/'supporting/data';E.write_csv(data/'drift_decomposition.csv',rows);E.write_csv(data/'within_tenure_msd.csv',msdrows);E.write_csv(data/'detector_free_mle.csv',mles)
    return dict(conditions_with_retained_tracks=sum(r.get('flicker_frames')==5 and r.get('MSD_range_fraction')==.25 for r in rows),missing_conditions=missing,
                raw_track_sources=raw_receipts,note='Primary raw, unsmoothed round1 detector tracks reused at saved profile cadence; complete new compact trajectories pending. Pooled independent histories have seed labels and missing separators; no seed-boundary displacements.')


def simulate_null(n,C,Ds,lam,T,dt,seed,bernoulli=False):
    """Exact sampled continuous Brownian drift with independent uniform resets.

    Brownian displacement after a reset leaves the reset position uniform, so
    draw uniform final phase whenever >=1 event occurs in a sampling interval.
    This is exact for the transition law, including lambda*dt>1.
    """
    rng=np.random.default_rng(seed);steps=int(round(T/dt));theta=np.empty((n,steps+1));theta[:,0]=rng.uniform(-np.pi,np.pi,n)
    Dtheta=(2*np.pi/C)**2*Ds
    noise=rng.normal(0,np.sqrt(2*Dtheta*dt),(n,steps))
    reset=rng.random((n,steps)) < (lam*dt if bernoulli else -np.expm1(-lam*dt)) if lam else None
    for i in range(steps):
        theta[:,i+1]=theta[:,i]+noise[:,i]
        if reset is not None:
            mask=reset[:,i];theta[mask,i+1]=rng.uniform(-np.pi,np.pi,mask.sum())
    return theta


def null_moments(theta,dt,lags):
    result=np.zeros((len(theta),3,len(lags)))
    for j,lag in enumerate(lags):
        step=int(round(lag/dt))
        if step==0:
            result[:,0,j]=theta.shape[1];result[:,1,j]=theta.shape[1];result[:,2,j]=theta.shape[1];continue
        difference=theta[:,step:]-theta[:,:-step]
        result[:,0,j]=np.cos(difference).sum(1);result[:,1,j]=difference.shape[1];result[:,2,j]=difference.shape[1]
    return result


def nulls(out,only='all'):
    data=out/'supporting/data';data.mkdir(parents=True,exist_ok=True)
    previous=out/'supporting/metadata/null_analysis_receipt.json';summary=json.loads(previous.read_text()) if previous.exists() else {}
    rates=[];acf=[];exemplars={}
    T=24000.;dt=20.;target,lags=E.lag_grid('chemical',T+dt,dt)
    if only in ('all','A'):
        index=[]
        for kind in NULLS:
            for C in (20.,60.):
                Ds=.002 if kind=='drift' else 0.;lam=0. if kind=='drift' else (1/7000)*(C/40 if kind=='perimeter' else 1)
                # Same seed per null across sizes, including the deterministic displayed history.
                seed=E.stable_seed('A:'+kind);theta=simulate_null(1000,C,Ds,lam,T,dt,seed,bernoulli=True)
                moments=null_moments(theta,dt,lags);s=E.summarize(moments,np.full(1000,C),'A:'+kind+str(C),target,lags)
                ktheory=(2*np.pi/C)**2*Ds+(-np.log1p(-lam*dt)/dt if lam else 0)
                j=E.analysis_lag(s)
                rates.append(dict(model=kind,C_label=C,Ds=Ds,lambda_s=lam,analytic_k=ktheory,continuum_k=(2*np.pi/C)**2*Ds+lam,n_seeds=1000,seed=seed,
                                  observed_k=s['k'][j] if j is not None else math.nan,analysis_lag=s['lags'][j] if j is not None else math.nan))
                curve_lags=np.arange(0,T+dt,1000.);curve=E.summarize(null_moments(theta,dt,curve_lags),np.full(1000,C),'A:'+kind+str(C),curve_lags,curve_lags)
                for j,lag in enumerate(curve_lags):
                    lo,hi=E.percentile(curve['Cboot'][:,j]);acf.append(dict(model=kind,C_label=C,lag=lag,target_lag=lag,C=curve['C'][j],C_ci_lower=lo,C_ci_upper=hi,n_seeds=1000,analytic_C=np.exp(-ktheory*lag)))
                x=(np.arange(240)+.5)*C/240;position=np.mod(theta[0],2*np.pi)*C/(2*np.pi)
                distance=(x[None,:]-position[:,None]+C/2)%C-C/2
                prefix=f'{kind}_C{C:g}__';exemplars[prefix+'time']=np.arange(len(position))*dt;exemplars[prefix+'field']=np.exp(-.5*(distance/1.5)**2).astype('float32');exemplars[prefix+'circumference']=np.array(C)
                name=f'null_A_{kind}_C{C:g}.npz';np.savez_compressed(data/name,time=np.arange(len(position))*dt,field=exemplars[prefix+'field'],circumference=np.array(C));index.append(dict(file=name,null=kind,model=kind,C_label=C,C_mean=C,seed=seed))
        E.write_csv(data/'null_A_rates.csv',rates);E.write_csv(data/'null_A_acf.csv',acf);E.write_csv(data/'null_A_kymographs.csv',index)
        summary['A']=dict(replicates_per_null_size=1000,point_width=1.5,T=T,dt=dt,seed_policy='fixed seed per null shared at C20/C60; first trajectory displayed',
                          transition='Literal Bernoulli uniform relocation probability lambda*dt; otherwise Gaussian drift.',
                          finite_dt_uniform_jump_rate='-log(1-lambda*dt)/dt',acf_grid='0,1000,...,24000; actual1000replicate mean endpoint cosines')
    if only in ('all','B'):
        groups={};contexts={};checks=[]
        fs=np.geomspace(1e-4,1,41)
        for kind in NULLS:
            for C in E.CHEMICAL_SIZES:
                for f in fs:
                    Ds=8.8*f*f;lam=0. if kind=='drift' else 135*f**4*(C/40 if kind=='perimeter' else 1)
                    analytic=(2*np.pi/C)**2*Ds+lam
                    # Entire requested lag range out of bounds: retain simulation,
                    # but estimate precision at the nearest positive rate rather
                    # than inventing an eligible local-alpha point.
                    moments=[];n=0;SE=math.inf;pointC=None;j=None;batch=512
                    while True:
                        theta=simulate_null(batch,C,Ds,lam,T,dt,E.stable_seed(('B',kind,C,float(f),n)))
                        moments.append(null_moments(theta,dt,lags));n+=batch
                        mm=np.concatenate(moments);num=mm[:,0];den=mm[:,1];pointC=num.sum(0)/den.sum(0)
                        eligible=(pointC>.03)&(pointC<.995)
                        # Precision at closest-to-half admissible lag (same shared estimator).
                        ix=np.flatnonzero(eligible)
                        if not len(ix):
                            # There is no eligible k to plot. Guard against a noise-created
                            # eligible positive correlation by requiring 4096 histories
                            # when the analytic curve is below the lower C threshold.
                            if np.exp(-analytic*lags[0])<.03 and n<4096:
                                batch=4096-n;continue
                            break
                        j=int(ix[np.argmin(abs(pointC[ix]-.5))])
                        weights=np.random.default_rng(E.stable_seed(('B',kind,C,float(f)))).multinomial(n,np.full(n,1/n),size=E.DRAW_COUNT)
                        cb=(weights@num[:,j])/(weights@den[:,j])
                        with np.errstate(invalid='ignore',divide='ignore'):logk=np.log(-np.log(cb)/lags[j])
                        SE=float(np.std(logk[np.isfinite(logk)],ddof=1)) if np.isfinite(logk).mean()>=.975 else math.inf
                        if SE<=.01:break
                        # Use the observed precision to size the next independent
                        # batch; still verify the actual 2000-draw bootstrap.
                        batch=min(16384,max(512,int(np.ceil(n*(1.1*(SE/.01)**2-1))))) if np.isfinite(SE) else max(512,n)
                        # No artificial replicate cap; independent batches continue to precision.
                    s=E.summarize(mm,np.full(n,C),('B',kind,C,float(f)),target,lags)
                    groups[(kind,C,float(f))]=s;contexts[(kind,C,float(f))]=dict(source='simulated_null',forced_combined=False)
                    j=E.analysis_lag(s)
                    checks.append(dict(model=kind,C_label=C,f=f,n_seeds=n,analytic_k=analytic,
                                       selected_lag=s['lags'][j] if j is not None else math.nan,observed_k=s['k'][j] if j is not None else math.nan,
                                       relative_deviation=s['k'][j]/analytic-1 if j is not None else math.nan,SE_log_k=SE,
                                       status='precision_met' if j is not None and SE<=.01 else 'unavailable_no_eligible_lag'))
                print('null B',kind,C,flush=True)
        # Null labels use chemical five-size neighborhoods explicitly.
        result=export_analysis(groups,contexts,data,prefix='null_B_');E.write_csv(data/'null_B_analytic_agreement.csv',checks)
        finite=[r for r in checks if np.isfinite(r['relative_deviation'])]
        summary['B']=dict(conditions=len(checks),max_replicates=max(r['n_seeds'] for r in checks),eligible_conditions=len(finite),
                          largest_absolute_relative_deviation=max(finite,key=lambda r:abs(r['relative_deviation'])) if finite else None,
                          precision_unmet=[r for r in checks if r['status']=='precision_met' and r['SE_log_k']>.01],
                          omitted_no_eligible_lag=sum(r['status']!='precision_met' for r in checks))
    if 'B' in summary:summary['B']['transition']='Exact continuous Brownian + Poisson uniform resets sampled atdt20; reset probability1-exp(-lambda*dt).'
    summary.update(generated_utc=datetime.now(timezone.utc).isoformat(),transition={'A':'Bernoulli lambda*dt','B':'Exact Poisson 1-exp(-lambda*dt)'},
                   reason='For B, lambda*dt exceeds1 at strong inputs; exact uniform-reset transition is defined at all strengths. A retains the requested Bernoulli discretization.',
                   bootstrap_resamples=2000,analytic_rates_plotted=False)
    meta=out/'supporting/metadata';meta.mkdir(parents=True,exist_ok=True);(meta/'null_analysis_receipt.json').write_text(json.dumps(summary,indent=2,allow_nan=True)+'\n')
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=DEFAULT_OUT);p.add_argument('--source',choices=['round2','compact','none'],default='round2');p.add_argument('--campaign',type=Path);p.add_argument('--nulls',choices=['none','all','A','B'],default='none');p.add_argument('--skip-tracked',action='store_true');p.add_argument('--forcing-candidates',type=Path);a=p.parse_args()
    if a.forcing_candidates:
        if a.campaign is None:p.error('--forcing-candidates requires --campaign')
        forcing_candidates(a.campaign,a.forcing_candidates);return
    if a.source=='round2':fallback(a.output,not a.skip_tracked)
    if a.source=='compact':
        if a.campaign is None:p.error('--source compact requires --campaign')
        compact(a.output,a.campaign,not a.skip_tracked)
    if a.nulls!='none':nulls(a.output,a.nulls)


if __name__=='__main__':main()
