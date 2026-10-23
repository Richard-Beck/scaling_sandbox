"""Rescore existing Figure 4 histories; no simulation or detector-free refit.

Run --reanalyze in SLURM. Default is a portable redraw of the revised records.
The preceding round-3 exports and all compact simulation evidence are immutable.
"""
from __future__ import annotations
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import math
from multiprocessing import get_context
from pathlib import Path
import shutil
import numpy as np
import figure4_round3_analysis_20261022 as A
import figure4_round3_estimators_20261022 as E
import figure4_round3_render_20261022 as R
import rebuild_figure4_round3_20261022 as B

OUT=R.ROOT/'figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4'
CAMPAIGN=R.ROOT/'data/runs/signature_figure4_round3_compact_20261022'


def key(row):return row['model'],float(row['C_label']),float(row['f'])


def initialize(out):
    data=out/'supporting/data';meta=out/'supporting/metadata'
    data.mkdir(parents=True,exist_ok=True);meta.mkdir(parents=True,exist_ok=True)
    if (meta/'baseline_round3_sources.json').exists():return
    sources={}
    for p in (R.OUT/'supporting/data').iterdir():
        if p.is_file():shutil.copy2(p,data/p.name);sources[str(p.relative_to(R.REPO))]=R.sha(p)
    for name in ['analysis_receipt.json','null_analysis_receipt.json','campaign_manifest.json',
                 'campaign_progress_receipt.json','campaign_final_allocation.json','campaign_final_precision.json',
                 'model_D_provenance.json']:
        p=R.OUT/'supporting/metadata'/name
        if p.exists():shutil.copy2(p,meta/name);sources[str(p.relative_to(R.REPO))]=R.sha(p)
    R.write_json(meta/'baseline_round3_sources.json',dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        baseline=str(R.OUT.relative_to(R.REPO)),source_sha256=sources))


def condition(task):
    """One condition: exactly the existing cohort, tenures, MSD and RNG seeds."""
    identity,records,old=task
    seeds=[];moments=[];means=[]
    for row in records:
        path=Path(row['compact'])
        if R.sha(path)!=row['compact_sha256']:raise ValueError('Changed compact evidence: '+str(path))
        with np.load(path,allow_pickle=False) as z:
            t=z['time'].astype(float);C=z['circumference'].astype(float)
            P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float)
            offsets=z['region_offset'];ids=z['region_id'];heading=z['region_heading']
            maps=[{int(ids[j]):float(heading[j]) for j in range(int(a),int(b)) if np.isfinite(heading[j])}
                  for a,b in zip(offsets[:-1],offsets[1:])]
            dt=float(np.median(np.diff(t)));target,lags=E.lag_grid(identity[0],t[-1]-t[0],dt)
            moments.append(E.seed_moments(P,abs(P)>=.1,dt,lags));means.append(float(C.mean()))
            seeds.append(dict(time=t,P=P,C=C,ids=z['dominant_track_id'].astype(int),headings=maps,cadence=dt))
    summary=E.summarize(moments,means,('compact',identity),target,lags)
    base={k:old[k] for k in ['model','C_label','f','C_mean','n_seeds','arm','view','source',
          'forced_combined','valid_fraction','multiple_front_fraction','multiple_front_fraction_ci_lower','multiple_front_fraction_ci_upper']}
    ids=[];maps=[];cs=[];ps=[]
    for seed in seeds:
        ids.extend([-1]+list(seed['ids']));maps.extend([{}]+seed['headings'])
        cs.extend([seed['C'][0]]+list(seed['C']));ps.extend([complex(np.nan)]+list(seed['P']))
    time=np.arange(len(ids))*dt;j=E.analysis_lag(summary)
    kobs=summary['k'][j] if j is not None else math.nan
    if np.isfinite(kobs) and not np.isclose(kobs,float(old['kobs']),rtol=1e-13,atol=1e-16):raise ValueError('Observed rate changed')
    results=[];draws=None
    for hold,fraction in [(5,.25),(2,.25),(10,.25),(5,.125),(5,.5)]:
        measure,_=E.drift_jump(time,np.asarray(cs),ids,maps,hold=hold,range_fraction=fraction,P=np.asarray(ps))
        measure.pop('jump_angles',None);ci={}
        if hold==5 and fraction==.25:
            for name in ['Ds','lambda_s','n_switches','n_scored_switches','mean_cos_jump','valid_time','MSD_max_lag']:
                if not np.isclose(measure[name],float(old[name]),rtol=1e-12,atol=1e-14,equal_nan=True):
                    raise ValueError(f'Unintended change to {name}: {identity}')
            ci=A.tracked_uncertainty(seeds,measure,summary,return_draws=True);draws=ci.pop('_draws')
        predicted=(2*np.pi/summary['circumference'])**2*measure['Ds']+measure['lambda_eff']
        results.append(dict(base,**measure,**ci,kobs=kobs,kpred=predicted,
                            kpred_over_kobs=predicted/kobs if kobs>0 else math.nan,
                            analysis_lag=summary['lags'][j] if j is not None else math.nan))
    return identity,results,draws


def enrich(out):
    """Join fixed detector-free fits, normalization and regime-specific checks."""
    data=out/'supporting/data';rows=R.read_csv(data/'drift_decomposition.csv')
    main=[r for r in rows if float(r['flicker_frames'])==5 and float(r['MSD_range_fraction'])==.25]
    tracked={key(r):r for r in main};mles=R.read_csv(data/'detector_free_mle.csv');fits=[]
    for row in mles:
        t=tracked[key(row)];kobs=float(t['kobs']);D=R.n(row,'Dtheta');eff=R.n(row,'lambda_eff')
        predicted=D+eff if row['status']=='resolved' else math.nan
        fits.append(dict(row,kobs=kobs,kpred=predicted,kpred_over_kobs=predicted/kobs if kobs>0 else math.nan,
                         analysis_lag=t['analysis_lag'],method='detector_free_wrapped_normal_uniform'))
    E.write_csv(data/'detector_free_decomposition.csv',fits)
    for row in rows:
        model,size,f=key(row);ref=tracked[(model,1. if model=='spring_ring' else 20.,f)]
        denominator=float(ref['Ds']);value=float(row['Ds'])/denominator if denominator>0 else math.nan
        row.update(Ds_over_smallest=value,Ds_reference_C_label=ref['C_label'],Ds_reference=ref['Ds'])
        if row in main:
            if size==float(ref['C_label']):lo=hi=1.
            else:
                with np.load(data/'bootstrap'/filename(key(row))) as z,np.load(data/'bootstrap'/filename(key(ref))) as rz:
                    with np.errstate(divide='ignore',invalid='ignore'):boot=z['Ds']/rz['Ds']
                lo,hi=E.percentile(boot,True)
            row.update(Ds_over_smallest_ci_lower=lo,Ds_over_smallest_ci_upper=hi)
    E.write_csv(data/'drift_decomposition.csv',rows)
    # The mean jump loss is conditional on scorable events. Exact worst-case
    # angle-loss bounds (0..2 for every missing event) show that limitation.
    coverage=[]
    for row in main:
        total=int(row['n_switches']);scored=int(row['n_scored_switches']);exposure=float(row['valid_time'])
        loss=scored*(1-float(row['mean_cos_jump'])) if scored else 0.
        coverage.append(dict(model=row['model'],C_label=row['C_label'],f=row['f'],n_switches=total,n_scored_switches=scored,
                             jump_scored_fraction=scored/total if total else math.nan,
                             lambda_eff_missing_angle_lower=loss/exposure if exposure>0 else math.nan,
                             lambda_eff_missing_angle_upper=(loss+2*(total-scored))/exposure if exposure>0 else math.nan,
                             lambda_eff=row['lambda_eff']))
    E.write_csv(data/'jump_scoring_coverage.csv',coverage)
    alpha={key(r):R.n(r,'alpha') for r in R.read_csv(data/'local_alpha.csv')};audit=[]
    for method,rr in [('scorable_events',main),('detector_free',[r for r in fits if int(r['detector_free_step'])==1])]:
        for regime,predicate in [('all',lambda r:True),('multiple_front_majority',lambda r:float(r['multiple_front_fraction'])>.5),
                                 ('alpha_below_2',lambda r:alpha.get(key(r),math.nan)<2),
                                 ('multiple_front_majority_and_alpha_below_2',lambda r:float(r['multiple_front_fraction'])>.5 and alpha.get(key(r),math.nan)<2)]:
            selected=[r for r in rr if predicate(r)];estimable=[r for r in selected if np.isfinite(float(r['kpred_over_kobs']))]
            audit.append(dict(method=method,regime=regime,conditions=len(selected),estimable=len(estimable),
                              within_15_percent=sum(abs(float(r['kpred_over_kobs'])-1)<=.15 for r in estimable),
                              positive_jump=sum(float(r['lambda_eff'])>0 for r in estimable)))
    E.write_csv(data/'closure_by_regime.csv',audit)
    detector={key(r):r for r in fits if int(r['detector_free_step'])==1}
    zero=[r for r in main if int(r['n_switches'])==0]
    receipt=dict(generated_utc=datetime.now(timezone.utc).isoformat(),no_new_simulations=True,
                 conditions=len(main),counted_replacements=sum(int(r['n_switches']) for r in main),
                 scored_replacements=sum(int(r['n_scored_switches']) for r in main),
                 jump_finite=sum(bool(np.isfinite(float(r['lambda_eff']))) for r in main),
                 zero_tracked_switch_conditions=len(zero),
                 zero_tracked_switch_positive_detector_free=sum(R.n(detector[key(r)],'lambda_eff')>0 for r in zero),
                 closure_by_regime=audit,
                 angle_rule='Original per-event five-frame P windows; condition mean uses all scorable events, with no all-events-scored gate.',
                 assumption='The scorable-event angle distribution estimates the angle distribution of all counted events. Coverage and missing-angle bounds are reported; this assumption is not a validation of missingness at random.',
                 detector_free='Existing resolved wrapped-normal plus uniform fits reused without refitting; shortest and twice cadence retained. It is a model-dependent diagnostic, not an independent jump detector.',
                 normalized_drift='Ds / Ds at smallest designed size, separately by model and input; reference 20 or Spring 1. Independent size ensembles retain deterministic whole-seed bootstrap draws.')
    R.write_json(out/'supporting/metadata/scoring_revision_receipt.json',receipt)
    return receipt


def filename(identity):return f'{identity[0]}_C{identity[1]:g}_f{identity[2]:g}.npz'


def reanalyze(out,campaign):
    initialize(out);records=A.read_csv(campaign/'results/compact_index.csv');per={}
    for row in records:per.setdefault(row['condition_id'],[]).append(row)
    manifest=json.loads((campaign/'manifest.json').read_text());old={key(r):r for r in R.read_csv(R.OUT/'supporting/data/drift_decomposition.csv') if float(r['flicker_frames'])==5 and float(r['MSD_range_fraction'])==.25}
    tasks=[]
    for c in manifest['conditions']:
        identity=c['model'],float(c['size']),float(c['fraction']);rr=sorted(per[c['condition_id']],key=lambda r:int(r['replicate']))
        assert len(rr)==int(old[identity]['n_seeds']) and {int(r['replicate']) for r in rr}==set(range(len(rr)))
        tasks.append((identity,rr,old[identity]))
    workers=A._compact_tracked_worker_limit();rows=[];bootstrap=out/'supporting/data/bootstrap';bootstrap.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
        for i,(identity,rr,draws) in enumerate(pool.map(condition,tasks,chunksize=1),1):
            rows.extend(rr);np.savez_compressed(bootstrap/filename(identity),**draws)
            if i%10==0:print(f'Rescored {i}/{len(tasks)} conditions',flush=True)
    E.write_csv(out/'supporting/data/drift_decomposition.csv',rows)
    receipt=enrich(out)
    receipt.update(workers=workers,compact_index_sha256=R.sha(campaign/'results/compact_index.csv'),
                   source_sha256={str(p.relative_to(R.REPO)):R.sha(p) for p in [Path(__file__),Path(E.__file__),Path(A.__file__)]})
    R.write_json(out/'supporting/metadata/scoring_revision_receipt.json',receipt)
    print(json.dumps(receipt,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--reanalyze',action='store_true')
    parser.add_argument('--output',type=Path,default=OUT);parser.add_argument('--campaign',type=Path,default=CAMPAIGN)
    parser.add_argument('--analysis-only',action='store_true');args=parser.parse_args()
    if args.reanalyze:reanalyze(args.output,args.campaign)
    if not args.analysis_only:
        if not (args.output/'supporting/metadata/scoring_revision_receipt.json').exists():raise FileNotFoundError('Run --reanalyze in SLURM first')
        B.rebuild(args.output)


if __name__=='__main__':main()
