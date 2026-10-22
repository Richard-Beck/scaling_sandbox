"""Round-2 Figure 4 interim, existing scalar outputs only; no integrations.

Reuses the accepted observation selectors and exact condition bootstrap seeds.
Requested lags not on a view's saved-time grid remain explicitly unavailable.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import signature_figure4_analysis_20261020 as A
from signature_figure4_e0_20261020 import VIEWS, select

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/signature_figure4_final_statistics_20261020'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4'
MODELS=('wave_pinning','goryachev','holmes_model3','spring_ring')
COLORS=dict(wave_pinning='#2878b5',goryachev='#c83b3b',holmes_model3='#36944c',spring_ring='#c59a25')
NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes',spring_ring='Spring')

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for x in iter(lambda:f.read(4<<20),b''): h.update(x)
    return h.hexdigest()

def read(p): return json.loads(Path(p).read_text())
def rows(p):
    with open(p,newline='') as f:return list(csv.DictReader(f))
def write(p,values):
    fields=list(dict.fromkeys(k for r in values for k in r))
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(values)
def num(x):return float(x) if x not in ('',None) else math.nan
def yes(x):return x is True or x=='True'
def forced_pair(a,b):
    if yes(a) or yes(b):return True
    return False if a in (False,'False') and b in (False,'False') else None
def q(x):
    x=np.asarray(x); x=x[~np.isnan(x)]
    return tuple(np.quantile(x,[.025,.975],method='inverted_cdf')) if len(x) else (math.nan,math.nan)
def rate(c,lag):
    c=np.asarray(c)
    with np.errstate(divide='ignore',invalid='ignore'):return np.where(c>0,-np.log(np.minimum(c,1))/lag,np.nan)

def seed_lags(P,valid,cadence,lags):
    """Pair counts and cosine sums within the declared view's endpoints."""
    out=np.full((3,len(lags)),np.nan)
    for j,lag in enumerate(lags):
        step=int(round(lag/cadence))
        if step<1 or not math.isclose(step*cadence,lag,abs_tol=1e-7) or step>=len(P):continue
        mask=valid[:-step]&valid[step:]; count=int(mask.sum())
        out[:,j]=[np.real(P[:-step][mask]*np.conj(P[step:][mask])/(abs(P[:-step][mask])*abs(P[step:][mask]))).sum(),count,len(P)-step]
    return out

def crossing(f,c):
    """First downward 0.5 crossing; missing/forced gaps give interval bounds."""
    f=np.asarray(f);c=np.asarray(c);missing=False;last=None
    for i,value in enumerate(c):
        if not np.isfinite(value):missing=True;continue
        if value<=.5:
            if i==0:return dict(status='left_censored',estimate=None,lower=0.,upper=float(f[0]))
            if last is None or missing:return dict(status='bounded_missing_or_forced_gap',estimate=None,lower=float(f[last]) if last is not None else 0.,upper=float(f[i]))
            w=(c[last]-.5)/(c[last]-value);x=float(np.exp(np.log(f[last])+w*np.log(f[i]/f[last])))
            return dict(status='resolved',estimate=x,lower=x,upper=x)
        last=i
    return dict(status='right_censored' if not missing else 'unavailable_or_excluded_inputs',estimate=None,lower=float(f[-1]) if not missing else 0.,upper=math.inf)

def condition(group):
    lags=np.array([250,500,1000,2000,4000,8000],float)/(10 if group['model']=='spring_ring' else 1)
    key=group.get('bootstrap_pairing_key',f"{group['model']}:{group['C_label']}:{group['f']}")
    seed=2026102004+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)
    weights=A.bootstrap_weights(len(group['noise']),2000,seed)
    samples={v:[] for v in VIEWS};circ={v:[] for v in VIEWS};receipts=[]
    for r in group['noise']:
        p=Path(r['scalar_path']);digest=sha(p)
        if digest!=r['scalar_sha256']:raise ValueError('Native scalar hash changed '+str(p))
        receipts.append(dict(path=str(p),sha256=digest,id=r['id'],seed=r['seed']))
        with np.load(p,allow_pickle=False) as z:
            t=z['time'];P=z['P_real']+1j*z['P_imag'];C=np.broadcast_to(z['circumference'],t.shape)
            if not np.allclose(np.diff(t),r['scalar_cadence']):raise ValueError('Scalar cadence changed')
            for v in VIEWS:
                ix=select(t,v,r['native_horizon']);pview=P[ix];dt=float(t[ix[1]]-t[ix[0]])
                samples[v].append(seed_lags(pview,abs(pview)>=.1,dt,lags));circ[v].append(C[ix].mean())
    result={}
    for v in VIEWS:
        a=np.asarray(samples[v]);n=a[:,0,:];d=a[:,1,:];possible=a[:,2,:]
        with np.errstate(divide='ignore',invalid='ignore'):
            c=n.sum(0)/d.sum(0);boot=(weights@n)/(weights@d);cov=d.sum(0)/possible.sum(0)
        cb=weights@np.asarray(circ[v])/len(circ[v]);k=rate(c,lags);kb=rate(boot,lags)
        result[v]=dict(lags=lags,c=c,boot=boot,cov=cov,k=k,kboot=kb,Cmean=np.mean(circ[v]),Cboot=cb,
            pairs=d.sum(0),possible=possible.sum(0),seed_num=n,seed_count=d,seed_possible=possible)
    return result,receipts

def ratio_record(s,l,sm,lm,j):
    usable=bool(min(s['cov'][j],l['cov'][j])>=.5 and (q(1-s['boot'][:,j])[0]>0 or q(1-l['boot'][:,j])[0]>0))
    ks,kl=s['k'][j],l['k'][j]
    with np.errstate(divide='ignore',invalid='ignore'):b=s['kboot'][:,j]/l['kboot'][:,j]
    valid=float(np.mean(~np.isnan(b)))
    ci=q(b) if valid>=.975 else (0.,math.inf)
    geom=(l['Cmean']/s['Cmean'])**2; gb=(l['Cboot']/s['Cboot'])**2
    gci=q(b/gb) if valid>=.975 else (0.,math.inf)
    point=ks/kl if kl>0 and np.isfinite(ks) else math.nan
    return dict(lag=float(s['lags'][j]),usable=usable,rho_k=point,rho_k_lower=point,rho_k_upper=point,
        rho_k_ci_lower=ci[0],rho_k_ci_upper=ci[1],rho_k_ci_upper_unbounded=math.isinf(ci[1]),
        rho_geo=point/geom,rho_geo_ci_lower=gci[0],rho_geo_ci_upper=gci[1],
        geometry_expectation=geom,coverage_small=s['cov'][j],coverage_large=l['cov'][j],
        C_small=s['c'][j],C_large=l['c'][j],k_small=ks,k_large=kl,
        bootstrap_valid_fraction=valid,status='resolved' if np.isfinite(point) else 'unavailable_nonpositive_C_or_rate')

def secondary(s,l):
    values=[]
    for a in (s,l):
        mask=(a['c']>=.2)&(a['c']<=.98)&(a['cov']>=.5)
        x=a['lags'][mask]
        if len(x)<2:values.append((math.nan,np.full(2000,np.nan),len(x),math.nan));continue
        center=x-x.mean();y=-np.log(a['c'][mask]);slope=float(y@center/(center@center))
        with np.errstate(invalid='ignore',divide='ignore'):bs=-np.log(a['boot'][:,mask])@center/(center@center)
        bs[bs<=0]=np.nan
        values.append((slope,bs,len(x),float(y.mean()-slope*x.mean())))
    a,b=values
    with np.errstate(divide='ignore',invalid='ignore'):ratio=a[1]/b[1]
    return dict(secondary_rho_k=a[0]/b[0] if a[0]>0 and b[0]>0 else math.nan,secondary_ci_lower=q(ratio)[0],secondary_ci_upper=q(ratio)[1],
        secondary_small_slope=a[0],secondary_large_slope=b[0],secondary_small_lags=a[2],secondary_large_lags=b[2],
        secondary_small_intercept=a[3],secondary_large_intercept=b[3])

def run(out):
    for d in ('vectors','data','metadata','viewers'):(out/'supporting'/d).mkdir(parents=True,exist_ok=True)
    data=out/'supporting/data';meta=out/'supporting/metadata'
    complete=read(RUN/'results/final_analysis/complete.json');old=RUN/'results/final_analysis/condition_metrics.csv'
    if sha(old)!=complete['outputs'][old.name]:raise ValueError('Accepted condition metrics changed')
    metrics=rows(old);spec=read(complete['input_spec_path']);groups={};sources=[]
    for stage in ('E2','E3','E4'):
        root=Path(spec['stages'][stage]['root']);arms=read(root/'results/arms.json')
        for a in arms['arms']:
            if float(a['sigma'])!=1:continue
            p=Path(a['path'])/'analysis_manifest.json'
            if sha(p)!=a['analysis_manifest_sha256']:raise ValueError('Manifest changed')
            sources.append(dict(path=str(p),sha256=sha(p)))
            for g in read(p)['groups']:groups[(stage,g['model'],float(g['C_label']),float(g['f']))]=g
    lookup={};cache={};receipts=[];lagrows=[];bootout={};seedrows=[]
    for r in metrics:
        arm=r['arm'];model=r['model'];C=float(r['C_label']);f=float(r['f']);view=r['view']
        origin=r.get('common_observation_arm') or r['selected_observation_arm']
        stage='E4' if r['selected_observation_arm']=='preselected_E4' else 'E3' if origin=='common_new_E3_size' else 'E2'
        key=(stage,model,C,f)
        if key not in cache:
            cache[key],receipt=condition(groups[key]);receipts.extend(receipt)
            print('condition',len(cache),key,flush=True)
        a=cache[key][view];lookup[(arm,model,C,f,view)]=(r,a)
        base=dict(arm=arm,model=model,C_label=C,f=f,view=view,forced_combined=r['forced_combined'],regime=r['regime'],n_seeds=r['n_seeds'])
        for j,lag in enumerate(a['lags']):
            resolved=math.isfinite(a['c'][j]);cvalid=float(np.isfinite(a['boot'][:,j]).mean())
            ci=q(a['boot'][:,j]) if cvalid>=.975 or not resolved else (-1.,1.)
            ki=(float(rate(ci[1],lag)),float(rate(ci[0],lag)) if ci[0]>0 else math.inf) if resolved else (math.nan,math.nan)
            lagrows.append(dict(base,lag=lag,C=a['c'][j],C_ci_lower=ci[0],C_ci_upper=ci[1],C_valid_bootstrap_fraction=cvalid,one_minus_C_ci_lower=1-ci[1],one_minus_C_ci_upper=1-ci[0],
                k=a['k'][j],k_ci_lower=ki[0],k_ci_upper=ki[1],k_valid_bootstrap_fraction=float(np.isfinite(a['kboot'][:,j]).mean()),coverage=a['cov'][j],eligible_pairs=a['pairs'][j],possible_pairs=a['possible'][j],
                status='available' if resolved else 'unavailable_saved_grid_or_window_or_coverage'))
        bkey=f'{arm}|{model}|{C}|{f}|{view}';bootout[bkey+'|C']=a['boot'];bootout[bkey+'|circumference']=a['Cboot']
        for i,case in enumerate(groups[key]['noise']):
            for j,lag in enumerate(a['lags']):seedrows.append(dict(base,seed=case['seed'],case_id=case['id'],lag=lag,cosine_sum=a['seed_num'][i,j],eligible_pairs=a['seed_count'][i,j],possible_pairs=a['seed_possible'][i,j]))
    ratios=[];alllags=[]
    for arm in ('common','longest'):
        for model in MODELS:
            sizes=sorted({k[2] for k in lookup if k[:2]==(arm,model)});sm,lg=sizes[0],sizes[-1]
            for view in VIEWS:
                fs=sorted({k[3] for k in lookup if k[:2]==(arm,model) and k[2] in (sm,lg)})
                for f in fs:
                    sk=(arm,model,sm,f,view);lk=(arm,model,lg,f,view)
                    base=dict(arm=arm,model=model,view=view,f=f,small_C_label=sm,large_C_label=lg)
                    if sk not in lookup or lk not in lookup:
                        ratios.append(dict(base,evaluable=False,status='unavailable_unmatched_designed_input'));continue
                    sr,s=lookup[sk];lr,l=lookup[lk]
                    base.update(forced_combined=forced_pair(sr['forced_combined'],lr['forced_combined']),large_multiple_front_fraction=lr['multiple_fronts_fraction'],large_regime=lr['regime'],small_regime=sr['regime'])
                    rr=[ratio_record(s,l,sr,lr,j) for j in range(6)];usable=[j for j in range(6) if rr[j]['usable']]
                    candidates=[j for j in usable if s['c'][j]<=.98 and s['c'][j]>0 and l['c'][j]>0]
                    drift=[rr[j]['rho_k'] for j in usable if np.isfinite(rr[j]['rho_k']) and rr[j]['rho_k']>0]
                    driftfold=max(drift)/min(drift) if len(drift)>1 else math.nan
                    for x in rr:alllags.append(dict(base,**x,non_exponential=driftfold>1.5,lag_drift_fold=driftfold))
                    if not candidates:
                        ratios.append(dict(base,evaluable=False,status='unavailable_no_qualifying_primary_lag',non_exponential=driftfold>1.5,**secondary(s,l)));continue
                    j=candidates[0];chosen=rr[j].copy();chosen.update(primary_lag=chosen.pop('lag'),evaluable=True,non_exponential=driftfold>1.5,lag_drift_fold=driftfold)
                    if all(q(1-l['boot'][:,u])[0]<=0 for u in usable):
                        upper=q(l['kboot'][:,j])[1];lower=q(s['kboot'][:,j])[0]
                        bound=s['k'][j]/upper if upper>0 else math.inf;conservative=lower/upper if upper>0 else math.inf
                        chosen.update(status='one_sided_lower_bound',rho_k=math.nan,rho_k_lower=bound,rho_k_upper=math.inf,rho_k_ci_lower=conservative,rho_k_ci_upper=math.inf,rho_k_ci_upper_unbounded=True,
                            rho_geo=math.nan,rho_geo_ci_lower=conservative/chosen['geometry_expectation'],rho_geo_ci_upper=math.inf)
                    ratios.append(dict(base,**chosen,**secondary(s,l)))
    tolerances=[]
    for arm in ('common','longest'):
        for model in MODELS:
            for C in sorted({k[2] for k in lookup if k[:2]==(arm,model)}):
                for view in VIEWS:
                    pairs=sorted([(k[3],v) for k,v in lookup.items() if k[:3]==(arm,model,C) and k[4]==view])
                    for inclusive in (False,True):
                        fs=np.array([p[0] for p in pairs]);cs=np.array([p[1][1]['c'][3] if (inclusive or p[1][0]['forced_combined']=='False') and p[1][1]['cov'][3]>=.5 else np.nan for p in pairs])
                        boots=np.array([p[1][1]['boot'][:,3] if np.isfinite(cs[i]) else np.full(2000,np.nan) for i,p in enumerate(pairs)])
                        est=crossing(fs,cs);br=[crossing(fs,boots[:,b]) for b in range(2000)]
                        lo=float(np.quantile([b['lower'] for b in br],.025,method='inverted_cdf'));hi=float(np.quantile([b['upper'] for b in br],.975,method='inverted_cdf'))
                        tolerances.append(dict(arm=arm,model=model,C_label=C,view=view,forced_inclusive=inclusive,reference_lag=200 if model=='spring_ring' else 2000,
                            fhalf=est['estimate'],fhalf_lower=est['lower'],fhalf_upper=est['upper'],fhalf_ci_lower=lo,fhalf_ci_upper=hi,fhalf_ci_upper_unbounded=math.isinf(hi),status=est['status'],
                            designed_inputs=json.dumps(fs.tolist()),excluded_or_unavailable_inputs=json.dumps(fs[~np.isfinite(cs)].tolist()),bootstrap_resamples=2000))
    write(data/'condition_lag_metrics.csv',lagrows);write(data/'seed_lag_summaries.csv',seedrows);write(data/'size_ratios.csv',ratios);write(data/'all_lag_size_ratios.csv',alllags);write(data/'noise_tolerance.csv',tolerances)
    np.savez_compressed(data/'whole_seed_bootstrap.npz',**bootout)
    write(data/'old_tau_condition_metrics.csv',metrics)
    receipt=dict(generated_utc=datetime.now(timezone.utc).isoformat(),existing_outputs_only=True,no_simulations=True,bootstrap_resamples=2000,condition_groups=len(cache),condition_view_rows=len(metrics),lag_rows=len(lagrows),ratio_rows=len(ratios),
        input_condition_metrics=dict(path=str(old),sha256=sha(old)),manifests=sources,scalar_sources=list({r['path']:r for r in receipts}.values()),
        qualifications=['All requested lags use exact saved endpoints within each declared view. Off-grid lags are unavailable; no interpolation in time. All within-view separations are used, without the old half-window maximum.',
        'Whole-seed multinomial resampling, exact original condition key and random seed; same weights across lags and views. Independent size ensembles.',
        'Primary selection is fixed from observed curves, not reselected within bootstrap draws. Pointwise 95% percentile intervals use inverted empirical CDF.',
        'Rate defined only for positive autocorrelation. Roundoff C>1 clips to1. Zero-rate denominator produces an unbounded ratio interval.',
        'If more than2.5% of ratio bootstrap draws have undefined rate ratios, its CI is conservatively expanded to [0,infinity] rather than conditioning on defined draws.',
        'One-sided bound uses observed k_small / upper95 k_large; conservative lower endpoint uses lower95 k_small / upper95 k_large.',
        'Secondary slopes use ordinary least-squares with free intercept, observed eligible .2<=C<=.98 lags; no extrapolated primary estimator.',
        'Noise tolerance uses first downcrossing of absolute C=.5 at fixed reference lag, log-input interpolation; missing/forced gaps remain bounded.',
        'Geometric normalization is relative to single-front expectation using actual mean circumference; it is not site diffusivity.',
        'deBelly absent from the original four-model Figure4 design; no new arm invented.'])
    owned=('condition_lag_metrics.csv','seed_lag_summaries.csv','size_ratios.csv','all_lag_size_ratios.csv','noise_tolerance.csv','whole_seed_bootstrap.npz','old_tau_condition_metrics.csv')
    receipt['outputs']={name:sha(data/name) for name in owned}
    (meta/'analysis_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return ratios,alllags,lagrows,tolerances

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUT);args=p.parse_args();run(args.output)
