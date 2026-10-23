"""Read-only full-campaign analysis; never integrates a noisy model.

Certifiable births are a conservative subset, not an estimate of all births.
The compact exporter discarded masks and parents, which prevents exact splits.
"""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math, os
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np

REPO=Path(__file__).resolve().parents[4]
BASE=REPO/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
OLD=BASE/'revision_round4_20261023/figure4'
OUT=BASE/'revision_round5_20261023/figure4'
CAMPAIGN=REPO/'polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022'
COMPONENTS=('motion','coexisting_amount','birth_loss','moment_residual')

def read_csv(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def write_csv(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None:fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def wrapped(x):return (x+np.pi)%(2*np.pi)-np.pi

def certified_no_overlap(h,w,C,old_h,old_w,old_C):
    """Enclose each contiguous <half-ring component about its weighted centroid.

    Positive weights place its circular centroid in the component's shortest
    support arc. A radius of full cell-quadrature width (not half-width) encloses
    every occupied cell, including half-cell boundaries and float32 rounding.
    Material cells live on a fixed periodic angular grid, even in Spring.
    """
    if not len(old_h):return True
    if not np.isfinite(h+w+C) or C<=0 or w>=C/2:return False
    if np.any(~np.isfinite(old_h+old_w+old_C)) or np.any(old_w>=old_C/2):return False
    radius=2*np.pi*w/C
    old_radius=2*np.pi*old_w/old_C
    return bool(np.all(np.abs(wrapped(h-old_h))>radius+old_radius+2e-6))

def dense(z):
    off=z['region_offset'];n=len(off)-1;counts=np.diff(off);k=max(1,int(counts.max()))
    ids=np.full((n,k),-1,dtype=np.int32);A=np.zeros((n,k));e=np.zeros((n,k),complex)
    finite=np.ones(n,bool)
    for j in range(k):
        good=counts>j;ix=off[:-1][good]+j
        ids[good,j]=z['region_id'][ix];A[good,j]=z['region_amount'][ix]
        e[good,j]=np.exp(1j*z['region_heading'][ix].astype(float))
        finite[good]&=np.isfinite(A[good,j])&np.isfinite(e[good,j])&z['region_amount_available'][ix]
    return ids,A,e,finite

def decomposition(z,lag):
    """Exact endpoint Shapley moment changes and signed heading-loss projection."""
    t=z['time'].astype(float);dt=float(np.median(np.diff(t)))
    steps=int(round(lag/dt)) if np.isfinite(lag) else 0
    result=dict(lag_steps=steps,analysis_lag=lag,potential_pairs=max(0,len(t)-steps) if steps else 0,
                eligible_pairs=0,low_coherence_excluded=0,amount_or_centroid_excluded=0,
                empty_hybrid_excluded=0,total_heading_loss=0.,max_complex_closure_error=0.,
                max_heading_closure_error=0.,**{c+'_loss':0. for c in COMPONENTS})
    if steps<1 or steps>=len(t):return result
    ids,A,e,finite=dense(z)
    i0,i1=ids[:-steps],ids[steps:];a0,a1=A[:-steps],A[steps:];e0,e1=e[:-steps],e[steps:]
    # One-to-one track IDs; compare each left slot with each right slot.
    match=(i0[:,:,None]==i1[:,None,:])&(i0[:,:,None]>=0)
    common0=match.any(2);common1=match.any(1)
    right_slot=match.argmax(2)
    ar=np.take_along_axis(a1,right_slot,axis=1);er=np.take_along_axis(e1,right_slot,axis=1)
    ac0=np.where(common0,a0,0);ac1=np.where(common0,ar,0)
    ec0=np.where(common0,e0,0);ec1=np.where(common0,er,0)
    lostA=np.where((i0>=0)&~common0,a0,0);bornA=np.where((i1>=0)&~common1,a1,0)
    lostZ=np.sum(lostA*e0,axis=1);bornZ=np.sum(bornA*e1,axis=1)
    lostW=lostA.sum(1);bornW=bornA.sum(1)
    states={};valid_hybrid=np.ones(len(i0),bool)
    for m,a,b in itertools.product((0,1),repeat=3):
        amounts=ac1 if a else ac0;positions=ec1 if m else ec0
        W=amounts.sum(1)+(bornW if b else lostW)
        Z=(amounts*positions).sum(1)+(bornZ if b else lostZ)
        valid_hybrid&=np.isfinite(W)&(W>1e-12)
        states[(m,a,b)]=np.divide(Z,W,out=np.full(Z.shape,np.nan+1j*np.nan),where=W>1e-12)
    contributions=[]
    for factor in range(3):
        value=np.zeros(len(i0),complex)
        others=[j for j in range(3) if j!=factor]
        for x,y in itertools.product((0,1),repeat=2):
            state=[0,0,0];state[others[0]]=x;state[others[1]]=y
            before=tuple(state);state[factor]=1;after=tuple(state)
            value+=(1/3 if x==y else 1/6)*(states[after]-states[before])
        contributions.append(value)
    P=(z['M1_real'].astype(float)+1j*z['M1_imag'].astype(float))/z['M0'].astype(float)
    p0,p1=P[:-steps],P[steps:];r0,r1=abs(p0),abs(p1)
    coherence=np.isfinite(p0)&np.isfinite(p1)&(r0>=.1)&(r1>=.1)
    amounts=finite[:-steps]&finite[steps:]
    ok=coherence&amounts&valid_hybrid
    result['low_coherence_excluded']=int((~coherence).sum())
    result['amount_or_centroid_excluded']=int((coherence&~amounts).sum())
    result['empty_hybrid_excluded']=int((coherence&amounts&~valid_hybrid).sum())
    result['eligible_pairs']=int(ok.sum())
    if not ok.any():return result
    residual=p1-p0-(states[(1,1,1)]-states[(0,0,0)])
    contributions.append(residual)
    deltaH=p1[ok]/r1[ok]-p0[ok]/r0[ok]
    losses=[np.real(c[ok]*np.conjugate(deltaH))/(r0[ok]+r1[ok]) for c in contributions]
    actual=.5*abs(deltaH)**2
    result['total_heading_loss']=float(actual.sum())
    for c,value in zip(COMPONENTS,losses):result[c+'_loss']=float(value.sum())
    result['max_complex_closure_error']=float(np.max(abs(sum(c[ok] for c in contributions)-(p1-p0)[ok])))
    result['max_heading_closure_error']=float(np.max(abs(sum(losses)-actual)))
    return result

def events(z):
    t=z['time'].astype(float);off=z['region_offset'];ids=z['region_id'];head=z['region_heading'].astype(float)
    widths=z['region_width'].astype(float);amount=z['region_amount'].astype(float);C=z['circumference'].astype(float)
    frame=np.repeat(np.arange(len(t)),np.diff(off));u,first,count=np.unique(ids,return_index=True,return_counts=True)
    ordering=np.argsort(ids,kind='stable');groupoff=np.r_[0,np.cumsum(count)]
    # Every native maximum-overlap identity ends permanently when absent.
    result=[];dominant=z['dominant_track_id'];contrast=z['contrast'].astype(float)
    accepted={};starts=np.r_[0,np.flatnonzero(np.diff(dominant)!=0)+1];ends=np.r_[starts[1:],len(dominant)]
    for start,end in zip(starts,ends):
        identity=int(dominant[start])
        if identity>=0 and end-start>=5:accepted.setdefault(identity,[]).append(int(start))
    for group,(identity,index,n) in enumerate(zip(u,first,count)):
        f=int(frame[index])
        if f==0:continue
        loc=ordering[groupoff[group]:groupoff[group+1]];frames=frame[loc]
        persistent=len(frames)>=5 and np.array_equal(frames[:5],np.arange(f,f+5))
        if not persistent:continue
        old=slice(int(off[f-1]),int(off[f]));oldC=np.full(off[f]-off[f-1],C[f-1])
        certified=certified_no_overlap(head[index],widths[index],C[f],head[old],widths[old],oldC)
        incumbent=int(dominant[f-1]);secondary=incumbent>=0
        oldloc=np.flatnonzero(ids[old]==incumbent)
        distance=halfwidth=math.nan
        if len(oldloc)==1 and np.isfinite(contrast[f-1]) and contrast[f-1]>1e-12:
            j=int(off[f-1])+int(oldloc[0])
            halfwidth=amount[j]/(2*contrast[f-1])
            # Physical distance uses birth-frame actual C; amount/height is a
            # previous-frame physical width (Spring changes negligibly per step).
            distance=abs(wrapped(head[index]-head[j]))*C[f]/(2*np.pi)-halfwidth
        after=[start for start in accepted.get(int(identity),[]) if start>=f]
        takeover=bool(after)
        # Dominance holds count all frames; ambiguous/missing incumbent IDs do
        # not invent successes. Event observation is finite-horizon.
        binid=0 if np.isfinite(distance) and distance<=0 else 1 if np.isfinite(distance) and distance<=C[f]/6 else 2 if np.isfinite(distance) else -1
        result.append(dict(track_id=int(identity),birth_frame=f,birth_time=t[f],persistent_frames=int(n),
            certified_no_parent=certified,secondary=secondary,incumbent_id=incumbent,
            takeover_5frames=takeover,takeover_time=t[int(after[0])+4] if takeover else math.nan,
            distance=distance,distance_over_C=distance/C[f],distance_bin=binid,
            incumbent_halfwidth=halfwidth,followup=t[-1]-t[f],right_censored=bool(frames[-1]==len(t)-1)))
    return result

def worker(index):
    manifest=json.loads((CAMPAIGN/'manifest.json').read_text());condition=manifest['conditions'][index]
    records=[r for r in read_csv(CAMPAIGN/'results/compact_index.csv') if r['condition_id']==condition['condition_id']]
    lags=read_csv(OLD/'supporting/data/analysis_lag_metrics.csv')
    lagrow=next(r for r in lags if r['model']==condition['model'] and float(r['C_label'])==condition['size'] and float(r['f'])==condition['fraction'])
    lag=float(lagrow.get('analysis_lag','nan'))
    seedrows=[];eventrows=[]
    for row in records:
        with np.load(row['compact'],allow_pickle=False) as z:
            ev=events(z);dec=decomposition(z,lag)
            exposure=float(np.sum(.5*(z['circumference'][:-1]+z['circumference'][1:])*np.diff(z['time'])))
            base=dict(condition_id=condition['condition_id'],model=condition['model'],C_label=condition['size'],f=condition['fraction'],replicate=int(row['replicate']))
            secondary=[e for e in ev if e['secondary']]
            cert=[e for e in secondary if e['certified_no_parent']]
            seedrows.append(dict(base,perimeter_time_exposure=exposure,persistent_incident_tracks=len(ev),
                persistent_secondary_onsets=len(secondary),certified_births=len(cert),
                unresolved_secondary_onsets=len(secondary)-len(cert),certified_takeovers=sum(e['takeover_5frames'] for e in cert),**dec))
            eventrows.extend(dict(base,**e) for e in secondary)
    dest=OUT/'supporting/condition_analysis'/f'c{index:03d}'
    write_csv(dest/'seeds.csv',seedrows)
    write_csv(dest/'events.csv',eventrows,fields=list(seedrows[0])[:5]+list(events_schema()))
    (dest/'complete.json').write_text(json.dumps(dict(condition=condition,n_seeds=len(records),events=len(eventrows)),indent=2))
    print(condition['condition_id'],len(records),len(eventrows),flush=True)

def events_schema():
    return ('track_id','birth_frame','birth_time','persistent_frames','certified_no_parent','secondary','incumbent_id','takeover_5frames','takeover_time','distance','distance_over_C','distance_bin','incumbent_halfwidth','followup','right_censored')

def collect():
    manifest=json.loads((CAMPAIGN/'manifest.json').read_text());summaries=[];bins=[];seedall=[];seedbins=[]
    for c in manifest['conditions']:
        dest=OUT/'supporting/condition_analysis'/f"c{c['condition_index']:03d}"
        if not (dest/'complete.json').exists():raise ValueError('Missing '+str(dest))
        seeds=read_csv(dest/'seeds.csv');ev=read_csv(dest/'events.csv');seedall.extend(seeds)
        row={k:seeds[0][k] for k in ('condition_id','model','C_label','f')};row['n_seeds']=len(seeds)
        keys=('perimeter_time_exposure','persistent_incident_tracks','persistent_secondary_onsets','certified_births','unresolved_secondary_onsets','certified_takeovers','potential_pairs','eligible_pairs','low_coherence_excluded','amount_or_centroid_excluded','empty_hybrid_excluded','total_heading_loss')+tuple(c+'_loss' for c in COMPONENTS)
        for key in keys:row[key]=sum(float(s[key]) for s in seeds)
        row['analysis_lag']=float(seeds[0]['analysis_lag'])
        row['certified_birth_rate_lower_bound']=row['certified_births']/row['perimeter_time_exposure']
        row['persistent_onset_rate_upper_bound']=row['persistent_secondary_onsets']/row['perimeter_time_exposure']
        row['certification_fraction']=row['certified_births']/row['persistent_secondary_onsets'] if row['persistent_secondary_onsets'] else math.nan
        row['pair_coverage']=row['eligible_pairs']/row['potential_pairs'] if row['potential_pairs'] else math.nan
        row['max_complex_closure_error']=max(float(s['max_complex_closure_error']) for s in seeds)
        row['max_heading_closure_error']=max(float(s['max_heading_closure_error']) for s in seeds)
        for component in COMPONENTS:row[component+'_fraction']=row[component+'_loss']/row['total_heading_loss'] if row['total_heading_loss']>1e-10 else math.nan
        weights=np.random.default_rng(20261023+c['condition_index']).multinomial(len(seeds),np.full(len(seeds),1/len(seeds)),size=2000)
        births=np.array([float(s['certified_births']) for s in seeds]);exposure=np.array([float(s['perimeter_time_exposure']) for s in seeds])
        lower,upper=np.percentile((weights@births)/(weights@exposure),[2.5,97.5])
        row.update(certified_birth_rate_ci_lower=lower,certified_birth_rate_ci_upper=upper)
        for component in COMPONENTS:
            denominator=weights@np.array([float(s['total_heading_loss']) for s in seeds]);numerator=weights@np.array([float(s[component+'_loss']) for s in seeds])
            frac=np.divide(numerator,denominator,out=np.full(2000,np.nan),where=denominator>1e-10)
            finite=frac[np.isfinite(frac)];ci=np.percentile(finite,[2.5,97.5]) if len(finite) else (math.nan,math.nan)
            row[component+'_fraction_ci_lower'],row[component+'_fraction_ci_upper']=ci
        summaries.append(row)
        for b in range(3):
            members=[e for e in ev if e['certified_no_parent']=='True' and int(e['distance_bin'])==b]
            horizon=100. if row['model']=='spring_ring' else 1000.
            common=[e for e in members if float(e['followup'])>=horizon]
            count=len(common);success=sum(e['takeover_5frames']=='True' and float(e['takeover_time'])-float(e['birth_time'])<=horizon for e in common)
            bins.append(dict({k:row[k] for k in ('condition_id','model','C_label','f')},distance_bin=b,
                births=count,takeovers=success,probability=success/count if count else math.nan,
                common_followup=horizon,all_certified_births=len(members),
                excluded_short_followup=len(members)-count,right_censored_tracks=sum(e['right_censored']=='True' for e in members),
                median_available_followup=float(np.median([float(e['followup']) for e in members])) if members else math.nan,
                cohort='certified_no_overlap_subset',full_cohort_probability_status='unidentifiable_missing_parent_masks'))
            per={int(s['replicate']):[0,0] for s in seeds}
            for e in common:
                per[int(e['replicate'])][0]+=1
                per[int(e['replicate'])][1]+=int(e['takeover_5frames']=='True' and float(e['takeover_time'])-float(e['birth_time'])<=horizon)
            seedbins.extend(dict({k:row[k] for k in ('condition_id','model','C_label','f')},replicate=rep,distance_bin=b,births=v[0],takeovers=v[1]) for rep,v in per.items())
    data=OUT/'supporting/data';write_csv(data/'birth_heading_condition_metrics.csv',summaries)
    write_csv(data/'takeover_by_condition_distance.csv',bins);write_csv(data/'birth_heading_seed_metrics.csv',seedall)
    write_csv(data/'takeover_seed_counts.csv',seedbins)
    # Pool non-drift conditions by distance and size, retaining exact f-specific CSV.
    pooled=[]
    for model in sorted({r['model'] for r in bins}):
        for size in sorted({float(r['C_label']) for r in bins if r['model']==model}):
            for b in range(3):
                relevant=[r for r in bins if r['model']==model and float(r['C_label'])==size and int(r['distance_bin'])==b and float(r['f'])>.01]
                n=sum(r['births'] for r in relevant);k=sum(r['takeovers'] for r in relevant)
                p=k/n if n else math.nan
                relevant_seeds=[s for s in seedbins if s['model']==model and float(s['C_label'])==size and s['distance_bin']==b and float(s['f'])>.01]
                weights=np.random.default_rng(20261023+b).multinomial(len(relevant_seeds),np.full(len(relevant_seeds),1/len(relevant_seeds)),size=2000)
                bootn=weights@np.array([s['births'] for s in relevant_seeds]);bootk=weights@np.array([s['takeovers'] for s in relevant_seeds])
                boot=np.divide(bootk,bootn,out=np.full(2000,np.nan),where=bootn>0);finite=boot[np.isfinite(boot)]
                lo,hi=np.percentile(finite,[2.5,97.5]) if len(finite) else (math.nan,math.nan)
                pooled.append(dict(model=model,C_label=size,distance_bin=b,births=n,takeovers=k,probability=p,
                    ci_lower=lo,ci_upper=hi,pooled_input_rule='f>0.01',common_followup=100. if model=='spring_ring' else 1000.,
                    all_certified_births=sum(r['all_certified_births'] for r in relevant),
                    excluded_short_followup=sum(r['excluded_short_followup'] for r in relevant),
                    right_censored_tracks=sum(r['right_censored_tracks'] for r in relevant),
                    cohort='certified_no_overlap_subset',interval='whole_independent_history_bootstrap_2000'))
    write_csv(data/'takeover_by_size_distance.csv',pooled)
    receipt=dict(conditions=len(summaries),seeds=len(seedall),certified_births=sum(r['certified_births'] for r in summaries),
        persistent_secondary_onsets=sum(r['persistent_secondary_onsets'] for r in summaries),
        max_heading_closure_error=max(r['max_heading_closure_error'] for r in summaries),
        max_complex_closure_error=max(r['max_complex_closure_error'] for r in summaries),
        requested_exact_births='partially_identified_lower_bounds_only',splits='unidentifiable_missing_parent_masks',
        no_noise_simulations=True)
    (OUT/'supporting/metadata').mkdir(parents=True,exist_ok=True)
    (OUT/'supporting/metadata/round5_analysis_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['worker','collect']);parser.add_argument('--index',type=int)
    args=parser.parse_args();worker(args.index) if args.mode=='worker' else collect()
