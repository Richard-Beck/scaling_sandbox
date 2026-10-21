"""Independent direct native-field checks of additive literal A3 MSDs.

No producer coordinate, FFT or aggregation functions are imported. Every
native retained sample supplies its own qualification mask and literal
coordinate. Pair counts are checked at every lag; selected short/interior/last
lag MSDs and 2,000 whole-seed intervals are checked by direct subtraction.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import warnings
import numpy as np


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(4<<20),b''):h.update(b)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def close(a,b,label):
    if not np.allclose(a,b,rtol=2e-8,atol=2e-10,equal_nan=True):raise AssertionError(label)


def route(signal,C,valid,policy='split',tolerance=1e-12):
    x=np.full(len(C),np.nan);angle=0.;cuts=np.zeros(len(C),bool);ties=np.zeros(len(C),bool)
    for i in range(len(C)):
        if not valid[i]:continue
        phase=math.atan2(signal[i].imag,signal[i].real)
        if i and valid[i-1]:
            previous=math.atan2(signal[i-1].imag,signal[i-1].real);step=(phase-previous+math.pi)%(2*math.pi)-math.pi
            ties[i]=abs(abs(step)-math.pi)<=tolerance;cuts[i]=ties[i] and policy=='split'
        if i==0 or not valid[i-1] or cuts[i]:angle=phase%(2*math.pi)
        else:
            if ties[i]:step=math.pi if policy=='plus_pi' else -math.pi
            angle+=step
        x[i]=angle*C[i]/(2*math.pi)
    return x,cuts,ties


def coordinate(signal,C,valid,policy='split'):
    return route(signal,C,valid,policy)[0]


def tracks(row):
    for k in ('scalar','profile','basal_profile','source_status'):
        if digest(row[k+'_path'])!=row[k+'_sha256']:raise ValueError('Native input hash changed')
    with np.load(row['scalar_path'],allow_pickle=False) as z:
        t=z['time'];M=z['M1_real']+1j*z['M1_imag'];R=np.hypot(z['P_real'],z['P_imag']);C=np.broadcast_to(z['circumference'],t.shape).copy()
        contrast=z['contrast'];regions=z['current_half_contrast_regions'];v=(R>=.1)&(regions==1)&(contrast>=.1*row['basal_contrast'])
    with np.load(row['profile_path'],allow_pickle=False) as z:
        pt=z['time'];pc=np.broadcast_to(z['circumference'],pt.shape).copy();field=z['primary'];q=z['material_positions'] if 'material_positions' in z else z['x']/pc[0]
        peak=np.mod(q[field.argmax(axis=1)],1.);contrast=field.max(axis=1)-field.min(axis=1)
    ix=np.searchsorted(t,pt);close(t[ix],pt,'Exact native profile alignment')
    pv=(R[ix]>=.1)&(regions[ix]==1)&(contrast>=.1*row['basal_contrast'])
    return {'heading_literal_actualC':(t,M,C,v),'profile_peak_literal_actualC':(pt,np.cos(2*np.pi*peak)+1j*np.sin(2*np.pi*peak),pc,pv)}


def track_labels(valid,cuts=None):
    labels=np.full(len(valid),-1,dtype=np.int64);current=-1
    for i,good in enumerate(valid):
        if not good:continue
        if i==0 or not valid[i-1] or (cuts is not None and cuts[i]):current+=1
        labels[i]=current
    return labels


def direct(values,valid,lags,cuts=None):
    labels=track_labels(valid,cuts);nums=[];counts=[]
    for lag in lags:
        good=(labels[lag:]==labels[:-lag])&(labels[lag:]>=0)
        delta=values[lag:]-values[:-lag]
        nums.append(np.square(delta[good]).sum());counts.append(good.sum())
    return np.asarray(nums),np.asarray(counts)


def all_counts(valid,maxlag,cuts=None):
    labels=track_labels(valid,cuts);lengths=np.bincount(labels[labels>=0])
    histogram=np.bincount(lengths,minlength=len(valid)+1);sizes=np.arange(len(histogram))
    tail_count=np.cumsum(histogram[::-1])[::-1];tail_length=np.cumsum((histogram*sizes)[::-1])[::-1]
    lag=np.arange(1,maxlag+1)
    return tail_length[lag+1]-lag*tail_count[lag+1]


def divide(a,b):return np.divide(a,b,out=np.full(np.broadcast_shapes(np.shape(a),np.shape(b)),np.nan),where=np.asarray(b)>0)


def audit(root,task):
    m=read(root/'manifest.json')
    if str(Path(__file__).resolve()) not in m['source_sha256']:raise ValueError('Execute the actual frozen independent site auditor')
    for p,h in {**m['source_sha256'],**m['input_receipts']}.items():
        if digest(p)!=h:raise ValueError('Frozen additive site closure changed')
    g=m['groups'][task];folder=root/f'workers/{task:04d}';r=read(folder/'result.json')
    if not r['passed'] or r['task']!=task or r['manifest_sha256']!=digest(root/'manifest.json') or digest(folder/'arrays.npz')!=r['arrays_sha256']:raise ValueError('Literal site worker identity changed')
    key=g.get('bootstrap_pairing_key',f"{g['model']}:{g['C_label']}:{g['f']}");seed=2026102004+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)
    n=len(g['noise']);w=np.random.default_rng(seed).multinomial(n,np.full(n,1./n),size=2000)
    per={view:{method:{policy:[] for policy in m['antipodal_policies']} for method in m['methods']} for view in m['views']};frames=0;coordinate_samples=0
    recorded={(x['id'],x['view'],x['method']):x for x in r['case_ambiguities']}
    for row in g['noise']:
        source=tracks(row)
        for view,(lo,hi,step) in m['views'].items():
            for method,(t,signal,C,valid) in source.items():
                ix=np.flatnonzero((t>=lo*row['native_horizon']-1e-8)&(t<=hi*row['native_horizon']+1e-8))[::step]
                mask=valid[ix];lag=(len(ix)-1)//2
                chosen=np.array(sorted(set(range(1,min(10,lag)+1))|{max(1,lag//4),max(1,lag//2),max(1,3*lag//4),lag}),int)
                _,split,ties=route(signal[ix],C[ix],mask,'split',m['antipodal_atol_radians']);base=all_counts(mask,lag);unambiguous=all_counts(mask,lag,split)
                declared=recorded[(row['id'],view,method)]
                for key,actual in dict(qualified_frames=int(mask.sum()),total_frames=len(ix),qualified_adjacent_steps=int(np.sum(mask[1:]&mask[:-1])),antipodal_steps=int(ties.sum()),tracks_before_split=len(set(track_labels(mask)[mask])),tracks_after_split=len(set(track_labels(mask,split)[mask]))).items():
                    if declared[key]!=actual:raise AssertionError('Native antipodal/qualification count '+key)
                for policy in m['antipodal_policies']:
                    x,cuts,_=route(signal[ix],C[ix],mask,policy,m['antipodal_atol_radians']);nums,counts=direct(x,mask,chosen,cuts);full=all_counts(mask,lag,cuts)
                    close(counts,full[chosen-1],'Direct eligibility counts')
                    per[view][method][policy].append((nums,counts,full,len(ix)-np.arange(1,lag+1),chosen,float(t[ix[1]]-t[ix[0]]),base,unambiguous))
                frames+=len(ix);coordinate_samples+=int(mask.sum())
    checks=[]
    with np.load(folder/'arrays.npz',allow_pickle=False) as z:
        for view in m['views']:
            for method in m['methods']:
                for policy in m['antipodal_policies']:
                    rows=per[view][method][policy];chosen=rows[0][4];nums=np.array([r[0] for r in rows]);counts=np.array([r[1] for r in rows]);full=np.sum([r[2] for r in rows],axis=0);possible=np.sum([r[3] for r in rows],axis=0);base=np.sum([r[6] for r in rows],axis=0);unambiguous=np.sum([r[7] for r in rows],axis=0)
                    prefix=view+'__'+method+'__'+policy+'__';close(z[prefix+'lag'],np.arange(1,len(full)+1)*rows[0][5],'All actual lag values')
                    close(z[prefix+'eligible_pairs'],full,'Every contiguous-section pair count');close(z[prefix+'possible_pairs'],possible,'Every possible pair count');close(z[prefix+'coverage'],divide(full,possible),'Every lag coverage')
                    for name,value in dict(eligible_pairs_without_ambiguity_split=base,coverage_without_ambiguity_split=divide(base,possible),unambiguous_pairs=unambiguous,pairs_lost_to_ambiguity_split=base-unambiguous,coverage_lost_to_ambiguity_split=divide(base-unambiguous,possible),route_assumed_pairs=full-unambiguous).items():close(z[prefix+name],value,'Ambiguity coverage '+name)
                    close(z[prefix+'MSD'][chosen-1],divide(nums.sum(axis=0),counts.sum(axis=0)),'Direct native coordinate MSD')
                    bootstrap=divide(w@nums,w@counts)
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore',RuntimeWarning);ci=np.nanquantile(bootstrap,[.025,.975],axis=0)
                    close(z[prefix+'MSD_ci_lower'][chosen-1],ci[0],'Whole-seed lower MSD interval');close(z[prefix+'MSD_ci_upper'][chosen-1],ci[1],'Whole-seed upper MSD interval')
                    checks.append(dict(view=view,method=method,antipodal_policy=policy,every_lag_pair_coverage=True,direct_subtraction_lags=chosen.tolist()))
    out=root/f'audits/{task:04d}.json'
    if out.exists():raise FileExistsError('Independent site receipts are immutable')
    out.write_text(json.dumps(dict(passed=True,task=task,source_version=2,manifest_sha256=digest(root/'manifest.json'),result_sha256=digest(folder/'result.json'),arrays_sha256=digest(folder/'arrays.npz'),native_frames=frames,qualified_coordinate_samples=coordinate_samples,bootstrap_resamples=2000,whole_seeds=n,checks=checks,all_qualified_native_frames_preserved=True,scope='Independent native coordinates/masks and antipodal counts at all saved samples; every-lag split/route pair coverage; direct subtraction and bootstrap at all declared short/interior/last audit lags'),indent=2)+'\n')
    return dict(passed=True,task=task,native_frames=frames)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=int,required=True);args=p.parse_args();print(json.dumps(audit(args.out.resolve(),args.task),indent=2))
