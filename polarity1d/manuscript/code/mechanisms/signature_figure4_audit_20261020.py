"""Independent direct-sum checks of E0 native scalar estimands.

Uses direct correlation rather than the producer FFT, and explicit contiguous
increments rather than its coordinate/FFT-MSD functions. No producer numerical
helper is imported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime, timezone
import warnings
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for b in iter(lambda:stream.read(4<<20),b''):h.update(b)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def declared_audit_sources(root,manifest):
    declared=dict(manifest['source_sha256'])
    legacy=root/'audit_freeze.json'
    if legacy.is_file():
        record=read(legacy);sources=record.get('source_sha256')
        if isinstance(sources,dict):declared.update(sources)
        elif record.get('source_path') and isinstance(sources,str):declared[record['source_path']]=sources
    return declared


def close(actual,expected,label):
    if not np.allclose(actual,expected,rtol=2e-8,atol=2e-10,equal_nan=True):
        raise AssertionError(label)


def divide(a,b):
    return np.divide(a,b,out=np.full(np.broadcast_shapes(np.shape(a),np.shape(b)),np.nan),where=np.asarray(b)>0)


def coordinates(P,C,valid):
    n=len(P);theta=np.zeros(n);arc=np.zeros(n)
    for t in range(1,n):
        if valid[t] and valid[t-1]:
            increment=float(np.angle(P[t]*P[t-1].conjugate()))
            theta[t]=theta[t-1]+increment
            arc[t]=arc[t-1]+increment*(C[t]+C[t-1])/(4*np.pi)
    return theta,arc,theta*C/(2*np.pi)


def direct_msd(coordinate,valid,lags):
    prefix=np.r_[0,np.cumsum(~valid)]
    numerators=[];counts=[]
    for lag in lags:
        eligible=(prefix[lag+1:]-prefix[:-lag-1])==0
        numerators.append(float(np.sum((coordinate[lag:]-coordinate[:-lag])[eligible]**2)))
        counts.append(int(np.sum(eligible)))
    return np.array(numerators),np.array(counts)


def intervals(curves,dt):
    lower=[];upper=[];codes=[]
    for curve in curves:
        hits=np.flatnonzero(np.isfinite(curve)&(curve<=math.exp(-1)))
        if len(hits):
            i=int(hits[0]);prior=1. if i==0 else curve[i-1]
            if not np.isfinite(curve[:i+1]).all():lo,hi,code=0.,math.inf,3
            else:
                point=i*dt+(prior-math.exp(-1))/(prior-curve[i])*dt
                lo,hi,code=(point,point,0) if point>=2*dt else (0.,2*dt,1)
        elif np.isfinite(curve).all():lo,hi,code=len(curve)*dt,math.inf,2
        else:lo,hi,code=0.,math.inf,3
        lower.append(lo);upper.append(hi);codes.append(code)
    return np.array(lower),np.array(upper),np.array(codes)


def audit(root,task):
    manifest=read(root/'analysis_manifest.json');group=manifest['groups'][task]
    source=str(Path(__file__).resolve())
    if declared_audit_sources(root,manifest).get(source)!=sha(source):raise ValueError('The actual independent audit must belong to the frozen source closure')
    w=root/f'analysis_workers/{task:04d}';result=read(w/'result.json')
    if sha(w/'arrays.npz')!=result['arrays_sha256']:raise ValueError('Worker arrays changed')
    lookup={(r['id'],r['view']):r for r in result['case_metrics']}
    source=[]
    for row in group['noise']+[group['mean'],group['zero']]:
        if sha(row['scalar_path'])!=row['scalar_sha256']:raise ValueError('Native scalar changed')
        with np.load(row['scalar_path'],allow_pickle=False) as z:
            source.append((row,{k:z[k] for k in ('time','P_real','P_imag','circumference','contrast','current_half_contrast_regions')}))
    key=group.get('bootstrap_pairing_key',f"{group['model']}:{group['C_label']}:{group['f']}")
    seed=2026102004+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)
    n=len(group['noise']);weights=np.random.default_rng(seed).multinomial(n,np.full(n,1/n),size=2000)
    checked_lags=0;checked_pairs=0
    with np.load(w/'arrays.npz',allow_pickle=False) as arrays:
        for view,(start,stop,step) in manifest['views'].items():
            nums=[];counts=[];angular_nums=[];angular_counts=[];site_nums={x:[] for x in ('heading_arc','heading_direct')};site_counts={x:[] for x in site_nums}
            fractions=[];Cmeans=[];mean_forced=None
            for index,(row,z) in enumerate(source):
                ix=np.flatnonzero((z['time']>=start*row['native_horizon']-1e-8)&(z['time']<=stop*row['native_horizon']+1e-8))[::step]
                P=(z['P_real']+1j*z['P_imag'])[ix];C=np.broadcast_to(z['circumference'],z['time'].shape)[ix]
                contrast=z['contrast'][ix];raw=z['current_half_contrast_regions'][ix];valid=abs(P)>=.1
                low=contrast<.1*row['basal_contrast'];single=(~low)&(raw==1);multi=(~low)&(raw>=2)
                frac=np.array([single.mean(),multi.mean(),low.mean(),(~(single|multi|low)).mean()])
                native=lookup[(row['id'],view)]
                for j,name in enumerate(('single_front','multiple_fronts','low_contrast','unclassified')):close(native[name+'_fraction'],frac[j],'native regime '+name)
                close(native['R_valid_fraction'],valid.mean(),'native R coverage')
                close(native['actual_C_mean'],C.mean(),'actual circumference')
                forced=np.mean((~valid)|low)
                close(native['forced_union_fraction'],forced,'matched mean/zero forcing union')
                if index==n:mean_forced=forced>.5
                dt=float(z['time'][ix[1]]-z['time'][ix[0]]);lag=(len(ix)-1)//2
                unit=np.zeros(len(P),complex);unit[valid]=P[valid]/abs(P[valid])
                direct_num=np.correlate(unit,unit,'full')[len(P):len(P)+lag].real
                mask=valid.astype(np.int64);direct_count=np.correlate(mask,mask,'full')[len(P):len(P)+lag]
                if 'per_case_tau_lower' in native:
                    lo,hi,_=intervals(divide(direct_num,direct_count)[None,:],dt)
                    close(native['per_case_tau_lower'],lo[0],'Individual/control persistence lower')
                    close(math.inf if native['per_case_tau_upper'] is None else native['per_case_tau_upper'],hi[0],'Individual/control persistence upper')
                if index>=n:continue
                nums.append(direct_num);counts.append(direct_count)
                chosen=np.array(sorted(set(range(1,11))|{lag}),int)
                theta,_,_=coordinates(P,C,valid)
                an,ac=direct_msd(theta,valid,chosen);angular_nums.append(an);angular_counts.append(ac)
                sitevalid=valid&single
                _,arc,direct=coordinates(P,C,sitevalid)
                for method,coord in (('heading_arc',arc),('heading_direct',direct)):
                    sn,sc=direct_msd(coord,sitevalid,chosen);site_nums[method].append(sn);site_counts[method].append(sc)
                fractions.append(frac);Cmeans.append(C.mean())
                checked_pairs+=len(P);checked_lags+=lag
            pre=view+'__';nums=np.array(nums);counts=np.array(counts);possible=len(P)-np.arange(1,nums.shape[1]+1)
            curve=divide(nums.sum(axis=0),counts.sum(axis=0));coverage=divide(counts.sum(axis=0),n*possible)
            close(arrays[pre+'C'],curve,'Every heading-cosine lag, direct correlation')
            close(arrays[pre+'eligible_pairs'],counts.sum(axis=0),'Every eligible heading pair count')
            close(arrays[pre+'coverage'],coverage,'Every endpoint coverage lag')
            boot=divide(weights@nums,weights@counts)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning);ci=np.nanquantile(boot,[.025,.975],axis=0)
            close(arrays[pre+'C_ci_lower'],ci[0],'Whole-seed lower cosine CI')
            close(arrays[pre+'C_ci_upper'],ci[1],'Whole-seed upper cosine CI')
            blo,bhi,bc=intervals(boot,dt)
            close(arrays[pre+'tau_boot_lower'],blo,'Bootstrap persistence lower intervals')
            close(arrays[pre+'tau_boot_upper'],bhi,'Bootstrap persistence upper intervals')
            close(arrays[pre+'tau_boot_code'],bc,'Bootstrap persistence statuses')
            an=np.array(angular_nums);ac=np.array(angular_counts)
            msd=divide(an.sum(axis=0),ac.sum(axis=0))
            close(arrays[pre+'angular_MSD'][chosen-1],msd,'Contiguous angular MSD direct displacement')
            close(arrays[pre+'angular_MSD_coverage'][chosen-1],divide(ac.sum(axis=0),n*possible[chosen-1]),'Continuous-track coverage')
            x=np.arange(2,11)*dt;x=x-x.mean();db=(divide(weights@an,weights@ac)[:,1:10]@x)/(2*float(x@x))
            close(arrays[pre+'Dtheta_boot'],db,'Whole-seed short-lag signed mobility')
            cb=weights@np.array(Cmeans)/n
            close(arrays[pre+'Cmean_boot'],cb,'Whole-seed actual circumference')
            close(arrays[pre+'Dx_boot'],db*(cb/(2*np.pi))**2,'Actual-geometry physical conversion')
            for method in site_nums:
                sn=np.array(site_nums[method]);sc=np.array(site_counts[method]);sp=pre+method+'__'
                close(arrays[sp+'MSD'][chosen-1],divide(sn.sum(axis=0),sc.sum(axis=0)),method+' direct MSD')
                close(arrays[sp+'coverage'][chosen-1],divide(sc.sum(axis=0),n*possible[chosen-1]),method+' single-qualified-track coverage')
            metric=result['metrics'][view]
            if metric['forced']!=mean_forced:raise AssertionError('Matched-mean forced classification')
            close(arrays[pre+'regime_boot'],weights@np.array(fractions)/n,'Whole-seed regime fractions')
    out=root/'estimand_audit_workers';out.mkdir(exist_ok=True)
    receipt=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,task=task,
        analysis_manifest_sha256=sha(root/'analysis_manifest.json'),producer_result_sha256=sha(w/'result.json'),
        audit_source_sha256=sha(__file__),noisy_realizations=n,views=7,all_heading_lags_checked=checked_lags,
        retained_scalar_frames_rechecked=checked_pairs,method='Independent direct correlation, explicit contiguous increments, native mean-control union; all 2000 bootstrap intervals rechecked')
    write(out/f'{task:04d}.json',receipt)
    return receipt


def collect(root):
    m=read(root/'analysis_manifest.json');rows=[]
    for group in m['groups']:
        p=root/f"estimand_audit_workers/{group['task']:04d}.json";r=read(p)
        if not r['passed'] or r['analysis_manifest_sha256']!=sha(root/'analysis_manifest.json'):raise ValueError('Incomplete independent audit')
        if r['audit_source_sha256'] not in declared_audit_sources(root,m).values():raise ValueError('Independent audit used an undeclared executable source')
        if r['producer_result_sha256']!=sha(root/f"analysis_workers/{group['task']:04d}/result.json"):raise ValueError('Producer output changed after independent audit')
        rows.append(dict(task=group['task'],receipt_sha256=sha(p),audit_source_sha256=r['audit_source_sha256']))
    result=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,complete=True,condition_groups=len(rows),
        noisy_realizations=sum(len(g['noise']) for g in m['groups']),all_seven_views=True,bootstrap_resamples=2000,
        analysis_manifest_sha256=sha(root/'analysis_manifest.json'),inventory_sha256=sha(root/'inventory.json'),receipts=rows)
    write(root/'results/analysis/independent_estimand_audit.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=('worker','collect'));p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=int)
    args=p.parse_args();result=audit(args.out.resolve(),args.task) if args.command=='worker' else collect(args.out.resolve())
    print(json.dumps(result,indent=2,allow_nan=False))
