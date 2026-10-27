"""Per-seed second-order archives and diagnostics for the Poisson ring assay.

This module does not delete data. Exact bootstrap half-crossing brackets and
cumulative integrals accompany a compact lag grid. Frame histories remain
until a separately verified retirement operation is authorized and performed.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.fft import fft,ifft,next_fast_len
from scipy.integrate import cumulative_trapezoid


def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');q.replace(p)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def table(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def correlations(x):
    """C_ij[l]=mean_t x_i[t+l] conj(x_j[t]), with actual pair counts."""
    x=np.asarray(x,np.complex128);N=x.shape[1];L=N//2+1
    f=fft(x,next_fast_len(2*N-1),axis=1)
    pairs=np.arange(N,N-L,-1)
    c=np.empty((len(x),len(x),L),complex)
    for i in range(len(x)):
        for j in range(len(x)):c[i,j]=ifft(f[i]*f[j].conj())[:L]/pairs
    return c


def half(curve,lags):
    if curve[1]<=.5:return None,'below_frame_resolution',math.log(2)/lags[1]
    hits=np.flatnonzero(curve[1:]<=.5)+1
    if not len(hits):return None,'beyond_analysis_lag',math.log(2)/lags[-1]
    k=int(hits[0]);lag=lags[k-1]+(lags[k]-lags[k-1])*(curve[k-1]-.5)/(curve[k-1]-curve[k])
    return float(math.log(2)/lag),'resolved',None


def bootstrap_curves(raw,counts):
    curves=counts@raw
    curves/=curves[:,0,None]
    return curves


def lag_grid(primary,counts):
    """Retain every original bracket used by any four-seed bootstrap draw."""
    L=primary.shape[1]
    keep={0,1,L-1}
    keep.update(np.geomspace(1,L-1,min(512,L-1)).astype(int).tolist())
    curves=bootstrap_curves(primary,counts)
    for curve in curves:
        for threshold in (.5,math.exp(-1),0.):
            hits=np.flatnonzero(curve[1:]<=threshold)+1
            if len(hits):
                k=int(hits[0]);keep.update((k-1,k))
    return np.array(sorted(keep),int)


def centered(c,x):
    N=x.shape[1];L=c.shape[-1];ell=np.arange(L);pairs=N-ell
    prefix=np.column_stack([np.zeros(len(x),complex),np.cumsum(x,axis=1)])
    late=(prefix[:,-1,None]-prefix[:,ell])/pairs
    early=prefix[:,N-ell]/pairs
    mu=x.mean(1)
    return c-late[:,None,:]*mu.conj()[None,:,None]-mu[:,None,None]*early.conj()[None,:,:]+mu[:,None,None]*mu.conj()[None,:,None]


def summaries(raw,indices,cadence,counts,weights):
    dense_lags=np.arange(raw.shape[1])*cadence
    grid_lags=indices*cadence
    full=bootstrap_curves(raw,counts);thin=bootstrap_curves(raw[:,indices],counts)
    records=[]
    for a,b in zip(full,thin):
        x=half(a,dense_lags);y=half(b,grid_lags)
        assert x[1]==y[1]
        assert np.isclose(x[0] if x[0] is not None else x[2],y[0] if y[0] is not None else y[2],rtol=1e-9,atol=1e-14)
        records.append(x)
    vals=np.array([r[0] if r[0] is not None else r[2] for r in records])
    kinds=np.array(['resolved' if r[1]=='resolved' else 'upper_bound' if r[1]=='beyond_analysis_lag' else 'lower_bound' for r in records])
    order=np.argsort(vals,kind='stable');cdf=np.cumsum(weights[order])
    ends=[order[min(np.searchsorted(cdf,q),len(order)-1)] for q in (.025,.975)]
    pooled=half(raw.mean(0)/raw[:,0].mean(),dense_lags)
    return dict(gamma=pooled[0],status=pooled[1],bound=pooled[2],
                gamma_ci_low=float(vals[ends[0]]),gamma_ci_high=float(vals[ends[1]]),
                ci_low_status=str(kinds[ends[0]]),ci_high_status=str(kinds[ends[1]])),vals,kinds


def load(O,slot,mode,species):return O.series(slot,mode,species)


def archive_condition(O,ci):
    m=O.verify();c=m['conditions'][ci];j=m['jobs'][c['slots'][0]];counts,weights=O.bootstrap_counts()
    levels=read(O.RUN/'expected_levels.json');level=int(levels[str(c['slots'][0])])
    dest=O.RUN/'compact_statistics'/f'{ci:03d}';dest.mkdir(parents=True,exist_ok=True)
    seed_arrays=[{} for _ in c['slots']];summary=[];transport=[];coherence=[]
    N=len(load(O,c['slots'][0],0,'u'));windows={'full':(0,N),'first_half':(0,N//2),'second_half':(N//2,N)}
    spring=c['model']=='spring_ring';E=O.engine();cad=j['frame_interval']
    basal=next(b for b in m['basals'] if b['model']==j['model'] and b['size']==j['size'])
    with np.load(O.RUN/basal['input']) as z:state=z['state']
    cells=(len(state)-1)//2 if spring else state.shape[1]
    if spring:diff=(.01,10.)
    else:diff=E.native('chemical').ring.native_model(j['model']).diffusivities[:2]
    circles=[]
    for slot in c['slots']:
        circles.append(np.concatenate([np.load(O.segment(slot,k)/'circumference.npy') for k in range(level+1)]) if spring else np.full(N,j['size']))
    # Preserve heading strength and increment distributions before losing time order.
    for si,slot in enumerate(c['slots']):
        u0=load(O,slot,0,'u').real;u1=load(O,slot,1,'u')
        P=np.divide(u1.conj(),u0,out=np.full(N,np.nan+0j),where=u0>0)
        for window,(start,stop) in windows.items():
            p=P[start:stop];strength=abs(p);valid=np.isfinite(p)&(strength>0)
            h=np.zeros(len(p),complex);h[valid]=p[valid]/strength[valid]
            hc=correlations(h[None,:])[0,0]
            pairs=np.arange(len(h),len(h)-len(hc),-1)
            good_pairs=correlations(valid[None,:].astype(float))[0,0].real*pairs
            nums=hc*pairs
            idx=np.unique(np.r_[0,1,np.geomspace(1,len(hc)-1,min(512,len(hc)-1)).astype(int)])
            prefix=window+'_heading'
            seed_arrays[si][prefix+'_lag']=idx*cad
            seed_arrays[si][prefix+'_sum']=nums[idx]
            seed_arrays[si][prefix+'_pairs']=np.rint(good_pairs[idx]).astype('int64')
            seed_arrays[si][prefix+'_mean']=np.mean(h[valid]) if valid.any() else np.nan+0j
            bins=np.linspace(-math.pi,math.pi,129);hist=[];eligible_hist=[];pairs_out=[];large=[]
            for lag in (j['tau']/4,j['tau'],4*j['tau'],16*j['tau']):
                k=int(round(lag/cad));ok=valid[k:]&valid[:-k]
                theta=np.angle(h[k:]*h[:-k].conj());high=ok&(strength[k:]>=.1)&(strength[:-k]>=.1)
                hist.append(np.histogram(theta[ok],bins)[0]);eligible_hist.append(np.histogram(theta[high],bins)[0])
                pairs_out.append([int(ok.sum()),int(high.sum())]);large.append([int(np.count_nonzero(abs(theta[ok])>math.pi/2)),int(np.count_nonzero(abs(theta[high])>math.pi/2))])
            seed_arrays[si][prefix+'_increment_edges']=bins
            seed_arrays[si][prefix+'_increment_hist']=np.array(hist)
            seed_arrays[si][prefix+'_increment_high_coherence_hist']=np.array(eligible_hist)
            seed_arrays[si][prefix+'_increment_pairs']=np.array(pairs_out)
            seed_arrays[si][prefix+'_large_increment_counts']=np.array(large)
            seed_arrays[si][prefix+'_increment_lags']=np.array([j['tau']/4,j['tau'],4*j['tau'],16*j['tau']])
            finite=strength[np.isfinite(strength)];q=np.quantile(finite,[0,.01,.05,.25,.5,.75,.95,.99,1])
            seed_arrays[si][window+'_coherence_hist']=np.histogram(finite,np.linspace(0,1,129))[0]
            coherence.append(dict(condition=ci,slot=slot,window=window,frames=len(p),valid_headings=int(valid.sum()),
                                  mean=float(finite.mean()),q0=float(q[0]),q01=float(q[1]),q05=float(q[2]),q25=float(q[3]),
                                  median=float(q[4]),q75=float(q[5]),q95=float(q[6]),q99=float(q[7]),q100=float(q[8]),
                                  fraction_below_005=float(np.mean(finite<.05)),fraction_below_01=float(np.mean(finite<.1))))
    for n in range(j['maximum_mode']+1):
        uv=[np.stack([load(O,s,n,'u'),load(O,s,n,'v')]) for s in c['slots']]
        for window,(start,stop) in windows.items():
            x=[a[:,start:stop] for a in uv];corr=np.array([correlations(a) for a in x])
            raw=corr[:,0,0].real;idx=lag_grid(raw,counts)
            stats,boot,bkind=summaries(raw,idx,cad,counts,weights)
            summary.append(dict(condition=ci,model=c['model'],size=c['size'],rate_level=c['rate_level'],amplitude_level=c['amplitude_level'],window=window,mode=n,frames=stop-start,analysis_lag=(corr.shape[-1]-1)*cad,**stats))
            for si,slot in enumerate(c['slots']):
                key=f'{window}_n{n:02d}';data=seed_arrays[si];C=circles[si][start:stop]
                data[key+'_lag']=idx*cad;data[key+'_pairs']=stop-start-idx
                data[key+'_C']=corr[si][:,:,idx]
                data[key+'_connected_C']=centered(corr[si],x[si])[:,:,idx]
                data[key+'_mean']=x[si].mean(1)
                data[key+'_powers']=np.array([corr[si,0,0,0].real,corr[si,1,1,0].real,corr[si,:,:,0].sum().real])
                raw3=np.stack([corr[si,0,0].real,corr[si,1,1].real,corr[si].sum((0,1)).real])
                data[key+'_raw_cumulative_integral']=cumulative_trapezoid(raw3,dx=cad,axis=1,initial=0)[:,idx]
                data[key+'_primary_bootstrap_rates']=boot;data[key+'_primary_bootstrap_status']=bkind.astype('U16')
                kappa=4*(cells/C)**2*math.sin(math.pi*n/cells)**2
                M=C*(x[si][0]+x[si][1]);flux=-C*kappa*(diff[0]*x[si][0]+diff[1]*x[si][1])
                delta=np.diff(M);pred=.5*cad*(flux[:-1]+flux[1:]);res=delta-pred
                transport.append(dict(condition=ci,slot=slot,window=window,mode=n,
                                      residual_rms=float(np.sqrt(np.mean(abs(res)**2))),increment_rms=float(np.sqrt(np.mean(abs(delta)**2))),
                                      predicted_increment_rms=float(np.sqrt(np.mean(abs(pred)**2))),residual_max=float(abs(res).max()),
                                      integrated_residual_abs=float(abs(res.sum())),
                                      mode0_relative_mass_error=float(np.max(abs(M/M[0]-1))) if n==0 else None))
    for si,slot in enumerate(c['slots']):
        seed_arrays[si]['bootstrap_counts']=counts.astype('int8');seed_arrays[si]['bootstrap_weights']=weights
        seed_arrays[si]['coherence_hist_edges']=np.linspace(0,1,129)
        seed_arrays[si]['mean_circumference']=np.array(circles[si].mean())
        np.savez_compressed(dest/f'seed_{slot:03d}.npz',**seed_arrays[si])
    table(dest/'rates_by_window.csv',summary);table(dest/'transport_balance.csv',transport);table(dest/'coherence.csv',coherence)
    # Compare the archived reconstruction with the original independently produced rates.
    with (O.RUN/'analysis'/f'{ci:03d}'/'rates.csv').open() as f:original={int(r['mode']):r for r in csv.DictReader(f)}
    checked=0
    for r in summary:
        if r['window']!='full':continue
        old=original[r['mode']];assert r['status']==old['status']
        for name in ('gamma','bound','gamma_ci_low','gamma_ci_high'):
            a=r[name];b=float(old[name]) if old[name] else None
            assert (a is None and b is None) or (a is not None and b is not None and np.isclose(a,b,rtol=2e-8,atol=1e-13)),(ci,r['mode'],name,a,b)
        assert r['ci_low_status']==old['ci_low_status'] and r['ci_high_status']==old['ci_high_status'];checked+=1
    # Read archives back and independently reconstruct every primary bootstrap curve.
    with np.load(dest/f"seed_{c['slots'][0]:03d}.npz") as z:keys=[k[:-2] for k in z.files if k.endswith('_C') and not k.endswith('connected_C')]
    for key in keys:
        arrays=[];lags=None
        for slot in c['slots']:
            with np.load(dest/f'seed_{slot:03d}.npz') as z:
                arrays.append(z[key+'_C'][0,0].real.copy());lags=z[key+'_lag'].copy()
        curves=bootstrap_curves(np.array(arrays),counts)
        with np.load(dest/f"seed_{c['slots'][0]:03d}.npz") as z:
            expected=z[key+'_primary_bootstrap_rates'];kinds=z[key+'_primary_bootstrap_status']
        for k,curve in enumerate(curves):
            result=half(curve,lags);value=result[0] if result[0] is not None else result[2]
            assert np.isclose(value,expected[k],rtol=1e-9,atol=1e-14)
    sources={}
    for slot in c['slots']:
        for k in range(level+1):
            for p in sorted(O.segment(slot,k).glob('*.npy')):
                sources[str(p.relative_to(O.RUN))]=dict(sha256=sha(p),logical_bytes=p.stat().st_size,allocated_bytes=p.stat().st_blocks*512)
    outputs={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(dest.iterdir()) if p.suffix in ('.npz','.csv')}
    write(dest/'validation.json',dict(generated_utc=now(),passed=True,condition=ci,primary_modes_verified=checked,
          all35_bootstrap_crossings_verified_all_windows=True,seed_ids=c['slots'],extension_level=level,
          source_histories=sources,output_files=outputs,raw_histories_deleted=False,
          limits='Retains second-order statistics on a declared lag grid and selected heading increments; arbitrary new trajectory event detectors cannot be recovered.'))


def check_archive(O):
    m=O.verify();records=[];levels=read(O.RUN/'expected_levels.json')
    for c in m['conditions']:
        folder=O.RUN/'compact_statistics'/f"{c['index']:03d}";v=read(folder/'validation.json')
        assert v['passed']
        assert all(v['extension_level']==int(levels[str(s)]) for s in c['slots'])
        for name,d in v['output_files'].items():assert sha(folder/name)==d['sha256']
        records.append(dict(condition=c['index'],modes=v['primary_modes_verified'],files=len(v['output_files'])))
    write(O.RUN/'results/compact_archive_validation.json',dict(generated_utc=now(),passed=True,conditions=len(records),
          records=records,raw_histories_deleted=False,retirement_ready=True))
    return records


def checks():
    rng=np.random.default_rng(4);x=rng.normal(size=(2,19))+1j*rng.normal(size=(2,19));c=correlations(x)
    d=centered(c,x);z=x-x.mean(1)[:,None]
    for i in range(2):
        for j in range(2):
            direct=np.array([np.mean(x[i,k:]*x[j,:19-k].conj()) for k in range(c.shape[-1])])
            conn=np.array([np.mean(z[i,k:]*z[j,:19-k].conj()) for k in range(c.shape[-1])])
            assert np.allclose(c[i,j],direct,atol=1e-12) and np.allclose(d[i,j],conn,atol=1e-12)
    counts=np.array([a for a in __import__('itertools').product(range(5),repeat=4) if sum(a)==4])
    raw=np.exp(-np.arange(100)[None,:]/np.array([5.,8.,20.,80.])[:,None]);idx=lag_grid(raw,counts)
    weights=np.array([math.factorial(4)/math.prod(math.factorial(int(v)) for v in a)/4**4 for a in counts])
    summaries(raw,idx,.5,counts,weights)
    print(json.dumps(dict(complex_cross_correlations='pass',centered_correlations='pass',bootstrap_bracket_retention='pass')))


if __name__=='__main__':checks()
