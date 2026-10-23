"""Detector-free Fourier evidence, inherited lag grid and whole-seed bootstrap."""
from common import *
import numpy as np
import figure4_round3_estimators_20261022 as E

def fourier(p):
    n=p.shape[1];theta=2*np.pi*(np.arange(n)+.5)/n
    return p@np.exp(1j*np.outer(theta,np.arange(min(33,n//2+1))))/n

def sufficient(t,a,model,Pabs=None):
    t=np.asarray(t,float);a=np.asarray(a,complex)
    cadence=float(np.min(np.diff(t)));target,lags=E.lag_grid(model,t[-1]-t[0],cadence)
    # Higher modes can turn over before the original heading grid's first
    #250-unit lag. Extend that SAME sqrt2 grid down to saved cadence, retaining
    # every original heading lag and the exact-endpoint/nearest-half rule.
    low=25. if model=='spring_ring' else 250.
    lower=low/2.**(np.arange(1,40)/2.);lower=lower[lower>=cadence/np.sqrt(2)]
    extra=np.maximum(1,np.ceil(lower/cadence-.5-1e-10))*cadence
    mapping={float(l):float(tar) for tar,l in zip(target,lags)}
    for tar,l in zip(lower,extra):mapping.setdefault(float(l),float(tar))
    lags=np.array(sorted(mapping));target=np.array([mapping[l] for l in lags])
    nm=a.shape[1];nl=len(lags);n=np.zeros((nl,nm));d=n.copy();counts=np.zeros(nl)
    phase_n=np.zeros(nl);phase_d=np.zeros(nl)
    if Pabs is None:Pabs=np.abs(a[:,1]/a[:,0])
    for j,lag in enumerate(lags):
        i0=np.arange(len(t));i1=np.searchsorted(t,t+lag)
        ok=i1<len(t);i0=i0[ok];i1=i1[ok]
        ok=np.isclose(t[i1]-t[i0],lag,atol=cadence*1e-5);i0=i0[ok];i1=i1[ok]
        if not len(i0):continue
        x,y=a[i0],a[i1]
        n[j]=np.real(y*np.conjugate(x)).sum(0)
        d[j]=.5*(abs(x)**2+abs(y)**2).sum(0);counts[j]=len(i0)
        valid=(Pabs[i0]>=.1)&(Pabs[i1]>=.1)
        phase_n[j]=np.cos(np.angle(y[valid,1]*np.conj(x[valid,1]))).sum();phase_d[j]=valid.sum()
    return dict(target=target,lags=lags,num=n,den=d,counts=counts,phase_num=phase_n,phase_den=phase_d,
                spectrum=np.mean(abs(a)**2,axis=0),frames=len(t))

def aggregate_modes(seeds,c,identity):
    # Same nonuniform exact-endpoint lag grid for all exemplar seeds.
    lags=seeds[0]['lags'];nm=seeds[0]['num'].shape[1];ns=len(seeds)
    nums=np.array([s['num'] for s in seeds]);dens=np.array([s['den'] for s in seeds])
    counts=np.array([s['counts'] for s in seeds]);pn=np.array([s['phase_num'] for s in seeds]);pd=np.array([s['phase_den'] for s in seeds])
    with np.errstate(divide='ignore',invalid='ignore'):corr=nums.sum(0)/dens.sum(0);phase=pn.sum(0)/pd.sum(0)
    w=np.random.default_rng(E.stable_seed(identity)).multinomial(ns,np.full(ns,1/ns),size=2000)
    spec=np.mean([s['spectrum'] for s in seeds],axis=0);S=spec/spec[0]
    C=c.get('C_mean',c['size']);base=dict(condition_id=c['condition_id'],model=c['model'],C_label=c['size'],f=c['fraction'],C_mean=C,n_seeds=ns)
    rates=[];curves=[]
    for n in range(nm):
        cv=corr[:,n];valid=np.flatnonzero((cv>.03)&(cv<.995)&(counts.sum(0)>0))
        gamma=se=lo=hi=lag=np.nan;status='unresolved_no_admissible_lag'
        if len(valid):
            j=valid[np.argmin(abs(cv[valid]-.5))];lag=lags[j];gamma=-np.log(cv[j])/lag
            with np.errstate(divide='ignore',invalid='ignore'):
                boot=-np.log((w@nums[:,j,n])/(w@dens[:,j,n]))/lag;lb=np.log(boot)
            good=np.isfinite(lb)
            if good.mean()>=.975:
                se=float(np.std(lb[good],ddof=1));lo,hi=np.quantile(boot[good],[.025,.975]);status='resolved'
            else:status='bootstrap_invalid_fraction_above_0.025'
        crossed=bool(np.any(cv>=.5) and np.any(cv<=.5))
        finite=np.isfinite(cv);bound=np.nan;direction='none'
        if finite.any() and not crossed:
            if np.all(cv[finite]>.5):bound=np.log(2)/np.max(lags[finite]);direction='upper'
            elif np.all(cv[finite]<.5):bound=np.log(2)/np.min(lags[finite]);direction='lower'
        rates.append(dict(base,n=n,q=2*np.pi*n/C,S=S[n],gamma=gamma,SE_log_gamma=se,gamma_ci_lower=lo,gamma_ci_upper=hi,
             analysis_lag=lag,status=status,half_crossing_observed=bool(np.any(cv>=.5) and np.any(cv<=.5)),
             crossing_rate_bound=bound,crossing_bound_direction=direction,n_star=int(np.argmax(S[1:])+1),lambda_star=C/(int(np.argmax(S[1:])+1))))
        for j,l in enumerate(lags):curves.append(dict(base,n=n,q=2*np.pi*n/C,lag=l,correlation=cv[j],pairs=int(counts[:,j].sum()),phase_only_a1=phase[j] if n==1 else np.nan,phase_eligible_pairs=int(pd[:,j].sum()) if n==1 else 0))
    return rates,curves

def model_test(rates,alpha):
    rows=[];collapse=[]
    # Fixed from round5, before inspecting targeted mode evidence.
    transition=dict(wave_pinning=.03,holmes_model3=.02,spring_ring=.005,goryachev=.01)
    for r in rates:
        if num(r,'n')!=1 or not np.isfinite(num(r,'gamma')) or num(r,'gamma')<=0 or r.get('status')!='resolved' or str(r.get('half_crossing_observed'))!='True':continue
        model=r['model'];f=num(r,'f');C=num(r,'C_label');q=num(r,'q')
        candidates=[s for s in rates if s['model']==model and num(s,'f')==f and num(s,'C_label')!=C and num(s,'n')>=1 and num(s,'gamma')>0 and s.get('status')=='resolved' and str(s.get('half_crossing_observed'))=='True']
        # Nearest three physical wavenumbers from OTHER sizes prevent reading
        # the same size trend out twice. One-sided predictions stay flagged.
        candidates=sorted(candidates,key=lambda s:abs(np.log(num(s,'q')/q)))[:8]
        if len(candidates)<3:continue
        x=np.log([num(s,'q') for s in candidates]);y=np.log([num(s,'gamma') for s in candidates])
        if np.ptp(x)<1e-5:continue
        slope,offset=np.polyfit(x,y,1);pred=np.exp(offset+slope*np.log(q));err=np.log(num(r,'gamma')/pred)
        aa=next((a for a in alpha if a['model']==model and num(a,'f')==f and num(a,'C_label')==C),{})
        rows.append(dict(model=model,f=f,C_label=C,q=q,gamma=num(r,'gamma'),gamma_other_sizes_prediction=pred,
            log_prediction_error=err,local_mode_slope=slope,alpha=num(aa,'alpha'),
            extrapolated=not (x.min()<=np.log(q)<=x.max()),other_size_mode_points=len(candidates),
            input_stratum='transition_or_stronger' if f>=transition[model] else 'below_transition',transition_input=transition[model]))
    for model in MODELS:
        rr=[r for r in rows if r['model']==model and r['f']>=transition[model]]
        kval=[r for r in rr if np.isfinite(r['alpha'])]
        collapse.append(dict(model=model,transition_input=transition[model],stratum='round5 transition input and all larger designed inputs; fixed before target results',n_comparisons=len(rr),log_gamma_RMSE=float(np.sqrt(np.mean([r['log_prediction_error']**2 for r in rr]))) if rr else np.nan,
             n_alpha_comparisons=len(kval),alpha_slope_RMSE=float(np.sqrt(np.mean([(r['alpha']-r['local_mode_slope'])**2 for r in kval]))) if kval else np.nan))
    return rows,collapse
