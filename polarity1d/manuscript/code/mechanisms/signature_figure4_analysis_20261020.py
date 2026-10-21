"""Finite-window heading and site estimands for Figure 4 E0.

Every time series is one whole realization. Heading cosines permit valid
endpoints across gaps; unwrapped angular and site MSDs never bridge a gap.
No diffusion or exponential law is assumed by the empirical persistence.
"""
from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.optimize import minimize_scalar

THRESHOLD=math.exp(-1.)


def sections(valid):
    valid=np.asarray(valid,dtype=bool)
    changes=np.diff(np.r_[False,valid,False].astype(np.int8))
    return list(zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1)))


def _correlation(values):
    values=np.asarray(values)
    n=len(values);length=1<<(2*n-1).bit_length()
    transformed=np.fft.fft(values,length)
    return np.fft.ifft(transformed*np.conj(transformed)).real[:n]


def heading_cosine(P,valid,max_lag):
    P=np.asarray(P,dtype=complex);valid=np.asarray(valid,dtype=bool)
    if P.ndim!=1 or valid.shape!=P.shape or max_lag>=len(P) or max_lag<1:
        raise ValueError('One series, matching validity and a supported lag range are required')
    unit=np.zeros(len(P),dtype=complex)
    if np.any(abs(P[valid])<=0):raise ValueError('A valid heading has zero magnitude')
    unit[valid]=P[valid]/abs(P[valid])
    numerator=_correlation(unit)[1:max_lag+1]
    count=np.rint(_correlation(valid.astype(float))[1:max_lag+1]).astype(np.int64)
    possible=np.arange(len(P)-1,len(P)-max_lag-1,-1,dtype=np.int64)
    numerator[count==0]=0.
    return numerator,count,possible


def section_msd(values,valid,max_lag):
    """Real-coordinate MSD numerators; both points stay in one valid section."""
    values=np.asarray(values,dtype=float);valid=np.asarray(valid,dtype=bool)
    if values.ndim!=1 or valid.shape!=values.shape:raise ValueError('One coordinate and matching validity required')
    numerator=np.zeros(max_lag);count=np.zeros(max_lag,dtype=np.int64)
    for start,stop in sections(valid):
        x=values[start:stop]
        if len(x)<2:continue
        if not np.isfinite(x).all():raise ValueError('Valid coordinate is not finite')
        x=x-x[0];m=min(max_lag,len(x)-1);lags=np.arange(1,m+1)
        squared=np.r_[0.,np.cumsum(x*x)]
        terms=squared[len(x)-lags]+squared[-1]-squared[lags]-2*_correlation(x)[lags]
        tolerance=1e-10*max(1.,float(squared[-1]))
        if np.min(terms)<-tolerance:raise FloatingPointError('MSD cancellation exceeds the declared roundoff allowance')
        numerator[:m]+=np.maximum(terms,0.);count[:m]+=len(x)-lags
    possible=np.maximum(len(values)-np.arange(1,max_lag+1),0)
    return numerator,count,possible


def angular_coordinate(P,valid):
    values=np.full(len(P),np.nan)
    for start,stop in sections(valid):values[start:stop]=np.unwrap(np.angle(P[start:stop]))
    return values


def physical_coordinate(P,circumference,valid,*,direct=False):
    """Reset the origin per track; primary arc increments use midpoint C."""
    P=np.asarray(P,dtype=complex);C=np.asarray(circumference,dtype=float)
    if C.shape!=P.shape or np.any(C<=0):raise ValueError('Positive actual circumference at every scalar frame required')
    values=np.full(len(P),np.nan)
    for start,stop in sections(valid):
        theta=np.unwrap(np.angle(P[start:stop]));theta-=theta[0]
        if direct:values[start:stop]=theta*C[start:stop]/(2*np.pi)
        else:
            steps=np.diff(theta)*.5*(C[start+1:stop]+C[start:stop-1])/(2*np.pi)
            values[start:stop]=np.r_[0.,np.cumsum(steps)]
    return values


def ratio(numerator,denominator):
    numerator=np.asarray(numerator,dtype=float);denominator=np.asarray(denominator,dtype=float)
    return np.divide(numerator,denominator,out=np.full(np.broadcast_shapes(numerator.shape,denominator.shape),np.nan),where=denominator>0)


def persistence(curve,coverage,cadence):
    curve=np.asarray(curve,dtype=float);coverage=np.asarray(coverage,dtype=float)
    if curve.ndim!=1 or len(curve)!=len(coverage) or cadence<=0:raise ValueError('Finite cadence and one matching curve required')
    hits=np.flatnonzero(np.isfinite(curve)&(curve<=THRESHOLD))
    if len(hits):
        i=int(hits[0]);previous=1. if i==0 else curve[i-1]
        lower=float(i*cadence);upper=float((i+1)*cadence)
        if not np.isfinite(previous) or not np.isfinite(curve[:i+1]).all():
            return dict(status='unavailable_lag_gap',estimate=None,lower=0.,upper=math.inf,
                        coverage=float(coverage[i]),crossing_bracket=[lower,upper])
        estimate=lower+(previous-THRESHOLD)/(previous-curve[i])*cadence
        cov=float(coverage[i]) if i==0 else float(min(coverage[i-1],coverage[i]))
        unresolved=estimate<2*cadence
        return dict(status='unresolved_below_two_cadences' if unresolved else 'resolved',
            estimate=None if unresolved else float(estimate),descriptive_crossing_estimate=float(estimate),
            lower=0. if unresolved else float(estimate),upper=2*cadence if unresolved else float(estimate),
            coverage=cov,crossing_bracket=[lower,upper])
    if len(curve) and np.isfinite(curve).all():
        return dict(status='right_censored',estimate=None,lower=len(curve)*cadence,upper=math.inf,
                    coverage=float(coverage[-1]),crossing_bracket=None)
    return dict(status='unavailable_no_contiguous_lag_curve',estimate=None,lower=0.,upper=math.inf,
                coverage=0.,crossing_bracket=None)


def persistence_draws(curves,cadence):
    """Bootstrap intervals keep right censoring and subcadence resolution."""
    curves=np.asarray(curves,dtype=float);finite=np.isfinite(curves)
    hits=finite&(curves<=THRESHOLD);observed=hits.any(axis=1);index=hits.argmax(axis=1)
    row=np.arange(len(curves));previous=curves[row,np.maximum(index-1,0)].copy();previous[index==0]=1.
    current=curves[row,index]
    with np.errstate(divide='ignore',invalid='ignore'):
        estimate=index*cadence+(previous-THRESHOLD)/(previous-current)*cadence
    finite_prefix=np.cumsum(~finite,axis=1)[row,index]==0
    valid_crossing=observed&np.isfinite(previous)&finite_prefix
    resolved=valid_crossing&(estimate>=2*cadence)
    unresolved=valid_crossing&~resolved
    censored=~observed&finite.all(axis=1)
    lower=np.zeros(len(curves));upper=np.full(len(curves),np.inf)
    lower[resolved]=estimate[resolved];upper[resolved]=estimate[resolved]
    upper[unresolved]=2*cadence;lower[censored]=curves.shape[1]*cadence
    return lower,upper,np.where(resolved,0,np.where(unresolved,1,np.where(censored,2,3)))


def interval_quantiles(lower,upper):
    lower=np.asarray(lower,dtype=float);upper=np.asarray(upper,dtype=float)
    return float(np.quantile(lower,.025,method='inverted_cdf')),float(np.quantile(upper,.975,method='inverted_cdf'))


def interval_ratio(numerator_lower,numerator_upper,denominator_lower,denominator_upper):
    a,b,c,d=np.broadcast_arrays(numerator_lower,numerator_upper,denominator_lower,denominator_upper)
    with np.errstate(divide='ignore',invalid='ignore'):
        low=np.where(np.isinf(d),0.,a/d)
        high=np.where((c<=0)|np.isinf(b),np.inf,b/c)
    return low,high


def signed_diffusivity(msd,cadence):
    """Free-intercept initial slope over exactly 2--10 saved cadence steps."""
    msd=np.asarray(msd,dtype=float);x=np.arange(2,11,dtype=float)*cadence
    y=msd[...,1:10]
    if y.shape[-1]!=9:raise ValueError('At least ten cadence lags required')
    center=x-x.mean();slope=np.sum(y*center,axis=-1)/float(center@center)
    valid=np.isfinite(y).all(axis=-1)
    return np.where(valid,slope/2,np.nan)


def exponential_comparison(curve,counts,cadence):
    curve=np.asarray(curve,dtype=float);counts=np.asarray(counts,dtype=float)
    time=np.arange(1,len(curve)+1)*cadence;valid=np.isfinite(curve)&(counts>0)
    if valid.sum()<2:return dict(status='unavailable',tau=None)
    x,y,w=time[valid],curve[valid],counts[valid];w=w/w.sum()
    low,high=math.log(cadence/1000),math.log(time[-1]*1e6)
    def loss(value):return float(np.sum(w*(np.exp(-x/math.exp(value))-y)**2))
    fitted=minimize_scalar(loss,bounds=(low,high),method='bounded')
    tau=math.exp(float(fitted.x));at_bound=abs(fitted.x-high)<.02 or abs(fitted.x-low)<.02
    return dict(status='fit_at_comparison_bound' if at_bound else 'comparison_fit',tau=tau,
                weighted_squared_error=loss(float(fitted.x)),used_lags=int(valid.sum()),empirical_bound_replaced=False)


def regime(contrast,raw_regions,basal_contrast):
    contrast=np.asarray(contrast,dtype=float);raw=np.asarray(raw_regions,dtype=int)
    low=contrast<.1*basal_contrast
    single=(~low)&(raw==1);multiple=(~low)&(raw>=2)
    other=~(low|single|multiple)
    fractions=np.array([single.mean(),multiple.mean(),low.mean(),other.mean()])
    names=('single_front','multiple_fronts','low_contrast','unclassified')
    majority=np.flatnonzero(fractions>.5)
    return fractions,names[int(majority[0])] if len(majority) else 'mixed'


def conservative_fhalf(fractions,tau_lower,tau_upper):
    """The first crossing of half the weakest-input value, with finite bounds."""
    f=np.asarray(fractions,dtype=float);lo=np.asarray(tau_lower,dtype=float);hi=np.asarray(tau_upper,dtype=float)
    if len(f)<2 or np.any(f<=0) or np.any(np.diff(f)<=0):raise ValueError('Ordered positive input strengths required')
    if not np.isfinite(hi[0]) or hi[0]!=lo[0]:
        below=np.flatnonzero(hi[1:]<.5*lo[0])+1
        return dict(status='bounded_weak_input_censored' if len(below) else 'inconclusive_weak_input_censored',
            estimate=None,lower=float(f[0]),upper=float(f[below[0]]) if len(below) else math.inf)
    target=.5*lo[0]
    for i in range(1,len(f)):
        if hi[i]<=target:
            if hi[i]!=lo[i] or hi[i-1]!=lo[i-1]:
                return dict(status='bounded_crossing',estimate=None,lower=float(f[i-1]),upper=float(f[i]))
            weight=(lo[i-1]-target)/(lo[i-1]-lo[i])
            value=float(math.exp(math.log(f[i-1])+weight*(math.log(f[i])-math.log(f[i-1]))))
            return dict(status='resolved',estimate=value,lower=value,upper=value)
        if lo[i]<=target<hi[i]:return dict(status='inconclusive_interval_overlap',estimate=None,lower=float(f[i-1]),upper=math.inf)
    return dict(status='right_censored_in_strength',estimate=None,lower=float(f[-1]),upper=math.inf)


def bootstrap_weights(n,draws,seed):
    if n<1 or draws<1:raise ValueError('At least one whole seed and one draw required')
    return np.random.default_rng(seed).multinomial(n,np.full(n,1./n),size=draws)


def finite_quantiles(values,axis=0):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        return np.nanquantile(values,[.025,.975],axis=axis)
