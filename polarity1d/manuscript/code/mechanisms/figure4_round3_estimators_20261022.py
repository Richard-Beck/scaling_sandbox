"""Shared point/whole-seed estimators for Figure 4 round 3.

The fixed five-size neighborhood is a design property, including when data are
missing. Lag selection uses observed curves once and stays fixed in bootstrap.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.fft import rfft,irfft,next_fast_len

CHEMICAL_SIZES = (20., 30., 40., 50., 60.)
SPRING_SIZES = (1., 1.5, 2., 2.5, 3.)
DRAW_COUNT = 2000


def write_csv(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields); writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(',', ':')) if isinstance(v, (list, dict, tuple)) else v for k, v in row.items()})


def stable_seed(identity):
    return 2026102204 + int(hashlib.sha256(str(identity).encode()).hexdigest()[:8], 16)


def lag_grid(model='chemical', horizon=math.inf, cadence=None):
    """sqrt(2) target grid; round to the closest native endpoint, ties down.

    Returned target and actual lags retain the cadence discretization explicitly.
    No interpolated time samples are used.
    """
    low = 25. if model == 'spring_ring' else 250.
    target = low * 2. ** (np.arange(13) / 2.)
    actual = target.copy() if cadence is None else np.maximum(1, np.ceil(target / cadence - .5 - 1e-10)) * cadence
    keep = actual <= horizon
    _, unique = np.unique(actual[keep], return_index=True)
    return target[keep][unique], actual[keep][unique]


def seed_moments(P, valid, cadence, lags):
    P = np.asarray(P, complex); valid = np.asarray(valid, bool)
    result = np.zeros((3, len(lags)))
    for j, lag in enumerate(lags):
        step = int(round(lag / cadence))
        if step < 1 or step >= len(P) or not np.isclose(step * cadence, lag):
            result[:, j] = np.nan; continue
        mask = valid[:-step] & valid[step:] & np.isfinite(P[:-step]) & np.isfinite(P[step:])
        phase = np.angle(P[:-step][mask] * np.conj(P[step:][mask]))
        result[:, j] = np.cos(phase).sum(), mask.sum(), len(P) - step
    return result


def percentile(values, conservative=False):
    values = np.asarray(values, float)
    finite = np.isfinite(values)
    if not finite.any(): return math.nan, math.nan
    if conservative and finite.mean() < .975: return -math.inf, math.inf
    return tuple(float(x) for x in np.quantile(values[finite], [.025, .975], method='inverted_cdf'))


def summarize(moments, circumferences, identity, target_lags, lags, draws=DRAW_COUNT):
    moments = np.asarray(moments, float); n = len(moments)
    weights = np.random.default_rng(stable_seed(identity)).multinomial(n, np.full(n, 1 / n), size=draws)
    numerator, count, possible = moments[:, 0], moments[:, 1], moments[:, 2]
    with np.errstate(invalid='ignore', divide='ignore'):
        C = numerator.sum(0) / count.sum(0)
        Cboot = (weights @ numerator) / (weights @ count)
        coverage = count.sum(0) / possible.sum(0)
        k = -np.log(C) / lags
        kboot = -np.log(Cboot) / lags
    return dict(target_lags=np.asarray(target_lags), lags=np.asarray(lags), C=C, Cboot=Cboot, coverage=coverage,
                k=k, kboot=kboot, circumference=float(np.mean(circumferences)),
                circumference_boot=weights @ np.asarray(circumferences) / n, n_seeds=n,
                eligible_pairs=count.sum(0), possible_pairs=possible.sum(0), seed_moments=moments,bootstrap_seed=stable_seed(identity))


def fixed_neighborhood(size, model='chemical'):
    sizes = SPRING_SIZES if model == 'spring_ring' else CHEMICAL_SIZES
    i = sizes.index(float(size))
    return sizes[max(0, i - 1):min(len(sizes), i + 2)]


def analysis_lag(summary):
    ok = (summary['C'] > .03) & (summary['C'] < .995)
    ix = np.flatnonzero(ok)
    if not len(ix): return None
    return int(ix[np.argmin(abs(summary['C'][ix]-.5))])


def local_alpha(groups, size, model='chemical'):
    neighborhood = fixed_neighborhood(size, model)
    base = dict(C_label=size, neighborhood=neighborhood, n_neighbors=len(neighborhood))
    absent = [c for c in neighborhood if c not in groups]
    if absent: return dict(base, status='unavailable_missing_fixed_neighbor', missing_sizes=absent)
    ss = [groups[c] for c in neighborhood]
    # Actual endpoint grids must match, not merely their target lag labels.
    shared = sorted(set.intersection(*(set(s['lags']) for s in ss)))
    candidates = []
    for lag in shared:
        indices = [int(np.flatnonzero(s['lags'] == lag)[0]) for s in ss]
        C = np.array([s['C'][i] for s, i in zip(ss, indices)])
        cov = np.array([s['coverage'][i] for s, i in zip(ss, indices)])
        if np.all((C > .03) & (C < .995) & (cov >= .5)):
            score = float(np.max(abs(np.log(-np.log(C)) - math.log(math.log(2)))))
            candidates.append((score, lag, indices))
    if not candidates: return dict(base, status='unavailable_no_shared_lag_with_C_0.03_0.995_and_coverage_0.5')
    score, lag, indices = min(candidates, key=lambda z: (z[0], z[1]))
    x = np.log([s['circumference'] for s in ss]); center = x - x.mean()
    y = np.log([s['k'][i] for s, i in zip(ss, indices)])
    alpha = -float(y @ center / (center @ center))
    xb = np.log(np.array([s['circumference_boot'] for s in ss])); xb -= xb.mean(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        yb = np.log(np.array([s['kboot'][:, i] for s, i in zip(ss, indices)]))
        boot = -np.sum(xb * yb, axis=0) / np.sum(xb ** 2, axis=0)
    lo, hi = percentile(boot, conservative=True)
    own = groups[float(size)]; j = int(np.flatnonzero(own['lags'] == lag)[0])
    klo, khi = percentile(own['kboot'][:, j], conservative=True)
    return dict(base, status='resolved', alpha=alpha, alpha_ci_lower=lo, alpha_ci_upper=hi,
                bootstrap_valid_fraction=float(np.isfinite(boot).mean()), lag=lag, lag_selection_score=score,
                C_mean=own['circumference'], C=own['C'][j], coverage=own['coverage'][j], k=own['k'][j],
                k_ci_lower=klo, k_ci_upper=khi, n_seeds=own['n_seeds'])


def crossing(fs, C):
    """Round-2 first downcrossing and missing-gap bounds, forced included."""
    fs = np.asarray(fs, float); C = np.asarray(C, float)
    last = None; gap = False
    for i, value in enumerate(C):
        if not np.isfinite(value): gap = True; continue
        if value <= .5:
            if i == 0: return dict(status='left_censored', estimate=math.nan, lower=0., upper=float(fs[0]))
            if last is None or gap:
                return dict(status='bounded_missing_gap', estimate=math.nan, lower=float(fs[last]) if last is not None else 0., upper=float(fs[i]))
            w = (C[last] - .5) / (C[last] - value)
            value = float(np.exp(np.log(fs[last]) + w * np.log(fs[i] / fs[last])))
            return dict(status='resolved', estimate=value, lower=value, upper=value)
        last = i
    return dict(status='right_censored' if not gap else 'unavailable_inputs', estimate=math.nan,
                lower=float(fs[-1]) if len(fs) and not gap else 0., upper=math.inf)


def tolerance(fs, C, Cboot):
    point = crossing(fs, C); bounds = [crossing(fs, c) for c in np.asarray(Cboot).T]
    lower = np.array([b['lower'] for b in bounds]); upper = np.array([b['upper'] for b in bounds])
    return dict(fhalf=point['estimate'], fhalf_lower=point['lower'], fhalf_upper=point['upper'], status=point['status'],
                fhalf_ci_lower=float(np.quantile(lower, .025, method='inverted_cdf')),
                fhalf_ci_upper=float(np.quantile(upper, .975, method='inverted_cdf')),
                bootstrap_valid_fraction=float(np.mean([b['status'] == 'resolved' for b in bounds])),
                bootstrap_lower=lower, bootstrap_upper=upper)


def tolerance_fold(small, large):
    with np.errstate(divide='ignore', invalid='ignore'):
        lo = np.where(np.isinf(small['bootstrap_upper']), 0., large['bootstrap_lower'] / small['bootstrap_upper'])
        hi = np.where(small['bootstrap_lower'] <= 0, np.inf, large['bootstrap_upper'] / small['bootstrap_lower'])
    point = large['fhalf'] / small['fhalf'] if np.isfinite(large['fhalf']) and small['fhalf'] > 0 else math.nan
    return dict(fhalf_fold=point, fhalf_fold_ci_lower=float(np.quantile(lo, .025, method='inverted_cdf')),
                fhalf_fold_ci_upper=float(np.quantile(hi, .975, method='inverted_cdf')),
                status='resolved' if np.isfinite(point) else 'bounded_or_unavailable',
                bootstrap_valid_fraction=float(np.mean(np.isfinite(lo) & np.isfinite(hi))))


def smooth_flicker(ids, hold=5):
    """A new observed dominant must persist hold frames; no missing-ID bridge."""
    ids = np.asarray(ids, int); result = np.full(len(ids), -1, int)
    accepted = -1; start = 0
    while start < len(ids):
        stop = start + 1
        while stop < len(ids) and ids[stop] == ids[start]: stop += 1
        identity = int(ids[start])
        if identity < 0: accepted = -1
        elif accepted < 0 or identity == accepted or stop - start >= hold: accepted = identity
        result[start:stop] = accepted
        start = stop
    return result


def tenure_sections(ids):
    ids = np.asarray(ids)
    starts = np.r_[0, np.flatnonzero(ids[1:] != ids[:-1]) + 1]
    stops = np.r_[starts[1:], len(ids)]
    return [(int(a), int(b)) for a, b in zip(starts, stops) if ids[a] >= 0]


def drift_jump(time, circumference, dominant_ids, track_headings, *, hold=5, range_fraction=.25, P=None,max_lag=None):
    """track_headings[frame] maps retained identity to its native region heading.

    Dominance tenures and switch exposure use the persistence-filtered IDs.
    Flicker suppression retains the accepted identity's centroid, rather than
    substituting the instantaneous winner's coordinate. A missing accepted
    centroid splits only its displacement segment; it does not erase an
    observed dominant-ID interval or change its tenure length.
    """
    time = np.asarray(time, float); C = np.broadcast_to(circumference, time.shape)
    dt = float(np.median(np.diff(time))); dominant = smooth_flicker(dominant_ids, hold)
    tenures=tenure_sections(dominant);durations=[(b-a)*dt for a,b in tenures]
    heading = np.array([track_headings[i].get(int(identity), math.nan) for i, identity in enumerate(dominant)])
    ids=dominant.copy()
    ids[~np.isfinite(heading)] = -1
    sections = tenure_sections(ids)
    valid_time=float(np.sum((dominant[:-1]>=0)&(dominant[1:]>=0))*dt)
    base = dict(flicker_frames=hold, MSD_range_fraction=range_fraction, n_tenures=len(tenures),n_centroid_segments=len(sections),
                source_cadence=dt,flicker_hold_time=hold*dt,jump_heading_mean_window=5*dt,
                valid_time=valid_time, n_switches=0,n_switches_with_missing_adjacent_centroid=0,
                Ds=math.nan, MSD_exponent=math.nan, lambda_s=math.nan, mean_cos_jump=math.nan,
                lambda_eff=math.nan, status='unavailable_no_valid_tenures')
    if not durations: return base, []
    maxstep = int(np.floor((np.median(durations)*range_fraction if max_lag is None else max_lag)/dt))
    numerator = np.zeros(maxstep); count = np.zeros(maxstep)
    for a, b in sections:
        theta = np.unwrap(heading[a:b]); x = np.r_[0., np.cumsum(np.diff(theta) * .5 * (C[a+1:b] + C[a:b-1]) / (2*np.pi))]
        n=b-a;m=min(maxstep,n-1)
        if m<1:continue
        if m<16:
            for step in range(1,m+1):
                d=x[step:]-x[:-step];numerator[step-1]+=d@d;count[step-1]+=len(d)
        else:
            # Exact non-circular MSD sums from linear FFT autocorrelation and
            # prefix square sums. Translation centering controls cancellation.
            x=x-x.mean();prefix=np.r_[0.,np.cumsum(x*x)];fft=rfft(x,next_fast_len(2*n-1));cross=irfft(fft*np.conj(fft),next_fast_len(2*n-1))[:m+1]
            steps=np.arange(1,m+1);squares=prefix[n]-prefix[steps]+prefix[n-steps]
            numerator[:m]+=np.maximum(squares-2*cross[1:m+1],0.);count[:m]+=n-steps
    lag = np.arange(1, maxstep+1) * dt
    with np.errstate(invalid='ignore', divide='ignore'): MSD = numerator / count
    ok = np.isfinite(MSD) & (count > 0)
    Ds = float(np.dot(lag[ok], MSD[ok]) / (2 * np.dot(lag[ok], lag[ok]))) if ok.sum() >= 2 else math.nan
    power = ok & (MSD > 0)
    exponent = float(np.polyfit(np.log(lag[power]), np.log(MSD[power]), 1)[0]) if power.sum() >= 3 else math.nan
    jump_heading=heading if P is None else np.angle(np.asarray(P,complex))
    jump_valid=np.isfinite(jump_heading) if P is None else np.isfinite(P)&(abs(np.asarray(P))>=.1)
    jumps = [];counted_switches=0;missing_adjacent_centroid=0
    for (a, b), (c, d) in zip(tenures[:-1], tenures[1:]):
        if b != c or dominant[a] == dominant[c]:continue
        counted_switches+=1
        if not np.isfinite(heading[b-1]) or not np.isfinite(heading[c]):missing_adjacent_centroid+=1
        if b<5 or c+5>len(heading):continue
        if P is None and (b-a<5 or d-c<5):continue
        if not np.all(jump_valid[b-5:b]) or not np.all(jump_valid[c:c+5]):continue
        before = np.mean(np.exp(1j * jump_heading[b-5:b])); after = np.mean(np.exp(1j * jump_heading[c:c+5]))
        if abs(before) > 1e-12 and abs(after) > 1e-12: jumps.append(float(np.angle(after * np.conj(before))))
    valid_time = base['valid_time']
    lam = counted_switches/valid_time if valid_time > 0 else math.nan
    mean_cos = float(np.cos(jumps).mean()) if jumps else math.nan
    base.update(Ds=Ds, MSD_exponent=exponent, lambda_s=lam, mean_cos_jump=mean_cos,
                lambda_eff=lam*(1-mean_cos) if jumps and len(jumps)==counted_switches else 0. if counted_switches==0 and valid_time>0 else math.nan,
                n_switches=counted_switches,n_switches_with_missing_adjacent_centroid=missing_adjacent_centroid,n_scored_switches=len(jumps),jump_heading_source='whole_P' if P is not None else 'tracked_region',
                median_tenure=float(np.median(durations)), MSD_max_lag=maxstep*dt,
                status='resolved' if np.isfinite(Ds) and len(jumps)==counted_switches else 'unavailable_MSD_range' if not np.isfinite(Ds) else 'unavailable_jump_heading_windows',
                jump_angles=jumps)
    return base, [dict(lag=float(t), MSD=float(m), eligible_pairs=int(n)) for t, m, n in zip(lag, MSD, count)]


def wrapped_mixture_mle(P, circumference, cadence, step=1,seed_labels=None):
    P=np.asarray(P,complex); valid=(abs(P)>=.1)&np.isfinite(P)
    mask=valid[:-step]&valid[step:]
    if seed_labels is not None:
        labels=np.asarray(seed_labels);mask &= (labels[:-step]==labels[step:])&(labels[:-step]>=0)
    delta=np.angle(P[step:][mask]*np.conj(P[:-step][mask]))
    out=dict(detector_free_step=step,delta=cadence*step,n_increments=len(delta),status='unavailable_too_few_valid_increments')
    if len(delta)<20:return out
    # Fourier wrapped-normal expansion is stable for broad distributions;
    # narrow normals use image sums to avoid cancellation.
    dt=cadence*step
    def loss(p):
        D,lam=np.exp(p);variance=2*D*dt;weight=np.exp(-lam*dt)
        if variance<1:
            offsets=delta[:,None]+2*np.pi*np.arange(-1,2)
            components=np.exp(-offsets**2/(2*variance)) / np.sqrt(2*np.pi*variance)
            density=components.sum(1)
            derivative_D=(components*(-.5+offsets**2/(2*variance))).sum(1)
        else:
            modes=np.arange(1,int(np.ceil(np.sqrt(-2*np.log(1e-15)/variance)))+1)
            components=np.cos(delta[:,None]*modes)*np.exp(-.5*variance*modes**2)
            density=(1+2*components.sum(1))/(2*np.pi)
            derivative_D=2*(components*(-.5*variance*modes**2)).sum(1)/(2*np.pi)
        mixed=np.maximum(weight*density+(1-weight)/(2*np.pi),1e-300)
        gradient=np.array([-np.sum(weight*derivative_D/mixed),np.sum(lam*dt*weight*(density-1/(2*np.pi))/mixed)])
        return -float(np.log(mixed).sum()),gradient
    guesses=[]
    for scale in (.05,.3,.8):
        start=np.log([max(np.var(delta)/(2*dt)*scale,1e-14),max(.01/dt,1e-14)])
        fit=minimize(loss,start,method='L-BFGS-B',jac=True,bounds=[(math.log(1e-16),math.log(100/dt)),(math.log(1e-16),math.log(30/dt))])
        guesses.append(fit)
    fit=min(guesses,key=lambda r:r.fun);D,lam=np.exp(fit.x)
    lower=np.array([math.log(1e-16)]*2);upper=np.log([100/dt,30/dt]);at_boundary=bool(np.any(np.isclose(fit.x,lower,atol=.05)|np.isclose(fit.x,upper,atol=.05)))
    out.update(Dtheta=float(D),Ds=float(D*circumference**2/(4*np.pi**2)),lambda_eff=float(lam),
               negative_log_likelihood=float(fit.fun),status='boundary_solution' if fit.success and at_boundary else 'resolved' if fit.success else 'optimizer_not_converged',fit_at_boundary=at_boundary)
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['precision']);parser.add_argument('--input',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();source=args.input/'results/compact_index.csv' if args.input.is_dir() else args.input
    if source.suffix=='.csv':
        with source.open(newline='') as stream:records=list(csv.DictReader(stream))
    else:
        records=json.loads(source.read_text());records=records.get('seeds',records) if isinstance(records,dict) else records
    groups={}
    for row in records:
        if row.get('status','complete').lower() not in ['complete','completed','success','passed']:continue
        groups.setdefault((row['model'],float(row.get('C_label',row.get('size'))),float(row.get('f',row.get('fraction')))),[]).append(row)
    output=[]
    for key,rows in groups.items():
        values=[];circ=[]
        for row in rows:
            with np.load(row.get('compact_path',row.get('compact',row.get('path'))),allow_pickle=False) as z:
                t=z['time'].astype(float);P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float)
                target,lags=lag_grid(key[0],t[-1]-t[0],float(np.median(np.diff(t))))
                values.append(seed_moments(P,abs(P)>=.1,float(np.median(np.diff(t))),lags));circ.append(float(z['circumference'].mean()))
        s=summarize(values,circ,('compact',key),target,lags);j=analysis_lag(s)
        with np.errstate(invalid='ignore',divide='ignore'):
            logboot=np.log(s['kboot'][:,j]) if j is not None else np.full(DRAW_COUNT,np.nan)
        valid=float(np.isfinite(logboot).mean());se=float(np.std(logboot[np.isfinite(logboot)],ddof=1)) if valid>=.975 else None
        output.append(dict(condition_id=rows[0].get('condition_id',f'{key[0]}:C{key[1]:g}:f{key[2]:g}'),model=key[0],C_label=key[1],f=key[2],n_seeds=len(rows),se_log_k=se,analysis_lag=float(lags[j]) if j is not None else None,
                           reason='resolved' if se is not None else 'unavailable_no_analysis_lag' if j is None else 'bootstrap_valid_fraction_below_0.975',bootstrap_valid_fraction=valid,
                           target_precision_met=se is not None and se<=.07,additional_seeds=0 if se is None or se<=.07 or len(rows)>=384 else min(96,384-len(rows))))
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(dict(conditions=output),indent=2)+'\n')


if __name__=='__main__':main()
