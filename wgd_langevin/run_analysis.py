"""Reproducible summary-data compatibility analysis; no raw-track likelihood."""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import simulate, metrics, msd

ROOT = Path(__file__).resolve().parent
GROUPS = ('reference', 'polyploid')
COLORS = ('#3569a8', '#c45642')


def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def read_data():
    with (ROOT / 'data/observations.csv').open() as f:
        data = list(csv.DictReader(f))
    for row in data:
        for key in ('value', 'low', 'high', 'cadence_h', 'duration_h', 'lag_h'):
            row[key] = float(row[key]) if row[key] else None
    return data


def surface(rows, tau, unit, speeds, cv=0, tolerance=1):
    v = speeds[None, :]
    score = np.zeros((len(tau), len(speeds)))
    compatible = np.ones_like(score, dtype=bool)
    for r in rows:
        name = r['metric']
        if name == 'speed': pred = unit['speed'][:, None] * v
        elif name == 'path_length': pred = unit['speed'][:, None] * v * r['duration_h']
        elif name == 'directness': pred = unit['directness'][:, None] + np.zeros_like(v)
        elif name == 'msd': pred = msd(r['lag_h'], v, tau[:, None], cv=cv)
        else: raise ValueError(name)
        width = max(np.log(r['high'] / r['low']) / 2, 0.001) * tolerance
        score += (np.log(pred / r['value']) / width)**2
        # Symmetric expansion in log-space of original reading intervals.
        lo = r['value'] * (r['low'] / r['value'])**tolerance
        hi = r['value'] * (r['high'] / r['value'])**tolerance
        compatible &= (pred >= lo) & (pred <= hi)
    return score, compatible


def best(score, tau, speeds):
    i, j = np.unravel_index(np.argmin(score), score.shape)
    return {'tau_h': float(tau[i]), 'Dr_per_h': float(1 / tau[i]),
            'propulsion_um_h': float(speeds[j]), 'score': float(score[i, j]),
            'tau_index': int(i), 'speed_index': int(j)}


def interval(mask, tau, speeds):
    i, j = np.where(mask)
    if not len(i): return {'compatible_points': 0}
    return {'compatible_points': int(len(i)), 'tau_min_h': float(tau[i].min()),
            'tau_max_h': float(tau[i].max()), 'speed_min_um_h': float(speeds[j].min()),
            'speed_max_um_h': float(speeds[j].max()),
            'touches_tau_lower_bound': bool(i.min() == 0),
            'touches_tau_upper_bound': bool(i.max() == len(tau)-1)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cells', type=int, default=2048)
    parser.add_argument('--grid', type=int, default=81)
    parser.add_argument('--substeps', type=int, default=32)
    parser.add_argument('--seed', type=int, default=20260921)
    args = parser.parse_args()
    out = ROOT / 'results'; out.mkdir(exist_ok=True)
    data = read_data()
    studies = list(dict.fromkeys(r['study'] for r in data))
    speeds = np.geomspace(0.1, 500, 1401)
    profiles, comparisons, checks, fitted, sensitivities = [], [], [], [], []
    fits = {}; banks = {}; report = {}
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for column, study in enumerate(studies):
        rows = [r for r in data if r['study'] == study]
        cadence, duration = rows[0]['cadence_h'], rows[0]['duration_h']
        knots = np.geomspace(cadence / 64, 48, args.grid)
        tau = np.geomspace(knots[0], knots[-1], 321)
        unit_raw = {key: [] for key in ('speed', 'directness', 'endpoint', 'directness_sd')}
        print(f'{study}: simulating {len(knots)} persistence values, {args.cells} cells', flush=True)
        for t in knots:
            p = simulate(1, t, cadence, duration, args.cells, max(args.substeps, int(np.ceil(8 * cadence / t))), args.seed)
            m = metrics(p, cadence)
            for key in unit_raw:
                unit_raw[key].append(float(m['directness'].std(ddof=1)) if key=='directness_sd' else float(m[key].mean()))
        unit = {k: np.interp(np.log(tau), np.log(knots), val) for k, val in unit_raw.items()}
        banks[study] = (tau, unit, cadence, duration)
        write_csv(out / f'{study}_simulation_bank.csv', [dict(tau_h=t, **{k: unit_raw[k][i] for k in unit_raw}) for i, t in enumerate(knots)])
        surfaces, masks = [], []
        fits[study] = {}; report[study] = {}
        for group, color in zip(GROUPS, COLORS):
            targets = [r for r in rows if r['role'] == 'fit' and r['group'] == group]
            score, mask = surface(targets, tau, unit, speeds)
            surfaces.append(score); masks.append(mask)
            fit = best(score, tau, speeds); fits[study][group] = fit
            bounds = interval(mask, tau, speeds)
            report[study][group] = {'best': fit, 'compatibility_range': bounds, 'persistence_informative_observations': not study.startswith('Jemaa')}
            fitted.append(dict(study=study, group=group, **fit, **bounds))
            for factor, cv in ((2, 0), (1, .35), (1, .7)):
                sc, ma = surface(targets, tau, unit, speeds, cv=cv, tolerance=factor)
                sensitivities.append(dict(study=study, group=group, tolerance_multiplier=factor,
                                          speed_cv=cv, **best(sc, tau, speeds), **interval(ma, tau, speeds)))
            bestj = score.argmin(axis=1)
            for i, t in enumerate(tau):
                profiles.append(dict(study=study, group=group, tau_h=t, Dr_per_h=1/t,
                                     best_speed_um_h=speeds[bestj[i]], score=score[i, bestj[i]],
                                     any_compatible=bool(mask[i].any())))
            axes[0, column].plot(tau, speeds[bestj], color=color, label=group)
            ii, jj = np.where(mask)
            if len(ii):
                lo = np.array([speeds[mask[i]].min() if mask[i].any() else np.nan for i in range(len(tau))])
                hi = np.array([speeds[mask[i]].max() if mask[i].any() else np.nan for i in range(len(tau))])
                axes[0, column].fill_between(tau, lo, hi, color=color, alpha=.3)
            axes[1, column].plot(tau, score.min(axis=1), color=color)
            # Model-based finite-track sampling sensitivity, not an empirical CI.
            if study == 'Mallin2023':
                ntracks = 50 if group == 'reference' else 49
                sc_sampling, ma_sampling = surface([r for r in targets if r['metric'] != 'directness'], tau, unit, speeds)
                direct = next(r for r in targets if r['metric']=='directness')
                sampling_se = unit['directness_sd'] / np.sqrt(ntracks)
                lo = direct['low'] - 1.96 * sampling_se
                hi = direct['high'] + 1.96 * sampling_se
                ma_sampling &= (unit['directness'][:,None] >= lo[:,None]) & (unit['directness'][:,None] <= hi[:,None])
                report[study][group]['model_based_sampling_sensitivity'] = interval(ma_sampling,tau,speeds)
            # Predict excluded measurements; no adjustment to accommodate checks.
            v, t, ti = fit['propulsion_um_h'], fit['tau_h'], fit['tau_index']
            for r in rows:
                if r['group'] != group or r['role'] not in ('check', 'diagnostic'): continue
                if r['metric'] == 'msd': pred = float(msd(r['lag_h'], v, t))
                elif r['metric'] == 'endpoint': pred = float(unit['endpoint'][ti] * v)
                elif r['metric'] == 'msd_slope':
                    # Window sensitivity is reported separately; this window is an assumption.
                    times = np.geomspace(cadence, duration / 2, 30)
                    pred = float(np.polyfit(np.log(times), np.log(msd(times, v, t)), 1)[0])
                else: continue
                checks.append(dict(id=r['id'], study=study, group=group, metric=r['metric'],
                                   predicted=pred, low=r['low'], high=r['high'],
                                   inside_reading_interval=bool(r['low'] <= pred <= r['high']),
                                   status='window_assumed_diagnostic' if r['role']=='diagnostic' else 'same_experiment_predictive_check'))
        a, b = surfaces
        # No AIC/p-values: these are distance scores on approximate, correlated summaries.
        null_i, null_j = np.unravel_index(np.argmin(a+b), a.shape)
        speed_i = int(np.argmin(a.min(axis=1)+b.min(axis=1)))
        persist_j = int(np.argmin(a.min(axis=0)+b.min(axis=0)))
        models = {'shared': ((null_i, null_j), (null_i, null_j)),
                  'speed_only': ((speed_i, int(a[speed_i].argmin())), (speed_i, int(b[speed_i].argmin()))),
                  'persistence_only': ((int(a[:,persist_j].argmin()), persist_j), (int(b[:,persist_j].argmin()), persist_j)),
                  'speed_and_persistence': (np.unravel_index(a.argmin(), a.shape), np.unravel_index(b.argmin(), b.shape))}
        feasible = {'shared': bool((masks[0] & masks[1]).any()),
                    'speed_only': bool((masks[0].any(axis=1) & masks[1].any(axis=1)).any()),
                    'persistence_only': bool((masks[0].any(axis=0) & masks[1].any(axis=0)).any()),
                    'speed_and_persistence': bool(masks[0].any() and masks[1].any())}
        for name, ((i,j),(k,l)) in models.items():
            comparisons.append(dict(study=study, model=name, distance_score=float(a[i,j]+b[k,l]),
                                    compatible_with_all_fit_intervals=feasible[name],
                                    reference_speed=speeds[j], polyploid_speed=speeds[l],
                                    reference_tau=tau[i], polyploid_tau=tau[k],
                                    speed_ratio=speeds[l]/speeds[j], tau_ratio=tau[k]/tau[i]))
        axes[0, column].set(title=study, xscale='log', yscale='log', ylabel='Propulsion (µm/h)', xlabel='Persistence time (h)')
        axes[1, column].set(xscale='log', yscale='symlog', ylim=(0, 10000), ylabel='Summary distance score', xlabel='Persistence time (h)')
        axes[0, column].legend(fontsize=8)
    fig.suptitle('Lines: conditional best fit; shading: compatible with all reading intervals\nSpeed-only studies have broad parameter ambiguity; scores are not likelihoods')
    fig.tight_layout(); fig.savefig(out/'parameter_profiles.png', dpi=160); plt.close(fig)

    # Crossover is a prediction, not used in fitting.
    xfits = fits['Xuan2018']
    def difference(t):
        return float(msd(t, xfits['polyploid']['propulsion_um_h'], xfits['polyploid']['tau_h']) -
                     msd(t, xfits['reference']['propulsion_um_h'], xfits['reference']['tau_h']))
    ts = np.geomspace(.001, 100, 500)
    cross = [brentq(difference, t0, t1) for t0,t1 in zip(ts[:-1], ts[1:]) if difference(t0)*difference(t1)<0]
    checks.append(dict(id='xuan_crossover', study='Xuan2018', group='both', metric='crossover',
                       predicted=cross[0] if cross else None, low=2, high=3,
                       inside_reading_interval=bool(cross and 2<=cross[0]<=3), status='same_experiment_predictive_check'))

    # Observation/prediction plots.
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    for group, color in zip(GROUPS, COLORS):
        fit = fits['Xuan2018'][group]
        ts = np.geomspace(1/6, 12, 150)
        ax[0,0].plot(ts, msd(ts, fit['propulsion_um_h'], fit['tau_h']), color=color, label=group)
        for r in data:
            if r['study']=='Xuan2018' and r['metric']=='msd' and r['group']==group:
                ax[0,0].errorbar(r['lag_h'], r['value'], yerr=[[r['value']-r['low']], [r['high']-r['value']]], fmt='o' if r['role']=='fit' else 's', color=color)
        ti = fit['tau_index']; unit = banks['Xuan2018'][1]
        ax[0,1].bar(group, unit['speed'][ti]*fit['propulsion_um_h'], color=color, alpha=.65)
        r = next(r for r in data if r['study']=='Xuan2018' and r['metric']=='speed' and r['group']==group)
        ax[0,1].errorbar(group, r['value'], yerr=[[r['value']-r['low']], [r['high']-r['value']]], fmt='ko')
        for row, metric in enumerate(('path_length','endpoint')):
            fit = fits['Mallin2023'][group]; ti=fit['tau_index']; unit=banks['Mallin2023'][1]
            pred = fit['propulsion_um_h'] * (unit['speed'][ti]*24 if metric=='path_length' else unit['endpoint'][ti])
            ax[row,2].bar(group, pred, color=color, alpha=.65)
            r=next(r for r in data if r['study']=='Mallin2023' and r['metric']==metric and r['group']==group)
            ax[row,2].errorbar(group,r['value'],yerr=[[r['value']-r['low']],[r['high']-r['value']]],fmt='ko')
    ax[0,0].axvspan(2,3,color='grey',alpha=.12,label='reported crossover range')
    ax[0,0].set(xscale='log',yscale='log',xlabel='Lag (h)',ylabel='MSD (µm²)',title='Xuan: MSD fit and early-lag check'); ax[0,0].legend(fontsize=8)
    ax[0,1].set(ylabel='Sampled speed (µm/h)',title='Xuan: speed fit')
    ax[0,2].set(ylabel='Accumulated distance (µm)',title='Mallin: fit')
    ax[1,2].set(ylabel='Endpoint distance (µm)',title='Mallin: excluded-metric check')
    # Show illustrative speed-only explanation at fixed 2h persistence, not arbitrary global best.
    for col, study in enumerate(('Jemaa2023','Jemaa2017')):
        tau, unit, cadence, duration=banks[study]; t=2.0
        att=np.interp(np.log(t),np.log(tau),unit['speed'])
        for group,color in zip(GROUPS,COLORS):
            r=next(r for r in data if r['study']==study and r['group']==group and r['role']=='fit')
            ax[1,col].bar(group,r['value'],color=color,alpha=.65)
            ax[1,col].errorbar(group,r['value'],yerr=[[r['value']-r['low']],[r['high']-r['value']]],fmt='ko')
        ax[1,col].set(title=f'{study}: sampled speed fit',ylabel='Sampled speed (µm/h)')
    fig.suptitle('Bars/curves: homogeneous active-particle model; points: approximate literature values\nError bars are reading/sensitivity intervals, not confidence intervals')
    fig.tight_layout(); fig.savefig(out/'observations_vs_model.png',dpi=160); plt.close(fig)

    # Tail scenarios, expressly not fitted: unpublished sample sizes prevent tail likelihood.
    fig, axes = plt.subplots(1, 2, figsize=(11,4))
    tail_rows=[]
    for col,study in enumerate(('Jemaa2023','Jemaa2017')):
        tau,unit,cadence,duration=banks[study]
        att=np.interp(np.log(2),np.log(tau),unit['speed'])
        for cv in (0,.35,.7):
            target=next(r['value'] for r in data if r['study']==study and r['group']=='polyploid' and r['role']=='fit')
            p=simulate(target/att,2,cadence,duration,n=2048,substeps=args.substeps,seed=args.seed+1,cv=cv)
            vals=metrics(p,cadence)['speed']; ordered=np.sort(vals)
            axes[col].plot(ordered,np.arange(1,len(vals)+1)/len(vals),label=f'CV={cv}')
            tail_rows.append(dict(study=study,group='polyploid',assumed_tau_h=2,assumed_cv=cv,mean_sampled_speed=float(vals.mean()),p95=float(np.quantile(vals,.95)),p99=float(np.quantile(vals,.99)),status='illustrative_not_fitted'))
        axes[col].set(title=study,xlabel='Per-track sampled speed (µm/h)',ylabel='Cumulative probability'); axes[col].legend()
    fig.suptitle('Heterogeneity scenarios at fixed 2h persistence — CV is not identified by available tail summaries')
    fig.tight_layout(); fig.savefig(out/'heterogeneity_sensitivity.png',dpi=160); plt.close(fig)
    # Slope-window sensitivity: no invented equivalence between exponent and tau.
    slopes=[]
    for group in GROUPS:
        fit=xfits[group]
        for lo,hi in ((1/6,1),(1/6,6),(1,6),(2,10),(1/6,12)):
            t=np.geomspace(lo,hi,30)
            slopes.append(dict(group=group,lag_min_h=lo,lag_max_h=hi,slope=float(np.polyfit(np.log(t),np.log(msd(t,fit['propulsion_um_h'],fit['tau_h'])),1)[0])))
    for name, rows in [('parameter_fits',fitted),('parameter_profiles',profiles),('model_comparison',comparisons),('predictive_checks',checks),('sensitivity',sensitivities),('heterogeneity_scenarios',tail_rows),('slope_window_sensitivity',slopes)]:
        # Different compatibility masks can yield different keys.
        keys=list(dict.fromkeys(k for r in rows for k in r))
        write_csv(out/f'{name}.csv',[{k:r.get(k,'') for k in keys} for r in rows])
    report['run']={'arguments':vars(args),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                   'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'model.py',ROOT/'run_analysis.py',ROOT/'data/observations.csv')},
                   'interpretation':'Reading-interval compatibility, not confidence intervals. Dt=0, localization_sd=0; CV=0 primary fit.'}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Finished. Results:',out,flush=True)

if __name__=='__main__': main()
