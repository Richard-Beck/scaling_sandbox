"""Deterministic reduction against the retained De Belly Figure 2 observables.

No De Belly simulations are rerun. Calibration uses only A/C trajectories;
countercue trajectories and duration labels are never used by the optimizer.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numba
from numba import njit
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / 'polarity1d/results/data/figure2'
CAMPAIGNS = ROOT / 'polarity1d/campaigns'
PREFIX = 'debelly_2026_standalone_figure2_julia_activation_defaults_'
NAMES = ['tau_p', 'tau_v', 'alpha', 'beta', 'ps', 'vs', 'cue_gain']
# All seven optimization coordinates are logarithmic.
LOW = np.array([1., 1., 1e-6, 1e-6, .003, .003, 1e-5])
HIGH = np.array([20000., 20000., 2., 2., 2., 2., 2.])
LOG = np.array([True, True, True, True, True, True, True])


def encode(values):
    out = np.asarray(values, dtype=float).copy()
    out[LOG] = np.log(out[LOG])
    return out


def decode(values):
    out = np.asarray(values, dtype=float).copy()
    out[LOG] = np.exp(out[LOG])
    return out


def parameter_dict(theta):
    result = dict(zip(NAMES, map(float, theta)))
    result['sigma_p'] = result['sigma_v'] = 0.
    return result


@njit(cache=True)
def forcing(t, pulses):
    result = 0.
    for start, stop, taper, amplitude in pulses:
        if start <= t < stop:
            weight = 1.
            if taper > 0 and t > stop - taper:
                weight = .5 * (1. + np.cos(np.pi * (t - stop + taper) / taper))
            result += amplitude * weight
    return result


@njit(cache=True)
def rhs(t, p, v, theta, pulses):
    tp, tv, alpha, beta, ps, vs, gain = theta
    return (-p/tp + beta*v/(1.+abs(v)/vs) + gain*forcing(t, pulses),
            -v/tv + alpha*p/(1.+abs(p)/ps))


@njit(cache=True)
def integrate(theta, times, pulses, initial, max_step=.5):
    """RK4 aligned with samples and all cue events; one-sided event evaluation."""
    result = np.empty((len(times), 2))
    p, v = initial
    result[0, 0], result[0, 1] = p, v
    t = times[0]
    cap = min(max_step, .2*theta[0], .2*theta[1])
    for i in range(1, len(times)):
        target = times[i]
        while t < target - 1e-9:
            h = min(cap, target-t)
            for start, stop, taper, amplitude in pulses:
                for event in (start, stop-taper, stop):
                    if event > t + 1e-9:
                        h = min(h, event-t)
            a, b = rhs(t+1e-10, p, v, theta, pulses)
            c, d = rhs(t+h/2, p+h*a/2, v+h*b/2, theta, pulses)
            e, f = rhs(t+h/2, p+h*c/2, v+h*d/2, theta, pulses)
            g, j = rhs(t+h-1e-10, p+h*e, v+h*f, theta, pulses)
            p += h*(a+2*c+2*e+g)/6
            v += h*(b+2*d+2*f+j)/6
            t += h
        result[i, 0], result[i, 1] = p, v
    return result


def pulses_from_task(task):
    protocol = task['protocol']
    configs = protocol['pulses'] if protocol['name'] == 'finite_pulse_sequence' else [protocol]
    pulses = []
    for c in configs:
        asym = c.get('asymmetry', 0.)
        signed = sum(-(1+asym/2) if side == 'left' else (1-asym/2) for side in c['sides'])
        # Exact right-half minus left-half mean of the cosine spatial stimulus.
        amplitude = c['amplitude'] * c['width_fraction'] * signed
        pulses.append([c['start'], c['stop'], c.get('taper', 0.), amplitude])
    return np.asarray(pulses, dtype=float)


def load_data():
    manifests = {}
    paths = []
    for panel, suffix in [('A', 'panel_a_10000'), ('C', 'base'), ('E', 'linear50_response')]:
        path = CAMPAIGNS / (PREFIX + suffix + '.json')
        paths.append(path)
        manifests[panel] = {t['task_id']: t for t in json.loads(path.read_text())['tasks']}
    data = pd.read_csv(DATA / 'figure2_v2_polarity_trajectories.csv.gz')
    cases = []
    for (panel, task_id), frame in data.groupby(['panel', 'task_id'], sort=True):
        wide = frame.pivot(index='time_s', columns='field', values='x50_over_L').sort_index()
        task = manifests[panel][task_id]
        values = np.column_stack([wide.rac.to_numpy()-.5, .5-wide.mca_bound.to_numpy()])
        cases.append(dict(panel=panel, task_id=task_id, variant=frame.variant.iloc[0],
                          length=float(frame.length_um.iloc[0]), times=wide.index.to_numpy(),
                          values=values, task=task, pulses=pulses_from_task(task)))
    return cases, manifests, paths


def training_cases(cases, length, variant):
    a_variant = {'beta_0x': 'beta_zero', 'beta_1x': 'full'}.get(variant)
    return [c for c in cases if c['length'] == length and
            ((c['panel'] == 'C' and c['variant'] == variant) or
             (c['panel'] == 'A' and c['variant'] == a_variant))]


def objective_cases(cases):
    prepared = []
    for case in cases:
        t = case['times']
        last_cue = case['pulses'][:, 1].max()
        # Resolve early kinetics while avoiding thousands of redundant tail points.
        early = np.flatnonzero(t <= last_cue + 500.)
        late = np.flatnonzero(t > last_cue + 500.)
        early = early[::max(1, round(5./np.median(np.diff(t))))]
        late = late[::max(1, round(50./np.median(np.diff(t))))]
        idx = np.unique(np.r_[early, late, len(t)-1])
        weights = np.empty(len(idx))
        first = t[idx] <= last_cue + 500.
        for mask in [first, ~first]:
            weights[mask] = 1. / np.sqrt(max(1, mask.sum()) * 2 * len(cases))
        scale = np.maximum(.02, np.max(np.abs(case['values']), axis=0))
        prepared.append((t[idx], case['values'][idx], case['pulses'], weights, scale))
    return prepared


def fit(cases, starts, max_nfev, out, label):
    prepared = objective_cases(cases)
    def residual(z):
        theta = decode(z)
        pieces = []
        for t, y, pulses, weight, scale in prepared:
            prediction = integrate(theta, t, pulses, np.zeros(2))
            pieces.append(((prediction-y)/scale * weight[:, None]).ravel())
        return np.concatenate(pieces)
    candidates = []
    for index, start in enumerate(starts):
        tick = time.monotonic()
        result = least_squares(residual, np.clip(encode(start), encode(LOW)+1e-8, encode(HIGH)-1e-8),
                               bounds=(encode(LOW), encode(HIGH)), max_nfev=max_nfev,
                               ftol=1e-7, xtol=1e-7, gtol=1e-7, diff_step=1e-4)
        candidates.append(dict(start=index, cost=float(result.cost), nfev=result.nfev,
                               success=bool(result.success), message=result.message,
                               parameters=parameter_dict(decode(result.x)), theta=decode(result.x).tolist(),
                               wall_seconds=time.monotonic()-tick,
                               active_bounds={NAMES[i]: int(a) for i, a in enumerate(result.active_mask) if a},
                               jacobian_singular_values=np.linalg.svd(result.jac, compute_uv=False).tolist()))
        (out / f'fit_{label}.json').write_text(json.dumps(candidates, indent=2)+'\n')
        print(f'{label}: start {index+1}/{len(starts)}, cost={result.cost:.6g}, nfev={result.nfev}, {time.monotonic()-tick:.1f}s', flush=True)
    best = min(candidates, key=lambda c: c['cost'])
    return np.asarray(best['theta']), best


def settling(times, values, after):
    tail = times >= times[-1]-1000.
    final = float(np.median(values[tail]))
    if np.max(np.abs(values[tail]-final)) > .002:
        return None, final
    remains = np.logical_and.accumulate((np.abs(values-final) <= .002)[::-1])[::-1]
    idx = np.flatnonzero(remains & (times >= after) & (times <= times[-1]-1000.))
    return (float(times[idx[0]]) if len(idx) else None), final


def reversal(times, values, after, reference):
    st, final = settling(times, values, after)
    success = st is not None and reference*final < 0 and abs(final+reference) <= .05 and min(abs(reference), abs(final)) >= .02
    return bool(success), st, final


def evaluate_case(case, theta, step=.5):
    y = integrate(theta, case['times'], case['pulses'], np.zeros(2), step)
    stop = case['pulses'][:, 1].max()
    rows = []
    for j, field in enumerate(['rac', 'mca']):
        st, final = settling(case['times'], y[:, j], stop)
        refst, reffinal = settling(case['times'], case['values'][:, j], stop)
        rows.append(dict(panel=case['panel'], task_id=case['task_id'], variant=case['variant'],
                         length_um=case['length'], field=field,
                         rmse=float(np.sqrt(np.mean((y[:, j]-case['values'][:, j])**2))),
                         max_error=float(np.max(np.abs(y[:, j]-case['values'][:, j]))),
                         final_error=abs(final-reffinal), predicted_final=final, reference_final=reffinal,
                         predicted_settling_s=st, reference_settling_s=refst,
                         physical_coordinate_valid=bool(np.max(np.abs(y[:, j])) <= .5)))
    return y, rows


def plot_cases(selected, theta, path, title, countercue=False):
    fig, axes = plt.subplots(len(selected), 2, figsize=(11, 2.9*len(selected)), squeeze=False, layout='constrained')
    for i, case in enumerate(selected):
        y, _ = evaluate_case(case, theta)
        shift = case['pulses'][-1, 0] if countercue else 0.
        for j, field in enumerate(['Rac', 'Bound MCA']):
            ax = axes[i, j]
            sign = 1 if j == 0 else -1
            ax.plot(case['times']-shift, .5+sign*case['values'][:, j], color='black', lw=1.4, label='De Belly')
            ax.plot(case['times']-shift, .5+sign*y[:, j], color='#d55e00', lw=1.3, label='Reduced model')
            for start, stop, taper, amplitude in case['pulses']:
                ax.axvspan(start-shift, stop-shift, color='#e8be54', alpha=.2)
            ax.axhline(.5, lw=.5, color='gray', ls=':')
            ax.set_ylabel(f'{field} x50/L')
            ax.set_xlabel('Time from countercue (s)' if countercue else 'Time (s)')
            if countercue:
                duration = case['pulses'][-1, 1]-shift
                ax.set_xlim(-1000, 5000)
                ax.set_title(f'{duration:g}-s countercue; L = {case["length"]:g} um')
            else:
                ax.set_title(f'Panel {case["panel"]}; L = {case["length"]:g} um; {case["variant"]}')
            ax.grid(alpha=.15)
            ax.legend(fontsize=8)
    fig.suptitle(title)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def run_sweeps(manifests, model_map, base_theta, out):
    reference = pd.read_csv(DATA/'figure2_v2_countercue_repolarization_audit.csv').set_index('task_id')
    rows = []
    tasks = [t for t in manifests['E'].values() if t['assay'] == 'countercue_threshold_sweep']
    for i, task in enumerate(tasks):
        length = float(task['simulation']['length'])
        variant = task['coupling_variant']
        # Frozen base transfer covers every length and all three feedback settings.
        models = [('frozen_reference', base_theta)]
        if (length, variant) in model_map:
            models.append(('condition_calibrated', model_map[length, variant]))
        for model_name, theta in models:
            end = task['simulation']['end_time']
            times = np.unique(np.r_[np.arange(0., end, 5.), end])
            pulses = pulses_from_task(task)
            y = integrate(theta, times, pulses, np.zeros(2))
            audit = task['polarity_reversal_audit']
            pre_index = max(0, np.searchsorted(times, audit['opposite_start'], side='left')-1)
            # The reference-conditioned score uses the exact retained target amplitude;
            # an additional self-conditioned score diagnoses amplitude mismatch.
            ref = reference.loc[task['task_id']]
            row = dict(model=model_name, task_id=task['task_id'], length_um=length, variant=variant,
                       duration_s=task['countercue_duration_s'])
            for j, field in enumerate(['rac', 'mca']):
                sign = 1 if j == 0 else -1
                refpre = sign*(audit[f'pre_{field}_reference']-.5)
                success, st, final = reversal(times, y[:, j], audit['opposite_stop'], refpre)
                selfsuccess, _, _ = reversal(times, y[:, j], audit['opposite_stop'], y[pre_index, j])
                row.update({field+'_predicted_success': success, field+'_self_conditioned_success': selfsuccess,
                            field+'_reference_success': bool(ref[field+'_reversal_success']),
                            field+'_predicted_final': final, field+'_predicted_settling_s': st,
                            field+'_pre_error': float(abs(y[pre_index, j]-refpre))})
            rows.append(row)
        if (i+1) % 500 == 0:
            print(f'Reversal sweep {i+1}/{len(tasks)}', flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(out/'reversal_sweep.csv', index=False)
    thresholds = []
    for (model, length, variant), frame in result.groupby(['model', 'length_um', 'variant']):
        for field in ['rac', 'mca']:
            a = frame.loc[frame[field+'_predicted_success'], 'duration_s']
            b = frame.loc[frame[field+'_reference_success'], 'duration_s']
            thresholds.append(dict(model=model, length_um=length, variant=variant, field=field,
                                   predicted_min_s=float(a.min()) if len(a) else np.nan,
                                   reference_min_s=float(b.min()) if len(b) else np.nan,
                                   agreement=float(np.mean(frame[field+'_predicted_success'] == frame[field+'_reference_success']))))
    threshold = pd.DataFrame(thresholds)
    threshold.to_csv(out/'reversal_thresholds.csv', index=False)
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True, sharey=True, layout='constrained')
    for i, variant in enumerate(['beta_0x', 'beta_1x', 'beta_2x']):
        for j, field in enumerate(['rac', 'mca']):
            ax = axes[i, j]
            subset = threshold[(threshold.variant == variant) & (threshold.field == field)]
            frozen = subset[subset.model == 'frozen_reference'].sort_values('length_um')
            calibrated = subset[subset.model == 'condition_calibrated'].sort_values('length_um')
            ax.plot(frozen.length_um, frozen.reference_min_s, 'ko-', ms=3, label='De Belly')
            ax.plot(frozen.length_um, frozen.predicted_min_s, 's--', color='#d55e00', ms=3, label='Frozen 29-um fit')
            ax.scatter(calibrated.length_um, calibrated.predicted_min_s, marker='D', s=30, color='#0072b2', label='Condition fit')
            # No observed success is explicitly represented above the sweep range.
            for frame, col, marker, color in [(frozen, 'reference_min_s', 'x', 'black'), (frozen, 'predicted_min_s', '+', '#d55e00'), (calibrated, 'predicted_min_s', 'D', '#0072b2')]:
                missing = frame[frame[col].isna()]
                ax.scatter(missing.length_um, np.full(len(missing), 1080), marker=marker, s=20, color=color)
            ax.set_title(f'{variant}: {field}')
            ax.set_xlabel('Length (um)')
            ax.set_ylabel('Minimum successful cue (s)')
            ax.set_yticks([0, 250, 500, 750, 1000, 1080], ['0', '250', '500', '750', '1000', 'No success'])
            ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle('Held-out reversal sweep: 50 durations from 10 to 1000 s\nNo success includes unpolarized, unsettled, or failed target matching')
    fig.savefig(out/'reversal_thresholds.png', dpi=170)
    plt.close(fig)
    return threshold


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE/'results_ungated')
    parser.add_argument('--starts', type=int, default=8)
    parser.add_argument('--max-nfev', type=int, default=180)
    parser.add_argument('--seed', type=int, default=20260924)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    tick = time.monotonic()
    cases, manifests, paths = load_data()
    train = training_cases(cases, 29., 'beta_1x')
    assert len(train) == 2
    rng = np.random.default_rng(args.seed)
    initial = np.array([25., 1000., .003, .05, .1, .1, .008])
    starts = [initial]
    for _ in range(args.starts-1):
        z = encode(initial) + rng.normal(0., 1.2, 7)
        starts.append(decode(np.clip(z, encode(LOW)+.01, encode(HIGH)-.01)))
    basepath = out/'fit_reference.json'
    if args.resume and basepath.exists():
        best = min(json.loads(basepath.read_text()), key=lambda r: r['cost'])
        theta = np.array(best['theta'])
    else:
        theta, best = fit(train, starts, args.max_nfev, out, 'reference')
    print('Reference parameters:', parameter_dict(theta), flush=True)
    plot_cases(train, theta, out/'calibration.png', 'Calibration only: shared parameters fitted to Rac and MCA')
    heldout = [c for c in cases if c['panel'] == 'E' and c['length'] == 29.]
    plot_cases(heldout, theta, out/'heldout_countercues.png', 'Held-out countercues: no refitting', True)
    # Step-halving against a tolerance far smaller than the figure's 0.002 criterion.
    errors = []
    for c in train + heldout:
        coarse = integrate(theta, c['times'], c['pulses'], np.zeros(2), .5)
        fine = integrate(theta, c['times'], c['pulses'], np.zeros(2), .25)
        errors.append(dict(task_id=c['task_id'], max_abs_difference=float(np.max(np.abs(coarse-fine)))))
    (out/'step_check.json').write_text(json.dumps(errors, indent=2)+'\n')
    if max(e['max_abs_difference'] for e in errors) > .0002:
        raise RuntimeError('Step-halving discrepancy exceeds 0.0002; refine integration')
    metrics = []
    for case in cases:
        _, rows = evaluate_case(case, theta)
        for row in rows:
            row['model'] = 'frozen_reference'
            row['used_for_calibration'] = case['task_id'] in {c['task_id'] for c in train} and case['panel'] in ['A', 'C']
        metrics.extend(rows)
    # Limited conditional refits quantify whether changing parameters fixes transfer.
    # No response/sweep data enter these calibrations.
    model_map = {(29., 'beta_1x'): theta}
    parameter_rows = [dict(length_um=29., variant='beta_1x', cost=best['cost'], **parameter_dict(theta))]
    for length, variant in [(17., 'beta_1x'), (21., 'beta_1x'), (25., 'beta_1x'), (33., 'beta_1x'), (39., 'beta_1x'), (29., 'beta_0x'), (29., 'beta_2x')]:
        label = f'L{length:g}_{variant}'
        local = training_cases(cases, length, variant)
        path = out/f'fit_{label}.json'
        if args.resume and path.exists():
            localbest = min(json.loads(path.read_text()), key=lambda r: r['cost'])
            localtheta = np.array(localbest['theta'])
        else:
            localtheta, localbest = fit(local, [theta], args.max_nfev, out, label)
        model_map[length, variant] = localtheta
        parameter_rows.append(dict(length_um=length, variant=variant, cost=localbest['cost'], **parameter_dict(localtheta)))
        local_ids = {(c['panel'], c['task_id']) for c in local}
        for case in cases:
            if case['length'] == length and (case['variant'] == variant or (case['panel'], case['task_id']) in local_ids):
                _, rows = evaluate_case(case, localtheta)
                for row in rows:
                    row['model'] = 'condition_calibrated'
                    row['used_for_calibration'] = case['panel'] in ['A', 'C']
                metrics.extend(rows)
    metric = pd.DataFrame(metrics)
    metric.to_csv(out/'trajectory_metrics.csv', index=False)
    pd.DataFrame(parameter_rows).to_csv(out/'parameter_mapping.csv', index=False)
    thresholds = run_sweeps(manifests, model_map, theta, out)
    # Compare final polarity transfer, including independently calibrated conditions.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for j, field in enumerate(['rac', 'mca']):
        ax = axes[j]
        frozen = metric[(metric.panel == 'C') & (metric.variant == 'beta_1x') & (metric.field == field) & (metric.model == 'frozen_reference')].sort_values('length_um')
        local = metric[(metric.panel == 'C') & (metric.variant == 'beta_1x') & (metric.field == field) & (metric.model == 'condition_calibrated')].sort_values('length_um')
        sign = 1 if j == 0 else -1
        ax.plot(frozen.length_um, .5+sign*frozen.reference_final, 'ko-', label='De Belly')
        ax.plot(frozen.length_um, .5+sign*frozen.predicted_final, '--', color='#d55e00', label='Frozen 29-um fit')
        ax.scatter(local.length_um, .5+sign*local.predicted_final, color='#0072b2', marker='D', label='Length-specific fit')
        ax.set(xlabel='Length (um)', ylabel=f'{field} final x50/L', title='Panel C, beta 1x')
        ax.legend(fontsize=8)
    fig.savefig(out/'length_transfer.png', dpi=170)
    plt.close(fig)
    provenance = dict(generated_utc=datetime.now(timezone.utc).isoformat(), wall_seconds=time.monotonic()-tick,
                      seed=args.seed, train_length_um=29, train_task_ids=[c['task_id'] for c in train],
                      observable='p=Rac_x50/L-0.5; v=0.5-bound_MCA_x50/L',
                      coupling='ungated', noise='zero', solver='event-aligned RK4, max step 0.5 s',
                      optimizer='bounded multistart least_squares; seven log positive parameters; no gate',
                      bounds={name:[float(lo),float(hi)] for name,lo,hi in zip(NAMES,LOW,HIGH)},
                      versions=dict(numpy=np.__version__, scipy=scipy.__version__, numba=numba.__version__, pandas=pd.__version__),
                      source_hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                     [Path(__file__), *paths, DATA/'figure2_v2_polarity_trajectories.csv.gz', DATA/'figure2_v2_countercue_repolarization_audit.csv']})
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (out/'best_parameters.json').write_text(json.dumps(parameter_dict(theta), indent=2)+'\n')
    selected = metric[(metric.model == 'frozen_reference') & (metric.length_um == 29) & ((metric.panel == 'E') | metric.used_for_calibration)]
    print(selected[['panel','field','rmse','final_error']].to_string(index=False), flush=True)
    print(thresholds[(thresholds.length_um == 29) & (thresholds.variant == 'beta_1x')].to_string(index=False), flush=True)
    print(f'Completed in {time.monotonic()-tick:.1f}s. Results: {out}', flush=True)


if __name__ == '__main__':
    main()
