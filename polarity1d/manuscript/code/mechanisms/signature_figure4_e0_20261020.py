"""Frozen, whole-realization E0 analysis. This module performs no integrations."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone

import numpy as np

import signature_figure4_analysis_20261020 as A

VIEWS = {
    'main': (.5, 1., 1),
    'analysis_first_half': (.5, .75, 1),
    'analysis_second_half': (.75, 1., 1),
    'late_quarter': (.875, 1., 1),
    'cadence_x2': (.5, 1., 2),
    'cadence_x4': (.5, 1., 4),
    'cadence_x8': (.5, 1., 8),
}
REGIMES = ('single_front', 'multiple_fronts', 'low_contrast', 'unclassified')
BOOTSTRAPS = 2000


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def require_hash(path, expected):
    actual = digest(path)
    if actual != expected:
        raise ValueError(f'Immutable input hash changed: {path}: {actual} != {expected}')


def source_check(root):
    m = read_json(root / 'analysis_manifest.json')
    for path, expected in m['source_sha256'].items():
        require_hash(path, expected)
    require_hash(root / 'inventory.json', m['inventory_sha256'])
    require_hash(root / 'inputs/estimands.json', m['estimands_sha256'])
    require_hash(root / 'results/native_identity_verification.json', m['native_identity_sha256'])
    return m


def prepare(root):
    root = root.resolve()
    if (root / 'analysis_manifest.json').exists():
        raise FileExistsError('Analysis preparation is immutable; use its existing manifest')
    inv = read_json(root / 'inventory.json')
    verify = read_json(root / 'results/native_identity_verification.json')
    if not verify['passed'] or verify['inventory_sha256'] != digest(root / 'inventory.json'):
        raise ValueError('Complete native-file identity verification is required')
    if inv['counts'] != {'mean': 110, 'noise': 840, 'zero': 12}:
        raise ValueError('The E0 scientific inventory changed')
    keys = sorted({(r['model'], r['C_label'], r['f']) for r in inv['cases'] if r['kind'] == 'noise'})
    if len(keys) != 110:
        raise ValueError('All 110 independent noisy conditions are required')
    groups = []
    for model, C, f in keys:
        noise = [r for r in inv['cases'] if (r['model'], r['C_label'], r['f'], r['kind']) == (model, C, f, 'noise')]
        means = [r for r in inv['cases'] if (r['model'], r['C_label'], r['f'], r['kind']) == (model, C, f, 'mean')]
        zeros = [r for r in inv['cases'] if (r['model'], r['C_label'], r['kind']) == (model, C, 'zero')]
        if len(means) != 1 or len(zeros) != 1 or len(noise) not in (2, 12):
            raise ValueError(f'Incomplete noisy/control condition {(model, C, f)}')
        groups.append(dict(task=len(groups), model=model, C_label=C, f=f,
                           noise=sorted(noise, key=lambda r: r['seed']), mean=means[0], zero=zeros[0]))
    frozen = root / 'frozen_sources/analysis'
    frozen.mkdir(parents=True)
    for path in (Path(__file__), Path(A.__file__)):
        shutil.copy2(path, frozen / path.name)
    sources = {str(p.resolve()): digest(p) for p in frozen.glob('*.py')}
    manifest = dict(created_utc=stamp(), scope='E0 A1-A4/A7 only; no new integrations',
        inventory_sha256=digest(root / 'inventory.json'), estimands_sha256=digest(root / 'inputs/estimands.json'),
        native_identity_sha256=digest(root / 'results/native_identity_verification.json'),
        source_sha256=sources, bootstrap_resamples=BOOTSTRAPS, views=VIEWS, groups=groups)
    write_json(root / 'analysis_manifest.json', manifest)
    script = root / 'submit_analysis.slurm'
    script.write_text('\n'.join([
        '#!/bin/bash', '#SBATCH --job-name=f4_e0', '#SBATCH --partition=red', '#SBATCH --qos=small',
        '#SBATCH --cpus-per-task=1', '#SBATCH --mem=4G', '#SBATCH --time=01:00:00',
        '#SBATCH --array=0-109', f'#SBATCH --output={root}/logs/analysis_%A_%a.log',
        'set -euo pipefail', 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1',
        f'/home/4473331/.conda/envs/scaling_sandbox/bin/python -B {frozen / Path(__file__).name} worker --out {root} --task "$SLURM_ARRAY_TASK_ID"', ''
    ]))
    (root / 'logs').mkdir(exist_ok=True)
    (root / 'analysis_workers').mkdir(exist_ok=True)
    return dict(manifest_sha256=digest(root / 'analysis_manifest.json'), groups=len(groups), submission=str(script))


def load_case(row):
    for prefix in ('scalar', 'profile', 'basal_profile', 'source_status'):
        require_hash(row[f'{prefix}_path'], row[f'{prefix}_sha256'])
    with np.load(row['scalar_path'], allow_pickle=False) as z:
        time = z['time'].astype(float)
        P = z['P_real'].astype(float) + 1j * z['P_imag']
        C = np.broadcast_to(z['circumference'], time.shape).copy().astype(float)
        contrast = z['contrast'].astype(float)
        regions = z['current_half_contrast_regions'].astype(int)
        M = z['M1_real'] + 1j * z['M1_imag']
    if not np.allclose(P / P[0], M / M[0], rtol=1e-10, atol=1e-10):
        raise ValueError(f'P is inconsistent with the native first moment in {row["id"]}')
    if len(time) != len(P) or len(P) < 50 or np.any(C <= 0):
        raise ValueError('Malformed scalar trajectory')
    if not all(np.isfinite(a).all() for a in (time, P, C, contrast)):
        raise ValueError('Nonfinite native scalar samples')
    if not np.allclose(np.diff(time), row['scalar_cadence'], rtol=1e-9, atol=1e-8):
        raise ValueError('Scalar cadence differs from the inventory')
    if not math.isclose(time[-1], row['native_horizon'], rel_tol=1e-9, abs_tol=1e-8):
        raise ValueError('Incomplete native scalar horizon')
    with np.load(row['profile_path'], allow_pickle=False) as z:
        pt = z['time'].astype(float)
        primary = z['primary'].astype(float)
        pc = np.broadcast_to(z['circumference'], pt.shape).copy().astype(float)
        if 'material_positions' in z:
            q = z['material_positions'].astype(float)
        else:
            q = z['x'].astype(float) / float(pc[0])
        peak = q[np.argmax(primary, axis=1)]
        pcontrast = np.ptp(primary, axis=1)
    if not np.allclose(np.diff(pt), row['profile_cadence'], rtol=1e-9, atol=1e-8):
        raise ValueError('Profile cadence differs from the inventory')
    pi = np.rint((pt - time[0]) / row['scalar_cadence']).astype(int)
    if not np.allclose(time[pi], pt, rtol=0, atol=1e-7):
        raise ValueError('Profiles are not exactly aligned to their scalar trajectory')
    return dict(row=row, time=time, P=P, C=C, contrast=contrast, regions=regions,
                pt=pt, pc=pc, peak=peak, pcontrast=pcontrast, pi=pi)


def select(time, view, horizon):
    start, stop, step = VIEWS[view]
    return np.flatnonzero((time >= horizon * start - 1e-8) & (time <= horizon * stop + 1e-8))[::step]


def case_view(case, view):
    row = case['row']
    ix = select(case['time'], view, row['native_horizon'])
    time = case['time'][ix]
    cadence = float(np.diff(time)[0])
    P, C, contrast, regions = [case[k][ix] for k in ('P', 'C', 'contrast', 'regions')]
    valid = abs(P) >= .1
    single = valid & (regions == 1) & (contrast >= .1 * row['basal_contrast'])
    lag = (len(time) - 1) // 2
    cn, cc, possible = A.heading_cosine(P, valid, lag)
    mn, mc, _ = A.section_msd(A.angular_coordinate(P, valid), valid, lag)
    site = {}
    for name, direct in (('heading_arc', False), ('heading_direct', True)):
        site[name] = A.section_msd(A.physical_coordinate(P, C, single, direct=direct), single, lag)
    fractions, label = A.regime(contrast, regions, row['basal_contrast'])
    forced_fraction = float(np.mean((~valid) | (contrast < .1 * row['basal_contrast'])))
    metrics = dict(id=row['id'], model=row['model'], C_label=row['C_label'], f=row['f'],
        kind=row['kind'], seed=row['seed'], view=view, cadence=cadence, frames=len(ix),
        native_sigma=row.get('sigma',1. if row['kind']=='noise' else 0.),native_horizon=row['native_horizon'],
        source_experiment=row.get('origin','unspecified'),
        actual_C_mean=float(C.mean()), actual_C_min=float(C.min()), actual_C_max=float(C.max()),
        R_valid_fraction=float(valid.mean()), single_site_fraction=float(single.mean()),
        forced_union_fraction=forced_fraction, regime=label, basal_contrast=row['basal_contrast'])
    metrics.update({name + '_fraction': float(fractions[i]) for i, name in enumerate(REGIMES)})
    individual=A.persistence(A.ratio(cn,cc),A.ratio(cc,possible),cadence)
    metrics.update(per_case_tau_status=individual['status'],per_case_tau=individual['estimate'],
        per_case_tau_lower=individual['lower'],per_case_tau_upper=individual['upper'],
        per_case_tau_upper_unbounded=math.isinf(individual['upper']),per_case_tau_coverage=individual['coverage'])
    return dict(metrics=metrics, cn=cn, cc=cc, possible=possible, mn=mn, mc=mc, site=site, fractions=fractions)


def peak_view(case, view):
    row = case['row']
    ix = select(case['pt'], view, row['native_horizon'])
    pi = case['pi'][ix]
    valid = (abs(case['P'][pi]) >= .1) & (case['regions'][pi] == 1) & (case['pcontrast'][ix] >= .1 * row['basal_contrast'])
    peakP = np.exp(2j * np.pi * case['peak'][ix])
    pc = case['pc'][ix]
    coordinate = A.physical_coordinate(peakP, pc, valid)
    return (*A.section_msd(coordinate, valid, (len(ix) - 1) // 2), float(np.diff(case['pt'][ix])[0]))


def pooled_bootstrap(numbers, counts, weights):
    return A.ratio(np.sum(numbers, axis=0), np.sum(counts, axis=0)), A.ratio(weights @ numbers, weights @ counts)


def analyze_group(root, task):
    m = source_check(root)
    group = m['groups'][task]
    out = root / f'analysis_workers/{task:04d}'
    out.mkdir(exist_ok=True)
    if (out / 'result.json').exists():
        raise FileExistsError(f'Immutable completed worker already exists: {out}')
    cases = [load_case(r) for r in group['noise']]
    controls = [load_case(group[k]) for k in ('mean', 'zero')]
    key = group.get('bootstrap_pairing_key',f"{group['model']}:{group['C_label']}:{group['f']}")
    seed = 2026102004 + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
    weights = A.bootstrap_weights(len(cases), BOOTSTRAPS, seed)
    arrays, results, case_rows = {}, {}, []
    for view in VIEWS:
        values = [case_view(c, view) for c in cases]
        control_values = [case_view(c, view) for c in controls]
        case_rows.extend(v['metrics'] for v in values + control_values)
        cadence = values[0]['metrics']['cadence']
        nums, counts = np.array([v['cn'] for v in values]), np.array([v['cc'] for v in values])
        possible = np.array([v['possible'] for v in values])
        curve, boot = pooled_bootstrap(nums, counts, weights)
        cov = A.ratio(counts.sum(axis=0), possible.sum(axis=0))
        cin = A.finite_quantiles(boot)
        persistence = A.persistence(curve, cov, cadence)
        lower, upper, code = A.persistence_draws(boot, cadence)
        tci = A.interval_quantiles(lower, upper)
        angular_nums = np.array([v['mn'] for v in values])
        angular_counts = np.array([v['mc'] for v in values])
        msd, msdboot = pooled_bootstrap(angular_nums, angular_counts, weights)
        msdci = A.finite_quantiles(msdboot)
        dtheta = float(A.signed_diffusivity(msd, cadence))
        db = A.signed_diffusivity(msdboot, cadence)
        dci = A.finite_quantiles(db)
        C = np.array([v['metrics']['actual_C_mean'] for v in values])
        cb = weights @ C / len(cases)
        Cmean = float(C.mean())
        dx = dtheta * (Cmean / (2 * np.pi)) ** 2
        dxb = db * (cb / (2 * np.pi)) ** 2
        dxci = A.finite_quantiles(dxb)
        fracs = np.array([v['fractions'] for v in values])
        pooled_fracs = fracs.mean(axis=0)
        fboot = weights @ fracs / len(cases)
        fci = A.finite_quantiles(fboot)
        label = REGIMES[int(np.argmax(pooled_fracs))] if max(pooled_fracs) > .5 else 'mixed'
        mean = control_values[0]['metrics']
        metric = dict(model=group['model'], C_label=group['C_label'], f=group['f'], view=view,
            experiment=m.get('stage','E0'),sigma=group.get('sigma',1.),native_horizon=cases[0]['row']['native_horizon'],
            native_profile_cadence=cases[0]['row']['profile_cadence'],
            source_experiments=sorted({c['row'].get('origin','unspecified') for c in cases}),
            R_threshold=.1,relative_contrast_threshold=.1,raw_region_contour=.5,replication_unit='whole independent noisy realization',
            n_seeds=len(cases), bootstrap_resamples=BOOTSTRAPS, limited_two_seed_ensemble=len(cases) == 2,
            cadence=cadence, window_start_fraction=VIEWS[view][0], window_end_fraction=VIEWS[view][1],
            actual_C_mean=Cmean, tau_status=persistence['status'], tau=persistence['estimate'],
            tau_lower=persistence['lower'], tau_upper=persistence['upper'],
            tau_ci_lower=tci[0], tau_ci_upper=tci[1], tau_coverage=persistence['coverage'],
            tau_crossing_bracket=persistence['crossing_bracket'],
            tau_bootstrap_resolved_fraction=float(np.mean(code == 0)),
            tau_bootstrap_right_censored_fraction=float(np.mean(code == 2)),
            exponential_comparison=A.exponential_comparison(curve, counts.sum(axis=0), cadence),
            Dtheta=dtheta, Dtheta_ci_lower=float(dci[0]), Dtheta_ci_upper=float(dci[1]),
            Dtheta_status='positive_short_lag_slope' if dtheta > 0 else 'unresolved_nonpositive_or_missing_slope',
            Dtheta_valid_bootstrap_fraction=float(np.isfinite(db).mean()),
            angular_MSD_coverage_min_fit=float(A.ratio(angular_counts.sum(axis=0), possible.sum(axis=0))[1:10].min()),
            Dx=dx, Dx_ci_lower=float(dxci[0]), Dx_ci_upper=float(dxci[1]), regime=label,
            forced=mean['forced_union_fraction'] > .5, mean_forced_union_fraction=mean['forced_union_fraction'],
            forced_E1='unavailable_until_sigma_sweep', R_valid_fraction=float(np.mean([v['metrics']['R_valid_fraction'] for v in values])))
        for j, name in enumerate(REGIMES):
            metric.update({name + '_fraction': float(pooled_fracs[j]), name + '_ci_lower': float(fci[0, j]), name + '_ci_upper': float(fci[1, j])})
        prefix = view + '__'
        for name, value in dict(lag=np.arange(1, len(curve) + 1) * cadence, C=curve, C_ci_lower=cin[0], C_ci_upper=cin[1],
                coverage=cov, eligible_pairs=counts.sum(axis=0), possible_pairs=possible.sum(axis=0),
                angular_MSD=msd, angular_MSD_ci_lower=msdci[0], angular_MSD_ci_upper=msdci[1],
                angular_MSD_coverage=A.ratio(angular_counts.sum(axis=0), possible.sum(axis=0)),
                tau_boot_lower=lower, tau_boot_upper=upper, tau_boot_code=code,
                Dtheta_boot=db, Dx_boot=dxb, Cmean_boot=cb, regime_boot=fboot).items():
            arrays[prefix + name] = value
        for method in ('heading_arc', 'heading_direct', 'profile_peak'):
            if method == 'profile_peak':
                tracks = [peak_view(c, view) for c in cases]
                sn, sc, sp = [np.array([t[j] for t in tracks]) for j in range(3)]
                sdt = tracks[0][3]
            else:
                sn, sc, sp = [np.array([v['site'][method][j] for v in values]) for j in range(3)]
                sdt = cadence
            smsd, sb = pooled_bootstrap(sn, sc, weights)
            sci = A.finite_quantiles(sb)
            for name, value in dict(lag=np.arange(1, len(smsd) + 1) * sdt, MSD=smsd,
                    MSD_ci_lower=sci[0], MSD_ci_upper=sci[1], coverage=A.ratio(sc.sum(axis=0), sp.sum(axis=0)),
                    eligible_pairs=sc.sum(axis=0), possible_pairs=sp.sum(axis=0)).items():
                arrays[prefix + method + '__' + name] = value
        results[view] = metric
    np.savez_compressed(out / 'arrays.npz', **arrays)
    write_json(out / 'result.json', dict(generated_utc=stamp(), passed=True, task=task,
        analysis_manifest_sha256=digest(root / 'analysis_manifest.json'), metrics=results, case_metrics=case_rows,
        arrays_sha256=digest(out / 'arrays.npz'), input_ids=[c['row']['id'] for c in cases + controls],
        source_sha256=m['source_sha256'], no_native_integrations=True))
    return dict(task=task, passed=True, conditions=1, views=len(VIEWS), noisy_seeds=len(cases))


def csv_write(path, rows):
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(clean(v), separators=(',', ':')) if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})


def nfloat(value, infinity=False):
    return math.inf if value is None and infinity else math.nan if value is None else float(value)


def condition_row(metric):
    row = dict(metric)
    exp = row.pop('exponential_comparison')
    row.update(exponential_tau=exp.get('tau'), exponential_status=exp['status'], exponential_weighted_squared_error=exp.get('weighted_squared_error'))
    row['tau_upper_unbounded'] = row['tau_upper'] is None
    row['tau_ci_upper_unbounded'] = row['tau_ci_upper'] is None
    return row


def ratio_row(small, large, sa, la, view):
    p = dict(model=small['model'], f=small['f'], view=view, small_C_label=small['C_label'], large_C_label=large['C_label'],
             n_seeds_small=small['n_seeds'], n_seeds_large=large['n_seeds'], forced=small['forced'] or large['forced'])
    lower, upper = A.interval_ratio(nfloat(large['tau_lower']), nfloat(large['tau_upper'], True), nfloat(small['tau_lower']), nfloat(small['tau_upper'], True))
    blo, bhi = A.interval_ratio(la[view + '__tau_boot_lower'], la[view + '__tau_boot_upper'], sa[view + '__tau_boot_lower'], sa[view + '__tau_boot_upper'])
    ci = A.interval_quantiles(blo, bhi)
    geom = (large['actual_C_mean'] / small['actual_C_mean']) ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        mobility = sa[view + '__Dtheta_boot'] / la[view + '__Dtheta_boot']
    mobility[(sa[view + '__Dtheta_boot'] <= 0) | (la[view + '__Dtheta_boot'] <= 0)] = np.nan
    mci = A.finite_quantiles(mobility)
    mobility_point = small['Dtheta'] / large['Dtheta'] if small['Dtheta'] is not None and large['Dtheta'] is not None and min(small['Dtheta'], large['Dtheta']) > 0 else math.nan
    dxratio = A.ratio(la[view + '__Dx_boot'], sa[view + '__Dx_boot'])
    dxratio[(sa[view + '__Dx_boot'] <= 0) | (la[view + '__Dx_boot'] <= 0)] = np.nan
    dxci = A.finite_quantiles(dxratio)
    point = large['tau'] / small['tau'] if large['tau'] is not None and small['tau'] is not None else None
    evaluable = min(small['tau_coverage'], large['tau_coverage']) >= .5 and all(r['tau_status'] in ('resolved', 'right_censored') for r in (small, large))
    geometry_boot=(la[view + '__Cmean_boot']/sa[view + '__Cmean_boot'])**2
    excess_ci=A.interval_quantiles(blo/geometry_boot,bhi/geometry_boot)
    p.update(tau_ratio=point, tau_ratio_lower=float(lower), tau_ratio_upper=float(upper), tau_ratio_ci_lower=ci[0], tau_ratio_ci_upper=ci[1],
        tau_ratio_upper_unbounded=math.isinf(float(upper)), tau_ratio_ci_upper_unbounded=math.isinf(ci[1]),
        evaluable=evaluable, coverage_small=small['tau_coverage'], coverage_large=large['tau_coverage'],
        tau_status_small=small['tau_status'], tau_status_large=large['tau_status'],
        Dtheta_mobility_proxy=mobility_point, Dtheta_mobility_proxy_ci_lower=float(mci[0]), Dtheta_mobility_proxy_ci_upper=float(mci[1]),
        mobility_proxy_valid_bootstrap_fraction=float(np.isfinite(mobility).mean()),
        reference_geometry_ratio=9., actual_geometry_ratio=geom,
        excess_persistence=point / geom if point is not None else None,
        excess_persistence_ci_lower=excess_ci[0], excess_persistence_ci_upper=excess_ci[1],
        Dx_large_over_small=large['Dx'] / small['Dx'] if small['Dx'] is not None and large['Dx'] is not None and min(small['Dx'], large['Dx']) > 0 else math.nan,
        Dx_ratio_ci_lower=float(dxci[0]), Dx_ratio_ci_upper=float(dxci[1]), large_multiple_front_fraction=large['multiple_fronts_fraction'])
    p['geometric_qualifier'] = ('more_than_geometry' if ci[0] > geom or dxci[1] < 1 else
        'less_than_geometry' if dxci[0] > 1 else 'compatible_with_flat_Dx' if dxci[0] <= 1 <= dxci[1] else 'inconclusive')
    return p


def fhalf_row(metrics, arrays, model, C, view, inclusive, required_weak_f=.0001):
    available = sorted([x for x in metrics if x['model'] == model and x['C_label'] == C and x['view'] == view], key=lambda x: x['f'])
    selected=[x for x in available if inclusive or not x['forced']]
    if len(selected)<2 or selected[0]['f']!=available[0]['f'] or available[0]['f']!=required_weak_f:
        return dict(model=model,C_label=C,view=view,forced_inclusive=inclusive,weakest_f=available[0]['f'],
            weakest_n_seeds=available[0]['n_seeds'],fhalf_status='unavailable_missing_declared_weak_input' if available[0]['f']!=required_weak_f else 'unavailable_after_forced_exclusion',fhalf=None,
            fhalf_lower=required_weak_f,fhalf_upper=math.inf,fhalf_ci_lower=required_weak_f,fhalf_ci_upper=math.inf,
            fhalf_upper_unbounded=True,fhalf_ci_upper_unbounded=True,weakest_tau_status=available[0]['tau_status'],
            fhalf_bootstrap_resolved_fraction=0.,excluded_forced_strengths=[x['f'] for x in available if x['forced']])
    f = np.array([x['f'] for x in selected])
    point = A.conservative_fhalf(f, [nfloat(x['tau_lower']) for x in selected], [nfloat(x['tau_upper'], True) for x in selected])
    los = np.array([arrays[(model, C, x['f'])][view + '__tau_boot_lower'] for x in selected])
    his = np.array([arrays[(model, C, x['f'])][view + '__tau_boot_upper'] for x in selected])
    draws = [A.conservative_fhalf(f, los[:, b], his[:, b]) for b in range(BOOTSTRAPS)]
    ci = A.interval_quantiles([d['lower'] for d in draws], [d['upper'] for d in draws])
    return dict(model=model, C_label=C, view=view, forced_inclusive=inclusive, weakest_f=float(f[0]), weakest_n_seeds=selected[0]['n_seeds'],
        fhalf_status=point['status'], fhalf=point['estimate'], fhalf_lower=point['lower'], fhalf_upper=point['upper'],
        fhalf_ci_lower=ci[0], fhalf_ci_upper=ci[1], fhalf_upper_unbounded=math.isinf(point['upper']),
        fhalf_ci_upper_unbounded=math.isinf(ci[1]), weakest_tau_status=selected[0]['tau_status'],
        fhalf_bootstrap_resolved_fraction=float(np.mean([d['status'] == 'resolved' for d in draws])),
        excluded_forced_strengths=[] if inclusive else [x['f'] for x in available if x['forced']])


def decisions(ratios, fhalfs):
    rows = []
    for model in sorted({x['model'] for x in ratios}):
        for view in VIEWS:
            for inclusive in (False, True):
                candidates = [r for r in ratios if r['model'] == model and r['view'] == view and (inclusive or not r['forced'])]
                evaluated = [r for r in candidates if r['evaluable']]
                gains = [r for r in evaluated if r['tau_ratio_ci_lower'] > 1]
                reversals = [r for r in evaluated if r['tau_ratio_ci_upper'] < 1]
                smallC = min(x['C_label'] for x in fhalfs if x['model'] == model)
                largeC = max(x['C_label'] for x in fhalfs if x['model'] == model)
                fs = next(x for x in fhalfs if (x['model'], x['C_label'], x['view'],x['forced_inclusive']) == (model, smallC, view,inclusive))
                fl = next(x for x in fhalfs if (x['model'], x['C_label'], x['view'],x['forced_inclusive']) == (model, largeC, view,inclusive))
                tolerance_forced = any(r['forced'] for r in ratios if r['model'] == model and r['view'] == view)
                tolerance_gain = fl['fhalf_ci_lower'] > fs['fhalf_ci_upper']
                h1 = 'pass_provisional' if evaluated and len(gains) == len(evaluated) and tolerance_gain else 'fail' if reversals else 'inconclusive'
                peak_drop = []
                for later in evaluated:
                    if not (later['tau_ratio_ci_lower'] <= 1 <= later['tau_ratio_ci_upper']):
                        continue
                    for earlier in evaluated:
                        if earlier['f'] >= later['f'] or earlier['tau_ratio'] is None or later['tau_ratio'] is None:
                            continue
                        span = [r for r in candidates if earlier['f'] <= r['f'] <= later['f']]
                        if earlier['tau_ratio'] >= 2 * later['tau_ratio'] and all(r['large_multiple_front_fraction'] > .5 for r in span):
                            peak_drop.append((earlier['f'], later['f']))
                h2 = 'pass_model_criterion_provisional' if reversals or peak_drop else 'inconclusive'
                rows.append(dict(model=model, view=view, forced_inclusive=inclusive, evaluable_conditions=len(evaluated),
                    excluded_forced_conditions=sum(r['forced'] for r in ratios if r['model'] == model and r['view'] == view) if not inclusive else 0,
                    H1=h1, gain_CI_excludes_one_count=len(gains), reversal_CI_excludes_one_count=len(reversals),
                    evaluable_f=[r['f'] for r in evaluated], reversal_f=[r['f'] for r in reversals],
                    fhalf_increase_nonoverlap=tolerance_gain, fhalf_small_ci=[fs['fhalf_ci_lower'], fs['fhalf_ci_upper']],
                    fhalf_large_ci=[fl['fhalf_ci_lower'], fl['fhalf_ci_upper']], fhalf_forced_contribution=tolerance_forced,
                    H2=h2, peak_to_unity_multiple_front_ranges=peak_drop,
                    H2_overall='pending_A6_and_transition_refinement', H3='pending_E1_and_A5_fingerprints',
                    caveat='E0 only; empirical censoring and limited two-seed outer strengths retained'))
    return rows


def collect(root):
    m = source_check(root)
    collector_freeze=read_json(root/'collector_freeze_v3.json')
    for path,expected in collector_freeze['source_sha256'].items():
        require_hash(path,expected)
    if str(Path(__file__).resolve()) not in collector_freeze['source_sha256']:
        raise ValueError('The actual export collector must belong to its immutable source closure')
    out = root / 'results/analysis'
    if out.exists():
        raise FileExistsError('Do not overwrite an immutable completed analysis export')
    metrics, cases, heading, site, arrays, receipts = [], {}, [], [], {}, []
    for group in m['groups']:
        worker = root / f"analysis_workers/{group['task']:04d}"
        r = read_json(worker / 'result.json')
        if not r['passed'] or r['task'] != group['task'] or r['analysis_manifest_sha256'] != digest(root / 'analysis_manifest.json'):
            raise ValueError('A worker is incomplete or belongs to another immutable manifest')
        require_hash(worker / 'arrays.npz', r['arrays_sha256'])
        key = group['model'], group['C_label'], group['f']
        with np.load(worker / 'arrays.npz', allow_pickle=False) as storage:
            z={k:storage[k] for k in storage.files}
            arrays[key] = {k: z[k] for k in z if '__tau_boot_' in k or k.endswith(('__Dtheta_boot', '__Dx_boot', '__Cmean_boot'))}
            for view, metric in r['metrics'].items():
                metrics.append(metric)
                prefix = view + '__'
                for j, lag in enumerate(z[prefix + 'lag']):
                    row = dict(model=key[0], C_label=key[1], f=key[2], view=view, lag=float(lag), forced=metric['forced'])
                    for name in ('C', 'C_ci_lower', 'C_ci_upper', 'coverage', 'eligible_pairs', 'possible_pairs', 'angular_MSD', 'angular_MSD_ci_lower', 'angular_MSD_ci_upper', 'angular_MSD_coverage'):
                        row[name] = z[prefix + name][j].item()
                    heading.append(row)
                for method in ('heading_arc', 'heading_direct', 'profile_peak'):
                    sp = prefix + method + '__'
                    for j, lag in enumerate(z[sp + 'lag']):
                        row = dict(model=key[0], C_label=key[1], f=key[2], view=view, method=method,
                            lag=float(lag), forced=metric['forced'], cadence=float(z[sp + 'lag'][0]),
                            profile_cadence_limited=method == 'profile_peak')
                        for name in ('MSD', 'MSD_ci_lower', 'MSD_ci_upper', 'coverage', 'eligible_pairs', 'possible_pairs'):
                            row[name] = z[sp + name][j].item()
                        site.append(row)
        for case in r['case_metrics']:
            casekey=(case['id'],case['view'])
            if casekey in cases and cases[casekey]!=case:
                raise ValueError('A reused zero control has conflicting metrics')
            cases[casekey]=case
        receipts.append(dict(task=group['task'], result_sha256=digest(worker / 'result.json'), arrays_sha256=r['arrays_sha256']))
    group_count=len(m['groups'])
    if len(metrics) != group_count * len(VIEWS):
        raise ValueError('Missing condition/view work is an error, not censoring')
    ratios, fhalfs = [], []
    for model in sorted({x['model'] for x in metrics}):
        sizes = sorted({x['C_label'] for x in metrics if x['model'] == model})
        for view in VIEWS:
            sm = {x['f']: x for x in metrics if (x['model'], x['C_label'], x['view']) == (model, sizes[0], view)}
            lg = {x['f']: x for x in metrics if (x['model'], x['C_label'], x['view']) == (model, sizes[-1], view)}
            for f in sorted(sm.keys() & lg.keys()):
                ratios.append(ratio_row(sm[f], lg[f], arrays[(model, sizes[0], f)], arrays[(model, sizes[-1], f)], view))
            for C in sizes:
                for inclusive in (False, True):
                    fhalfs.append(fhalf_row(metrics, arrays, model, C, view,inclusive,m.get('required_fhalf_baseline',.0001)))
    decision = decisions(ratios, fhalfs)
    out.mkdir(parents=True)
    census=m.get('census_counts',dict(noise=840,mean=110,zero=12))
    if len(cases)!=sum(census.values())*len(VIEWS):
        raise ValueError('Every independent native realization and control must appear exactly once per view')
    tables = dict(condition_metrics=[condition_row(x) for x in metrics], case_metrics=list(cases.values()),
        heading_curves=heading, site_curves=site, size_ratios=ratios, noise_tolerance=fhalfs, decisions=decision)
    for name, rows in tables.items():
        csv_write(out / (name + '.csv'), rows)
    outputs = {name + '.csv': digest(out / (name + '.csv')) for name in tables}
    status = dict(generated_utc=stamp(), passed=True, complete=True, declared_condition_groups_complete=True,
        condition_groups=group_count, noisy_realizations=census['noise'], mean_controls=census['mean'], zero_controls=census['zero'],
        inventory_sha256=m['inventory_sha256'], estimands_sha256=m['estimands_sha256'],
        native_identity_sha256=m['native_identity_sha256'], analysis_manifest_sha256=digest(root / 'analysis_manifest.json'),
        source_sha256=m['source_sha256'], collector_source_sha256=collector_freeze['source_sha256'],
        collector_freeze_sha256=digest(root/'collector_freeze_v3.json'),outputs=outputs, rows={k: len(v) for k, v in tables.items()},
        worker_receipts=receipts, scope=m.get('scope','E0 A1-A4/A7 interim; A5/A6/E1/E2/E3/E4 pending'),
        stage=m.get('stage','E0'),sigma=m.get('sigma',1.),
        horizon_extension='observed_E4' if m.get('stage')=='E4' else 'unavailable_in_this_arm', no_native_integrations=True)
    write_json(out / 'status.json', status)
    return dict(passed=True, groups=group_count, output=str(out), status_sha256=digest(out / 'status.json'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=('prepare', 'worker', 'collect'))
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=int)
    args = p.parse_args()
    if args.command == 'prepare':
        result = prepare(args.out)
    elif args.command == 'worker':
        result = analyze_group(args.out.resolve(), args.task)
    else:
        result = collect(args.out.resolve())
    print(json.dumps(clean(result), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
