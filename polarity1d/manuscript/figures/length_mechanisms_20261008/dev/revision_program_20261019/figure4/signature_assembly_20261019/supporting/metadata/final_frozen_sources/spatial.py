"""Read-only F4 window sensitivity and instantaneous periodic region widths.

No integration or event/establishment/lifetime classifier. Widths use saved
native primary profiles at their actual cadence and physical circumference.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.util
import json
import os
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'polarity1d/manuscript').is_dir())
MS = ROOT / 'polarity1d/manuscript'
RUN = MS / 'data/runs/revision_noise_strength_broad_20261019'
OUT = MS / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/postburn_signature_candidate_20261019'
COLLECTOR = RUN / 'results/monitor/collector_v2.py'
COLLECTOR_SHA = '79747ebe6d831571e302a1fec528b7c6f58b223d7c9d7a8ea827b22a134fdd83'
METRICS = ('multiregion_fraction', 'qualified_multiregion_fraction', 'weak_R_fraction',
           'low_contrast_fraction', 'mean_R', 'mean_A', 'mean_Q', 'mean_contrast_ratio',
           'mean_raw_q1', 'mean_q2', 'factorization_maximum_error')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, (np.floating, float)) and not np.isfinite(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def write(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')


def width_cells(mask):
    """Descending widths of circular connected True components, in cells."""
    n = len(mask)
    if not np.any(mask):
        return np.empty(0, dtype=int)
    if np.all(mask):
        return np.array([n], dtype=int)
    # Cut at a below-threshold cell so no True component can cross the new seam.
    seam = int(np.flatnonzero(~mask)[0])
    shifted = np.roll(mask, -seam)
    diff = np.diff(np.r_[False, shifted, False].astype(int))
    return np.sort(np.flatnonzero(diff == -1) - np.flatnonzero(diff == 1))[::-1]


def windows(T):
    return [('full', 0., T), ('first_half', 0., T/2), ('last_half', T/2, T),
            ('quarter_1', 0., T/4), ('quarter_2', T/4, T/2),
            ('quarter_3', T/2, 3*T/4), ('quarter_4', 3*T/4, T),
            ('tail_after_quarter_1', T/4, T)]


def profile_samples(path, scalar):
    with np.load(path) as z:
        t = z['time'].copy(); u = z['primary'].copy()
        C = np.broadcast_to(z['circumference'], len(t)).copy()
    contrast = np.ptp(u, axis=1)
    valid = contrast > 1e-12
    above = u >= u.min(axis=1)[:, None] + .5 * contrast[:, None]
    counts = np.sum(above & ~np.roll(above, 1, axis=1), axis=1)
    counts[above.all(axis=1)] = 1; counts[~valid] = 0
    ix = np.searchsorted(scalar['time'], t)
    np.testing.assert_allclose(scalar['time'][ix], t, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(counts, scalar['current_half_contrast_regions'][ix])
    largest = np.zeros(len(t)); second = np.full(len(t), np.nan)
    total = np.zeros(len(t))
    for k in np.flatnonzero(valid):
        widths = width_cells(above[k])
        assert len(widths) == counts[k]
        largest[k] = widths[0] / u.shape[1]
        total[k] = widths.sum() / u.shape[1]
        if len(widths) > 1:
            second[k] = widths[1] / u.shape[1]
    return dict(time=t, circumference=C, contrast_ratio=contrast/contrast[0],
                regions=counts, largest_fraction=largest, second_fraction=second,
                total_above_fraction=total, largest_physical=largest*C,
                second_physical=second*C, qualified=valid & (contrast/contrast[0] >= .1))


def summarize_spatial(samples, label, lo, hi):
    keep = (samples['time'] >= lo) & (samples['time'] <= hi)
    qualified = samples['qualified'][keep]
    counts = samples['regions'][keep]
    one = qualified & (counts == 1); multiple = qualified & (counts > 1)
    result = dict(window=label, first_time=lo, last_time=hi,
        profile_samples=int(keep.sum()), profile_cadence=float(np.diff(samples['time'])[0]),
        qualified_profile_samples=int(qualified.sum()),
        qualified_profile_fraction=float(qualified.mean()),
        qualified_one_region_fraction=float(one.mean()),
        qualified_two_region_fraction=float(np.mean(qualified & (counts == 2))),
        qualified_three_or_more_fraction=float(np.mean(qualified & (counts >= 3))),
        qualified_multiregion_fraction_profile_cadence=float(multiple.mean()),
        raw_multiregion_fraction_profile_cadence=float(np.mean(counts > 1)),
        qualified_multiple_profile_samples=int(multiple.sum()),
        maximum_regions=int(counts.max()), mean_raw_region_count=float(counts.mean()),
        mean_qualified_region_count=float(counts[qualified].mean()) if qualified.any() else np.nan)
    for key in ('largest_fraction', 'largest_physical', 'total_above_fraction'):
        values = samples[key][keep][qualified]
        result['mean_' + key] = float(values.mean()) if len(values) else np.nan
        for quantile in (.1, .5, .9):
            result[f'q{int(quantile*100):02d}_' + key] = float(np.quantile(values, quantile)) if len(values) else np.nan
    for key in ('second_fraction', 'second_physical'):
        values = samples[key][keep][multiple]
        result['mean_' + key] = float(values.mean()) if len(values) else np.nan
        for quantile in (.1, .5, .9):
            result[f'q{int(quantile*100):02d}_' + key] = float(np.quantile(values, quantile)) if len(values) else np.nan
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    parser.add_argument('--run', type=Path, default=RUN)
    parser.add_argument('--expected-cases', type=int, choices=(576,936), default=576)
    args = parser.parse_args(); out = args.out; run = args.run
    data = out / 'supporting/data'; metadata = out / 'supporting/metadata'
    start = time.monotonic()
    assert sha(COLLECTOR) == COLLECTOR_SHA
    spec = importlib.util.spec_from_file_location('accepted_scalar_descriptors', COLLECTOR)
    collector = importlib.util.module_from_spec(spec); spec.loader.exec_module(collector)
    cachepath = run / 'results/partial/case_cache.json'
    cache = json.loads(cachepath.read_text()); assert len(cache) == args.expected_cases
    status = json.loads((run / 'results/partial/status.json').read_text())
    assert status['producer'] == 'terminal'
    provenance = []; spatial = []; scalar = []; memory = []; samples = {}
    for index, (case_id, case) in enumerate(cache.items()):
        old = case['descriptors'][0]
        base = {key: old[key] for key in ('case_id', 'model', 'platform', 'size', 'fraction', 'kind', 'seed', 'legacy')}
        for label in ('trace', 'profile'):
            assert sha(case[label]) == case[label + '_sha256']
        r = collector.trace(case['trace'])
        T = float(r['time'][-1])
        lags = [5., 40., 100., 500., 2000., 8000.] if old['platform'] == 'chemical' else [.5, 4., 10., 50., 200., 800.]
        s = profile_samples(case['profile'], r)
        for label, lo, hi in windows(T):
            descriptor, mm = collector.describe(r, label, lo, hi, lags)
            scalar.append(dict(**base, **{key: descriptor[key] for key in
                ('window', 'first_time', 'last_time', 'observations', 'scalar_cadence') + METRICS}))
            memory.extend(dict(**base, **value) for value in mm)
            spatial.append(dict(**base, **summarize_spatial(s, label, lo, hi)))
        for key, array in s.items():
            samples[case_id + '__' + key] = array
        provenance.append(dict(**base, trace=case['trace'], trace_sha256=case['trace_sha256'],
            profile=case['profile'], profile_sha256=case['profile_sha256'],
            scalar_frames=len(r['time']), profile_frames=len(s['time'])))
        if index % 48 == 0:
            print(json.dumps(dict(completed_cases=index+1, elapsed_seconds=time.monotonic()-start)), flush=True)
    pd.DataFrame(provenance).to_csv(data / 'postburn_native_source_index.csv', index=False)
    spatial = pd.DataFrame(spatial); spatial.to_csv(data / 'postburn_spatial_case_windows.csv', index=False)
    scalar = pd.DataFrame(scalar); scalar.to_csv(data / 'postburn_scalar_case_windows.csv', index=False)
    memory = pd.DataFrame(memory); memory.to_csv(data / 'postburn_memory_case_windows.csv', index=False)
    np.savez_compressed(data / 'postburn_region_samples.npz', **samples)
    keys = ['model', 'size', 'fraction', 'kind', 'window']
    descriptor_keys = [c for c in spatial if c not in list(provenance[0]) + ['window']]
    groups = []
    for key, group in spatial.groupby(keys, sort=True):
        result = dict(zip(keys, key), seeds=len(group))
        for metric in descriptor_keys:
            values = group[metric].dropna()
            result[metric] = values.mean() if len(values) else np.nan
            result[metric + '_seed_min'] = values.min() if len(values) else np.nan
            result[metric + '_seed_max'] = values.max() if len(values) else np.nan
            result[metric + '_contributing_seeds'] = len(values)
        groups.append(result)
    pd.DataFrame(groups).to_csv(data / 'postburn_spatial_condition_windows.csv', index=False)
    result = dict(generated_at=now(), elapsed_seconds=time.monotonic()-start,
        source=str(Path(__file__).resolve().relative_to(ROOT)), source_sha256=sha(__file__),
        cache_sha256=sha(cachepath), accepted_scalar_source=str(COLLECTOR.relative_to(ROOT)),
        accepted_scalar_source_sha256=COLLECTOR_SHA, native_integrations=0,
        cases=len(cache), scalar_frames=sum(p['scalar_frames'] for p in provenance),
        profile_frames=sum(p['profile_frames'] for p in provenance),
        exact_native_archives_verified=2*len(cache), windows=[w[0] for w in windows(1)],
        region_definition='Periodic connected grid cells at or above current min+0.5*(max-min); flat contrast<=1e-12 gives zero; seam joined before width',
        amplitude_qualification='Current contrast>=0.1*initial basal contrast',
        width_weighting='One largest width per qualified frame; one second-largest per qualified multiregion frame; condition means average available seed descriptors equally; coverage/contributors exported',
        physical_width='number of cells/N times actual C(t); material fraction widths exported separately',
        no_birth_establishment_lifetime_or_event_rate=True,
        stationarity_test='Finite-window sensitivity only; no equilibrium certification',
        outputs={str(p.relative_to(out)):sha(p) for p in data.glob('postburn_*') if p.is_file()})
    write(metadata / 'SPATIAL_WINDOW_ANALYSIS.json', result)
    print(json.dumps(clean(result), indent=2), flush=True)


if __name__ == '__main__':
    main()
