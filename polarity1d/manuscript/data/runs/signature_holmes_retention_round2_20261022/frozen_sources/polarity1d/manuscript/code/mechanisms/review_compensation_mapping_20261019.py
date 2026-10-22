"""Independent read-only replay of Figure 2 background-aware readouts.

Runs no solver and imports no producer measurement routines. Native dynamics
and full-state archive identities are covered by the preceding native audits.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import os
import time
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

ROOT = Path('/home/4473331/projects/scaling_sandbox/polarity1d/manuscript')
OUT = ROOT/'data/runs/review_compensation_mapping_20261019'
MODELS = ['wave_pinning', 'goryachev', 'holmes_model3', 'spring', 'debelly']
FLOORS = ['current_minimum', 'current_sample_quantile_0p05', 'frozen_initial_minimum']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def table(path):
    return list(csv.DictReader(Path(path).open()))


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp.json')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def entry(times, fraction, loss):
    good = fraction <= 1-loss
    hits = np.flatnonzero(good)
    if len(hits) == 0:
        return dict(status='no_sampled_entry_through_finite_horizon',
                    lower=float(times[-1]), upper=None)
    first = int(hits[0])
    out = dict(status='already_present' if first == 0 else 'sampled_entry',
               lower=float(times[max(0, first-1)]), upper=float(times[first]))
    if not good[-1]:
        out.update(retained_status='not_retained_through_last_sample',
                   retained_lower=None, retained_upper=None)
    else:
        failures = np.flatnonzero(~good)
        onset = int(failures[-1]+1) if len(failures) else 0
        out.update(retained_status='sampled_retained_entry' if onset else 'present_all_samples',
                   retained_lower=float(times[max(0, onset-1)]), retained_upper=float(times[onset]))
    return out


def independent_readouts(native, case, region):
    variable = {'debelly':'rac', 'holmes_model3':'rac_active'}.get(case['model'], 'active')
    with np.load(native) as archive:
        times = archive['time'].copy()
        active = archive['field_'+variable].copy()
        lengths = archive['field_length'].copy() if case['model'] == 'spring' else np.ones(len(times))
    basal = region['initial_active']
    q = region['material_fraction']
    mask = region['initial_core_mask'].astype(bool)
    threshold = (float(basal.min())+float(basal.max()))/2
    last = next((i for i, u in enumerate(basal) if u < threshold), len(basal))
    edge = 1. if last == len(basal) else float(np.interp(threshold, [basal[last], basal[last-1]], [q[last], q[last-1]]))
    initial_length = float(region['initial_physical_length']) if case['model'] == 'spring' else 1.
    initial_physical_worker_weights = region['initial_weights']*mask

    if case['model'] == 'spring':
        def core(values):
            keep = q <= edge
            positions = q[keep]
            selected = values[..., keep]
            if positions[-1] != edge:
                right = int(np.searchsorted(q, edge))
                fraction = (edge-q[right-1])/(q[right]-q[right-1])
                boundary = values[..., right-1]*(1-fraction)+values[..., right]*fraction
                selected = np.concatenate([selected, boundary[..., None]], axis=-1)
                positions = np.append(positions, edge)
            return np.trapezoid(selected, positions, axis=-1)
    else:
        boundaries = np.linspace(0, case['length'], len(basal)+1)
        portions = np.maximum(0, np.minimum(boundaries[1:], edge*case['length'])-boundaries[:-1])
        def core(values):
            return np.sum(values*portions, axis=-1)

    def worker(values):
        return np.sum(values*initial_physical_worker_weights, axis=-1)/initial_length

    floors = [np.min(active, axis=1), np.percentile(active, 5, axis=1), np.full(len(times), basal.min())]
    result = dict(time=times, current_global_minimum=floors[0],
                  current_global_sample_quantile_0p05=floors[1],
                  global_active_contrast=np.max(active, axis=1)-np.min(active, axis=1),
                  fixed_initial_core_peak=np.max(active[:, mask], axis=1))
    for roi, integral in [('organization_core', core), ('worker_core', worker)]:
        raw = integral(active)*lengths
        result[roi+'_raw_amount'] = raw
        result[roi+'_raw_fraction'] = raw/(float(integral(basal))*initial_length)
        result[roi+'_extent'] = np.full(len(times), float(integral(np.ones_like(basal))))*lengths
        for rule, floor in zip(FLOORS, floors):
            concentration = np.maximum(active-floor[:, None], 0)
            initial_floor = np.percentile(basal, 5) if rule == FLOORS[1] else basal.min()
            denominator = float(integral(np.maximum(basal-initial_floor, 0)))*initial_length
            amount = integral(concentration)*lengths
            result[roi+'_'+rule+'_excess_amount'] = amount
            result[roi+'_'+rule+'_excess_fraction'] = amount/denominator
    for rule, floor in zip(FLOORS, floors):
        peak = np.maximum(result['fixed_initial_core_peak']-floor, 0)
        initial_floor = np.percentile(basal, 5) if rule == FLOORS[1] else basal.min()
        result['worker_core_'+rule+'_peak_contrast'] = peak
        result['worker_core_'+rule+'_peak_contrast_fraction'] = peak/(basal[mask].max()-initial_floor)
    return result


def audit(task, limit=None):
    started = time.monotonic()
    label = ['baseline', 'strength'][task//5]
    model = MODELS[task%5]
    run = ROOT/f'data/runs/revision_compensation_excess_{label}_20261019'
    manifest = load(run/'manifest.json')
    result = load(run/'workers'/model/'result.json')
    assert result['manifest_sha256'] == sha(run/'manifest.json')
    events = table(run/'workers'/model/'sampled_loss_events.csv')
    endpoints = table(run/'workers'/model/'pulse_endpoints.csv')
    key = lambda r:(int(r['case_slot']), r['phase'], float(r['pulse_duration']))
    event_map = {}
    for row in events:
        extended = key(row)+(row['roi'], row['background_rule'], row['observable'], float(row['loss_fraction']))
        assert extended not in event_map
        event_map[extended] = row
    endpoint_map = {key(row):row for row in endpoints}
    checked = dict(phases=0, frames=0, arrays=0, values=0, endpoints=0, events=0)
    maximum = 0.
    hashes = {}

    def compare(actual, expected):
        nonlocal maximum
        a = np.asarray(actual, dtype=float)
        b = np.asarray(expected, dtype=float)
        assert a.shape == b.shape, (a.shape, b.shape)
        discrepancy = float(np.max(np.abs(a-b), initial=0))/max(1., float(np.max(np.abs(b), initial=0)))
        maximum = max(maximum, discrepancy)
        assert discrepancy < 2e-11, (model, context, discrepancy)
        checked['values'] += a.size

    refs = [ref for ref in manifest['cases'] if ref['case']['model'] == model]
    if limit is not None:
        refs = refs[:limit]
    for ref in refs:
        case = ref['case']
        region_path = Path(ref['regional_input'])
        assert sha(region_path) == ref['regional_input_sha256']
        assert sha(ref['result_path']) == ref['result_sha256']
        with np.load(region_path) as archive:
            region = {k:archive[k].copy() for k in archive.files}
        for index, phase in enumerate(ref['phases']):
            native = Path(phase['archive'])
            context = str(native)
            digest = sha(native)
            assert digest == phase['archive_sha256']
            assert digest == result['native_source_hashes'][str(native)]
            hashes[str(native)] = digest
            expected = independent_readouts(native, case, region)
            derived = run/'workers'/model/f'{case["slot"]:04d}_{index:02d}_{phase["phase"]}.npz'
            assert sha(derived) == result['derived_artifact_hashes'][str(derived.relative_to(run))]
            with np.load(derived) as archive:
                assert set(archive.files) == set(expected)
                for field, values in expected.items():
                    context = str(native)+' '+field
                    compare(archive[field], values)
                    checked['arrays'] += 1
            checked['phases'] += 1
            checked['frames'] += len(expected['time'])
            durations = case['pulse_durations'] if phase['phase'] == 'cue_on' else [phase['pulse_duration']]
            for duration in durations:
                use = expected['time'] <= duration if phase['phase'] == 'cue_on' else np.ones(len(expected['time']), dtype=bool)
                sample = {k:v[use] for k,v in expected.items()}
                base = (case['slot'], phase['phase'], float(duration))
                endpoint = endpoint_map[base]
                assert endpoint['native_archive'] == str(native) and endpoint['native_archive_sha256'] == digest
                for field, values in sample.items():
                    if field == 'time':
                        continue
                    context = (base, field)
                    compare(float(endpoint[field+'_first']), values[0])
                    compare(float(endpoint[field+'_last']), values[-1])
                for flag in ['scalar_orientation_success', 'scalar_full_polarity_success', 'native_field_tail_pass']:
                    assert endpoint[flag] == ('' if phase.get(flag) is None else str(phase[flag]))
                checked['endpoints'] += 1
                for roi in ['organization_core', 'worker_core']:
                    for floor, observable, losses, field in [
                        *[(rule, 'excess_fraction', [.1,.25,.5,.75,.9], roi+'_'+rule+'_excess_fraction') for rule in FLOORS],
                        ('none_raw_activity', 'raw_fraction', [.1,.25,.5], roi+'_raw_fraction')]:
                        for loss in losses:
                            event = event_map[base+(roi, floor, observable, loss)]
                            rebuilt = entry(sample['time'], sample[field], loss)
                            for name, value in rebuilt.items():
                                context = (base, roi, floor, loss, name)
                                if isinstance(value, (int,float)):
                                    compare(float(event[name]), value)
                                else:
                                    assert event[name] == ('' if value is None else value), context
                            checked['events'] += 1
    if limit is None:
        assert checked['endpoints'] == len(endpoints)
        assert checked['events'] == len(events)
    report = dict(passed=True, generated_at=datetime.now(timezone.utc).isoformat(),
                  dataset=label, model=model, manifest_sha256=sha(run/'manifest.json'),
                  producer_result_sha256=sha(run/'workers'/model/'result.json'),
                  reviewer_source_sha256=sha(Path(__file__)), cases=len(refs),
                  checked=checked, maximum_scaled_error=maximum,
                  native_archive_sha256=hashes, wall_seconds=time.monotonic()-started,
                  new_native_integrations=0, limit=limit)
    path = OUT/'results'/f'{task:02d}.json'
    dump(path, report)
    print(json.dumps({k:v for k,v in report.items() if k != 'native_archive_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', type=int, required=True)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    audit(args.task, args.limit)
