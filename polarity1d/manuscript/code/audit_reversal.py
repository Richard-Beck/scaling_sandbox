"""Independent audit of corrected matched-reversal inputs and native trial arrays.

Does not import the campaign runner, gate functions, or model runtime. Run final
audits on SLURM alongside collection; --partial permits a read-only progress check.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
MODELS = {'wave_pinning', 'goryachev', 'debelly', 'holmes_model3', 'spring'}
CONDITIONS = {'fixed_10pct', 'fixed_20pct', 'relative_20pct'}
GATES = ('orientation', 'full_polarity')


def read(path):
    return json.loads(path.read_text())


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, message):
    require(np.isfinite(actual) and np.isfinite(expected)
            and np.isclose(actual, expected, rtol=1e-10, atol=1e-12), message)


def checkpoint_fields(z, case):
    """Duplicate only the documented field-coordinate convention, not dynamics."""
    fields = {name[6:]:z[name][-1] for name in z.files if name.startswith('field_')}
    if case['model'] == 'debelly':
        names = ['rac', 'rho', 'mca_bound', 'mca_unbound', 'displacement',
                 'velocity', 'tension', 'boundary_separation']
        fields = {name:fields[name] for name in names}
        fields['displacement'] = fields['displacement'] - np.mean(fields['displacement'])
        fields['physical_length'] = case['length'] + fields.pop('boundary_separation')
    return fields


def audit_inputs(run, manifest, package):
    provenance = read(run / 'provenance.json')
    require(sha(run / 'manifest.json') == provenance['manifest_sha256'], 'Manifest hash changed')
    for name, digest in provenance['sources'].items():
        require(sha(package / name) == digest, f'Source/input hash changed: {name}')
        require(sha(run / 'frozen_sources' / name) == digest, f'Frozen source hash changed: {name}')
    require(set(manifest['models']) == MODELS, 'Unexpected active models')
    require(set(manifest['excluded_models']) == {'legi', 'otsuji'}, 'Wrong excluded models')
    require(manifest['amplitude_factor'] == 1., 'Wrong amplitude factor')
    require(manifest['recovery_horizon'] == 10000. and manifest['gate_window'] == 500.
            and manifest['sample_interval'] == 2.5, 'Wrong observation clocks')
    cases = manifest['cases']
    require(len(cases) == 330, 'Expected 330 preparation conditions')
    require([c['case_index'] for c in cases] == list(range(len(cases))), 'Case indices inconsistent')
    reference = package / 'data/establishment_basal'
    original = read(reference / 'manifest.json')
    originals = {c['case_index']:c for c in original['cases']}
    summaries = {int(c['case_index']):c for c in csv.DictReader(
        (reference / 'results/endpoints.csv').open())}
    expected = {i for i,c in originals.items()
                if c['model'] in MODELS and c['amplitude_factor'] == 1.}
    require({c['establishment_case_index'] for c in cases} == expected, 'Preparation membership changed')
    grids = read(package / 'data/inputs/reversal_duration_grids.json')['durations_by_model']
    initial = {}
    for c in cases:
        old = originals[c['establishment_case_index']]
        summary = summaries[c['establishment_case_index']]
        for key in ('model','length','cells','condition','width','width_mode',
                    'amplitude','amplitude_factor','width_fraction_at_reference','fixed_scales'):
            require(c[key] == old[key], f'Preparation input changed for case {c["case_index"]}: {key}')
        require(c['condition'] in CONDITIONS, 'Unknown cue condition')
        require(c['durations'] == grids[c['model']], 'Duration grid differs from frozen input')
        require(c['durations'][0] == 0. and all(a < b for a,b in zip(
            c['durations'], c['durations'][1:])), 'Grid must be unique, sorted, nonnegative')
        require(all(np.isfinite(s) and s > 0 for s in c['fixed_scales'].values()), 'Invalid field scale')
        if summary['status'] != 'polarized':
            require(not c['eligible'] and c['checkpoint'] is None, 'Unqualified endpoint used')
            continue
        checkpoint = package / c['checkpoint']
        require(sha(checkpoint) == c['checkpoint_sha256'] == summary['selected_checkpoint_sha256'],
                'Selected checkpoint digest changed')
        require(checkpoint.resolve() == (reference / summary['selected_checkpoint']).resolve(),
                'Selected checkpoint replaced')
        with np.load(checkpoint, allow_pickle=False) as z:
            ori, contrast = float(z['orientation'][-1]), float(z['contrast'][-1])
            close(c['initial_orientation'], ori, 'Initial orientation differs from checkpoint')
            close(c['initial_contrast'], contrast, 'Initial contrast differs from checkpoint')
            require(c['eligible'] == bool(np.isfinite(ori) and ori < -.05 and contrast >= .05),
                    'Incorrect eligibility gate')
            initial[c['case_index']] = checkpoint_fields(z, c)
    return cases, initial


def audit_trial(row, arrays, case, initial):
    duration = row['duration']
    require(np.isfinite(duration) and duration >= 0., 'Invalid pulse duration')
    t = arrays['recovery_time']
    o = arrays['recovery_orientation']
    contrast = arrays['recovery_contrast']
    require(len(t) == len(o) == len(contrast) == 4001, 'Recovery sample coverage incomplete')
    require(np.array_equal(t, np.arange(4001) * 2.5), 'Recovery sampling times changed')
    require(np.isfinite(o).all() and np.isfinite(contrast).all(), 'Nonfinite scalar histories')
    tail = t >= 9500.
    require(int(tail.sum()) == 201, 'Incomplete final-500 sample set')
    y, q = o[tail], contrast[tail]
    orientation = bool(np.min(y) > .05 and np.min(q) >= .05 and np.ptp(y) < .005)
    error = float(np.max(np.abs(y / abs(case['initial_orientation']) - 1.)))
    full = bool(orientation and error <= .01)
    require(row['orientation_success'] == orientation, 'Orientation gate differs from arrays')
    require(row['full_polarity_success'] == full, 'Scalar-strength gate differs from arrays')
    close(row['maximum_relative_polarity_error_last_500'], error, 'Strength error differs')
    close(row['minimum_contrast_last_500'], float(np.min(q)), 'Tail contrast differs')
    close(row['orientation_range_last_500'], float(np.ptp(y)), 'Orientation range differs')
    close(row['final_orientation'], float(y[-1]), 'Final orientation differs')
    pulse_t = arrays['pulse_time']
    require(len(pulse_t) >= 1 and pulse_t[0] == 0. and np.isclose(pulse_t[-1], duration),
            'Pulse timing differs')
    close(float(arrays['pulse_orientation'][0]), case['initial_orientation'],
          'Pulse does not restart original orientation')
    close(float(arrays['pulse_contrast'][0]), case['initial_contrast'],
          'Pulse does not restart original contrast')
    maximum, excursion = 0., 0.
    fields = {name[len('fieldtail_'):-len('_min')] for name in arrays
              if name.startswith('fieldtail_') and name.endswith('_min')}
    require(fields == set(initial), 'Monitored native-field coverage changed')
    for name in fields:
        scale = case['fixed_scales'][name]
        lower, upper, previous, last = [arrays[f'fieldtail_{name}_{key}']
                                      for key in ('min','max','penultimate','last')]
        require(lower.shape == upper.shape == previous.shape == last.shape == initial[name].shape,
                'Field-tail shape changed')
        require(all(np.isfinite(v).all() for v in (lower,upper,previous,last)), 'Nonfinite field statistic')
        require(np.all(lower <= upper) and np.all(previous >= lower) and np.all(previous <= upper)
                and np.all(last >= lower) and np.all(last <= upper), 'Inconsistent tail extrema')
        variation = float(np.max(upper - lower)) / scale
        drift = float(np.max(np.abs(last - previous))) / 2.5 * 500. / scale
        maximum = max(maximum, variation, drift)
        close(row['field_tail_metrics'][name]['range_scaled'], variation, 'Field range differs')
        close(row['field_tail_metrics'][name]['projected_drift_500'], drift, 'Field drift differs')
        if duration == 0.:
            a, b, start = [arrays[f'no_cue_{name}_{key}'] for key in ('min','max','initial')]
            require(np.array_equal(start, initial[name]), 'No-cue control initial field differs from checkpoint')
            require(a.shape == b.shape == start.shape and np.isfinite(a).all() and np.isfinite(b).all()
                    and np.all(a <= b), 'Invalid no-cue whole-history extrema')
            excursion = max(excursion, float(max(np.max(np.abs(a-start)), np.max(np.abs(b-start)))) / scale)
    close(row['field_tail_max_scaled'], maximum, 'Field flag maximum differs')
    require(row['field_tail_pass'] == bool(maximum <= .001), 'Field-tail flag differs')
    final = [v for name,v in arrays.items() if name == 'finalstate' or name.startswith('finalstate_')]
    require(bool(final) and all(np.isfinite(v).all() for v in final), 'Missing/nonfinite final native state')
    if duration == 0.:
        close(row['no_cue_maximum_scaled_excursion'], excursion, 'No-cue full-history excursion differs')
        left = bool(np.max(y) < -.05 and np.min(q) >= .05 and np.ptp(y) < .005)
        valid = bool(excursion <= .01 and not orientation and left)
        require(row['no_cue_left_tail_valid'] == left, 'No-cue left-state gate differs')
        require(row['no_cue_preparation_valid'] == valid, 'No-cue preparation validity differs')


def audit(run, package, partial=False, inputs_only=False, validate_collected=True):
    manifest = read(run / 'manifest.json')
    cases, initial = audit_inputs(run, manifest, package)
    report = dict(scope='inputs' if inputs_only else 'all retained trials', cases=len(cases),
                  eligible_cases=sum(c['eligible'] for c in cases), trials=0,
                  pending_chunks=0, execution_errors=0, status_counts={}, errors=[], selected_readouts=[])
    if inputs_only:
        return dict(report, passed=True, complete=False)
    groups = {c['case_index']:[] for c in cases}
    for path in sorted((run / 'work/jobs').glob('r*.json')):
        if not re.fullmatch(r'r\d+\.json', path.name):
            continue
        require(sha(path) == read(path.with_suffix('.sha256.json'))['sha256'], 'Frozen job hash changed')
        round_id = int(path.stem[1:])
        for job in read(path):
            case = cases[job['case_index']]
            require(case['eligible'], 'Scheduled ineligible case')
            chunk = run / 'work/chunks' / f'r{round_id}_{job["slot"]:05d}.json'
            if chunk.with_suffix('.error.json').exists():
                report['execution_errors'] += 1
                continue
            if not chunk.exists():
                report['pending_chunks'] += 1
                continue
            data = read(chunk)
            require(data['case_index'] == job['case_index'] and data['slot'] == job['slot']
                    and data['round'] == round_id, 'Chunk does not match scheduled job')
            require(data['checkpoint_sha256'] == case['checkpoint_sha256'], 'Trial checkpoint changed')
            require([t['duration'] for t in data['trials']] == job['durations'], 'Trial coverage differs')
            trace = chunk.with_suffix('.npz')
            require(sha(trace) == data['trace_sha256'], 'Native trial archive hash changed')
            with np.load(trace, allow_pickle=False) as z:
                for row in data['trials']:
                    require(run / row['trace'] == trace, 'Wrong trial trace path')
                    prefix = row['trace_prefix']
                    arrays = {key[len(prefix):]:z[key] for key in z.files if key.startswith(prefix)}
                    audit_trial(row, arrays, case, initial[case['case_index']])
                    groups[case['case_index']].append(row)
                    report['trials'] += 1
    complete = not report['pending_chunks'] and not report['execution_errors']
    counts = Counter()
    published = {}
    result_file = run / 'results/records.json.gz'
    if result_file.exists():
        import gzip
        with gzip.open(result_file, 'rt') as stream:
            published = {c['case_index']:c for c in json.load(stream)['cases']}
    for case in cases:
        ts = sorted(groups[case['case_index']], key=lambda t:t['duration'])
        require(len(ts) == len({t['duration'] for t in ts}), 'Duplicate durations')
        zero = next((t for t in ts if t['duration'] == 0.), None)
        valid = bool(case['eligible'] and zero and zero['no_cue_preparation_valid'])
        if not partial:
            if case['eligible']:
                require(zero is not None, 'No-cue control missing')
            if valid:
                require(set(case['durations']) <= {t['duration'] for t in ts}, 'Base duration coverage incomplete')
        for gate in GATES:
            passes = [t for t in ts if t['duration'] > 0. and t[gate+'_success']]
            selected = passes[0] if valid and passes else None
            status = 'ineligible' if not case['eligible'] else ('preparation_invalid' if zero and not valid
                      else 'success' if selected else 'censored' if complete and not partial else 'pending')
            counts[status] += 1
            for a,b in zip(ts, ts[1:]):
                if valid and not a[gate+'_success'] and b[gate+'_success']:
                    resolved = (b['duration'] - a['duration']) / b['duration'] <= .01000000001
                    if not resolved and not partial:
                        raise ValueError('Unresolved sampled fail-to-pass boundary')
            report['selected_readouts'].append(dict(case_index=case['case_index'], definition=gate,
                status=status, endpoint=selected['duration'] if selected else None,
                selected_tail_field_pass_0_1_percent=selected['field_tail_pass'] if selected else None))
            if not partial and validate_collected:
                require(case['case_index'] in published, 'Collected case missing')
                q = published[case['case_index']]['result'][gate]
                require(q['status'] == status, 'Collected status differs from trial evidence')
                require(q['endpoint'] == (selected['duration'] if selected else None),
                        'Collected endpoint is not shortest tested passing duration')
                if selected:
                    require(q['selected_tail_field_pass_0_1_percent'] == selected['field_tail_pass'],
                            'Collected selected field flag differs')
    report['status_counts'] = dict(counts)
    return dict(report, passed=complete or partial, complete=complete and not partial)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', '--campaign', type=Path, default=PACKAGE / 'data/runs/corrected_reversal_20261003')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--partial', action='store_true')
    parser.add_argument('--inputs-only', action='store_true')
    parser.add_argument('--before-publication', action='store_true',
                        help='Audit complete trial evidence and export recomputed endpoints before final tables exist.')
    args = parser.parse_args()
    try:
        result = audit(args.run.resolve(), PACKAGE, args.partial, args.inputs_only, not args.before_publication)
    except Exception as error:
        result = dict(passed=False, complete=False, errors=[str(error)])
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'selected_readouts'}, indent=2))
    if not result['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
