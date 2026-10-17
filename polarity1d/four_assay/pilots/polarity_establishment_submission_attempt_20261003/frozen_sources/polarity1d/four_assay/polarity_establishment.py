"""Independent native-rest -> physical patch -> verified release duration campaign."""
from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import traceback
from collections import Counter

import numpy as np

import runner as r
import settling

CODE = Path(__file__).resolve().parent
DEFAULT = CODE / 'pilots/polarity_establishment_20261003'
FACTORS = [.1, .25, .5, 1., 2., 4., 10.]
REFERENCES = {'wave_pinning': 21., 'otsuji': 21., 'goryachev': 21.,
              'legi': 17., 'debelly': 27., 'holmes_model3': 30., 'spring': 1.6}


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def read(path):
    return json.loads(path.read_text())


class Patch:
    """Fixed native physical width/peak; right is the spatial mirror of left."""
    def __init__(self, amplitude, width, side='left', stop=None):
        self.amplitude, self.width, self.side, self.stop = amplitude, width, side, stop
        self.event_times = () if stop is None else (float(stop),)

    def __call__(self, t, x, length):
        x = np.asarray(x, float)
        if self.stop is not None and t >= self.stop:
            return np.zeros_like(x)
        distance = x if self.side == 'left' else length - x
        q = distance / self.width
        return np.where((q >= 0) & (q <= 1),
                        self.amplitude * .5 * (1 + np.cos(np.pi * np.clip(q, 0, 1))), 0.)


def sources():
    hashes = r.sources()
    for path in [Path(__file__), CODE / 'polarity_establishment_plot.py']:
        hashes[str(path.relative_to(r.ROOT))] = r.sha(path)
    return hashes


def verify(out):
    provenance = read(out / 'provenance.json')
    if sources() != provenance['sources']:
        raise RuntimeError('Campaign sources changed after freezing')
    if r.sha(out / 'manifest.json') != provenance['manifest_sha256']:
        raise RuntimeError('Campaign manifest changed')
    return read(out / 'manifest.json')


def make_manifest():
    scale_rows = read(CODE / r.PROTOCOL['steady_state']['fixed_scales'])
    cases = []
    for task in r.tasks():
        if task['assay'] != 'finite_countercue':
            continue
        model, length = task['model'], task['length']
        scales = {}
        for row in scale_rows:
            if row['model'] == model and row['length'] == length:
                for field, value in row['fixed_scales'].items():
                    scales[field] = max(scales.get(field, 1e-8), value)
        lo, hi = (.002, 50.) if model == 'spring' else ((.5, 2000.) if model == 'debelly' else (.05, 5120.))
        durations = [0.] + [float(v) for v in np.geomspace(lo, hi, 17)]
        for fraction in [.1, .2]:
            for factor in FACTORS:
                cases.append(dict(case_index=len(cases), model=model, length=length,
                                  cells=task['cells'], width=REFERENCES[model] * fraction,
                                  width_fraction_at_reference=fraction, reference_length=REFERENCES[model],
                                  amplitude_factor=factor,
                                  amplitude=r.PROTOCOL['normalization']['amplitude_by_model'][model] * factor,
                                  fixed_scales=scales, base_durations=durations))
    return dict(created_at=r.now(), cases=cases, amplitude_factors=FACTORS,
                steady_state=r.PROTOCOL['steady_state'],
                protocol='Native resting input -> fixed physical left cosine patch -> zero cue. Independent native start for every duration; no pre-settling.',
                mirror_rule='Identical A and width for future opposite-end countercue; right(x)=left(L-x).',
                normalization='Fixed peak A and physical width across lengths within each model. Native units; Spring units are not microns.',
                short_domain_policy='Restrict the unchanged patch to the actual domain; do not rescale width or peak.',
                scale_rule='Per-field maximum of frozen v11 conditioning/release scales at this model/length, constant across A, width and duration.',
                selection=dict(relative_strength_tolerance=.01, minimum_strength=.05,
                               minimum_contrast=.05, absolute_polarity_arrival_floor=1e-6,
                               objective='Maximum cue-aligned final polarity (-orientation), among verified endpoints with contrast >=0.05. This avoids optimizing normalized polarity of a vanishing signal.',
                               duration='Shortest tested duration (including zero if no cue is needed) within 1% of strongest verified polarity and with contrast >=0.05.',
                               total_time='Selected duration plus first persistent post-release 1% polarity arrival; excludes convergence confirmation/holdout time.'),
                refinement=dict(strength_rounds=2, top_durations=3, upper_extensions=2,
                                extension_factor=4., duration_cap=60000., threshold_steps=8,
                                relative_duration_resolution=.01,
                                monotonicity='No global monotonicity assumption; refine every sampled transition into the best-strength band.'),
                endpoint_diagnostic=dict(tolerance=.01, method='Greedy complete-link clusters of full monitored native profiles on frozen scales; no endpoint-count-specific refinement.',
                                         interpretation='Apparent number among verified sampled endpoints, including zero-duration control. Continuous families and numerical tolerance can change the count.'),
                outputs='One final collection: cue duration, total establishment time, selected final polarity, apparent endpoint count.')


def persistent_arrival(time, values, target, tolerance, dwell=1000.):
    if not np.isfinite(target):
        return None
    outside = np.flatnonzero(~np.isfinite(values) | (np.abs(values - target) > tolerance))
    index = outside[-1] + 1 if len(outside) else 0
    if index >= len(time) or time[-1] - time[index] < dwell:
        return None
    return float(time[index])


def release(case, state, manifest):
    rule = manifest['steady_state']
    key, length, cells = case['model'], case['length'], case['cells']
    scales, trace, elapsed, consecutive, checks = case['fixed_scales'], None, 0., 0, []
    verified = False
    # Resolve rapid post-withdrawal changes without sampling every entire
    # 60000-unit continuation at this fine interval (especially Spring).
    early = r.run(key, length, state, r.Cue('zero', key), 10.,
                  interval=.05 if key == 'spring' else .25, n=cells)
    if not all(np.all(np.isfinite(v)) for v in settling.monitored(early, key, length).values()):
        raise RuntimeError('Nonfinite native field in early release')
    state, trace, elapsed = early.state, early, 10.
    while elapsed + rule['chunk'] <= rule['cap']:
        tr = r.run(key, length, state, r.Cue('zero', key), rule['chunk'],
                   interval=rule['sample_interval'], n=cells)
        state = tr.state
        fields = settling.monitored(tr, key, length)
        if not all(np.all(np.isfinite(v)) for v in fields.values()):
            raise RuntimeError('Nonfinite native field in release')
        metric = settling.diagnostics(tr, key, length, scales, rule)
        trace = r.concatenate(trace, tr, elapsed)
        elapsed += rule['chunk']
        checks.append(dict(time=elapsed, kind='convergence', **metric))
        consecutive = consecutive + 1 if metric['steady'] else 0
        if consecutive < rule['consecutive_windows']:
            continue
        if elapsed + rule['holdout'] > rule['cap']:
            break
        target = {k: v[-1].copy() for k, v in settling.monitored(trace, key, length).items()}
        hold = r.run(key, length, state, r.Cue('zero', key), rule['holdout'],
                     interval=rule['sample_interval'], n=cells)
        if not all(np.all(np.isfinite(v)) for v in settling.monitored(hold, key, length).values()):
            raise RuntimeError('Nonfinite native field in holdout')
        state = hold.state
        excursion = float(settling.deviations(settling.monitored(hold, key, length), target, scales).max())
        final = settling.diagnostics(hold, key, length, scales, rule)
        trace = r.concatenate(trace, hold, elapsed)
        elapsed += rule['holdout']
        verified = excursion <= rule['holdout_tolerance'] and final['steady']
        checks.append(dict(time=elapsed, kind='holdout', maximum_scaled_excursion=excursion,
                           verified=verified, **final))
        if verified:
            break
        consecutive = 0
    return trace, dict(verified=bool(verified), elapsed=elapsed, checks=checks)


def token(duration):
    return format(float(duration), '.12g').replace('.', 'p').replace('+', '')


def trial(out, case, duration, manifest):
    folder = out / 'work/cases' / f"case_{case['case_index']:04d}"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / ('duration_' + token(duration))
    if path.with_suffix('.json').exists():
        return read(path.with_suffix('.json'))
    key, length, cells = case['model'], case['length'], case['cells']
    state = r.initial(key, length, cells)
    pulse = None
    if duration > 0:
        # Apply the pulse through its left limit, then restart under zero input.
        pulse = r.run(key, length, state, Patch(case['amplitude'], case['width']),
                      duration, interval=min(5., duration / 20), n=cells)
        state = pulse.state
    tr, convergence = release(case, state, manifest)
    fields = settling.monitored(tr, key, length)
    target = {k: v[-1].copy() for k, v in fields.items()}
    distances = settling.deviations(fields, target, case['fixed_scales'])
    final = float(tr.orientation[-1])
    valid_polarity = bool(np.isfinite(final))
    tolerance = max(abs(final) * .01, manifest['selection']['absolute_polarity_arrival_floor'])
    arrival = persistent_arrival(tr.time, tr.orientation, final, tolerance) if convergence['verified'] else None
    field_arrival = persistent_arrival(tr.time, distances, 0., .01) if convergence['verified'] else None
    strength = max(0., -final) if valid_polarity else 0.
    checkpoint = path.with_suffix('.checkpoint.npz')
    # A compact full-native restart record, independent of profile clustering.
    last = type(tr)(**{**vars(tr), 'time': tr.time[-1:]})
    for name in ['active', 'orientation', 'contrast', 'left', 'right']:
        setattr(last, name, getattr(tr, name)[-1:])
    last.fields = {k: v[-1:] for k, v in tr.fields.items()}
    r.save_checkpoint(checkpoint, key, tr.state, last)
    sample = np.unique(np.linspace(0, len(tr.time) - 1, min(61, len(tr.time))).astype(int))
    arrays = dict(release_time=tr.time, release_orientation=tr.orientation,
                  release_contrast=tr.contrast, maximum_scaled_field_distance=distances,
                  profile_time=tr.time[sample], x=tr.x,
                  **{'profile_' + k: v[sample] for k, v in fields.items()},
                  **{'endpoint_' + k: v for k, v in target.items()})
    if pulse is not None:
        pidx = np.unique(np.linspace(0, len(pulse.time) - 1, min(41, len(pulse.time))).astype(int))
        arrays.update(pulse_time=pulse.time[pidx], pulse_orientation=pulse.orientation[pidx],
                      pulse_contrast=pulse.contrast[pidx], pulse_active=pulse.active[pidx])
    np.savez_compressed(path.with_suffix('.npz'), **arrays)
    record = dict(duration=float(duration), verified=convergence['verified'],
                  final_orientation=final if valid_polarity else None,
                  cue_aligned_strength=strength, final_contrast=float(tr.contrast[-1]),
                  polarity_arrival=arrival, field_arrival=field_arrival,
                  total_establishment_time=float(duration + arrival) if arrival is not None else None,
                  polarity_arrival_tolerance=tolerance if valid_polarity else None,
                  convergence=convergence, checkpoint=str(checkpoint.relative_to(out)),
                  checkpoint_sha256=r.sha(checkpoint), trace=str(path.with_suffix('.npz').relative_to(out)))
    write(path.with_suffix('.json'), record)
    print(case['case_index'], duration, record['verified'], strength, flush=True)
    return record


def qualifying(trials, manifest):
    valid = [t for t in trials if t['verified'] and t['final_orientation'] is not None]
    if not valid:
        return [], None
    rule = manifest['selection']
    observable = [t for t in valid if t['final_contrast'] >= rule['minimum_contrast']]
    best = max(t['cue_aligned_strength'] for t in (observable or valid))
    passed = [t for t in valid if best >= rule['minimum_strength']
              and t['cue_aligned_strength'] >= best * (1 - rule['relative_strength_tolerance'])
              and t['final_contrast'] >= rule['minimum_contrast'] and t['polarity_arrival'] is not None]
    return sorted(passed, key=lambda t: t['duration']), best


def strength_candidates(trials, manifest):
    ordered = sorted(trials, key=lambda t: t['duration'])
    ranked = sorted([t for t in ordered if t['verified'] and t['duration'] > 0
                     and t['final_contrast'] >= manifest['selection']['minimum_contrast']],
                    key=lambda t: (-t['cue_aligned_strength'], t['duration']))
    if not ranked or ranked[0]['cue_aligned_strength'] < manifest['selection']['minimum_strength']:
        return []
    candidates = set()
    for t in ranked[:manifest['refinement']['top_durations']]:
        i = ordered.index(t)
        for j in [i - 1, i + 1]:
            if 0 <= j < len(ordered):
                a, b = sorted([t['duration'], ordered[j]['duration']])
                if a > 0 and b / a > 1.01:
                    candidates.add(float(np.sqrt(a * b)))
    return sorted(candidates)


def threshold_candidates(trials, manifest):
    passed, _ = qualifying(trials, manifest)
    if not passed:
        return []
    pass_d = {t['duration'] for t in passed}
    ordered = sorted(trials, key=lambda t: t['duration'])
    candidates = []
    for a, b in zip(ordered, ordered[1:]):
        if b['duration'] in pass_d and a['duration'] not in pass_d:
            if (b['duration'] - a['duration']) / b['duration'] > .01:
                candidates.append(float(np.sqrt(a['duration'] * b['duration'])) if a['duration'] else b['duration'] / 2)
    return candidates


def endpoint_count(out, case, trials, tolerance):
    profiles, clusters, assignments = [], [], []
    for t in sorted([t for t in trials if t['verified']], key=lambda t: t['duration']):
        with np.load(out / t['trace']) as z:
            vector = np.concatenate([np.ravel(z['endpoint_' + name]) / scale
                                     for name, scale in sorted(case['fixed_scales'].items())])
        profiles.append(vector)
        index = len(profiles) - 1
        match = next((i for i, members in enumerate(clusters)
                      if all(np.max(abs(vector - profiles[j])) <= tolerance for j in members)), None)
        if match is None:
            match = len(clusters)
            clusters.append([])
        clusters[match].append(index)
        assignments.append(dict(duration=t['duration'], endpoint_cluster=match))
    return len(clusters), assignments


def worker(out, slot):
    manifest = verify(out)
    case = manifest['cases'][slot]
    trials = []
    def evaluate(duration):
        t = trial(out, case, duration, manifest)
        if not any(v['duration'] == duration for v in trials):
            trials.append(t)
        return t
    for duration in case['base_durations']:
        evaluate(duration)
    # Bounded upper exploration only when strongest polarity is still increasing.
    for _ in range(manifest['refinement']['upper_extensions']):
        ordered = sorted([t for t in trials if t['verified']
                          and t['final_contrast'] >= manifest['selection']['minimum_contrast']], key=lambda t: t['duration'])
        if len(ordered) < 2:
            break
        last = ordered[-1]
        previous_best = max(t['cue_aligned_strength'] for t in ordered[:-1])
        if last['cue_aligned_strength'] < .05 or last['cue_aligned_strength'] <= previous_best * 1.01:
            break
        duration = min(last['duration'] * 4, manifest['refinement']['duration_cap'])
        if duration <= max(t['duration'] for t in trials):
            break
        evaluate(duration)
    for _ in range(manifest['refinement']['strength_rounds']):
        for duration in strength_candidates(trials, manifest):
            evaluate(duration)
    for _ in range(manifest['refinement']['threshold_steps']):
        candidates = threshold_candidates(trials, manifest)
        if not candidates:
            break
        for duration in candidates:
            evaluate(duration)
    verify(out)
    trials.sort(key=lambda t: t['duration'])
    passed, best = qualifying(trials, manifest)
    selected = passed[0] if passed else None
    count, assignments = endpoint_count(out, case, trials, manifest['endpoint_diagnostic']['tolerance'])
    upper = trials[-1]
    earlier = [t['cue_aligned_strength'] for t in trials[:-1] if t['verified']
               and t['final_contrast'] >= manifest['selection']['minimum_contrast']]
    upper_unresolved = bool(not upper['verified'] or (earlier and upper['cue_aligned_strength'] >= .05
                            and upper['cue_aligned_strength'] > max(earlier) * 1.01))
    lower = max((t['duration'] for t in trials if selected and t['duration'] < selected['duration']), default=0.)
    bracket = ((selected['duration'] - lower) / selected['duration'] if selected['duration'] else 0.) if selected else None
    failures = [t for t in trials if not t['verified']]
    status = ('polarized' if selected else ('unsettled' if best is None else
              'unpolarized' if best < .05 else 'no_qualifying_endpoint'))
    row = dict(**{k: case[k] for k in ['case_index', 'model', 'length', 'cells', 'width',
                                      'width_fraction_at_reference', 'amplitude_factor', 'amplitude']},
               status=status, best_observed_strength=best,
               cue_duration=selected['duration'] if selected else None,
               total_establishment_time=selected['total_establishment_time'] if selected else None,
               selected_final_polarity=selected['cue_aligned_strength'] if selected else best,
               selected_final_contrast=selected['final_contrast'] if selected else None,
               post_release_arrival=selected['polarity_arrival'] if selected else None,
               apparent_endpoint_count=count if count else None,
               tested_durations=len(trials), unsettled_trials=len(failures),
               upper_optimum_unresolved=upper_unresolved,
               selected_relative_duration_bracket=bracket,
               duration_boundary_unresolved=bool(selected and bracket > .010000001),
               selected_checkpoint=selected['checkpoint'] if selected else None,
               selected_checkpoint_sha256=selected['checkpoint_sha256'] if selected else None)
    write(out / 'work/cases' / f'case_{slot:04d}' / 'result.json',
          dict(generated_at=r.now(), summary=row, trials=trials, endpoint_assignments=assignments))


def submit(out, action, array=None, dependency=None):
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '
    command += shlex.join([sys.executable, str(Path(__file__).resolve()), action, '--out', str(out)])
    if array is not None:
        command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    opts = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
            '--mem=3G', '--time=12:00:00' if action in ['worker', 'monitor'] else '--time=00:30:00',
            '--job-name=polarity-establish-' + action, '--chdir=' + str(r.ROOT),
            '--output=' + str(out / 'work/logs' / (action + '_%A_%a.out')),
            '--error=' + str(out / 'work/logs' / (action + '_%A_%a.err'))]
    if array is not None:
        opts.append('--array=' + array)
    if dependency:
        opts.append('--dependency=afterany:' + dependency)
    return subprocess.check_output(opts + ['--wrap', command], text=True).strip().split(';')[0]


def initialize(out):
    out.mkdir(parents=True, exist_ok=False)
    for folder in ['work/cases', 'work/logs', 'results', 'frozen_sources']:
        (out / folder).mkdir(parents=True, exist_ok=True)
    manifest = make_manifest()
    write(out / 'manifest.json', manifest)
    hashes = sources()
    for rel in hashes:
        src, dest = r.ROOT / rel, out / 'frozen_sources' / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
    write(out / 'provenance.json', dict(generated_at=r.now(), sources=hashes,
                                      manifest_sha256=r.sha(out / 'manifest.json')))
    collect(out, final=False)


def launch(out, resume=False):
    if not resume:
        initialize(out)
    manifest = verify(out)
    missing = [c['case_index'] for c in manifest['cases']
               if not (out / 'work/cases' / f"case_{c['case_index']:04d}" / 'result.json').exists()]
    if not missing:
        collect(out, final=True)
        return
    # Split arrays to stay below site MaxArraySize; no concurrency throttle.
    limit = int(subprocess.check_output(['scontrol', 'show', 'config'], text=True).split('MaxArraySize')[1].split('=')[1].split()[0])
    jobs = []
    for start in range(0, len(missing), min(900, limit - 1)):
        if max(missing) >= limit:
            raise RuntimeError('Case indices exceed site MaxArraySize; use array mapping before submission')
        jobs.append(submit(out, 'worker', array=','.join(map(str, missing[start:start + min(900, limit - 1)]))))
    final = submit(out, 'finalize', dependency=':'.join(jobs))
    monitor = submit(out, 'monitor')
    write(out / 'launch.json', dict(generated_at=r.now(), arrays=jobs, finalize_job=final,
                                    monitor_job=monitor, submitted_cases=len(missing), stage='running'))
    print(json.dumps(read(out / 'launch.json')), flush=True)


def collect(out, final=False):
    import fcntl
    # Serialize monitor/finalizer writes so a progress snapshot cannot overwrite
    # final artifacts or mark a finished producer as running again.
    with (out / 'work/collection.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        status_path = out / 'results/status.json'
        if not final and status_path.exists():
            existing = read(status_path)
            if existing['producer'] == 'finished':
                return existing
        return _collect_locked(out, final)


def _collect_locked(out, final=False):
    manifest = verify(out)
    completed, errors = [], []
    for c in manifest['cases']:
        folder = out / 'work/cases' / f"case_{c['case_index']:04d}"
        if (folder / 'result.json').exists():
            completed.append(read(folder / 'result.json'))
        elif (folder / 'error.json').exists():
            errors.append(read(folder / 'error.json'))
    status = dict(generated_at=r.now(), complete=len(completed) == len(manifest['cases']),
                  producer='finished' if final else 'running', cases=len(manifest['cases']),
                  completed_cases=len(completed), execution_errors=len(errors),
                  status_counts=dict(Counter(c['summary']['status'] for c in completed)),
                  heatmaps='Generated once, at complete finalization.')
    if not final and (out / 'launch.json').exists():
        jobs = read(out / 'launch.json')['arrays']
        query = subprocess.run(['squeue', '-h', '-j', ','.join(jobs), '-o', '%i'],
                               text=True, capture_output=True)
        if query.returncode == 0 and not query.stdout.strip():
            status['producer'] = 'workers stopped; awaiting finalization'
    write(out / 'results/status.json', status)
    rows = [c['summary'] for c in completed]
    if rows:
        temp = out / 'results/endpoints.csv.tmp'
        r.csvwrite(temp, rows)
        temp.replace(out / 'results/endpoints.csv')
    plot_html = ''
    if final and status['complete']:
        for c in completed:
            passed, best = qualifying(c['trials'], manifest)
            selected = passed[0] if passed else None
            summary = c['summary']
            if summary['best_observed_strength'] != best or summary['cue_duration'] != (selected['duration'] if selected else None):
                raise RuntimeError('Strongest-endpoint / shortest-duration selection audit failed')
            if selected and summary['selected_final_polarity'] != selected['cue_aligned_strength']:
                raise RuntimeError('Selected polarity audit failed')
            for t in c['trials']:
                if r.sha(out / t['checkpoint']) != t['checkpoint_sha256']:
                    raise RuntimeError('Trial checkpoint hash mismatch')
                with np.load(out / t['trace']) as z:
                    expected = persistent_arrival(z['release_time'], z['release_orientation'],
                                                  t['final_orientation'] if t['final_orientation'] is not None else np.nan,
                                                  t['polarity_arrival_tolerance'] or 1e-6) if t['verified'] else None
                    if expected != t['polarity_arrival']:
                        raise RuntimeError('Polarity arrival audit failed')
        with gzip.open(out / 'results/records.json.gz', 'wt') as f:
            json.dump(dict(generated_at=r.now(), manifest=manifest, cases=completed), f, allow_nan=False)
        write(out / 'results/validation.json', dict(generated_at=r.now(), passed=True,
              complete_case_coverage=True, checkpoint_hashes_verified=True,
              persistent_polarity_arrivals_audited=True, source_hashes_verified=True,
              strongest_endpoint_and_shortest_duration_audited=True,
              cases=len(completed), trials=sum(len(c['trials']) for c in completed)))
        subprocess.run([sys.executable, str(CODE / 'polarity_establishment_plot.py'), '--out', str(out)], check=True)
        plot_html = ''.join(f'<p><a href="{name}.pdf">PDF</a> · <a href="{name}.svg">SVG</a></p><img width="1500" src="{name}.png">'
                            for name in ['cue_duration', 'total_establishment_time', 'final_polarity', 'apparent_endpoint_count'])
    elif final:
        status['heatmaps'] = 'Not generated: incomplete execution. Use resume; missing cases are not biological censoring.'
        write(out / 'results/status.json', status)
    html = f'''<!doctype html><html><head><meta charset="utf-8"><title>Polarity establishment</title>
    <style>body{{font-family:sans-serif;max-width:1550px;margin:30px auto}}img{{max-width:100%}}</style></head>
    <body><h1>Polarity establishment</h1><p>Data generated {status['generated_at']}; producer {status['producer']}.</p>
    <p>{status['completed_cases']} / {status['cases']} combinations completed; {status['execution_errors']} execution errors.</p>
    <p>{status['heatmaps']}</p><p><a href="status.json">Status</a> · <a href="endpoints.csv">Current endpoint table</a> ·
    <a href="../manifest.json">Protocol and definitions</a></p>{plot_html}</body></html>'''
    temp = out / 'results/index.html.tmp'
    temp.write_text(html)
    temp.replace(out / 'results/index.html')
    if final and (out / 'launch.json').exists():
        launch_record = read(out / 'launch.json')
        launch_record.update(stage='complete' if status['complete'] else 'execution_incomplete',
                             finished_at=r.now())
        write(out / 'launch.json', launch_record)
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['init', 'launch', 'resume', 'worker', 'status', 'finalize', 'monitor'])
    parser.add_argument('--out', type=Path, default=DEFAULT)
    parser.add_argument('--slot', type=int)
    args = parser.parse_args()
    out = args.out.resolve()
    if args.action == 'init':
        initialize(out)
    elif args.action in ['launch', 'resume']:
        launch(out, resume=args.action == 'resume')
    elif args.action == 'worker':
        try:
            worker(out, args.slot)
        except Exception:
            write(out / 'work/cases' / f'case_{args.slot:04d}' / 'error.json',
                  dict(case_index=args.slot, generated_at=r.now(), traceback=traceback.format_exc()))
            raise
    elif args.action == 'monitor':
        import time
        for _ in range(1300):
            current = read(out / 'results/status.json')
            if current['producer'] == 'finished':
                return
            launch_path = out / 'launch.json'
            if launch_path.exists():
                jobs = read(launch_path)['arrays']
                active = subprocess.check_output(['squeue', '-h', '-j', ','.join(jobs), '-o', '%i'], text=True).strip()
                if not active:
                    break
            collect(out)
            time.sleep(30)
        else:
            # Continue reporting even if queued workers outlive this monitor.
            submit(out, 'monitor')
    elif args.action == 'status':
        print(json.dumps(read(out / 'results/status.json')))
    else:
        status = collect(out, final=True)
        print(json.dumps(status))
        if not status['complete']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
