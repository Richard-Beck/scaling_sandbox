"""Reference-only continuation and retrospectively measured preparation arrival times.

Does not change assay inputs, gates, results, or response experiments.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, os, shlex, shutil, subprocess
import numpy as np
import runner as r

BASE = Path(__file__).resolve().parent
CAMPAIGN = BASE / 'reference' / 'continuation_20261001'
CHUNK = 1000.
LIMIT = 50000.
TOL = 1e-4
ARRIVAL = .01


def fields(tr, key, length):
    """Physical fields on the native material grid; remove rigid translation only."""
    out = {k: np.asarray(v) for k, v in tr.fields.items()}
    if key == 'debelly':
        out = {k: out[k] for k in ['rac', 'rho', 'mca_bound', 'mca_unbound',
                                  'displacement', 'velocity', 'tension', 'boundary_separation']}
        out['displacement'] = out['displacement'] - out['displacement'].mean(axis=1)[:, None]
        out['physical_length'] = length + out.pop('boundary_separation')
    return out


def scales(history, endpoint, key, length):
    result = {k: max(float(np.max(abs(v))), float(np.max(abs(endpoint[k]))), 1e-8)
              for k, v in history.items()}
    for k in ['displacement', 'physical_length', 'length']:
        if k in result:
            result[k] = length
    return result


def check(tr, key, length, scale):
    out = {}
    for k, y in fields(tr, key, length).items():
        variation = float(np.max(np.ptp(y, axis=0))) / scale[k]
        drift = float(np.max(abs(y[-1] - y[-2]))) / (tr.time[-1] - tr.time[-2]) * CHUNK / scale[k]
        out[k] = dict(range=variation, projected_drift=drift)
    return dict(fields=out, passed=all(max(v.values()) <= TOL for v in out.values()))


def distance(profile, reference, scale):
    return {k: float(np.max(abs(profile[k] - reference[k]))) / scale[k] for k in scale}


def verify():
    manifest = json.loads((CAMPAIGN / 'manifest.json').read_text())
    for path, sha in manifest['source_hashes'].items():
        if r.sha(Path(path)) != sha:
            raise RuntimeError('Frozen source changed: ' + path)
    return manifest


def setup():
    if CAMPAIGN.exists():
        raise RuntimeError('Campaign already exists')
    for name in ['inputs', 'profiles', 'histories', 'records', 'logs']:
        (CAMPAIGN / name).mkdir(parents=True)
    source = BASE / 'results/records.json.gz'
    records = json.load(gzip.open(source, 'rt'))
    preps = {(p['model'], p['length'], p['assay']): p for p in records['preparations']}
    tasks = []
    for task in r.tasks():
        if task['assay'] != 'gradient_reversal':
            continue
        key, length = task['model'], task['length']
        pair = [preps[key, length, a] for a in ['gradient_reversal', 'finite_countercue']]
        paths = []
        for p in pair:
            filename = f"task_{p['task_index']:04d}.npz"
            origin = BASE / 'work/checkpoints' / filename if key == 'spring' else (
                BASE.parent / 'archive/four_assay_normalized_v9_20261001/work/checkpoints' / filename)
            if r.sha(origin) != p['checkpoint_sha256']:
                raise RuntimeError('Checkpoint provenance mismatch: ' + str(origin))
            target = CAMPAIGN / 'inputs' / filename
            shutil.copy2(origin, target)
            paths.append(dict(path=str(target), sha256=r.sha(target), origin=str(origin)))
        for phase_index, event in enumerate(['conditioning', 'release']):
            p = pair[phase_index]
            if len(p['phases']) <= phase_index:
                raise RuntimeError('Missing old phase')
            tasks.append(dict(slot=len(tasks), model=key, length=length, cells=task['cells'],
                              event=event, phase=p['phases'][phase_index], checkpoint=paths[phase_index],
                              conditioning_checkpoint=paths[0]))
    hashes = r.sources()
    absolute = {str(r.ROOT / p): sha for p, sha in hashes.items()}
    absolute[str(Path(__file__).resolve())] = r.sha(Path(__file__).resolve())
    absolute[str(source)] = r.sha(source)
    r.write(CAMPAIGN / 'manifest.json', dict(created_at=r.now(), tasks=tasks, source_hashes=absolute,
        rule=dict(chunk=CHUNK, additional_cap=LIMIT, tolerance=TOL, consecutive_windows=2,
                  sample_interval=5, arrival_tolerance=ARRIVAL),
        purpose='Reference analysis only; original assay checkpoints remain unchanged'))
    for p in absolute:
        target = CAMPAIGN / 'sources' / Path(p).relative_to(r.ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
    python = json.loads((BASE / 'environment.json').read_text())['python']
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; ' + shlex.join(
        [python, str(Path(__file__).resolve()), 'worker']) + ' --slot "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(['sbatch', '--parsable', '--qos=small', '--ntasks=1',
        '--cpus-per-task=1', '--mem=3G', '--time=01:00:00', '--array=0-307',
        '--job-name=polarity-ss', f'--chdir={r.ROOT}',
        f'--output={CAMPAIGN}/logs/%A_%a.out', f'--error={CAMPAIGN}/logs/%A_%a.err',
        '--wrap', command], text=True).strip()
    r.write(CAMPAIGN / 'launch.json', dict(job=job, submitted_at=r.now()))
    print(job)


def worker(slot):
    manifest = verify()
    task = manifest['tasks'][slot]
    key, length, n = task['model'], task['length'], task['cells']
    for p in [task['checkpoint'], task['conditioning_checkpoint']]:
        if r.sha(Path(p['path'])) != p['sha256']:
            raise RuntimeError('Frozen input changed')
    old_state, old_trace = r.read_checkpoint(task['checkpoint']['path'], key)
    old_time = float(task['phase']['elapsed'])
    cue = r.Cue('left' if task['event'] == 'conditioning' else 'zero', key)
    state = r.initial(key, length, n) if task['event'] == 'conditioning' else r.read_checkpoint(
        task['conditioning_checkpoint']['path'], key)[0]
    # Replay the exact original restart schedule, retaining its previously discarded history.
    elapsed = 0.
    history = None
    while elapsed < old_time:
        duration = min(1000. if history is None else 100., old_time - elapsed)
        tr = r.run(key, length, state, cue, duration, interval=5, n=n)
        state = tr.state
        history = r.concatenate(history, tr, elapsed)
        elapsed += duration
    original = {k: v[-1] for k, v in fields(old_trace, key, length).items()}
    scale = scales(fields(history, key, length), original, key, length)
    replay_error = distance({k: v[-1] for k, v in fields(history, key, length).items()}, original, scale)
    # Continue the actual saved state, not a replay approximation or a different branch.
    state = old_state
    checks = []
    consecutive = 0
    elapsed = 0.
    while elapsed < LIMIT:
        tr = r.run(key, length, state, cue, CHUNK, interval=5, n=n)
        state = tr.state
        elapsed += CHUNK
        metric = check(tr, key, length, scale)
        consecutive = consecutive + 1 if metric['passed'] else 0
        checks.append(dict(additional_time=elapsed, **metric))
        history = r.concatenate(history, tr, old_time + elapsed - CHUNK)
        print(key, length, task['event'], elapsed, metric['passed'], flush=True)
        if consecutive >= 2:
            break
    converged = consecutive >= 2
    profile = {k: v[-1].copy() for k, v in fields(tr, key, length).items()}
    discrepancy = distance(original, profile, scale)
    h = fields(history, key, length)
    errors = np.zeros(len(history.time))
    for k, y in h.items():
        d = abs(y - profile[k]).reshape(len(y), -1).max(axis=1) / scale[k]
        errors = np.maximum(errors, d)
    outside = np.flatnonzero(errors > ARRIVAL)
    arrival = (float(history.time[outside[-1] + 1]) if len(outside) and outside[-1] + 1 < len(errors)
               else 0. if not len(outside) else None) if converged else None
    # Require a complete confirmation window after entry and accurate numerical replay.
    if arrival is not None and (history.time[-1] - arrival < CHUNK or max(replay_error.values()) > TOL):
        arrival = None
    stem = f'{slot:03d}_{key}_{task["event"]}_L{length:g}'
    r.save_checkpoint(CAMPAIGN / 'profiles' / (stem + '.npz'), key, state, tr)
    np.savez_compressed(CAMPAIGN / 'histories' / (stem + '.npz'), time=history.time, x=history.x,
                        maximum_scaled_distance=errors, **{'field_' + k: v for k, v in h.items()})
    record = dict(slot=slot, model=key, length=length, event=task['event'],
        converged=converged, old_settled=task['phase']['settled'],
        old_time=task['phase']['time_to_settle'], old_elapsed=old_time,
        continuation_time=elapsed, arrival_time=arrival,
        old_minus_arrival=old_time - arrival if arrival is not None else None,
        maximum_scaled_endpoint_difference=max(discrepancy.values()),
        endpoint_difference_by_field=discrepancy, fixed_scales=scale,
        maximum_scaled_replay_error=max(replay_error.values()), replay_error_by_field=replay_error,
        checks=checks, profile='profiles/' + stem + '.npz', history='histories/' + stem + '.npz',
        generated_at=r.now(), slurm_job=os.environ.get('SLURM_JOB_ID'),
        reference_status='verified_continuation' if converged else 'unresolved_at_cap')
    r.write(CAMPAIGN / 'records' / (stem + '.json'), record)


def collect():
    manifest = verify()
    records = [json.loads(p.read_text()) for p in sorted((CAMPAIGN / 'records').glob('*.json'))]
    if sorted(p['slot'] for p in records) != list(range(len(manifest['tasks']))):
        raise RuntimeError(f'Incomplete: {len(records)}/{len(manifest["tasks"])}')
    rows = [{k: v for k, v in p.items() if not isinstance(v, (list, dict))} for p in records]
    with (CAMPAIGN / 'comparison.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summaries = []
    for key in r.MODELS:
        for event in ['conditioning', 'release']:
            rr = [p for p in records if p['model'] == key and p['event'] == event]
            verified = [p for p in rr if p['converged']]
            timed = [p for p in verified if p['arrival_time'] is not None]
            summaries.append(dict(model=key, event=event, total=len(rr), converged=len(verified),
                maximum_endpoint_difference=max((p['maximum_scaled_endpoint_difference'] for p in verified), default=None),
                old_time_range=[min(p['old_elapsed'] for p in rr), max(p['old_elapsed'] for p in rr)],
                arrival_time_range=[min(p['arrival_time'] for p in timed), max(p['arrival_time'] for p in timed)] if timed else None,
                median_time_reduction=float(np.median([p['old_minus_arrival'] for p in timed])) if timed else None))
    r.write(CAMPAIGN / 'summary.json', dict(generated_at=r.now(), complete=True, events=len(records),
        verified=sum(p['converged'] for p in records), summaries=summaries,
        maximum_replay_error=max(p['maximum_scaled_replay_error'] for p in records),
        unresolved=[{k: p[k] for k in ['model', 'length', 'event', 'old_settled', 'continuation_time']} for p in records if not p['converged']]))
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['setup', 'worker', 'collect'])
    parser.add_argument('--slot', type=int)
    args = parser.parse_args()
    {'setup': setup, 'worker': lambda: worker(args.slot), 'collect': collect}[args.action]()
