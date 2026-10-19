"""Read-only collection of the submitted scoped F4 replication.

Reuse all576 audited original cases, append completed new native cases, and
write atomic timestamped tables/HTML. No native submission or resubmission.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import collections
import hashlib
import importlib.util
import json
import os
import re
import signal
import subprocess
import sys
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
sys.dont_write_bytecode = True
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'polarity1d/manuscript').is_dir())
MS = ROOT / 'polarity1d/manuscript'
PILOT = MS / 'data/runs/revision_noise_strength_broad_20261019'
RUN = MS / 'data/runs/revision_noise_strength_replication_scoped_20261019'
SOURCE = PILOT / 'results/monitor/collector_v2.py'
SOURCE_SHA = '79747ebe6d831571e302a1fec528b7c6f58b223d7c9d7a8ea827b22a134fdd83'
MANIFEST_SHA = '57d461508929b4d5957b408c4a67db5616f0c0713f576562ce4834e5516ca3e9'
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


assert sha(SOURCE) == SOURCE_SHA
spec = importlib.util.spec_from_file_location('accepted_f4_readouts', SOURCE)
base = importlib.util.module_from_spec(spec); spec.loader.exec_module(base)


def scheduler(run, count):
    handle = base.read(run / 'submission.json')['job_id']
    states = {}; error = None
    try:
        command = ['sacct', '-j', handle, '--format=JobID,State,ExitCode,ElapsedRaw',
                   '--noheader', '--parsable2']
        response = subprocess.run(command, capture_output=True, text=True, timeout=30, check=True)
        for line in response.stdout.splitlines():
            fields = line.split('|')
            match = re.fullmatch(re.escape(handle) + r'_(\d+)', fields[0])
            expanded = re.fullmatch(re.escape(handle) + r'_\[(\d+)-(\d+)\]', fields[0])
            slots = [int(match[1])] if match else range(int(expanded[1]), int(expanded[2])+1) if expanded else []
            for slot in slots:
                if slot < count:
                    states[slot] = dict(state=fields[1].split()[0], exit_code=fields[2], elapsed_seconds=int(fields[3]))
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    return dict(queried_at=now(), job_id=handle, states=states, error=error)


def memory_groups(rows):
    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[tuple(row[k] for k in ('model', 'size', 'fraction', 'kind', 'window', 'lag'))].append(row)
    output = []
    for key, group in sorted(grouped.items()):
        eligible = sum(r['pairs_eligible'] for r in group)
        total = sum(r['pairs_available'] for r in group)
        reported = eligible >= 100 and eligible/total >= .01
        cosine = sum(r['direction_cosine_sum'] for r in group)/eligible if eligible else np.nan
        supported = [r['direction_correlation'] for r in group if r['support']]
        energy = np.sqrt(sum(r['energy_a'] for r in group)*sum(r['energy_b'] for r in group))
        output.append(dict(zip(('model', 'size', 'fraction', 'kind', 'window', 'lag'), key),
            cases=len(group), direction_pairs=eligible, all_pairs=total, direction_coverage=eligible/total,
            seeds_with_pairs=sum(r['pairs_eligible'] > 0 for r in group),
            seeds_with_declared_support=sum(r['support'] for r in group), group_support=reported,
            direction_cosine_pooled=cosine, direction_cosine_reported=cosine if reported else np.nan,
            supported_seed_cosine_min=min(supported) if supported else np.nan,
            supported_seed_cosine_max=max(supported) if supported else np.nan,
            raw_vector_correlation=sum(r['raw_product_sum'] for r in group)/energy if energy else np.nan))
    return output


def report(out, status, conditions):
    lines = ['<!doctype html><html><head><meta charset="utf-8"><title>F4 scoped replication: observed progress</title>',
        '<style>body{font:16px system-ui;max-width:1250px;margin:2rem auto;padding:0 1rem;color:#263640}table{border-collapse:collapse;width:100%;font-size:14px;margin:1rem 0 2rem}th,td{padding:.4rem;border-bottom:1px solid #d8dfe4;text-align:right}td:first-child,th:first-child{text-align:left}.notice{padding:1rem;background:#fff1d2;line-height:1.5}</style></head><body>',
        '<h1>Figure 4: common-strength replication progress</h1>',
        f"<p>Data generated {status['generated_at']} · native last observed {status['data_observed_at']}: {status['producer']} · collector {status['collector_state']} (job {status['collector_job_id']}).</p>",
        f"<p>{status['new_completed_outputs']}/360 new native cases complete; {status['scientific_cases_available']}/936 scientific cases available including576 exact retained cases. Native failures: {status['failures']}.</p>",
        '<p class="notice">This is a saved observation, updated by a separate read-only collector while it runs. Open or refresh manually; no timed refresh. New f=0.003/0.03/0.1 groups approach12 seeds as cases finish; the other broad endpoints remain2. Read actual n, rather than treating partial groups as final frequencies. The selected Spring upper-f run-q90 gate remains failed. Only the separately accepted ten-metric scope supports this replication; durations are sampled-cadence descriptors, without continuous lifetime or event-rate claims.</p>',
        '<p><a href="status.json">Observed status</a> · <a href="condition_window_summaries.csv">All condition windows</a> · <a href="memory_all_conditions_windows_lags.csv">Direction/coverage</a> · <a href="case_window_descriptors.csv">Case descriptors</a> · <a href="worker_inventory.csv">Every worker</a></p>',
        '<p>RAW regions: connected arcs above current min plus half current contrast. Qualified: multiple raw regions AND contrast≥10% of basal contrast. Low contrast: below that10% threshold. Fractions retain all saved scalar times. Mean positive input and fluctuation amplitude scale together.</p>']
    for model, name in base.NAMES.items():
        for size in (1., 2., 3.) if model == 'spring_ring' else (20., 40., 60.):
            lines.append(f'<h2>{name}, {"C₀" if model == "spring_ring" else "C"}={size:g}</h2><table><tr><th>f</th><th>Actual/target n</th><th>RAW multi</th><th>Qualified multi</th><th>Low contrast</th><th>Mean R</th><th>Mean A</th><th>Mean Q</th></tr>')
            for r in conditions:
                if r['model'] == model and r['size'] == size and r['kind'] == 'noise' and r['window'] == 'last_half':
                    values = ''.join(f'<td>{100*r[k]:.2f}%</td>' for k in ('multiregion_fraction', 'qualified_multiregion_fraction', 'low_contrast_fraction'))
                    values += ''.join(f'<td>{r[k]:.3f}</td>' for k in ('mean_R', 'mean_A', 'mean_Q'))
                    lines.append(f"<tr><td>{r['fraction']:g}</td><td>{r['available_cases']}/{r['target_cases']}</td>{values}</tr>")
            lines.append('</table>')
    lines.append('<p>Displayed table uses second-half windows: chemicals16000–32000, Spring1600–3200. This is not proof of stationarity. Full/half tables and native profile provenance remain downloadable. Activity direction is not cell migration; native laws/clocks are uncalibrated across models.</p></body></html>')
    target = out / 'index.html'; temporary = target.with_suffix('.html.tmp')
    temporary.write_text('\n'.join(lines)); temporary.replace(target)


def collect(run, final_state=None):
    run = Path(run).resolve(); assert sha(run / 'manifest.json') == MANIFEST_SHA
    manifest = base.read(run / 'manifest.json'); out = run / 'results/partial'
    out.mkdir(parents=True, exist_ok=True)
    path = out / 'case_cache.json'
    cache = base.read(path) if path.exists() else base.read(PILOT / 'results/partial/case_cache.json')
    assert len(cache) >= 576
    retained = base.rows(run / manifest['reuse']['index'])
    expected = {(r['trace_sha256'], r['primary_sha256']) for r in retained}
    actual = {(r['trace_sha256'], r['profile_sha256']) for r in cache.values()
              if not r['descriptors'][0]['case_id'].startswith('replication_')}
    assert len(expected) == 576 and actual == expected
    observed = scheduler(run, len(manifest['jobs'])); inventory = []
    changed = not path.exists()
    for job in manifest['jobs']:
        folder = run / 'workers' / f"{job['slot']:04d}"
        native = base.read(folder / 'status.json') if (folder / 'status.json').exists() else {}
        state = observed['states'].get(job['slot'], {})
        inventory.append(dict(slot=job['slot'], model=job['model'], size=job['size'], fraction=job['fraction'],
            seed=job['seed'], replicate=job['replicate'], native_status=native.get('status', 'not_started'),
            producer=native.get('producer', 'not_started'), complete=native.get('complete', False),
            scheduler_state=state.get('state', 'not_reported'), scheduler_exit_code=state.get('exit_code'),
            elapsed_seconds=native.get('elapsed_seconds', state.get('elapsed_seconds')),
            error=native.get('error', native.get('message'))))
        if not native.get('complete'):
            continue
        assert native['manifest_sha256'] == MANIFEST_SHA and native['job'] == job
        case_id = f"replication_{job['slot']:04d}"
        hashes = native['output_sha256']
        if case_id not in cache:
            assert sha(folder / 'trace.npz') == hashes['trace.npz']
            assert sha(folder / 'primary.npz') == hashes['primary.npz']
            identity = dict(case_id=case_id, model=job['model'], platform=job['platform'],
                size=job['size'], fraction=job['fraction'], kind='noise', seed=job['seed'], legacy=False)
            item = base.summarize_case(identity, folder / 'trace.npz', folder / 'primary.npz')
            item.update(trace=str(folder / 'trace.npz'), trace_sha256=hashes['trace.npz'],
                        profile=str(folder / 'primary.npz'), profile_sha256=hashes['primary.npz'])
            cache[case_id] = base.clean(item); changed = True
        else:
            assert cache[case_id]['trace_sha256'] == hashes['trace.npz'] and cache[case_id]['profile_sha256'] == hashes['primary.npz']
    if changed:
        base.write(path, cache)
    descriptors = [r for case in cache.values() for r in case['descriptors']]
    memory = [r for case in cache.values() for r in case['memory']]
    excluded = {'size', 'fraction', 'seed', 'legacy', 'observations', 'scalar_cadence', 'first_time', 'last_time'}
    metrics = sorted({k for r in descriptors for k, v in r.items()
                      if isinstance(v, (float, int)) and not isinstance(v, bool) and k not in excluded})
    grouped = collections.defaultdict(list)
    for r in descriptors:
        grouped[tuple(r[k] for k in ('model', 'size', 'fraction', 'kind', 'window'))].append(r)
    conditions = []
    for key, group in sorted(grouped.items()):
        model, size, f, kind, window = key
        target = 12 if kind == 'noise' and f in (.001, .003, .01, .03, .1) else 2 if kind == 'noise' else 1
        row = dict(zip(('model', 'size', 'fraction', 'kind', 'window'), key),
            available_cases=len(group), target_cases=target,
            coverage='complete' if len(group) == target else 'partial',
            inference='independent_seed_observations' if kind == 'noise' else 'deterministic_control')
        for metric in metrics:
            values = [r[metric] for r in group if r.get(metric) is not None]
            row[metric] = float(np.mean(values)) if values else np.nan
            row[metric + '_seed_min'] = min(values) if values else np.nan
            row[metric + '_seed_max'] = max(values) if values else np.nan
        conditions.append(row)
    base.csvwrite(out / 'case_window_descriptors.csv', descriptors)
    base.csvwrite(out / 'case_window_memory.csv', memory)
    base.csvwrite(out / 'condition_window_summaries.csv', conditions)
    base.csvwrite(out / 'memory_all_conditions_windows_lags.csv', memory_groups(memory))
    base.csvwrite(out / 'worker_inventory.csv', inventory)
    failed_states = {'FAILED', 'TIMEOUT', 'OUT_OF_MEMORY', 'CANCELLED', 'NODE_FAIL', 'PREEMPTED'}
    failures = [r for r in inventory if r['native_status'] == 'failed' or r['scheduler_state'] in failed_states]
    terminal = all(r['producer'] == 'terminal' or r['scheduler_state'] in failed_states for r in inventory)
    terminal = terminal and not any(r['scheduler_state'] in {'RUNNING', 'PENDING', 'COMPLETING', 'CONFIGURING'} for r in inventory)
    pending = all(r['native_status'] == 'not_started' and r['scheduler_state'] == 'PENDING' for r in inventory)
    status = dict(generated_at=now(), data_observed_at=observed['queried_at'],
        producer='terminal_with_failures' if terminal and failures else 'terminal' if terminal else 'pending' if pending else 'running',
        collector_state=final_state or ('stopped_native_terminal' if terminal else 'running' if os.environ.get('SLURM_JOB_ID') else 'one_shot'),
        collector_job_id=os.environ.get('SLURM_JOB_ID'), complete=terminal and not failures,
        scientific_cases_available=len(cache), scientific_cases_target=936, retained_cases=576,
        new_completed_outputs=sum(r['complete'] for r in inventory), requested_new_outputs=360,
        failures=len(failures), scheduler=observed, manifest_sha256=MANIFEST_SHA,
        collector_sha256=sha(__file__), accepted_readout_source_sha256=SOURCE_SHA,
        qualification='Separate root-accepted ten-metric scope; original Spring full gate remains failed for sampled run-q90',
        no_native_submission=True)
    report(out, status, conditions)
    base.write(out / 'status.json', status)
    return status


def main():
    global STOP
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, default=RUN)
    p.add_argument('--loop', action='store_true')
    p.add_argument('--interval', type=float, default=90.)
    args = p.parse_args()
    assert not args.loop or os.environ.get('SLURM_JOB_ID'), 'Long collection requires SLURM'
    def stop(*_):
        global STOP
        STOP = True
    signal.signal(signal.SIGTERM, stop); signal.signal(signal.SIGUSR1, stop)
    try:
        while True:
            result = collect(args.run, None if args.loop else 'one_shot_complete')
            print(json.dumps({k: v for k, v in result.items() if k != 'scheduler'}), flush=True)
            if result['producer'].startswith('terminal') or not args.loop or STOP:
                break
            for _ in range(max(1, int(args.interval))):
                if STOP:
                    break
                time.sleep(1)
        if STOP and not result['producer'].startswith('terminal'):
            collect(args.run, 'stopped_signal')
    except BaseException as exc:
        path = args.run / 'results/partial/status.json'
        result = base.read(path) if path.exists() else dict(generated_at=now(), producer='not_observed')
        result.update(collector_state='stopped_error', collector_stopped_at=now(),
                      collector_error=f'{type(exc).__name__}: {exc}')
        base.write(path, result)
        raise


if __name__ == '__main__':
    main()
