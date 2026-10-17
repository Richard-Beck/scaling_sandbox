"""Isolated paired pilot: linear gradient cues for reversal and finite countercue.

Reuse accepted v8 physics, preparation gates and complete duration grids.
Only train/counter spatial cues are mapped to the existing left/right gradient.
The accepted runner, protocol and results are never rewritten.
"""
from pathlib import Path
import argparse
from collections import Counter
import copy
import csv
import gzip
import json
import os
import shlex
import subprocess
import sys

import numpy as np

SCRIPT = Path(__file__).resolve()
CODE = SCRIPT.parents[1]
sys.path.insert(0, str(CODE))
import runner as r

LENGTHS = {
    'wave_pinning': [3., 21., 39.],
    'otsuji': [3., 21., 39.],
    'goryachev': [3., 13., 15.],
    'legi': [1., 31., 33.],
    'holmes_model3': [10., 14., 24., 34., 50.],
    'debelly': [3., 9., 11., 18., 19., 29., 39.],
    'spring': [.6, 1.1, 1.6, 2.1, 2.6],
}
CHANGED = {'holmes_model3', 'debelly', 'spring'}
ORIGINAL_CUE = r.Cue
ORIGINAL_SOURCES = r.sources


class NormalizedCue(ORIGINAL_CUE):
    def __init__(self, kind, model, stop=None):
        super().__init__({'train': 'left', 'counter': 'right'}.get(kind, kind), model, stop)


def configure(work):
    r.HERE = work.resolve()
    r.Cue = NormalizedCue
    r.PROTOCOL = copy.deepcopy(r.PROTOCOL)
    for task in r.PROTOCOL['tasks']:
        included = (task['assay'] in ['gradient_reversal', 'finite_countercue']
                    and task['length'] in LENGTHS[task['model']])
        task['action'] = 'steady_state_rerun' if included else 'pilot_excluded'
    r.sources = sources


def selected():
    return [t for t in r.tasks() if t['action'] == 'steady_state_rerun']


def sources():
    files = [SCRIPT, r.HERE / 'pilot_protocol.json',
             CODE / 'results' / 'records.json.gz', CODE / 'results' / 'endpoints.csv']
    return {**ORIGINAL_SOURCES(), **{str(p.relative_to(r.ROOT)): r.sha(p) for p in files}}


def submit(action, array=None, dependency=None):
    executable = json.loads((CODE / 'environment.json').read_text())['python']
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; ' + shlex.join(
        [executable, str(SCRIPT), action, '--work', str(r.HERE)])
    if array is not None:
        command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
            '--mem=3G', '--time=02:00:00', '--job-name=cue-pilot-' + action,
            '--chdir=' + str(r.ROOT), '--output=' + str(r.HERE / 'logs' / (action + '_%A_%a.out')),
            '--error=' + str(r.HERE / 'logs' / (action + '_%A_%a.err'))]
    if array is not None:
        args.append('--array=' + array)
    if dependency:
        args.append('--dependency=afterok:' + dependency)
    return subprocess.check_output(args + ['--wrap', command], text=True).strip().split(';')[0]


def cue_checks():
    checks = []
    q = np.linspace(0, 1, 201)
    for model in r.MODELS:
        for kind, target in [('train', 'left'), ('counter', 'right')]:
            a = NormalizedCue(kind, model)(0, q, 1)
            b = ORIGINAL_CUE(target, model)(0, q, 1)
            assert np.array_equal(a, b)
            original = ORIGINAL_CUE(kind, model)(0, q, 1)
            if model not in CHANGED:
                assert np.array_equal(original, a)
            checks.append(dict(model=model, cue=kind, matches_gradient=True,
                               differs_from_baseline=not np.array_equal(original, a),
                               peak_amplitude=float(a.max())))
        assert np.array_equal(NormalizedCue('counter', model, stop=2)(2, q, 1), np.zeros_like(q))
    return checks


def launch():
    if (r.HERE / 'frozen_provenance.json').exists():
        raise RuntimeError('Pilot already initialized; use another workspace')
    for name in ['logs', 'preparations', 'checkpoints', 'data', 'trials', 'results']:
        (r.HERE / name).mkdir(parents=True, exist_ok=True)
    r.write(r.HERE / 'pilot_protocol.json', dict(
        generated_at=r.now(), tasks=selected(), lengths=LENGTHS, cue_checks=cue_checks(),
        normalization='Same existing linear-gradient shape and amplitude within each model; train->left and counter->right',
        changed_models=sorted(CHANGED), shared=r.PROTOCOL['shared'],
        duration_grids=r.PROTOCOL['duration_grids'],
        success='Every sampled baseline-successful endpoint remains successful; new ineligible/unsettled/censored readouts are losses; gains do not offset losses',
        retained='Native initialization, Spring spontaneous gradient preparation, finite no-cue preparation, all eligibility/settling gates, crossfade, horizons, dense duration grids and endpoint definitions'))
    hashes=sources()
    r.write(r.HERE / 'frozen_provenance.json', dict(created_at=r.now(), source_hashes=hashes,
            provenance_id=r.hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()))
    preparation=submit('prepare', array=f'0-{len(selected())-1}')
    controller=submit('stage2', dependency=preparation)
    record=dict(stage='preparations_running', preparation_job_id=preparation,
                controller_job_id=controller, case_count=len(selected()), launched_at=r.now())
    r.write(r.HERE / 'launch.json', record)
    print(json.dumps(record))


def stage2():
    r.verify()
    jobs=[]
    for task in selected():
        idx=task['task_index']; path=r.HERE / 'preparations' / f'task_{idx:04d}.json'
        if not path.exists():
            raise RuntimeError(f'Missing preparation: {idx}')
        prep=json.loads(path.read_text())
        if not prep['eligible']:
            r.write(r.HERE / 'data' / f'task_{idx:04d}.json', dict(**task,
                    status=prep['status'], eligible=False, endpoint=None))
            continue
        base=dict(task_index=idx, checkpoint_sha256=prep['checkpoint_sha256'])
        if task['assay']=='finite_countercue':
            durations=r.PROTOCOL['duration_grids'][task['model']]['refined_durations']
            for i in range(0, len(durations), 10):
                jobs.append(dict(**base, durations=durations[i:i+10]))
        else:
            jobs.append(base)
    for i, job in enumerate(jobs):
        job['slot']=i
    r.write(r.HERE / 'response_manifest.json', jobs)
    r.write(r.HERE / 'response_manifest_sha256.json', dict(sha256=r.sha(r.HERE / 'response_manifest.json')))
    response=submit('respond', array=f'0-{len(jobs)-1}') if jobs else None
    final=submit('summarize', dependency=response)
    record=json.loads((r.HERE / 'launch.json').read_text())
    record.update(stage='responses_running', response_job_id=response,
                  response_jobs=len(jobs), summarize_job_id=final)
    r.write(r.HERE / 'launch.json', record)
    print(json.dumps(record))


def summarize():
    r.verify()
    with gzip.open(CODE / 'results' / 'records.json.gz', 'rt') as f:
        baseline=json.load(f)
    base_preps={p['task_index']:p for p in baseline['preparations']}
    base_cases={int(p['task_index']):p for p in csv.DictReader((CODE / 'results' / 'endpoints.csv').open())}
    jobs=json.loads((r.HERE / 'response_manifest.json').read_text())
    results=[]; comparisons=[]; control_errors=[]
    for task in selected():
        idx=task['task_index']; prep=json.loads((r.HERE / 'preparations' / f'task_{idx:04d}.json').read_text())
        if task['assay']=='finite_countercue' and prep['eligible']:
            planned=[j for j in jobs if j['task_index']==idx]
            chunks=[json.loads((r.HERE / 'trials' / f'chunk_{j["slot"]:05d}.json').read_text()) for j in planned]
            tests=sorted([t for c in chunks for t in c['tests']], key=lambda t:t['duration'])
            expected=r.PROTOCOL['duration_grids'][task['model']]['refined_durations']
            if [t['duration'] for t in tests]!=expected:
                raise RuntimeError(f'Incomplete duration sweep: {idx}')
            passing=[t['duration'] for t in tests if t['success']]
            row=dict(**task, status='success' if passing else 'censored', eligible=True,
                     endpoint=min(passing) if passing else None, tests=tests, tested_count=len(tests))
        else:
            row=json.loads((r.HERE / 'data' / f'task_{idx:04d}.json').read_text())
        results.append(row)
        old=base_cases[idx]; bp=base_preps[idx]
        baseline_endpoint=float(old['endpoint']) if old['endpoint'] else None
        lost=old['mode']=='success' and row['status']!='success'
        gained=old['mode']!='success' and row['status']=='success'
        original_phases={p['cue']:p for p in bp['phases']}
        new_unsettled=[p['cue'] for p in prep['phases'] if not p['settled'] and original_phases.get(p['cue'],{}).get('settled')]
        comp=dict(model=task['model'], assay=task['assay'], length=task['length'], task_index=idx,
                  changed_model=task['model'] in CHANGED, baseline_status=old['mode'], pilot_status=row['status'],
                  baseline_endpoint=baseline_endpoint, pilot_endpoint=row['endpoint'], lost_readout=lost, gained_readout=gained,
                  baseline_preparation_time=bp['total_preparation_time'], pilot_preparation_time=prep['total_preparation_time'],
                  baseline_preparation_status=bp['status'], pilot_preparation_status=prep['status'],
                  newly_unsettled_phases=';'.join(new_unsettled),
                  endpoint_ratio=(row['endpoint']/baseline_endpoint if baseline_endpoint and row['endpoint'] is not None else None))
        comparisons.append(comp)
        # All gradient protocols and the four unchanged models are exact controls.
        if task['assay']=='gradient_reversal' or task['model'] not in CHANGED:
            if old['mode']!=row['status'] or baseline_endpoint!=row['endpoint'] or bp['total_preparation_time']!=prep['total_preparation_time']:
                control_errors.append(comp)
    r.csvwrite(r.HERE / 'comparison.csv', comparisons)
    r.write(r.HERE / 'pilot_records.json', dict(generated_at=r.now(), cases=results))
    groups=[]
    for model in r.MODELS:
        for assay in ['gradient_reversal','finite_countercue']:
            rows=[c for c in comparisons if c['model']==model and c['assay']==assay]
            groups.append(dict(model=model, assay=assay, n_cases=len(rows),
                 baseline_success=sum(c['baseline_status']=='success' for c in rows),
                 pilot_success=sum(c['pilot_status']=='success' for c in rows),
                 losses=sum(c['lost_readout'] for c in rows), gains=sum(c['gained_readout'] for c in rows),
                 baseline_status_counts=dict(Counter(c['baseline_status'] for c in rows)),
                 pilot_status_counts=dict(Counter(c['pilot_status'] for c in rows))))
    summary=dict(generated_at=r.now(), complete=True, passed=not control_errors and not any(c['lost_readout'] for c in comparisons),
                 total_cases=len(comparisons), baseline_success=sum(c['baseline_status']=='success' for c in comparisons),
                 pilot_success=sum(c['pilot_status']=='success' for c in comparisons),
                 lost_readouts=[c for c in comparisons if c['lost_readout']],
                 gained_readouts=[c for c in comparisons if c['gained_readout']],
                 control_mismatch_count=len(control_errors), control_mismatches=control_errors,
                 full_duration_tests=sum(len(c.get('tests',[])) for c in results), groups=groups,
                 scope='Sampled lengths only; complete original duration grid at every eligible finite preparation')
    r.write(r.HERE / 'summary.json', summary)
    record=json.loads((r.HERE / 'launch.json').read_text())
    record.update(stage='finished', finished_at=r.now(), passed=summary['passed'])
    r.write(r.HERE / 'launch.json', record)
    print(json.dumps(summary))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['launch','prepare','stage2','respond','summarize','status','check'])
    ap.add_argument('--work', type=Path, default=CODE / 'pilots' / 'linear_shared_20261001')
    ap.add_argument('--slot', type=int)
    args=ap.parse_args(); configure(args.work)
    if args.action=='launch': launch()
    elif args.action=='check': print(json.dumps(dict(cases=len(selected()), cue_checks=cue_checks())))
    elif args.action=='prepare': r.run_preparation(args.slot)
    elif args.action=='stage2': stage2()
    elif args.action=='respond':
        if r.sha(r.HERE / 'response_manifest.json')!=json.loads((r.HERE / 'response_manifest_sha256.json').read_text())['sha256']:
            raise RuntimeError('Response manifest changed')
        r.run_response(args.slot)
    elif args.action=='summarize': summarize()
    else: print((r.HERE / 'launch.json').read_text())


if __name__=='__main__':
    main()
