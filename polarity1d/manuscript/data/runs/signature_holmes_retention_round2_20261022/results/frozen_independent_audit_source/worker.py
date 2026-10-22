"""Independent reconstruction of the frozen Holmes retention campaign.

This reader runs no integration and never changes a frozen input or worker.
One SLURM audit worker per trajectory keeps full-field checks bounded in RAM.
"""
from pathlib import Path
from datetime import datetime, timezone
from types import SimpleNamespace
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
RUN = ROOT / 'polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022'


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def verify_freeze():
    manifest = read(RUN / 'manifest.json')
    assert sha(RUN / 'manifest.json') == read(RUN / 'freeze.json')['manifest_sha256']
    assert sha(RUN / 'worker.py') == manifest['worker_sha256']
    for relative, digest in {**manifest['source_sha256'], **manifest['input_sha256']}.items():
        assert sha(RUN / relative) == digest, relative
    return manifest


def worker(task):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Full native-array reconstruction requires SLURM')
    manifest = verify_freeze()
    frozen = RUN / 'frozen_sources/polarity1d/manuscript/code/mechanisms'
    sys.path.insert(0, str(frozen))
    import signature_analysis_20261020 as analysis
    import signature_native_20261020 as native
    import signature_readouts_20261020 as readouts
    case = manifest['cases'][task]
    result_path = RUN / 'workers' / f'{task:04d}' / 'result.json'
    result = read(result_path)
    assert result['case'] == case
    assert result['manifest_sha256'] == sha(RUN / 'manifest.json')
    folder = result_path.parent / 'trial'
    actual_files = {str(p.relative_to(RUN)) for p in folder.rglob('*') if p.is_file()}
    assert actual_files == set(result['artifacts'])
    for relative, digest in result['artifacts'].items():
        assert sha(RUN / relative) == digest, relative
    mapping = analysis._trace_mapping(folder / 'trace.npz', case['model'])
    basal = analysis._trace_mapping(RUN / case['checkpoint'], case['model'])
    regenerated = readouts.analyze_trace(mapping, analysis._readout_parameters(case, basal), basal)
    measurement = read(folder / 'readouts.json')
    analysis._same_numeric_tree(regenerated, measurement)
    endpoints = analysis.pulse_gap_endpoints(regenerated['series'], result['trial']['events'])
    analysis._same_numeric_tree(endpoints, result['endpoints'])
    state_fields = tuple('field_' + name for name in native.model_for(case['model']).state_names)
    phase_audit = analysis.audit_native_phase_join(folder / 'trace.npz', result['trial']['events'],
                                                 mapping, state_fields=state_fields)
    assert phase_audit['passed']
    # Full initial native state must equal the exact copied basal, not only Rac.
    with np.load(folder / 'phase_000.npz', allow_pickle=False) as trace, np.load(RUN / case['checkpoint'], allow_pickle=False) as start:
        restart_checks = {}
        for field in state_fields:
            wanted, actual = start[field][-1], trace[field][0]
            allowance = 32 * np.finfo(float).eps * np.maximum(1., np.abs(wanted))
            assert np.all(np.abs(actual - wanted) <= allowance), field
            restart_checks[field] = float(np.max(np.abs(actual - wanted)))
    # Reconstruct the independent full-field tail diagnostic from final phase.
    final_phase = Path(result['trial']['events'][-1]['archive']['path'])
    with np.load(final_phase, allow_pickle=False) as data:
        trace = SimpleNamespace(time=data['time'], fields={name[6:]: data[name] for name in state_fields},
                                model=case['model'], length=case['native_length'])
        tail = native.native_field_tail_gate(trace, baseline_scales=case['fixed_scales'])
    # Times in the standalone final phase are local; all other quantities agree.
    observed_tail = read(folder / 'native_tail.json')
    for key in ('tail_start', 'tail_stop'):
        tail[key] += result['trial']['events'][-1]['start']
    analysis._same_numeric_tree(tail, observed_tail)
    noise = result['forcing']
    if noise:
        assert noise['hard_bound_pass']
        maximum = noise['actual_queries']['maximum']
        if maximum is None:
            assert 'continuation' in result and 'hard_bound_scope' in noise
        else:
            assert maximum <= .04 + 1e-14
    row = dict(task=task, case=case, trial=result['trial'], endpoints=endpoints, forcing=noise,
               result_path=str(result_path), result_sha256=sha(result_path), campaign_root=str(RUN),
               execution_manifest_sha256=sha(RUN/'manifest.json'),
               worker_wall_seconds=result['wall_seconds'], scientific_artifacts=len(actual_files),
               restart_maximum_absolute_error=restart_checks, all_artifact_hashes_verified=True,
               all_readouts_independently_recomputed=True, all_native_phases_agree=True,
               field_tail_independently_recomputed=True, phase_audit=phase_audit, passed=True)
    write(RUN / 'results/independent_audit/workers' / f'{task:04d}.json', row)
    print(json.dumps(dict(task=task, passed=True)))


def collect(continuation=None):
    manifest = verify_freeze()
    scheduler = subprocess.run(['sacct', '--jobs=21279666', '--noheader', '--parsable2', '-X',
                                '-o', 'JobID,State,ExitCode'], check=True, text=True,
                               capture_output=True)
    states = [line.split('|') for line in scheduler.stdout.splitlines() if line.strip()]
    assert len(states) == 32
    if continuation is None:
        assert all(row[1] == 'COMPLETED' and row[2] == '0:0' for row in states)
    else:
        additional = read(continuation / 'manifest.json')
        assert sha(continuation/'manifest.json')==read(continuation/'freeze.json')['manifest_sha256']
        assert additional['continuation']['parent_manifest_sha256']==sha(RUN/'manifest.json')
        assert additional['cases']==manifest['cases']
        assert not subprocess.run(['squeue','--jobs=21279666','--noheader','--format=%i'],check=True,text=True,capture_output=True).stdout.strip()
    rows = []
    for case in manifest['cases']:
        path = RUN / 'results/independent_audit/workers' / f"{case['task']:04d}.json"
        if not path.exists() and continuation is not None:
            path = continuation / 'results/independent_audit/workers' / f"{case['task']:04d}.json"
        row = read(path)
        assert row['passed'] and row['case'] == case
        assert sha(row['result_path']) == row['result_sha256']
        rows.append(row)
    assert len(rows) == 32
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True, workers=32,
                  manifest_sha256=sha(RUN / 'manifest.json'),
                  audit_source_sha256=sha(Path(__file__)), all_native_phase_join_audits_pass=True,
                  all_readout_recomputations_pass=True, all_restart_checks_pass=True,
                  all_field_tail_recomputations_pass=True, all_artifact_hashes_verified=True,
                  producer_state='stopped', trajectories=rows)
    report['scheduler_observation'] = dict(job_id='21279666', records=states,
                                           observed_utc=report['generated_utc'])
    if continuation is not None:
        report['continuation_campaign']=dict(path=str(continuation),manifest_sha256=sha(continuation/'manifest.json'),
                                             original_phases_preserved=True,forcing_realization_unchanged=True)
    write(RUN / 'results/independent_native_readout_audit.json', report)
    allocated = sum(path.stat().st_blocks * 512 for path in RUN.rglob('*') if path.is_file())
    write(RUN / 'results/status.json', dict(generated_utc=report['generated_utc'], state='complete',
          producer_state='stopped', completed=32, total=32, independently_audited=32,
          scientific_audit_sha256=sha(RUN / 'results/independent_native_readout_audit.json'),
          allocated_bytes=allocated, storage_scope='Allocated file blocks under campaign root',
          native_job_id='21279666', audit_job_ids=['21279728', '21279740']))
    print(json.dumps(dict(passed=True, workers=32)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['worker', 'collect'])
    parser.add_argument('--task', type=int)
    parser.add_argument('--run',type=Path,default=RUN)
    parser.add_argument('--continuation-run',type=Path)
    args = parser.parse_args()
    RUN=args.run.resolve()
    worker(args.task) if args.command == 'worker' else collect(args.continuation_run.resolve() if args.continuation_run else None)
