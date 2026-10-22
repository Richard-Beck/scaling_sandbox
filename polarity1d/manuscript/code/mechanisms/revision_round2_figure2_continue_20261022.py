"""Terminal-only exact phase continuation contingency for Holmes gap noise.

No new integration is launched by inspect or prepare. A fresh immutable run is
prepared only after the original array is terminal. Completed phase archives
are referenced unchanged; only unsaved phases are integrated. The original
native adapter already restarts its BDF solver at every phase boundary.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'polarity1d/manuscript').is_dir())
PARENT = ROOT / 'polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022'
DEFAULT = ROOT / 'polarity1d/manuscript/data/runs/signature_holmes_retention_continuation_round2_20261022'
PYTHON = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'


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
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def frozen_imports(run):
    sys.path.insert(0, str(run / 'frozen_sources/polarity1d/manuscript/code/mechanisms'))
    import signature_campaign_20261020 as campaign
    import signature_noise_campaign_20261020 as noise_campaign
    from signature_noise_20261020 import LinearGapNoise, ShiftedInput
    from signature_analysis_20261020 import pulse_gap_endpoints
    native, inputs, readouts, _ = campaign.modules()
    return campaign, noise_campaign, LinearGapNoise, ShiftedInput, pulse_gap_endpoints, native, inputs, readouts


def verify(run):
    manifest = read(run / 'manifest.json')
    assert sha(run / 'manifest.json') == read(run / 'freeze.json')['manifest_sha256']
    assert sha(run / 'worker.py') == manifest['worker_sha256']
    for relative, digest in {**manifest['input_sha256'], **manifest['source_sha256']}.items():
        assert sha(run / relative) == digest, relative
    return manifest


def parent_terminal():
    result = subprocess.run(['squeue', '--jobs=21279666', '--noheader', '--format=%i'],
                            capture_output=True, text=True, check=True)
    if result.stdout.strip():
        raise RuntimeError('Original array is live: refuse continuation preparation or integration')
    history = subprocess.run(['sacct', '--jobs=21279666', '--noheader', '--parsable2', '-X',
                              '-o', 'JobID,State,ExitCode'], capture_output=True, text=True, check=True)
    records = [line.split('|') for line in history.stdout.splitlines() if line.strip()]
    terminal = {'COMPLETED', 'FAILED', 'TIMEOUT', 'CANCELLED', 'OUT_OF_MEMORY', 'NODE_FAIL',
                'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
    assert len(records) == 32 and all(r[1].split()[0] in terminal for r in records)
    return records


def inspect():
    manifest = verify(PARENT)
    rows = []
    for case in manifest['cases']:
        if case['kind'] != 'noise':
            continue
        phases = sorted((PARENT / 'workers' / f"{case['task']:04d}" / 'trial').glob('phase_[0-9][0-9][0-9].npz'))
        details = []
        for path in phases:
            with np.load(path, allow_pickle=False) as data:
                details.append(dict(path=str(path), duration=float(data['time'][-1]),
                                    complete_native_final_state_shape=list(data['state'].shape),
                                    field_names=[key for key in data if key.startswith('field_')],
                                    exact_restart_flag=json.loads(str(data['restart_check_json']))['passed']))
        rows.append(dict(task=case['task'], seed=case['seed'], completed_phases=details))
    report = dict(continuation_possible=True, integration_launched=False, cases=rows,
                  forcing_reconstruction='Identical PCG64 seed, whole-train end time and array shape; identical frozen spatial/temporal lfilter and variance-normalized interpolation. No draws during RHS evaluation.',
                  solver_restart='Exact saved full state at a boundary where the original adapter already restarts BDF.',
                  activation_guard='Prepare and worker both require all original array tasks terminal; no duplicate live integration.')
    write(PARENT / 'results/continuation_preflight.json', report)
    print(json.dumps(dict(continuation_possible=True, cases=len(rows), integration_launched=False)))


def prepare(out):
    states = parent_terminal()
    old = verify(PARENT)
    missing = [case['task'] for case in old['cases']
               if not (PARENT / 'workers' / f"{case['task']:04d}" / 'result.json').is_file()]
    if not missing:
        raise RuntimeError('Original campaign is complete; no continuation is necessary')
    if out.exists():
        raise FileExistsError(out)
    for directory in ('inputs', 'workers', 'logs', 'results'):
        (out / directory).mkdir(parents=True, exist_ok=True)
    shutil.copytree(PARENT / 'frozen_sources', out / 'frozen_sources')
    for relative in old['input_sha256']:
        shutil.copy2(PARENT / relative, out / relative)
    shutil.copy2(Path(__file__), out / 'worker.py')
    reused = {}
    for task in missing:
        source = PARENT / 'workers' / f'{task:04d}' / 'trial'
        files = sorted(source.glob('phase_[0-9][0-9][0-9].npz'))
        assert [path.name for path in files] == [f'phase_{i:03d}.npz' for i in range(len(files))]
        reused[str(task)] = [dict(path=str(path.resolve()), sha256=sha(path)) for path in files]
    manifest = dict(old, worker_sha256=sha(out / 'worker.py'), status='prepared_continuation_not_submitted',
                    continuation=dict(parent=str(PARENT), parent_manifest_sha256=sha(PARENT / 'manifest.json'),
                                      planned_tasks=missing, reused_phases=reused, terminal_scheduler_records=states,
                                      forcing_identity='Exact regeneration from original full-field PCG64 seed and shape'))
    write(out / 'manifest.json', manifest)
    write(out / 'freeze.json', dict(manifest_sha256=sha(out / 'manifest.json')))
    worker = out / 'worker.sh'
    worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '
                      + PYTHON + ' -B ' + str(out / 'worker.py') + ' worker --out ' + str(out)
                      + ' --task "$SLURM_ARRAY_TASK_ID"\n')
    worker.chmod(0o750)
    command = ['sbatch', '--parsable', '--partition=red', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
               '--mem=8G', '--time=04:00:00', '--array=' + ','.join(map(str, missing)),
               '--job-name=F2-Holmes-continuation', '--output=' + str(out / 'logs/%A_%a.out'),
               '--error=' + str(out / 'logs/%A_%a.err'), str(worker)]
    import shlex
    (out / 'submission_command.txt').write_text(shlex.join(command) + '\n')
    print(json.dumps(dict(prepared=str(out), tasks=missing, submitted=False)))


def worker(out, task):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Native continuation requires SLURM')
    parent_terminal()
    manifest = verify(out)
    assert sha(PARENT / 'manifest.json') == manifest['continuation']['parent_manifest_sha256']
    assert task in manifest['continuation']['planned_tasks']
    assert not (PARENT / 'workers' / f'{task:04d}' / 'result.json').exists()
    case = next(case for case in manifest['cases'] if case['task'] == task)
    campaign, nc, Noise, Shift, endpoints_fn, native, inputs, readouts = frozen_imports(out)
    started = time.monotonic()
    folder = out / 'workers' / f'{task:04d}'
    if (folder / 'result.json').exists():
        raise FileExistsError(folder)
    trial_folder = folder / 'trial'
    trial_folder.mkdir(parents=True, exist_ok=True)
    state, basal = native.read_checkpoint(out / case['checkpoint'], case['model'])
    schedule = inputs.make_train(case['calibration_duration'], 1.5, 4, case['gap'])
    latent = Noise(case['native_length'], 201, schedule.end_time, case['seed'], case['mean'],
                   ell=1.9, tau=10, knot=2.5, sigma=1, kind=case['kind'], transform='bounded_tanh')
    forcing = nc.AppliedInputAudit(latent, 1.)
    forcing_audit = latent.audit()
    cue = inputs.AreaCue(case['cue_area'], case['cue_width'], grid_kind='cell')
    zero = inputs.AreaCue(0., case['cue_width'], grid_kind='cell')
    phases = [(f'{i:02d}_{kind}', start, end, cue if kind == 'on' else Shift(forcing, start))
              for i, (kind, start, end) in enumerate(schedule.segments)]
    phases.append(('recovery', schedule.end_time, schedule.end_time + case['recovery_duration'], zero))
    reused = manifest['continuation']['reused_phases'][str(task)]
    traces, events = [], []
    for index, (name, start, stop, stimulus) in enumerate(phases):
        archive = trial_folder / f'phase_{index:03d}.npz'
        if index < len(reused):
            record = reused[index]
            assert sha(record['path']) == record['sha256']
            archive.symlink_to(record['path'])
            _, trace = native.read_checkpoint(archive, case['model'])
            assert trace.time[0] == 0 and trace.time[-1] == stop - start
            trace.restart_check = native.restart_error(state, trace)
            assert trace.restart_check['passed']
            artifact = dict(path=str(archive), sha256=record['sha256'], frames=len(trace.time),
                            fields=list(trace.fields), full_native_fields=True)
        else:
            trace = native.run_native(case['model'], case['native_length'], state, stimulus,
                                      native.sample_times(case['model'], stop - start), cells=201, max_step=1.25)
            artifact = native.save_phase(archive, trace)
        event = dict(index=index, name=name, start=start, stop=stop, frames=len(trace.time),
                     zero_duration=False, restart_check=trace.restart_check, archive=artifact)
        if hasattr(stimulus, 'profile_diagnostics'):
            event['cue_initial_final'] = [stimulus.profile_diagnostics(trace.x, float(trace.physical_length[i])) for i in (0, -1)]
        events.append(event)
        traces.append(trace)
        state = trace.state
    trace = native.concatenate(traces)
    native.save_phase(trial_folder / 'trace.npz', trace)
    measurement = campaign._analyze(trace, basal, case, readouts)
    tail = native.native_field_tail_gate(trace, baseline_scales=case['fixed_scales'])
    write(trial_folder / 'readouts.json', measurement)
    write(trial_folder / 'native_tail.json', tail)
    phase_endpoints = []
    times = np.asarray(measurement['series']['time'])
    for event in events:
        for label in ('start', 'stop'):
            indices = np.flatnonzero(np.isclose(times, event[label], atol=1e-8, rtol=0))
            assert len(indices) == 1
            i = int(indices[0])
            phase_endpoints.append(dict(phase=event['name'], endpoint=label, time=event[label], scalar_sample_index=i,
                                        readouts={key: values[i] for key, values in measurement['series'].items()}))
    on = [[event['start'], event['stop']] for event in events if event['name'].endswith('_on')]
    recovery = [events[-1]['start'], events[-1]['stop']]
    write(trial_folder / 'phase_evidence.json', dict(on_intervals=on, recovery_interval=recovery, endpoints=phase_endpoints))
    trial = dict(final_d=measurement['final_d'], readout_eligible=measurement['readout_eligible'],
                 ineligibility_reason=measurement['ineligibility_reason'], field_tail_pass=tail['field_tail_pass'],
                 field_tail=tail, schedule=schedule.as_dict(), events=events, trial_directory=str(trial_folder.relative_to(out)),
                 phase_on_intervals=on, recovery_interval=recovery, phase_endpoint_readouts=phase_endpoints,
                 readouts='readouts.json', full_trace_sha256=sha(trial_folder / 'trace.npz'))
    forcing_audit['actual_queries'] = forcing.audit()
    forcing_audit['actual_queries']['scope'] += '; continuation phases only; original completed-phase queries are unavailable'
    forcing_audit['hard_bound_pass'] = forcing.maximum <= .04 + 1e-14 if forcing.values else True
    forcing_audit['hard_bound_scope'] = 'Analytic bounded-tanh map establishes [0,.04] over the entire original realization; resumed queries separately verified.'
    assert forcing_audit['hard_bound_pass']
    artifacts = {str(path.relative_to(out)): sha(path) for path in trial_folder.rglob('*') if path.is_file()}
    write(folder / 'result.json', dict(case=case, trial=trial, endpoints=endpoints_fn(measurement['series'], events),
          forcing=forcing_audit, wall_seconds=time.monotonic() - started, manifest_sha256=sha(out / 'manifest.json'),
          artifacts=artifacts, continuation=manifest['continuation']['reused_phases'][str(task)]))
    print(json.dumps(dict(task=task, complete=True, phases_reused=len(reused))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['inspect', 'prepare', 'worker'])
    parser.add_argument('--out', type=Path, default=DEFAULT)
    parser.add_argument('--task', type=int)
    args = parser.parse_args()
    if args.command == 'inspect':
        inspect()
    elif args.command == 'prepare':
        prepare(args.out.resolve())
    else:
        worker(args.out.resolve(), args.task)
