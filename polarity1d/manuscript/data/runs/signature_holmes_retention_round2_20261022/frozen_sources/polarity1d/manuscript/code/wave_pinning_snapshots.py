"""Continuous right-cue wave-pinning snapshots, with exact basal preparations."""
from __future__ import annotations

import argparse
import copy
import csv
import fcntl
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import traceback

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE / 'core'))
import numpy as np
import native as n
import polarity_establishment as e

DEFAULT = n.DATA / 'runs/wave_pinning_continuous_snapshots_20261005'
BASE = n.DATA / 'runs/corrected_reversal_20261003/manifest.json'
TEMPLATE = CODE / 'wave_pinning_snapshots.html'


def read(path):
    return json.loads(path.read_text())


def sources():
    return {**e.sources(), **{str(p.relative_to(n.ROOT)): n.sha(p)
            for p in [Path(__file__), TEMPLATE, BASE]}}


def verify(out):
    provenance = read(out / 'provenance.json')
    if sources() != provenance['sources']:
        raise RuntimeError('Frozen snapshot sources changed')
    if n.sha(out / 'manifest.json') != provenance['manifest_sha256']:
        raise RuntimeError('Snapshot manifest changed')
    return read(out / 'manifest.json')


def initialize(out):
    n.protect_retained_output(out)
    out.mkdir(parents=True, exist_ok=False)
    for folder in ['results', 'work/logs', 'work/cases', 'frozen_sources']:
        (out / folder).mkdir(parents=True)
    baseline = read(BASE)
    conditions = ['fixed_20pct', 'relative_20pct']
    cases, new = [], []
    for length in [15., 30.]:
        for condition in conditions:
            # L=30 is absent from the retained grid. Only convergence scales
            # come from L=31; no native state/profile is reused or interpolated.
            original = next(c for c in baseline['cases'] if c['model'] == 'wave_pinning'
                            and c['length'] == (15. if length == 15. else 31.)
                            and c['condition'] == condition)
            c = copy.deepcopy(original)
            c.pop('durations')
            c.update(case_index=len(cases), length=length,
                     width=3. if condition == 'fixed_20pct' else .2 * length)
            if length == 15.:
                c['preparation_source'] = 'Exact retained corrected basal-release checkpoint'
                if n.sha(n.ROOT / c['checkpoint']) != c['checkpoint_sha256']:
                    raise RuntimeError('Retained checkpoint hash mismatch')
            else:
                c.update(checkpoint=None, checkpoint_sha256=None,
                         establishment_case_index=None, initial_orientation=None,
                         initial_contrast=None, preparation_status='pending', eligible=False,
                         preparation_source='Fresh independent native-rest finite-left-patch/basal-release duration search',
                         scale_source='Frozen L=31 wave-pinning scales used only for convergence diagnostics')
                prep = copy.deepcopy(c)
                prep['case_index'] = len(new)
                c['new_preparation_case_index'] = prep['case_index']
                new.append(prep)
            cases.append(c)
    prep_out = out / 'preparation'
    for folder in ['work/cases', 'results']:
        (prep_out / folder).mkdir(parents=True)
    prep_manifest = read(n.DATA / 'establishment_basal/manifest.json')
    prep_manifest.update(created_at=n.now(), cases=new, amplitude_factors=[1.],
                         scale_rule='L=31 frozen active/inactive scales for new L=30 convergence only; states start independently from native rest.')
    n.write(prep_out / 'manifest.json', prep_manifest)
    n.write(prep_out / 'provenance.json', dict(sources=e.sources(), manifest_sha256=n.sha(prep_out / 'manifest.json')))
    manifest = dict(created_at=n.now(), model='wave_pinning', cases=cases,
                    time=list(range(0, 301, 5)), end_time=300., sample_interval=5.,
                    cue='Continuous right raised-cosine patch at nominal peak 0.1, stop=None, including t=300; intrinsic basal reactions remain active.',
                    fixed_width=3., relative_fraction=.2,
                    initial_state='Exact qualifying corrected basal-release checkpoint; fresh matched preparation at L=30.',
                    convention='At t=0 the field is the pre-cue state; the right cue switches on at this instant. No cue withdrawal or recovery in the snapshot window.')
    n.write(out / 'manifest.json', manifest)
    hashes = sources()
    for rel in hashes:
        dest = out / 'frozen_sources' / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(n.ROOT / rel, dest)
    n.write(out / 'provenance.json', dict(generated_at=n.now(), sources=hashes,
                                         manifest_sha256=n.sha(out / 'manifest.json')))
    publish(out)


def worker(out, slot):
    manifest = verify(out)
    c = copy.deepcopy(manifest['cases'][slot])
    if c['length'] == 30.:
        prep_out = out / 'preparation'
        index = c['new_preparation_case_index']
        e.worker(prep_out, index)
        summary = read(prep_out / 'work/cases' / f'case_{index:04d}' / 'result.json')['summary']
        if summary['status'] != 'polarized' or summary['duration_boundary_unresolved'] or summary['upper_optimum_unresolved']:
            raise RuntimeError('New preparation failed qualification or has an unresolved selection')
        c.update(checkpoint=str((prep_out / summary['selected_checkpoint']).relative_to(n.ROOT)),
                 checkpoint_sha256=summary['selected_checkpoint_sha256'],
                 establishment_summary=summary, preparation_status='polarized', eligible=True)
    checkpoint = n.ROOT / c['checkpoint']
    if n.sha(checkpoint) != c['checkpoint_sha256']:
        raise RuntimeError('Input checkpoint hash mismatch')
    state, prep = n.read_checkpoint(checkpoint, 'wave_pinning')
    c.update(initial_orientation=float(prep.orientation[-1]), initial_contrast=float(prep.contrast[-1]))
    if not (c['initial_orientation'] < -.05 and c['initial_contrast'] >= .05):
        raise RuntimeError('Initial polarity qualification failed')
    cue = e.Patch(c['amplitude'], c['width'], side='right', stop=None,
                  width_mode=c['width_mode'], fraction=.2, background=c['background_stimulus'])
    tr = n.run('wave_pinning', c['length'], state, cue, 300., interval=5., n=c['cells'])
    stimulus = np.stack([cue(t, tr.x, c['length']) for t in tr.time])
    active, inactive = tr.fields['active'], tr.fields['inactive']
    if not np.array_equal(tr.time, np.arange(0, 301, 5)):
        raise RuntimeError('Incorrect sample clocks')
    onset_error = float(np.max(abs(np.vstack([active[0], inactive[0]]) - state)))
    if onset_error > 8 * np.finfo(float).eps * max(1., float(np.max(abs(state)))):
        raise RuntimeError('Snapshot onset differs from exact restart')
    # solve_ivp evaluates dense output at t=0, which can add an ulp. The
    # onset snapshot is defined by the exact checkpoint passed to the solver.
    active[0], inactive[0] = state[0], state[1]
    if not all(np.isfinite(v).all() and v.min() >= -1e-8 for v in [active, inactive]):
        raise RuntimeError('Invalid native fields')
    if not np.array_equal(stimulus, np.broadcast_to(stimulus[0], stimulus.shape)):
        raise RuntimeError('Cue is not continuous across the window')
    # Grid1D uses equal finite-volume cells; sum*dx is the native mass.
    mass = (active + inactive).sum(axis=1) * c['length'] / c['cells']
    mass_error = float(np.max(abs(mass / mass[0] - 1)))
    if mass_error > 1e-7:
        raise RuntimeError(f'Native mass conservation failed: {mass_error}')
    path = out / 'work/cases' / f'case_{slot:04d}'
    np.savez_compressed(path.with_suffix('.npz'), time=tr.time, x=tr.x,
                        active=active, inactive=inactive, stimulus=stimulus,
                        orientation=tr.orientation, total_mass=mass)
    with path.with_suffix('.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['time_s', 'x_um', 'active', 'inactive', 'extra_cue'])
        for i, t in enumerate(tr.time):
            writer.writerows(zip(np.full(len(tr.x), t), tr.x, active[i], inactive[i], stimulus[i]))
    c.update(generated_at=n.now(), data=str(path.with_suffix('.npz').relative_to(out)),
             data_sha256=n.sha(path.with_suffix('.npz')),
             csv=str(path.with_suffix('.csv').relative_to(out)),
             samples=len(tr.time), maximum_relative_mass_error=mass_error,
             solver_onset_roundoff=onset_error,
             final_orientation=float(tr.orientation[-1]), validation_passed=True)
    verify(out)
    n.write(path.with_suffix('.json'), c)
    publish(out)
    print(json.dumps({k: c[k] for k in ['length', 'condition', 'samples', 'final_orientation', 'maximum_relative_mass_error']}), flush=True)


def publish(out):
    with (out / 'results/.publish.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        manifest = verify(out)
        payload = dict(generated_at=n.now(), complete=False, producer='running', cases=[], time=manifest['time'])
        errors = []
        for c in manifest['cases']:
            path = out / 'work/cases' / f"case_{c['case_index']:04d}.json"
            error = path.with_suffix('.error.json')
            if error.exists():
                errors.append(read(error))
            if not path.exists():
                payload['cases'].append(dict(**c, pending=True))
                continue
            meta = read(path)
            with np.load(out / meta['data']) as z:
                meta.update({k: z[k].tolist() for k in ['time', 'x', 'active', 'inactive', 'stimulus', 'orientation', 'total_mass']})
            payload['cases'].append(meta)
        payload['complete'] = all(not c.get('pending') for c in payload['cases']) and not errors
        payload['producer'] = 'finished' if payload['complete'] else 'failed' if errors else 'running'
        payload['errors'] = errors
        compact = json.dumps(payload, separators=(',', ':'), allow_nan=False).replace('</', '<\\/')
        html = TEMPLATE.read_text().replace('__SNAPSHOT_DATA__', compact)
        temporary = out / 'results/index.html.tmp'
        temporary.write_text(html)
        temporary.replace(out / 'results/index.html')
        n.write(out / 'results/status.json', {k: v for k, v in payload.items() if k != 'cases'} | dict(
            completed_cases=sum(not c.get('pending') for c in payload['cases']), requested_cases=4))
        if payload['complete']:
            comparison = {}
            for field in ['active', 'inactive']:
                comparison[field] = float(np.max(abs(np.asarray(payload['cases'][0][field]) - np.asarray(payload['cases'][1][field]))))
            if max(comparison.values()) > 1e-10:
                raise RuntimeError('Equivalent L=15 width cases disagree')
            n.write(out / 'results/validation.json', dict(generated_at=payload['generated_at'], passed=True,
                cases=4, snapshots_per_case=61, exact_checkpoint_start=True,
                continuous_cue_through_300=True, finite_nonnegative_fields=True,
                maximum_relative_mass_error=max(c['maximum_relative_mass_error'] for c in payload['cases']),
                L15_fixed_relative_maximum_absolute_difference=comparison,
                input_checkpoint_hashes={f"L{c['length']:g}_{c['condition']}": c['checkpoint_sha256'] for c in payload['cases']}))


def launch(out):
    initialize(out)
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '
    command += shlex.join([sys.executable, str(Path(__file__).resolve()), 'worker', '--out', str(out)])
    command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(['sbatch', '--parsable', '--qos=small', '--ntasks=1',
        '--cpus-per-task=1', '--mem=3G', '--time=02:00:00', '--array=0-3',
        '--job-name=wp-snapshots', '--chdir=' + str(n.ROOT),
        '--output=' + str(out / 'work/logs/%A_%a.out'),
        '--error=' + str(out / 'work/logs/%A_%a.err'), '--wrap', command], text=True).strip().split(';')[0]
    n.write(out / 'launch.json', dict(generated_at=n.now(), job_id=job, slots=4, qos='small'))
    print(job)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['launch', 'worker', 'publish'])
    parser.add_argument('--out', type=Path, default=DEFAULT)
    parser.add_argument('--slot', type=int)
    args = parser.parse_args()
    out = args.out.resolve()
    if args.action == 'launch':
        launch(out)
    elif args.action == 'publish':
        publish(out)
    else:
        try:
            worker(out, args.slot)
        except Exception:
            n.write(out / 'work/cases' / f'case_{args.slot:04d}.error.json', dict(
                generated_at=n.now(), slot=args.slot, traceback=traceback.format_exc()))
            publish(out)
            raise
