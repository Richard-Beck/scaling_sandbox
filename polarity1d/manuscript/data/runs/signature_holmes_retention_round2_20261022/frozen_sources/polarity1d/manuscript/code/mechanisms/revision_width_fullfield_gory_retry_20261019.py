"""Four Goryachev retries after saved-output roundoff assertion failure; freeze/preflight never submit.

Scientific integration is the unchanged native.run path. Dense output retains
both native species throughout cue and release; legacy 2.5-unit gate samples
remain a strict subset of the recovery output. No threshold search is run.
"""
from pathlib import Path
import argparse
import copy
import csv
import json
import shutil
import sys
import traceback
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code/core'))
import native as n
import polarity_establishment as e
import settling

SOURCE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
DEFAULT = ROOT / 'data/runs/revision_width_fullfield_gory_retry_20261019'
DESIGN = [('goryachev', 9., 74.), ('goryachev', 25., 815.)]
CUE_INTERVAL = .1
OFF_INTERVAL = .25
OFF_DURATION = 10000.


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n.write(path, data)


def csvwrite(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp.csv')
    n.csvwrite(temp, rows)
    temp.replace(path)


def cue(case):
    return e.Patch(case['response_amplitude'], case['width'], side='right',
                   width_mode=case['width_mode'], fraction=.2,
                   background=case['background_stimulus'])


def regional_inputs(case, state, x):
    u = state[0]
    half = float(u.min() + .5 * np.ptp(u))
    above = u > half
    peak = int(np.argmax(u))
    left = peak
    right = peak
    while left and above[left - 1]:
        left -= 1
    while right < len(u) - 1 and above[right + 1]:
        right += 1
    core = np.zeros(len(u), bool)
    core[left:right + 1] = True
    width = .2 * case['length'] if case['width_mode'] == 'relative' else case['width']
    footprint = (x >= case['length'] - width) & (x <= case['length'])
    signal = cue(case)(0., x, case['length'])
    return dict(x=x, initial_state=state, incumbent_core_mask=core,
                cue_footprint_mask=footprint, signal=signal,
                initial_halfcontrast_threshold=np.array(half))


def initialize(out):
    if out.exists():
        raise ValueError('Use a fresh campaign destination')
    original = read(SOURCE / 'manifest.json')
    if not read(SOURCE / 'results/status.json')['complete']:
        raise RuntimeError('Original campaign is incomplete')
    with (SOURCE / 'results/endpoints.csv').open() as f:
        endpoints = list(csv.DictReader(f))
    cases = []
    for model, length, duration in DESIGN:
        selected = sorted([c for c in original['cases'] if c['model'] == model
                           and c['length'] == length and c['condition'] in
                           ('fixed_physical', 'fixed_relative')], key=lambda c: c['width'])
        assert len(selected) == 2
        assert len({c['checkpoint_sha256'] for c in selected}) == 1
        for old in selected:
            c = {k: copy.deepcopy(v) for k, v in old.items()
                 if k not in ('base_durations', 'durations', 'reused_positive_trials')}
            assert c['eligible'] and c['cells'] == 201
            assert c['response_amplitude'] == (.1 if model == 'wave_pinning' else 6.4)
            gate_rows = [r for r in endpoints if int(r['case_index']) == c['case_index']]
            assert len(gate_rows) == 2
            assert all(r['status'] == 'success' and r['selected_field_tail_pass'] == 'True'
                       and r['nonmonotonic_success_pattern'] == 'False' for r in gate_rows)
            c.update(slot=len(cases), source_case_index=old['case_index'],
                     pair_id=f'{model}_L{length:g}', cue_duration=duration,
                     width_label='narrow' if old == selected[0] else 'broad',
                     original_thresholds={r['definition']: dict(
                         lower_failing_duration=float(r['lower_failing_duration']),
                         passing_duration=float(r['endpoint']),
                         field_tail_pass=True) for r in gate_rows})
            cases.append(c)
        # Selected pulse lies strictly between the previously disjoint passing brackets.
        low, high = sorted([c['original_thresholds']['full_polarity'] for c in cases[-2:]],
                           key=lambda r: r['passing_duration'])
        assert low['passing_duration'] < duration < high['lower_failing_duration']
    for folder in ('inputs/basals', 'inputs/regions', 'workers', 'logs', 'results', 'frozen_sources'):
        (out / folder).mkdir(parents=True, exist_ok=True)
    frozen = out / 'frozen_sources'
    source_hashes = n.sources()
    for path in (Path(__file__), ROOT / 'code/three_test_reversal.py'):
        source_hashes[str(path.relative_to(ROOT))] = n.sha(path)
    for rel, digest in source_hashes.items():
        destination = frozen / rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, destination)
        assert n.sha(destination) == digest
    checks = []
    for c in cases:
        source = ROOT / c['checkpoint']
        assert n.sha(source) == c['checkpoint_sha256']
        dest = out / 'inputs/basals' / f'{c["pair_id"]}.npz'
        if not dest.exists():
            shutil.copy2(source, dest)
        assert n.sha(dest) == c['checkpoint_sha256']
        c['frozen_checkpoint'] = str(dest.relative_to(out))
        state, tr = n.read_checkpoint(dest, c['model'])
        assert state.shape == (2, 201) and np.isfinite(state).all() and state.min() >= 0
        np.testing.assert_array_equal(tr.fields['active'][-1], state[0])
        np.testing.assert_array_equal(tr.fields['inactive'][-1], state[1])
        np.testing.assert_allclose(tr.x, (np.arange(201) + .5) * c['length'] / 201,
                                   rtol=0, atol=4e-15)
        z = regional_inputs(c, state, tr.x)
        path = out / 'inputs/regions' / f'{c["slot"]:02d}.npz'
        np.savez_compressed(path, **z)
        c['regional_input'] = str(path.relative_to(out))
        c['regional_input_sha256'] = n.sha(path)
        signal = z['signal']
        model = n.model_for(c['model'])
        response = model.reaction(0, tr.x, state, signal) - model.reaction(
            0, tr.x, state, np.zeros_like(signal))
        expected = signal * state[1]
        error = float(np.max(abs(response - np.vstack((expected, -expected)))))
        assert error < 1e-12
        assert not np.any(z['incumbent_core_mask'] & z['cue_footprint_mask'])
        dx = c['length'] / c['cells']
        c.update(sampled_peak=float(signal.max()), sampled_spatial_integral=float(signal.sum() * dx),
                 cue_roi_cells=int(z['cue_footprint_mask'].sum()),
                 initial_core_cells=int(z['incumbent_core_mask'].sum()))
        checks.append(dict(slot=c['slot'], native_channel='stimulus * inactive, added to active and subtracted from inactive',
                           pointwise_native_reaction_difference_error=error,
                           checkpoint_exact=True, core_cue_overlap_cells=0,
                           sampled_peak=c['sampled_peak'], sampled_spatial_integral=c['sampled_spatial_integral']))
    frames = sum(round(c['cue_duration'] / CUE_INTERVAL) + 1 + round(OFF_DURATION / OFF_INTERVAL) + 1
                 for c in cases)
    manifest = dict(created_at=n.now(), experiment='four_goryachev_fullfield_validation_retry',
                    retry_of='data/runs/revision_width_fullfield_20261019',
                    retry_reason='Original exact-input integration aborted only because saved t=0 BDF interpolation was required to be bitwise identical; scientific solver and inputs are unchanged.',
                    original_campaign=str(SOURCE.relative_to(ROOT)), cases=cases,
                    source_manifest_sha256=n.sha(SOURCE / 'manifest.json'),
                    source_endpoints_sha256=n.sha(SOURCE / 'results/endpoints.csv'),
                    cue_sample_interval=CUE_INTERVAL, recovery_sample_interval=OFF_INTERVAL,
                    recovery_duration=OFF_DURATION, legacy_gate_interval=2.5, legacy_gate_tail=500.,
                    integration='Unmodified native.run; BDF, rtol=2e-6, atol=1e-9, max_step=10 with original short-phase cap; two independent cue/off calls with exact withdrawal restart.',
                    frozen_regions='Initial connected above-halfcontrast region containing the incumbent maximum; fixed actual cue footprint. Native cell-centre quadrature dx=L/201.',
                    primary_readouts='Native incumbent-core amount/peak and cue-footprint challenger amount/peak versus time; fractional regional depletion is an activity change, not identified molecular damage.',
                    limitations=['Only eight selected protocols; no global monotonicity assumption.',
                                 'Fixed requested peak, different spatial integral; no pure width versus dose attribution.',
                                 'Native channel and clocks are not calibrated optical actuators or reporters.',
                                 'N=201 inherited discretization; no spatial convergence claim.',
                                 'Scalar polarity success and native field-tail pass are separate from regional geography.'],
                    resources=dict(cpus_per_task=1, memory='4G', walltime='02:00:00', qos='small',
                                   tasks=4, array_throttle=None, saved_frames=frames,
                                   active_inactive_uncompressed_bytes=frames * 2 * 201 * 8,
                                   estimate_basis='Array-size estimate only; no benchmark/integration run. About 2 h reservation per independent task is conservative for these native two-species trajectories.'))
    write(out / 'manifest.json', manifest)
    input_hashes = {str(p.relative_to(out)): n.sha(p) for p in sorted((out / 'inputs').rglob('*.npz'))}
    write(out / 'provenance.json', dict(created_at=n.now(), manifest_sha256=n.sha(out / 'manifest.json'),
         frozen_sources=source_hashes, frozen_inputs=input_hashes))
    write(out / 'results/preflight.json', dict(generated_at=n.now(), passed=True, integration_executed=False,
         jobs_submitted=False, cases=checks, exact_basal_pairs=2, total_conditions=4,
         planned_fullfield_frames=frames, legacy_grid_is_subset=True))
    write(out / 'results/latest_status.json', dict(generated_at=n.now(), producer='frozen_not_submitted',
         complete=False, completed=0, total=4, pending=4, errors=0))
    rows = []
    for c in cases:
        r = dict(slot=c['slot'], pair_id=c['pair_id'], model=c['model'], length=c['length'],
                 width=c['width'], width_mode=c['width_mode'], label=c['width_label'],
                 requested_peak=c['response_amplitude'], sampled_peak=c['sampled_peak'],
                 sampled_spatial_integral=c['sampled_spatial_integral'], pulse_duration=c['cue_duration'],
                 off_duration=OFF_DURATION, source_case_index=c['source_case_index'],
                 checkpoint_sha256=c['checkpoint_sha256'], cue_roi_cells=c['cue_roi_cells'],
                 initial_core_cells=c['initial_core_cells'])
        r.update(c['original_thresholds']['full_polarity'])
        rows.append(r)
    csvwrite(out / 'results/condition_manifest.csv', rows)
    frozen_worker = out / 'frozen_sources/code/mechanisms' / Path(__file__).name
    shell = '#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n'
    shell += f'exec /home/4473331/.conda/envs/scaling_sandbox/bin/python "{frozen_worker}" worker --out "{out}" --slot "$SLURM_ARRAY_TASK_ID"\n'
    (out / 'worker.sh').write_text(shell)
    (out / 'submission_command.txt').write_text(
        f'sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=02:00:00 --array=0-3 '
        f'--job-name=F2-width-fullfield --output={out}/logs/%A_%a.out --error={out}/logs/%A_%a.err {out}/worker.sh\n')
    (out / 'README.md').write_text('# Four Goryachev full-field validation retries\n\nFrozen and not submitted. '
        'See manifest.json, results/condition_manifest.csv and results/preflight.json. '
        'The worker uses the frozen scientific source tree and exact copied basals. '
        'submission_command.txt is a proposal for root review, not a launch record.\n\n'
        'Both native fields are saved every 0.1 cue-on unit and every 0.25 cue-free unit through 10,000; '
        'the original 2.5-unit recovery grid is a subset. Original scalar gates and full-field tail diagnostics '
        'are computed on that original grid, separately from dense regional curves. Fixed initial core and actual '
        'cue footprint are in inputs/regions. No interpolation of trajectories or spatial refinement is performed.\n')
    print(json.dumps(dict(out=str(out), manifest_sha256=n.sha(out / 'manifest.json'),
                         preflight_sha256=n.sha(out / 'results/preflight.json'),
                         conditions=len(cases), frames=frames), indent=2))


def verify(out):
    p = read(out / 'provenance.json')
    assert n.sha(out / 'manifest.json') == p['manifest_sha256']
    for rel, digest in p['frozen_sources'].items():
        assert n.sha(out / 'frozen_sources' / rel) == digest, rel
    for rel, digest in p['frozen_inputs'].items():
        assert n.sha(out / rel) == digest, rel
    return read(out / 'manifest.json')


def save_phase(path, tr):
    arrays = dict(time=tr.time, x=tr.x, orientation=tr.orientation, contrast=tr.contrast,
                  finalstate=tr.state)
    arrays.update({'field_' + k: v for k, v in tr.fields.items()})
    temp = path.with_suffix('.tmp.npz')
    np.savez_compressed(temp, **arrays)
    temp.replace(path)


def legacy_gates(tr, case):
    indices = np.arange(0, len(tr.time), round(2.5 / OFF_INTERVAL))
    time = tr.time[indices]
    np.testing.assert_array_equal(time, np.arange(0, OFF_DURATION + 2.5, 2.5))
    tail = time >= OFF_DURATION - 500.
    o = tr.orientation[indices][tail]
    contrast = tr.contrast[indices][tail]
    orientation = bool(o.min() > .05 and contrast.min() >= .05 and np.ptp(o) < .005)
    error = float(np.max(abs(o / abs(case['initial_orientation']) - 1.)))
    metrics = {}
    maximum = 0.
    for name, y in tr.fields.items():
        a = y[indices][tail]
        scale = case['fixed_scales'][name]
        change = float(np.max(np.ptp(a, axis=0)) / scale)
        drift = float(np.max(abs(a[-1] - a[-2])) / 2.5 * 500. / scale)
        metrics[name] = dict(range_scaled=change, projected_drift_500=drift)
        maximum = max(maximum, change, drift)
    return dict(orientation_success=orientation, full_polarity_success=bool(orientation and error <= .01),
                final_orientation=float(o[-1]), minimum_contrast_last_500=float(contrast.min()),
                orientation_range_last_500=float(np.ptp(o)), maximum_relative_polarity_error_last_500=error,
                field_tail_pass=bool(maximum <= .001), field_tail_max_scaled=maximum,
                field_tail_metrics=metrics, gate_sample_interval=2.5, gate_window=500.)


def regional_rows(tr, c, z, phase):
    dx = c['length'] / c['cells']
    u, v = tr.fields['active'], tr.fields['inactive']
    core, roi = z['incumbent_core_mask'], z['cue_footprint_mask']
    initial_amount = float(z['initial_state'][0, core].sum() * dx)
    initial_total = float(z['initial_state'].sum() * dx)
    core_amount = u[:, core].sum(axis=1) * dx
    challenger_amount = u[:, roi].sum(axis=1) * dx
    amount = (u + v).sum(axis=1) * dx
    rows = []
    for i, time in enumerate(tr.time):
        rows.append(dict(slot=c['slot'], pair_id=c['pair_id'], width_label=c['width_label'],
                         phase=phase, time_phase_native=float(time),
                         time_since_cue_onset_native=float(time + (c['cue_duration'] if phase == 'off' else 0.)),
                         incumbent_core_active_amount=float(core_amount[i]),
                         incumbent_core_active_fraction=float(core_amount[i] / initial_amount),
                         incumbent_core_regional_depletion=float(1. - core_amount[i] / initial_amount),
                         incumbent_core_active_peak=float(u[i, core].max()),
                         cue_footprint_active_amount=float(challenger_amount[i]),
                         cue_footprint_active_peak=float(u[i, roi].max()),
                         cue_footprint_inactive_mean=float(v[i, roi].mean()),
                         total_active_amount=float(u[i].sum() * dx),
                         total_active_plus_inactive_amount=float(amount[i]),
                         total_amount_relative_error=float(abs(amount[i] / initial_total - 1.)),
                         scalar_orientation=float(tr.orientation[i]), contrast=float(tr.contrast[i])))
    return rows


def saved_initial_roundoff(actual, expected):
    # Exact checkpoint/restart arrays are passed to native.run. SciPy's saved
    # t=0 dense-output interpolation can change tiny tails at roundoff level.
    scale = max(1., float(np.max(np.abs(expected))))
    error = float(np.max(np.abs(actual - expected)))
    bound = 32 * np.finfo(float).eps * scale
    assert error <= bound, (error, bound)
    return dict(maximum_absolute_error=error, scale=scale, scaled_error=error / scale,
                permitted_scaled_roundoff=32 * np.finfo(float).eps,
                interpretation='Saved initial output roundoff; exact input identity is checked separately.')


def worker(out, slot):
    m = verify(out)
    c = m['cases'][slot]
    folder = out / 'workers' / f'{slot:02d}'
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / 'result.json').exists():
        raise RuntimeError('Worker already completed; refusing to overwrite native evidence')
    write(folder / 'status.json', dict(generated_at=n.now(), state='running', slot=slot))
    state, _ = n.read_checkpoint(out / c['frozen_checkpoint'], c['model'])
    with np.load(out / c['regional_input']) as zz:
        z = {k: zz[k].copy() for k in zz.files}
    np.testing.assert_array_equal(state, z['initial_state'])
    pulse = n.run(c['model'], c['length'], state, cue(c), c['cue_duration'],
                  interval=CUE_INTERVAL, n=c['cells'])
    cue_roundoff = saved_initial_roundoff(np.vstack([pulse.fields[k][0] for k in ('active', 'inactive')]), state)
    save_phase(folder / 'cue.npz', pulse)
    pulse_rows = regional_rows(pulse, c, z, 'cue_on')
    restart = pulse.state.copy()
    del pulse
    write(folder / 'status.json', dict(generated_at=n.now(), state='recovery', slot=slot,
                                      cue_sha256=n.sha(folder / 'cue.npz')))
    recovery = n.run(c['model'], c['length'], restart, e.withdrawal(c), OFF_DURATION,
                     interval=OFF_INTERVAL, n=c['cells'])
    restart_roundoff = saved_initial_roundoff(np.vstack([recovery.fields[k][0] for k in ('active', 'inactive')]), restart)
    save_phase(folder / 'recovery.npz', recovery)
    rows = pulse_rows + regional_rows(recovery, c, z, 'off')
    csvwrite(folder / 'regional_readouts.csv', rows)
    diagnostics = legacy_gates(recovery, c)
    finite = all(np.isfinite(y).all() for y in recovery.fields.values())
    minimum = min(float(y.min()) for y in recovery.fields.values())
    result = dict(generated_at=n.now(), slot=slot, state='complete', case=c,
                  manifest_sha256=n.sha(out / 'manifest.json'),
                  legacy_diagnostics=diagnostics,
                  maximum_total_amount_relative_error=max(r['total_amount_relative_error'] for r in rows),
                  finite=finite, minimum_saved_recovery_field=minimum,
                  cue_input_state_exact=True, withdrawal_restart_input_exact=True,
                  saved_initial_output_roundoff=dict(cue=cue_roundoff, recovery=restart_roundoff),
                  cue_sha256=n.sha(folder / 'cue.npz'), recovery_sha256=n.sha(folder / 'recovery.npz'),
                  regional_readouts_sha256=n.sha(folder / 'regional_readouts.csv'),
                  saved_readouts=len(rows))
    assert finite and minimum >= -1e-7
    write(folder / 'result.json', result)
    write(folder / 'status.json', dict(generated_at=n.now(), state='complete', slot=slot))


def collect(out):
    m = verify(out)
    rows = []
    for c in m['cases']:
        folder = out / 'workers' / f'{c["slot"]:02d}'
        p = folder / 'result.json'
        r = dict(slot=c['slot'], pair_id=c['pair_id'], width_label=c['width_label'], state='not_started')
        if p.exists():
            result = read(p)
            assert result['case'] == c
            for key, name in (('cue_sha256', 'cue.npz'), ('recovery_sha256', 'recovery.npz'),
                              ('regional_readouts_sha256', 'regional_readouts.csv')):
                assert n.sha(folder / name) == result[key]
            r.update(state='complete', **{k: result['legacy_diagnostics'][k] for k in
                ('orientation_success', 'full_polarity_success', 'field_tail_pass')})
        elif (folder / 'error.json').exists():
            r['state'] = 'error'
        elif (folder / 'status.json').exists():
            r['state'] = read(folder / 'status.json')['state']
        rows.append(r)
    write(out / 'results/condition_status.json', rows)
    done = sum(r['state'] == 'complete' for r in rows)
    errors = sum(r['state'] == 'error' for r in rows)
    write(out / 'results/latest_status.json', dict(generated_at=n.now(), complete=done == 4,
          producer='finished' if done == 4 else 'execution_errors' if errors else 'partial_or_not_started',
          completed=done, pending=4 - done - errors, errors=errors, total=4))
    print(json.dumps(read(out / 'results/latest_status.json'), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('init', 'preflight', 'worker', 'collect'))
    parser.add_argument('--out', type=Path, default=DEFAULT)
    parser.add_argument('--slot', type=int)
    args = parser.parse_args()
    out = args.out.resolve()
    if args.action == 'init':
        initialize(out)
    elif args.action == 'preflight':
        m = verify(out)
        print(json.dumps(dict(passed=True, conditions=len(m['cases']), integration_executed=False,
                             manifest_sha256=n.sha(out / 'manifest.json')), indent=2))
    elif args.action == 'collect':
        collect(out)
    else:
        if args.slot is None or not 0 <= args.slot < 4:
            raise ValueError('slot must be 0..3')
        try:
            worker(out, args.slot)
        except Exception:
            write(out / 'workers' / f'{args.slot:02d}' / 'error.json',
                  dict(generated_at=n.now(), state='error', traceback=traceback.format_exc()))
            raise


if __name__ == '__main__':
    main()
