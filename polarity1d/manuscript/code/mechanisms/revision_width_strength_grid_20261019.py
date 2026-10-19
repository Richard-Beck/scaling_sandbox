"""Bounded 72-protocol two-model strength x physical-width x pulse experiment.

No automatic launch. The native scientific solver is unchanged; explicit saved
times add early-release resolution without an extra integration restart.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import copy
import csv
import json
import shutil
import sys
import traceback
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code/mechanisms'))
import revision_width_fullfield_gory_retry_20261019 as h
n, e = h.n, h.e
from polarity1d.contract import output_times

DEFAULT = ROOT / 'data/runs/revision_width_strength_grid_20261019'
ORIGINAL = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
FIRST = ROOT / 'data/runs/revision_width_fullfield_20261019'
RETRY = ROOT / 'data/runs/revision_width_fullfield_gory_retry_20261019'
COHORTS = [('wave_pinning', 9.), ('wave_pinning', 19.), ('goryachev', 9.), ('goryachev', 25.)]
DURATIONS = dict(wave_pinning=[16., 60., 240.], goryachev=[74., 400., 1600.])
WIDTHS = [3., 5.]
MULTIPLIERS = [.5, 1., 2.]
OFF = 10000.


def off_times():
    return np.unique(np.r_[output_times(OFF, 2.5), output_times(100., .25)])


def first_results():
    available = {}
    for root in (FIRST, RETRY):
        for p in sorted((root / 'workers').glob('*/result.json')):
            r = h.read(p)
            c = r['case']
            key = (c['model'], c['length'], c['width'], c['response_amplitude'], c['cue_duration'])
            available[key] = dict(result_path=str(p.relative_to(ROOT)), result_sha256=n.sha(p),
                cue_path=str((p.parent / 'cue.npz').relative_to(ROOT)), cue_sha256=r['cue_sha256'],
                recovery_path=str((p.parent / 'recovery.npz').relative_to(ROOT)), recovery_sha256=r['recovery_sha256'])
    assert len(available) == 8, 'Finish the eight native-field protocols before freezing the next stage'
    return available


def initialize(out):
    if out.exists():
        raise ValueError('Use a fresh output root')
    original = h.read(ORIGINAL / 'manifest.json')
    available = first_results()
    for directory in ('inputs/basals', 'inputs/regions', 'frozen_sources', 'workers', 'results', 'logs'):
        (out / directory).mkdir(parents=True, exist_ok=True)
    sources = n.sources()
    for p in (Path(__file__), Path(h.__file__)):
        sources[str(p.relative_to(ROOT))] = n.sha(p)
    for rel, digest in sources.items():
        destination = out / 'frozen_sources' / rel
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, destination)
        assert n.sha(destination) == digest
    evidence = ROOT / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/width_fullfield_experiment_20261019/supporting/data/actual_fullfield_summary.csv'
    shutil.copy2(evidence, out / 'inputs/first_eight_summary.csv')
    cases, jobs, checks = [], [], []
    for model, length in COHORTS:
        base = next(c for c in original['cases'] if c['model'] == model and c['length'] == length
                    and c['condition'] == 'fixed_physical')
        assert base['eligible'] and base['cells'] == 201
        src = ROOT / base['checkpoint']
        assert n.sha(src) == base['checkpoint_sha256']
        dest = out / 'inputs/basals' / f'{model}_L{length:g}.npz'
        shutil.copy2(src, dest)
        state, tr = n.read_checkpoint(dest, model)
        for width in WIDTHS:
            for multiplier in MULTIPLIERS:
                for duration in DURATIONS[model]:
                    c = {k: copy.deepcopy(v) for k, v in base.items()
                         if k not in ('base_durations', 'durations', 'reused_positive_trials')}
                    c.update(slot=len(cases), case_index=len(cases), source_case_index=base['case_index'],
                             condition='physical_width_strength_duration_grid', width_mode='physical',
                             width=width, native_peak_multiplier=multiplier,
                             baseline_requested_peak=base['response_amplitude'],
                             response_amplitude=base['response_amplitude'] * multiplier,
                             response_factor=base['response_factor'] * multiplier,
                             cue_duration=duration, pair_id=f'{model}_L{length:g}',
                             width_label='narrow' if width == WIDTHS[0] else 'broad',
                             frozen_checkpoint=str(dest.relative_to(out)))
                    z = h.regional_inputs(c, state, tr.x)
                    assert not np.any(z['incumbent_core_mask'] & z['cue_footprint_mask'])
                    region = out / 'inputs/regions' / f'{c["slot"]:03d}.npz'
                    np.savez_compressed(region, **z)
                    c['regional_input'] = str(region.relative_to(out))
                    c['regional_input_sha256'] = n.sha(region)
                    c['sampled_peak'] = float(z['signal'].max())
                    c['sampled_spatial_integral'] = float(z['signal'].sum() * length / 201)
                    key = (model, length, width, c['response_amplitude'], duration)
                    c['exact_reuse'] = available.get(key)
                    if c['exact_reuse'] is None:
                        jobs.append(dict(task_index=len(jobs), case_slot=c['slot']))
                    expected = z['signal'] * state[1]
                    m = n.model_for(model)
                    delta = m.reaction(0, tr.x, state, z['signal']) - m.reaction(0, tr.x, state, np.zeros(201))
                    error = float(np.max(abs(delta - np.vstack((expected, -expected)))))
                    assert error < 1e-12
                    checks.append(dict(slot=c['slot'], channel_max_absolute_error=error,
                                       cue_core_overlap_cells=0, exact_checkpoint=True,
                                       new_integration=c['exact_reuse'] is None))
                    cases.append(c)
    assert len(cases) == 72 and len(jobs) == 69
    canonical = output_times(OFF, 2.5)
    np.testing.assert_array_equal(off_times()[np.searchsorted(off_times(), canonical)], canonical)
    frames = sum(len(output_times(c['cue_duration'], .1)) + len(off_times())
                 for c in cases if c['exact_reuse'] is None)
    manifest = dict(created_at=n.now(), experiment='bounded_two_model_strength_width_duration_stage',
        cases=cases, jobs=jobs, models=['wave_pinning', 'goryachev'], models_not_in_this_stage=['debelly', 'holmes_model3', 'spring'],
        source_manifest=str((ORIGINAL / 'manifest.json').relative_to(ROOT)), source_manifest_sha256=n.sha(ORIGINAL / 'manifest.json'),
        independent_controls=dict(physical_widths=WIDTHS, native_peak_multipliers=MULTIPLIERS,
                                  pulse_durations_by_model=DURATIONS, size_cohorts=COHORTS),
        motivation='The eight full-field trials show opposite width-fate effects at GoryachevL9 versusL25, while a broader cue improves WP outcomes. This stage varies strength, width and pulse independently within each model/size; pulses are identical across the two sizes within each model.',
        first_eight_summary_sha256=n.sha(out / 'inputs/first_eight_summary.csv'),
        integration='Original native BDF scientific closure, rtol2e-6, atol1e-9, max_step10 with original phase cap. Two phases only; exact withdrawal restart. Only t_eval is changed for nonuniform saved recovery times.',
        sampling=dict(cue_interval=.1, recovery_base_interval=2.5, recovery_early_interval=.25,
                      recovery_early_end=100., recovery_end=OFF, canonical_gates_interval=2.5, tail_window=500.),
        primary_measurements='Fixed initial incumbent-core native activity and actual cue-support challenger activity versus time, cue-on sampled depletion latencies, post-withdrawal regional fate. Scalar-strength and full-field-tail diagnostics remain separate.',
        limits=['Two-model stage, not five-model functional coverage.', 'Native actuator is stimulus times inactive; no optical calibration.',
                'Fixed peak does not fix integrated dose as width changes.', 'InheritedN201; no spatial convergence claim.',
                'Sampled regional activity loss is not identified molecular damage or unique causal attribution.',
                'Three pulse durations sample responses; no monotonicity assumption or minimum-duration estimate.'],
        resources=dict(new_jobs=len(jobs), exact_reused_protocols=3, cpus_per_task=1, memory='4G', qos='small',
                       walltime='02:00:00', array_throttle=None, new_saved_frames=frames,
                       new_two_field_uncompressed_bytes=frames * 2 * 201 * 8))
    h.write(out / 'manifest.json', manifest)
    inputs = {str(p.relative_to(out)): n.sha(p) for p in sorted((out / 'inputs').rglob('*')) if p.is_file()}
    h.write(out / 'provenance.json', dict(created_at=n.now(), manifest_sha256=n.sha(out / 'manifest.json'), frozen_sources=sources, frozen_inputs=inputs))
    h.write(out / 'results/preflight.json', dict(generated_at=n.now(), passed=True, conditions=72, new_jobs=69,
         exact_reuses=3, integration_executed=False, native_channel_checks=checks, canonical_gate_times_exact_subset=True))
    h.write(out / 'results/latest_status.json', dict(generated_at=n.now(), producer='frozen_not_submitted',
         complete=False, total=72, completed=3, exact_reuses=3, new_completed=0, pending=69, errors=0))
    rows = [dict(slot=c['slot'], model=c['model'], length=c['length'], width=c['width'],
                 peak_multiplier=c['native_peak_multiplier'], requested_peak=c['response_amplitude'],
                 sampled_peak=c['sampled_peak'], sampled_integral=c['sampled_spatial_integral'],
                 pulse_duration=c['cue_duration'], off_duration=OFF,
                 exact_reuse=bool(c['exact_reuse']), checkpoint_sha256=c['checkpoint_sha256']) for c in cases]
    h.csvwrite(out / 'results/condition_manifest.csv', rows)
    script = out / 'frozen_sources/code/mechanisms' / Path(__file__).name
    (out / 'worker.sh').write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n'
        f'exec /home/4473331/.conda/envs/scaling_sandbox/bin/python "{script}" worker --out "{out}" --task "$SLURM_ARRAY_TASK_ID"\n')
    (out / 'submission_command.txt').write_text(f'sbatch --parsable --qos=small --ntasks=1 --cpus-per-task=1 --mem=4G --time=02:00:00 --array=0-68 --job-name=F2-width-strength --output={out}/logs/%A_%a.out --error={out}/logs/%A_%a.err {out}/worker.sh\n')
    (out / 'README.md').write_text('# Bounded strength × width × pulse grid\n\n72 two-model protocols;69 new native integrations and3 exact full-field reuses. Frozen, not submitted. '
        'The eight first-stage histories motivated the grid; Gory1600 extends beyond the established large-cell baseline threshold. '
        'Widths3/5 and peak multipliers0.5/1/2 are independently crossed with the same three pulse durations at both sizes within each model. '
        'This is not five-model functional coverage.\n\nSee manifest.json, results/condition_manifest.csv, results/preflight.json and submission_command.txt. '
        'The original scientific BDF solver accepts explicit output times; off samples preserve every canonical2.5-unit gate sample and add0.25-unit samples through100, without extra restarts. '
        'Exact reused native archives remain at their original paths and hashes.\n')
    print(json.dumps(dict(out=str(out), manifest_sha256=n.sha(out / 'manifest.json'), conditions=72,
                         new_jobs=69, raw_two_field_bytes=frames * 2 * 201 * 8), indent=2))


def verify(out):
    m = h.verify(out)
    package = out.parents[2]
    for c in m['cases']:
        r = c['exact_reuse']
        if r:
            for path, digest in (('result_path', 'result_sha256'), ('cue_path', 'cue_sha256'), ('recovery_path', 'recovery_sha256')):
                assert n.sha(package / r[path]) == r[digest]
    return m


def native_sampled_run(case, state, stimulus, times):
    # Same native._sim scientific kernel and flags; explicit t_eval instead of
    # uniform output_interval. No intermediate restarts or changed max_step.
    state = np.asarray(state, dtype=float)
    model = n.model_for(case['model'])
    if model.nonnegative:
        if state.min() < -1e-8:
            raise ValueError('Invalid native nonnegative restart')
        state = np.maximum(state, 0.)
    tr = n.simulate(model, initial_state=state, length=float(case['length']), cells=201,
                    protocol=stimulus, times=times, end_time=float(times[-1]),
                    max_step=10., backend='bdf')
    fields = {name: tr.fields[name] for name in model.state_names}
    active = np.maximum(n.rd._active(tr), 0.)
    return SimpleNamespace(time=tr.time, x=tr.x, fields=fields, active=active,
        orientation=n.rd._orientation(tr), contrast=np.ptp(active, axis=1),
        state=np.vstack([fields[name][-1] for name in model.state_names]))


def gates(tr, case):
    target = output_times(OFF, 2.5)
    index = np.searchsorted(tr.time, target)
    np.testing.assert_array_equal(tr.time[index], target)
    tail = target >= OFF - 500.
    o, contrast = tr.orientation[index][tail], tr.contrast[index][tail]
    orientation = bool(o.min() > .05 and contrast.min() >= .05 and np.ptp(o) < .005)
    error = float(np.max(abs(o / abs(case['initial_orientation']) - 1.)))
    maximum, details = 0., {}
    for name, values in tr.fields.items():
        a = values[index][tail]
        scale = case['fixed_scales'][name]
        change = float(np.max(np.ptp(a, axis=0)) / scale)
        drift = float(np.max(abs(a[-1] - a[-2])) / 2.5 * 500. / scale)
        maximum = max(maximum, change, drift)
        details[name] = dict(range_scaled=change, projected_drift_500=drift)
    return dict(orientation_success=orientation, full_polarity_success=bool(orientation and error <= .01),
                final_orientation=float(o[-1]), maximum_relative_polarity_error_last_500=error,
                field_tail_pass=bool(maximum <= .001), field_tail_max_scaled=maximum,
                field_tail_metrics=details)


def worker(out, task):
    m = verify(out)
    c = m['cases'][m['jobs'][task]['case_slot']]
    folder = out / 'workers' / f'{c["slot"]:03d}'
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / 'result.json').exists():
        raise RuntimeError('Refusing to overwrite completed native evidence')
    h.write(folder / 'status.json', dict(generated_at=n.now(), state='cue_on', task=task))
    state, _ = n.read_checkpoint(out / c['frozen_checkpoint'], c['model'])
    with np.load(out / c['regional_input']) as zz:
        z = {k: zz[k].copy() for k in zz.files}
    np.testing.assert_array_equal(state, z['initial_state'])
    pulse = native_sampled_run(c, state, h.cue(c), output_times(c['cue_duration'], .1))
    initial_roundoff = h.saved_initial_roundoff(np.vstack([pulse.fields[k][0] for k in ('active', 'inactive')]), state)
    h.save_phase(folder / 'cue.npz', pulse)
    rows = h.regional_rows(pulse, c, z, 'cue_on')
    restart = pulse.state.copy()
    del pulse
    h.write(folder / 'status.json', dict(generated_at=n.now(), state='recovery', task=task))
    recovery = native_sampled_run(c, restart, e.withdrawal(c), off_times())
    restart_roundoff = h.saved_initial_roundoff(np.vstack([recovery.fields[k][0] for k in ('active', 'inactive')]), restart)
    h.save_phase(folder / 'recovery.npz', recovery)
    rows += h.regional_rows(recovery, c, z, 'off')
    h.csvwrite(folder / 'regional_readouts.csv', rows)
    result = dict(generated_at=n.now(), state='complete', case=c, task_index=task,
                  manifest_sha256=n.sha(out / 'manifest.json'), legacy_diagnostics=gates(recovery, c),
                  maximum_total_amount_relative_error=max(r['total_amount_relative_error'] for r in rows),
                  exact_input_identity=True, saved_initial_output_roundoff=dict(cue=initial_roundoff, recovery=restart_roundoff),
                  cue_sha256=n.sha(folder / 'cue.npz'), recovery_sha256=n.sha(folder / 'recovery.npz'),
                  regional_readouts_sha256=n.sha(folder / 'regional_readouts.csv'))
    h.write(folder / 'result.json', result)
    h.write(folder / 'status.json', dict(generated_at=n.now(), state='complete', task=task))


def collect(out):
    m = verify(out)
    package = out.parents[2]
    rows = []
    for c in m['cases']:
        folder = out / 'workers' / f'{c["slot"]:03d}'
        reuse = c['exact_reuse']
        r = dict(slot=c['slot'], model=c['model'], length=c['length'], width=c['width'],
                 peak_multiplier=c['native_peak_multiplier'], pulse_duration=c['cue_duration'],
                 source_kind='exact_fullfield_reuse' if reuse else 'new_native', state='not_started')
        p = package / reuse['result_path'] if reuse else folder / 'result.json'
        if p.exists():
            result = h.read(p)
            if not reuse:
                assert result['case'] == c
                for key, filename in (('cue_sha256', 'cue.npz'), ('recovery_sha256', 'recovery.npz'), ('regional_readouts_sha256', 'regional_readouts.csv')):
                    assert n.sha(folder / filename) == result[key]
            r.update(state='complete', **{k: result['legacy_diagnostics'][k] for k in
                     ('orientation_success', 'full_polarity_success', 'field_tail_pass')})
            r['actual_cue_path'] = reuse['cue_path'] if reuse else str((folder / 'cue.npz').relative_to(package))
            r['actual_recovery_path'] = reuse['recovery_path'] if reuse else str((folder / 'recovery.npz').relative_to(package))
        elif (folder / 'error.json').exists():
            r['state'] = 'error'
        elif (folder / 'status.json').exists():
            r['state'] = h.read(folder / 'status.json')['state']
        rows.append(r)
    h.write(out / 'results/condition_status.json', rows)
    done = sum(r['state'] == 'complete' for r in rows)
    errors = sum(r['state'] == 'error' for r in rows)
    h.write(out / 'results/latest_status.json', dict(generated_at=n.now(), total=72, complete=done == 72,
         completed=done, exact_reuses=3, new_completed=done - 3, pending=72 - done - errors, errors=errors,
         producer='finished' if done == 72 else 'execution_errors' if errors else 'partial_or_not_started'))
    print(json.dumps(h.read(out / 'results/latest_status.json'), indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('init', 'preflight', 'worker', 'collect'))
    parser.add_argument('--out', type=Path, default=DEFAULT)
    parser.add_argument('--task', type=int)
    args = parser.parse_args()
    out = args.out.resolve()
    if args.action == 'init':
        initialize(out)
    elif args.action == 'preflight':
        m = verify(out)
        print(json.dumps(dict(passed=True, protocols=len(m['cases']), new_tasks=len(m['jobs']),
                             manifest_sha256=n.sha(out / 'manifest.json')), indent=2))
    elif args.action == 'collect':
        collect(out)
    else:
        if args.task is None or not 0 <= args.task < 69:
            raise ValueError('task must be0..68')
        try:
            worker(out, args.task)
        except Exception:
            slot = h.read(out / 'manifest.json')['jobs'][args.task]['case_slot']
            h.write(out / 'workers' / f'{slot:03d}' / 'error.json',
                    dict(generated_at=n.now(), state='error', traceback=traceback.format_exc()))
            raise


if __name__ == '__main__':
    main()
