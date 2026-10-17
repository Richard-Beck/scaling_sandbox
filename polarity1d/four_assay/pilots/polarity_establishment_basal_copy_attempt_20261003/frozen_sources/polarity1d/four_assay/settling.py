"""One continuation/arrival ruleset for both events in the preparation chain."""
import json
import numpy as np


def monitored(tr, key, length):
    out = {k: np.asarray(v) for k, v in tr.fields.items()}
    if key == 'debelly':
        out = {k: out[k] for k in ['rac', 'rho', 'mca_bound', 'mca_unbound',
                                  'displacement', 'velocity', 'tension', 'boundary_separation']}
        out['displacement'] = out['displacement'] - out['displacement'].mean(axis=1)[:, None]
        out['physical_length'] = length + out.pop('boundary_separation')
    return out


def deviations(fields, target, scales):
    count = len(next(iter(fields.values())))
    distance = np.zeros(count)
    for name, y in fields.items():
        error = np.abs(y - target[name]).reshape(count, -1).max(axis=1) / scales[name]
        distance = np.maximum(distance, error)
    return distance


def diagnostics(tr, key, length, scales, rule):
    keep = tr.time >= tr.time[-1] - rule['chunk']
    t = tr.time[keep]
    metrics = {}
    for name, y in monitored(tr, key, length).items():
        y = y[keep]
        variation = float(np.ptp(y, axis=0).max()) / scales[name]
        drift = float(abs(y[-1] - y[-2]).max()) / (t[-1] - t[-2]) * rule['chunk'] / scales[name]
        metrics[name] = dict(relative_range=variation, projected_relative_drift=drift)
    return dict(steady=all(max(v.values()) <= rule['tolerance'] for v in metrics.values()), metrics=metrics)


def establish(r, key, length, state, kind, n, history_path):
    rule = r.protocol()['steady_state']
    event = 'conditioning' if kind == 'left' else 'release'
    rows = json.loads((r.CODE / rule['fixed_scales']).read_text())
    scales = next(p['fixed_scales'] for p in rows if p['model'] == key and p['length'] == length and p['event'] == event)
    cue = r.Cue(kind, key)
    trace = None
    elapsed = 0.
    checks = []
    consecutive = 0
    verified = False
    while elapsed + rule['chunk'] <= rule['cap']:
        tr = r.run(key, length, state, cue, rule['chunk'], interval=rule['sample_interval'], n=n)
        state = tr.state
        metric = diagnostics(tr, key, length, scales, rule)
        trace = r.concatenate(trace, tr, elapsed)
        elapsed += rule['chunk']
        checks.append(dict(time=elapsed, kind='convergence', **metric))
        consecutive = consecutive + 1 if metric['steady'] else 0
        print(key, length, event, elapsed, metric['steady'], flush=True)
        if consecutive < rule['consecutive_windows']:
            continue
        if elapsed + rule['holdout'] > rule['cap']:
            break
        target = {k: v[-1].copy() for k, v in monitored(trace, key, length).items()}
        hold = r.run(key, length, state, cue, rule['holdout'], interval=rule['sample_interval'], n=n)
        state = hold.state
        excursion = float(deviations(monitored(hold, key, length), target, scales).max())
        final = diagnostics(hold, key, length, scales, rule)
        trace = r.concatenate(trace, hold, elapsed)
        elapsed += rule['holdout']
        verified = excursion <= rule['holdout_tolerance'] and final['steady']
        checks.append(dict(time=elapsed, kind='holdout', maximum_scaled_excursion=excursion,
                           verified=verified, **final))
        if verified:
            break
        consecutive = 0
    target = {k: v[-1].copy() for k, v in monitored(trace, key, length).items()}
    errors = deviations(monitored(trace, key, length), target, scales)
    outside = np.flatnonzero(errors > rule['arrival_tolerance'])
    arrival = (float(trace.time[outside[-1] + 1]) if len(outside) and outside[-1] + 1 < len(errors)
               else 0. if not len(outside) else None) if verified else None
    if arrival is not None and trace.time[-1] - arrival < rule['chunk']:
        verified = False
        arrival = None
    np.savez_compressed(history_path, time=trace.time, x=trace.x, maximum_scaled_distance=errors,
                        **{'field_' + k: v for k, v in monitored(trace, key, length).items()})
    # Save only the final verification window in the restart checkpoint.
    keep = trace.time >= elapsed - rule['chunk']
    trace.time = trace.time[keep]
    for name in ['active', 'orientation', 'contrast', 'left', 'right']:
        setattr(trace, name, getattr(trace, name)[keep])
    trace.fields = {k: v[keep] for k, v in trace.fields.items()}
    phase = dict(cue=kind, event=event, settled=verified, time_to_settle=arrival,
                 arrival_time=arrival, elapsed=elapsed, cap=rule['cap'],
                 confirmation_window=rule['chunk'], fixed_scales=scales, checks=checks,
                 arrival_tolerance=rule['arrival_tolerance'], history=str(history_path.relative_to(r.ROOT)),
                 time_definition='First sampled entry into 1% neighborhood that persists through verification',
                 reference_definition='Final endpoint of accepted continuation plus 10000-unit holdout')
    return state, trace, phase


def prepare_pair(r, slot):
    pid = r.verify()
    pairs = [t for t in r.tasks() if t['assay'] == 'gradient_reversal']
    task = pairs[slot]
    key, length, n = task['model'], task['length'], task['cells']
    gradient = task
    finite = next(t for t in r.tasks() if t['model'] == key and t['length'] == length and t['assay'] == 'finite_countercue')
    gidx, fidx = gradient['task_index'], finite['task_index']
    cpath = r.HERE / 'checkpoints' / f'task_{gidx:04d}.npz'
    zpath = r.HERE / 'checkpoints' / f'task_{fidx:04d}.npz'
    state, tr, conditioning = establish(r, key, length, r.initial(key, length, n), 'left', n,
        r.HERE / 'histories' / f'conditioning_{gidx:04d}.npz')
    r.save_checkpoint(cpath, key, state, tr)
    conditioning.update(reference_checkpoint=str(cpath.relative_to(r.ROOT)),
                        reference_checkpoint_sha256=r.sha(cpath), input_definition='Native resting initial state')

    def save(task, phases, path, trace):
        settled = all(p['settled'] for p in phases)
        eligible, criteria = r.left_eligible(trace, key, finite=task['assay'] == 'finite_countercue') if settled else (False, {})
        record = dict(**task, phases=phases, settled=settled, eligible=eligible,
                      eligibility_diagnostics=criteria, status='qualified' if eligible else 'ineligible' if settled else 'unsettled',
                      total_preparation_time=sum(p['elapsed'] for p in phases),
                      reported_preparation_time=sum(p['time_to_settle'] for p in phases) if settled else None,
                      cap_per_phase=r.protocol()['steady_state']['cap'], relative_tolerance=r.protocol()['steady_state']['tolerance'],
                      provenance_id=pid, checkpoint=str(path.relative_to(r.ROOT)), checkpoint_sha256=r.sha(path),
                      generated_at=r.now(), slurm_job=r.os.environ.get('SLURM_ARRAY_JOB_ID'), slurm_slot=slot)
        r.write(r.HERE / 'preparations' / f'task_{task["task_index"]:04d}.json', record)

    save(gradient, [conditioning], cpath, tr)
    if conditioning['settled']:
        state, tr, release = establish(r, key, length, state, 'zero', n,
            r.HERE / 'histories' / f'release_{fidx:04d}.npz')
        r.save_checkpoint(zpath, key, state, tr)
        release.update(reference_checkpoint=str(zpath.relative_to(r.ROOT)), reference_checkpoint_sha256=r.sha(zpath),
                       input_definition='Fully settled shared conditioning', input_checkpoint_sha256=r.sha(cpath))
        save(finite, [conditioning, release], zpath, tr)
    else:
        r.save_checkpoint(zpath, key, state, tr)
        save(finite, [conditioning], zpath, tr)
    r.verify()
