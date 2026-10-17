"""Independent long-holdout check and a compact review of continuation references."""
from pathlib import Path
import argparse, csv, json, os, shlex, subprocess
import numpy as np
import steady_states as s
import runner as r

C = s.CAMPAIGN


def validate(slot):
    manifest = s.verify()
    frozen = json.loads((C / 'validation_manifest.json').read_text())
    if r.sha(Path(__file__).resolve()) != frozen['review_source_sha256']:
        raise RuntimeError('Validation source changed')
    task = manifest['tasks'][slot]
    record = next(p for p in (C / 'records').glob(f'{slot:03d}_*.json'))
    item = json.loads(record.read_text())
    if not item['converged']:
        raise RuntimeError('Cannot validate unresolved reference')
    key, length = item['model'], item['length']
    state, old = r.read_checkpoint(C / item['profile'], key)
    cue = r.Cue('left' if item['event'] == 'conditioning' else 'zero', key)
    tr = r.run(key, length, state, cue, 10000., interval=5, n=task['cells'])
    before = {k: v[-1] for k, v in s.fields(old, key, length).items()}
    after = {k: v[-1] for k, v in s.fields(tr, key, length).items()}
    error = s.distance(before, after, item['fixed_scales'])
    maximum_excursion = max(float(np.max(abs(y - before[k]))) / item['fixed_scales'][k]
                            for k, y in s.fields(tr, key, length).items())
    # Recompute arrival against the later endpoint to quantify target sensitivity.
    with np.load(C / item['history']) as z:
        t = z['time']
        d = np.zeros(len(t))
        for k, scale in item['fixed_scales'].items():
            y = z['field_' + k]
            d = np.maximum(d, abs(y - after[k]).reshape(len(t), -1).max(axis=1) / scale)
        outside = np.flatnonzero(d > s.ARRIVAL)
        arrival = float(t[outside[-1] + 1]) if len(outside) else 0.
    stem = Path(item['profile']).stem
    keep = tr.time >= tr.time[-1] - s.CHUNK
    tr.time = tr.time[keep]
    for name in ['active', 'orientation', 'contrast', 'left', 'right']:
        setattr(tr, name, getattr(tr, name)[keep])
    tr.fields = {k: v[keep] for k, v in tr.fields.items()}
    r.save_checkpoint(C / 'validation' / (stem + '.npz'), key, tr.state, tr)
    final = s.check(tr, key, length, item['fixed_scales'])
    r.write(C / 'validation' / (stem + '.json'), dict(slot=slot, model=key, length=length,
        event=item['event'], holdout_time=10000, endpoint_error_by_field=error,
        maximum_scaled_endpoint_shift=max(error.values()), maximum_scaled_excursion=maximum_excursion,
        accepted=maximum_excursion <= .001 and final['passed'], final_window=final,
        later_reference_arrival=arrival, arrival_shift=arrival - item['arrival_time'],
        generated_at=r.now(), slurm_job=os.environ.get('SLURM_JOB_ID')))


def launch():
    s.verify()
    (C / 'validation').mkdir()
    r.write(C / 'validation_manifest.json', dict(created_at=r.now(), holdout_time=10000,
        excursion_tolerance=.001, final_window_tolerance=s.TOL,
        review_source_sha256=r.sha(Path(__file__).resolve())))
    python = json.loads((s.BASE / 'environment.json').read_text())['python']
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; ' + shlex.join(
        [python, str(Path(__file__).resolve()), 'validate']) + ' --slot "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(['sbatch', '--parsable', '--qos=small', '--ntasks=1',
        '--cpus-per-task=1', '--mem=3G', '--time=01:00:00', '--array=0-307',
        '--job-name=polarity-ss-check', f'--chdir={r.ROOT}',
        f'--output={C}/logs/validation_%A_%a.out', f'--error={C}/logs/validation_%A_%a.err',
        '--wrap', command], text=True).strip()
    r.write(C / 'validation_launch.json', dict(job=job, submitted_at=r.now()))
    print(job)


def report():
    s.verify()
    items = [json.loads(p.read_text()) for p in sorted((C / 'records').glob('*.json'))]
    validations = [json.loads(p.read_text()) for p in sorted((C / 'validation').glob('*.json'))]
    if sorted(p['slot'] for p in validations) != list(range(308)):
        raise RuntimeError('Incomplete validation')
    summary = dict(generated_at=r.now(), total=len(validations), passed=sum(p['accepted'] for p in validations),
        maximum_scaled_excursion=max(p['maximum_scaled_excursion'] for p in validations),
        maximum_absolute_arrival_shift=max(abs(p['arrival_shift']) for p in validations),
        failures=[p for p in validations if not p['accepted']])
    r.write(C / 'validation_summary.json', summary)
    rows = [{k: v for k, v in p.items() if not isinstance(v, (list, dict))} for p in validations]
    with (C / 'validation.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 7, figsize=(17, 6), sharey='row')
    colors = ['#3974a8', '#c17534']
    for j, key in enumerate(r.MODELS):
        for e, event in enumerate(['conditioning', 'release']):
            rr = sorted([p for p in items if p['model'] == key and p['event'] == event], key=lambda p: p['length'])
            x = [p['length'] for p in rr]
            axs[0, j].plot(x, [p['old_time'] if p['old_settled'] else np.nan for p in rr],
                           color=colors[e], ls='--', alpha=.6)
            axs[0, j].plot(x, [p['arrival_time'] for p in rr], color=colors[e], label=event)
            axs[1, j].plot(x, [max(p['maximum_scaled_endpoint_difference'] * 100, 1e-12) for p in rr], color=colors[e])
            for p in rr:
                if not p['old_settled']:
                    axs[0, j].scatter(p['length'], p['old_elapsed'], color=colors[e], marker='x')
        axs[0, j].set_title(key.replace('_', ' '), fontsize=10)
        axs[1, j].set_yscale('log')
        axs[1, j].axhline(1., color='grey', lw=.8, ls=':')
        axs[1, j].set_xlabel('Native length')
        for ax in axs[:, j]:
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=8)
    axs[0, 0].set_ylabel('Native time')
    axs[1, 0].set_ylabel('Endpoint difference (% of fixed scale)')
    axs[0, -1].legend(fontsize=8)
    fig.suptitle('Solid: 1% reference arrival; dashed: original confirmation time; x: original cap', fontsize=11)
    fig.tight_layout()
    fig.savefig(C / 'comparison.png', dpi=180)
    fig.savefig(C / 'comparison.pdf')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('action', choices=['launch', 'validate', 'report'])
    ap.add_argument('--slot', type=int)
    args = ap.parse_args()
    {'launch': launch, 'validate': lambda: validate(args.slot), 'report': report}[args.action]()
