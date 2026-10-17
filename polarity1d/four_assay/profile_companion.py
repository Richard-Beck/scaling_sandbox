"""Rac/active profiles accompanying the revision-11 four-component timing figure.

Uses retained histories and response traces; reconstructs only successful finite
pulses from verified release checkpoints, with campaign pulse-end validation.
Run on SLURM (QoS small) with the existing scaling_sandbox Python environment.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, gzip, hashlib, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[1]
sys.path.insert(0, str(CODE))
import runner as r

MODELS = r.MODELS
LABELS = ['Wave-pinning', 'Otsuji', 'Goryachev', 'LEGI', 'deBelly', 'Holmes 3', 'Spring']
METRICS = ['conditioning', 'release', 'gradient_reversal', 'finite_countercue']
TITLES = ['Conditioning: 1% arrival', 'Release: 1% arrival',
          'Gradient reversal response', 'Countercue pulse duration']
COLORS = ['#176d9c', '#db7b24', '#27804b']
STYLES = [':', '--', '-']
FRACTIONS = np.array([0., .5, .99])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=CODE / 'results')
    parser.add_argument('--length-policy', choices=['nearby', 'exact'], default='nearby')
    args = parser.parse_args()
    out = args.results.resolve()
    status = json.loads((out / 'status.json').read_text())
    assert status['complete'] and status['source_revision'] == r.PROTOCOL['revision'] == 11
    records = json.load(gzip.open(out / 'records.json.gz', 'rt'))
    selections = {m: [7., 17., 27.] for m in MODELS}
    if args.length_policy == 'nearby':
        selections.update(holmes_model3=[10., 16., 26.], spring=[.7, 1.7, 2.6])
    sources = {}
    def source(path):
        path = path.resolve()
        sources[str(path.relative_to(ROOT))] = sha(path)
        return path
    for p in [Path(__file__), out / 'records.json.gz', out / 'status.json', CODE / 'protocol.json',
              CODE / 'runner.py', *sorted((CODE / 'models').glob('*.py')),
              *sorted((ROOT / 'polarity1d/src/polarity1d').glob('*.py'))]:
        source(p)
    def resolve(name):
        return source(ROOT / name)
    def select(items, model, length, assay):
        return next((x for x in items if x['model'] == model and x['length'] == length
                     and x['assay'] == assay), None)
    def checkpoint(prep, model):
        path = resolve(prep['checkpoint'])
        assert sha(path) == prep['checkpoint_sha256']
        return r.read_checkpoint(path, model)
    def field_name(model):
        return {'debelly': 'rac', 'holmes_model3': 'rac_active', 'legi': 'response'}.get(model, 'active')
    def interpolate(time, values, times):
        assert time[0] <= times.min() <= times.max() <= time[-1]
        return np.array([np.interp(times, time, col) for col in values.T]).T

    field_dir = out / 'length_response_profiles_fields'
    field_dir.mkdir(exist_ok=True)
    panels, metadata, checks = {}, [], []
    for model in MODELS:
        for length in selections[model]:
            gp = select(records['preparations'], model, length, 'gradient_reversal')
            fp = select(records['preparations'], model, length, 'finite_countercue')
            if gp is None:
                for metric in METRICS:
                    rec = dict(model=model, length=length, metric=metric, status='length unavailable',
                               indicated_time=None, times=[], source=None, field=field_name(model))
                    metadata.append(rec)
                    panels[(model, length, metric)] = (rec, None)
                continue
            gstate, gprep = checkpoint(gp, model)
            fstate, fprep = checkpoint(fp, model)
            assert fp['phases'][0]['history'] == gp['phases'][0]['history']
            assert fp['phases'][1]['input_checkpoint_sha256'] == gp['checkpoint_sha256']
            for metric in METRICS:
                endpoint = None
                if metric in ['conditioning', 'release']:
                    phase = gp['phases'][0] if metric == 'conditioning' else fp['phases'][1]
                    origin = resolve(phase['history'])
                    endpoint = float(phase['time_to_settle'])
                    times = endpoint * FRACTIONS
                    with np.load(origin) as z:
                        t, x = z['time'], z['x']
                        u = z['field_' + field_name(model)]
                        profiles = interpolate(t, u, times)
                        start = gprep.active[-1] if metric == 'release' else None
                        if start is not None:
                            np.testing.assert_allclose(profiles[0], start, rtol=1e-9, atol=1e-10)
                    mode = 'success'
                else:
                    case = select(records['cases'], model, length, metric)
                    mode = case['status']
                    prep = gprep if metric == 'gradient_reversal' else fprep
                    if mode != 'success':
                        # No measured clock exists. Show only the actual initial profile.
                        times, x, profiles = np.array([0.]), prep.x, prep.active[-1:]
                        origin = resolve((gp if metric == 'gradient_reversal' else fp)['checkpoint'])
                    elif metric == 'gradient_reversal':
                        origin = resolve(case['diagnostic_trace'])
                        endpoint = float(case['endpoint'])
                        times = endpoint * FRACTIONS
                        with np.load(origin) as z:
                            x = z['x']
                            profiles = interpolate(z['time'], z['active'], times)
                        np.testing.assert_allclose(profiles[0], prep.active[-1], rtol=1e-9, atol=1e-10)
                    else:
                        endpoint = float(case['endpoint'])
                        times = endpoint * FRACTIONS
                        origin = resolve(fp['checkpoint'])
                        # Match original integration/output settings, then interpolate.
                        pulse = r.run(model, length, fstate, r.Cue('counter', model, stop=endpoint),
                                      endpoint, interval=min(2.5, endpoint / 20), n=fp['cells'])
                        x = pulse.x
                        profiles = interpolate(pulse.time, pulse.active, times)
                        np.testing.assert_allclose(profiles[0], prep.active[-1], rtol=1e-9, atol=1e-10)
                        # Validate this reconstruction against the frozen trial's full pulse-end profile.
                        matched = False
                        for path in sorted((CODE / 'work/trials').glob('*.json')):
                            chunk = json.loads(path.read_text())
                            if chunk.get('task_index') != case['task_index']:
                                continue
                            diag = ROOT / chunk['diagnostic_trace']
                            with np.load(diag) as z:
                                indices = np.flatnonzero(np.isclose(z['durations'], endpoint, rtol=0, atol=1e-10))
                                if not len(indices):
                                    continue
                                expected = z['pulse_final_active'][indices[0]]
                                np.testing.assert_allclose(pulse.active[-1], expected, rtol=1e-7, atol=1e-8)
                                checks.append(dict(model=model, length=length, duration=endpoint,
                                                   maximum_pulse_end_error=float(np.max(abs(pulse.active[-1] - expected)))))
                            source(path); source(diag)
                            matched = True
                            break
                        if not matched:
                            raise RuntimeError(f'Missing original pulse-end diagnostic: {model} {length}')
                assert profiles.shape == (len(times), len(x)) and np.isfinite(profiles).all()
                # Native concentrations; normalized spatial coordinates allow lengths to be overlaid.
                q = x / length
                assert 0 <= q[0] < q[-1] <= 1 and np.all(np.diff(q) > 0)
                rec = dict(model=model, length=length, metric=metric, status=mode,
                           indicated_time=endpoint, times=times.tolist(), source=str(origin.relative_to(ROOT)),
                           field=field_name(model))
                metadata.append(rec)
                panels[(model, length, metric)] = (rec, (q, profiles))
                np.savez_compressed(field_dir / f'{model}_{length:g}_{metric}.npz',
                                    position_fraction=q, nominal_x=x, times=times, profile=profiles,
                                    indicated_time=np.nan if endpoint is None else endpoint)
                print(f'{model} L={length:g} {metric}: {mode}, T={endpoint}', flush=True)

    plt.rcParams.update({'font.size':8, 'axes.spines.top':False, 'axes.spines.right':False,
                         'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig, axes = plt.subplots(7, 4, figsize=(20, 19), sharex=True)
    for i, (model, label) in enumerate(zip(MODELS, LABELS)):
        rowvalues = []
        for j, metric in enumerate(METRICS):
            ax = axes[i, j]
            annotations = []
            for color, length in zip(COLORS, selections[model]):
                rec, data = panels[(model, length, metric)]
                endpoint = rec['indicated_time']
                text = f'L={length:g}: ' + (f'T={endpoint:g}' if endpoint is not None else rec['status'] + '; initial only' if data else rec['status'])
                annotations.append((color, text))
                if data:
                    q, profiles = data
                    rowvalues.extend(profiles.ravel())
                    # Draw initial last so coincident zero-clock snapshots retain their dotted signature.
                    for k in reversed(range(len(profiles))):
                        ax.plot(q, profiles[k], color=color, ls=STYLES[k], lw=1.4, alpha=.9)
            for k, (color, text) in enumerate(annotations):
                ax.text(.98, .98-k*.09, text, transform=ax.transAxes, ha='right', va='top',
                        color=color, fontsize=7, bbox=dict(facecolor='white', alpha=.8, edgecolor='none', pad=.5))
            ax.set_xlim(0, 1)
            ax.set_xticks([0, .5, 1])
            ax.tick_params(labelsize=7)
            ax.set_ylabel(('Rac' if model in ['debelly', 'holmes_model3'] else 'Response' if model == 'legi' else 'Active') + ' (native)')
            if i == 0:
                ax.set_title(TITLES[j], fontsize=10)
            if i == 6:
                ax.set_xlabel('Position / length (material coordinate for Spring)', fontsize=8)
        if rowvalues:
            lo, hi = min(rowvalues), max(rowvalues)
            span = max(hi-lo, abs(hi)*.05, 1e-6)
            for ax in axes[i]:
                ax.set_ylim(max(0, lo-.05*span), hi+.38*span)
        axes[i, 0].set_ylabel(label + '\n' + axes[i, 0].get_ylabel(), fontsize=9)
    stamp = datetime.now(timezone.utc).isoformat()
    fig.suptitle('Full spatial Rac / active profiles at initial, 50% and 99% of each measured time\n'
                 'Native concentrations; common y-axis scale within each model. T uses that length’s timing-figure readout.', fontsize=12)
    handles = [Line2D([0], [0], color=c, lw=2, label=f'Length choice {i+1}') for i,c in enumerate(COLORS)]
    handles += [Line2D([0], [0], color='.25', ls=s, lw=1.6, label=t)
                for s,t in zip(STYLES, ['Initial (t=0)', 't=0.50 T', 't=0.99 T'])]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5,.025), ncol=6, frameon=False, fontsize=9)
    fig.text(.5,.015, 'Censored/ineligible: initial only. T=0: all three profiles coincide. Countercue: during the shortest passing pulse; recovery excluded.',
             ha='center', fontsize=8)
    fig.text(.5,.005, f'Generated {stamp} · results {status["generated_at"]} · producer finished', ha='center', fontsize=7)
    fig.subplots_adjust(left=.07, right=.99, bottom=.065, top=.94, wspace=.34, hspace=.43)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(out / f'length_response_profiles.{ext}', dpi=200)
    plt.close(fig)
    with (out / 'length_response_profiles_times.csv').open('w') as f:
        rows = [dict(model=rec['model'], length=rec['length'], metric=rec['metric'], status=rec['status'],
                     indicated_time=rec['indicated_time'], initial=rec['times'][0] if rec['times'] else None,
                     half=rec['times'][1] if len(rec['times']) > 1 else None,
                     ninety_nine_percent=rec['times'][2] if len(rec['times']) > 1 else None) for rec in metadata]
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    provenance = dict(generated_at=stamp, results_generated_at=status['generated_at'], source_revision=11,
                      selections=selections, length_policy=args.length_policy, panels=metadata,
                      finite_reconstruction_checks=checks, sources=sources,
                      sampling='Pointwise linear interpolation of native histories at exactly 0, 0.50 T, 0.99 T; finite pulses reconstructed with original solver and sampling, validated against original pulse-end profiles.',
                      concentrations='Native rac for deBelly, rac_active for Holmes 3, response for LEGI, active for other models including Spring.',
                      spatial_axis='x / nominal length; Spring material fraction. No concentration normalization.',
                      unavailable_clocks='Only actual initial profile shown; no invented endpoint or horizon substituted.',
                      outputs={p.name:sha(p) for p in [*field_dir.glob('*.npz'),
                               *[out / f'length_response_profiles.{ext}' for ext in ['png','pdf','svg']],
                               out / 'length_response_profiles_times.csv']})
    (out / 'length_response_profiles_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Saved length_response_profiles.png/pdf/svg', flush=True)

if __name__ == '__main__':
    main()
