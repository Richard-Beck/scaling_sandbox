"""Read-only native-field review of the eight fixed-pulse width contrasts."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'data/runs/revision_width_fullfield_20261019'
RETRY = ROOT / 'data/runs/revision_width_fullfield_gory_retry_20261019'
OUT = ROOT / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/width_fullfield_experiment_20261019'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load(p):
    return json.loads(p.read_text())


def csvwrite(p, rows):
    with p.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_cases():
    data = []
    sources = {}
    for run in (RUN, RETRY):
        if not (run / 'manifest.json').exists():
            continue
        m = load(run / 'manifest.json')
        sources[str((run / 'manifest.json').relative_to(ROOT))] = sha(run / 'manifest.json')
        for c in m['cases']:
            if run == RUN and c['model'] == 'goryachev':
                continue
            folder = run / 'workers' / f'{c["slot"]:02d}'
            result_path = folder / 'result.json'
            if not result_path.exists():
                continue
            result = load(result_path)
            assert result['case'] == c
            phases = {}
            for phase, filename in (('cue_on', 'cue.npz'), ('off', 'recovery.npz')):
                p = folder / filename
                digest = sha(p)
                assert digest == result['cue_sha256' if phase == 'cue_on' else 'recovery_sha256']
                sources[str(p.relative_to(ROOT))] = digest
                with np.load(p) as z:
                    phases[phase] = {k: z[k].copy() for k in z.files}
            with np.load(run / c['regional_input']) as z:
                masks = {k: z[k].copy() for k in z.files}
            sources[str((run / c['regional_input']).relative_to(ROOT))] = c['regional_input_sha256']
            data.append(dict(case=c, result=result, phases=phases, masks=masks,
                             source_run=str(run.relative_to(ROOT))))
    return data, sources


def metrics(d):
    c, z = d['case'], d['masks']
    core, roi = z['incumbent_core_mask'], z['cue_footprint_mask']
    dx = c['length'] / c['cells']
    initial = float(z['initial_state'][0, core].sum() * dx)
    curves = {}
    for phase, p in d['phases'].items():
        u, v = p['field_active'], p['field_inactive']
        curves[phase] = dict(time=p['time'], incumbent=u[:, core].sum(axis=1) * dx / initial,
                             challenger_peak=u[:, roi].max(axis=1),
                             inactive_mean=v.mean(axis=1))
    d['curves'] = curves
    on, off = curves['cue_on'], curves['off']
    final = d['phases']['off']['field_active'][-1]
    above = final > float(z['initial_halfcontrast_threshold'])
    regions = int(np.count_nonzero(np.diff(np.r_[False, above, False].astype(int)) == 1))
    row = dict(model=c['model'], length=c['length'], width=c['width'], width_label=c['width_label'],
               pulse_duration=c['cue_duration'], requested_peak=c['response_amplitude'],
               sampled_peak=c['sampled_peak'], sampled_spatial_integral=c['sampled_spatial_integral'],
               initial_core_active_amount=initial,
               incumbent_fraction_at_withdrawal=float(on['incumbent'][-1]),
               minimum_incumbent_fraction_during_cue=float(on['incumbent'].min()),
               cue_time_minimum_incumbent=float(on['time'][np.argmin(on['incumbent'])]),
               incumbent_fraction_final=float(off['incumbent'][-1]),
               challenger_peak_at_withdrawal=float(on['challenger_peak'][-1]),
               challenger_peak_final=float(off['challenger_peak'][-1]),
               final_above_initial_halfcontrast_regions=regions,
               final_native_active_peak=float(final.max()),
               source_run=d['source_run'], source_slot=c['slot'])
    for fraction in (.9, .75, .5, .25):
        hit = np.flatnonzero(on['incumbent'] <= fraction)
        row[f'first_cue_on_time_core_fraction_le_{fraction:g}'] = float(on['time'][hit[0]]) if hit.size else ''
    row.update({k: d['result']['legacy_diagnostics'][k] for k in
                ('orientation_success', 'full_polarity_success', 'field_tail_pass')})
    row['maximum_total_amount_relative_error'] = d['result']['maximum_total_amount_relative_error']
    return row


def page(model, pairs, rows):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.4), gridspec_kw=dict(width_ratios=(1.05, 1.05, 1.15)))
    fig.subplots_adjust(left=.065, right=.94, top=.84, bottom=.24, wspace=.56, hspace=.6)
    label = 'Wave-pinning' if model == 'wave_pinning' else 'Goryachev'
    fig.text(.065, .958, f'{label}: what changes when the same-peak cue is broader?', fontsize=16, weight='bold')
    fig.text(.065, .918, 'New full-field experiments from exact shared basals; two widths at the same pulse within each size.', fontsize=11)
    fig.text(.065, .887, 'Solid = original-core activity / initial amount (left axis). Dashed = native active peak inside actual cue support (right axis).', fontsize=10)
    colors = dict(narrow='#2171b5', broad='#983ca0')
    letters = 'ABCDEF'
    for row_index, pair in enumerate(pairs):
        c0 = pair[0]['case']
        for col in (0, 1):
            ax = axes[row_index, col]
            right = ax.twinx()
            phase = 'cue_on' if col == 0 else 'off'
            for d in pair:
                c, p = d['case'], d['curves'][phase]
                # Plot a subset only; all exact saved fields and primary curves remain archived.
                stride = max(1, len(p['time']) // 1600)
                ids = np.unique(np.r_[np.arange(0, len(p['time']), stride), len(p['time']) - 1])
                if phase == 'off':
                    ids = np.unique(np.r_[np.flatnonzero(p['time'] <= 100),
                                          np.searchsorted(p['time'], np.geomspace(100, 10000, 1200)),
                                          len(p['time']) - 1])
                color = colors[c['width_label']]
                ax.plot(p['time'][ids], p['incumbent'][ids], color=color, lw=1.7,
                        label=f'w = {c["width"]:g}')
                right.plot(p['time'][ids], p['challenger_peak'][ids], '--', color=color, lw=1.4)
            ax.set_ylim(-.025, 1.08)
            right.set_ylim(bottom=0)
            if col:
                ax.set_xscale('symlog', linthresh=10)
                ax.set_xlim(0, 10000)
                ax.set_xticks([0, 10, 100, 1000, 10000], ['0', '10', '100', '1000', '10000'])
                ax.set_xlabel('Time after withdrawal (native; linear ≤10, log >10)', fontsize=8.5)
                title = 'Cue-free continuation'
            else:
                ax.set_xlim(0, c0['cue_duration'])
                ax.set_xlabel('Time since cue onset (native)', fontsize=9)
                title = f'L = {c0["length"]:g}, pulse = {c0["cue_duration"]:g}'
                ax.legend(frameon=False, fontsize=8, loc='upper right')
            ax.set_ylabel('Original-core active amount / initial', fontsize=8.5)
            right.set_ylabel('Cue-ROI active peak (native)', fontsize=8.5)
            ax.set_title(f'{letters[row_index*3+col]}  {title}', loc='left', fontsize=11, weight='bold')
            ax.tick_params(labelsize=8)
            right.tick_params(labelsize=8)
            ax.grid(alpha=.15)
        ax = axes[row_index, 2]
        x = pair[0]['masks']['x']
        ax.plot(x, pair[0]['masks']['initial_state'][0], color='.5', lw=1.5, label='initial')
        for d in pair:
            c = d['case']
            color = colors[c['width_label']]
            ax.plot(x, d['phases']['cue_on']['field_active'][-1], color=color, lw=1.7,
                    label=f'w={c["width"]:g}, withdrawal')
            ax.plot(x, d['phases']['off']['field_active'][-1], ':', color=color, lw=2,
                    label=f'w={c["width"]:g}, +10000 off')
        ax.set_ylim(bottom=0)
        ax.set_xlim(0, c0['length'])
        ax.set_xlabel('Native interval position', fontsize=9)
        ax.set_ylabel('Native active concentration', fontsize=9)
        ax.set_title(f'{letters[row_index*3+2]}  Actual activity profiles', loc='left', fontsize=11, weight='bold')
        ax.tick_params(labelsize=8)
        ax.legend(frameon=False, fontsize=7, loc='upper right')
    selected = [r for r in rows if r['model'] == model]
    statements = []
    for length in sorted({r['length'] for r in selected}):
        rr = sorted([r for r in selected if r['length'] == length], key=lambda r: r['width'])
        statements.append(f'L{length:g}: ' + '; '.join(
            f'w{r["width"]:g}: scalar-strength {"PASS" if r["full_polarity_success"] else "FAIL"}, '
            f'field-tail {"PASS" if r["field_tail_pass"] else "FAIL"}' for r in rr))
    fig.text(.065, .117, '\n'.join(statements), fontsize=9, linespacing=1.7)
    fig.text(.065, .057, 'Regional activity loss is not identified molecular damage. Width also changes integrated input; native stimulus × inactive is not a calibrated optical actuator.', fontsize=8.5)
    fig.text(.065, .029, 'Same pulse within each width pair; pulses differ across the two sizes. N = 201; no spatial-convergence claim. Scalar gates and regional profiles are separate observations.', fontsize=8.5)
    return fig


def main():
    for directory in ('supporting/data', 'supporting/vectors', 'supporting/metadata', 'supporting/viewers'):
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    data, sources = read_cases()
    if not data:
        raise RuntimeError('No completed native fields')
    rows = [metrics(d) for d in data]
    csvwrite(OUT / 'supporting/data/actual_fullfield_summary.csv', rows)
    figures = []
    with PdfPages(OUT / 'supporting/vectors/figure2_width_fullfield_review.pdf') as pdf:
        for model in ('wave_pinning', 'goryachev'):
            dd = [d for d in data if d['case']['model'] == model]
            sizes = sorted({d['case']['length'] for d in dd})
            pairs = [sorted([d for d in dd if d['case']['length'] == size], key=lambda d:d['case']['width']) for size in sizes]
            if len(pairs) != 2 or any(len(p) != 2 for p in pairs):
                continue
            fig = page(model, pairs, rows)
            stem = f'figure2_{model}_fullfield'
            fig.savefig(OUT / f'{stem}.png', dpi=150)
            fig.savefig(OUT / 'supporting/vectors' / f'{stem}.svg')
            pdf.savefig(fig)
            plt.close(fig)
            figures.append(stem)
    (OUT / 'README.md').write_text('# Fixed-pulse cue-width experiments: actual native-field review\n\n'
        '[Current review PDF](supporting/vectors/figure2_width_fullfield_review.pdf)\n\n'
        + '\n\n'.join(f'![{stem}]({stem}.png)' for stem in figures) + '\n\n'
        'These are new full-field measurements from eight selected fixed-pulse protocols, not a finished Figure 2. '
        'Only completed pairs are shown. The original all-five incumbent organization remains the foundation: '
        '[unchanged prior organization](../finite_response_candidate_20261019/figure2_organization.png).\n\n'
        'In WP, the broader cue leaves less activity in the original core at withdrawal: 0.789→0.665 of its initial '
        'amount at L9 and 0.671→0.605 at L19. Both broader cues end with retained right polarity; both narrower cues '
        'return to the incumbent. The L19 challenger peaks at withdrawal are almost identical (1.4599/1.4615), '
        'despite the different eventual outcomes.\n\n'
        'In Goryachev, broadening lowers the withdrawal challenger peak at both sizes (12.025→7.433 at L9; '
        '22.014→15.519 at L25) and also lowers original-core activity (0.449→0.412; 0.446→0.320). '
        'The narrower cue wins at L9, whereas the broader cue wins at L25. Thus neither withdrawal challenger height '
        'nor incumbent-region depletion alone explains all four fates. The cue-free histories show the subsequent '
        'growth/return directly. These are conditional native-model observations, not a unique mechanism identification.\n\n'
        'Pulses match within each width pair but differ between size cohorts. The separately frozen next stage crosses '
        'physical width, native peak and pulse duration independently, using the same pulse grid at both sizes within '
        'each model; it is explicitly a two-model stage.\n\n'
        '[Exact numerical summary](supporting/data/actual_fullfield_summary.csv) · '
        '[Methods and limits](EVIDENCE.md) · [Source hashes](supporting/metadata/validation.json)\n')
    (OUT / 'EVIDENCE.md').write_text('# Evidence and readout scope\n\n'
        'The cue is the original right-end raised-cosine patch, with requested peak 0.1 (WP) or 6.4 (Goryachev), '
        'and unchanged stimulus × inactive activation transfer. Each pair starts from an identical copied, hash-checked '
        'corrected basal; support is physical in these fixed-length models even where inherited metadata says relative. '
        'The pulse is identical within each pair. It differs between size cohorts and is therefore not a controlled '
        'cross-size pulse comparison. Width changes integrated input as well as footprint.\n\n'
        'Primary plotted reads are native active concentration maxima within the actual cue footprint and native '
        'active amount in the fixed initial incumbent region, divided by that region’s initial amount. The region is '
        'the connected component above the initial background-plus-half-range threshold containing the initial maximum. '
        'It is a fixed spatial ROI, not a tracked site identity or identified molecular damage. Full native active and '
        'inactive arrays are retained every 0.1 cue unit and 0.25 off unit through 10000. Plotting may subsample lines; '
        'metrics use every saved frame. First cue-on depletion crossings in the table are sampled threshold entries, '
        'not persistent state changes, optical latency calibrations or pulse-duration thresholds.\n\n'
        'Original scalar orientation and retained-strength gates use the unchanged 2.5-unit grid over the final 500 '
        'off units. Native full-field tail stationarity is reported separately; neither scalar success nor one regional '
        'curve proves complete field replacement. The inherited N=201 grid is not a spatial convergence study.\n\n'
        'The first four WP protocols directly replay exact old scalar trials. Gory74/815 were not exact prior tested '
        'pulses. Their outcomes must be measured without assuming global monotonicity. The initial Gory run failed only '
        'because an exact-output assertion rejected BDF interpolation roundoff below 5e-22; the explicitly linked retry '
        'preserves exact solver inputs and records returned-frame absolute/scaled roundoff separately. No model, cue, '
        'parameter or scientific integration rule changed.\n')
    sources[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    artifacts = {str(p.relative_to(OUT)): sha(p) for p in OUT.rglob('*') if p.is_file()
                 and p.name != 'validation.json'}
    (OUT / 'supporting/metadata/validation.json').write_text(json.dumps(dict(
        completed_native_conditions=len(data), requested_conditions=8,
        status='complete' if len(data) == 8 else 'partial',
        source_hashes=sources, artifact_hashes=artifacts), indent=2) + '\n')
    print(json.dumps(dict(output=str(OUT), completed=len(data), pages=len(figures)), indent=2))


if __name__ == '__main__':
    main()
