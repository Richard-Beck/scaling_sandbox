"""Polarity-establishment report: all four metrics grouped within each model."""
import argparse
import csv
import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

MODELS = ['wave_pinning', 'otsuji', 'goryachev', 'legi', 'debelly', 'holmes_model3', 'spring']
LABELS = ['Wave-pinning', 'Otsuji', 'Goryachev', 'LEGI', 'deBelly', 'Holmes 3', 'Spring']
SETS = [('cue_duration', 'cue_duration', 'Cue duration', 'native time', True),
        ('total_establishment_time', 'total_establishment_time', 'Total establishment time', 'native time', True),
        ('final_polarity', 'selected_final_polarity', 'Final cue-aligned polarity', 'polarity score', False),
        ('apparent_endpoint_count', 'apparent_endpoint_count', 'Apparent number of settled endpoints', 'count', False)]


def edges(values, logarithmic=False):
    v = np.log(values) if logarithmic else np.asarray(values, float)
    if len(v) == 1:
        e = np.array([v[0] - .5, v[0] + .5])
    else:
        mid = (v[:-1] + v[1:]) / 2
        e = np.r_[v[0] - (mid[0] - v[0]), mid, v[-1] + (v[-1] - mid[-1])]
    return np.exp(e) if logarithmic else e


def plot(out, figures=None):
    result = out / 'results'
    status = json.loads((result / 'status.json').read_text())
    if not status['complete']:
        raise RuntimeError('Final collection requires complete execution')
    with (result / 'endpoints.csv').open() as f:
        rows = list(csv.DictReader(f))
    endpoint_hash = hashlib.sha256((result / 'endpoints.csv').read_bytes()).hexdigest()
    manifest = json.loads((out / 'manifest.json').read_text())
    destination = Path(figures).resolve() if figures else result
    destination.mkdir(parents=True,exist_ok=True)
    conditions = manifest.get('conditions', [dict(condition='fixed_10pct', width_mode='physical', fraction=.1),
                                            dict(condition='fixed_20pct', width_mode='physical', fraction=.2)])
    combined = PdfPages(destination / 'polarity_establishment_by_model.pdf')
    for model, label in zip(MODELS, LABELS):
        filename = 'model_' + model
        fig, axes = plt.subplots(4, len(conditions), figsize=(7.5 * len(conditions), 15), constrained_layout=True)
        fig.get_layout_engine().set(rect=(0., .13, 1., .81))
        model_rows = [r for r in rows if r['model'] == model]
        for i, (_, field, title, unit, log) in enumerate(SETS):
            vals = [float(r[field]) for r in model_rows if r[field] and (not log or float(r[field]) > 0)]
            if field == 'selected_final_polarity':
                norm = Normalize(0, 1)
            elif log:
                low, high = (min(vals), max(vals)) if vals else (.01, 1.)
                norm = LogNorm(low, high if high > low else low * 10)
            else:
                norm = Normalize(1, max(2, max(vals)) if vals else 2)
            for j, condition in enumerate(conditions):
                ax = axes[i, j]
                fraction = condition['fraction']
                chosen = [r for r in model_rows if r.get('condition', 'fixed_10pct' if float(r['width_fraction_at_reference']) == .1 else 'fixed_20pct') == condition['condition']]
                lengths = sorted({float(r['length']) for r in chosen})
                amplitudes = sorted({float(r['amplitude']) for r in chosen})
                matrix = np.full((len(amplitudes), len(lengths)), np.nan)
                for row in chosen:
                    x, y = float(row['length']), float(row['amplitude'])
                    if row[field] and (not log or float(row[field]) > 0):
                        matrix[amplitudes.index(y), lengths.index(x)] = float(row[field])
                cmap = plt.get_cmap('viridis').copy()
                cmap.set_bad('#dedede')
                length_edges, amplitude_edges = edges(lengths), edges(amplitudes, True)
                mesh = ax.pcolormesh(length_edges, amplitude_edges,
                                     np.ma.masked_invalid(matrix), cmap=cmap, norm=norm, shading='flat')
                ax.set_yscale('log')
                ax.set_yticks(amplitudes, [f'{v:g}' for v in amplitudes])
                ax.yaxis.set_minor_formatter(NullFormatter())
                ax.set_xticks(lengths[::max(1, len(lengths) // 8)])
                ax.set_xlabel('L (native length)' + ('; not μm' if model == 'spring' else ''))
                ax.set_ylabel('A (native peak amplitude)')
                width_label = ('Relative width 20% of instantaneous L' if model == 'spring' else 'Relative width 20% of L') if condition['width_mode'] == 'relative' else f"Fixed width {float(chosen[0]['width']):g} ({fraction:.0%} of reference L)"
                ax.set_title(f"{title}\n{width_label}", fontsize=11)
                for row in chosen:
                    x, y = float(row['length']), float(row['amplitude'])
                    xi, yi = lengths.index(x), amplitudes.index(y)
                    dx = length_edges[xi + 1] - length_edges[xi]
                    dy = np.log(amplitude_edges[yi + 1] / amplitude_edges[yi])
                    marker = {'unpolarized': 'x', 'unsettled': '+', 'no_qualifying_endpoint': 's'}.get(row['status'])
                    if marker:
                        ax.plot(x, y, marker=marker, color='black', ms=4, mfc='none', linestyle='none')
                    if row['upper_optimum_unresolved'] == 'True' or row['unsettled_trials'] != '0':
                        ax.plot(x + .22 * dx, y * np.exp(.20 * dy), marker='^', color='white', mec='black', ms=4, linestyle='none')
                    if row['duration_boundary_unresolved'] == 'True':
                        ax.plot(x - .22 * dx, y * np.exp(.20 * dy), marker='o', color='black', mfc='none', ms=5, linestyle='none')
                    if row['cue_duration'] and float(row['cue_duration']) == 0:
                        ax.plot(x, y, marker='D', color='black', mfc='white', ms=5, linestyle='none')
                cb = fig.colorbar(mesh, ax=ax, label=unit, fraction=.035, pad=.02)
                if field == 'apparent_endpoint_count':
                    cb.locator = MaxNLocator(integer=True)
                    cb.update_ticks()
        fig.suptitle(label + ': polarity establishment\nRest → added left patch → verified basal release', fontsize=17, y=.995)
        handles = [Line2D([], [], marker=m, color='black', mfc='white' if m in ['^', 'D'] else 'none', linestyle='none', label=l)
                   for m, l in [('x', 'Unpolarized'), ('+', 'No verified endpoint'),
                                ('s', 'No qualifying polarized trial'), ('^', 'Upper optimum or some releases unresolved'),
                                ('o', 'Cue-duration threshold unresolved to 1%'), ('D', 'No polarizing cue required (D=0)')]]
        fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .068), ncol=3, fontsize=9)
        fig.text(.5, .01, 'Strongest observed cue-aligned endpoint, then shortest tested duration within 1% of its strength.\n'
                      'Counts: approximate complete-link profile clusters; zero-duration control included; no extra search for counts.\n'
                      'Colors share a scale between all conditions within each model; native time units differ between models.\n'
                      'Generated ' + status['generated_at'] + '; producer finished.', fontsize=9, ha='center', va='bottom')
        for extension in ['png', 'pdf', 'svg']:
            fig.savefig(destination / f'{filename}.{extension}', dpi=160)
        combined.savefig(fig)
        plt.close(fig)
    combined.close()
    if hashlib.sha256((result / 'endpoints.csv').read_bytes()).hexdigest() != endpoint_hash:
        raise RuntimeError('Endpoint table changed while rendering')
    provenance = dict(rendered_at=datetime.now(timezone.utc).isoformat(),
                      data_generated_at=status['generated_at'],
                      endpoint_table_sha256=endpoint_hash,
                      renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      presentation=f'Seven model figures: four metric rows by {len(conditions)} condition columns. Identical combination flags across metrics; concurrent flags offset within cells.',
                      conditions=conditions, reference_lengths=manifest.get('reference_lengths'),
                      artifact_hashes={f'{name}.{ext}': hashlib.sha256((destination / f'{name}.{ext}').read_bytes()).hexdigest()
                                       for name in ['model_' + m for m in MODELS] for ext in ['png', 'pdf', 'svg']})
    provenance['artifact_hashes']['polarity_establishment_by_model.pdf'] = hashlib.sha256((destination / 'polarity_establishment_by_model.pdf').read_bytes()).hexdigest()
    (destination / 'rendering_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    import shutil
    if destination != result:
        for name in ['endpoints.csv','status.json']:
            shutil.copy2(result/name,destination/name)
        shutil.copy2(out/'manifest.json',destination/'manifest.json')
    navigation = ' · '.join(f'<a href="#{model}">{label}</a>' for model, label in zip(MODELS, LABELS))
    sections = ''.join(f'<section id="{model}"><h2>{label}</h2><p><a href="model_{model}.pdf">PDF</a> · '
                       f'<a href="model_{model}.svg">SVG</a> · <a href="model_{model}.png">PNG</a></p>'
                       f'<img src="model_{model}.png" alt="{label}: all four polarity establishment metrics across {len(conditions)} conditions"></section>'
                       for model, label in zip(MODELS, LABELS))
    html = f'''<!doctype html><html><head><meta charset="utf-8"><title>Polarity establishment by model</title>
    <style>body{{font-family:sans-serif;max-width:1550px;margin:30px auto;padding:0 20px}}img{{width:100%;height:auto}}
    section{{margin-top:45px}}nav{{line-height:2}}a{{color:#185a99}}</style></head><body>
    <h1>Polarity establishment by model</h1><p>Data generated {status['generated_at']}; producer {status['producer']}.
    Report rendered {provenance['rendered_at']}.</p><p>{status['completed_cases']} / {status['cases']} combinations completed;
    {status['execution_errors']} execution errors.</p><p>Each model groups cue duration, total establishment time,
    final polarity and apparent endpoint count in four rows, with {len(conditions)} conditions in columns. Symbols are identical across metrics.</p>
    <p>Columns: {', '.join('20% of cell length' if c['width_mode']=='relative' else f"fixed {c['fraction']:.0%} of reference length" for c in conditions)}.
    Reference length {manifest.get('standardized_reference_length', 'recorded per model')}; Spring native reference {manifest.get('spring_reference_length', 'recorded per model')}.</p>
    <p><a href="polarity_establishment_by_model.pdf">Download all seven models as one PDF</a> ·
    <a href="endpoints.csv">Endpoint table</a> · <a href="status.json">Status</a> ·
    <a href="manifest.json">Protocol</a> · <a href="rendering_provenance.json">Rendering provenance</a></p>
    <p>Withdrawal removes the added spatial cue and restores basal conditions. LEGI uses uniform S=0.1, starts at A=I=0.1 and R=0.5, and returns to homogeneous response. Missing timing cells indicate no qualifying persistent polarized endpoint.</p><nav>{navigation}</nav>{sections}</body></html>'''
    temp = destination / 'index.html.tmp'
    temp.write_text(html)
    temp.replace(destination / 'index.html')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--figures',type=Path)
    args=parser.parse_args()
    import native as r
    destination=args.figures or (r.ROOT/'figures/establishment_basal')
    plot(args.out.resolve(), destination)
