"""The single final collection of four polarity-establishment heatmap sets."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, MaxNLocator
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


def plot(out):
    result = out / 'results'
    status = json.loads((result / 'status.json').read_text())
    if not status['complete']:
        raise RuntimeError('Final collection requires complete execution')
    with (result / 'endpoints.csv').open() as f:
        rows = list(csv.DictReader(f))
    for filename, field, title, unit, log in SETS:
        fig, axes = plt.subplots(7, 2, figsize=(15, 22), constrained_layout=True)
        fig.get_layout_engine().set(rect=(0., .09, 1., .87))
        for i, (model, label) in enumerate(zip(MODELS, LABELS)):
            model_rows = [r for r in rows if r['model'] == model]
            vals = [float(r[field]) for r in model_rows if r[field] and (not log or float(r[field]) > 0)]
            if field == 'selected_final_polarity':
                norm = Normalize(0, 1)
            elif log:
                low, high = (min(vals), max(vals)) if vals else (.01, 1.)
                norm = LogNorm(low, high if high > low else low * 10)
            else:
                norm = Normalize(1, max(2, max(vals)) if vals else 2)
            for j, fraction in enumerate([.1, .2]):
                ax = axes[i, j]
                chosen = [r for r in model_rows if float(r['width_fraction_at_reference']) == fraction]
                lengths = sorted({float(r['length']) for r in chosen})
                amplitudes = sorted({float(r['amplitude']) for r in chosen})
                matrix = np.full((len(amplitudes), len(lengths)), np.nan)
                for row in chosen:
                    x, y = float(row['length']), float(row['amplitude'])
                    if row[field] and (not log or float(row[field]) > 0):
                        matrix[amplitudes.index(y), lengths.index(x)] = float(row[field])
                cmap = plt.get_cmap('viridis').copy()
                cmap.set_bad('#dedede')
                mesh = ax.pcolormesh(edges(lengths), edges(amplitudes, True),
                                     np.ma.masked_invalid(matrix), cmap=cmap, norm=norm, shading='flat')
                ax.set_yscale('log')
                ax.set_yticks(amplitudes, [f'{v:g}' for v in amplitudes])
                ax.yaxis.set_minor_formatter(NullFormatter())
                ax.set_xticks(lengths[::max(1, len(lengths) // 8)])
                ax.set_xlabel('L (native length)' + ('; not μm' if model == 'spring' else ''))
                ax.set_ylabel('A (native peak amplitude)')
                ax.set_title(f"{label} · fixed width {float(chosen[0]['width']):g} ({fraction:.0%} of reference L)")
                for row in chosen:
                    x, y = float(row['length']), float(row['amplitude'])
                    marker = {'unpolarized': 'x', 'unsettled': '+', 'no_qualifying_endpoint': 's'}.get(row['status'])
                    if marker:
                        ax.plot(x, y, marker=marker, color='black', ms=4, mfc='none', linestyle='none')
                    if row['upper_optimum_unresolved'] == 'True' or row['unsettled_trials'] != '0':
                        ax.plot(x, y, marker='^', color='white', mec='black', ms=4, linestyle='none')
                    if log and row['duration_boundary_unresolved'] == 'True':
                        ax.plot(x, y, marker='o', color='black', mfc='none', ms=5, linestyle='none')
                    if field == 'cue_duration' and row[field] and float(row[field]) == 0:
                        ax.plot(x, y, marker='D', color='black', mfc='white', ms=5, linestyle='none')
                cb = fig.colorbar(mesh, ax=ax, label=unit, fraction=.035, pad=.02)
                if field == 'apparent_endpoint_count':
                    cb.locator = MaxNLocator(integer=True)
                    cb.update_ticks()
        fig.suptitle(title + '\nNative rest → fixed left patch → verified zero-cue release', fontsize=17, y=.995)
        handles = [Line2D([], [], marker=m, color='black', mfc='none', linestyle='none', label=l)
                   for m, l in [('x', 'Unpolarized'), ('+', 'No verified endpoint'),
                                ('s', 'No qualifying polarized trial'), ('^', 'Upper optimum or some releases unresolved'),
                                ('o', 'Duration boundary >1% bracket'), ('D', 'No polarizing cue required (D=0)')]]
        fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .047), ncol=3, fontsize=9)
        fig.text(.5, .01, 'Strongest observed cue-aligned endpoint, then shortest tested duration within 1% of its strength.\n'
                      'Counts: approximate complete-link profile clusters; zero-duration control included; no extra search for counts.\n'
                      'Colors share a scale between widths within each model; native time units differ between models.\n'
                      'Generated ' + status['generated_at'] + '; producer finished.', fontsize=9, ha='center', va='bottom')
        for extension in ['png', 'pdf', 'svg']:
            fig.savefig(result / f'{filename}.{extension}', dpi=160)
        plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    plot(parser.parse_args().out.resolve())
