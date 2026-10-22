"""Condense report pages 2–3 from audited pre-cue metrics; no integration."""
from pathlib import Path
import csv
import hashlib
import json
from datetime import datetime, timezone
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/runs/initial_polarity_length_corrected_20261008'
OUT = ROOT / 'figures/initial_polarity_length'
MODELS = ['wave_pinning', 'goryachev', 'debelly', 'holmes_model3', 'spring']
LABELS = ['Wave-pinning', 'Goryachev', 'deBelly', 'Holmes 3', 'Spring']
BLUE, ORANGE = '#0072B2', '#D55E00'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    source = DATA / 'steady_state_metrics.csv'
    provenance = json.loads((DATA / 'report_provenance.json').read_text())
    assert sha(source) == provenance['summary_csv_sha256']
    for key in ['all_checkpoint_hashes_pass', 'all_three_geometry_checkpoint_identities_pass',
                'independent_active_moments_pass', 'pole_amount_background_identity_pass']:
        assert provenance[key], key
    rows = list(csv.DictReader(source.open()))
    assert len(rows) == 97
    for r in rows:
        assert sha(ROOT / r['checkpoint']) == r['checkpoint_sha256']
        for key in ['domain_length', 'old_pole_width', 'old_pole_mean_active',
                    'old_pole_peak_active', 'old_pole_total_active', 'old_pole_excess_active']:
            assert np.isfinite(float(r[key])) and float(r[key]) > 0
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .65, 'xtick.major.width': .65,
                         'ytick.major.width': .65, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    fig, axes = plt.subplots(4, 5, figsize=(10.6, 8.0), sharex='col')
    fig.subplots_adjust(left=.15, right=.985, bottom=.16, top=.82,
                        wspace=.40, hspace=.35)
    fig.text(.15, .968, 'Longer cells carry more active material in the pre-existing pole',
             fontsize=13, weight='bold', va='top')
    fig.text(.15, .932, 'Pole growth is mainly through width, except in Goryachev, where concentration increases.',
             fontsize=9, color='#404a53', va='top')
    metrics = [
        ('old_pole_total_active', 'old_pole_excess_active', 'Active material', 'concentration × length'),
        ('old_pole_width', None, 'Pole width', 'native length'),
        ('old_pole_mean_active', 'old_pole_peak_active', 'Pole concentration', 'native concentration'),
        ('opposite_inactive', None, 'Opposite-end inactive', 'native concentration')]
    legends = [
        [Line2D([], [], color=BLUE, lw=1.5, marker='o', ms=3, label='Total'),
         Line2D([], [], color=ORANGE, lw=1.5, ls='--', label='Above background')],
        [Line2D([], [], color=BLUE, lw=1.5, label='Connected half-contrast width')],
        [Line2D([], [], color=BLUE, lw=1.5, marker='o', ms=3, label='Mean'),
         Line2D([], [], color=ORANGE, lw=1.5, ls='--', label='Peak')],
        [Line2D([], [], color=BLUE, lw=1.5, label='At the opposite end')]]
    counts = []
    for col, (model, label) in enumerate(zip(MODELS, LABELS)):
        rs = sorted([r for r in rows if r['model'] == model], key=lambda r: float(r['domain_length']))
        counts.append(len(rs))
        x = np.array([float(r['domain_length']) for r in rs])
        excluded = provenance['excluded_lengths'][model]
        ax = axes[0, col]
        ax.text(.5, 1.34, label, transform=ax.transAxes, ha='center', fontsize=10, weight='bold')
        ax.text(.5, 1.17, 'Taller pole' if model == 'goryachev' else 'Wider pole',
                transform=ax.transAxes, ha='center', fontsize=8, color='#53616b')
        for row, (key, extra, _, _) in enumerate(metrics):
            ax = axes[row, col]
            if key == 'opposite_inactive' and model == 'debelly':
                assert all(r[key] == '' for r in rs)
                ax.text(.5, .5, 'No explicit\ninactive Rac pool', transform=ax.transAxes,
                        ha='center', va='center', fontsize=8, color='#626b73', linespacing=1.5)
                ax.set_yticks([])
            else:
                y = np.array([float(r[key]) for r in rs])
                ax.plot(x, y, color=BLUE, lw=1.5, marker='o', ms=2.4, mew=0, zorder=3)
                ymax = y.max()
                if extra:
                    yy = np.array([float(r[extra]) for r in rs])
                    ax.plot(x, yy, color=ORANGE, lw=1.4, ls='--', zorder=4)
                    ymax = max(ymax, yy.max())
                ax.set_ylim(0, ymax * 1.13)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
                ax.grid(axis='y', color='#e7eaed', lw=.55, zorder=0)
            ax.tick_params(labelsize=7.5, pad=2, length=3)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
            if excluded:
                ax.scatter(excluded, [-.045] * len(excluded), marker='x', color='#8b9095',
                           s=13, lw=.75, transform=ax.get_xaxis_transform(), clip_on=False)
            lo = min([x.min()] + excluded)
            span = x.max() - lo
            ax.set_xlim(lo - .035 * span, x.max() + .035 * span)
        axes[-1, col].tick_params(axis='x', pad=6)
        axes[-1, col].set_xlabel('Resting length L0' if model == 'spring' else 'Cell length L', fontsize=8, labelpad=5)
    for row, (_, _, title, unit) in enumerate(metrics):
        box = axes[row, 0].get_position()
        fig.text(.025, box.y1, 'ABCD'[row], fontsize=12, weight='bold', va='top')
        fig.text(.025, box.y1 - .034, title.replace('Opposite-end ', 'Opposite-end\n'),
                 fontsize=8.5, weight='bold', va='top', linespacing=1.4)
        fig.text(.025, box.y1 - (.087 if row == 3 else .067), unit.replace(' × ', ' ×\n').replace('native concentration', 'native\nconcentration'),
                 fontsize=7.5, color='#59636c', va='top', linespacing=1.4)
        if row in (0, 2):
            fig.legend(handles=legends[row], loc='upper left', bbox_to_anchor=(.025, box.y1 - .113),
                       bbox_transform=fig.transFigure, frameon=False, fontsize=6.6,
                       handlelength=1.8, handletextpad=.4, borderaxespad=0,
                       labelspacing=.45)
    fig.text(.15, .093, 'Verified basal states before the reversal cue; 97 eligible preparations.  ×  No eligible polarized state.',
             fontsize=8, color='#47515a')
    fig.text(.15, .066, 'Pole: boundary-connected region above half contrast. Material and mean use the same region; background = field minimum.',
             fontsize=7.6, color='#47515a')
    fig.text(.15, .041, 'All axes use model-native units and separate scales. Spring width uses physical length; deBelly uses nominal solver coordinates.',
             fontsize=7.6, color='#47515a')
    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / 'initial_polarity_length'
    for suffix in ['png', 'pdf', 'svg']:
        path = stem.with_suffix('.' + suffix)
        # Preserve the original multipage PDF under its established name.
        if suffix != 'png':
            path = OUT / ('initial_polarity_length_manuscript.' + suffix)
        fig.savefig(path, dpi=600, facecolor='white')
    plt.close(fig)
    caption = ('Initial polarized-state burden increases with cell length through model-dependent pole geometry. '
               '(A) Integrated active material in the original pole, including background (solid blue) or after subtracting '
               'the minimum field concentration (dashed orange). (B) Width of the boundary-connected region above halfway '
               'between the field minimum and maximum. (C) Mean concentration in that region (solid blue) and peak '
               'concentration (dashed orange). (D) Inactive concentration at the opposite end; deBelly has no explicit '
               'inactive Rac field. Longer cells have wider poles in Wave-pinning, deBelly, Holmes 3 and Spring, whereas '
               'Goryachev retains a nearly constant pole width and increases pole concentration. Each point is an exact '
               'verified basal-release checkpoint shared by the three reversal tests, prepared with a nominal-amplitude '
               'left patch of 20% support. Lines connect sampled lengths; gray crosses mark ineligible preparations. '
               'All 97 eligible states are included (20, 21, 12, 23 and 21 per model, respectively). Native scales differ '
               'between models. Spring is plotted against resting length L0, with width and material evaluated on its '
               'physical domain. deBelly uses nominal chemical-solver coordinates. Integrated material is a one-dimensional '
               'concentration integral, rather than a calibrated molecule count. These starting-state measurements '
               'describe the material present before reversal; they do not by themselves establish a reversal-time mechanism.\n')
    (OUT / 'initial_polarity_length_caption.txt').write_text(caption)
    review = ROOT.parents[1] / 'current manuscript figs/initial_polarity_length.png'
    shutil.copy2(stem.with_suffix('.png'), review)
    record = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer='finished',
                  source_metrics_sha256=sha(source), source_report_provenance_sha256=sha(DATA/'report_provenance.json'),
                  renderer_sha256=sha(Path(__file__)), source_checkpoint_hashes_pass=True,
                  preparations=len(rows), models=MODELS, counts=counts, excluded_lengths=provenance['excluded_lengths'],
                  outputs={str(p.relative_to(ROOT)): sha(p) for p in [stem.with_suffix('.png'),
                      OUT/'initial_polarity_length_manuscript.pdf', OUT/'initial_polarity_length_manuscript.svg',
                      OUT/'initial_polarity_length_caption.txt']},
                  manuscript_copy=str(review), manuscript_copy_sha256=sha(review))
    (DATA/'manuscript_figure_provenance.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'png': str(stem.with_suffix('.png')), 'review_copy': str(review), 'checkpoints': len(rows)}, indent=2))

if __name__ == '__main__':
    main()
