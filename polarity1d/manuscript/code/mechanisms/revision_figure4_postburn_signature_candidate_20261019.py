"""Timestamped F4 candidate from existing post-burn-in observations only.

No native integration, event classifier or change to accepted report/main art.
The first export consumes the independently audited last-half memory table.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import shutil

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[4]
MS = ROOT / 'polarity1d/manuscript'
RUN = MS / 'data/runs/revision_noise_strength_broad_20261019'
FIG = MS / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
COMPLETE = FIG / 'broad_strength_campaign_20261019/review_complete_scientific'
DEFAULT = FIG / 'postburn_signature_candidate_20261019'
MODELS = dict(wave_pinning='Wave-pinning', goryachev='Goryachev',
              holmes_model3='Holmes 3', spring_ring='Spring ring')
STRENGTHS = (.001, .003, .01, .03, .1)
COLORS = ('#4055a8', '#0d8a93', '#66923d', '#d2882c', '#ba445b')


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def export(fig, out, stem):
    fig.savefig(out / (stem + '.png'), dpi=190)
    for extension in ('pdf', 'svg'):
        fig.savefig(out / 'supporting/vectors' / (stem + '.' + extension))


def memory(out):
    source = COMPLETE / 'supporting/data/memory_all_conditions_windows_lags.csv'
    assert sha(source) == 'f9c0bc752b303a4f25ab445d5671d4183c6d8d6c713db73e42506b26ee2021ca'
    df = pd.read_csv(source)
    selected = df[(df.kind == 'noise') & (df.window == 'last_half')
        & df.fraction.isin(STRENGTHS)
        & (((df.model == 'spring_ring') & (df.lag == 800))
           | ((df.model != 'spring_ring') & (df.lag == 8000)))].copy()
    assert len(selected) == 60
    manifest = json.loads((RUN / 'manifest.json').read_text())
    physical = {}
    for b in manifest['basals']:
        with np.load(RUN / b['input']) as z:
            # Scalar native C is the final state entry of the Spring basal.
            C = float(z['state'][-1]) if b['model'] == 'spring_ring' else float(b['size'])
        physical[(b['model'], float(b['size']))] = C
    selected['basal_circumference'] = [physical[(r.model, r['size'])] for _, r in selected.iterrows()]
    selected['analysis_first_time'] = np.where(selected.model == 'spring_ring', 1600., 16000.)
    selected['analysis_last_time'] = np.where(selected.model == 'spring_ring', 3200., 32000.)
    selected.to_csv(out / 'supporting/data/postburn_memory_coverage_points.csv', index=False)
    shutil.copy2(source, out / 'supporting/data/source_memory_all_windows.csv')
    fig, axes = plt.subplots(2, 4, figsize=(13.2, 7.1), sharey='row')
    fig.subplots_adjust(left=.072, right=.98, bottom=.23, top=.72, wspace=.24, hspace=.25)
    fig.suptitle('Second-half directional memory across size and input regimes',
                 x=.072, y=.97, ha='left', fontsize=17)
    fig.text(.072, .928, 'Discard first half: chemicals t=16000–32000; Spring t=1600–3200. These finite windows are not assumed stationary.', fontsize=10.5)
    handles = [Line2D([], [], color=c, lw=1.8, label=f'f={f:g}') for f, c in zip(STRENGTHS, COLORS)]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.066, .906),
               ncol=5, frameon=False, fontsize=10, handlelength=2.8, columnspacing=2.1)
    fig.text(.072, .815, '■ original strengths: n=12 seeds     ○ new strengths: n=2 seeds     Thin whiskers: supported-seed range, not a confidence interval.', fontsize=9.6)
    for column, (model, label) in enumerate(MODELS.items()):
        top, bottom = axes[:, column]
        for f, color in zip(STRENGTHS, COLORS):
            d = selected[(selected.model == model) & (selected.fraction == f)].sort_values('size')
            x = d.basal_circumference.to_numpy()
            y = d.direction_cosine_reported.to_numpy()
            n = int(d.cases.iloc[0])
            marker = 's' if n == 12 else 'o'
            face = color if n == 12 else 'white'
            top.plot(x, y, color=color, marker=marker, mfc=face, mec=color, ms=5, lw=1.6)
            top.vlines(x, d.supported_seed_cosine_min, d.supported_seed_cosine_max,
                       color=color, alpha=.32, lw=.8)
            missing = ~np.isfinite(y)
            top.plot(x[missing], np.full(missing.sum(), -.31), linestyle='none',
                     marker='x', color=color, ms=7, clip_on=False)
            bottom.plot(x, d.direction_coverage, color=color, marker=marker,
                        mfc=face, mec=color, ms=5, lw=1.6)
        for ax in (top, bottom):
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(axis='y', color='#e5e9eb', lw=.7)
            ax.tick_params(labelsize=9)
            d = selected[selected.model == model].drop_duplicates('size').sort_values('size')
            ticks = d.basal_circumference.to_numpy()
            ax.set_xticks(ticks, [f'{v:.3f}' if model == 'spring_ring' else f'{v:g}' for v in ticks])
            ax.set_xlim(ticks[0] - .09 * np.ptp(ticks), ticks[-1] + .09 * np.ptp(ticks))
        top.set_title(f'{chr(65 + column)}  {label}\nlag={800 if model == "spring_ring" else 8000} native units', loc='left', fontsize=10.5)
        top.axhline(0, color='#89969b', lw=.65)
        bottom.set_xlabel('Basal physical circumference' if model == 'spring_ring' else 'Circumference', fontsize=10)
        top.set_ylim(-1.05, 1.065)
        bottom.set_ylim(-.035, 1.05)
        bottom.set_yticks([0, .5, 1])
    axes[0, 0].set_ylabel('Conditional direction cosine', fontsize=10)
    axes[1, 0].set_ylabel('Usable-pair fraction', fontsize=10)
    fig.text(.072, .132, 'Both endpoint R≥0.1. Direction is reported only with ≥100 eligible pairs and ≥1% coverage; × marks unavailable, while coverage remains shown.', fontsize=9.5)
    fig.text(.072, .101, 'Heading is the activity-vector direction, not cell movement. Whole positive input mean and fluctuation amplitude scale together; native laws/clocks differ by model.', fontsize=9.5)
    fig.text(.072, .07, 'Spring uses actual basal circumference for this size axis and actual C(t) in physical activity moments. Connecting segments guide the eye between three sizes.', fontsize=9.5)
    export(fig, out, 'postburn_memory_coverage')
    plt.close(fig)
    curves = df[(df.kind == 'noise') & (df.window == 'last_half')
        & df.fraction.isin(STRENGTHS)].copy()
    curves.to_csv(out / 'supporting/data/postburn_temporal_curve_points.csv', index=False)
    with PdfPages(out / 'supporting/vectors/postburn_temporal_curves_with_coverage.pdf') as pdf:
        for metric, stem, title, ylabel, ylim in (
            ('direction_cosine_reported', 'postburn_temporal_memory',
             'Second-half temporal direction memory', 'Conditional direction cosine', (-.4, 1.05)),
            ('direction_coverage', 'postburn_temporal_coverage',
             'Usable-pair coverage of the same temporal memory estimates', 'Usable-pair fraction', (-.035, 1.05))):
            fig, axes = plt.subplots(4, 3, figsize=(12.8, 13.3), sharey=True)
            fig.subplots_adjust(left=.075, right=.98, bottom=.115, top=.865, hspace=.65, wspace=.2)
            fig.suptitle(title, x=.075, y=.975, ha='left', fontsize=18)
            fig.text(.075, .946, 'Six stored native lag checkpoints; discard first half. Colors/seed counts match the circumference summary.', fontsize=10.5)
            fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.068, .933),
                       ncol=5, frameon=False, fontsize=10, handlelength=2.8, columnspacing=2.1)
            for row, (model, label) in enumerate(MODELS.items()):
                for column, size in enumerate((1., 2., 3.) if model == 'spring_ring' else (20., 40., 60.)):
                    ax = axes[row, column]
                    for f, color in zip(STRENGTHS, COLORS):
                        d = curves[(curves.model == model) & (curves['size'] == size)
                            & (curves.fraction == f)].sort_values('lag')
                        n = int(d.cases.iloc[0])
                        ax.plot(d.lag, d[metric], color=color, marker='s' if n == 12 else 'o',
                                mfc=color if n == 12 else 'white', ms=3.8, lw=1.35)
                    ax.set_xscale('log'); ax.set_ylim(*ylim)
                    ax.spines[['top', 'right']].set_visible(False)
                    ax.grid(axis='y', color='#e5e9eb', lw=.7)
                    ax.tick_params(labelsize=9)
                    ax.set_title(f'{label} · {"C₀" if model == "spring_ring" else "C"}={size:g}',
                                 loc='left', fontsize=10.5)
                    ax.set_xlabel('Lag (native time units)', fontsize=9.5)
                    if column == 0:
                        ax.set_ylabel(ylabel, fontsize=9.5)
            fig.text(.075, .065, 'Model clocks and forcing laws differ. These finite-window curves describe the activity-vector heading, not a cell-migration autocorrelation.', fontsize=9.3)
            fig.text(.075, .04, 'Use with its paired coverage sheet: weak-resultant endpoints are excluded from direction; sparse estimates remain unavailable. Connecting segments guide the eye.', fontsize=9.3)
            export(fig, out, stem); pdf.savefig(fig); plt.close(fig)
    return dict(memory_source=str(source.relative_to(ROOT)), memory_source_sha256=sha(source),
        points=60, temporal_curve_points=len(curves), primary_window='last half after discarding first half', stationarity_assumed=False,
        analysis_windows=dict(chemical=[16000, 32000], spring=[1600, 3200]),
        chemical_lag=8000, spring_lag=800, strengths=list(STRENGTHS),
        uncertainty='Observed supported-seed ranges only, not CIs; n=2 is exploratory',
        unavailable_points=int(selected.direction_cosine_reported.isna().sum()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=DEFAULT)
    a = p.parse_args()
    for sub in ('data', 'metadata', 'vectors'):
        (a.out / 'supporting' / sub).mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    result = memory(a.out)
    result.update(generated_at=now(), source=str(Path(__file__).relative_to(ROOT)), source_sha256=sha(__file__),
        native_integrations=0, original_report_changed=False, main_figure_replaced=False)
    result['outputs'] = {str(path.relative_to(a.out)): sha(path) for path in sorted(a.out.rglob('*'))
                         if path.is_file() and path.name != 'MEMORY_RENDER_RECEIPT.json'}
    write(a.out / 'supporting/metadata/MEMORY_RENDER_RECEIPT.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
