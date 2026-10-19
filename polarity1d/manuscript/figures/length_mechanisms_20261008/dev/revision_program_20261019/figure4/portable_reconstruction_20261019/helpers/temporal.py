"""Comprehensive observed-window F4 temporal supplement, all model/size/f.

Direct read-only plotting of existing audited sufficient-statistic aggregates.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS = ROOT / 'polarity1d/manuscript'
FIG = MS / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
SOURCE = FIG / 'broad_strength_campaign_20261019/review_complete_scientific/supporting/data/memory_all_conditions_windows_lags.csv'
OUT = FIG / 'postburn_signature_candidate_20261019'
MODELS = dict(wave_pinning='Wave-pinning', goryachev='Goryachev', holmes_model3='Holmes 3', spring_ring='Spring ring')
GRID = (.0001, .0003, .001, .003, .01, .03, .1, .3, 1.)
COLORS = ('#28206f', '#324991', '#4055a8', '#0d8a93', '#66923d', '#d2882c', '#ba445b', '#8c467c', '#5f3579')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=SOURCE); p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(); out = args.out
    for sub in ('vectors', 'data', 'metadata', 'viewers'):
        (out / 'supporting' / sub).mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.source); data = df[df.kind == 'noise'].copy()
    assert len(data) == 1944 and set(data.fraction) == set(GRID)
    data.to_csv(out / 'supporting/data/full_temporal_supplement_points.csv', index=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    panels = []; artists = []
    with PdfPages(out / 'supporting/vectors/all_models_sizes_strengths_temporal_windows.pdf') as pdf:
        for window, window_label, limits in (
            ('full', 'Whole observed run', 'chemical0–32000; Spring0–3200'),
            ('first_half', 'First half', 'chemical0–16000; Spring0–1600'),
            ('last_half', 'Second half after discarding the first half', 'chemical16000–32000; Spring1600–3200')):
            for metric, label in (
                ('direction_cosine_reported', 'Conditional direction memory'),
                ('direction_coverage', 'Usable direction-pair coverage'),
                ('raw_vector_correlation', 'Uncentered first-moment correlation')):
                fig, axes = plt.subplots(4, 3, figsize=(12.8, 13.4), sharey=True)
                fig.subplots_adjust(left=.075, right=.98, top=.83, bottom=.115, hspace=.66, wspace=.2)
                fig.suptitle(label + ': all nine strengths', x=.075, y=.973, ha='left', fontsize=18)
                fig.text(.075, .94, window_label + ' · ' + limits, fontsize=11)
                handles = [Line2D([], [], color=c, lw=1.8, label=f'f={f:g}') for f, c in zip(GRID, COLORS)]
                fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.068, .922),
                           ncol=5, frameon=False, fontsize=9.5, handlelength=2.4, columnspacing=1.7)
                fig.text(.075, .852, 'Filled squares: n=12. Open circles: n=2 exploratory levels. The table retains exact n, pairs and contributing seeds.', fontsize=9.5)
                for row, (model, title) in enumerate(MODELS.items()):
                    for column, size in enumerate((1., 2., 3.) if model == 'spring_ring' else (20., 40., 60.)):
                        ax = axes[row, column]
                        for f, color in zip(GRID, COLORS):
                            d = data[(data.model == model) & (data['size'] == size)
                                & (data.window == window) & (data.fraction == f)].sort_values('lag')
                            assert len(d) == 6
                            n = int(d.cases.iloc[0])
                            line, = ax.plot(d.lag, d[metric], color=color, marker='s' if n == 12 else 'o',
                                mfc=color if n == 12 else 'white', ms=3.5, lw=1.2)
                            artists.append(dict(window=window, metric=metric, model=model, size=size, fraction=f,
                                seeds=n, lag=d.lag.tolist(), plotted_y=[float(v) if np.isfinite(v) else None for v in d[metric]], color=color))
                        ax.set_xscale('log'); ax.set_ylim(-.04, 1.05) if metric == 'direction_coverage' else ax.set_ylim(-1.05, 1.05)
                        ax.spines[['top', 'right']].set_visible(False)
                        ax.grid(axis='y', color='#e5e9eb', lw=.7); ax.tick_params(labelsize=8.5)
                        ax.set_title(f'{title} · {"C₀" if model == "spring_ring" else "C"}={size:g}', loc='left', fontsize=10.5)
                        ax.set_xlabel('Lag (native time units)', fontsize=9)
                        if column == 0:
                            ax.set_ylabel('Usable-pair fraction' if metric == 'direction_coverage' else 'Correlation', fontsize=9)
                fig.text(.075, .074, 'Direction requires both endpoint R≥0.1 and reports groups only with≥100 eligible pairs and≥1% coverage. Read direction with its coverage page.', fontsize=9.2)
                fig.text(.075, .048, 'Uncentered vector correlation uses all raw first-moment pairs and lag-energy normalization; it also reflects magnitude/forced adaptation. It is descriptive here.', fontsize=9.2)
                fig.text(.075, .022, 'No stationarity or migration-velocity law is assumed. Chemical and Spring clocks/noise laws differ. All six lag checkpoints are retained; lines guide the eye.', fontsize=9.2)
                stem = f'full_temporal_{window}_{metric}'
                fig.savefig(out / 'supporting/viewers' / (stem + '.png'), dpi=140)
                pdf.savefig(fig); plt.close(fig)
                panels.append(stem)
    receipt = dict(generated_at=datetime.now(timezone.utc).isoformat(), source=str(Path(__file__).relative_to(ROOT)),
        source_sha256=sha(__file__), aggregate_source=portable_path(args.source), aggregate_sha256=sha(args.source),
        observed_noise_groups=108, points=1944, metrics=3, pdf_pages=9, panels=panels,
        native_integrations=0, noise_only=True, deterministic_controls_retained_in_source=True,
        artist_count=len(artists), artists=artists)
    (out / 'supporting/metadata/FULL_TEMPORAL_RENDER_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k != 'artists'}, indent=2))


def portable_path(path):
    path=Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


if __name__ == '__main__':
    main()
