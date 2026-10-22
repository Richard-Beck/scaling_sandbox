"""Figure 12 presentation revision with distinct length and process palettes.

Read the validated 13 October plotted tables verbatim. No native solves,
new scientific calculations, parameter changes, or narrative changes occur.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shutil

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
PRIOR = REPORT / 'revision_20261013/conclusion'
OUT = REPORT / 'revision_20261014/conclusion'
STEM = 'Figure_12_long_cell_reversal'
LENGTH_COLOR = {31: '#882255', 37: '#88CCEE', 39: '#999933'}
EXAMPLES = [31, 37, 39]
TERMS = ['velocity_gradient', 'profile_motion', 'binding_exchange',
         'uniform_replacement', 'bound_diffusion']
TERM_LABEL = {'velocity_gradient': 'Flow gradient', 'profile_motion': 'Profile advection',
              'binding_exchange': 'Net binding', 'uniform_replacement': 'Uniform replacement',
              'bound_diffusion': 'Bound diffusion'}
TERM_COLOR = {'velocity_gradient': '#D55E00', 'profile_motion': '#0072B2',
              'binding_exchange': '#009E73', 'uniform_replacement': '#CC79A7',
              'bound_diffusion': '#7F7F7F'}
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(source(path).read_text())


def render(traces, native, windows):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.0))
    fig.subplots_adjust(left=.085, right=.985, top=.955, bottom=.14, hspace=.42, wspace=.28)
    for ax, letter in zip(axes.flat, 'ABCD'):
        ax.text(-.13, 1.025, letter, transform=ax.transAxes, fontsize=15, weight='bold')
        ax.grid(axis='y', color='#DFE3E5', lw=.65, alpha=.7)
        ax.tick_params(labelsize=10)
    for L in EXAMPLES:
        q, rate, c = traces[traces.length_um == L], native[native.length_um == L], LENGTH_COLOR[L]
        for key, ls in [('right_bound_outflux', '-'), ('left_bound_outflux', '--')]:
            axes[0, 0].plot(rate.time_end_s, rate[key] * 1e3, color=c, lw=1.9, ls=ls)
        for ax, key in [(axes[0, 1], 'old_endpoint_Rac'), (axes[1, 0], 'front_bound_MCA')]:
            ax.plot(q.time_s, q[key], color=c, lw=2.05)
            ax.plot(q.time_s.iloc[-1], q[key].iloc[-1], 'o', ms=4.8, mfc='white', mec=c, mew=1.3)
    handles = [Line2D([], [], color=LENGTH_COLOR[L], lw=2, label=f'{L} µm') for L in EXAMPLES]
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        legend = ax.legend(handles=handles, ncol=3, loc='upper right', frameon=False,
            fontsize=9.5, columnspacing=1.1, handlelength=1.5, title='Cell length', title_fontsize=9.5)
        if ax is axes[0, 0]:
            ax.add_artist(legend)
            ax.legend(handles=[Line2D([], [], color='#555555', ls=ls, lw=1.8, label=label)
                for ls, label in [('--', 'Old left end'), ('-', 'New right end')]],
                loc='upper left', bbox_to_anchor=(0., .82), ncol=2, frameon=False, fontsize=9.5)
        ax.set(xlim=(0, 835), xticks=[0, 200, 400, 600, 800], xlabel='Time under relative-20% cue (s)')
    axes[0, 0].set(ylabel='Bound-MCA loss at each end\n(10⁻³ conc. µm/s)', ylim=(-.15, 16.))
    axes[0, 1].set(ylabel='Rac at the old left endpoint', ylim=(.69, .96))
    ax = axes[1, 0]
    ax.set(ylabel='Mean bound MCA in rightmost 20%', ylim=(.038, .106))
    ax.axvspan(60, 180, color='#CCCCCC', alpha=.27, lw=0)
    ax.axvspan(507, 700, facecolor='#E7E7E7', edgecolor='#BBBBBB', alpha=.35, hatch='///', lw=0)
    q = traces[traces.length_um == 39]
    ax.plot(507, q.front_bound_MCA.min(), 'o', ms=5, color=LENGTH_COLOR[39])
    ax.text(512, .0414, '507 s', color=LENGTH_COLOR[39], fontsize=9.5, va='top')
    ax = axes[1, 1]
    xpos = np.array([0., 1., 2., 3.45])
    positive, negative = np.zeros(4), np.zeros(4)
    for term in TERMS:
        values = windows[term + '_conc_s'].to_numpy() * 1e4
        ax.bar(xpos, values, width=.62, bottom=np.where(values >= 0, positive, negative),
            color=TERM_COLOR[term], edgecolor='white', linewidth=.35, label=TERM_LABEL[term], zorder=3)
        positive += np.maximum(values, 0)
        negative += np.minimum(values, 0)
    ax.plot(xpos, windows.direct_observed_rate_conc_s.to_numpy() * 1e4,
        'D', color='#202020', ms=5.1, label='Net change', zorder=5)
    ax.axhline(0, color='#565E62', lw=.85)
    ax.set(xticks=xpos, xticklabels=['31 µm', '37 µm', '39 µm', '39 µm'],
        ylim=(-8.1, 6.8), ylabel='Bound-MCA change in rightmost 20%\n(10⁻⁴ conc./s; positive = accumulation)')
    ax.text(1, -.16, 'Common 60–180 s', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=9.6)
    ax.text(3.45, -.16, 'Rebound\n507–700 s', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=9.6)
    ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(.5, 1.28), frameon=False,
        fontsize=8.8, columnspacing=1., handlelength=1.3, handletextpad=.4, title='Bound-MCA contributions', title_fontsize=9)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / (STEM + '.' + ext), dpi=220)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source(Path(__file__))
    prior = read(PRIOR / (STEM + '.manifest.json'))
    assert prior['passed'] and prior['producer'] == 'finished'
    for name, digest in prior['artifacts'].items():
        assert sha(source(PRIOR / name)) == digest, name
    for name, digest in prior['sources'].items():
        assert sha(source(ROOT / name)) == digest, name
    names = ['windows', 'native_interval_budgets', 'all_length_persistence', 'traces']
    unchanged = {}
    for name in names:
        src, dst = PRIOR / f'{STEM}_{name}.csv', OUT / f'{STEM}_{name}.csv'
        shutil.copy2(src, dst)
        assert sha(src) == sha(dst)
        unchanged[name] = sha(dst)
    traces = pd.read_csv(OUT / (STEM + '_traces.csv'))
    native = pd.read_csv(OUT / (STEM + '_native_interval_budgets.csv'))
    windows = pd.read_csv(OUT / (STEM + '_windows.csv'))
    render(traces, native, windows)
    chapter = read(PRIOR / 'chapter.json')
    previous_paragraphs = chapter['paragraphs'].copy()
    chapter.update(image=f'revision_20261014/conclusion/{STEM}.png',
        manifest=f'revision_20261014/conclusion/{STEM}.manifest.json')
    chapter['caption'] = chapter['caption'].replace('Shading identifies the common 60–180-s comparison and descriptive 39-µm rebound interval, 507–700 s.',
        'Plain gray shading identifies the common 60–180-s comparison; hatched gray shading identifies the descriptive 39-µm rebound interval, 507–700 s. Local keys distinguish cell-length colors in A–C from process colors in D.')
    chapter['sources'] = [[url.replace('revision_20261013/conclusion/', 'revision_20261014/conclusion/'), label]
                          for url, label in chapter['sources']]
    assert chapter['paragraphs'] == previous_paragraphs
    for name in ['METHODS.md', STEM + '_prose.md']:
        shutil.copy2(PRIOR / name, OUT / name)
    (OUT / (STEM + '_caption.txt')).write_text(chapter['caption'] + '\n')
    (OUT / 'chapter.json').write_text(json.dumps(chapter, indent=2, ensure_ascii=False) + '\n')
    audit = dict(passed=True, generated_at=datetime.now(timezone.utc).isoformat(), producer='finished',
        presentation_only=True, new_integrations=0, scientific_paragraphs_unchanged=True,
        data_byte_identity=unchanged, source_hashes_passed=True, length_colors=LENGTH_COLOR,
        process_colors=TERM_COLOR, process_palette_matches_figure9=True,
        local_keys=['A: length and endpoint style', 'B: length', 'C: length', 'D: bound-MCA contribution'],
        categorical_length_colors_distinct_from_process_colors=True,
        inherited_numerical_audit=prior['audit'])
    (OUT / (STEM + '_audit.json')).write_text(json.dumps(audit, indent=2) + '\n')
    artifacts = {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and not p.name.endswith('.manifest.json')}
    manifest = dict(generated_at=audit['generated_at'], passed=True, producer='finished',
        renderer=str(Path(__file__).relative_to(ROOT)), renderer_sha256=sha(Path(__file__)),
        sources=SOURCES, artifacts=artifacts, audit=audit,
        reproduction='/home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/code/mechanisms/length_revision_20261014_conclusion.py')
    (OUT / (STEM + '.manifest.json')).write_text(json.dumps(manifest, indent=2) + '\n')
    print(OUT / (STEM + '.png'))


if __name__ == '__main__':
    main()
