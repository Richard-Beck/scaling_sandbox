"""Standalone Figure 3B layout proposal using the current manuscript fit."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
CAMPAIGN = ROOT / 'wgd_langevin/refined_ou_20261003'
FIGURE_PROVENANCE = CAMPAIGN / 'figures/main_03_migration_OU_provenance.json'
original = json.loads(FIGURE_PROVENANCE.read_text())
fit_path = CAMPAIGN / 'results' / original['selected_primary_fit']
assert hashlib.sha256(fit_path.read_bytes()).hexdigest() == original['source_sha256'][str(fit_path.relative_to(ROOT))]
fit = json.loads(fit_path.read_text())
sys.path.insert(0, str(CAMPAIGN / 'scripts'))
from figure_tracks import msd_expected

colors = ['#168b87', '#d87529']
navy = '#19394d'
curves = []
for group in (0, 1):
    rows = sorted((r for r in fit['rows'] if r['study'] == 'Xuan2018'
                   and r['metric'] == 'msd_curve' and int(r['group']) == group),
                  key=lambda r: r['lag_h'])
    curves.append((np.array([r['lag_h'] for r in rows]),
                   np.array([r['obs'] for r in rows])))
pars = np.exp(np.asarray(fit['z'][:20]).reshape(5, 2, 2))[0]
lower = max(t[0] for t, y in curves)
upper = min(t[-1] for t, y in curves)
native = np.unique(np.concatenate([t[(t >= lower) & (t <= upper)]
                                  for t, y in curves]))
lag = np.geomspace(lower, upper, 600)
observed = np.array([np.interp(lag, t, y) for t, y in curves])
predicted = np.array([msd_expected(lag, *p) for p in pars])
ratio_data = observed[1] / observed[0]
ratio_fit = predicted[1] / predicted[0]

def data_difference(t):
    return np.interp(t, *curves[1]) - np.interp(t, *curves[0])

def fit_difference(t):
    return float(msd_expected(t, *pars[1]) - msd_expected(t, *pars[0]))

def crossings(grid, fn):
    values = np.array([fn(t) for t in grid])
    return [float(brentq(fn, grid[i], grid[i+1]))
            for i in np.flatnonzero(values[:-1] * values[1:] < 0)]

data_crossings = crossings(native, data_difference)
fit_crossings = crossings(lag, fit_difference)
assert np.allclose(data_crossings, original['population_curve_crossovers_h'], atol=1e-9)
assert np.allclose(fit_crossings, original['fitted_ou_crossovers_h'], atol=1e-9)
for group in (0, 1):
    rows = [r for r in fit['rows'] if r['study'] == 'Xuan2018'
            and r['metric'] == 'msd_curve' and int(r['group']) == group]
    assert np.allclose(msd_expected(np.array([r['lag_h'] for r in rows]), *pars[group]),
                       [r['predicted'] for r in rows], rtol=1e-10)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig = plt.figure(figsize=(7.5, 6.2))
gs = fig.add_gridspec(2, 1, height_ratios=[2.4, 1.25], hspace=.12,
                      left=.14, right=.79, bottom=.15, top=.90)
ax = fig.add_subplot(gs[0])
ratio = fig.add_subplot(gs[1], sharex=ax)
for group, (color, (t, y)) in enumerate(zip(colors, curves)):
    # Thin full traces preserve the digitized shape; symbols are display samples.
    ax.plot(t, y, color=color, lw=.8, alpha=.45)
    ax.plot(lag, predicted[group], color=color, lw=1.8, ls='--')
    indices = np.unique(np.r_[np.arange(0, len(t), 12), len(t)-1])
    ax.plot(t[indices], y[indices], linestyle='none', marker='o', markersize=4.5,
            markerfacecolor='white', markeredgecolor=color, markeredgewidth=1.15)
    label = ['non-PGCC\n(WGD−)', 'PGCC\n(WGD+)'][group]
    ax.text(1.025, y[-1], label, transform=ax.get_yaxis_transform(),
            color=color, va='center', fontsize=10, fontweight='bold')
ax.set(xscale='log', yscale='log', ylabel='MSD (µm²)', ylim=(5, 11000))
ax.set_xlim(lower * .94, upper * 1.08)
ax.tick_params(axis='x', which='both', labelbottom=False)
ax.legend(handles=[Line2D([], [], color='#536470', marker='o', markerfacecolor='white',
                          lw=.8, markersize=4.5, label='Digitized population curve'),
                   Line2D([], [], color='#536470', lw=1.8, ls='--', label='OU fit')],
          loc='upper left', frameon=False, fontsize=8.5)

# Offset reference segments show slopes without asserting an observed power law.
guide_lag = np.array([.22, .44])
for slope in (1, 2):
    guide_y = 400 * (guide_lag / guide_lag[0])**slope
    ax.plot(guide_lag, guide_y, color='#89939b', lw=.9)
    ax.text(.46, guide_y[-1], f'slope {slope}', color='#78828a',
            va='center', fontsize=8)
ax.grid(axis='y', which='major', color='#e9edf0', lw=.7)

ratio.plot(lag, ratio_data, color=navy, lw=1.9)
ratio.plot(lag, ratio_fit, color='#74818b', lw=1.7, ls='--')
ratio.axhline(1, color='#475963', lw=.9, zorder=0)
for cross, color, marker, label, position in (
        (data_crossings[0], navy, 'o', 'Digitized ≈2.3 h', (3.8, .72)),
        (fit_crossings[0], '#74818b', 'D', 'OU fit ≈1.1 h', (.42, 1.88))):
    ratio.scatter(cross, 1, marker=marker, s=27, facecolor='white', edgecolor=color,
                  linewidth=1.1, zorder=5)
    ratio.annotate(label, xy=(cross, 1), xytext=position, fontsize=8.5,
                   color=color, va='center',
                   arrowprops={'arrowstyle': '-', 'color': color, 'lw': .8})
ratio.text(1.025, 1, 'Equal MSD', transform=ratio.get_yaxis_transform(),
           va='center', color='#475963', fontsize=9)
ratio.set(ylabel='MSD ratio\nWGD+ / WGD−', xlabel='Lag time (h; logarithmic scale)',
          ylim=(.25, 2.25), yticks=[.5, 1, 1.5, 2])
ratio.xaxis.set_major_locator(FixedLocator([.2, .5, 1, 2, 5, 10]))
ratio.xaxis.set_major_formatter(FixedFormatter(['0.2', '0.5', '1', '2', '5', '10']))
ratio.xaxis.set_minor_formatter(NullFormatter())
ratio.grid(axis='y', color='#edf0f2', lw=.6)
fig.text(.035, .92, 'B', fontsize=17, fontweight='bold', color=navy)
fig.text(.14, .94, 'Xuan 2018 · population MSD', fontsize=12, color=navy)
fig.text(.14, .035, 'Design preview · symbols sample published curves, not individual cells.\n'
         'Same digitized data and fitted parameters as the current manuscript.',
         fontsize=8, color='#6d7c86')
for ext in ('png', 'pdf', 'svg'):
    fig.savefig(OUT / f'figure3B_proposal.{ext}', dpi=220, facecolor='white')
plt.close(fig)

with (OUT / 'ratio_curves.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['lag_h', 'digitized_msd_ratio', 'fitted_ou_msd_ratio'])
    writer.writerows(zip(lag, ratio_data, ratio_fit))
sources = [fit_path, FIGURE_PROVENANCE, Path(__file__),
           CAMPAIGN / 'scripts/figure_tracks.py',
           ROOT / 'wgd_langevin/persistence_only/campaign.py']
provenance = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'status': 'Standalone design proposal; manuscript unchanged',
              'selected_primary_fit': fit_path.name,
              'population_curve_crossovers_h': data_crossings,
              'fitted_ou_crossovers_h': fit_crossings,
              'ratio_interpolation': 'Linear in native lag/MSD, matching manuscript crossover',
              'data_symbols': 'Every twelfth native raster sample, plus endpoints; not independent cells',
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}}
(OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(json.dumps({'output': str(OUT / 'figure3B_proposal.png'),
                  'data_crossovers_h': data_crossings,
                  'fit_crossovers_h': fit_crossings,
                  'verification': 'Fitted predictions and both crossovers reproduce current manuscript'}, indent=2))
