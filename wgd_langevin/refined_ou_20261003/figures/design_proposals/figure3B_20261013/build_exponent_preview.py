"""Points for digitized MSD, model lines, and a matched-window slope comparison."""
from datetime import datetime, timezone
from pathlib import Path
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

ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
CAMPAIGN = ROOT / 'wgd_langevin/refined_ou_20261003'
source_provenance = CAMPAIGN / 'figures/main_03_migration_OU_provenance.json'
original = json.loads(source_provenance.read_text())
fit_path = CAMPAIGN / 'results' / original['selected_primary_fit']
assert hashlib.sha256(fit_path.read_bytes()).hexdigest() == original['source_sha256'][str(fit_path.relative_to(ROOT))]
fit = json.loads(fit_path.read_text())
sys.path.insert(0, str(CAMPAIGN / 'scripts'))
from figure_tracks import msd_expected

colors = ['#168b87', '#d87529']
labels = ['non-PGCC (WGD−)', 'PGCC (WGD+)']
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
log_edges = np.linspace(np.log(lower), np.log(upper), 7)
halfwidth = (log_edges[1] - log_edges[0]) / 2
log_centers = (log_edges[:-1] + log_edges[1:]) / 2
lag = np.geomspace(lower, upper, 600)

def log_slope(log_t, log_y):
    x = log_t - np.mean(log_t)
    return float(np.dot(x, log_y - np.mean(log_y)) / np.dot(x, x))

def model_window_slope(center, group):
    x = center + np.linspace(-halfwidth, halfwidth, 151)
    return log_slope(x, np.log(msd_expected(np.exp(x), *pars[group])))

# A power law must produce a constant recovered exponent with this transform.
check_x = np.linspace(np.log(lower), np.log(upper), 151)
for exponent in (.75, 1, 1.5, 2):
    assert abs(log_slope(check_x, np.log(7.0) + exponent * check_x) - exponent) < 1e-12
for group in (0, 1):
    rows = [r for r in fit['rows'] if r['study'] == 'Xuan2018'
            and r['metric'] == 'msd_curve' and int(r['group']) == group]
    assert np.allclose(msd_expected(np.array([r['lag_h'] for r in rows]), *pars[group]),
                       [r['predicted'] for r in rows], rtol=1e-10)

slope_rows = []
for group, (t, y) in enumerate(curves):
    x = np.log(t)
    for index, center in enumerate(log_centers):
        mask = (x >= log_edges[index]) & (x < log_edges[index+1])
        if index == len(log_centers)-1:
            mask = (x >= log_edges[index]) & (x <= log_edges[index+1])
        assert mask.sum() >= 20
        slope_rows.append({'group': group, 'window': index+1,
                           'lag_low_h': float(np.exp(log_edges[index])),
                           'lag_center_h': float(np.exp(center)),
                           'lag_high_h': float(np.exp(log_edges[index+1])),
                           'raster_samples': int(mask.sum()),
                           'observed_window_exponent': log_slope(x[mask], np.log(y[mask])),
                           'model_window_exponent': model_window_slope(center, group)})
model_centers = np.linspace(log_centers[0], log_centers[-1], 300)
model_slopes = np.array([[model_window_slope(c, g) for c in model_centers]
                         for g in (0, 1)])
assert np.all((model_slopes >= 1) & (model_slopes <= 2))
assert np.all(np.diff(model_slopes, axis=1) < 0)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig = plt.figure(figsize=(7.5, 6.2))
gs = fig.add_gridspec(2, 1, height_ratios=[2.25, 1.35], hspace=.15,
                      left=.14, right=.79, bottom=.16, top=.90)
ax = fig.add_subplot(gs[0])
slope_ax = fig.add_subplot(gs[1], sharex=ax)
display_counts = []
for group, (color, (t, y)) in enumerate(zip(colors, curves)):
    # These are display samples of a raster line, not recovered experimental lags.
    # Their count is a layout choice and is deliberately sparse.
    indices = np.unique(np.rint(np.linspace(0, len(t)-1, 12)).astype(int))
    display_counts.append(len(indices))
    ax.plot(t[indices], y[indices], linestyle='none', marker='o', markersize=2.65,
            markerfacecolor=color, markeredgewidth=0, zorder=3)
    ax.plot(lag, msd_expected(lag, *pars[group]), color=color, lw=1.65, zorder=2)
    ax.text(1.025, y[-1], labels[group], transform=ax.get_yaxis_transform(),
            color=color, va='center', fontsize=9.4, fontweight='bold')
    group_slopes = [r for r in slope_rows if r['group'] == group]
    slope_ax.plot(np.exp(model_centers), model_slopes[group], color=color, lw=1.65)
    slope_ax.plot([r['lag_center_h'] for r in group_slopes],
                  [r['observed_window_exponent'] for r in group_slopes],
                  linestyle='none', marker='o', color=color, markersize=4.3,
                  markeredgecolor='white', markeredgewidth=.45, zorder=3)
ax.set(xscale='log', yscale='log', ylabel='MSD (µm²)', ylim=(5, 11000))
ax.set_xlim(lower * .94, upper * 1.08)
ax.tick_params(axis='x', which='both', labelbottom=False)
ax.legend(handles=[Line2D([], [], color='#536470', marker='o', linestyle='none',
                          markeredgewidth=0, markersize=3.2, label='Display samples of digitized curve'),
                   Line2D([], [], color='#536470', lw=1.65, label='OU model')],
          loc='upper left', frameon=False, fontsize=8.5)
ax.grid(axis='y', which='major', color='#e9edf0', lw=.7)

slope_ax.set(ylabel='Local exponent, alpha', xlabel='Lag time (h; logarithmic scale)',
             ylim=(.72, 2.10), yticks=[1, 1.5, 2])
slope_ax.xaxis.set_major_locator(FixedLocator([.2, .5, 1, 2, 5, 10]))
slope_ax.xaxis.set_major_formatter(FixedFormatter(['0.2', '0.5', '1', '2', '5', '10']))
slope_ax.xaxis.set_minor_formatter(NullFormatter())
slope_ax.grid(axis='y', color='#e9edf0', lw=.7)
for level, text in ((2, 'Straight-line\nmotion'), (1, 'Diffusive\nscaling')):
    slope_ax.text(1.025, level, text, transform=slope_ax.get_yaxis_transform(),
                  va='center', fontsize=8.3, color='#6d7c86')
fig.text(.035, .92, 'B', fontsize=17, fontweight='bold', color='#19394d')
fig.text(.14, .94, 'Xuan 2018 · population MSD', fontsize=12, color='#19394d')
fig.text(.14, .025, 'Dots sample published lines; original experimental lag locations are unknown.\n'
         'Lower dots: slopes we estimated from digitized curves, not author-reported exponents.\n'
         'A constant exponent would be horizontal; observed/model slope windows have equal width.',
         fontsize=7.7, color='#6d7c86')
for ext in ('png', 'pdf', 'svg'):
    fig.savefig(OUT / f'figure3B_exponent_proposal.{ext}', dpi=220, facecolor='white')
plt.close(fig)

with (OUT / 'exponent_window_estimates.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(slope_rows[0]))
    writer.writeheader()
    writer.writerows(slope_rows)
with (OUT / 'model_window_exponents.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['lag_center_h', 'non_pgcc_exponent', 'pgcc_exponent'])
    writer.writerows(zip(np.exp(model_centers), *model_slopes))
sources = [fit_path, source_provenance, Path(__file__),
           CAMPAIGN / 'scripts/figure_tracks.py', ROOT / 'wgd_langevin/persistence_only/campaign.py']
provenance = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'status': 'Separate exponent-focused design preview; manuscript and ratio preview unchanged',
              'selected_primary_fit': fit_path.name,
              'upper_display_counts': display_counts,
              'upper_display_rule': 'Twelve display samples per raster curve, including endpoints; count is a layout choice, not experimental lag count',
              'experimental_lag_grid': 'Unknown. Display locations are not recovered author-calculated lag locations.',
              'previous_preview_correction': 'Earlier 64/66-dot display was every fourth raster column, not 64/66 recovered experimental lags. Corrected sparse display, legend, and footnote explicitly identify curve sampling.',
              'slope_estimator': 'OLS of log(MSD) versus log(lag); descriptive curve shape only',
              'observed_slope_windows': 6,
              'window_log10_width': float(2 * halfwidth / np.log(10)),
              'model_slope_estimator': 'Moving log-linear fit using 151 equally spaced log lags over identical window width',
              'slope_display_extent': 'Only centers for which the entire model fitting window is in the measured lag interval',
              'uncertainty': 'Raster samples are correlated. Slopes depend on window choice; no biological confidence intervals inferred.',
              'caption': 'Upper: twelve display samples of each digitized published population MSD line, with current fitted OU expectations as solid lines. The display points are a visualization choice; neither their number nor their lag locations reconstruct the original experimental MSD lag grid, which is unknown. Lower: slopes we estimated in six equal log-lag windows of the digitized curves, with matching-window moving OU exponents as solid lines. These are not author-reported local exponents. Constant exponents would appear horizontal. All native raster samples within each slope window enter the estimate; raster samples are correlated and are not independent cell observations.',
              'validation': 'Current fit SHA-256 and all native MSD predictions match manuscript; known power-law exponents recovered to 1e-12; OU window exponents remain between 1 and 2 and decrease with lag.',
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}}
(OUT / 'exponent_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
print(json.dumps({'output': str(OUT / 'figure3B_exponent_proposal.png'),
                  'display_points_per_group': display_counts,
                  'window_width_decades': provenance['window_log10_width'],
                  'slope_rows': slope_rows,
                  'validation': provenance['validation']}, indent=2))
