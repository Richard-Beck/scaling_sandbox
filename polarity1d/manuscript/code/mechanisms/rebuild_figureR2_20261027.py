"""Six-case R2 display rebuild; no model tuning or event-log regeneration.

extract: fixed-prefix sample profiles and analytic nulls. render: portable CSVs
only. Missing Spring profiles and model statistics have bounded SLURM helpers.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('MPLCONFIGDIR', f'/tmp/r2-display-mpl-{os.getuid()}')
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from scipy.special import hyp2f1

ROOT = Path(__file__).resolve().parents[4]
MS = ROOT / 'polarity1d/manuscript'
OUT = MS / 'figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/figureR2'
DATA = OUT / 'supporting/data'
META = OUT / 'supporting/metadata'
SOURCE = MS / 'data/runs/poisson_ring_modes_20261027'
REPLAY = MS / 'data/runs/ring_r2_display_profiles_20261029'
R7 = MS / 'figures/length_mechanisms_20261008/dev/revision_round7_20261026/figure4/supporting/data'
CASES = ('relocation_constant', 'relocation_inverse_square', 'patch_field',
         'goryachev_perimeter', 'spring_perimeter', 'spring_cell')
LABELS = ('Per-perimeter relocation · constant success', 'Per-perimeter relocation · success ∝ 1/C²',
          'Patch field · uniform births', 'Goryachev · per-perimeter events',
          'Spring · per-perimeter events', 'Spring · per-cell events')
COLORS = {20: '#E69F00', 60: '#0072B2', 1: '#E69F00', 3: '#0072B2'}
PANELS = tuple('ABCDEFGHIJKL')

def now():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text())

def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def relative(path):
    return str(Path(path).relative_to(ROOT))

def profile_path(case, size):
    return DATA / f'{case}_C{size:g}_profile.csv.gz'

def export_profile(case, size, time, circumference, field):
    """Wide CSV: physical x_j(t)=(j+1/2) C(t)/n; signed activity retained."""
    frame = pd.DataFrame(field, columns=[f'activity_{i:03d}' for i in range(field.shape[1])])
    frame.insert(0, 'circumference', np.broadcast_to(circumference, len(time)))
    frame.insert(0, 'time', time)
    frame.to_csv(profile_path(case, size), index=False, float_format='%.8g')

def load_profile(case, size):
    frame = pd.read_csv(profile_path(case, size))
    return frame.time.to_numpy(), frame.circumference.to_numpy(), frame.filter(like='activity_').to_numpy()

def import_long(path, case, size):
    frame = pd.read_csv(path)
    matrix = frame.pivot(index='time', columns='x' if 'x' in frame else 'relative_arc_position', values='activity')
    export_profile(case, size, matrix.index.to_numpy(), size, matrix.to_numpy())

def nulls():
    """Reuse relocation records; draw a stationary immigration/death field."""
    sources = []
    for size in (20, 60):
        src = R7 / f'A_perimeter_C{size}.csv.gz'
        import_long(src, CASES[0], size)
        sources.append(dict(case_id=CASES[0], size=size, path=relative(src), sha256=sha(src)))
        # This immutable compact sample is reused, never redrawn or searched.
        if not profile_path(CASES[1], size).exists():
            raise FileNotFoundError('Retained inverse-success sample path is required')
        rng = np.random.default_rng(np.random.SeedSequence([20261029, 2, size]))
        r = 1 / 7000
        b = r / 6
        t = np.arange(0., 24000. + 20., 20.)
        x = (np.arange(240) + .5) * size / 240
        n = rng.poisson(size / 6)
        centres = list(rng.uniform(0, size, n))
        starts = [0.] * n
        stops = list(rng.exponential(1 / r, n))
        birth = 0.
        while True:
            birth += rng.exponential(1 / (b * size))
            if birth > t[-1]:
                break
            starts.append(birth)
            centres.append(rng.uniform(0, size))
            stops.append(birth + rng.exponential(1 / r))
        field = np.zeros((len(t), len(x)))
        for start, stop, centre in zip(starts, stops, centres):
            distance = (x - centre + size / 2) % size - size / 2
            field[(t >= start) & (t < stop)] += np.exp(-.5 * (distance / 1.5) ** 2)
        export_profile(CASES[2], size, t, size, field)
        sources.append(dict(case_id=CASES[2], size=size, seed=[20261029, 2, size],
                            initialization='Stationary Poisson(C/6) live patches', birth_density=b,
                            death_rate=r, gaussian_sd=1.5, initial_patches=int(n), all_patches=len(starts)))
    rows = []
    for case in CASES[:3]:
        for size in (20, 60):
            for lag in np.linspace(0, 24000, 481):
                if case == CASES[2]:
                    rho = np.exp(-lag / 7000)
                    value = np.pi / 4 * rho * hyp2f1(.5, .5, 2, rho * rho)
                else:
                    rate = (1 / 7000) * (size / 40 if case == CASES[0] else 40 / size)
                    value = np.exp(-rate * lag)
                rows.append(dict(case_id=case, size=size, lag=lag, correlation=value))
    pd.DataFrame(rows).to_csv(DATA / 'null_heading.csv', index=False)
    if not (META / 'null_sources.json').exists():
        write(META / 'null_sources.json', dict(generated_utc=now(), sources=sources))

def model_profiles():
    manifest = read(SOURCE / 'manifest.json')
    records = []
    for case, model, sizes, window in ((CASES[3], 'goryachev', (20, 60), 2000.),
                                       (CASES[4], 'spring_ring', (1, 3), 200.)):
        for size in sizes:
            job = next(j for j in manifest['jobs'] if j['model'] == model and j['size'] == size
                       and j['rate_level'] == 'low' and j['amplitude_level'] == 'high' and j['replicate'] == 0)
            folder = SOURCE / f'workers/{job["slot"]:03d}/segments/00'
            time = np.load(folder / 'time.npy', mmap_mode='r')
            keep = np.flatnonzero(time <= time[0] + window + 1e-8)
            keep = keep[np.unique(np.linspace(0, len(keep) - 1, 801).round().astype(int))]
            coefficients = np.load(folder / 'u_modes.npy', mmap_mode='r')[keep]
            q = (np.arange(240) + .5) / 240
            phase = np.exp(2j * np.pi * np.arange(coefficients.shape[1])[:, None] * q)
            field = coefficients[:, 0].real[:, None] + 2 * np.real(coefficients[:, 1:] @ phase[1:])
            circumference = np.load(folder / 'circumference.npy', mmap_mode='r')[keep] if model == 'spring_ring' else size
            export_profile(case, size, time[keep] - time[0], circumference, field)
            records.append(dict(case_id=case, size=size, seed=job['seed'], slot=job['slot'],
                                source=relative(folder), source_start=float(time[0]), duration=window,
                                mode_count=int(coefficients.shape[1]), reconstruction_points=240,
                                u_modes_sha256=sha(folder / 'u_modes.npy')))
    replay_check = read(REPLAY / 'validation.json')
    if not replay_check['passed']:
        raise ValueError('Missing-profile replay has not passed validation')
    for slot, size in ((0, 1), (1, 3)):
        folder = REPLAY / f'workers/{slot:03d}/segments/00'
        time = np.load(folder / 'time.npy')
        keep = np.unique(np.linspace(0, len(time) - 1, 801).round().astype(int))
        export_profile(CASES[5], size, time[keep] - time[0], np.load(folder / 'circumference.npy')[keep],
                       np.load(folder / 'u_profiles.npy', mmap_mode='r')[keep])
        records.append(dict(case_id=CASES[5], size=size, source=relative(folder),
                            duration=float(time[-1] - time[0]), native_full_profiles=True,
                            profile_sha256=sha(folder / 'u_profiles.npy')))
    write(META / 'profile_sources.json', dict(generated_utc=now(), records=records,
          selection='Fixed replicate 0; first 40 event decay times after burn, 801 uniformly spaced frames; no heading-based selection.',
          replay_validation=replay_check))

def extract():
    DATA.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)
    nulls()
    model_profiles()

def render(round_number):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    from PIL import Image
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.linewidth': .6,
                         'lines.linewidth': 1.25, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig = plt.figure(figsize=(180 / 25.4, 260 / 25.4))
    bottoms = [.842, .685, .528, .345, .188, .031]
    height = .085
    null = pd.read_csv(DATA / 'null_heading.csv')
    model = pd.read_csv(DATA / 'model_heading.csv')
    profiles = {(case, size): load_profile(case, size) for i, case in enumerate(CASES)
                for size in ((1, 3) if i >= 4 else (20, 60))}
    maxima = {'relocation': 1., 'patch_field': max(profiles[(CASES[2], s)][2].max() for s in (20, 60)),
              'goryachev': max(profiles[(CASES[3], s)][2].max() for s in (20, 60)),
              'spring': max(profiles[(case, s)][2].max() for case in CASES[4:] for s in (1, 3))}
    boxes = {}
    heading_records = []
    crossings = []
    for i, (case, label, bottom) in enumerate(zip(CASES, LABELS, bottoms)):
        sizes = (1, 3) if i >= 4 else (20, 60)
        pair = [profiles[(case, s)] for s in sizes]
        extents = [float(c.max()) for t, c, v in pair]
        left, right, gap = .10, .535, .020
        widths = (right - left - gap) * np.array(extents) / sum(extents)
        xpos = [left, left + widths[0] + gap]
        group = 'relocation' if i < 2 else 'patch_field' if i == 2 else 'goryachev' if i == 3 else 'spring'
        vmax = maxima[group]
        for j, (size, (time, circumference, values)) in enumerate(zip(sizes, pair)):
            ax = fig.add_axes([xpos[j], bottom, widths[j], height])
            n = values.shape[1]
            x = circumference[:, None] * (np.arange(n)[None, :] + .5) / n
            mesh = ax.pcolormesh(x, np.broadcast_to(time[:, None], x.shape), np.maximum(values, 0),
                                 shading='nearest', cmap='viridis', vmin=0, vmax=vmax, rasterized=True, linewidth=0)
            if i >= 3:
                phase = np.exp(2j*np.pi*(np.arange(n)+.5)/n)
                P = (values @ phase)/values.sum(axis=1)
                position = np.mod(np.angle(P), 2*np.pi)/(2*np.pi)*circumference
                good = np.isfinite(P) & (abs(P) >= .1)
                line = np.where(good, position, np.nan)
                # Do not connect across the periodic seam or excluded frames.
                seam = np.r_[False, abs(np.diff(position/circumference)) > .5]
                line[seam] = np.nan
                ax.plot(line, time, color='white', lw=.6)
                heading_records.extend(dict(case_id=case, size=size, time=t, arclength=p,
                                             coherence=c, plotted=bool(ok and not cut))
                                       for t, p, c, ok, cut in zip(time, position, abs(P), good, seam))
            ax.set(xlim=(0, extents[j]), ylim=(0, time[-1]), xticks=[0, sizes[j]],
                   ylabel='Time' if j == 0 else '')
            ax.set_yticks([0, time[-1]/2, time[-1]])
            if j:
                ax.tick_params(labelleft=False, left=False)
            ax.set_title(f'C{size}', fontsize=7, color='black', pad=3, fontweight='normal')
            ax.tick_params(labelsize=6.5, pad=2)
        cb = fig.colorbar(mesh, cax=fig.add_axes([.555, bottom, .010, height]), ticks=[0, vmax])
        cb.ax.set_yticklabels(['0', f'{vmax:.2g}'])
        cb.ax.tick_params(labelsize=6, length=2, pad=2)
        cb.set_label('Activity', fontsize=6.5, labelpad=1)
        cb.outline.set_linewidth(.4)
        if i >= 3:
            import matplotlib.patheffects as pe
            handle = Line2D([], [], color='white', lw=.7, label='Heading',
                            path_effects=[pe.Stroke(linewidth=2, foreground='.3'), pe.Normal()])
            fig.legend(handles=[handle], loc='lower right', bbox_to_anchor=(.56, bottom+height+.004),
                       fontsize=5.8, frameon=False, handlelength=1.4, borderaxespad=0)
        fig.text((left+right)/2, bottom-.025, 'Arclength', fontsize=7, ha='center')
        fig.text(.10, bottom + height + .027, label, fontsize=7.5, color='black', va='bottom')
        fig.text(.038, bottom + height + .013, PANELS[2*i], fontsize=11, color='black', fontweight='bold', va='bottom')
        ax = fig.add_axes([.715, bottom, .258, height])
        frame = null if i < 3 else model
        frame = frame[frame.case_id == case]
        for size in sizes:
            q = frame[frame['size'] == size].sort_values('lag')
            if i >= 3:
                q = q[q.lag > 0]
            style = '--' if size == sizes[-1] else '-'
            ax.plot(q.lag, q.correlation, color=COLORS[size], ls=style, zorder=3 if size == sizes[0] else 4)
            if i >= 3:
                ax.fill_between(q.lag, q.ci_low, q.ci_high, color=COLORS[size], alpha=.23, lw=0)
                below = np.flatnonzero(q.correlation.to_numpy() <= .5)
                if len(below) and below[0] > 0:
                    k = below[0]
                    a, b = q.iloc[k-1], q.iloc[k]
                    weight = (.5-a.correlation)/(b.correlation-a.correlation)
                    lag = np.exp(np.log(a.lag)+weight*np.log(b.lag/a.lag))
                    ax.plot(lag, .5, 'o', ms=3, color=COLORS[size], zorder=6)
                    crossings.append(dict(case_id=case, size=size, lag=lag, correlation=.5,
                                          interpolation='Intersection with displayed log-lag line segment'))
        ax.set(xlim=(0, frame.lag.max()), ylim=(-.3 if i >= 3 else 0, 1.04),
               xlabel='Lag (log scale)' if i >= 3 else 'Lag', ylabel='Heading correlation', yticks=[0, .5, 1])
        if i >= 3:
            from matplotlib.ticker import NullLocator
            ax.axhline(.5, color='.7', lw=.6, ls=':', zorder=0)
            ax.set_xscale('log')
            ax.set_xlim(frame.loc[frame.lag > 0, 'lag'].min(), frame.lag.max())
            ticks = [1, 100, 10000] if i == 3 else [.1, 10, 1000]
            ax.set_xticks(ticks, labels=[f'{x:g}' for x in ticks])
            ax.xaxis.set_minor_locator(NullLocator())
        else:
            ax.xaxis.set_major_locator(MaxNLocator(3, steps=[1, 2, 4, 5, 10]))
            ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.tick_params(labelsize=6.5, pad=2)
        handles = [Line2D([], [], color=COLORS[s], ls='--' if s == sizes[-1] else '-', label=f'C{s}') for s in sizes]
        if i >= 3:
            from matplotlib.patches import Patch
            from matplotlib.legend_handler import HandlerTuple
            handles.append(tuple(Patch(facecolor=COLORS[s], alpha=.23, edgecolor='none') for s in sizes))
        main_key = ax.legend(handles=handles, frameon=False, ncol=len(handles), loc='lower left', bbox_to_anchor=(-.03, 1.02),
                  labels=[f'C{s}' for s in sizes]+(['95% CI'] if i >= 3 else []),
                  handler_map={tuple: HandlerTuple(ndivide=None)} if i >= 3 else {},
                  fontsize=6.1, borderaxespad=0, handlelength=1.4, columnspacing=.7, handletextpad=.4)
        if i >= 3:
            ax.add_artist(main_key)
            ax.legend(handles=[Line2D([], [], color='.3', marker='o', ls='none', ms=3, label='Half-height')],
                      loc='lower left', fontsize=5.7, frameon=False, borderaxespad=.25,
                      handlelength=1, handletextpad=.3)
        fig.text(.715, bottom+height+.027,
                 'Many-patch approximation' if i == 2 else 'Closed form' if i < 3 else 'Simulation ensemble',
                 fontsize=7, color='black', va='bottom')
        fig.text(.663, bottom + height + .013, PANELS[2*i+1], fontsize=11, color='black', fontweight='bold', va='bottom')
        boxes[PANELS[2*i]] = [.025, bottom-.033, .625, bottom+height+.045]
        boxes[PANELS[2*i+1]] = [.64, bottom-.033, .995, bottom+height+.045]
    fig.add_artist(Line2D([.038, .975], [.477, .477], transform=fig.transFigure, color='.45', lw=.7))
    fig.text(.10, .483, 'Model events · low rate, high amplitude', fontsize=7, color='black')
    pd.DataFrame(heading_records).to_csv(DATA/'profile_heading.csv', index=False, float_format='%.8g')
    pd.DataFrame(crossings).to_csv(DATA/'half_height_markers.csv', index=False, float_format='%.8g')
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(OUT / f'figureR2_main.{ext}', dpi=300)
    svg = OUT/'figureR2_main.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)
    review = META / f'review_round_{round_number}'
    review.mkdir(exist_ok=True)
    (review/'legend.md').write_text((META/'main_legend.md').read_text().replace('(../../README.md)', '(../../../README.md)'))
    with Image.open(OUT / 'figureR2_main.png') as picture:
        picture.save(review / 'print_180mm.png', dpi=(300, 300))
        width, pixels_high = picture.size
        for letter, (x0, y0, x1, y1) in boxes.items():
            picture.crop((round(x0*width), max(0, round((1-y1)*pixels_high)),
                          round(x1*width), min(pixels_high, round((1-y0)*pixels_high)))).save(review / f'panel_{letter}.png', dpi=(300, 300))
    write(META / 'render.json', dict(rendered_utc=now(), print_width_mm=180, print_height_mm=260,
          round=round_number, renderer_sha256=sha(__file__), color_maxima=maxima,
          output_sha256={p.name: sha(p) for p in OUT.glob('figureR2_main.*')},
          data_sha256={p.name: sha(p) for p in sorted(DATA.glob('*.csv*'))}, panel_boxes=boxes))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['extract', 'render', 'all'])
    parser.add_argument('--round', type=int, default=1)
    args = parser.parse_args()
    if args.command in ('extract', 'all'):
        extract()
    if args.command in ('render', 'all'):
        render(args.round)
