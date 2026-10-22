"""Standalone draft MCA profiles and uniform-replacement means; no integration.

Reuse corrected independent-endpoint relative-20% trajectories and audited
native-step budgets. The canonical report and all original data stay intact.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import length_dev_figure_x_overview as overview
from length_dev_output_layout import dev_output_path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'dev/mca_length_drafts_20261014'
BUDGET = ROOT / 'data/runs/debelly_mca_speed_decomposition_20261010'
LENGTHS = [23, 31, 39]
TIMES = [0., 30., 180.]
COLORS = {23: '#332288', 31: '#882255', 39: '#999933'}
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def read(path):
    return json.loads(source(path).read_text())


def output(name):
    return dev_output_path(OUT, name)


def write_csv(name, rows):
    with output(name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, stem):
    for ext in ['.png', '.pdf', '.svg']:
        fig.savefig(output(stem + ext), dpi=180, facecolor='white')
    plt.close(fig)


def clean_axes(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(direction='out')
    ax.grid(axis='y', color='#e6e6e6', lw=.7)
    ax.set_axisbelow(True)


def profile_figure(profiles, normalized):
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.7), sharey=True)
    fig.subplots_adjust(left=.073, right=.99, bottom=.21, top=.72, wspace=.13)
    fig.text(.073, .955, 'DRAFT · Bound MCA at the same elapsed times', fontsize=17, weight='bold')
    fig.text(.073, .892, 'Right cue: peak 3.6, support 20% of each cell length; cue remains on at 30 and 180 s.', fontsize=11)
    all_values = np.concatenate([p['bound'] for p in profiles.values()])
    lo, hi = float(all_values.min()), float(all_values.max())
    pad = .075 * (hi - lo)
    for col, (ax, time) in enumerate(zip(axes, TIMES)):
        for length in LENGTHS:
            p = profiles[length, time]
            ax.plot(p['x'] / length if normalized else p['x'], p['bound'],
                    color=COLORS[length], lw=2.2, label=f'{length} µm')
        ax.set_title(['Basal state', '30 s after cue onset', '180 s after cue onset'][col], fontsize=13, pad=43)
        ax.legend(title='Cell length', ncol=3, loc='lower left', bbox_to_anchor=(0, 1.015),
                  frameon=False, fontsize=10, title_fontsize=10, columnspacing=1.0,
                  handlelength=1.6, borderaxespad=0)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(0, 1 if normalized else 39)
        ax.set_xticks([0, .25, .5, .75, 1] if normalized else [0, 10, 20, 30, 39])
        ax.set_xlabel('Reference position, x / L' if normalized else 'Reference position, x (µm)')
        clean_axes(ax)
    axes[0].set_ylabel('Bound MCA (native concentration)')
    note = ('Position is normalized; concentration is unchanged.' if normalized else
            'Curves end at their own nominal cell lengths; x is the undeformed reference coordinate.')
    fig.text(.073, .065, note + ' All panels use the same concentration limits.', fontsize=10)
    stem = 'bound_MCA_profiles_' + ('normalized_position' if normalized else 'physical_position')
    save_figure(fig, stem)
    return stem


def replacement_data():
    manifest = read(BUDGET / 'manifest.json')
    status = read(BUDGET / 'status.json')
    assert status['complete'] and status['audit_passed']
    assert manifest['geometry'] == 'fixed_relative'
    rows, checks = [], []
    for slot, length in enumerate(manifest['lengths_um']):
        base = BUDGET / 'workers' / f'{slot:03d}'
        meta = read(base.with_suffix('.json'))
        path = source(base.with_suffix('.csv'))
        assert meta['passed'] and sha(path) == meta['csv_sha256']
        p = meta['parameters']
        assert p['center_boundary_motion'] is False
        fraction_bound = p['mca_on'] / (p['mca_on'] + p['mca_off'])
        with path.open() as stream:
            native = list(csv.DictReader(stream))
        start = np.array([float(r['time_start_s']) for r in native])
        end = np.array([float(r['time_end_s']) for r in native])
        dt = np.array([float(r['step_s']) for r in native])
        rate = np.array([float(r['bound_uniform_source_conc_s']) for r in native])
        duration = float(native[0]['cue_duration_s'])
        assert duration > 180 and start[0] == 0 and end[-1] == duration
        assert all(int(r['length_um']) == length for r in native)
        np.testing.assert_allclose(start[1:], end[:-1], rtol=0, atol=1e-12)
        np.testing.assert_allclose(end - start, dt, rtol=0, atol=1e-12)
        np.testing.assert_allclose(rate, -np.array([float(r['front_depletion_uniform_replenishment_conc_s'])
                                                   for r in native]), rtol=1e-12, atol=1e-16)
        for window, stop in [('whole_cue', duration), ('first_180_s', 180.)]:
            weight = np.maximum(0., np.minimum(end, stop) - start)
            assert abs(weight.sum() - stop) < 1e-9
            mean_bound = float(np.sum(rate * weight) / stop)
            mean_total = mean_bound / fraction_bound
            rows.append(dict(length_um=length, window=window, window_start_s=0., window_end_s=stop,
                             cue_duration_s=duration, mean_bound_replacement_conc_s=mean_bound,
                             mean_unbound_replacement_conc_s=mean_total - mean_bound,
                             mean_total_replacement_conc_s=mean_total,
                             mean_total_replacement_amount_s=length * mean_total,
                             bound_fraction_of_uniform_correction=fraction_bound))
        checks.append(dict(length_um=length, intervals=len(native), duration_s=duration,
                           native_csv_hash_verified=True, continuous_time_coverage_verified=True,
                           uniform_rate_matches_local_budget=True,
                           original_numerical_errors=meta['checks'][0]['errors']))
    return rows, checks


def replacement_figure(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)
    fig.subplots_adjust(left=.088, right=.99, bottom=.22, top=.72, wspace=.16)
    fig.text(.088, .954, 'DRAFT · Mean uniform MCA replacement during the reversing cue', fontsize=16, weight='bold')
    fig.text(.088, .89, 'Relative-20% right cue · time-weighted native correction rates · 12 eligible lengths', fontsize=11)
    for ax, window, title in zip(axes, ['whole_cue', 'first_180_s'],
                                 ['Whole selected cue', 'Common first 180 s']):
        selected = [r for r in rows if r['window'] == window]
        lengths = [r['length_um'] for r in selected]
        for column, color, marker, label in [
                ('mean_total_replacement_conc_s', '#666666', 'o', 'Total MCA'),
                ('mean_bound_replacement_conc_s', '#CC79A7', 's', 'Bound MCA')]:
            ax.plot(lengths, [r[column] for r in selected], color=color, marker=marker,
                    markersize=5, lw=1.8, label=label)
        ax.set_title(title, fontsize=13, pad=43)
        ax.legend(ncol=2, loc='lower left', bbox_to_anchor=(0, 1.015), frameon=False,
                  fontsize=11, borderaxespad=0)
        ax.set_xlabel('Nominal cell length (µm)')
        ax.set_xticks([18, 23, 27, 31, 35, 39])
        ax.set_xlim(17, 40)
        ax.set_ylim(bottom=0)
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_powerlimits((-4, -4))
        ax.yaxis.set_major_formatter(formatter)
        clean_axes(ax)
    axes[0].set_ylabel('Mean replacement rate (native concentration / s)')
    fig.text(.088, .087, 'Whole-cue durations vary with length (226–809 s); every cue exceeds 180 s.', fontsize=10)
    fig.text(.088, .043, 'The uniform correction assigns 78.125% to bound MCA and 21.875% to unbound MCA.', fontsize=10)
    stem = 'mean_uniform_replacement_vs_length'
    save_figure(fig, stem)
    return stem


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    protected = {name: sha(REPORT / name) for name in ['index.html', 'length_mechanisms.pdf']}
    traces, _, profile_audits, _ = overview.load()
    SOURCES.update(overview.SOURCES)
    source(Path(__file__))
    source(Path(overview.__file__))
    source(Path(__file__).with_name('length_dev_output_layout.py'))
    profiles, records = {}, []
    for length in LENGTHS:
        trace = traces[length]
        for time in TIMES:
            index = np.flatnonzero(np.isclose(trace['time'], time, rtol=0, atol=1e-10))
            assert len(index) == 1 and time < trace['duration']
            bound = trace['fields']['mca_bound'][index[0]].copy()
            profiles[length, time] = dict(x=trace['x'].copy(), bound=bound)
            for x, b in zip(trace['x'], bound):
                records.append(dict(length_um=length, elapsed_time_s=time, x_um=float(x),
                                    x_over_L=float(x / length), bound_MCA=float(b),
                                    selected_cue_duration_s=trace['duration']))
    write_csv('bound_MCA_profile_samples.csv', records)
    stems = [profile_figure(profiles, True), profile_figure(profiles, False)]
    rates, rate_audits = replacement_data()
    write_csv('uniform_replacement_means.csv', rates)
    stems.append(replacement_figure(rates))
    generated = datetime.now(timezone.utc).isoformat()
    captions = [
        'Bound MCA at the exact saved basal, 30-s and 180-s samples for 23, 31 and 39 µm. '
        'The peak-3.6 raised-cosine reversing cue occupies the right 20% of nominal length. '
        'All curves share concentration limits; x/L changes position only. The basal sample '
        'is the exact qualifying checkpoint, with no reflection or additional settling.',
        'The same profiles in undeformed reference coordinates measured in µm. '
        'The common x range exposes the different cell lengths; each curve ends at its own native cell centres.',
        'Time-weighted mean of the native spatially uniform MCA concentration correction, '
        'not net binding or endpoint flux alone. Left: average over each length’s own shortest '
        'tested scalar-strength-passing cue. Right: average over the common 0–180-s cue window. '
        'Total MCA is the sum of bound and unbound correction rates; the fixed native resting '
        'partition assigns 78.125% to bound MCA. The correction closes the pre-correction '
        'total-pool deficit, including endpoint flux and the native sequential-exchange residual. '
        'Values are concentration per second; spatially integrated total amount per second '
        'is exported separately in the CSV. Connecting lines guide the eye; no fit is used.'
    ]
    for stem, caption in zip(stems, captions):
        output(stem + '_caption.txt').write_text(caption + '\n')
    sections = []
    for stem, caption in zip(stems, captions):
        sections.append(f'<section><h2>{stem.replace("_", " ")}</h2>'
                        f'<p><a href="supporting/vectors/{stem}.pdf">PDF</a> · '
                        f'<a href="supporting/vectors/{stem}.svg">SVG</a></p>'
                        f'<img src="{stem}.png" alt="{stem.replace("_", " ")}">'
                        f'<p>{caption}</p></section>')
    (OUT / 'index.html').write_text('<!doctype html><html><head><meta charset="utf-8">'
        '<title>Draft MCA comparisons</title><style>body{font:16px system-ui;max-width:1400px;'
        'margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.5}section{margin:40px 0}'
        '</style></head><body><h1>Draft MCA length comparisons</h1><p>Generated ' + generated +
        ' · producer finished · saved trajectories only.</p><p>'
        '<a href="supporting/data/bound_MCA_profile_samples.csv">Profile samples</a> · '
        '<a href="supporting/data/uniform_replacement_means.csv">Replacement means</a></p>' +
        ''.join(sections) + '</body></html>')
    assert protected == {name: sha(REPORT / name) for name in protected}
    audit = dict(generated_at=generated, producer='finished', passed=True,
                 profile_lengths_um=LENGTHS, exact_profile_times_s=TIMES,
                 profile_source_checks=profile_audits, replacement_checks=rate_audits,
                 protected_report_sha256=protected, source_sha256=SOURCES,
                 no_integration=True, no_profile_interpolation=True,
                 rate_average='sum(native per-step rate * overlap duration) / window duration')
    output('draft_mca_length_comparison_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    artifacts = [p for p in OUT.rglob('*') if p.is_file() and not p.name.endswith('.manifest.json')]
    output('draft_mca_length_comparison.manifest.json').write_text(json.dumps(dict(
        generated_at=generated, producer='finished', passed=True, artifact_root='../..',
        sha256={str(p.relative_to(OUT)): sha(p) for p in artifacts}, source_sha256=SOURCES,
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '
        '/home/4473331/.conda/envs/scaling_sandbox/bin/python -B '
        'polarity1d/manuscript/code/mechanisms/draft_mca_length_comparison_20261014.py'), indent=2) + '\n')
    print(json.dumps(dict(output=str(OUT), generated_at=generated, figures=stems,
                          profile_samples=len(records), rate_means=len(rates)), indent=2))


if __name__ == '__main__':
    main()
