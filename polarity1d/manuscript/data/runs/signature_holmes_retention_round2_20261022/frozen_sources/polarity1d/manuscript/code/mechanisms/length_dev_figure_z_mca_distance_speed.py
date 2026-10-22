"""Draft Figure Z from audited, strength-selected corrected native histories.

No integration or manuscript-report changes. Centroids use nominal coordinates.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'figures/length_mechanisms_20261008/dev/figure_z_mca_redistribution'
RUN = ROOT / 'data/runs/length_narrative_memory_20261009'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
GEOS = ['fixed_physical', 'fixed_relative', 'fixed_slope']
TITLES = ['Fixed patch: 3 µm', 'Relative patch: 20%', 'Constant slope: 0.24']
LENGTHS = [18, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39]
WINDOW = 60.0
STEM = 'Figure_Z_MCA_distance_speed'
SOURCES = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def register(p):
    SOURCES[str(p.relative_to(ROOT))] = sha(p)


def read(p):
    register(p)
    return json.loads(p.read_text())


def output(name):
    return layout.dev_output_path(OUT, name)


def write(name, value):
    output(name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def centroid(b, x):
    assert np.isfinite(b).all() and b.min() >= -1e-10
    assert np.all(b.sum(axis=-1) > 0)
    return b @ x / b.sum(axis=-1)


def measure():
    manifest = read(RUN / 'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(RUN / 'status.json')['audit_passed']
    assert read(BASE / 'results/independent_audit.json')['passed']
    assert read(ROOT / 'data/runs/length_mechanisms_20261008/narrative_20261009/memory/independent_audit.json')['passed']
    records_path = BASE / 'results/records.json.gz'
    register(records_path)
    with gzip.open(records_path, 'rt') as stream:
        records = {r['case_index']: r for r in json.load(stream)['cases']}
    rows, checks = [], []
    for job in manifest['jobs']:
        case = job['case']
        L, geo, D = int(case['length']), case['condition'], float(job['duration'])
        record = records[case['case_index']]
        assert D == record['result']['full_polarity']['endpoint']
        assert D == min(t['duration'] for t in record['tests'] if t['full_polarity_success'])
        assert D >= 2 * WINDOW
        path = ROOT / job['source'] if job['source'] else RUN / f"workers/{manifest['missing_jobs'].index(job):03d}.npz"
        meta = read(path.with_suffix('.json'))
        register(path)
        assert sha(path) == meta['archive_sha256']
        checkpoint_path = ROOT / case['checkpoint']
        register(checkpoint_path)
        assert sha(checkpoint_path) == case['checkpoint_sha256']
        with np.load(checkpoint_path) as checkpoint:
            initial_b = checkpoint['state_mca_bound'].copy()
            initial_rac = checkpoint['state_rac'].copy()
        with np.load(path) as z:
            x, t, b = z['x'], z['pulse_time'], z['pulse_field_mca_bound']
            np.testing.assert_allclose(x, (np.arange(len(x)) + .5) * L / len(x), atol=1e-12)
            np.testing.assert_allclose(b[0], initial_b, atol=1e-12)
            np.testing.assert_allclose(z['pulse_field_rac'][0], initial_rac, atol=1e-12)
            np.testing.assert_allclose(z['recovery_field_mca_bound'][0], b[-1], atol=1e-12)
            np.testing.assert_allclose(z['recovery_field_rac'][0], z['pulse_field_rac'][-1], atol=1e-12)
            assert abs(t[0]) < 1e-12 and abs(t[-1] - D) < 1e-9
            assert np.all(np.diff(t) > 0)
            c = centroid(b, x)
            tail = z['recovery_time'] >= 9500
            rac = z['recovery_field_rac'][tail]
            polarity = rac @ (2 * x / L - 1) / rac.sum(axis=-1)
            contrast = np.ptp(rac, axis=1)
            strength_error = float(np.max(abs(polarity / abs(case['initial_orientation']) - 1)))
            assert polarity.min() > .05 and contrast.min() >= .05 and np.ptp(polarity) < .005
            assert strength_error <= .01 and meta['metrics']['full_polarity_success']
            c0, cD = float(c[0]), float(c[-1])
            reflected = float(centroid(initial_b[::-1], x))
            np.testing.assert_allclose(reflected, L - c0, atol=1e-12)
            np.testing.assert_allclose(c, L / 2 * (1 + b @ (2 * x / L - 1) / b.sum(axis=-1)), atol=1e-12)
            assert c0 > L / 2
            delta = c0 - cD
            early = (c0 - float(np.interp(WINDOW, t, c))) / WINDOW
            late = (float(np.interp(D - WINDOW, t, c)) - cD) / WINDOW
            average = delta / D
            row = dict(condition=geo, length_um=L, cue_duration_s=D,
                       initial_centroid_um=c0, withdrawal_centroid_um=cD,
                       reflected_target_centroid_um=reflected,
                       initial_distance_to_midpoint_um=c0 - L / 2,
                       initial_distance_to_reflected_state_um=c0 - reflected,
                       centroid_shift_during_cue_um=delta,
                       centroid_path_length_during_cue_um=float(np.sum(abs(np.diff(c)))),
                       remaining_distance_to_midpoint_um=abs(cD - L / 2),
                       remaining_distance_to_reflected_state_um=abs(cD - reflected),
                       early_mean_speed_um_s=early, late_mean_speed_um_s=late,
                       cue_average_speed_um_s=average,
                       early_mean_speed_per_length_s=early / L,
                       late_mean_speed_per_length_s=late / L,
                       cue_average_speed_per_length_s=average / L,
                       cue_fraction_of_reflected_centroid_distance=delta / (c0 - reflected),
                       early_window_start_s=0, early_window_end_s=WINDOW,
                       late_window_start_s=D - WINDOW, late_window_end_s=D,
                       minimum_sampled_centroid_speed_um_s=float((-np.gradient(c, t)).min()),
                       strength_error_final500=strength_error,
                       native_field_tail_pass=meta['metrics']['field_tail_pass'],
                       source_archive=str(path.relative_to(ROOT)),
                       checkpoint_sha256=case['checkpoint_sha256'])
            rows.append(row)
            checks.append(dict(condition=geo, length_um=L,
                               minimum_strength_duration_verified=True,
                               archive_hash_verified=True, checkpoint_hash_verified=True,
                               checkpoint_and_withdrawal_continuity_verified=True,
                               scalar_gates_recomputed=True,
                               centroid_moment_identity_verified=True,
                               reflected_target_identity_verified=True,
                               maximum_sample_interval_s=float(np.diff(t).max())))
    assert len(rows) == 36
    for geo in GEOS:
        assert sorted(r['length_um'] for r in rows if r['condition'] == geo) == LENGTHS
    return rows, checks


def render(rows):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.2), sharex=True, sharey='row')
    fig.subplots_adjust(left=.078, right=.985, top=.81, bottom=.17, hspace=.38, wspace=.14)
    fig.text(.078, .972, 'Figure Z | Bound-MCA redistribution across cell lengths',
             fontsize=17, weight='bold', va='top')
    fig.text(.078, .927, 'The redistribution distance grows; cue-average centroid speed does not keep pace.',
             fontsize=11.5, va='top', color='#333333')
    distance_styles = [
        ('initial_distance_to_reflected_state_um', '#555555', '--', 's', 'Initial → reflected polarized centroid'),
        ('initial_distance_to_midpoint_um', '#8c6bb1', ':', '^', 'Initial → midpoint'),
        ('centroid_shift_during_cue_um', '#0072B2', '-', 'o', 'Shift during selected cue')]
    speed_styles = [
        ('early_mean_speed_um_s', '#009E73', '-', '^', 'Early: first 60 s'),
        ('late_mean_speed_um_s', '#D55E00', '-', 's', 'Late: final 60 s before withdrawal'),
        ('cue_average_speed_um_s', '#0072B2', '-', 'o', 'Cue average: onset → withdrawal')]
    for col, geo in enumerate(GEOS):
        rr = sorted((r for r in rows if r['condition'] == geo), key=lambda r: r['length_um'])
        lengths = np.array([r['length_um'] for r in rr])
        ref = next(r for r in rr if r['length_um'] == 23)
        for row, styles in enumerate([distance_styles, speed_styles]):
            ax = axes[row, col]
            ax.text(-.12, 1.035, chr(65 + row * 3 + col), transform=ax.transAxes,
                    fontsize=13, weight='bold')
            ax.axhline(0, color='#aaaaaa', lw=.7)
            ax.grid(axis='y', alpha=.16)
            for key, color, ls, marker, label in styles:
                values = np.array([r[key] for r in rr]) * (1000 if row else 1)
                ax.plot(lengths, values, color=color, ls=ls, marker=marker,
                        ms=4.3, lw=1.7, label=label)
            ax.set(xlim=(17, 40), xticks=[18, 23, 27, 31, 35, 39])
            if row == 0:
                ax.set_title(TITLES[col], fontsize=12, pad=15)
            else:
                ax.plot(lengths, ref['cue_average_speed_um_s'] * lengths / 23 * 1000,
                        color='#777777', ls=(0, (5, 3)), lw=1.2, zorder=1)
                ax.set_xlabel('Cell length (µm)')
                last = rr[-1]
                ratio = last['cue_average_speed_um_s'] / ref['cue_average_speed_um_s']
                ax.text(.04, .965, f'23 → 39 µm: distance ×{last["initial_distance_to_reflected_state_um"] / ref["initial_distance_to_reflected_state_um"]:.2f}'
                        f'\ncue-average speed ×{ratio:.2f}', transform=ax.transAxes,
                        va='top', fontsize=9, color='#333333',
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=2))
        axes[0, col].set_ylim(bottom=0)
    axes[0, 0].set_ylabel('Centroid distance / net shift (µm)')
    axes[1, 0].set_ylabel('Centroid speed toward new rear\n(10⁻³ µm/s)')
    fig.legend(handles=[Line2D([], [], color=c, ls=ls, marker=m, ms=4, label=label)
                        for _, c, ls, m, label in distance_styles],
               loc='upper center', bbox_to_anchor=(.53, .889), ncol=3,
               frameon=False, fontsize=9.3, columnspacing=2.1)
    speed_handles = [Line2D([], [], color=c, ls=ls, marker=m, ms=4, label=label)
                     for _, c, ls, m, label in speed_styles]
    speed_handles.append(Line2D([], [], color='#777777', ls=(0, (5, 3)), lw=1.2,
                               label='Speed proportional to length (cue average at 23 µm)'))
    fig.legend(handles=speed_handles, loc='lower center', bbox_to_anchor=(.53, .055),
               ncol=2, frameon=False, fontsize=9.3, columnspacing=2.6)
    fig.text(.078, .022,
             'Shortest tested 1%-strength-passing cues · nominal coordinates · centroid motion includes redistribution and exchange',
             fontsize=9, color='#555555')
    paths = []
    for ext in ['png', 'pdf', 'svg']:
        p = output(f'{STEM}.{ext}')
        fig.savefig(p, dpi=220, facecolor='white')
        paths.append(p)
    plt.close(fig)
    return paths


def main():
    register(Path(__file__))
    register(Path(layout.__file__))
    rows, checks = measure()
    table = output(f'{STEM}.csv')
    with table.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    paths = render(rows)
    caption = (
        'Figure Z. Bound-MCA redistribution distance and centroid speed across the twelve eligible corrected deBelly lengths. '
        'Columns compare a fixed 3-µm patch, a relative 20% patch, and constant physical slope 0.24. '
        'Every trajectory begins at its exact corrected relative-20% basal checkpoint with independent endpoint commands. '
        'A–C: net centroid shift toward the new rear during the shortest tested canonical 1%-strength-passing finite cue (blue); '
        'initial centroid distance to the nominal midpoint (purple); and initial centroid separation from the exact spatially reflected '
        'pre-cue bound-MCA profile (gray). The latter two curves are initial geometric benchmarks, not remaining distances at withdrawal. '
        'D–F: signed centroid displacement divided by elapsed time over the first 60 s (green), final 60 s before withdrawal (orange), '
        'and entire cue (blue). All windows lie within the cue-on phase, with linearly interpolated window boundaries; the full-cue mean '
        'equals the integral of signed centroid velocity divided by cue duration. Dashed gray speed curves show proportional-to-length '
        'scaling anchored to each geometry’s cue-average speed at L=23 µm. Insets compare L=23 and 39; the reflected-centroid distance '
        'increases 1.83-fold while cue-average speed changes by 0.90-, 0.55-, and 1.09-fold for fixed, relative, and slope cues respectively. '
        'These endpoint comparisons do not assume monotonic intermediate trends. Early and late rates differ, and selected '
        'cue durations vary with length; cue-average speed is a descriptive readout, not an independent predictor of selected duration. '
        'Centroids use the mass-normalized bound-MCA distribution in nominal reference coordinates; these speeds include local exchange '
        'and changes in bound amount, and do not measure tracked molecular travel, cortical flow, or laboratory-frame motion. '
        'The reflected state is a geometric reference, not an assertion that MCA has reached that profile at the scalar Rac reversal gate. '
        'Separate native-field stability flags, all cue durations, centroid path lengths and remaining distances at withdrawal are retained '
        'in the CSV; no length is removed because of its field-tail flag. Lines connect sampled lengths without a fitted scaling law.'
    )
    caption_path = output(f'{STEM}_caption.txt')
    caption_path.write_text(caption + '\n')
    summaries = {}
    for geo in GEOS:
        rr = [r for r in rows if r['condition'] == geo]
        short, long = (next(r for r in rr if r['length_um'] == L) for L in [23, 39])
        summaries[geo] = dict(
            length_ratio_39_over_23=39 / 23,
            target_distance_ratio_39_over_23=long['initial_distance_to_reflected_state_um'] / short['initial_distance_to_reflected_state_um'],
            **{key + '_ratio_39_over_23': long[key] / short[key] for key in
               ['early_mean_speed_um_s', 'late_mean_speed_um_s', 'cue_average_speed_um_s']},
            native_field_tail_passes=sum(r['native_field_tail_pass'] for r in rr))
    audit_name = f'{STEM}_audit.json'
    write(audit_name, dict(passed=True, generated_at=datetime.now(timezone.utc).isoformat(),
                           trials=36, checks=checks, new_integrations=0, summary=summaries))
    paths += [table, caption_path, output(audit_name)]
    write(f'{STEM}.manifest.json', dict(
        artifact_root='../..', producer='finished', generated_at=datetime.now(timezone.utc).isoformat(),
        passed=True, legend=caption, description=caption,
        sha256={layout.dev_artifact_key(OUT, p): sha(p) for p in paths},
        source_sha256=SOURCES, renderer_sha256=sha(Path(__file__)),
        methods=dict(centroid='sum(x * bound_MCA) / sum(bound_MCA), uniform finite volumes',
                     target='Exact pre-cue profile reflected on nominal reference coordinates',
                     positive_direction='Left, toward the new rear',
                     early_window_s=[0, WINDOW], late_window_s='[D-60, D]',
                     integrated_average='[centroid(0)-centroid(D)] / D; signed displacement, not path length',
                     coordinates='Nominal reference finite-volume centers',
                     cue_selection='Shortest tested scalar 1%-strength-passing duration; independent field-tail flag'),
        native_histories=36, new_integrations=0, current_report_changed=False,
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '
                      '/home/4473331/.conda/envs/scaling_sandbox/bin/python ' + str(Path(__file__))]))
    print(json.dumps(summaries, indent=2))
    print(output(f'{STEM}.png'))


if __name__ == '__main__':
    main()
