"""Figure X: annotated native deBelly reversal kymographs at three lengths.

Uses the relative-20% geometry of Figure 11S2 and the current shortest tested
scalar-strength-passing pulses. No integrations or manuscript-report writes.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / 'figures/length_mechanisms_20261008'
OUT = REPORT / 'dev/figure_x_reversal_overview'
RUN = ROOT / 'data/runs/length_narrative_memory_20261009'
BASE = ROOT / 'data/runs/debelly_corrected_three_test_20261008'
STEM = 'Figure_X_reversal_overview'
LENGTHS = [23, 31, 39]
FIELDS = ['rac', 'rho', 'mca_bound', 'velocity', 'tension']
LABELS = ['Rac', 'Rho', 'Bound MCA', 'Cortical flow', 'Membrane tension']
UNITS = ['Native conc.', 'Native conc.', 'Native conc.', 'µm/s', 'pN']
LIMIT = 1200.
PRE_CUE = 60.
SOURCES = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def register(p):
    SOURCES[str(p.relative_to(ROOT))] = sha(p)
    return p


def read(p):
    return json.loads(register(p).read_text())


def output(name):
    return layout.dev_output_path(OUT, name)


def write(name, value):
    output(name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def crossings(t, y):
    indices = np.flatnonzero(y[:-1] * y[1:] < 0)
    return [float(t[i] + (t[i + 1] - t[i]) * (-y[i]) / (y[i + 1] - y[i]))
            for i in indices]


def load():
    manifest = read(RUN / 'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(RUN / 'status.json')['audit_passed']
    assert read(BASE / 'results/independent_audit.json')['passed']
    records_path = register(BASE / 'results/records.json.gz')
    with gzip.open(records_path, 'rt') as stream:
        records = {r['case_index']: r for r in json.load(stream)['cases']}
    traces, measurements, audits, arrays = {}, [], [], {}
    for L in LENGTHS:
        job = next(j for j in manifest['jobs'] if j['case']['length'] == L
                   and j['case']['condition'] == 'fixed_relative')
        case, D = job['case'], float(job['duration'])
        record = records[case['case_index']]
        assert D == record['result']['full_polarity']['endpoint']
        assert D == min(trial['duration'] for trial in record['tests'] if trial['full_polarity_success'])
        p = ROOT / job['source'] if job['source'] else RUN / f"workers/{manifest['missing_jobs'].index(job):03d}.npz"
        register(p)
        meta = read(p.with_suffix('.json'))
        assert sha(p) == meta['archive_sha256']
        cp = register(ROOT / case['checkpoint'])
        assert sha(cp) == case['checkpoint_sha256']
        with np.load(cp) as z:
            initial = {key: z['state_' + key].copy() for key in FIELDS}
        with np.load(p) as z:
            x = z['x'].copy()
            tp, tr = z['pulse_time'], z['recovery_time']
            assert abs(tp[-1] - D) < 1e-9 and tr[-1] == 10000
            np.testing.assert_allclose(x, (np.arange(len(x)) + .5) * L / len(x), atol=1e-12)
            t = np.r_[tp, D + tr[1:]]
            assert t[0] == 0 and np.all(np.diff(t) > 0)
            stop = int(np.searchsorted(t, LIMIT)) + 1
            fields = {}
            for key in FIELDS:
                a, b = z['pulse_field_' + key], z['recovery_field_' + key]
                np.testing.assert_allclose(a[0], initial[key], atol=1e-12)
                np.testing.assert_allclose(b[0], a[-1], atol=1e-12)
                fields[key] = np.r_[a, b[1:]][:stop]
                assert np.isfinite(fields[key]).all()
                if key not in ['velocity', 'tension']:
                    assert fields[key].min() >= -1e-10
            tail = tr >= 9500
            r = z['recovery_field_rac'][tail]
            polarity = r @ (2 * x / L - 1) / r.sum(axis=-1)
            contrast = np.ptp(r, axis=1)
            strength_error = float(np.max(abs(polarity / abs(case['initial_orientation']) - 1)))
            assert polarity.min() > .05 and contrast.min() >= .05 and np.ptp(polarity) < .005
            assert strength_error <= .01 and meta['metrics']['full_polarity_success']
        time = t[:stop]
        pre_time = np.arange(-PRE_CUE, 0., .5)
        diag = {}
        for key in FIELDS:
            values = fields[key]
            # These descriptive summaries are distinct for chemical bias,
            # mean flow and the signed spatial moment of tension.
            if key == 'velocity':
                y = values.mean(axis=1)
            else:
                y = values @ (2 * x / L - 1) / np.abs(values).sum(axis=1)
                if key in ['rho', 'mca_bound', 'tension']:
                    y = -y
            zeros = crossings(time, y)
            diag[key] = dict(first_zero_s=zeros[0] if zeros else None,
                             all_zeros_in_display_s=zeros,
                             zero_crossings_during_cue_s=[v for v in zeros if v <= D])
            measurements.append(dict(length_um=L, field=key, cue_duration_s=D,
                diagnostic=('Mean signed velocity; positive to right' if key == 'velocity' else
                            'Signed normalized spatial first moment; positive to right' if key == 'rac' else
                            'Negative signed normalized spatial first moment; positive to left'),
                initial=float(y[0]), at_60_s=float(np.interp(60, time, y)),
                at_withdrawal=float(np.interp(D, time, y)),
                at_1200_s=float(np.interp(LIMIT, time, y)),
                first_zero_s=zeros[0] if zeros else None,
                all_zeros_in_display_s=json.dumps(zeros),
                native_field_tail_pass=meta['metrics']['field_tail_pass']))
            plateau = np.broadcast_to(initial[key], (len(pre_time), len(x)))
            arrays[f'L{L}_{key}'] = np.concatenate([plateau, values])
            np.testing.assert_array_equal(arrays[f'L{L}_{key}'][:len(pre_time)], plateau)
        display_time = np.r_[pre_time, time]
        arrays[f'L{L}_time'] = display_time
        arrays[f'L{L}_x'] = x
        arrays[f'L{L}_cue_duration'] = np.array(D)
        traces[L] = dict(x=x, time=display_time,
                         fields={key: arrays[f'L{L}_{key}'] for key in FIELDS},
                         duration=D, diagnostics=diag)
        audits.append(dict(length_um=L, archive=str(p.relative_to(ROOT)),
                           archive_hash_verified=True, checkpoint_sha256=case['checkpoint_sha256'],
                           exact_checkpoint_verified_for=FIELDS, phase_continuity_verified_for=FIELDS,
                           shortest_tested_strength_selection_verified=True,
                           independent_scalar_gate_pass=True, maximum_strength_error=strength_error,
                           native_field_tail_pass=meta['metrics']['field_tail_pass'],
                           retained_cue_output_intervals_s=sorted(set(np.round(np.diff(tp), 9))),
                           recovery_output_interval_s=2.5, diagnostics=diag,
                           pre_cue_display_s=[-PRE_CUE, 0],
                           pre_cue_initial_checkpoint_repeated_exactly=True))
    return traces, measurements, audits, arrays


def marker(ax, number, xy, offset):
    ax.annotate(str(number), xy=xy, xytext=offset, textcoords='offset points',
                ha='center', va='center', fontsize=9, weight='bold',
                bbox=dict(boxstyle='circle,pad=.25', facecolor='white', edgecolor='#222222'),
                arrowprops=dict(arrowstyle='->', color='white', lw=2.8), zorder=6)
    ax.annotate('', xy=xy, xytext=offset, textcoords='offset points',
                arrowprops=dict(arrowstyle='->', color='#222222', lw=.9), zorder=7)


def render(traces):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(11.8, 12.8))
    gs = fig.add_gridspec(5, 4, width_ratios=[23, 31, 39, 1.9],
                          left=.10, right=.93, bottom=.065, top=.965,
                          hspace=.24, wspace=.18)
    axes = {}
    for row, (key, label, unit) in enumerate(zip(FIELDS, LABELS, UNITS)):
        visible = [tr['fields'][key][tr['time'] <= LIMIT] for tr in traces.values()]
        vmin = min(float(a.min()) for a in visible)
        vmax = max(float(a.max()) for a in visible)
        if key == 'rac':
            # Allocate half the palette to the basal concentration range,
            # retaining the full cue-induced range in the upper half.
            norm, cmap = TwoSlopeNorm(vmin=vmin, vcenter=1., vmax=vmax), 'viridis'
        elif key in ['velocity', 'tension']:
            maximum = max(abs(vmin), abs(vmax))
            norm, cmap = Normalize(-maximum, maximum), 'RdBu_r'
        else:
            norm, cmap = Normalize(vmin, vmax), 'viridis'
        for col, L in enumerate(LENGTHS):
            ax = fig.add_subplot(gs[row, col])
            tr = traces[L]
            mesh = ax.pcolormesh(tr['x'], tr['time'], tr['fields'][key],
                                shading='nearest', cmap=cmap, norm=norm, rasterized=True)
            ax.set_yscale('symlog', linthresh=120, linscale=.65)
            ax.set(xlim=(0, L), ylim=(LIMIT, -PRE_CUE), xticks=[0, L / 2, L],
                   yticks=[-60, 0, 60, 120, 300, 600, 1200])
            ax.set_yticklabels(['−60', '0', '60', '120', '300', '600', '1200'])
            ax.tick_params(labelsize=9, direction='out', length=3)
            if col:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(label + '\nTime from cue onset (s)', fontsize=10, labelpad=8)
            ax.axhline(0, color='white', lw=1.5)
            ax.axhline(0, color='#444444', lw=.6)
            ax.axhline(60, color='white', ls=':', lw=1.5, alpha=.9)
            ax.axhline(tr['duration'], color='white', ls='--', lw=2)
            ax.axhline(tr['duration'], color='#222222', ls='--', lw=.7)
            ax.text(.025, .975, chr(65 + row * 3 + col), transform=ax.transAxes,
                    ha='left', va='top', fontsize=11, weight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=1.6))
            if row == 0:
                ax.set_title(f'L = {L} µm', fontsize=12, pad=12)
            if row == 4:
                ax.set_xlabel('Reference position x (µm)', fontsize=10)
            axes[row, L] = ax
        cax = fig.add_subplot(gs[row, 3])
        cb = fig.colorbar(mesh, cax=cax)
        if key == 'rac':
            cb.set_ticks([.6, .8, 1., 2., 4.])
        cb.ax.tick_params(labelsize=8, length=2)
        cb.set_label(unit, fontsize=9, labelpad=6)
        cb.ax.yaxis.get_offset_text().set_fontsize(8)
    marker(axes[0, 31], 1, (29, 30), (-35, -10))
    marker(axes[1, 31], 2, (26, 60), (-35, -10))
    marker(axes[2, 31], 3, (24, 60), (-25, -18))
    marker(axes[3, 39], 4, (29, 485), (-25, 28))
    marker(axes[4, 39], 5, (10, 530), (25, 28))
    paths = []
    for ext in ['png', 'pdf', 'svg']:
        p = output(f'{STEM}.{ext}')
        fig.savefig(p, dpi=200, facecolor='white')
        paths.append(p)
    plt.close(fig)
    return paths


def refresh_redraft_collection():
    """Keep the user's collected Figure X identical to the latest export."""
    destination = REPORT / 'dev/deBelly_redraft'
    collection_path = destination / 'supporting/metadata/deBelly_redraft_collection.manifest.json'
    if not collection_path.exists():
        return
    collection = json.loads(collection_path.read_text())
    manifest_path = output(f'{STEM}.manifest.json')
    manifest = json.loads(manifest_path.read_text())
    files = dict(manifest['sha256'])
    files[layout.dev_artifact_key(OUT, manifest_path)] = sha(manifest_path)
    collection['copies'] = [row for row in collection['copies'] if row['figure'] != 'X']
    for relative, digest in files.items():
        source, target = OUT / relative, destination / relative
        assert sha(source) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        assert sha(target) == digest
        collection['sha256'][relative] = digest
        collection['copies'].append(dict(figure='X', source=str(source.relative_to(ROOT)),
                                         destination=relative, sha256=digest))
    selected = collection['selected_figures']['X']
    selected.update(source_manifest=str(manifest_path.relative_to(ROOT)),
                    source_manifest_sha256=sha(manifest_path),
                    generated_at=manifest['generated_at'],
                    copied_manifest=layout.dev_artifact_key(OUT, manifest_path),
                    candidate_versions=[dict(manifest=str(manifest_path.relative_to(ROOT)),
                                             generated_at=manifest['generated_at'])])
    collection['generated_at'] = datetime.now(timezone.utc).isoformat()
    for relative, digest in collection['sha256'].items():
        assert sha(destination / relative) == digest
    collection_path.write_text(json.dumps(collection, indent=2, ensure_ascii=False) + '\n')


def main():
    before = {name: sha(REPORT / name) for name in ['index.html', 'length_mechanisms.pdf']}
    register(Path(__file__))
    register(Path(layout.__file__))
    traces, rows, checks, arrays = load()
    table = output(f'{STEM}_measurements.csv')
    with table.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fields_path = output(f'{STEM}_fields.npz')
    np.savez_compressed(fields_path, **arrays)
    paths = render(traces)
    caption = (
        'Figure X. Observable deBelly reversal dynamics precede the detailed mechanism analysis. Columns show L=23, 31 and 39 µm; '
        'rows show Rac (A–C), Rho (D–F), bound MCA (G–I), cortical velocity (J–L) and membrane tension (M–O). '
        'Each cell begins at its exact qualifying corrected relative-20% basal checkpoint, with independent endpoint commands. '
        'The right raised-cosine patch has peak 3.6 and support 20% of nominal length (4.6, 6.2 and 7.8 µm). '
        'The shortest tested canonical scalar 1%-strength-passing cues last 232.244898, 285.280612 and 809.081633 s, respectively. '
        'The solid line at time zero marks cue onset; dashed lines mark withdrawal; dotted lines mark the common 60-s comparison. '
        'The interval from −60 to 0 s displays the verified initial basal checkpoint repeated without change for visibility; it is not an additional simulated pre-cue history. '
        'No chemical or mechanical field is reset at withdrawal. '
        'Time increases downward and is linear from −60 through 120 s, then logarithmic, on a shared symlog axis through 1200 s after cue onset. '
        'Column widths scale with nominal cell length; positions are undeformed reference coordinates, not membrane positions in the laboratory frame. '
        'All recorded samples through the visible interval are retained, with one common absolute color scale per field across the three lengths. '
        'Chemical colors report native concentrations. Rac uses a continuous piecewise-linear color mapping: half the palette spans the displayed minimum to 1.0, '
        'and half spans 1.0 to the displayed maximum, revealing basal polarity while retaining the full cue-induced range without clipping. '
        'The 1.0 breakpoint is a display choice, not a biological threshold; colorbar ticks give native concentrations. '
        'Rho and bound MCA use linear scales; velocity and tension use zero-centered signed linear scales. '
        'Numbered arrows identify (1) early Rac gain at the challenged end, (2) early Rho reorientation, (3) persisting old-rear MCA bias, '
        '(4) delayed long-cell whole-cortex flow reversal, and (5) delayed long-cell tension redistribution. '
        'At 60 s, Rac and Rho first moments have reversed at all three lengths while bound-MCA centroids have not. '
        'The 39-µm case illustrates a mechanical exception to uniform fast reorientation: its mean velocity changes sign near 485 s and its tension spatial '
        'moment near 530 s. Bound MCA can cross and recross the midpoint; the overview does not assign a universal MCA reversal threshold. '
        'Rac, Rho and MCA diagnostic first moments are mass-normalized. Tension bias is its signed first moment divided by the spatial sum of absolute '
        'tension; flow timing refers to the spatial mean of signed cortical velocity. These descriptive sign changes are distinct observables, '
        'not interchangeable measures of settling or successful reversal. The weakening after withdrawal does not imply failed reversal: '
        'all three selected pulses satisfy the final-500-s scalar Rac gates after 10000 s of basal recovery, beyond this figure’s display. '
        'Separate native-field tail flags are retained in the measurements and audit. These histories motivate examination of MCA as a candidate slow '
        'memory variable; they do not establish a causal attribution. No new integrations were performed.'
    )
    cp = output(f'{STEM}_caption.txt')
    cp.write_text(caption + '\n')
    assert before == {name: sha(REPORT / name) for name in before}
    audit_name = f'{STEM}_audit.json'
    write(audit_name, dict(passed=True, generated_at=datetime.now(timezone.utc).isoformat(),
                           new_integrations=0, checks=checks, current_report_hashes=before,
                           current_report_unchanged=True))
    paths += [table, fields_path, cp, output(audit_name)]
    write(f'{STEM}.manifest.json', dict(
        artifact_root='../..', passed=True, producer='finished', generated_at=datetime.now(timezone.utc).isoformat(),
        description=caption, legend=caption,
        sha256={layout.dev_artifact_key(OUT, p): sha(p) for p in paths},
        source_sha256=SOURCES, renderer_sha256=sha(Path(__file__)),
        methods=dict(lengths_um=LENGTHS, condition='fixed_relative', cue_peak=3.6,
                     fields=FIELDS, displayed_time_s=[-PRE_CUE, LIMIT], time_axis='Symlog: linthresh=120 s, linscale=0.65',
                     pre_cue_display='Exact initial basal checkpoint repeated from -60 to 0 s; no added integration',
                     in_figure_text='Axis labels, units, lengths, panel letters and numbered arrows only',
                     color_scales='One common native concentration scale across lengths for each field; Rac piecewise linear with half the palette below 1.0 and half above; remaining fields linear, symmetric around zero for velocity/tension',
                     rac_color_breakpoint=1., rac_colorbar_ticks=[.6, .8, 1., 2., 4.],
                     annotations='Observed native field features; diagnostic crossings linearly interpolated between saved samples'),
        new_integrations=0, current_report_changed=False,
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '
                      '/home/4473331/.conda/envs/scaling_sandbox/bin/python ' + str(Path(__file__))]))
    refresh_redraft_collection()
    print(output(f'{STEM}.png'))
    print(json.dumps({L: traces[L]['diagnostics'] for L in LENGTHS}, indent=2))


if __name__ == '__main__':
    main()
