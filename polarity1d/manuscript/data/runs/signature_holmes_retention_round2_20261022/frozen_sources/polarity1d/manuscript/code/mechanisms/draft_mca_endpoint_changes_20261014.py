"""Native endpoint bound-MCA changes across length, during the reversing cue."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import ScalarFormatter

import draft_mca_length_comparison_20261014 as common
from length_dev_output_layout import dev_output_path

ROOT, REPORT, RUN = common.ROOT, common.REPORT, common.BUDGET
OUT = REPORT / 'dev/mca_endpoint_change_drafts_20261014'
TIMES = [30., 60., 120., 180.]
END_LABELS = ['Left endpoint (old front → new rear)', 'Right endpoint (cue end → new front)']
SOURCES = {}


def source(path):
    SOURCES[str(path.relative_to(ROOT))] = common.sha(path)
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


def save(fig, stem):
    for ext in ['.png', '.pdf', '.svg']:
        fig.savefig(output(stem + ext), dpi=180, facecolor='white')
    plt.close(fig)
    return stem


def setup(title, subtitle, height=5.2):
    fig, axes = plt.subplots(1, 2, figsize=(12.6, height), sharey=True)
    fig.subplots_adjust(left=.085, right=.99, bottom=.20, top=.68, wspace=.17)
    fig.text(.085, .95, 'DRAFT · ' + title, fontsize=16, weight='bold')
    fig.text(.085, .885, subtitle, fontsize=11)
    for ax, label in zip(axes, END_LABELS):
        ax.set_title(label, fontsize=12, pad=60)
        ax.set_xlabel('Nominal cell length (µm)')
        ax.set_xticks([18, 23, 27, 31, 35, 39])
        ax.set_xlim(17, 40)
        common.clean_axes(ax)
        ax.axhline(0, color='#aaaaaa', lw=.8)
    return fig, axes


def load():
    manifest = read(RUN / 'manifest.json')
    assert read(RUN / 'status.json')['audit_passed']
    assert manifest['geometry'] == 'fixed_relative'
    traces, rows, summary, audits = {}, [], [], []
    for slot, length in enumerate(manifest['lengths_um']):
        path = RUN / 'workers' / f'{slot:03d}.npz'
        meta = read(path.with_suffix('.json'))
        assert meta['passed'] and meta['parameters']['center_boundary_motion'] is False
        assert common.sha(source(path)) == meta['npz_sha256']
        check = meta['checks'][0]
        checkpoint = source(ROOT / check['checkpoint'])
        with np.load(path) as z:
            time, x = z[f'L{length}_time'].copy(), z[f'L{length}_x'].copy()
            bound = z[f'L{length}_bound'].copy()
        with np.load(checkpoint) as z:
            np.testing.assert_array_equal(bound[0], z['state_mca_bound'])
        duration = check['cue_duration_s']
        assert time[0] == 0 and time[-1] == duration and duration > max(TIMES)
        assert np.all(np.diff(time) > 0) and np.isfinite(bound).all()
        np.testing.assert_allclose(x[[0, -1]], [length/402, length-length/402], rtol=0, atol=1e-12)
        endpoint = bound[:, [0, -1]]
        delta = endpoint - endpoint[0]
        rate = np.diff(endpoint, axis=0) / np.diff(time)[:, None]
        np.testing.assert_allclose((rate * np.diff(time)[:, None]).sum(axis=0), delta[-1], rtol=0, atol=1e-14)
        traces[length] = dict(time=time, endpoint=endpoint, delta=delta, rate=rate, duration=duration)
        for i, t in enumerate(time):
            rows.append(dict(length_um=length, elapsed_time_s=float(t), cue_duration_s=duration,
                             left_x_um=float(x[0]), right_x_um=float(x[-1]),
                             left_bound_MCA=float(endpoint[i, 0]), right_bound_MCA=float(endpoint[i, 1]),
                             left_net_change_conc=float(delta[i, 0]), right_net_change_conc=float(delta[i, 1]),
                             left_net_rate_conc_s=float(rate[i-1, 0]) if i else '',
                             right_net_rate_conc_s=float(rate[i-1, 1]) if i else ''))
        for label, t in [(f'{int(t)}_s', t) for t in TIMES] + [('cue_end', duration)]:
            index = np.flatnonzero(np.isclose(time, t, rtol=0, atol=1e-10))
            assert len(index) == 1
            change = delta[index[0]]
            summary.append(dict(length_um=length, sample=label, elapsed_time_s=float(t), cue_duration_s=duration,
                                left_net_change_conc=float(change[0]), right_net_change_conc=float(change[1]),
                                left_mean_net_rate_conc_s=float(change[0]/t),
                                right_mean_net_rate_conc_s=float(change[1]/t)))
        audits.append(dict(length_um=length, source_hash_verified=True, exact_checkpoint_verified=True,
                           endpoint_native_cell_centres_um=x[[0, -1]].tolist(),
                           native_steps=len(time)-1, time_integration_matches_net_change=True,
                           dense_replay_verified_original_error=check['dense_replay_max_error']))
    return traces, rows, summary, audits


def snapshots(summary):
    fig, axes = setup('Endpoint bound-MCA change versus cell length',
                      'Net change = bound MCA at elapsed time minus bound MCA before the cue; relative-20% right cue.', 5.7)
    colors = plt.get_cmap('viridis')(np.linspace(.12, .85, len(TIMES)))
    for ax, endpoint in zip(axes, ['left', 'right']):
        for t, color in zip(TIMES, colors):
            selected = [r for r in summary if r['sample'] == f'{int(t)}_s']
            ax.plot([r['length_um'] for r in selected], [r[f'{endpoint}_net_change_conc'] for r in selected],
                    color=color, lw=1.8, marker='o', ms=4, label=f'{int(t)} s')
        selected = [r for r in summary if r['sample'] == 'cue_end']
        ax.plot([r['length_um'] for r in selected], [r[f'{endpoint}_net_change_conc'] for r in selected],
                color='#202020', ls='--', lw=1.8, marker='s', ms=4, label='Own cue end')
        ax.legend(title='Time from cue onset', loc='lower left', bbox_to_anchor=(0, 1.015), frameon=False,
                  ncol=3, fontsize=10, title_fontsize=10, columnspacing=1.1, handlelength=1.6, borderaxespad=0)
    axes[0].set_ylabel('Bound-MCA net change (native concentration)')
    fig.text(.085, .065, 'Common-time curves compare equal elapsed times; dashed curves use length-dependent cue ends (226–809 s).', fontsize=10)
    return save(fig, 'endpoint_bound_MCA_net_change_vs_length')


def means(summary):
    fig, axes = setup('Mean endpoint bound-MCA net change rate',
                      'Observed net rate averaged over the window; includes every contribution to bound MCA.')
    for ax, endpoint in zip(axes, ['left', 'right']):
        for sample, color, style, label in [('cue_end', '#202020', '--', 'Whole selected cue'),
                                           ('180_s', '#007F7F', '-', 'Common first 180 s')]:
            selected = [r for r in summary if r['sample'] == sample]
            ax.plot([r['length_um'] for r in selected], [r[f'{endpoint}_mean_net_rate_conc_s'] for r in selected],
                    color=color, ls=style, marker='o', ms=4, lw=1.8, label=label)
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.015), frameon=False, fontsize=10,
                  borderaxespad=0)
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_powerlimits((-4, -4))
        ax.yaxis.set_major_formatter(formatter)
    axes[0].set_ylabel('Mean bound-MCA net rate (native concentration / s)')
    fig.text(.085, .065, 'Mean net rate = (bound MCA at window end − starting bound MCA) / window duration. Positive means accumulation.', fontsize=10)
    return save(fig, 'endpoint_bound_MCA_mean_net_rate_vs_length')


def evolution(traces):
    lengths = list(traces)
    time_centres = np.arange(.5, max(v['duration'] for v in traces.values()), 1.)
    time_edges = np.r_[time_centres - .5, time_centres[-1] + .5]
    length_edges = np.r_[lengths[0]-(lengths[1]-lengths[0])/2,
                         (np.array(lengths[:-1])+lengths[1:])/2,
                         lengths[-1]+(lengths[-1]-lengths[-2])/2]
    values = np.full((2, len(time_centres), len(lengths)), np.nan)
    for column, length in enumerate(lengths):
        tr = traces[length]
        selected = time_centres <= tr['duration']
        for end in range(2):
            values[end, selected, column] = np.interp(time_centres[selected], tr['time'], tr['delta'][:, end])
    maximum = max(float(abs(v['delta']).max()) for v in traces.values())
    cmap = LinearSegmentedColormap.from_list('endpoint_change', ['#4C78A8', '#ffffff', '#E45756'])
    cmap.set_bad('#e4e4e4')
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.5), sharey=True)
    fig.subplots_adjust(left=.085, right=.87, bottom=.16, top=.76, wspace=.17)
    fig.text(.085, .95, 'DRAFT · Endpoint bound-MCA change throughout the cue', fontsize=16, weight='bold')
    fig.text(.085, .895, 'Change from each endpoint’s starting value · blue: depletion · red: accumulation', fontsize=11)
    for end, (ax, label) in enumerate(zip(axes, END_LABELS)):
        field = ax.pcolormesh(length_edges, time_edges, np.ma.masked_invalid(values[end]),
                             cmap=cmap, norm=TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum),
                             shading='flat', rasterized=True)
        ax.plot(lengths, [traces[length]['duration'] for length in lengths], color='#202020',
                marker='s', ms=4, ls='', label='Cue withdrawal')
        ax.set_title(label, fontsize=12, pad=12)
        ax.set_xlabel('Nominal cell length (µm)')
        ax.set_xticks(lengths)
        ax.tick_params(axis='x', labelsize=9)
        ax.set_ylim(0, max(time_edges))
        ax.legend(loc='upper left', fontsize=10, frameon=True, facecolor='white', framealpha=.95)
    axes[0].set_ylabel('Time from cue onset (s)')
    cax = fig.add_axes([.90, .16, .018, .60])
    fig.colorbar(field, cax=cax, label='Bound-MCA net change (native concentration)')
    fig.text(.085, .065, 'Each stripe is one simulated length. Gray is after cue withdrawal; no recovery is shown. Both ends share one color scale.', fontsize=10)
    return save(fig, 'endpoint_bound_MCA_full_cue_evolution')


def cumulative_outward_loss(traces):
    summary, intervals, audits = [], [], []
    for slot, (length, trace) in enumerate(traces.items()):
        base = RUN / 'workers' / f'{slot:03d}'
        meta = read(base.with_suffix('.json'))
        path = source(base.with_suffix('.csv'))
        assert common.sha(path) == meta['csv_sha256']
        with path.open() as stream:
            native = list(csv.DictReader(stream))
        start = np.array([float(r['time_start_s']) for r in native])
        end = np.array([float(r['time_end_s']) for r in native])
        dt = np.array([float(r['step_s']) for r in native])
        left = np.array([float(r['left_bound_outflux']) for r in native])
        right = np.array([float(r['right_bound_outflux']) for r in native])
        assert len(native) == len(trace['time'])-1
        assert np.min(left) >= 0 and np.min(right) >= 0
        np.testing.assert_array_equal(start, trace['time'][:-1])
        np.testing.assert_array_equal(end, trace['time'][1:])
        np.testing.assert_allclose(end-start, dt, rtol=0, atol=1e-12)
        # Endpoint faces use the implicit bound field before uniform closure,
        # and the actual commands for that native step, not stored interior v.
        uniform = np.array([float(r['bound_uniform_source_conc_s']) for r in native])
        preclosure = trace['endpoint'][1:] - (uniform*dt)[:, None]
        commands = np.array([[float(r['left_command_um_s']), float(r['right_command_um_s'])]
                             for r in native])
        reconstructed = commands * preclosure * np.array([-1., 1.])
        error = float(np.max(abs(reconstructed - np.column_stack([left, right]))))
        assert error < 1e-12
        cumulative_left, cumulative_right = np.cumsum(left*dt), np.cumsum(right*dt)
        assert np.all(np.diff(cumulative_left) >= 0) and np.all(np.diff(cumulative_right) >= 0)
        for i in range(len(native)):
            intervals.append(dict(length_um=length, time_start_s=float(start[i]), time_end_s=float(end[i]),
                step_s=float(dt[i]), left_bound_outflux_amount_s=float(left[i]),
                right_bound_outflux_amount_s=float(right[i]),
                cumulative_left_bound_loss_amount=float(cumulative_left[i]),
                cumulative_right_bound_loss_amount=float(cumulative_right[i]),
                cumulative_both_ends_bound_loss_amount=float(cumulative_left[i]+cumulative_right[i])))
        for sample, stop in [(f'{int(t)}_s', t) for t in TIMES] + [('cue_end', trace['duration'])]:
            weight = np.maximum(0., np.minimum(end, stop)-start)
            assert abs(weight.sum()-stop) < 1e-9
            lost_left, lost_right = float(left@weight), float(right@weight)
            summary.append(dict(length_um=length, sample=sample, elapsed_time_s=float(stop),
                cue_duration_s=trace['duration'], left_bound_loss_amount=lost_left,
                right_bound_loss_amount=lost_right, both_ends_bound_loss_amount=lost_left+lost_right))
        audits.append(dict(length_um=length, native_flux_csv_sha256_verified=True,
            reconstructed_endpoint_flux_max_error_amount_s=error,
            nonnegative_bound_outflux_verified=True, cumulative_loss_monotonic=True,
            native_time_coverage_verified=True))
    write_csv('outward_advected_bound_MCA_summary.csv', summary)
    write_csv('outward_advected_bound_MCA_native_intervals.csv', intervals)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.7))
    fig.subplots_adjust(left=.065, right=.99, bottom=.20, top=.66, wspace=.26)
    fig.text(.065, .95, 'DRAFT · Cumulative outward advected bound MCA versus cell length', fontsize=16, weight='bold')
    fig.text(.065, .885, 'Time integral of native endpoint outflux from cue onset · relative-20% right cue', fontsize=11)
    colors = plt.get_cmap('viridis')(np.linspace(.12, .85, len(TIMES)))
    for ax, route, title in zip(axes, ['left', 'right', 'both_ends'],
                                ['Lost through left end', 'Lost through right end', 'Total lost through both ends']):
        for t, color in zip(TIMES, colors):
            selected = [r for r in summary if r['sample'] == f'{int(t)}_s']
            ax.plot([r['length_um'] for r in selected], [r[f'{route}_bound_loss_amount'] for r in selected],
                    color=color, marker='o', ms=4, lw=1.8, label=f'{int(t)} s')
        selected = [r for r in summary if r['sample'] == 'cue_end']
        ax.plot([r['length_um'] for r in selected], [r[f'{route}_bound_loss_amount'] for r in selected],
                color='#202020', ls='--', marker='s', ms=4, lw=1.8, label='Own cue end')
        ax.set_title(title, fontsize=12, pad=65)
        ax.legend(title='Time from cue onset', loc='lower left', bbox_to_anchor=(0, 1.015),
                  frameon=False, ncol=3, fontsize=9, title_fontsize=10,
                  columnspacing=.9, handlelength=1.5, borderaxespad=0)
        ax.set_xlabel('Nominal cell length (µm)')
        ax.set_ylabel('Lost bound MCA (native concentration × µm)', fontsize=10)
        ax.set_xticks([18, 23, 27, 31, 35, 39])
        ax.set_xlim(17, 40)
        ax.set_ylim(bottom=0)
        common.clean_axes(ax)
    fig.text(.065, .088, 'Dashed curves use each cell’s own cue end (226–809 s). Each panel has its own vertical scale.', fontsize=10)
    fig.text(.065, .043, 'Cumulative loss measures outward throughput; uniform replacement can replenish the lost MCA.', fontsize=10)
    return save(fig, 'outward_advected_bound_MCA_vs_length'), audits


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'svg.fonttype': 'none'})
    protected = {name: common.sha(REPORT/name) for name in ['index.html', 'length_mechanisms.pdf']}
    source(Path(__file__))
    source(Path(common.__file__))
    traces, rows, summary, audits = load()
    write_csv('endpoint_bound_MCA_native_time_series.csv', rows)
    write_csv('endpoint_bound_MCA_change_summary.csv', summary)
    stems = [snapshots(summary), means(summary), evolution(traces)]
    loss_stem, loss_audits = cumulative_outward_loss(traces)
    stems.append(loss_stem)
    generated = datetime.now(timezone.utc).isoformat()
    captions = [
        'Observed bound-MCA concentration change at the first and last native finite-volume '
        'cell centres relative to their exact pre-cue values, for all twelve eligible lengths. '
        'Colored curves use exact native samples at common elapsed times. Dashed black curves '
        'use each length’s own selected cue duration, so they do not compare equal times. '
        'These are local endpoint changes, distinct from endpoint flux, net binding and whole-cell amount.',
        'The time-weighted mean observed endpoint net rate equals its concentration change '
        'divided by window duration. All physical and native numerical contributions are included. '
        'The whole-cue and common first-180-s windows match the uniform-replacement draft. '
        'Both endpoints share the same vertical limits; positive means accumulation.',
        'Full cue evolution of the same cumulative endpoint changes, with a common diverging '
        'color scale. Each stripe represents one simulated length; its spatial width is only '
        'a plotting bin. Native time series are linearly interpolated solely onto 1-s display '
        'centres; exact native values are retained in the CSV. Gray begins after the selected '
        'cue duration and contains no recovery data. Black squares mark actual cue withdrawal.'
        ,
        'Cumulative bound-MCA amount advected out through the left and right endpoint faces, '
        'and their sum, from cue onset to each indicated time. The native flux on each implicit '
        'step is multiplied by that step’s duration and summed; no trapezoidal approximation '
        'of stored endpoint concentrations or interior velocities is used. Endpoint donors '
        'are the pre-uniform-correction bound concentrations, with the cortical endpoint '
        'commands used in that step. Bound flux is outward or zero at both ends. Physical '
        'boundary diffusion is zero. All common-time samples precede withdrawal; the dashed '
        'curve uses each length’s own selected cue duration, including its shortened last step. '
        'The plotted quantity is amount in the native 1D concentration-times-length measure, '
        'not a concentration decrease at an endpoint. It includes bound MCA only; total here '
        'means the sum of losses through both ends. Pool replacement can replenish material '
        'that subsequently exits again, so cumulative throughput can exceed the initial bound '
        'inventory. Each panel uses its own vertical limits. The CSV exports exact native-step '
        'fluxes and cumulative losses for all twelve lengths.'
    ]
    sections = []
    for stem, caption in zip(stems, captions):
        output(stem+'_caption.txt').write_text(caption+'\n')
        sections.append(f'<section><h2>{stem.replace("_", " ")}</h2><p>'
                        f'<a href="supporting/vectors/{stem}.pdf">PDF</a> · '
                        f'<a href="supporting/vectors/{stem}.svg">SVG</a></p>'
                        f'<img src="{stem}.png" alt="{stem.replace("_", " ")}"><p>{caption}</p></section>')
    (OUT/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8">'
        '<title>Draft endpoint MCA changes</title><style>body{font:16px system-ui;max-width:1400px;'
        'margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.5}section{margin:40px 0}'
        '</style></head><body><h1>Draft endpoint bound-MCA changes</h1><p>Generated '+generated+
        ' · producer finished · corrected independent endpoints · right relative-20% cue.</p>'
        '<p>Endpoints are sampled at x=L/402 and x=L−L/402, the native cells nearest the boundaries.</p>'
        '<p><a href="supporting/data/endpoint_bound_MCA_native_time_series.csv">Native time series</a> · '
        '<a href="supporting/data/endpoint_bound_MCA_change_summary.csv">Changes and mean rates</a> · '
        '<a href="supporting/data/outward_advected_bound_MCA_summary.csv">Cumulative outward loss</a></p>'+
        ''.join(sections)+'</body></html>')
    assert protected == {name: common.sha(REPORT/name) for name in protected}
    output('endpoint_bound_MCA_audit.json').write_text(json.dumps(dict(generated_at=generated, passed=True,
        producer='finished', endpoint_checks=audits, outward_loss_checks=loss_audits,
        protected_report_sha256=protected,
        no_new_integration=True, source_sha256=SOURCES), indent=2)+'\n')
    artifacts = [p for p in OUT.rglob('*') if p.is_file() and not p.name.endswith('.manifest.json')]
    output('endpoint_bound_MCA.manifest.json').write_text(json.dumps(dict(generated_at=generated,
        passed=True, producer='finished', artifact_root='../..', source_sha256=SOURCES,
        sha256={str(p.relative_to(OUT)):common.sha(p) for p in artifacts},
        reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '
        '/home/4473331/.conda/envs/scaling_sandbox/bin/python -B '
        'polarity1d/manuscript/code/mechanisms/draft_mca_endpoint_changes_20261014.py'), indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT), figures=stems, native_samples=len(rows)), indent=2))


if __name__ == '__main__':
    main()
