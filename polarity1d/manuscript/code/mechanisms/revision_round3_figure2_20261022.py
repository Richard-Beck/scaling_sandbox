"""Figure 2 round 3: redraw immutable native evidence; never integrate or submit.

Run with the existing scaling_sandbox Python. Main, all-area width supplement,
endpoint audit, source data and exact preservation receipts are written together.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from pypdf import PdfReader, PdfWriter

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'polarity1d/manuscript').is_dir())
DEV = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
OUT = DEV / 'revision_round3_20261022/figure2'
OLD = DEV / 'signature_assays_20261020'
PREVIOUS = DEV / 'revision_round2_20261022/figure2'
RUNS = ROOT / 'polarity1d/manuscript/data/runs'
ORDER = ('goryachev', 'wave_pinning', 'holmes_model3', 'debelly', 'spring')
NAMES = dict(goryachev='Goryachev', wave_pinning='Wave-pinning', holmes_model3='Holmes', debelly='deBelly', spring='Spring')
COLORS = dict(goryachev='#c83b3b', wave_pinning='#2878b5', holmes_model3='#36944c', debelly='#82599c', spring='#c59a25')
LENGTHS = (19, 23, 27, 31, 35, 39)
WIDTHS = (.95, 1.9, 3.8, 7.6, 15.2)
AREAS = (.5, 1., 2.)
STYLES = ('-', (0, (5, 2)), (0, (3, 1, 1, 1)), (0, (1, 1.3)), (0, (2, 2)))
SOURCES = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return json.loads(path.read_text())


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def csvout(name, rows):
    path = OUT / 'supporting/data' / name
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def clean(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=2, width=.55, pad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(.55)


def native_suspect_audit(record):
    """Recompute pulse1/gap1/pulse2 from three native phase arrays, not D.

    This checks actual amounts, current-profile minimum, fractional half-cell
    integration and visibility. No solver or new numerical parameter is used.
    """
    source = record['source']
    analysis = read(record['audited_metadata']['path'])
    assert sha(record['audited_metadata']['path']) == source['worker_analysis_receipt']['sha256']
    result_path = Path(source['run_root']) / 'workers' / f"{source['task']:04d}" / 'result.json'
    assert sha(result_path) == source['worker_analysis_receipt']['result_sha256']
    result = read(result_path)
    folder = result_path.parent / 'trial'
    summary = []
    for index, endpoint in enumerate(record['endpoints'][:3]):
        path = folder / f'phase_{index:03d}.npz'
        expected = result['artifact_sha256'][str(path.relative_to(Path(source['run_root'])))]
        assert sha(path) == expected
        with np.load(path) as native:
            active = native['active'][-1]
            x = native['x']
        dx = float(x[1] - x[0])
        edges = np.r_[x[0] - dx / 2, x + dx / 2]
        midpoint = (edges[0] + edges[-1]) / 2
        right_weights = np.clip(edges[1:] - np.maximum(edges[:-1], midpoint), 0, dx)
        excess = active - np.min(active)
        total = float(np.sum(excess * np.diff(edges)))
        challenger = float(np.sum(excess * right_weights))
        raw_share = challenger / total
        contrast = float(np.ptp(active) / (np.max(active) + np.min(active)))
        saved = endpoint['end']
        for value, key in ((total, 'E_total'), (challenger, 'E_challenger_half'),
                           (raw_share, 'share_half_raw'), (contrast, 'contrast')):
            assert np.isclose(value, saved[key], rtol=1e-10, atol=1e-12), (key, value, saved[key])
        assert (contrast >= saved['contrast_threshold']) == saved['contrast_eligible']
        summary.append(dict(kind=endpoint['kind'], time=saved['time'],
            raw_share=raw_share, visible_share=saved['share_half'],
            contrast=contrast, contrast_threshold=saved['contrast_threshold'],
            challenger_amount=challenger, total_excess_amount=total,
            readout_eligible=saved['readout_eligible'], native_phase=str(path.relative_to(ROOT)),
            native_phase_sha256=expected))
    return dict(passed=True, source_task=source['task'], endpoint_checks=summary,
        finding='First gap is below the existing activity-visibility gate; no fully retained visible activity.',
        raw_gap_retention=summary[1]['raw_share'] / summary[0]['raw_share'],
        visible_gap_retention=0., pulse2_headroom_gain=(summary[2]['raw_share'] - summary[0]['raw_share']) / (1 - summary[0]['raw_share']),
        missing_value_policy='Only an observed endpoint explicitly below the existing activity gate is placed at x=0. Missing samples are never zero or forward-filled.')


def extract():
    coverage = read(OLD / 'supporting/metadata/coverage.json')
    assert coverage['final_export_eligible'] and not coverage['audit_errors']
    orgfile = Path(coverage['retained_organization']['path'])
    assert sha(orgfile) == coverage['retained_organization']['sha256']
    SOURCES[str(orgfile.relative_to(ROOT))] = sha(orgfile)
    organization = list(csv.DictReader(orgfile.open()))
    assert len(organization) == 97
    csvout('A_incumbent_organization.csv', organization)
    initiation = list(csv.DictReader((PREVIOUS / 'initiation_normalized.csv').open()))
    SOURCES[str((PREVIOUS / 'initiation_normalized.csv').relative_to(ROOT))] = sha(PREVIOUS / 'initiation_normalized.csv')
    definition = read(OLD / 'supporting/metadata/reporter_initiation_definition.json')
    csvout('B_initiation.csv', initiation)
    widths = []
    for name in ('signature_width_core_20261020', 'signature_width_holmes_spring_20261020'):
        folder = RUNS / name
        status = read(folder / 'results/analysis_20261020/status.json')
        assert status['declared_subset_complete'] and not status.get('errors') and not status['missing_tasks']
        assert sha(folder / 'manifest.json') == status['manifest_sha256']
        assert sha(folder / 'results/analysis_20261020/event_rows.json') == status['outputs']['event_rows.json']
        events = read(folder / 'results/analysis_20261020/event_rows.json')
        for row in events:
            assert row['label'] == 'trial'
            event = row['continuous_cue_on_events']['D0_5']
            assert event['status'] in ('observed', 'right_censored')
            widths.append(dict(model=row['model'], length=row['requested_reference_length'],
                width=row['requested_width_label'], area=row['integrated_area_factor'],
                time_lower=event['lower_elapsed'], time_upper=event['upper_elapsed'],
                time_estimate=event['estimated_elapsed'], lower_bound=event['right_censored_at'],
                status=event['status'], field_tail_pass=row['field_tail_pass'],
                run=name, task=row['task']))
    assert len(widths) == 450
    assert {(r['model'], r['length'], r['width'], r['area']) for r in widths} == {
        (m, length, width, area) for m in ORDER for length in LENGTHS for width in WIDTHS for area in AREAS}
    csvout('C_replacement_times_all450.csv', widths)
    holmes = [r for r in widths if r['model'] == 'holmes_model3' and r['width'] == 15.2]
    assert len(holmes) == 18 and all(r['status'] == 'observed' for r in holmes)
    write(OUT / 'supporting/metadata/Holmes_width_decision.json', dict(
        additional_runs_needed=False, widest_existing_width=15.2,
        finite_conditions=18, tested_lengths=list(LENGTHS), tested_areas=list(AREAS),
        reason='Existing width 15.2 yields a finite replacement time at every length and all three areas. Narrow-cue censoring does not prevent the per-width comparison.',
        no_kinetic_retuning=True))
    all_trains = read(OLD / 'supporting/data/figure2_S2_all_native_pulse_gap_endpoints.json')
    display_audit = read(OLD / 'supporting/metadata/S1_S2_display_audit.json')
    rel = 'supporting/data/figure2_S2_all_native_pulse_gap_endpoints.json'
    assert sha(OLD / rel) == display_audit['artifact_sha256'][rel]
    trains = [r for r in all_trains if r['case']['pulse_count'] == 4 and r['case']['total_on_factor_of_calibration'] == 1.5]
    additional = read(RUNS / 'signature_holmes_retention_round2_20261022/results/independent_native_readout_audit.json')
    assert additional['passed'] and additional['workers'] == 32
    for row in additional['trajectories']:
        case = row['case']
        if case['kind'] == 'train':
            assert row['passed'] and sha(row['result_path']) == row['result_sha256']
            trains.append(dict(case=dict(model='holmes_model3', requested_reference_length=case['length'],
                gap_multiple_of_pulse=case['gap']), endpoints=row['endpoints'],
                field_tail_pass=row['trial']['field_tail_pass'], source=dict(run_root=row['campaign_root'], task=row['task']),
                schedule=row['trial']['schedule']))
    assert len(trains) == 180
    suspect = next(r for r in trains if r['case']['model'] == 'debelly' and r['case']['requested_reference_length'] == 19 and r['case']['gap_multiple_of_pulse'] == 16)
    suspect_proof = native_suspect_audit(suspect)
    write(OUT / 'supporting/metadata/debelly_L19_gap16_native_audit.json', suspect_proof)
    endpoints = []
    for row in trains:
        case = row['case']
        on = [e for e in row['endpoints'] if e['kind'] == 'on']
        gap = next(e for e in row['endpoints'] if e['kind'] == 'gap')['end']
        first, second = on[0]['end'], on[1]['end']
        assert first['readout_eligible'] and second['readout_eligible']
        s1, s2 = first['share_half'], second['share_half']
        assert 0 < s1 < 1 and 0 <= s2 <= 1
        assert gap['status'] == 'observed_exact_endpoint'
        below_gate = gap['readout_eligible'] is False
        if below_gate:
            assert gap['share_half'] is None
            visible_gap_share = 0.
        else:
            assert gap['share_half'] is not None
            visible_gap_share = gap['share_half']
        endpoints.append(dict(model=case['model'], length=case['requested_reference_length'],
            gap=case['gap_multiple_of_pulse'], pulse1_share=s1, gap1_visible_share=visible_gap_share,
            gap1_raw_share=gap['share_half_raw'], pulse2_share=s2,
            activity_kept=visible_gap_share / s1,
            raw_share_retention=gap['share_half_raw'] / s1,
            headroom_gain=(s2 - s1) / (1 - s1),
            below_activity_gate=below_gate, gap_contrast=gap['contrast'],
            contrast_threshold=gap['contrast_threshold'], field_tail_pass=row['field_tail_pass'],
            target_length_calibrated=not(case['model'] == 'holmes_model3' and case['requested_reference_length'] >= 31),
            run=str(Path(row['source']['run_root']).relative_to(ROOT)), task=row['source']['task']))
    csvout('D_gap1_vs_pulse2_all180.csv', endpoints)
    examples = read(OLD / 'supporting/data/selected_train_endpoints.json')
    h39 = next(r for r in trains if r['case']['model'] == 'holmes_model3' and r['case']['requested_reference_length'] == 39 and r['case']['gap_multiple_of_pulse'] == 4)
    examples = [r for r in examples if r['model'] != 'holmes_model3'] + [dict(model='holmes_model3', **h39)]
    traces = []
    for row in examples:
        model = row['model']
        folder = Path(row['source']['run_root']) / 'workers' / f"{row['source']['task']:04d}" / 'trial'
        measurement = read(folder / 'readouts.json')
        result = read(folder.parent / 'result.json')
        assert sha(folder / 'readouts.json') == result['trial']['artifacts']['readouts']['sha256']
        series = measurement['series']
        pulse = row['schedule']['pulse_duration']
        stop = row['schedule']['segments'][-1][2]
        time = np.array(series['time'])
        # Positive-excess share remains mathematically defined below visibility;
        # keep raw shares in E and explicitly distinguish D's visibility gate.
        share = np.asarray(series['share_half_raw'], dtype=float)
        keep = np.ones(len(time), dtype=bool) if model == 'debelly' else time <= stop + pulse
        from signature_figure2_supplement_20261020 import display_indices
        selected = display_indices(time[keep], share[keep], endpoints=[e['end']['time'] for e in row['endpoints'] if e['kind'] in ('on', 'gap')], bins=300)
        indices = np.flatnonzero(keep)[selected]
        traces.extend(dict(model=model, time=float(time[i]), pulse_time=float(time[i] / pulse),
            challenger_share=float(share[i]), readout_eligible=series['D_half_eligible'][i],
            final_share=float(share[-1]), final_time=float(time[-1]), pulse_duration=pulse)
            for i in indices)
    csvout('E_L39_challenger_share_traces.csv', traces)
    return organization, initiation, widths, endpoints, traces, definition, suspect_proof


def width_axes(ax, widths, model, area, row_index, *, height_label=False):
    for width, style in zip(WIDTHS, STYLES):
        rows = sorted([r for r in widths if r['model'] == model and r['area'] == area and r['width'] == width], key=lambda r: r['length'])
        # Separate finite pieces; censored horizons are never joined as estimates.
        values = [r['time_estimate'] if r['status'] == 'observed' else np.nan for r in rows]
        ax.plot(LENGTHS, values, color=COLORS[model], ls=style, lw=.95)
        for r in rows:
            if r['status'] == 'right_censored':
                ax.plot(r['length'], r['lower_bound'], marker='^', mfc='white', mec=COLORS[model], ms=3, mew=.65, zorder=4)
                ax.annotate('', xy=(r['length'], r['lower_bound'] * 1.65), xytext=(r['length'], r['lower_bound'] * 1.02),
                    arrowprops=dict(arrowstyle='->', lw=.6, color=COLORS[model], shrinkA=0, shrinkB=0, mutation_scale=5))
    all_model = [r for r in widths if r['model'] == model]
    finite = [r['time_estimate'] for r in all_model if r['status'] == 'observed']
    highest = max([r['lower_bound'] or r['time_upper'] for r in all_model])
    ax.set(xlim=(17.7, 40.3), yscale='log', ylim=(min(finite) / 1.4, highest * 2.2), xticks=[19, 29, 39])
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=3))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.tick_params(labelbottom=row_index == 2)
    if row_index == 2:
        ax.set_xlabel('Length L', labelpad=1.5)
    clean(ax)


def schematic(ax, kind, *, width=False):
    """Deliberately dimensionless predictions, never labelled simulated data."""
    if width:
        x = np.linspace(0, 1, 60)
        if kind == 0:
            small = .32 - .08 * x
            large = .85 - .13 * x
        elif kind == 1:
            small = .12 + .78 * np.exp(-4 * x)
            large = .18 + 1.4 * np.exp(-3.8 * x)
        else:
            small = .36 - .04 * x
            large = .80 - .06 * x
        ax.plot(x, small, color='.42', lw=1.)
        ax.plot(x, large, color='.08', lw=1., ls=(0, (3, 1.6)))
        ax.set(xlim=(0, 1), ylim=(0, 1.2), xticks=[0, 1], xticklabels=['narrow', 'broad'], yticks=[])
        ax.set_xlabel('Width at fixed dose', labelpad=0)
    else:
        times = np.linspace(0, 8, 401)
        def pulse_line(large):
            values = np.full_like(times, .04)
            for j in range(4):
                on = (times >= 2*j) & (times < 2*j+1)
                off = (times >= 2*j+1) & (times < 2*j+2)
                phase = times[on] - 2*j
                if kind == 0:
                    delay = .5 if large else .15
                    values[on] = .04 + .64 * np.maximum(phase-delay, 0) / (1-delay)
                    values[off] = .04 + .64 * np.exp(-7*(times[off]-2*j-1))
                elif kind == 1:
                    base = .04 + j * (.13 if large else .2)
                    peak = base + (.13 if large else .2)
                    values[on] = base + (peak-base) * (1-np.exp(-7*phase))
                    values[off] = peak - .025*(times[off]-2*j-1)
                else:
                    peak = min(.96, .36 + j*(.12 if large else .24))
                    values[on] = .04 + (peak-.04) * (1-np.exp(-8*phase))
                    commits = j >= (3 if large else 2)
                    values[off] = peak if commits else .04 + (peak-.04)*np.exp(-8*(times[off]-2*j-1))
            return values
        for j in range(4):
            ax.axvspan(2*j, 2*j+1, color='.93', lw=0)
        ax.plot(times, pulse_line(False), color='.42', lw=1.)
        ax.plot(times, pulse_line(True), color='.08', lw=1., ls=(0, (3, 1.6)))
        ax.set(xlim=(0, 8), ylim=(0, 1.04), xticks=[], yticks=[])
    clean(ax)


def render(data):
    org, initiation, widths, endpoints, traces, _, _ = data
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.2,
        'axes.labelsize': 7.2, 'axes.titlesize': 7.5, 'xtick.labelsize': 6.8,
        'ytick.labelsize': 6.8, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'axes.titlecolor': 'black', 'text.color': 'black', 'axes.labelcolor': 'black'})
    fig = plt.figure(figsize=(7.2, 10.0))
    label_artists = []
    def label(letter, height):
        label_artists.append(fig.text(.035, height / 10, letter, fontsize=10,
            fontweight='bold', color='black', ha='left', va='top'))
    def facets(bottom, height, count=5, left=.63, right=7.08, gap=.24, share_y=None):
        width = (right-left-gap*(count-1))/count
        return [fig.add_axes([(left+j*(width+gap))/7.2, bottom/10, width/7.2, height/10],
            sharey=share_y if share_y is not None else None) for j in range(count)]
    def model_labels(axes):
        for ax, model in zip(axes, ORDER):
            ax.set_title(NAMES[model], color='black', pad=4)
    fig.legend([Line2D([], [], color=COLORS[m], lw=1.5) for m in ORDER],
        [NAMES[m] for m in ORDER], loc='upper center', bbox_to_anchor=(.55, .997),
        ncol=5, frameon=False, handlelength=1.5, columnspacing=1.2, fontsize=7.2)
    # A: exactly the preceding 97 states, quantities and normalization.
    label('A', 9.72)
    axes = facets(8.86, .72, count=3, gap=.47)
    for ax, metric, name in zip(axes, ('excess_amount', 'width', 'peak'), ('Amount', 'Extent', 'Peak')):
        for model in ORDER:
            rows = sorted([r for r in org if r['model'] == model], key=lambda r: float(r['size_native']))
            ax.plot([float(r['size_over_reference']) for r in rows], [float(r[metric+'_over_reference']) for r in rows], color=COLORS[model], lw=1.)
        ax.set(xlabel='Length / reference', ylabel=name+' / reference')
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=3, prune='both'))
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 2
        clean(ax)
    # B: shared normalized scale for four models; all six native bounds for deBelly.
    label('B', 8.50)
    axes = facets(7.61, .69)
    model_labels(axes)
    for ax, model in zip(axes, ORDER):
        rows = [r for r in initiation if r['model'] == model]
        if model == 'debelly':
            assert len(rows) == 6 and all(float(r['time_upper']) == .5 for r in rows)
            ax.plot(LENGTHS, [.5]*6, 'v', ms=3.5, mfc='white', mec=COLORS[model], mew=.75)
            for length in LENGTHS:
                ax.annotate('', xy=(length, .05), xytext=(length, .46), arrowprops=dict(arrowstyle='->', lw=.65, color=COLORS[model], mutation_scale=5))
            ax.set(ylim=(0, .6), yticks=[0, .5])
            ax.text(.5, .04, 'Native time', transform=ax.transAxes, ha='center', fontsize=6.5)
        else:
            ax.plot([float(r['length']) for r in rows], [float(r['ratio_estimate']) for r in rows], color=COLORS[model], lw=1.)
            for r in rows:
                value = float(r['ratio_estimate'])
                ax.errorbar(float(r['length']), value, yerr=[[value-float(r['ratio_lower'])], [float(r['ratio_upper'])-value]], fmt='none', color=COLORS[model], lw=.65, capsize=1)
            ax.set(yscale='log', ylim=(.75, 4.1), yticks=[1, 2, 4], yticklabels=['1', '2', '4'])
            ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.axhline(1, color='.7', ls=':', lw=.6)
            if model != 'goryachev':
                ax.tick_params(labelleft=False)
        ax.set(xlim=(17, 41), xticks=[19, 29, 39], xlabel='Length L')
        ax.xaxis.labelpad = 1.5
        clean(ax)
    axes[0].set_ylabel('Initiation / L19', labelpad=2)
    # C: all three areas and all five widths, without time ratios.
    label('C', 7.22)
    fig.legend([Line2D([], [], color='.15', lw=1, ls=s) for s in STYLES] +
        [Line2D([], [], color='.2', lw=0, marker='^', ms=3, mfc='white')],
        [f'w = {w:g}' for w in WIDTHS] + ['> horizon'], loc='upper center', bbox_to_anchor=(.56, .730),
        ncol=6, frameon=False, handlelength=1.8, columnspacing=.8, fontsize=6.8)
    for i, area in enumerate(AREAS):
        axes = facets(6.34 - i*.75, .55)
        if i == 0:
            model_labels(axes)
        for ax, model in zip(axes, ORDER):
            width_axes(ax, widths, model, area, i)
        axes[0].text(.035, .965, f'Area {area:g}×', transform=axes[0].transAxes,
            fontsize=6.5, ha='left', va='top', color='black')
    fig.text(.013, .590, 'Replacement time (native units)', fontsize=7.2, rotation=90, va='center')
    # D: other models behind each model; visibility-qualified left wall.
    label('D', 4.37)
    fig.text(.13, .437, 'Gap 1 retained share', fontsize=7.2)
    fig.text(.57, .437, 'L19', fontsize=6.7)
    for i, length in enumerate(LENGTHS):
        color = str(.72 - .12*i)
        fig.add_artist(Line2D([.615+i*.017], [.440], marker='o', ms=3.2, color=color, transform=fig.transFigure))
    fig.text(.724, .437, 'L39', fontsize=6.7)
    fig.text(.13, .422, 'Below activity gate: open circle at x = 0', fontsize=6.4)
    axes = facets(3.17, .82)
    model_labels(axes)
    for j, (ax, model) in enumerate(zip(axes, ORDER)):
        other = [r for r in endpoints if r['model'] != model]
        ax.scatter([r['activity_kept'] for r in other], [r['headroom_gain'] for r in other], s=6, c='.82', linewidths=0, zorder=1)
        for r in [r for r in endpoints if r['model'] == model]:
            light = .58*(39-r['length'])/20
            color = tuple((1-light)*c + light for c in to_rgb(COLORS[model]))
            ax.scatter(r['activity_kept'], r['headroom_gain'], s=11,
                facecolors='white' if r['below_activity_gate'] else [color],
                edgecolors=[color], linewidths=.75 if r['below_activity_gate'] else .15, zorder=3)
        ax.set(xlim=(-.055, 1.06), ylim=(-.05, 1.02), xticks=[0, .5, 1], yticks=[0, .5, 1])
        ax.axhline(0, lw=.5, color='.65', zorder=0)
        ax.axvline(0, lw=.5, color='.65', zorder=0)
        ax.set_xlabel('Gap-end / pulse 1', labelpad=1.5)
        if j:
            ax.tick_params(labelleft=False)
        clean(ax)
        if model == 'debelly':
            ax.annotate('L19, gap 16', xy=(0, .8594472898053895), xytext=(.14, .91), fontsize=5.9,
                arrowprops=dict(arrowstyle='-', lw=.5, color='.3'), ha='left')
    axes[0].set_ylabel('Pulse-2 gain / headroom', labelpad=2)
    # E: same selected native histories, on a common challenger-share scale.
    label('E', 2.78)
    axes = facets(1.62, .86)
    model_labels(axes)
    classes = dict(goryachev='lost', wave_pinning='kept', holmes_model3='kept', debelly='lost then commits', spring='lost then commits')
    for j, (ax, model) in enumerate(zip(axes, ORDER)):
        rows = [r for r in traces if r['model'] == model]
        ax.plot([r['pulse_time'] for r in rows], [r['challenger_share'] for r in rows], color=COLORS[model], lw=1.)
        for pulse in range(4):
            ax.axvspan(pulse*5, pulse*5+1, color=COLORS[model], alpha=.12, lw=0)
        ax.set(ylim=(-.02, 1.04), yticks=[0, .5, 1], xlabel='Time / pulse')
        xmax = max(r['pulse_time'] for r in rows)
        ax.set_xlim(0, xmax)
        ax.set_xticks([0, 16, 32] if model == 'debelly' else [0, 8, 16])
        ax.xaxis.labelpad = 1.5
        ax.set_title(NAMES[model], color='black', pad=15)
        ax.text(.5, 1.025, classes[model], transform=ax.transAxes, fontsize=6.4, va='bottom', ha='center', color='black')
        ax.text(.99, .025, 'L39', transform=ax.transAxes, fontsize=6.2, ha='right')
        if j:
            ax.tick_params(labelleft=False)
        clean(ax)
    axes[0].set_ylabel('Challenger share', labelpad=2)
    # F: three distinct prospective biosensor signatures under the same input.
    label('F', 1.25)
    fig.text(.13, .125, 'Predictions (schematic)', fontsize=7.2)
    fig.legend([Line2D([], [], color='.42', lw=1), Line2D([], [], color='.08', lw=1, ls=(0, (3, 1.6)))],
        ['WGD−', 'WGD+'], loc='upper right', bbox_to_anchor=(.99, .140), ncol=2,
        frameon=False, fontsize=6.7, handlelength=2, columnspacing=1.)
    upper = facets(.82, .27, count=3, gap=.38, right=6.92)
    lower = facets(.35, .27, count=3, gap=.38, right=6.92)
    for j, (ax, bx, name) in enumerate(zip(upper, lower,
            ('Initiation-limited', 'Saturating activation', 'Slow non-activity state'))):
        schematic(ax, j)
        schematic(bx, j, width=True)
        ax.set_title(name, fontsize=7.2, pad=3, color='black')
    upper[0].set_ylabel('Biosensor', labelpad=2, fontsize=6.7)
    lower[0].set_ylabel('Replace\ntime', labelpad=2, fontsize=6.7)
    # Same metrics for every letter, fixed position and black text, verified below.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    letter_checks = [dict(letter=a.get_text(), color=a.get_color(), size=a.get_fontsize(),
        weight=a.get_fontweight(), x=a.get_position()[0]) for a in label_artists]
    assert len(letter_checks) == 6 and len({(r['color'], r['size'], r['weight'], r['x']) for r in letter_checks}) == 1
    # Detect figure-edge clipping, allowing no text outside the printed canvas.
    clipped = []
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            if box.x0 < -.002 or box.x1 > 1.002 or box.y0 < -.002 or box.y1 > 1.002:
                clipped.append(dict(text=artist.get_text(), bounds=list(box.bounds)))
    write(OUT / 'supporting/metadata/display_audit.json', dict(
        page_inches=[7.2, 10], page_mm=[182.88, 254], labels=letter_checks,
        clipped_text=clipped, figure_title=None, facet_label_color='black'))
    assert not clipped, clipped
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(OUT / f'figure2_round3.{suffix}', dpi=300, facecolor='white')
    plt.close(fig)


def supplements(data):
    from matplotlib.backends.backend_pdf import PdfPages
    _, _, widths, endpoints, _, _, suspect = data
    added = OUT / 'supporting/vectors/round3_endpoint_and_width_audit.pdf'
    with PdfPages(added) as book:
        for area in AREAS:
            fig, axes = plt.subplots(1, 5, figsize=(10, 3))
            for ax, model in zip(axes, ORDER):
                width_axes(ax, widths, model, area, 2)
                ax.set_title(NAMES[model], color='black')
            axes[0].set_ylabel('Replacement time (native units)')
            fig.suptitle(f'Equal-area width series: area {area:g}×; all finite times and lower bounds', fontsize=11)
            fig.legend([Line2D([], [], color='.2', ls=s) for s in STYLES], [f'w = {w:g}' for w in WIDTHS], loc='lower center', ncol=5, frameon=False)
            fig.subplots_adjust(bottom=.24, top=.77, left=.07, right=.98, wspace=.4)
            book.savefig(fig)
            plt.close(fig)
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(.08, .94, 'deBelly L19, gap16: exact native endpoint audit', fontsize=14, weight='bold')
        lines = [
            'The first gap does not fully retain visible activity.',
            f"Pulse-1 end share: {suspect['endpoint_checks'][0]['raw_share']:.6f}.",
            f"Gap-1 raw share: {suspect['endpoint_checks'][1]['raw_share']:.6f}; raw ratio: {suspect['raw_gap_retention']:.6f}.",
            f"Gap contrast: {suspect['endpoint_checks'][1]['contrast']:.6f}, below the existing 0.05 gate.",
            'Eligible gap share is unavailable because activity is below that visibility gate.',
            'The main D places this known below-gate observation at x = 0 with an open circle.',
            'This is a visibility code, not an assertion of mathematically zero raw activity.',
            f"Pulse-2 end share: {suspect['endpoint_checks'][2]['raw_share']:.6f}.",
            f"Pulse-2 headroom gain: {suspect['pulse2_headroom_gain']:.6f}.",
            'Independent checks recompute all three endpoints from phase-native active fields.',
            'The native archives and audited endpoint table agree; no endpoint is forward-filled.',
            'D uses the first gap only. Later gaps must not substitute for it.',
            'All 180 D coordinates, raw shares, visibility and field-tail flags are in the CSV.',
            'Holmes L31–39 uses the existing L27-referenced train, with no target-length calibration.',
            'deBelly timing and endpoint-field qualifications remain in the retained supplement.',
        ]
        for i, line in enumerate(lines):
            fig.text(.08, .87-i*.043, line, fontsize=9)
        book.savefig(fig)
        plt.close(fig)
    parts = [added, PREVIOUS / 'figure2_round2.pdf', PREVIOUS / 'figure2_round2_supplement.pdf']
    writer = PdfWriter()
    for part in parts:
        writer.append(part)
    writer.write(OUT / 'figure2_round3_supplement.pdf')
    current = PdfReader(OUT / 'figure2_round3_supplement.pdf')
    offset = len(PdfReader(added).pages)
    retained = []
    for part in parts[1:]:
        original = PdfReader(part)
        for page in original.pages:
            assert page.get_contents().get_data() == current.pages[offset].get_contents().get_data()
            offset += 1
        retained.append(dict(path=str(part.relative_to(ROOT)), sha256=sha(part), pages=len(original.pages), decoded_page_identity=True))
    write(OUT / 'supporting/metadata/supplement_preservation.json', dict(
        passed=True, supplement_pages=len(current.pages), added_pages=len(PdfReader(added).pages),
        preserved=retained, superseded_panels='Complete round-2 main preserved, including width ratios, fate maps and class table; all previous supplements follow unchanged.'))


def content(data):
    _, _, widths, endpoints, _, _, suspect = data
    generated = datetime.now(timezone.utc).isoformat()
    legend = [
        'Enlarged cells resist replacement through distinct initiation, cue-width and pulse-memory signatures.',
        'A, Incumbent excess amount, extent and peak, normalized to each model’s reference, across 97 retained states.',
        'B, Initiation is the time from cue onset until maximum in-cue activity minus the profile minimum exceeds 10% of the model’s L19 incumbent excess peak. Area 1, width 3.8; times are normalized to L19. All six deBelly upper bounds are 0.5 native units (downward arrows); first-update sampling leaves its size dependence unresolved.',
        'C, Replacement time (first eligible challenger share ≥ 0.75) versus length, with a line per width and rows for areas 0.5, 1 and 2. Upward triangles/arrows are lower bounds: the cue never replaced within the horizon (6,000 Wave-pinning; 10,000 Goryachev/Holmes; 2,000 deBelly; 100 Spring, native units). Broad cues rescue Holmes across all lengths; narrow-cue censoring does not establish a growing penalty. Spring native widths equal the displayed widths divided by 38.',
        'D, Gap-1 activity kept versus pulse-2 gain as a fraction of headroom (Eq. 1). Four equal pulses, total ON 1.5 times the model-length commitment calibration; gaps 0–16 pulse durations. Holmes L31–39 uses the L27 calibration. Other models are grey; lightness indicates L. Open circles at x = 0 code activity below the existing visibility gate. At the left wall, deBelly/Spring gain on pulse 2 despite no visible activity; Goryachev never does. The deBelly L19/gap16 first-gap endpoint was checked against native fields.',
        'E, L39 gap4 challenger-share histories, labelled kept / lost / lost then commits; shading marks cue ON. deBelly includes full recovery; others show one post-train pulse duration. Continuous raw shares include below-gate values. deBelly’s later timing and fields remain numerically sensitive (supplement).',
        'F, Schematic biosensor predictions for the same pulse train and equal-dose width series in WGD− and enlarged WGD+ cells: delayed rises and resets; saturation with gap retention and strong width rescue; or gap collapse followed by commitment, with weaker width rescue. Predictions are prospective. Model kinetics are unchanged.'
    ]
    write(OUT / 'legend.json', dict(id='figure2', title='Replacement resistance in enlarged cells has distinct size and pulse-memory signatures',
        status='round3_complete_qualified', role='Separate size-dependent initiation and cue-width responses from observable pulse-memory classes.',
        assets=[dict(path=str((OUT / 'figure2_round3.png').relative_to(ROOT)), caption='Figure 2. Revised size and pulse-memory signatures.',
            page_size=[7.2, 10], minimum_printed_width_inches=7.2,
            native_pdf_page=str((OUT / 'figure2_round3.pdf').relative_to(ROOT)),
            page_title='Figure 2. Size and pulse-memory signatures')],
        legend=legend, equations=[dict(math=r'x=\frac{s_{\mathrm{gap}\,1}}{s_1},\qquad y=\frac{s_2-s_1}{1-s_1}', number='1',
            description='s₁ and s₂ are challenger shares at pulse ends 1 and 2; s_gap1 is the first-gap-end share. An observed endpoint below the existing visibility gate is coded as zero visible share in D; its raw value remains in the CSV.')], narrative=[legend[0],
            'The full width series shows that Holmes can replace at every tested length when the cue is broad enough; narrow-cue censoring cannot establish a growing width-penalty ratio.',
            'Pulse-2 headroom gain separates low visible retention with residual progress from the initiation-limited reset response. Saturating activation retains visible activity through gaps.'],
        measurements=['Measure first in-cue biosensor rise, then vary width at fixed integrated dose.',
            'Measure pulse-1, gap-1 and pulse-2 challenger shares under the same pulse train, followed by eventual commitment.'],
        limitations=['Commanded input areas and native clocks are not calibrated across models.',
            'deBelly initiation is unresolved within the first native update; later-gap timing and endpoint fields remain numerically sensitive.',
            'D codes known below-gate activity at x=0; the exact raw ratio remains in the plotted-value CSV.',
            'Holmes large-length trains have no target-length commitment calibration; visible retention differs from eventual replacement.',
            'The three classes describe selected deterministic regimes; WGD biosensor predictions are prospective.'],
        main_pdf=str((OUT / 'figure2_round3.pdf').relative_to(ROOT)), supplement_pdf=str((OUT / 'figure2_round3_supplement.pdf').relative_to(ROOT)),
        reader=str((OUT / 'README.md').relative_to(ROOT)),
        sources=[dict(label=label, path=str((OUT / path).relative_to(ROOT))) for label, path in [
            ('All 450 per-width replacement times', 'supporting/data/C_replacement_times_all450.csv'),
            ('All 180 gap1/pulse2 coordinates and raw shares', 'supporting/data/D_gap1_vs_pulse2_all180.csv'),
            ('deBelly L19/gap16 native-phase audit', 'supporting/metadata/debelly_L19_gap16_native_audit.json'),
            ('Preserved prior main panels and complete supplement', 'supporting/metadata/supplement_preservation.json')]],
        main_source_data=[dict(panel=p, path=str((OUT / 'supporting/data' / n).relative_to(ROOT)), sha256=sha(OUT / 'supporting/data' / n)) for p, n in [
            ('A', 'A_incumbent_organization.csv'), ('B', 'B_initiation.csv'), ('C', 'C_replacement_times_all450.csv'),
            ('D', 'D_gap1_vs_pulse2_all180.csv'), ('E', 'E_L39_challenger_share_traces.csv')]]))
    readme = f'''# Figure 2 — round 3

[Main PDF](figure2_round3.pdf) · [PNG](figure2_round3.png) · [SVG](figure2_round3.svg) · [Complete supplement](figure2_round3_supplement.pdf) · [Legend](legend.json)

Enlargement changes initiation and replacement responses. Pulse memory distinguishes initiation-limited, saturating-activation and slow non-activity-state classes. A–C show size responses; D–E distinguish temporal signatures; F gives prospective biosensor predictions.

The main is 7.2 × 10 inches (183 × 254 mm), with uniform black A–F labels and black facet labels. A uses the exact preceding 97 states. B includes all six deBelly upper bounds. C shows all 450 width-campaign observations, with separate lines per width and rows for areas 0.5, 1 and 2. All 18 Holmes width-15.2 conditions are finite; no extra width runs or kinetic retuning are needed. C replaces the unsupported claim that Holmes’s width penalty grows through L39.

D uses first-gap retention and pulse-2 headroom gain from the existing 180 endpoints, with other models grey and lightness showing length. The deBelly L19/gap16 audit recomputes the actual pulse1/gap1/pulse2 native fields: gap contrast is {suspect['endpoint_checks'][1]['contrast']:.6f}, below the existing 0.05 activity gate. Its raw retention is {suspect['raw_gap_retention']:.6f}, not full retention; its pulse-2 headroom gain is {suspect['pulse2_headroom_gain']:.6f}. The open point at x=0 codes known below-gate visible activity, not a mathematically zero raw share. Missing samples are never converted to zero. Raw and visibility-qualified values are both exported. E uses continuous raw positive-excess challenger share, including below-gate values.

F’s curves are explicitly qualitative assay predictions. The same pulse train and equal-dose width series produce different proposed signatures for each class. WGD-associated enlargement motivates length; simulated polarity, observed migration and conditional OU persistence remain separate evidence.

All previous supplements and the complete superseded round-2 main are appended with decoded PDF page contents verified unchanged. This preserves the fate maps, ratio panel and class tables. Existing deBelly smaller-step/field qualifications, Holmes large-L calibration limits and noisy-gap exceptions remain available there.

Sources and plotted values: [CSV directory](supporting/data/) · [native endpoint audit](supporting/metadata/debelly_L19_gap16_native_audit.json) · [width-run decision](supporting/metadata/Holmes_width_decision.json) · [display audit](supporting/metadata/display_audit.json) · [preservation receipt](supporting/metadata/supplement_preservation.json).

Generated UTC: {generated}. Producer complete. No automatic browser refresh. Redraw using the existing scaling_sandbox environment: `python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure2_20261022.py`.
'''
    (OUT / 'README.md').write_text(readme)
    body = '<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 2 round 3</title><style>body{max-width:1000px;margin:2rem auto;padding:0 1rem;font:16px/1.5 system-ui;color:#182330}img{max-width:100%}a{color:#166aa2}</style>'
    import html
    body += f'<h1>Figure 2 — round 3</h1><p>Generated {generated}. Producer complete.</p><p><a href="figure2_round3.pdf">Main PDF</a> · <a href="figure2_round3_supplement.pdf">Complete supplement</a> · <a href="README.md">Sources and qualifications</a></p><img src="figure2_round3.png" alt="Figure 2: size and pulse-memory signatures">'
    body += ''.join('<p>'+html.escape(p)+'</p>' for p in legend)
    body += '<div style="text-align:center;margin:1.5em 0"><span style="float:right">(1)</span><math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>x</mi><mo>=</mo><mfrac><msub><mi>s</mi><mrow><mtext>gap</mtext><mn>1</mn></mrow></msub><msub><mi>s</mi><mn>1</mn></msub></mfrac><mo>,</mo><mspace width="2em"/><mi>y</mi><mo>=</mo><mfrac><mrow><msub><mi>s</mi><mn>2</mn></msub><mo>−</mo><msub><mi>s</mi><mn>1</mn></msub></mrow><mrow><mn>1</mn><mo>−</mo><msub><mi>s</mi><mn>1</mn></msub></mrow></mfrac></math></div><p>s₁ and s₂ are pulse-end challenger shares. The first-gap share uses the existing visibility gate in D; raw shares remain in the CSV.</p>'
    (OUT / 'index.html').write_text(body + '</html>')
    write(OUT / 'supporting/metadata/source_sha256.json', SOURCES)
    write(OUT / 'supporting/metadata/build_receipt.json', dict(generated_utc=generated,
        native_runs_submitted=0, kinetics_retuned=False, source_states=97, width_conditions=len(widths),
        temporal_conditions=len(endpoints), deBelly_bounds=6, class2_name='Saturating activation',
        source_script_sha256=sha(Path(__file__)), completed=True))


def main():
    for directory in ('supporting/data', 'supporting/vectors', 'supporting/metadata', 'supporting/viewers'):
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    data = extract()
    render(data)
    supplements(data)
    content(data)
    print(json.dumps(dict(main=str(OUT / 'figure2_round3.pdf'), width_conditions=450, temporal_conditions=180)))


if __name__ == '__main__':
    main()
