"""Render the corrected manuscript protocol without historical reversal results.

Cue curves show the exact cosine-patch geometry, not simulated cell profiles.
Matched reversal results remain pending; no numerical length trend is invented.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
BLUE, ORANGE, GRAY = '#0072B2', '#D55E00', '#666666'
MODELS = ['wave_pinning', 'goryachev', 'debelly', 'holmes_model3', 'spring']
LABELS = dict(zip(MODELS, ['Wave-pinning', 'Goryachev', 'deBelly', 'Holmes 3', 'Spring']))
CONDITIONS = ['fixed_10pct', 'fixed_20pct', 'relative_20pct']
CONDITION_LABELS = ['Fixed width 1.5', 'Fixed width 3', 'Relative width 20%']
GATES = [('orientation', BLUE, 'o', 'Orientation gate'),
         ('full_polarity', ORANGE, 's', 'Polarity strength within 1%')]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def box(ax, x, y, width, height, title, text, color=BLUE):
    ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle='round,pad=0.006',
                               facecolor='white', edgecolor=color, linewidth=1.3))
    ax.text(x + width / 2, y + height - .052, title, ha='center', va='top',
            weight='bold', fontsize=9, color=color)
    ax.text(x + width / 2, y + height / 2 - .025, text, ha='center', va='center',
            fontsize=8, linespacing=1.55)


def arrow(ax, a, b, color=GRAY):
    ax.annotate('', xy=b, xytext=a, arrowprops=dict(arrowstyle='->', color=color, lw=1.4))


def cosine(q, width, end='left'):
    distance = q if end == 'left' else 1 - q
    z = distance / width
    return np.where((z >= 0) & (z <= 1), .5 * (1 + np.cos(np.pi * np.clip(z, 0, 1))), 0)


def render(data, output, campaign=None):
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = data / 'establishment_basal' / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    assert manifest['standardized_reference_length'] == 15
    assert manifest['spring_reference_length'] == 1
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(.065, .967, 'Basal preparation and matched finite reversal', fontsize=16, weight='bold')
    fig.text(.065, .936, 'Five models · consistent cosine patches · verified cue-free polarized starting states', fontsize=10)
    fig.text(.065, .908, 'Corrected preparation is complete; matched reversal measurements and their length trends are pending.',
             fontsize=9, color=ORANGE)

    fig.text(.025, .861, 'A', weight='bold', fontsize=12)
    fig.text(.065, .861, 'Establish polarity from independent native rest', weight='bold', fontsize=11)
    a = fig.add_axes([.065, .634, .90, .205]); a.set_axis_off()
    box(a, .01, .18, .25, .69, 'Native resting state', 'No sustained preconditioning\nNew resting start for each trial')
    box(a, .365, .18, .25, .69, 'Finite left-end patch', 'Search the applied duration\nHold peak and width rule fixed')
    box(a, .72, .18, .25, .69, 'Withdraw added cue', 'Restore native basal input\nVerify the full native endpoint')
    arrow(a, (.27, .52), (.35, .52)); arrow(a, (.63, .52), (.705, .52))
    a.text(.5, .035, 'Select the strongest observed persistent left polarity, then the shortest tested duration within 1% of that strength.',
           ha='center', fontsize=8.5)

    fig.text(.025, .608, 'B', weight='bold', fontsize=12)
    fig.text(.065, .608, 'Three cue widths; challenge is the spatial mirror of preparation', weight='bold', fontsize=11)
    q = np.linspace(0, 1, 501)
    for j, title in enumerate(['Fixed width 1.5', 'Fixed width 3', 'Width 20% of cell length']):
        ax = fig.add_axes([.065 + .317 * j, .435, .265, .149])
        for length, alpha, ls in [(15, 1, '-'), (30, .48, '--')]:
            width = ([1.5, 3][j] / length) if j < 2 else .2
            ax.plot(q, cosine(q, width, 'left'), color=BLUE, alpha=alpha, ls=ls, lw=1.5)
            ax.plot(q, cosine(q, width, 'right'), color=ORANGE, alpha=alpha, ls=ls, lw=1.5)
        ax.set(xlim=(0, 1), ylim=(-.025, 1.09), xticks=[0, .5, 1], yticks=[0, 1],
               xlabel='Position / instantaneous length', ylabel='Added cue / peak A' if j == 0 else '')
        ax.set_title(title, fontsize=10, loc='left')
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
    handles = [Line2D([0], [0], color=BLUE, label='Left preparation'),
               Line2D([0], [0], color=ORANGE, label='Right challenge'),
               Line2D([0], [0], color=GRAY, label='Length 15'),
               Line2D([0], [0], color=GRAY, alpha=.48, ls='--', label='Length 30')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .383), ncol=4,
               frameon=False, fontsize=8)
    fig.text(.065, .344, 'Peak A stays independent of length in every condition; mirror the same width and amplitude in the reversal trial.', fontsize=8)
    fig.text(.065, .323, 'Spring: native reference length 1, fixed widths 0.1 and 0.2; relative width follows its moving physical domain.', fontsize=8)

    fig.text(.025, .282, 'C', weight='bold', fontsize=12)
    fig.text(.065, .282, 'Reversal starts from the selected, verified basal endpoint', weight='bold', fontsize=11)
    c = fig.add_axes([.065, .130, .90, .129]); c.set_axis_off()
    box(c, .01, .12, .25, .77, 'Saved left-polarized state', 'Exact full native checkpoint\nNo reflection or fresh preparation')
    box(c, .365, .12, .25, .77, 'Finite right-end patch', 'Same peak and spatial width\nVary challenge duration', ORANGE)
    box(c, .72, .12, .25, .77, 'Basal recovery', 'Withdraw added patch\nTest persistent opposite polarity', ORANGE)
    arrow(c, (.27, .52), (.35, .52)); arrow(c, (.63, .52), (.705, .52), ORANGE)
    fig.text(.065, .097, 'Wave-pinning · Goryachev · deBelly · Holmes 3 · Spring', fontsize=9, weight='bold')
    fig.text(.065, .073, 'Otsuji and pure LEGI are excluded from the selected reversal suite; corrected LEGI withdrawal removes its polarity.', fontsize=8.5)
    fig.text(.065, .050, 'Unpolarized preparations are ineligible. A selected zero duration means no polarizing cue was required. Native units differ across models.', fontsize=8)
    fig.text(.065, .027, 'Cue curves illustrate the protocol only. No reversal-response trajectories, response slopes, or mechanisms are inferred here.', fontsize=8, color=GRAY)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(output / ('main_polarity.' + ext), dpi=240, facecolor='white')
    plt.close(fig)
    (output / 'CAPTIONS.md').write_text('''# Main figure: corrected basal preparation and matched finite reversal

A. Every establishment trial starts independently at the model's native resting
state. A finite left-end cosine patch is withdrawn to restore native basal
conditions. Full native fields are continued and checked against the recorded
convergence and holdout criteria. Among verified observable cue-aligned endpoints,
select the strongest observed polarity, then the shortest tested duration within
1% of its strength. This is a bounded sampled search, not a global optimum.

B. The reference length is 15 for the six physical-length models in the full
establishment comparison; Spring uses its native reference 1 as a comparison
convention, without a calibrated physical-unit conversion. Fixed widths are
1.5 and 3 (Spring 0.1 and 0.2); the relative condition uses 20% of the actual
instantaneous domain. Peak amplitude is independent of length. Drawn cue curves
show the exact cosine shape at lengths 15 and 30, normalized by peak A.
They are cue geometry, not simulated active-species profiles. Fixed patches
on shorter cells are restricted to the actual domain without width or peak
rescaling; no half-cell restriction from the previous pilot is imposed.

C. The planned matched finite-countercue experiment reuses the qualifying full
native establishment checkpoint and mirrors the cue spatially with identical
amplitude and width rule. The challenge is temporary; withdrawal returns to
native basal conditions before testing persistent opposite polarity. No
reflection, sustained linear conditioning, or old release state is used.
The selected five-model reversal scope is wave-pinning, Goryachev, deBelly,
Holmes 3 and Spring; Otsuji and LEGI are excluded by the user's scope decision.
Pure LEGI returns to a homogeneous state with positive
basal stimulus, so all 462 corrected establishment combinations are unpolarized
and provide no eligible starting state. A selected zero establishment duration
means spontaneous polarization or no cue needed; it remains explicitly recorded.

The corrected establishment comparison is complete. Matched reversal measurements,
their tested duration grid, final readout definitions and length-response results
must come from the new campaign, not the superseded historical width pilot.
This protocol-only main is a reviewable interim artifact and claims no response
trend or mechanistic explanation. Detailed establishment heatmaps are retained
as supplements; the old reversal-field report is a personal diagnostic.
''')
    derived = data / 'derived'; derived.mkdir(parents=True, exist_ok=True)
    provenance = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer='finished',
                      status='corrected_protocol_only; matched_reversal_pending',
                      integrated_new_trajectories=False,
                      input_hashes={'establishment_basal/manifest.json': digest(manifest_path)},
                      renderer_sha256=digest(Path(__file__)),
                      output_hashes={f'main_polarity.{ext}': digest(output / f'main_polarity.{ext}')
                                     for ext in ['png', 'pdf', 'svg']})
    (derived / 'figure_generation_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    if campaign is not None and render_complete_campaign(data, output, campaign):
        return
    print(json.dumps({'output': str(output), 'status': provenance['status']}, indent=2))


def read_endpoint_rows(path):
    """Keep the CSV schema explicit: no historical trial or source imports."""
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for name in ['length', 'width', 'amplitude', 'initial_orientation', 'endpoint',
                     'lower_failing_duration', 'relative_bracket_width', 'maximum_tested_duration']:
            value = row.get(name)
            row[name] = float(value) if value not in [None, '', 'None', 'null'] else None
        value = (row.get('selected_field_tail_pass') or '').lower()
        row['selected_field_tail_pass'] = True if value in ['true', '1'] else False if value in ['false', '0'] else None
        assert row['model'] in MODELS, row['model']
        assert row['condition'] in CONDITIONS, row['condition']
        assert row['definition'] in ['orientation', 'full_polarity']
        assert row['status'] in ['success', 'censored', 'ineligible', 'preparation_invalid'], 'Final endpoint table contains an incomplete outcome'
        if row['status'] == 'success':
            assert row['endpoint'] is not None and row['endpoint'] > 0
    paired = {(r['model'], r['length'], r['condition'], r['definition']): r for r in rows}
    assert len(paired) == len(rows), 'Duplicate endpoint rows'
    for row in rows:
        if row['definition'] == 'full_polarity' and row['status'] == 'success':
            original = paired[row['model'], row['length'], row['condition'], 'orientation']
            assert original['status'] == 'success', 'Full polarity must imply orientation gate'
            assert row['endpoint'] >= original['endpoint'] - 1e-9
    return rows


def fit_summary(rows):
    result = []
    for model in MODELS:
        for condition in CONDITIONS:
            for gate, _, _, _ in GATES:
                chosen = [r for r in rows if r['model'] == model and r['condition'] == condition and r['definition'] == gate]
                assert chosen, (model, condition, gate)
                scalar_passed = [r for r in chosen if r['status'] == 'success']
                passed = [r for r in scalar_passed if r['selected_field_tail_pass'] is True]
                row = dict(model=model, condition=condition, definition=gate,
                           n_total=len(chosen), n_success=len(scalar_passed), n_fit=len(passed),
                           n_censored=sum(r['status'] == 'censored' for r in chosen),
                           n_ineligible=sum(r['status'] == 'ineligible' for r in chosen),
                           n_preparation_invalid=sum(r['status'] == 'preparation_invalid' for r in chosen),
                           n_field_tail_failed=sum(r['selected_field_tail_pass'] is False for r in scalar_passed),
                           n_field_tail_unknown=sum(r['selected_field_tail_pass'] is None for r in scalar_passed),
                           fit_status='insufficient_successes', slope=None, intercept=None,
                           length_min=None, length_max=None, characteristic_length=None,
                           fitted_response_at_midpoint=None, normalized_slope=None, r_squared=None)
                if len(passed) > 3:
                    x = np.asarray([r['length'] for r in passed])
                    y = np.asarray([r['endpoint'] for r in passed])
                    slope, intercept = np.linalg.lstsq(np.column_stack([x, np.ones_like(x)]), y, rcond=None)[0]
                    if np.ptp(y) == 0:
                        slope, intercept = 0., y[0]
                    midpoint = (x.min() + x.max()) / 2
                    fitted = slope * midpoint + intercept
                    variation = np.sum((y - y.mean()) ** 2)
                    row.update(fit_status='fit' if fitted > 0 else 'nonpositive_midpoint',
                               slope=float(slope), intercept=float(intercept),
                               length_min=float(x.min()), length_max=float(x.max()),
                               characteristic_length=float(midpoint), fitted_response_at_midpoint=float(fitted),
                               normalized_slope=float(slope * midpoint / fitted) if fitted > 0 else None,
                               r_squared=float(1 - np.sum((y - slope * x - intercept) ** 2) / variation) if variation > 0 else None)
                result.append(row)
    return result


def save(fig, output, name, dpi=220):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(output / f'{name}.{ext}', dpi=dpi, facecolor='white')
    plt.close(fig)


def measured_row(fig, grid, index, letter, title, field, duration, arrival=None, failed=None):
    """Cue, recorded profiles, scalar clock; snapshots use actual stored samples."""
    cueax, profax, respax = [fig.add_subplot(grid[index, j]) for j in range(3)]
    fig.text(.025, cueax.get_position().y1 + .015, letter, weight='bold', fontsize=11)
    fig.text(.075, cueax.get_position().y1 + .015, title, weight='bold', fontsize=9)
    x, times, active = field['x'], field['time'], field['active']
    q = x / 15
    stop = float(field['orientation_time'][-1])
    end = 'left' if index == 0 else 'right'
    cue_times = np.unique(np.r_[np.linspace(0, duration, 21), duration,
                               np.geomspace(max(duration, .001), stop, 201)])
    values = np.asarray([cosine(q, .2, end) if t < duration else np.zeros_like(q) for t in cue_times])
    cueax.pcolormesh(cue_times, q, values.T, cmap='Greys', norm=matplotlib.colors.Normalize(0, 1),
                     shading='nearest', rasterized=True)
    cueax.set(xlim=(0, stop), ylim=(0, 1), yticks=[0, 1], ylabel='Position / length')
    targets = [0, duration, duration + 30, duration + 250, stop] if index == 0 else [0, duration, duration + 100, duration + 1000, stop]
    indices = sorted(set(int(np.argmin(np.abs(times - t))) for t in targets))
    chosen = active[indices]
    minimum, span = chosen.min(), np.ptp(chosen)
    colors = plt.cm.viridis(np.linspace(.08, .9, len(indices)))
    for k, color in zip(indices, colors):
        profax.plot(q, (active[k] - minimum) / span, color=color, lw=1.3)
        cueax.axvline(times[k], color=color, ls=':', lw=.7)
        respax.scatter(times[k], np.interp(times[k], field['orientation_time'], field['orientation']),
                       s=18, color=color, edgecolor='white', lw=.4, zorder=4)
    profax.set(xlim=(0, 1), ylim=(0, 1.08), xticks=[0, .5, 1], yticks=[0, 1], ylabel='Scaled active species')
    profax.legend(handles=[Line2D([0], [0], color=c, label=f'{times[k]:.3g}') for k, c in zip(indices, colors)],
                  ncol=len(indices), loc='upper center', bbox_to_anchor=(.5, -.26), frameon=False,
                  fontsize=6.7, handlelength=.8, handletextpad=.3, columnspacing=.7)
    respax.plot(field['orientation_time'], field['orientation'], color='.2', lw=1.1)
    respax.axhline(0, color='.7', ls=':', lw=.6)
    respax.set(ylim=(-.52, .64), ylabel='Orientation', xlim=(0, stop))
    if failed is not None:
        respax.plot(failed['time'], failed['orientation'], color='.5', ls='--', lw=1)
        respax.legend(handles=[Line2D([0], [0], color='.2', label=f'{duration:.3g} passes full gate'),
                               Line2D([0], [0], color='.5', ls='--', label=f'{failed["duration"]:.3g} fails full gate')],
                       frameon=False, fontsize=6.7, loc='upper left')
        respax.axvspan(stop - 500, stop, color='.92', zorder=0)
    if arrival is not None:
        respax.axvspan(arrival, stop, color='.92', zorder=0)
        respax.text(.98, .92, f'Cue duration {duration:.3g}\nTotal establishment {arrival:.3g}',
                    ha='right', va='top', transform=respax.transAxes, fontsize=7)
    for ax in [cueax, respax]:
        ax.set_xscale('symlog', linthresh=1 if index == 0 else max(1, duration), linscale=1)
        ax.axvline(duration, color='.4', ls='--', lw=.7)
        ax.set_xlabel('Time from patch onset', fontsize=8)
        ax.set_xticks([0, duration, 100, 1000, stop])
        ax.set_xticklabels(['0', f'{duration:.3g}', '100', '1,000', f'{stop:,.0f}'], fontsize=7)
    for ax in [cueax, profax, respax]:
        ax.tick_params(labelsize=7)
        ax.yaxis.label.set_size(8)


def preparation_example(data):
    base = data / 'establishment_basal'
    with (base / 'results/endpoints.csv').open() as handle:
        row = next(r for r in csv.DictReader(handle) if r['model'] == 'wave_pinning' and
                   float(r['length']) == 15 and r['condition'] == 'fixed_20pct' and float(r['amplitude_factor']) == 1)
    checkpoint = base / row['selected_checkpoint']
    assert digest(checkpoint) == row['selected_checkpoint_sha256']
    trace = checkpoint.with_name(checkpoint.name.replace('.checkpoint', ''))
    duration = float(row['cue_duration'])
    with np.load(trace) as z:
        field = dict(x=z['x'].copy(),
                     time=np.r_[z['pulse_time'][:-1], duration + z['profile_time']],
                     active=np.r_[z['pulse_active'][:-1], z['profile_active']],
                     orientation_time=np.r_[z['pulse_time'][:-1], duration + z['release_time']],
                     orientation=np.r_[z['pulse_orientation'][:-1], z['release_orientation']])
    return field, duration, float(row['total_establishment_time']), trace


def complete_main(data, output, rows, fits, example_path):
    prep, duration, arrival, prep_path = preparation_example(data)
    with np.load(example_path) as z:
        example = {k: z[k].copy() for k in z.files}
    field = {k: example[k] for k in ['time', 'x', 'active', 'orientation']}
    field['orientation_time'] = field['time']
    selected = float(example['passing_duration'])
    selected_row = next(r for r in rows if r['model'] == 'wave_pinning' and r['length'] == 15
                        and r['condition'] == 'fixed_20pct' and r['definition'] == 'full_polarity')
    assert selected_row['status'] == 'success' and abs(selected_row['endpoint'] - selected) < 1e-9
    failed = None
    if 'failed_time' in example:
        failed = dict(time=example['failed_time'], orientation=example['failed_orientation'], duration=float(example['failed_duration']))
    fig = plt.figure(figsize=(11, 10.5))
    fig.text(.075, .974, 'Matched basal polarity reversal and length response', fontsize=15, weight='bold')
    fig.text(.075, .949, 'Five models · nominal peak amplitudes · independent finite-patch establishment · mirrored finite challenge', fontsize=9)
    fig.text(.075, .924, 'A–B: recorded wave-pinning, length 15, width 3, peak 0.1 · native time units', fontsize=8.5)
    for x, label in zip([.075, .385, .695], ['Added cue: white 0 → black A', 'Active-species profiles vs x / L', 'Measured polarity']):
        fig.text(x, .893, label, weight='bold', fontsize=9)
    grid = fig.add_gridspec(2, 3, left=.075, right=.975, top=.865, bottom=.56, hspace=.65, wspace=.36)
    measured_row(fig, grid, 0, 'A', 'Finite preparation patch → basal verified endpoint', prep, duration, arrival)
    measured_row(fig, grid, 1, 'B', 'Mirrored finite countercue → basal recovery', field, selected, failed=failed)
    fig.text(.025, .484, 'C', weight='bold', fontsize=11)
    fig.text(.075, .484, 'The three length-dependent cue geometries', weight='bold', fontsize=9)
    for j, label in enumerate(CONDITION_LABELS):
        ax = fig.add_axes([.075, .365 - .086 * j, .27, .046])
        distance = np.linspace(0, .4, 400)
        for length, color, ls in [(15, '#333333', '-'), (30, '#73A6B8', '--')]:
            width = ([1.5, 3][j] / length) if j < 2 else .2
            ax.plot(distance, cosine(distance, width), color=color, ls=ls, lw=1.3)
        ax.set(xlim=(0, .4), ylim=(0, 1.1), yticks=[0, 1], xticks=[])
        ax.set_title(label, fontsize=7.5, loc='left', pad=2)
        ax.tick_params(labelsize=6.5)
        if j == 0:
            ax.legend(handles=[Line2D([0], [0], color='#333333', label='L=15'),
                               Line2D([0], [0], color='#73A6B8', ls='--', label='L=30')], frameon=False, fontsize=6.5)
        if j == 2:
            ax.set_xticks([0, .2, .4]); ax.set_xlabel('Distance from stimulated end / L', fontsize=7)
    fig.text(.075, .434, 'Added cue / peak A; peak stays fixed', fontsize=7)
    fig.text(.075, .119, 'Spring: native reference 1, fixed widths\n0.1 and 0.2; relative width follows\nits instantaneous physical domain.', fontsize=7, linespacing=1.5)
    fig.text(.410, .484, 'D', weight='bold', fontsize=11)
    fig.text(.485, .484, 'Descriptive length trends', weight='bold', fontsize=9)
    fig.text(.485, .448, 'OLS slopes; n = stable successes / scalar successes', fontsize=7.5)
    lookup = {(r['model'], r['condition'], r['definition']): r for r in fits}
    matrix = np.full((5, 6), np.nan)
    for i, model in enumerate(MODELS):
        for j, condition in enumerate(CONDITIONS):
            for k, (gate, _, _, _) in enumerate(GATES):
                value = lookup[model, condition, gate]['normalized_slope']
                if value is not None: matrix[i, 2*j+k] = value
    cmap = plt.colormaps['RdBu_r'].copy(); cmap.set_bad('#eeeeee')
    maximum = max(1, float(np.nanmax(np.abs(matrix)))) if np.any(np.isfinite(matrix)) else 1
    ax = fig.add_axes([.485, .195, .49, .198])
    im = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=-maximum, vmax=maximum,
                   aspect='auto', interpolation='nearest')
    for i, model in enumerate(MODELS):
        for j, condition in enumerate(CONDITIONS):
            for k, (gate, _, _, _) in enumerate(GATES):
                row = lookup[model, condition, gate]; value = row['normalized_slope']
                star = '*' if row['n_field_tail_failed'] or row['n_field_tail_unknown'] else ''
                counts = f'n={row["n_fit"]}/{row["n_success"]}'
                label = f'{value:+.2f}{star}\n{counts}' if value is not None else f'No fit\n{counts}'
                ax.text(2*j+k, i, label, ha='center', va='center', fontsize=7,
                        color='white' if value is not None and abs(value) > .65 * maximum else '#222222')
    ax.set_yticks(range(5), [LABELS[m] for m in MODELS], fontsize=8)
    ax.set_xticks(range(6), ['Orient.', 'Strength'] * 3, fontsize=7)
    ax.tick_params(length=0)
    for j, title in enumerate(['Fixed 1.5', 'Fixed 3', 'Relative 20%']):
        ax.text((2*j+1) / 6, 1.055, title, transform=ax.transAxes, ha='center', fontsize=7.5)
    for border in [1.5, 3.5]: ax.axvline(border, color='white', lw=2)
    cax = fig.add_axes([.485, .147, .25, .009])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.ax.tick_params(labelsize=6.5, length=2)
    fig.text(.755, .145, '* Some scalar successes excluded by\n   the native-field tail diagnostic.', fontsize=6.7, va='center')
    fig.text(.075, .095, 'Normalized slope = raw slope × valid-length midpoint / fitted duration there. Fits require scalar success AND passing native-field tail stability.', fontsize=7)
    fig.text(.075, .076, 'Gray: fewer than four stable successes or nonpositive fitted midpoint. Scalar-only excluded cases remain visible in counts and Supplement 1.', fontsize=7)
    fig.text(.075, .057, 'Full strength is a scalar polarity criterion, not complete pattern replacement. Duration minima and OLS trends do not establish scaling laws.', fontsize=7)
    fig.text(.075, .038, 'Preparation arrival excludes verification time; reversal readout excludes preparation and 10,000-unit basal recovery. Otsuji and LEGI are excluded.', fontsize=7)
    save(fig, output, 'main_polarity')
    return prep_path


def complete_supplement(rows, output, stamp):
    fig, axes = plt.subplots(5, 3, figsize=(14.5, 15.5))
    for i, model in enumerate(MODELS):
        for j, condition in enumerate(CONDITIONS):
            ax = axes[i, j]
            chosen = sorted([r for r in rows if r['model'] == model and r['condition'] == condition], key=lambda r: (r['length'], r['definition']))
            valid = [r['endpoint'] for r in chosen if r['status'] == 'success']
            caps = [r['maximum_tested_duration'] for r in chosen if r['maximum_tested_duration'] is not None]
            assert caps, 'Censoring limits must be supplied explicitly by the campaign'
            maximum = max(caps)
            for gate, color, marker, _ in GATES:
                gate_rows = [r for r in chosen if r['definition'] == gate]
                ax.plot([r['length'] for r in gate_rows], [r['endpoint'] if r['status'] == 'success' else np.nan for r in gate_rows],
                        color=color, ls='-' if gate == 'orientation' else '--', lw=1.3)
                for row in gate_rows:
                    if row['status'] == 'success':
                        ax.plot(row['length'], row['endpoint'], marker, ms=4, color=color,
                                mfc=color if row['selected_field_tail_pass'] is True else 'white')
                    elif row['status'] == 'censored':
                        ax.plot(row['length'], row['maximum_tested_duration'], '^', ms=4.5, mfc='white', mec=color)
                    elif gate == 'orientation':
                        ax.plot(row['length'], -.13, 'x' if row['status'] == 'ineligible' else '+', color='#999999',
                                transform=ax.get_xaxis_transform(), clip_on=False, ms=5)
            x = [r['length'] for r in chosen]; delta = max(x) - min(x)
            ax.set_xlim(min(x)-.03*delta, max(x)+.03*delta)
            ax.set_yscale('log')
            censored = any(r['status'] == 'censored' for r in chosen)
            ax.set_ylim(min(valid)/1.5 if valid else maximum/100,
                        max(max(valid) if valid else maximum, maximum if censored else 0)*1.6)
            if i == 0:
                ax.set_title(CONDITION_LABELS[j] + ('\n(Spring: 0.1)' if j == 0 else '\n(Spring: 0.2)' if j == 1 else ''), fontsize=10)
            if j == 0: ax.set_ylabel(LABELS[model]+'\nShortest passing tested pulse', fontsize=9)
            ax.set_xlabel('Native resting length L0' if model == 'spring' else 'Nominal cell length (native units)', fontsize=8, labelpad=20)
            ax.tick_params(labelsize=7); ax.grid(axis='y', color='#eeeeee', lw=.6)
            initial = [r for r in chosen if r['definition'] == 'orientation']
            ax.text(.98, .035, f'{sum(r["status"] == "ineligible" for r in initial)} ineligible; '
                    f'{sum(r["status"] == "preparation_invalid" for r in initial)} invalid preparation',
                    transform=ax.transAxes, ha='right', fontsize=6.5, color=GRAY)
    fig.suptitle('Supplement 1 · Corrected finite countercue: width and response definitions', fontsize=14, y=.98)
    fig.text(.5, .956, 'Nominal amplitude · qualifying finite-patch establishment checkpoint · identical mirrored challenge · basal recovery', ha='center', fontsize=9)
    handles = [Line2D([0], [0], color=c, marker=m, ls='-' if g == 'orientation' else '--', label=t) for g,c,m,t in GATES]
    handles.extend([Line2D([0],[0],color=GRAY,marker='^',mfc='white',ls='none',label='No qualifying tested pulse'),
                    Line2D([0],[0],color=GRAY,marker='x',ls='none',label='Ineligible preparation (below axis)'),
                    Line2D([0],[0],color=GRAY,marker='+',ls='none',label='Invalid preparation (below axis)'),
                    Line2D([0],[0],color=GRAY,marker='o',mfc='white',ls='none',label='Native-field tail diagnostic fails')])
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5,.937), ncol=3, fontsize=8, frameon=False)
    fig.text(.5, .034, 'Fixed widths use reference length 15 (Spring native reference 1). Relative width follows instantaneous physical length; peak is independent of length.\n'
             'Both scalar gates score the same tested durations. Full strength additionally requires polarity within 1% of the reversed initial magnitude in the final 500 units.\n'
             'Open markers distinguish native-field instability from scalar success. Censoring marks tested limits, not measured durations. Missing endpoints break lines.\n'
             'Native time units and y ranges differ by model. Otsuji and LEGI are excluded by scope. Source generation: '+stamp+'; producer finished.',
             ha='center', fontsize=7.5, linespacing=1.6)
    fig.subplots_adjust(left=.09,right=.985,top=.865,bottom=.135,hspace=.5,wspace=.22)
    save(fig, output, 'supplement_width_response', dpi=180)


def render_complete_campaign(data, output, campaign):
    results = campaign / 'results' if (campaign / 'results').exists() else campaign
    status_path = results / 'status.json'
    if not status_path.exists(): return False
    status = json.loads(status_path.read_text())
    if not status.get('complete', False): return False
    assert status.get('passed', True), 'Completed campaign failed its audits'
    validation_path = results / 'validation.json'
    validation = json.loads(validation_path.read_text())
    assert validation['passed'] and validation['complete_execution'], 'Final figures require passed execution and scientific readout audits'
    endpoint_path = results / 'endpoints.csv'
    example_path = results / 'explanation_wave_pinning.npz'
    assert example_path.exists(), 'Final main requires the recorded corrected example'
    rows = read_endpoint_rows(endpoint_path)
    fits = fit_summary(rows)
    summary_path = data / 'derived/corrected_length_trends.csv'
    with summary_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fits[0])); writer.writeheader(); writer.writerows(fits)
    prep_path = complete_main(data, output, rows, fits, example_path)
    complete_supplement(rows, output, status['generated_at'])
    caption = output / 'CAPTIONS.md'
    caption.write_text('''# Main figure: matched basal polarity reversal and length response

A. Recorded wave-pinning at length 15, fixed width 3 and nominal peak 0.1.
Independent native rest is followed by the selected finite left cosine patch,
withdrawal to native basal conditions, and verified polarized endpoint. The
establishment readout is applied duration plus first persistent post-release
1% polarity arrival; convergence confirmation and holdout are excluded from
that clock. The right challenge starts from the exact fully verified checkpoint.
Snapshots use actual recorded times. Profile curves share a minimum/range
within each row; colors match cue map lines and orientation markers.

B. The selected shortest tested right-patch duration passes the full scalar
polarity gate after 10000 native units of basal recovery. The dashed neighboring
lower-duration trial fails that full gate (it can still pass orientation alone).
Gray shading marks the final 500-unit readout interval. Applied duration excludes
preparation and recovery; native-field tail stability is an independent diagnostic.
The orientation gate requires every original 2.5-unit sample in that final
interval to have right orientation above 0.05 and active contrast at least 0.05,
with orientation range below 0.005. The full gate additionally requires every
sample's polarity within 1% of the magnitude of initial left polarity.
The native-field diagnostic requires each monitored field's maximum pointwise
tail range and its last sampled rate projected over 500 units to be at most
0.1% of its fixed field scale. These conditions measure sampled persistence,
not an exact mathematical equilibrium or identical complete pattern.

C. Identical left/right cosine patches use fixed widths 1.5 and 3 at reference
length 15, or 20% of instantaneous physical cell length. Peak remains fixed with
length. Nominal added peaks are 0.1 for wave-pinning and Goryachev, 3.6 for
deBelly, 1 for Holmes 3 and 4 for Spring, in their respective native units.
Spring uses native reference 1 and fixed widths 0.1 and 0.2, with the
relative patch following its moving domain; no physical-unit conversion is implied.
On short domains the unchanged patch is restricted to the actual domain without
renormalizing peak or width. There is no half-cell restriction from the old pilot.

D. Primary ordinary least-squares duration fits require scalar success AND the
selected native-field tail diagnostic to pass at more than three lengths, with
positive fitted duration at the midpoint of that valid range.
Normalized slope is raw slope times that midpoint divided by the fitted duration
there. Both gates are fit separately and can use different length ranges.
n displays native-field-stable scalar successes divided by all scalar successes;
asterisks mark at least one scalar success excluded by the native-field diagnostic.
Unknown diagnostics are also excluded and counted in the CSV. Gray denotes
unavailable fit rather than imputing a response or a trend.
Primary fits filter the selected shortest scalar readout; a longer stable pulse
does not replace an excluded endpoint. Preparation-invalid cases fail their
zero-added-cue control and are retained explicitly outside the fits.
The CSV retains raw slopes, intercepts, R-squared, length ranges, and all outcome
counts. Nonlinear sampled curves, selective fit membership and fitted slopes
do not establish scaling laws. Complete scalar polarity recovery does not mean
every native profile is replaced or converged.

# Supplement 1: width and scalar response definitions

Five selected models, three widths, all established nominal-amplitude lengths.
Blue circles/solid curves are the orientation gate. Orange squares/dashes add
the requirement that the final 500-unit polarity remain within 1% of the reversed
initial magnitude. Both definitions use the same tested trials and recovery
history. Triangles mark censoring at the explicit maximum tested duration;
they are not measured endpoints. Ineligible/invalid starting preparations are
marked below the axis. Open successful markers fail the separate native-field
tail diagnostic. Missing endpoints break lines; time units/y ranges differ by
model. Otsuji and LEGI are outside the user-selected five-model reversal scope.
All reversal results here come from the corrected basal-preparation campaign;
no historical conditioning or old width-pilot readouts are included.
''')
    provenance = dict(generated_at=datetime.now(timezone.utc).isoformat(),
                      producer='finished', source_generated_at=status['generated_at'],
                      status='corrected_basal_matched_reversal_complete',
                      input_hashes={str(p.relative_to(data)): digest(p) for p in
                                    [endpoint_path, status_path, validation_path, example_path, prep_path]},
                      renderer_sha256=digest(Path(__file__)), summary_sha256=digest(summary_path),
                      output_hashes={p.name:digest(p) for p in output.glob('*')
                                     if p.suffix in ['.png','.pdf','.svg'] and
                                     (p.stem == 'main_polarity' or p.stem == 'supplement_width_response')})
    (data / 'derived/figure_generation_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # The SLURM finalizer calls this renderer directly; publish the new links as
    # well as the figure bytes without invoking another rendering cycle.
    build_path = Path(__file__).with_name('build.py')
    if build_path.exists() and output.resolve() == (ROOT / 'figures').resolve():
        spec = importlib.util.spec_from_file_location('_polarity_figure_index', build_path)
        builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(builder)
        builder.index()
    print(json.dumps({'output':str(output), 'status':provenance['status'], 'summary_rows':len(fits)},indent=2))
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT / 'data')
    parser.add_argument('--output', type=Path, default=ROOT / 'figures')
    parser.add_argument('--campaign', type=Path, default=ROOT / 'data/runs/corrected_reversal_20261003')
    args = parser.parse_args()
    render(args.data.resolve(), args.output.resolve(), args.campaign.resolve())
