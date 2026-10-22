"""Selective review of performed native strength/width/pulse experiments."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/width_strength_grid_20261019'
PRIOR = OUT.parent / 'finite_response_candidate_20261019/figure2_organization.png'
DATA = OUT / 'supporting/data'


def readcsv(p):
    with p.open() as f:
        return list(csv.DictReader(f))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def get(rows, model, size, width, multiplier, duration):
    return next(r for r in rows if r['model'] == model and float(r['length']) == size
                and float(r['physical_width']) == width and float(r['native_peak_multiplier']) == multiplier
                and float(r['pulse_duration']) == duration)


def endpoint_templates(rows):
    run = ROOT / 'data/runs/revision_width_strength_grid_20261019'
    manifest = json.loads((run / 'manifest.json').read_text())
    matches = []
    for r in rows:
        c = manifest['cases'][int(r['slot'])]
        with np.load(ROOT / r['actual_recovery_path']) as z:
            final = z['finalstate']
        with np.load(run / c['regional_input']) as z:
            initial = z['initial_state']
        scales = np.array([c['fixed_scales'][k] for k in ('active', 'inactive')])[:, None]
        initial_error = float(np.max(abs(final - initial) / scales))
        reflected_error = float(np.max(abs(final - initial[:, ::-1]) / scales))
        reflected = reflected_error < initial_error
        assert reflected == (r['full_polarity_success'] == 'True')
        matches.append(dict(slot=r['slot'], model=r['model'], length=r['length'],
            physical_width=r['physical_width'], native_peak_multiplier=r['native_peak_multiplier'],
            pulse_duration=r['pulse_duration'], final_time_off=10000.,
            max_field_error_from_initial_fixed_scales=initial_error,
            max_field_error_from_reflected_initial_fixed_scales=reflected_error,
            closer_template='reflected_basal' if reflected else 'original_basal',
            closest_template_max_fixed_scale_error=min(initial_error, reflected_error),
            active_scale=c['fixed_scales']['active'], inactive_scale=c['fixed_scales']['inactive'],
            actual_recovery_path=r['actual_recovery_path'], checkpoint_sha256=c['checkpoint_sha256']))
    path = DATA / 'native_endpoint_template_matches.csv'
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(matches[0]))
        writer.writeheader()
        writer.writerows(matches)
    return path, max(r['closest_template_max_fixed_scale_error'] for r in matches)


def result_note(rows, events):
    run = ROOT / 'data/runs/revision_width_strength_grid_20261019'
    manifest = json.loads((run / 'manifest.json').read_text())
    prior_table = OUT.parent / 'finite_response_candidate_20261019/supporting/data/organization_all97.csv'
    organization = readcsv(prior_table)
    foundation_rows, sensitivity, channels = [], [], []
    for model in ('wave_pinning', 'goryachev'):
        sizes = [9., 19.] if model == 'wave_pinning' else [9., 25.]
        native, prior = [], []
        for size in sizes:
            c = next(c for c in manifest['cases'] if c['model'] == model and c['length'] == size)
            with np.load(run / c['regional_input']) as z:
                u, core = z['initial_state'][0], z['incumbent_core_mask']
                native.append(dict(amount=float(u[core].sum() * size / 201),
                                   width=float(core.sum() * size / 201), peak=float(u.max())))
            prior_row = next(o for o in organization if o['model'] == model and float(o['size_native']) == size)
            assert prior_row['checkpoint_sha256'] == c['checkpoint_sha256']
            prior.append(prior_row)
        f = dict(model=model, small_size=sizes[0], large_size=sizes[1], length_fold=sizes[1] / sizes[0])
        for key in ('amount', 'width', 'peak'):
            f['assay_core_' + key + '_small'] = native[0][key]
            f['assay_core_' + key + '_large'] = native[1][key]
            f['assay_core_' + key + '_large_over_small'] = native[1][key] / native[0][key]
        for name, key in (('excess_amount', 'excess_amount_native'), ('interpolated_width', 'width_native'), ('peak', 'peak_native')):
            f['prior_organization_' + name + '_large_over_small'] = float(prior[1][key]) / float(prior[0][key])
        foundation_rows.append(f)
        for size in sizes:
            for width in (3., 5.):
                for loss in (.1, .25, .5):
                    rr = dict(model=model, size=size, physical_width=width, regional_loss_fraction=loss)
                    for multiplier, prefix in ((.5, 'low_peak'), (1., 'baseline_peak'), (2., 'high_peak')):
                        source_row = max([r for r in rows if r['model'] == model and float(r['length']) == size
                                          and float(r['physical_width']) == width
                                          and float(r['native_peak_multiplier']) == multiplier], key=lambda r: float(r['pulse_duration']))
                        event = next(e for e in events if e['slot'] == source_row['slot'] and e['phase'] == 'cue_on'
                                     and e['regional_depletion_fraction'] and float(e['regional_depletion_fraction']) == loss)
                        assert event['status'] == 'sampled_transition'
                        rr[prefix + '_slot'] = source_row['slot']
                        rr[prefix + '_sampled_lower'] = float(event['lower'])
                        rr[prefix + '_sampled_upper'] = float(event['upper'])
                    rr['fourfold_peak_time_shortening_observed_ratio'] = rr['low_peak_sampled_upper'] / rr['high_peak_sampled_upper']
                    rr['adjacent_bracket_ratio_lower'] = rr['low_peak_sampled_lower'] / rr['high_peak_sampled_upper']
                    rr['adjacent_bracket_ratio_upper'] = rr['low_peak_sampled_upper'] / rr['high_peak_sampled_lower']
                    sensitivity.append(rr)
            for multiplier in (.5, 1., 2.):
                r = max([r for r in rows if r['model'] == model and float(r['length']) == size
                         and float(r['physical_width']) == 3.
                         and float(r['native_peak_multiplier']) == multiplier], key=lambda r: float(r['pulse_duration']))
                c = manifest['cases'][int(r['slot'])]
                with np.load(run / c['regional_input']) as z:
                    signal, roi = z['signal'], z['cue_footprint_mask']
                with np.load(ROOT / r['actual_cue_path']) as z:
                    time, inactive = z['time'], z['field_inactive']
                    early = time <= 10. + 1e-12
                    realized = np.sum(inactive[early] * signal, axis=1) * size / 201
                    t10 = int(np.argmin(abs(time - 10.)))
                    channels.append(dict(model=model, size=size, physical_width=3.,
                        native_peak_multiplier=multiplier, requested_peak=c['response_amplitude'],
                        native_window_start=0., native_window_end=10., sampled_interval=.1,
                        sampled_trapezoid_integrated_channel=float(np.trapezoid(realized, time[early])),
                        initial_spatially_integrated_channel=float(realized[0]),
                        channel_at10_over_initial=float(realized[-1] / realized[0]),
                        cue_roi_inactive_at10_over_initial=float(inactive[t10, roi].mean() / inactive[0, roi].mean()),
                        source_slot=r['slot'], actual_cue_path=r['actual_cue_path']))
    foundation_path = DATA / 'incumbent_size_ratios.csv'
    sensitivity_path = DATA / 'finite_loss_strength_sensitivity.csv'
    channel_path = DATA / 'native_channel_first10.csv'
    for path, values in ((foundation_path, foundation_rows), (sensitivity_path, sensitivity), (channel_path, channels)):
        with path.open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)
    ranges = []
    for model in ('wave_pinning', 'goryachev'):
        for width in (3., 5.):
            values = [r['fourfold_peak_time_shortening_observed_ratio'] for r in sensitivity
                      if r['model'] == model and r['physical_width'] == width]
            ranges.append(f'| {"WP" if model == "wave_pinning" else "Goryachev"} | {width:g} | {min(values):.3f}–{max(values):.3f}× |')
    (OUT / 'RESULTS.md').write_text('''# What the performed experiment adds

Larger basals contain more incumbent material, organized differently. In the unchanged prior measurements, WP L9→L19 increases incumbent excess amount **2.27×** and width **2.07×**, with only **1.02×** peak concentration. Goryachev L9→L25 increases excess amount **2.99×** and peak **2.89×**, while width changes only **1.03×**. The new assay's fixed native-cell core gives the same organizational reading: activity amount/width/peak ratios are **2.24/2.08/1.02** for WP and **3.02/1.04/2.89** for Goryachev. These are different size ranges; they are not matched growth factors. [Exact values and both region conventions](supporting/data/incumbent_size_ratios.csv).

Under the same cue grid within each model, the larger incumbent is harder to replace. At WP pulse 60, all six L9 peak/width combinations pass the retained scalar-strength gate, versus three of six L19 combinations; doubling the baseline peak at width 3, or widening the baseline cue to 5, makes L19 pass. At Goryachev pulse 400, all six L9 combinations pass and all six L25 combinations return. Pulse 1600 makes all six L25 combinations pass. The sparse pulse grid can miss shifted success boundaries, so this is not evidence that Goryachev strength is irrelevant or duration is its only control.

Finite weakening times add a continuous response beyond that coarse fate grid. Increasing native peak fourfold (0.5×→2× baseline) shortens time to substantial 10%, 25% or 50% original-region activity loss by these factors, over both tested sizes:

| Model | Physical width | Observed weakening-time shortening |
|---|---:|---:|
''' + '\n'.join(ranges) + '''

Thus the tested WP responses are more sensitive to peak strength in degree across both widths and all three loss thresholds. This is not yet a qualitative mechanism-specific signature. The native input operating points differ, the chosen size folds differ, and neither actuator nor clock is optically calibrated. [All 24 comparisons and adjacent-sample ratio brackets](supporting/data/finite_loss_strength_sensitivity.csv).

The weaker Goryachev strength sensitivity may reflect the tested dose range and saturation as well as the maintenance mechanism: matching fractional peak changes does not match effective actuation. At width 3 over the first 10 native cue units, a fourfold requested-peak increase produces only **1.20–1.21×** the sampled integral of the realized stimulus × inactive channel in Goryachev, versus **3.69×** in WP. Cue-region inactive concentration at time 10 is **6–19%** of initial in Goryachev and **90–97%** in WP. This is consistent with depletion-limited delivery at the chosen Goryachev operating point. It is a description of a sampled native reaction-channel integral, not net active gain, calibrated cross-model exposure or a unique causal explanation. [Twelve exact source/readout rows](supporting/data/native_channel_first10.csv).

The mechanistic connection is conditional: larger incumbents provide different native active distributions to the same type of inactive-to-active recruitment challenge. Stimulation adds stimulus × inactive locally, while each model's feedback, deactivation and transport continue to couple incumbent and challenger. The measured time needed to deplete a finite part of the original region therefore tests the coupled ability to weaken that incumbent; it is not an identified molecular-damage rate or a unique diagnosis of chemistry. The native withdrawal/final profiles show whether that weakening is followed by return or retention. No causal percentage is assigned to incumbent amount, height, width or substrate availability.

All 72 separately evaluated native field-tail tests pass. Complete final active/inactive fields also match either the original basal or its reflection to 5.283 × 10⁻⁸ of the inherited fixed species scales. This supports the finite-window return/reversal descriptions directly from fields, separately from the scalar gate. Native N = 201 resolution and finite 10000-unit withdrawal follow-up limit the inference. No new integrations are required for these comparisons.
''')
    return [foundation_path, sensitivity_path, channel_path, OUT / 'RESULTS.md'], prior_table


def foundation():
    fig = plt.figure(figsize=(15, 9.4))
    fig.text(.06, .945, 'Larger incumbents: the prior all-five foundation', fontsize=19, weight='bold')
    fig.text(.06, .902, 'Exact previously retained amount, width and peak comparison; unchanged artwork.', fontsize=12)
    ax = fig.add_axes([.025, .34, .95, .51])
    ax.imshow(plt.imread(PRIOR))
    ax.axis('off')
    fig.text(.065, .262, 'NEW: 72 full-field challenge protocols in WP and Goryachev', fontsize=17, weight='bold')
    fig.text(.065, .205, 'Exact polarized basal  →  independently vary native peak, physical support and pulse duration  →  observe during cue and release', fontsize=12)
    fig.text(.065, .152, 'At both sizes within each model: peak multipliers 0.5 / 1 / 2; physical widths 3 / 5; three common pulse durations.', fontsize=11)
    fig.text(.065, .112, 'Readouts: finite original-region activity loss, challenger growth, and retained field response after 10000 cue-free native units.', fontsize=11)
    fig.text(.065, .058, 'This is a two-model experimental stage, not all-five challenge coverage. Native stimulus × inactive is not a calibrated optogenetic actuator.', fontsize=9)
    fig.text(.065, .033, 'The new experiments illustrate consequences of the incumbent organization; they do not uniquely identify molecular mechanisms or allocate causal percentages.', fontsize=9)
    return fig


def model_page(model, rows, events):
    sizes = [9., 19.] if model == 'wave_pinning' else [9., 25.]
    durations = [16., 60., 240.] if model == 'wave_pinning' else [74., 400., 1600.]
    label = 'Wave-pinning' if model == 'wave_pinning' else 'Goryachev'
    fig = plt.figure(figsize=(15, 9.4))
    fig.text(.06, .945, f'{label}: performed strength × width × pulse experiment', fontsize=18, weight='bold')
    fig.text(.06, .905, f'36 protocols; exact shared basals at L{sizes[0]:g}/L{sizes[1]:g}; native requested peak = {0.1 if model == "wave_pinning" else 6.4:g} × multiplier.', fontsize=11)
    ax = fig.add_axes([.08, .535, .37, .285])
    combinations = [(width, multiplier) for width in (3., 5.) for multiplier in (.5, 1., 2.)]
    matrix = np.empty((6, 6))
    for i, (width, multiplier) in enumerate(combinations):
        for j, (size, duration) in enumerate((s, t) for s in sizes for t in durations):
            r = get(rows, model, size, width, multiplier, duration)
            matrix[i, j] = float(r['off_incumbent_fraction_last'])
    im = ax.imshow(np.clip(matrix, 0, 1), cmap='viridis', vmin=0, vmax=1, aspect='auto')
    ax.set_yticks(range(6), [f'w{w:g}, {a:g}×' for w, a in combinations], fontsize=9)
    ax.set_xticks(range(6), [f'{t:g}' for _ in sizes for t in durations], fontsize=9)
    ax.set_xlabel('Pulse duration (native)', fontsize=9)
    ax.text(.25, 1.025, f'L = {sizes[0]:g}', transform=ax.transAxes, ha='center', fontsize=10)
    ax.text(.75, 1.025, f'L = {sizes[1]:g}', transform=ax.transAxes, ha='center', fontsize=10)
    ax.axvline(2.5, color='white', lw=2)
    ax.axhline(2.5, color='white', lw=1)
    for i, (width, multiplier) in enumerate(combinations):
        for j, (size, duration) in enumerate((s, t) for s in sizes for t in durations):
            r = get(rows, model, size, width, multiplier, duration)
            if r['full_polarity_success'] == 'True':
                ax.plot(j, i, 'o', mfc='none', mec='white', ms=9, mew=1.7)
    ax.set_title('A  Original-core activity after 10000 off units', fontsize=11, loc='left', weight='bold', pad=27)
    cb = fig.colorbar(im, ax=ax, fraction=.025, pad=.025)
    cb.ax.tick_params(labelsize=8)
    cb.set_label('Active amount / initial original-core amount', fontsize=8.5)
    fig.text(.08, .475, '○ White ring: original retained scalar-strength gate passes. Color: measured regional activity.', fontsize=8.5)

    ax = fig.add_axes([.63, .535, .30, .285])
    colors = {sizes[0]: '#2171b5', sizes[1]: '#983ca0'}
    for size in sizes:
        for width in (3., 5.):
            lo, hi = [], []
            for multiplier in (.5, 1., 2.):
                r = get(rows, model, size, width, multiplier, durations[-1])
                event = next(e for e in events if e['slot'] == r['slot'] and e['phase'] == 'cue_on'
                             and e['regional_depletion_fraction'] == '0.25')
                assert event['status'] == 'sampled_transition'
                lo.append(float(event['lower']))
                hi.append(float(event['upper']))
            ax.errorbar([.5, 1., 2.], hi, yerr=[np.array(hi) - lo, np.zeros(3)],
                        color=colors[size], ls='-' if width == 3 else '--', marker='o', ms=4,
                        label=f'L{size:g}, w{width:g}', capsize=2)
    ax.set_xscale('log', base=2)
    ax.set_xticks([.5, 1., 2.], ['0.5', '1', '2'])
    ax.set_yscale('log')
    ax.set_xlabel('Native requested-peak multiplier', fontsize=9)
    ax.set_ylabel('Native time to 25% original-region activity loss', fontsize=9)
    ax.set_title('B  Finite cue-on weakening time', fontsize=11, loc='left', weight='bold')
    ax.tick_params(labelsize=8)
    ax.grid(alpha=.15)
    ax.legend(frameon=False, fontsize=8)
    fig.text(.63, .475, 'Bars: adjacent saved-time brackets (0.1 unit), not confidence intervals.', fontsize=8.5)

    ax = fig.add_axes([.08, .155, .37, .25])
    size, duration = sizes[1], durations[1]
    peakcolors = {.5: '#888888', 1.: '#983ca0', 2.: '#2171b5'}
    for width in (3., 5.):
        for multiplier in (.5, 1., 2.):
            r = get(rows, model, size, width, multiplier, duration)
            with np.load(ROOT / r['actual_cue_path']) as z:
                manifest = json.loads((ROOT / 'data/runs/revision_width_strength_grid_20261019/manifest.json').read_text())
                c = manifest['cases'][int(r['slot'])]
                with np.load(ROOT / 'data/runs/revision_width_strength_grid_20261019' / c['regional_input']) as zz:
                    core, initial = zz['incumbent_core_mask'], zz['initial_state']
                fraction = z['field_active'][:, core].sum(axis=1) / initial[0, core].sum()
                ax.plot(z['time'], fraction, color=peakcolors[multiplier], ls='-' if width == 3 else '--', lw=1.5)
    ax.set_xlim(0, duration)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel('Time since cue onset (native)', fontsize=9)
    ax.set_ylabel('Original-core active amount / initial', fontsize=9)
    ax.set_title(f'C  L{size:g}, common pulse {duration:g}: actual weakening', fontsize=11, loc='left', weight='bold')
    legend = [Line2D([], [], color=peakcolors[a], label=f'{a:g}× peak') for a in (.5, 1., 2.)]
    legend += [Line2D([], [], color='black', ls=ls, label=f'w={w:g}') for w, ls in ((3, '-'), (5, '--'))]
    ax.legend(handles=legend, frameon=False, fontsize=8, ncol=2)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=.15)

    ax = fig.add_axes([.63, .155, .30, .25])
    initial_drawn = False
    for multiplier in (1., 2.):
        r = get(rows, model, size, 3., multiplier, duration)
        with np.load(ROOT / r['actual_cue_path']) as z:
            x, u = z['x'], z['field_active'][-1]
            ax.plot(x, u, color=peakcolors[multiplier], lw=1.6, label=f'{multiplier:g}×, withdrawal')
        with np.load(ROOT / r['actual_recovery_path']) as z:
            ax.plot(x, z['field_active'][-1], color=peakcolors[multiplier], ls=':', lw=2,
                    label=f'{multiplier:g}×, +10000 off')
        if not initial_drawn:
            c = manifest['cases'][int(r['slot'])]
            with np.load(ROOT / 'data/runs/revision_width_strength_grid_20261019' / c['regional_input']) as z:
                ax.plot(x, z['initial_state'][0], color='.5', lw=1, label='initial')
            initial_drawn = True
    ax.set_xlim(0, size)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Native interval position', fontsize=9)
    ax.set_ylabel('Native active concentration', fontsize=9)
    ax.set_title(f'D  L{size:g}, w3: actual profiles', fontsize=11, loc='left', weight='bold')
    ax.legend(frameon=False, fontsize=7.5, loc='upper center')
    ax.tick_params(labelsize=8)

    if model == 'wave_pinning':
        statement = 'At pulse 60, doubling peak or widening to 5 overcomes the tested L19 scalar resistance; L9 passes every tested combination.'
    else:
        statement = 'At pulse 400, all L9 conditions reverse and all L25 conditions return; pulse 1600 reverses all L25 conditions. Sparse pulses can miss shifted boundaries.'
    fig.text(.065, .083, statement, fontsize=9.2, weight='bold')
    fig.text(.065, .048, 'All 36 native field-tail tests pass. Regional loss ≠ identified molecular damage; scalar strength ≠ complete field replacement. N = 201; no spatial-convergence claim.', fontsize=8.5)
    fig.text(.065, .025, 'Native stimulus × inactive and clocks are model-specific. Width changes integrated input. The 25% readout is descriptive; 10%/50% results remain in the source table.', fontsize=8.5)
    return fig


def main():
    rows = readcsv(DATA / 'finite_assay_readouts.csv')
    events = readcsv(DATA / 'regional_sampled_entries.csv')
    assert len(rows) == 72 and all(r['state'] == 'complete' for r in rows)
    assert all(r['field_tail_pass'] == 'True' for r in rows)
    template_path, template_error = endpoint_templates(rows)
    note_paths, prior_table = result_note(rows, events)
    paths = []
    pdf = OUT / 'supporting/vectors/figure2_performed_cue_grid_review.pdf'
    with PdfPages(pdf) as book:
        fig = foundation()
        path = OUT / 'figure2_prior_organization_foundation.png'
        fig.savefig(path, dpi=140)
        paths.append(path)
        book.savefig(fig)
        plt.close(fig)
        for model in ('wave_pinning', 'goryachev'):
            fig = model_page(model, rows, events)
            png = OUT / f'figure2_{model}_performed_grid.png'
            fig.savefig(png, dpi=150)
            fig.savefig(OUT / 'supporting/vectors' / f'figure2_{model}_performed_grid.svg')
            book.savefig(fig)
            paths.append(png)
            plt.close(fig)
    (OUT / 'REVIEW.md').write_text('# Performed cue experiments: fresh-read guide\n\n'
        '[Three-page review PDF](supporting/vectors/figure2_performed_cue_grid_review.pdf)\n\n'
        'Page 1 preserves the exact prior all-five organization artwork as the foundation. Pages 2–3 report the new '
        '72-protocol two-model experiment: actual retained original-region activity, finite weakening times and native fields. '
        'These are reviewable results, not an accepted final Figure2 composition.\n\n'
        + '\n\n'.join(f'![{p.stem}]({p.name})' for p in paths) + '\n\n'
        'A uses measured activity in a fixed initial spatial core, divided by that core’s initial amount. Rings identify '
        'the original retained scalar-strength gate separately. B uses the longest common pulse in each model to '
        'observe the first sampled 25% regional-depletion transition; bars are adjacent 0.1-unit saved brackets. '
        'C holds pulse fixed at the intermediate duration and shows all six peak/width combinations in the larger cohort. '
        'D displays actual native profiles for the baseline and doubled peak at fixed width 3.\n\n'
        'A separate direct full-field check compares each final active/inactive state with its exact basal and '
        f'with the spatial reflection of that basal. All 72 are within {template_error:.4g} of one template '
        'when each native species is divided by its inherited fixed scale, and the closer template agrees with '
        'the scalar outcome. This is a direct native-field comparison at 10000 cue-free units, not a deduction from '
        'the scalar gate or proof of asymptotic replacement. [All 72 template comparisons](supporting/data/native_endpoint_template_matches.csv).\n\n'
        '[Result-specific interpretation and robustness](RESULTS.md) · [Caption and measurement definitions](CAPTION.md) · [Complete native readouts and limits](README.md) · '
        '[Eight initial full-field width contrasts](../width_fullfield_experiment_20261019/README.md)\n')
    (OUT / 'CAPTION.md').write_text('''# Performed cue-strength, physical-width and pulse-duration experiments

The prior all-five incumbent amount/width/peak comparison is reproduced from unchanged artwork and exact retained measurements. The new experiment is explicitly a two-model stage: WP at L9/L19 and Goryachev at L9/L25. At each size, physical widths 3/5, native peak multipliers 0.5/1/2, and three pulse durations are independently crossed. Baseline requested peaks are 0.1 (WP) and 6.4 (Goryachev). Both sizes receive the same pulse grid within each model: 16/60/240 for WP, 74/400/1600 for Goryachev. The 72 protocols comprise 69 new native integrations and three exact full-field reuses from the preceding eight-pulse experiment. All start from exact copied corrected polarized basals and continue for 10000 native cue-free units after withdrawal.

**A. Retained original-region activity.** Color gives final native active amount in the fixed initial incumbent core, divided by that core's initial amount. The core is the connected region above the initial background-plus-half-range threshold containing the initial activity maximum; integration uses the native cell-centre quadrature dx = L/201. This is a fixed spatial ROI, not a tracked site identity. White rings mark the original scalar-strength criterion separately: reversed signed polarity with magnitude within 1% of its initial magnitude, together with the retained orientation, contrast and persistence gates over the last 500 cue-free units. All 72 native field-tail tests pass independently. Direct final active/inactive comparisons also reproduce either the exact initial basal or its reflection within 5.283e-8 of inherited species fixed scales; this finite-window field comparison is recorded separately from the scalar gate.

**B. Finite weakening time.** The ordinate is the first saved cue-on transition to a 25% decrease of active amount in the original core, measured during the longest common pulse (240 WP; 1600 Goryachev). Error bars enclose adjacent 0.1-unit samples, not confidence intervals or exact continuous-time first crossings. The 10% and 50% thresholds are retained in the complete event table. Regional depletion is not identified molecular damage, and this cue-on latency is not a minimum pulse duration.

**C. Actual regional response at a common pulse.** All six peak/width combinations are shown at the larger size, with pulse 60 (WP) or 400 (Goryachev). Both sizes were tested under those same pulses, as displayed in A. Each curve is native fixed-core active amount divided by its own initial amount; no response is divided by its observed maximum.

**D. Native profiles.** Initial, withdrawal and final activity profiles compare baseline and doubled requested peak at fixed physical width 3 and the common intermediate pulse. The vertical axes show actual native active concentration on a common zero-based scale within each panel.

The native stimulation channel in both models is stimulus × inactive concentration, transferred locally from inactive to active. Different intrinsic activation/deactivation laws and the coupled conserved pool remain active. These characterized native challenges are not matched optical actuators, reporter calibrations or unique molecular identifiers. Width changes integrated input as well as footprint. Three pulse durations can miss shifted success boundaries: the absence of a Goryachev strength-dependent fate change in this grid does not establish strength irrelevance. Native N=201 discretization and fixed-scale conventions are inherited; no independent spatial-convergence claim is made.

New trajectories save both native species every 0.1 cue-on unit and every 0.25 cue-free unit through 100, then every 2.5 units through 10000. The three reused archives retain 0.25-unit sampling throughout release. Event tables preserve adjacent-sample brackets and these cadence differences. All canonical 2.5-unit endpoint/gate samples are present. The original solver, parameters and two-phase exact withdrawal restart are unchanged; only output sampling differs. No new integration is performed by this review builder.
''')
    sources = {str(p.relative_to(ROOT)): sha(p) for p in
               (Path(__file__), PRIOR, prior_table, DATA / 'finite_assay_readouts.csv', DATA / 'regional_sampled_entries.csv')}
    artifacts = {str(p.relative_to(OUT)): sha(p) for p in paths + [pdf, OUT / 'REVIEW.md', OUT / 'CAPTION.md', template_path] + note_paths
                 + list((OUT / 'supporting/vectors').glob('*performed_grid.svg'))}
    (OUT / 'supporting/metadata/review_validation.json').write_text(json.dumps(dict(
        protocols=72, models=2, field_tail_passes=72, source_hashes=sources,
        artifact_hashes=artifacts), indent=2) + '\n')
    print(json.dumps(dict(pdf=str(pdf), pages=3, pngs=[str(p) for p in paths]), indent=2))


if __name__ == '__main__':
    main()
