"""Native amplitude-extension readouts with explicit finite censoring and channel quadrature checks."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'data/runs/revision_amplitude_operating_point_20261019'
OUT = ROOT / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/amplitude_operating_point_20261019'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_suffix('.tmp.json')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(p)


def csvwrite(p, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    temp = p.with_suffix('.tmp.csv')
    with temp.open('w') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    temp.replace(p)


def entries(time, condition):
    hits = np.flatnonzero(condition)
    if not len(hits):
        return dict(status='right_censored_at_observed_horizon', lower=float(time[-1]), upper='', retained_status='right_censored_at_observed_horizon',
                    retained_lower=float(time[-1]), retained_upper='')
    i = int(hits[0])
    failing = np.flatnonzero(~condition)
    retained = int(failing[-1] + 1) if len(failing) else 0
    result = dict(status='already_present' if i == 0 else 'sampled_transition',
                  lower=float(time[max(0, i - 1)]), upper=float(time[i]))
    if retained < len(time):
        result.update(retained_status='present_all_saved_frames' if retained == 0 else 'sampled_retained_entry',
                      retained_lower=float(time[max(0, retained - 1)]), retained_upper=float(time[retained]))
    else:
        result.update(retained_status='not_retained_through_last_frame', retained_lower='', retained_upper='')
    return result


def main():
    for directory in ('supporting/data', 'supporting/metadata', 'supporting/vectors', 'supporting/viewers'):
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    manifest = read(RUN / 'manifest.json')
    source_hashes = {str((RUN / 'manifest.json').relative_to(ROOT)): sha(RUN / 'manifest.json')}
    summary, events, profiles, channels = [], [], [], []
    completed = 0
    for c in manifest['cases']:
        row = dict(slot=c['slot'], model=c['model'], length=c['length'], physical_width=c['width'],
                   native_peak_multiplier=c['native_peak_multiplier'], requested_peak=c['response_amplitude'],
                   sampled_peak=c['sampled_peak'], sampled_spatial_integral=c['sampled_spatial_integral'],
                   pulse_duration=c['cue_duration'], off_duration=10000., state='not_started',
                   source_kind='exact_fullfield_reuse' if c['exact_reuse'] else 'new_native',
                   checkpoint_sha256=c['checkpoint_sha256'])
        folder = RUN / 'workers' / f'{c["slot"]:03d}'
        reuse = c['exact_reuse']
        result_path = ROOT / reuse['result_path'] if reuse else folder / 'result.json'
        if not result_path.exists():
            if (folder / 'error.json').exists():
                row['state'] = 'error'
            elif (folder / 'status.json').exists():
                row['state'] = read(folder / 'status.json')['state']
            summary.append(row)
            continue
        result = read(result_path)
        if reuse:
            assert sha(result_path) == reuse['result_sha256']
        source_hashes[str(result_path.relative_to(ROOT))] = sha(result_path)
        row['state'] = 'complete'
        completed += 1
        row['recovery_cadence'] = '0.25_through100_then2.5'
        row['cue_cadence'] = 'legacy_uniform0.1' if reuse else 'fine0.001_to1_0.005_to10_then0.1_plus_initiallog'
        with np.load(RUN / c['regional_input']) as zz:
            masks = {k: zz[k].copy() for k in zz.files}
        source_hashes[str((RUN / c['regional_input']).relative_to(ROOT))] = sha(RUN / c['regional_input'])
        core, roi = masks['incumbent_core_mask'], masks['cue_footprint_mask']
        initial = masks['initial_state']
        dx = c['length'] / 201
        initial_core = float(initial[0, core].sum() * dx)
        row['initial_core_active_amount'] = initial_core
        row['initial_core_physical_width'] = float(core.sum()*dx)
        row['initial_native_active_peak'] = float(initial[0].max())
        threshold = float(masks['initial_halfcontrast_threshold'])
        for phase, filename in (('cue_on', 'cue.npz'), ('off', 'recovery.npz')):
            p = ROOT / reuse['cue_path' if phase == 'cue_on' else 'recovery_path'] if reuse else folder / filename
            digest = sha(p)
            assert digest == result['cue_sha256' if phase == 'cue_on' else 'recovery_sha256']
            source_hashes[str(p.relative_to(ROOT))] = digest
            row['actual_cue_path' if phase == 'cue_on' else 'actual_recovery_path'] = str(p.relative_to(ROOT))
            row['actual_cue_sha256' if phase == 'cue_on' else 'actual_recovery_sha256'] = digest
            with np.load(p) as z:
                t, u, v = z['time'], z['field_active'], z['field_inactive']
                if phase == 'cue_on':
                    with np.load(RUN / c['regional_input']) as region:
                        signal = region['signal']
                    for window in manifest['channel_quadrature_qualification']['windows']:
                        found = np.flatnonzero(np.isclose(t, window, rtol=0, atol=1e-12))
                        channel_row = dict(slot=c['slot'], model=c['model'], length=c['length'],
                            native_peak_multiplier=c['native_peak_multiplier'], window_end=window,
                            source_kind=row['source_kind'], primary_quadrature_qualified=False)
                        if not len(found):
                            channel_row.update(status='window_endpoint_not_retained', fine_integral='',
                                               decimated_integral='', relative_discrepancy='')
                        else:
                            last = int(found[0]); times = t[:last+1]
                            transfer = np.sum(v[:last+1] * signal, axis=1) * dx
                            fine = float(np.trapezoid(transfer, times))
                            ids = np.unique(np.r_[np.arange(0, last+1, 2), last])
                            coarse = float(np.trapezoid(transfer[ids], times[ids]))
                            discrepancy = abs(fine-coarse)/max(abs(fine), 1e-300)
                            qualified = bool(not reuse and discrepancy <= manifest['channel_quadrature_qualification']['relative_fine_vs_decimated_tolerance'])
                            channel_row.update(status='legacy_sampling_unqualified' if reuse else 'passes_saved_quadrature_check' if qualified else 'fails_saved_quadrature_check',
                                fine_integral=fine, decimated_integral=coarse, relative_discrepancy=discrepancy,
                                primary_quadrature_qualified=qualified, saved_frames=last+1,
                                cue_roi_inactive_end_over_initial=float(v[last, roi].mean()/v[0, roi].mean()))
                        channels.append(channel_row)
                fraction = u[:, core].sum(axis=1) * dx / initial_core
                peak = u[:, roi].max(axis=1)
                row[f'{phase}_incumbent_fraction_last'] = float(fraction[-1])
                row[f'{phase}_minimum_incumbent_fraction'] = float(fraction.min())
                row[f'{phase}_cue_peak_last'] = float(peak[-1])
                row[f'{phase}_cue_active_amount_last'] = float(u[-1, roi].sum() * dx)
                row[f'{phase}_inactive_mean_last'] = float(v[-1].mean())
                row[f'{phase}_saved_frames'] = len(t)
                for loss in (.1, .25, .5):
                    events.append(dict(slot=c['slot'], model=c['model'], length=c['length'], phase=phase,
                        observable='fixed_original_core_active_fraction', condition='less_than_or_equal',
                        threshold=1 - loss, regional_depletion_fraction=loss,
                        time_units='native_phase_local', **entries(t, fraction <= 1 - loss)))
                events.append(dict(slot=c['slot'], model=c['model'], length=c['length'], phase=phase,
                    observable='actual_cue_footprint_native_active_peak', condition='greater_than_or_equal',
                    threshold=threshold, regional_depletion_fraction='',
                    time_units='native_phase_local', **entries(t, peak >= threshold)))
                if phase == 'off':
                    above = u[-1] > threshold
                    row['final_regions_above_initial_halfcontrast'] = int(np.count_nonzero(
                        np.diff(np.r_[False, above, False].astype(int)) == 1))
                    final = np.vstack((u[-1], v[-1]))
                    scales = np.array([c['fixed_scales'][k] for k in ('active', 'inactive')])[:, None]
                    row['final_initial_template_max_scaled_error'] = float(np.max(abs(final-initial)/scales))
                    row['final_reflected_template_max_scaled_error'] = float(np.max(abs(final-initial[:, ::-1])/scales))
                for i, x in enumerate(masks['x']):
                    profiles.append(dict(slot=c['slot'], phase='withdrawal' if phase == 'cue_on' else 'final_off10000',
                        x=float(x), active=float(u[-1, i]), inactive=float(v[-1, i]),
                        initial_active=float(initial[0, i]), initial_inactive=float(initial[1, i]),
                        incumbent_core=bool(core[i]), cue_footprint=bool(roi[i])))
        row.update({k: result['legacy_diagnostics'][k] for k in
                    ('orientation_success', 'full_polarity_success', 'field_tail_pass', 'final_orientation')})
        row['maximum_total_amount_relative_error'] = result['maximum_total_amount_relative_error']
        summary.append(row)
    csvwrite(OUT / 'supporting/data/finite_assay_readouts.csv', summary)
    csvwrite(OUT / 'supporting/data/regional_sampled_entries.csv', events)
    csvwrite(OUT / 'supporting/data/withdrawal_final_native_profiles.csv', profiles)
    csvwrite(OUT / 'supporting/data/realized_channel_quadrature.csv', channels)
    now = datetime.now(timezone.utc).isoformat()
    (OUT / 'README.md').write_text('# Native amplitude operating-point extension: actual readouts\n\n'
        f'Generated {now}. Complete scientific protocols: {completed}/32, including twelve exact full-field reuses. '
        'This is a two-model stage (WP/Goryachev), not five-model functional coverage.\n\n'
        '[All 32 protocol/readout rows](supporting/data/finite_assay_readouts.csv) · '
        '[Regional sampled-entry brackets](supporting/data/regional_sampled_entries.csv) · '
        '[Native withdrawal/final profiles](supporting/data/withdrawal_final_native_profiles.csv) · '
        '[Saved-output channel quadrature checks](supporting/data/realized_channel_quadrature.csv) · '
        '[Source hashes/status](supporting/metadata/latest_status.json)\n\n'
        'At fixed physical width 3, peak multipliers 0.01 / 0.03 / 0.1 / 0.3 / 3 extend the exact retained 0.5 / 1 / 2 conditions. '
        'Both sizes receive WP 240 / Gory 1600 pulses. The unchanged native input is '
        'stimulus × inactive. Exact corrected '
        'basals and N = 201 native grids are retained; no optical-actuator or spatial-convergence claim is made.\n\n'
        'Core activity uses the fixed initial connected above-halfcontrast region. Activity loss is regional depletion, '
        'not identified molecular damage. Challenger peaks use the actual cue footprint. Original scalar-strength '
        'and field-tail flags remain separate from regional field observations. Final whole-field deviations '
        'from the original and reflected basal templates are diagnostics on the inherited fixed species scales; '
        'they do not replace either gate.\n\n'
        'Entry brackets enclose adjacent saved samples at the first observed false→true transition; they do not '
        'exclude an unobserved earlier excursion between samples. Retained entry means all subsequent saved frames '
        'satisfy the condition. Cue-on timing is not pulse duration. Compare brackets, not unqualified '
        'point latencies. All canonical 2.5-unit gate samples are present. All twelve reused histories use the same '
        '0.25-through-100 / 2.5-afterwards off cadence as the new histories. The cue sampling differs: reused cases retain uniform 0.1 units, '
        'while new cases retain dense early times. Fine-versus-decimated integral comparisons quantify saved-output '
        'quadrature sensitivity on new cases, not solver convergence. Legacy channel integrals remain unqualified.\n')
    artifacts = {str(p.relative_to(OUT)): sha(p) for p in (OUT / 'supporting/data').glob('*.csv')}
    source_hashes[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    write(OUT / 'supporting/metadata/latest_status.json', dict(generated_at=now, complete=completed == 32,
        completed=completed, total=32, producer_status=read(RUN / 'results/latest_status.json'),
        source_hashes=source_hashes, artifact_hashes=artifacts))
    print(json.dumps(dict(completed=completed, events=len(events), profile_rows=len(profiles), channel_records=len(channels), output=str(OUT)), indent=2))


if __name__ == '__main__':
    main()
