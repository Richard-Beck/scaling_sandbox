"""Export representative trajectories, compare ungated results, and write the report."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from campaign import (HERE, NAMES, integrate, load_data, pulses_from_task,
                      reversal, training_cases)


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |',
                      '| '+' | '.join(['---']*len(headers))+' |'] +
                     ['| '+' | '.join(map(str, r))+' |' for r in rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE/'results_ungated')
    args = parser.parse_args()
    out = args.out
    cases, manifests, _ = load_data()
    candidates = json.loads((out/'fit_reference.json').read_text())
    best = min(candidates, key=lambda r: r['cost'])
    theta = np.array(best['theta'])
    train = training_cases(cases, 29., 'beta_1x')
    heldout = [c for c in cases if c['panel'] == 'E' and c['length'] == 29.]
    arrays = {}
    metrics = pd.read_csv(out/'trajectory_metrics.csv')
    summary_rows = []
    for c in train + heldout:
        y = integrate(theta, c['times'], c['pulses'], np.zeros(2))
        key = c['panel']+'_'+c['task_id']
        arrays[key+'_time_s'] = c['times']
        arrays[key+'_reference_x50_over_L'] = np.column_stack([.5+c['values'][:, 0], .5-c['values'][:, 1]])
        arrays[key+'_prediction_x50_over_L'] = np.column_stack([.5+y[:, 0], .5-y[:, 1]])
        row = metrics[(metrics.task_id == c['task_id']) & (metrics.panel == c['panel']) & (metrics.model == 'frozen_reference')].set_index('field')
        label = 'Calibration '+c['panel'] if c['panel'] != 'E' else f'Held-out {c["pulses"][-1, 1]-c["pulses"][-1, 0]:g}-s countercue'
        summary_rows.append([label, f'{row.loc["rac", "rmse"]:.5f}', f'{row.loc["mca", "rmse"]:.5f}',
                             f'{row.loc["rac", "final_error"]:.5f}', f'{row.loc["mca", "final_error"]:.5f}'])
    np.savez_compressed(out/'representative_trajectories.npz', **arrays)

    # Additional step checks cover all condition fits and all reference-length
    # 1x duration sweep labels, not just the displayed reference-fit curves.
    step_rows, label_changes = [], []
    for path in sorted(out.glob('fit_*.json')):
        best_local = min(json.loads(path.read_text()), key=lambda r: r['cost'])
        localtheta = np.array(best_local['theta'])
        if path.name == 'fit_reference.json':
            length, variant = 29., 'beta_1x'
        else:
            length_text, variant = path.stem.removeprefix('fit_L').split('_', 1)
            length = float(length_text)
        selected = training_cases(cases, length, variant)
        if variant == 'beta_1x':
            selected += [c for c in cases if c['panel'] == 'E' and c['length'] == length]
        for c in selected:
            coarse = integrate(localtheta, c['times'], c['pulses'], np.zeros(2), .5)
            fine = integrate(localtheta, c['times'], c['pulses'], np.zeros(2), .25)
            step_rows.append(dict(fit=path.stem, task_id=c['task_id'], max_abs_difference=float(np.max(np.abs(coarse-fine)))))
    for task in manifests['E'].values():
        if task['assay'] != 'countercue_threshold_sweep' or task['simulation']['length'] != 29 or task['coupling_variant'] != 'beta_1x':
            continue
        end = task['simulation']['end_time']
        t = np.unique(np.r_[np.arange(0., end, 5.), end])
        pulse = pulses_from_task(task)
        a = integrate(theta, t, pulse, np.zeros(2), .5)
        b = integrate(theta, t, pulse, np.zeros(2), .25)
        audit = task['polarity_reversal_audit']
        for j, field in enumerate(['rac', 'mca']):
            ref = (1 if j == 0 else -1)*(audit[f'pre_{field}_reference']-.5)
            if reversal(t, a[:, j], audit['opposite_stop'], ref)[0] != reversal(t, b[:, j], audit['opposite_stop'], ref)[0]:
                label_changes.append(dict(task_id=task['task_id'], field=field))
    check = dict(trajectory_checks=step_rows, reference_sweep_label_changes=label_changes,
                 maximum_trajectory_difference=max(r['max_abs_difference'] for r in step_rows))
    (out/'extended_numerical_checks.json').write_text(json.dumps(check, indent=2)+'\n')
    assert check['maximum_trajectory_difference'] < .0002
    assert not label_changes

    thresholds = pd.read_csv(out/'reversal_thresholds.csv')
    selected = thresholds[(thresholds.model == 'frozen_reference') & (thresholds.length_um == 29) & (thresholds.variant == 'beta_1x')]
    def number(v):
        return 'No success in 10–1000 s' if pd.isna(v) else f'{v:.2f}'
    threshold_rows = [[r.field, number(r.reference_min_s), number(r.predicted_min_s)] for r in selected.itertuples()]
    parameter_rows = [[name, f'{best["parameters"][name]:.6g}'] for name in NAMES]
    previous = pd.read_csv(HERE/'results/trajectory_metrics.csv')
    comparison = []
    for label, frame in [('Gated', previous), ('Ungated', metrics)]:
        for panel in ['A', 'C', 'E']:
            subset = frame[(frame.model == 'frozen_reference') & (frame.length_um == 29) & (frame.panel == panel)]
            subset = subset[subset.variant == ('full' if panel == 'A' else 'beta_1x')]
            comparison.append([label, panel, f'{subset[subset.field == "rac"].rmse.mean():.5f}',
                               f'{subset[subset.field == "mca"].rmse.mean():.5f}'])
    near = [r for r in candidates if r['cost'] <= 1.01*best['cost']]
    provenance = json.loads((out/'provenance.json').read_text())
    report = f'''# Ungated De Belly matching

Generated {datetime.now(timezone.utc).isoformat()}.

The logistic gate has been removed entirely from both feedback terms. Seven
parameters were recalibrated; there are no k or z0 parameters. The previous
gated outputs remain in `../results/` for comparison.

For this run, removing the gate restores the held-out 600-s reversal and gives
thresholds near De Belly's reference value, while making the establishment fit
less precise. The fitted parameters remain an approximation: transient cue
amplitudes and MCA relaxation are not matched exactly.

    dp/dt = -p/tau_p + beta*v/(1+abs(v)/vs) + cue_gain*u(t)
    dv/dt = -v/tau_v + alpha*p/(1+abs(p)/ps)
    p = Rac_x50/L - 0.5
    v = 0.5 - bound_MCA_x50/L

Time is seconds, states are centered relative half-mass coordinates, and both
noise strengths are zero. The signed projected cue, initial conditions, training
data, optimizer bounds on remaining parameters, phase weighting, and validation
criteria are unchanged. Only A/full and C/1x at 29 um enter the reference fit;
all countercues and reversal labels are held out.

## Reference parameters

{table(['Parameter', 'Value'], parameter_rows)}

{len(near)} of {len(candidates)} starts reached within 1% of the best cost
({best['cost']:.7g}). These are effective parameters, without confidence intervals.
tau_p and tau_v are seconds. ps and vs use the centered polarity units;
alpha and beta convert state changes per second. The cue gain converts the
spatially projected activation input to p change per second.

## Calibration and held-out trajectories

{table(['Assay', 'Rac RMSE', 'MCA RMSE', 'Rac final error', 'MCA final error'], summary_rows)}

Errors are in x50/L units. The RMS errors use all retained samples and differ
from the phase-balanced optimization objective. For panel E, the following
comparison averages the three separate trajectory RMS errors equally.

{table(['Model', 'Panel', 'Rac mean RMSE', 'MCA mean RMSE'], comparison)}

![Calibration](calibration.png)
![Held-out countercues](heldout_countercues.png)

## Held-out threshold at 29 um, 1x

{table(['Field', 'De Belly minimum cue (s)', 'Ungated minimum cue (s)'], threshold_rows)}

These are grid minima over the original 50 durations, not continuous threshold
estimates. The exact retained figure criteria are used: 0.002 settling tolerance,
1,000-s confirmation, 0.05 reflected-state tolerance, and at least 0.02 polarity
before and after reversal. Primary labels use the retained De Belly conditioning
reference; self-conditioned diagnostic labels are also saved.

![All reversal thresholds](reversal_thresholds.png)

## Transfer and numerical validation

All 3,000 original duration-sweep conditions were evaluated with frozen reference
parameters; eight length/feedback pairs also have conditional fits or the reference
fit. Condition-specific fits use establishment data only, with a single start.
They are exploratory mappings, not a fitted length-scaling law. The ungated model
has no explicit cell-length input; a common parameter set cannot produce length
dependence when projected cues are identical.

![Length transfer](length_transfer.png)

Step-halving across {len(step_rows)} calibration/validation trajectories changes
the coordinates by at most {check['maximum_trajectory_difference']:.3g}; none of
the 100 reference-length field-level sweep labels changed. The core campaign
ran in {provenance['wall_seconds']:.1f} s. No target coordinates are clipped.

## Saved outputs

`best_parameters.json`, `parameter_mapping.csv`, `fit_*.json`,
`trajectory_metrics.csv`, `reversal_sweep.csv`, `reversal_thresholds.csv`,
`representative_trajectories.npz`, `provenance.json`, and numerical-check JSONs
provide the parameters, predictions, data identities, and checks. The NPZ stores
the two calibration and three held-out 29-um trajectories in original x50/L
coordinates (columns Rac, bound MCA).

The separately regenerated 2D tracks use the user's specified parameters,
not these fitted parameters: tau_p=1, tau_v=10, alpha=beta=vs=ps=1, duration
1000, dt=0.005, sample interval 1, sigma_p=sigma_v=0.15 or 0.45. See
[2D overview](../../results_dimensional_ungated/trajectories_overview.png).
'''
    (out/'REPORT.md').write_text(report)
    (out/'report_provenance.json').write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    print(selected.to_string(index=False))
    print(f'Report: {out / "REPORT.md"}')
    print(f'Max step-halving error: {check["maximum_trajectory_difference"]:.3g}; changed sweep labels: {len(label_changes)}')


if __name__ == '__main__':
    main()
