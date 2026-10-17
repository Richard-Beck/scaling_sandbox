"""Export representative trajectories, diagnose the gate, and write the report."""
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
    parser.add_argument('--out', type=Path, default=HERE/'results')
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

    longest = max(heldout, key=lambda c: c['pulses'][-1, 1]-c['pulses'][-1, 0])
    y = integrate(theta, longest['times'], longest['pulses'], np.zeros(2))
    gate = 1/(1+np.exp(np.clip(theta[7]-theta[6]*y[:, 0]*y[:, 1], -700, 700)))
    onset = longest['pulses'][-1, 0]
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, layout='constrained')
    for j, label in enumerate(['p = Rac_x50/L - 0.5', 'v = 0.5 - MCA_x50/L']):
        axes[j].plot(longest['times']-onset, longest['values'][:, j], color='black', label='De Belly')
        axes[j].plot(longest['times']-onset, y[:, j], color='#d55e00', label='Reduced model')
        axes[j].set_ylabel(label)
        axes[j].axhline(0, color='gray', lw=.5)
        axes[j].legend(fontsize=8)
    axes[2].plot(longest['times']-onset, gate, color='#0072b2')
    axes[2].set(ylabel='Reduced gate S(p*v)', xlabel='Time from countercue onset (s)', ylim=(-.03, 1.03), xlim=(-300, 2000))
    for ax in axes:
        ax.axvspan(0, 600, color='#e8be54', alpha=.2)
        ax.grid(alpha=.15)
    fig.suptitle('Why the fitted model resists reversal: its shared gate closes during disagreement')
    fig.savefig(out/'gate_diagnostic.png', dpi=170)
    plt.close(fig)

    provenance = json.loads((out/'provenance.json').read_text())
    mapping = pd.read_csv(out/'parameter_mapping.csv')
    thresholds = pd.read_csv(out/'reversal_thresholds.csv')
    reference_thresholds = thresholds[(thresholds.model == 'frozen_reference') & (thresholds.length_um == 29) & (thresholds.variant == 'beta_1x')]
    threshold_text = ', '.join(f'{r.field}: {r.reference_min_s:.2f} s' for r in reference_thresholds.itertuples())
    parameter_rows = [[name, f'{best["parameters"][name]:.6g}'] for name in ['tau_p','tau_v','alpha','beta','ps','vs','k','z0','cue_gain']]
    local_rows = [[f'{r.length_um:g}', r.variant, f'{r.tau_p:.2f}', f'{r.tau_v:.2f}', f'{r.cost:.6g}'] for r in mapping.itertuples()]
    near = [r for r in candidates if r['cost'] <= 1.01*best['cost']]
    report = f'''# De Belly matching pilot

Generated {datetime.now(timezone.utc).isoformat()}.

The two-variable model closely reproduces polarity establishment and persistence
at the reference length, but the best establishment-only fit fails the held-out
reversal test. This is a useful local parameter mapping, not a validated surrogate
for all De Belly dynamics.

The exact retained Figure 2 half-mass observables were used:

    p = Rac_x50/L - 0.5
    v = 0.5 - bound_MCA_x50/L

No additional amplitude normalization, no experimental velocity calibration, and
no stochastic noise were introduced. Time is in seconds. The original signed cue
envelopes enter the p equation through one calibrated gain. All original De Belly
trajectories and campaign manifests remain unchanged.

## Reference calibration at 29 um

Only panel A/full and panel C/beta_1x were fitted. {len(candidates)} starts were tried;
{len(near)} reached costs within 1% of the best. Best weighted least-squares cost:
{best['cost']:.8g}. Different higher-cost local optima were also found. The
parameters are effective coefficients, not unique mechanistic measurements;
this pilot provides no statistical confidence intervals.

{table(['Parameter', 'Value'], parameter_rows)}

tau_p and tau_v have units seconds; ps and vs use centered relative half-mass
coordinates; alpha and beta convert between these states per second; k has
inverse-product-state units and z0 has product-state units. Large k reflects
the small numerical scale of the centered polarity product. cue_gain converts
the projected activation input to Rac-polarity change per second.

The asymptotic polarity amplitudes are reproduced closely. The full sampled
trajectory errors below are in x50/L units; the optimizer's phase-balanced loss
is different from these unweighted RMS errors.

{table(['Assay', 'Rac RMSE', 'MCA RMSE', 'Rac final error', 'MCA final error'], summary_rows)}

![Calibration](calibration.png)

## Held-out countercues

No countercue observations were used for fitting or candidate selection. The
150-s and 300-s cues return to the original direction in both models, although
the reduced model underestimates MCA redistribution and recovery delay. At 600 s,
De Belly reverses both fields, while the reduced model returns to the original
direction. In the 50-duration sweep at 29 um, the retained minimum successful
durations are {threshold_text}. The fitted reduced model has no successful
duration from 10 to 1000 s, even under its self-conditioned reversal criterion.
These grid minima are not continuous threshold estimates.

![Held-out trajectories](heldout_countercues.png)

The fitted gate is almost shut when the countercue flips Rac while MCA still
points in the old direction. The alpha-driven term that could redistribute MCA
then becomes very small, leaving mainly decay on the roughly {theta[1]:.0f}-s MCA memory
timescale. This explains the observed failure of this fitted parameter set. It
does not prove that every possible parameterization of the model must fail.

![Gate diagnostic](gate_diagnostic.png)

## Length and feedback transfer

An unchanged two-state model has no explicit length input and therefore cannot
recover De Belly's length-dependent polarization transition with otherwise
identical projected cues. Exploratory condition-specific refits reproduce the
final polarity more closely, including the unpolarized 17-um condition. Each
used only establishment data and one start from the reference fit. No simple
length-scaling law is inferred from these independent fits.

{table(['Length (um)', 'De Belly variant', 'tau_p (s)', 'tau_v (s)', 'Calibration cost'], local_rows)}

The 0x and 2x labels modify De Belly's rho_mechanical_feedback, not the reduced
model's beta by definition. The 2x conditional fit has only panel C available.
All these conditional fits still fail to reproduce successful tested reversals
in the polarized cases. Changing parameters by condition improves establishment
but has not solved the countercue mismatch.

![Length transfer](length_transfer.png)
![Reversal thresholds](reversal_thresholds.png)

## Numerical checks and interpretation

Four unit checks passed: uncoupled decay, cue projection, exact handling of cue
events/reflection, and the reversal/confirmation criteria. Across {len(step_rows)}
reference and conditional trajectory checks, halving the integration step from
0.5 to 0.25 s changes either coordinate by at most
{check['maximum_trajectory_difference']:.3g}. None of the 100 field-level
reference-length duration-sweep labels changed. Numerical error is much smaller
than the figure's 0.002 settling tolerance.

The campaign evaluated all 3,000 original sweep conditions with frozen reference
parameters, plus condition-specific predictions for eight length/feedback pairs.
The core campaign completed in {provenance['wall_seconds']:.1f} s
using one CPU. Numerical error does not explain the reversal mismatch.

The present parameters are suitable as an approximation for establishment and
withdrawal at 29 um. They should not be used as a validated reversal mapping.
A next model comparison could relax the shared gate on the p-to-MCA pathway or
include a cue-dependent MCA-driving contribution. Neither modification was
made in this campaign. An alternative is to fit reversal data too, explicitly
treating that as a new calibration and retaining other countercues for validation.

## Files

- `best_parameters.json`: reference parameters, both noise amplitudes zero.
- `parameter_mapping.csv`: all eight conditional parameter sets.
- `fit_*.json`: optimizer candidates, convergence, bounds, Jacobian singular values.
- `representative_trajectories.npz`: time and paired reference/predicted x50/L
  arrays for the two calibration and three held-out reference-length assays;
  coordinate order is Rac, bound MCA.
- `trajectory_metrics.csv`: full sampled trajectory and final-state errors.
- `reversal_sweep.csv`: every duration prediction and reference label, including
  a self-conditioned diagnostic and conditioning errors.
- `reversal_thresholds.csv`: grid minimum successful durations; empty means no
  success among tested durations, including unpolarized/ineligible cases.
- `provenance.json`, `step_check.json`, `extended_numerical_checks.json`: source
  hashes, environments, calibration IDs, and numerical validation.

See [campaign instructions](../README.md) for complete calibration and criterion
definitions. Rebuild this report with `python summarize.py` in the existing
scaling_sandbox Conda environment after completing the campaign.
'''
    (out/'REPORT.md').write_text(report)
    (out/'report_provenance.json').write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    print(f'Report: {out / "REPORT.md"}')
    print(f'Max step-halving error: {check["maximum_trajectory_difference"]:.3g}; changed sweep labels: {len(label_changes)}')


if __name__ == '__main__':
    main()
