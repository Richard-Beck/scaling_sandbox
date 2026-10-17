# Ungated De Belly matching

Generated 2026-09-24T18:50:29.602068+00:00.

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

| Parameter | Value |
| --- | --- |
| tau_p | 37.0196 |
| tau_v | 1021.69 |
| alpha | 0.00243249 |
| beta | 0.0239989 |
| ps | 0.161934 |
| vs | 0.188299 |
| cue_gain | 0.0100614 |

8 of 8 starts reached within 1% of the best cost
(0.003000106). These are effective parameters, without confidence intervals.
tau_p and tau_v are seconds. ps and vs use the centered polarity units;
alpha and beta convert state changes per second. The cue gain converts the
spatially projected activation input to p change per second.

## Calibration and held-out trajectories

| Assay | Rac RMSE | MCA RMSE | Rac final error | MCA final error |
| --- | --- | --- | --- | --- |
| Calibration A | 0.00537 | 0.00872 | 0.00097 | 0.00376 |
| Calibration C | 0.00436 | 0.00661 | 0.00135 | 0.00480 |
| Held-out 150-s countercue | 0.00630 | 0.01025 | 0.00119 | 0.00435 |
| Held-out 300-s countercue | 0.00951 | 0.01193 | 0.00082 | 0.00334 |
| Held-out 600-s countercue | 0.01044 | 0.01018 | 0.00114 | 0.00422 |

Errors are in x50/L units. The RMS errors use all retained samples and differ
from the phase-balanced optimization objective. For panel E, the following
comparison averages the three separate trajectory RMS errors equally.

| Model | Panel | Rac mean RMSE | MCA mean RMSE |
| --- | --- | --- | --- |
| Gated | A | 0.00169 | 0.00048 |
| Gated | C | 0.00119 | 0.00195 |
| Gated | E | 0.03955 | 0.06724 |
| Ungated | A | 0.00537 | 0.00872 |
| Ungated | C | 0.00436 | 0.00661 |
| Ungated | E | 0.00875 | 0.01079 |

![Calibration](calibration.png)
![Held-out countercues](heldout_countercues.png)

## Held-out threshold at 29 um, 1x

| Field | De Belly minimum cue (s) | Ungated minimum cue (s) |
| --- | --- | --- |
| rac | 353.47 | 373.67 |
| mca | 353.47 | 393.88 |

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

Step-halving across 33 calibration/validation trajectories changes
the coordinates by at most 8.04e-10; none of
the 100 reference-length field-level sweep labels changed. The core campaign
ran in 29.6 s. No target coordinates are clipped.

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
