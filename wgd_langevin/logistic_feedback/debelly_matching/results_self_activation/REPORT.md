# De Belly matching with polarity self-activation

Generated 2026-09-24T19:04:21.286045+00:00.

The ungated model now includes saturating polarity self-activation with n=2.
Nine parameters were recalibrated, including a>=0 and K>0. The prior ungated
outputs remain in `../results_ungated/`; no countercue data enter this fit.

    dp/dt = -p/tau_p + beta*v/(1+abs(v)/vs) + a*p*abs(p)/(K*K+p*p) + cue_gain*u(t)
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
| tau_p | 23.2632 |
| tau_v | 1060.71 |
| alpha | 0.00288497 |
| beta | 0.0274606 |
| ps | 0.0835263 |
| vs | 0.047697 |
| cue_gain | 0.00927808 |
| a | 0.0032455 |
| K | 0.0572758 |

8 of 8 starts reached within 1% of the best cost
(0.001311623). These are effective parameters, without confidence intervals.
tau_p and tau_v are seconds. ps and vs use the centered polarity units;
alpha and beta convert state changes per second. The cue gain converts the
spatially projected activation input to p change per second.

## Calibration and held-out trajectories

| Assay | Rac RMSE | MCA RMSE | Rac final error | MCA final error |
| --- | --- | --- | --- | --- |
| Calibration A | 0.00508 | 0.00364 | 0.00077 | 0.00153 |
| Calibration C | 0.00129 | 0.00382 | 0.00078 | 0.00161 |
| Held-out 150-s countercue | 0.01619 | 0.00239 | 0.00077 | 0.00158 |
| Held-out 300-s countercue | 0.01831 | 0.00674 | 0.00077 | 0.00152 |
| Held-out 600-s countercue | 0.02184 | 0.00503 | 0.00077 | 0.00154 |

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
| Self-activation | A | 0.00508 | 0.00364 |
| Self-activation | C | 0.00129 | 0.00382 |
| Self-activation | E | 0.01878 | 0.00472 |

![Calibration](calibration.png)
![Held-out countercues](heldout_countercues.png)

## Held-out threshold at 29 um, 1x

| Field | De Belly minimum cue (s) | Self-activation minimum cue (s) |
| --- | --- | --- |
| rac | 353.47 | 434.29 |
| mca | 353.47 | 434.29 |

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
the coordinates by at most 4.83e-08; none of
the 100 reference-length field-level sweep labels changed. The core campaign
ran in 36.2 s. No target coordinates are clipped.

## Saved outputs

`best_parameters.json`, `parameter_mapping.csv`, `fit_*.json`,
`trajectory_metrics.csv`, `reversal_sweep.csv`, `reversal_thresholds.csv`,
`representative_trajectories.npz`, `provenance.json`, and numerical-check JSONs
provide the parameters, predictions, data identities, and checks. The NPZ stores
the two calibration and three held-out 29-um trajectories in original x50/L
coordinates (columns Rac, bound MCA).

The regenerated 2D tracks use this new fitted parameter set, duration 10000 s,
sampling every 1 s and relative noise fractions 0.1/0.3/0.6. See
[2D overview](../../results_fitted_2d_self_activation/trajectories_overview.png).
