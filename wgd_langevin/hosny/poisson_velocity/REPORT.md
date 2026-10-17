# Joint Poisson fit: results and interpretation

The primary model jointly fits the untreated/control baseline and both knockdowns,
with HIF allowed only to reduce kick directionality and RhoA allowed only to reduce
velocity memory. It does not fully explain the reported distributions. The most
clear-cut discrepancy is persistent oxygen preference after RhoA KD, especially
in PACCs. Adding RhoA directionality as a jointly fitted parameter approximately
halves the PACC error. This is a descriptive model comparison, not a molecular
causal identification or a statistical significance result.

The standalone figure is [results/hosny_poisson_joint.png](results/hosny_poisson_joint.png).
The plotted parameters, tracks and distribution bands are always from the primary
model; panel C alone compares independently refitted expanded models. The plotted
scores are evaluated with fresh, larger simulated populations.

## What is fitted

For each PACC class: baseline tau, Poisson rate, Rayleigh kick scale, and von Mises
concentration; an HIF concentration no greater than baseline; a RhoA tau no greater
than baseline. Six jointly optimized quantities per class. Normoxia has three
parameters, with uniform directions fixed. Untreated and both siNTC controls are
three separate source constraints on the same baseline. Knockdown data inform
shared baseline parameters during the same optimization. All cells within each
condition share parameters; there is no permanent cell-to-cell heterogeneity.

The objective includes all three quartiles of sampled path length, directness and
parallel displacement/path. Endpoints and perpendicular fractions are held out.
The experimental cadence is 0.5 h for 16 h. Parameter scale is profiled analytically
against every relevant path constraint at each evaluation, so this is still a joint
fit. See README and metadata for exact loss, bounds and numerical procedures.

## Results from independent 32,768-cell simulations

| Joint model | non-PACC error | PACC error |
|---|---:|---:|
| Primary: HIF concentration only; RhoA tau only | 0.0489 | 0.0508 |
| Also allow HIF tau | 0.0413 | 0.0475 |
| Also allow RhoA concentration | 0.0403 | 0.0253 |
| Allow both additions | 0.0336 | 0.0218 |
| All three conditions have free parameter sets | 0.0250 | 0.0106 |

The last column does not result from holding the primary baseline fixed. Every row
refits its baseline and both KDs jointly. Added HIF tau and RhoA concentration are
unrestricted within the search bounds; the fully free model also releases primary
inequalities and rate/amplitude sharing. Parameter counts are 6, 7, 7, 8, and 12 per
PACC class. Scores have no complexity penalty and are not likelihoods.

Under the primary fit, non-PACC baseline tau is about 0.040 h and RhoA tau 0.033 h;
PACC baseline tau is about 0.326 h and RhoA tau 0.271 h. HIF concentration falls from
about 0.148 to 0.013 in non-PACCs and 0.310 to 0.076 in PACCs. Baseline rates are
about 0.71/h and 0.60/h. These are representative numerical optima, not identified
physical constants. Full values and jump scales are in panel A and fits.csv.

HIF directionality reduction captures loss of directional bias but leaves the
non-PACC path distribution too high. Adding HIF memory reduces that discrepancy.
RhoA memory reduction lowers path distance while retaining the directional kick
preference, leaving too much forward movement and directness. The RhoA-direction
expansion estimates near-zero RhoA concentration in non-PACCs and about 0.024 in
PACCs, compared with refitted baseline values around 0.229 and 0.512 respectively.

Poisson intermittency permits substantially broader path distributions than the
previous Gaussian OU fit, but important mismatches remain, especially normoxic and
non-PACC path IQRs. The displayed bands are population distributions, not confidence
intervals or variability between replicate experiments.

## What is weakly identified

The original 0.04 h lower memory bound was reached, prompting an extension to
0.005 h. Rate and memory profiles jointly refit all remaining shared and KD
parameters. For short tau, kick amplitude can increase approximately inversely
with tau, preserving displacement per kick while changing instantaneous peak speeds.
Thus the sampled trajectories do not establish the huge short-lived speeds implied
by some point estimates. The normoxic optimum moves to about 0.01 h after expanding
the range, with almost unchanged fit quality. Non-PACC profiles are also very flat
at short memory times; PACC profiles admit a broad range of alternatives.

In particular, do not interpret the primary PACC/non-PACC tau ratio as a robustly
identified biological persistence ratio. The alternative normalization and weighting
fits move tau substantially (sensitivity.json). The supplied notes leave the FMI
denominator unresolved; this materially affects inferred parameters.

Three optimizer seeds and polishing were used. All reported optimizer runs converged,
and nested bank losses pass ordering checks. Independent loss differs from bank loss
by approximately 0.0009–0.0051, depending on model. Small differences between nearby
optima or profile points are below this numerical scale. No claim of a proven global
minimum, formal confidence bounds, or precise discrimination among near-equal fits
is warranted. Parameter profiles are diagnostics, not inferential confidence regions.

## Other observations and predictions

Endpoint and perpendicular distributions are compared with every corresponding
source summary, including controls. Both knockdowns reduce the predicted fraction
of positive/high-endpoint PACCs, matching the qualitative note: about 21% baseline
versus 12% after either KD. The cutoff is the shared baseline PACC endpoint Q75;
positive means forward displacement >0. No numeric experimental fractions were
supplied, so this is a qualitative check only.

Rayleigh tests on simulated samples of 60 cells check the observed pattern of
significant baseline/control directionality and nonsignificant KD directionality.
The primary model often retains significant directionality after RhoA KD. These
checks assume independent cells; the three biological experiments in the notes
cannot be modeled without their separate data. Individual experimental p-values
are retained beside predicted p-value ranges in rayleigh_checks.csv.

MSD curves and cue removal (set concentration to zero, preserve other parameters)
are generated predictions, not observations fitted or verified against supplied
MSD data. Sample trajectories use the same first-hour prefixes as the 16 h tracks.

The homogeneous stationary model does not address radial oxygen response, temporal
adaptation, cell area or area-based retrospective PACC classification. The supplied
radial means have unresolved population composition; temporal and area correlations
lack numerical targets. Those observations have not been silently treated as fitted.

## Validation and provenance

All nine simulator checks passed, including analytic MSD/mean/covariance and exactly
matching trajectories when the observation grid is refined with the same events.
The final audit also verified parameter sharing/inequalities, nested objective
ordering and source hashes. No external papers or online data were consulted.

SLURM jobs 20511931 (initial bank/fit) and 20511938 (expanded memory range and final
fit) completed on the red partition, small QoS, four CPUs. Final analysis took about
115 seconds using the existing/extended simulation bank; PNG rendering followed.
All source, seeds, outputs and reproduction instructions are retained in this folder.
