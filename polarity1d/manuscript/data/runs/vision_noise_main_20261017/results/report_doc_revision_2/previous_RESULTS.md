# Replicated noisy rings: polarity memory, magnitude and local geometry

**All 243 prescribed integrations completed and all 216 noisy realizations are
included.** The experiment contains 12 independent realizations in each of the
18 model/circumference/strength groups, plus 18 deterministic mean and nine
zero-input controls. Every realization reaches 32000 native time units. This
report supplies the noise component of the manuscript mechanism-and-measurement
checkpoint; it does not assign a cell velocity or calibrate migration persistence.
Generated 2026-10-17T21:00:53.085715+00:00.

**Effective executed protocol:** $T=32000$; scalar and site-descriptor cadence 5;
saved-primary-profile cadence 40; shared maximum lag 8000. The operative
[frozen configuration](../../../../data/runs/vision_noise_main_20261017/configuration.json)
and [frozen analysis plan](../../../../data/runs/vision_noise_main_20261017/plan_freeze.json)
supersede prospective timing/cadence wording in the design-stage notes.
The [frozen scientific design](../../../../data/runs/vision_noise_main_20261017/protocol/scientific_design.md)
records the prospective rationale; it is not a live producer-status page.

At the higher predeclared input, **larger rings retain more directional memory
in Wave-pinning and Holmes while largely maintaining their resultant
magnitude**. At lag 8000, the direction cosine rises from 0.593
[0.481, 0.693] at C20 to 0.972 [0.967, 0.978] at C60 in Wave-pinning, and from
0.315 [0.147, 0.463] to 0.738 [0.614, 0.836] in Holmes. These are pointwise 95%
whole-seed bootstrap intervals conditional on the fixed model and input law;
more than 99% of directional pairs remain eligible in these comparisons.

**Angular memory is not a universal measure of slower physical site motion.**
Wave-pinning's eligible physical RMS displacement at lag 8000 decreases from
3.25 [2.71, 3.74] to 2.25 [2.03, 2.48] native length units between C20 and C60
(endpoint coverage 99.2% and 94.6%). Holmes instead gives 4.57 [3.88, 5.17]
and 5.54 [4.25, 6.70], with overlapping intervals and coverage falling from
92.1% to 37.3%. Geometry and endpoint conditioning therefore remain essential
to interpreting Holmes. Goryachev shows a separate regime: its C20 high-mean
control loses polarity, whereas C40/C60 retain it. That comparison cannot be
presented as noise acting around the same intact polarized state. Its
conservative single-peak gate also leaves insufficient physical-displacement
coverage; an isolated one-pair value is not a population drift estimate.

The lower strength produces weak lag decorrelation over this horizon. **Sixteen
of eighteen full-window noisy conditions do not cross $1/e$ by lag 8000.**
The two sampled crossings are Holmes C20 at high input and the distinct
high-input Goryachev C20 regime. The raw curves and signed restricted areas
therefore support a finite-window polarity-memory result, rather than fitted
unlimited or stationary persistence times, calibrated biological times or cell
trajectories.

## Question and declared comparison

Does enlargement change measurable whole-boundary polarity memory under local
fluctuating activation? Wave-pinning, Goryachev and Holmes begin from their exact
qualifying C20/C40/C60 basals at spatial step 0.05. The positive input law,
preparation, kinetic parameters, two strengths and 12 independent seeds per
condition were frozen before the main producer launched. Both strengths are
retained, including conditions in which the deterministic mean alters the
polarized state.

The native external activation input is

$$
I(s,t)=fA_m\exp[Z(s,t)-1/2],\quad
f\in\{0.001,0.01\},\quad
A_m\in\{0.1,6.4,1\}.\tag{1}
$$

The Gaussian field has unit local ensemble variance, physical spatial
covariance width 1.5 and OU knot correlation time 50. Knots are separated by
5 and use analytically variance-normalized linear interpolation. No
per-realization mean removal, amplitude normalization or state clipping is
introduced. This is the declared interpolated driving law, not an exact
continuous-time OU path. Local noise statistics are fixed within each model;
total input and whole-ring angular input modes change with circumference.
The three native input/time scales are not calibrated equal biological exposures.

## Realized forcing and whole-ring modes

The [216 realized-input audits](supporting/data/vision_noise_input_audits.csv)
link every worker noise-audit JSON and per-knot mode trace with exact hashes.
They retain realized local means/variances, Gaussian and lognormal spatial and
temporal covariances, integrated zeroth-mode variance and first-angular-mode
second moment alongside finite-ring expectations. The table averages the 12
realized audit summaries within each declared condition; it never rescales
any input realization. Expected normalized local mean and variance are 1 and
$e-1$, respectively. Zeroth-mode variance and first-mode second moment below
are divided by $\mu^2$; each entry is observed / theoretical.

| Model | C | f | Local mean / μ | Local variance / μ² | Zeroth-mode variance / μ² | First-mode second moment / μ² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Wave-pinning | 20 | 0.001 | 1.01 | 1.805 | 113.5 / 110.9 | 107 / 101.3 |
| Wave-pinning | 40 | 0.001 | 1.002 | 1.748 | 222.1 / 221.8 | 217.4 / 216.8 |
| Wave-pinning | 60 | 0.001 | 1.007 | 1.735 | 331.9 / 332.8 | 337.8 / 329.4 |
| Wave-pinning | 20 | 0.01 | 0.9996 | 1.674 | 109.4 / 110.9 | 97.83 / 101.3 |
| Wave-pinning | 40 | 0.01 | 1.003 | 1.735 | 229 / 221.8 | 216.5 / 216.8 |
| Wave-pinning | 60 | 0.01 | 0.9984 | 1.75 | 342.2 / 332.8 | 332.3 / 329.4 |
| Goryachev | 20 | 0.001 | 1.014 | 1.815 | 115 / 110.9 | 108.1 / 101.3 |
| Goryachev | 40 | 0.001 | 1.004 | 1.811 | 234.2 / 221.8 | 226.5 / 216.8 |
| Goryachev | 60 | 0.001 | 0.9977 | 1.692 | 328.6 / 332.8 | 323.2 / 329.4 |
| Goryachev | 20 | 0.01 | 1.009 | 1.737 | 113.1 / 110.9 | 101.6 / 101.3 |
| Goryachev | 40 | 0.01 | 1.007 | 1.767 | 229.7 / 221.8 | 225 / 216.8 |
| Goryachev | 60 | 0.01 | 0.9987 | 1.705 | 323.5 / 332.8 | 326.3 / 329.4 |
| Holmes 3 | 20 | 0.001 | 0.999 | 1.746 | 115.7 / 110.9 | 102.5 / 101.3 |
| Holmes 3 | 40 | 0.001 | 0.9971 | 1.68 | 213.9 / 221.8 | 212 / 216.8 |
| Holmes 3 | 60 | 0.001 | 1.005 | 1.724 | 329.5 / 332.8 | 333.1 / 329.4 |
| Holmes 3 | 20 | 0.01 | 0.9998 | 1.747 | 112.3 / 110.9 | 102.6 / 101.3 |
| Holmes 3 | 40 | 0.01 | 1.008 | 1.764 | 223.1 / 221.8 | 224.3 / 216.8 |
| Holmes 3 | 60 | 0.01 | 1.001 | 1.717 | 333.7 / 332.8 | 329.2 / 329.4 |

Thus equal local covariance is not equal whole-ring angular forcing or equal
total input at different circumferences. The complete technical archive verifier
regenerated every realization from its frozen source/seed; finite realized
moment deviations remain observations rather than normalization targets.

## Initial memory and magnitude

![Replicated initial-axis memory](vision_noise_initial_memory.png)

Solid curves are the independent-realization mean; shading is the 2.5–97.5%
percentile interval from 1000 complete-realization bootstrap resamples. Dashed
curves are deterministic mean-input controls. The same full range is used
across panels and contains every displayed mean and interval. No individual
seed or strength is removed because it gives an inconvenient result.

![Replicated vector magnitude](vision_noise_magnitude.png)

For the native active field $u(s,t)$, the measured moments and aligned vector are

$$
M_k(t)=\int_0^C u(s,t)e^{2\pi i k s/C}\,\mathrm ds,
\qquad P(t)=\frac{M_1(t)}{M_1(0)},\qquad
q_k(t)=\frac{|M_k(t)|}{M_0(t)}.\tag{2}
$$

The activity vector is initially $(1,0)$ and its relative magnitude is $R=|P|$;
R can exceed one because its denominator is the initial moment. Initial-axis memory is $\langle\operatorname{Re}P(t)\rangle$;
it can change through direction, magnitude or both. A small resultant can
coexist with active opposing sites and cannot by itself establish homogeneous
local activity. Total activity, $q_1$, $q_2$, and peak/region diagnostics remain
in the seedwise tables. Uniform mean input can act on a spatially nonuniform
inactive field, so its direct contribution to the active moment need not vanish.

## Lag memory and direction

![Replicated vector and direction lag correlations](vision_noise_correlations.png)

The solid statistic is the raw uncentered complex-vector lag correlation: the
ratio of ensemble-averaged within-realization endpoint cross-products to the
square root of both endpoint energies. It retains weak resultants and subtracts
no temporal or ensemble mean vector. Each dashed curve is the conditional
unit-direction cosine with both endpoint magnitudes at least 0.1. These answer
different questions; an amplitude-memory change is not automatically slower
turning. Confidence bands use whole seeds, never time origins as independent
replicates. Complete direction intervals and magnitude gates 0.05/0.2 are in
[the direction table](supporting/data/vision_noise_direction_correlations.csv).

![Direction eligibility and physical geometry](vision_noise_coverage_geometry.png)

Direction-pair coverage is shown explicitly. Physical displacement uses producer
scalar descriptors every 5 units, independently validated against all saved
40-unit primary profiles. Both endpoints must have one current half-contrast
region, at most one peak with native prominence at least 0.05, contrast at
least 0.05 and R at least 0.1. Missing eligible pairs remain unavailable. The
same physical-displacement scale is used throughout; angular persistence can
increase merely because the same physical displacement subtends a smaller
angle on a larger ring. These finite-cadence endpoint comparisons do not prove
continued site identity, define a diffusion coefficient or describe cell travel.
The tables preserve all eligibility fractions and bootstrap-defined fractions.
The absolute prominence gate is deliberately conservative: small secondary
peaks can fail it while a clear dominant half-contrast incumbent remains.
Missing eligible pairs do not establish instability, and every detected
0.05-prominence peak is not labelled a functional competitor or protrusion.
Raw dominant-site descriptors and region counts remain available. Goryachev's
physical RMS fragments are not displayed as a population curve: at high input
and lag 8000, C20 has only one eligible pair among 57612, while C40/C60 have
none. The collapsed one-pair interval is not a precision estimate. Complete
computed values and bootstrap-defined fractions remain in the unmodified CSV;
this restriction does not prevent observing its strong dominant incumbent.

## Restricted memory and adaptation

![Full and half-window restricted integrals](vision_noise_restricted_integrals.png)

The integral is the raw signed area through lag 8000, including negative
correlations. It is an observation-window descriptor in native time units,
not a fitted exponential persistence time. Full/first-half/last-half windows
use the same maximum lag and no chosen burn-in. Differences between them expose
continuing adaptation; agreement alone does not prove stationarity. A curve
that does not cross $1/e$ within this horizon does not supply a resolved
exponential time. Weak endpoint-energy statuses remain explicit in
[the complete integral table](supporting/data/vision_noise_persistence.csv).

| Model | C | f | Final mean Re P [95% interval] | Restricted area [95% interval] | Direction cosine at lag 8000 | Eligible direction pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Wave-pinning | 20 | 0.001 | 0.9861 [0.9737, 0.9964] | 7972 [7965, 7977] | 0.9925 | 1 |
| Wave-pinning | 40 | 0.001 | 0.9966 [0.9925, 0.9994] | 7996 [7995, 7998] | 0.9991 | 1 |
| Wave-pinning | 60 | 0.001 | 0.9981 [0.9955, 0.9996] | 7999 [7999, 7999] | 0.9997 | 1 |
| Wave-pinning | 20 | 0.01 | 0.2875 [-0.09834, 0.6248] | 6257 [5848, 6621] | 0.5931 | 1 |
| Wave-pinning | 40 | 0.01 | 0.6505 [0.4306, 0.8338] | 7627 [7533, 7718] | 0.9055 | 1 |
| Wave-pinning | 60 | 0.01 | 0.9118 [0.8542, 0.9617] | 7885 [7867, 7904] | 0.9723 | 1 |
| Goryachev | 20 | 0.001 | 0.9569 [0.9462, 0.965] | 7996 [7996, 7997] | 0.9998 | 1 |
| Goryachev | 40 | 0.001 | 0.981 [0.9791, 0.9829] | 8000 [8000, 8000] | 1 | 1 |
| Goryachev | 60 | 0.001 | 0.9858 [0.9843, 0.9871] | 8000 [8000, 8000] | 1 | 1 |
| Goryachev | 20 | 0.01 | -0.02096 [-0.1454, 0.1191] | 220 [158.1, 271.5] | 0.02724 | 0.7893 |
| Goryachev | 40 | 0.01 | 0.7943 [0.778, 0.8117] | 7945 [7940, 7951] | 0.998 | 1 |
| Goryachev | 60 | 0.01 | 0.8412 [0.8255, 0.8558] | 7992 [7992, 7993] | 0.9997 | 1 |
| Holmes 3 | 20 | 0.001 | 0.9131 [0.8726, 0.9512] | 7895 [7861, 7923] | 0.9738 | 1 |
| Holmes 3 | 40 | 0.001 | 0.9846 [0.9781, 0.9908] | 7979 [7976, 7982] | 0.9952 | 1 |
| Holmes 3 | 60 | 0.001 | 0.9936 [0.9899, 0.9971] | 7991 [7987, 7993] | 0.9976 | 1 |
| Holmes 3 | 20 | 0.01 | 0.0446 [-0.4139, 0.4882] | 4304 [3581, 4961] | 0.3155 | 0.9983 |
| Holmes 3 | 40 | 0.01 | 0.2764 [-0.08693, 0.6536] | 6332 [5974, 6753] | 0.5825 | 0.9967 |
| Holmes 3 | 60 | 0.01 | 0.3342 [-0.02332, 0.6471] | 7103 [6722, 7353] | 0.7381 | 0.9918 |

## Deterministic mean regimes

The mean controls determine whether noise is probing the original kind of
polarized state or a different forced-state regime. In particular the pilot's
high-input Goryachev C20 mean destroyed its polarized resultant. The complete
main controls are retained here rather than counting that regime as a clean
noise-specific directional-persistence comparison.

| Model | C | f | Mean-control final R | Mean-control final total activity / initial |
| --- | ---: | ---: | ---: | ---: |
| Wave-pinning | 20 | 0.001 | 1 | 1.001 |
| Wave-pinning | 40 | 0.001 | 1 | 1.001 |
| Wave-pinning | 60 | 0.001 | 1 | 1.001 |
| Wave-pinning | 20 | 0.01 | 0.9999 | 1.011 |
| Wave-pinning | 40 | 0.01 | 1.001 | 1.011 |
| Wave-pinning | 60 | 0.01 | 1.001 | 1.011 |
| Goryachev | 20 | 0.001 | 0.9562 | 0.9966 |
| Goryachev | 40 | 0.001 | 0.9798 | 0.9993 |
| Goryachev | 60 | 0.001 | 0.9853 | 0.9997 |
| Goryachev | 20 | 0.01 | 3.241e-17 | 0.9288 |
| Goryachev | 40 | 0.01 | 0.743 | 0.9882 |
| Goryachev | 60 | 0.01 | 0.8361 | 0.9962 |
| Holmes 3 | 20 | 0.001 | 1.001 | 1.002 |
| Holmes 3 | 40 | 0.001 | 1.001 | 1.002 |
| Holmes 3 | 60 | 0.001 | 1.001 | 1.002 |
| Holmes 3 | 20 | 0.01 | 1.013 | 1.023 |
| Holmes 3 | 40 | 0.01 | 1.013 | 1.023 |
| Holmes 3 | 60 | 0.01 | 1.013 | 1.023 |

## Local histories selected by the frozen rule

![Every condition medoid](vision_noise_selected_histories.png)

One medoid is shown for every noisy model/C/strength condition. The frozen
descriptors are final initial-axis memory, minimum magnitude, time-mean
magnitude and the profile fraction with multiple current half-contrast regions.
The sum of within-group median deviations, scaled by descriptor IQR with a
unit zero-IQR fallback, determines the medoid; ties use case identifier.
These are descriptive representatives, not the histories most favourable to a
size effect. Each model row has one raw native activity colour scale. The
same angular axis spans different physical lengths at different C. Exact seeds,
archives and hashes are retained in [the selection table](supporting/data/vision_noise_representative_histories.csv).

![Lowest-resultant diagnostic](vision_noise_lowest_resultant_diagnostic.png)

The additional diagnostic is the globally lowest sampled noisy resultant,
selected independently of which size trend it suggests. It occurs in
Goryachev at C=20,
f=0.01, seed 202610180093,
$t=20635$: $R=0.0010749$,
$q_1=0.0011594$, $q_2=0.28259$ and native
contrast 4.7001. The profile snapshot is the nearest
saved time, $t=20640$; it is not passed off as the
exact scalar-event state. The matched mean control and raw local profile
separate weak whole-ring asymmetry from local activity loss or cancellation.
No dominant-site switch alone is labelled incumbent replacement.

## Corresponding experimental measurement

A whole-boundary activity reporter could provide the same first angular moment,
its magnitude and its temporal correlations. Local site locations, intensities,
region counts and physical displacement would distinguish persistent asymmetry
from local wandering, multisite cancellation or loss of contrast. This is a
polarity-memory readout motivated by the enlargement/persistence hypothesis;
it is distinct from measured migration trajectories and conditional OU velocity
persistence estimates in the manuscript's empirical component.

## Complete execution and numerical evidence

The independent analysis verified 22
frozen source files, 18 original
provenance/input files and 1944
worker-output hashes. Every exact basal start is checked. At each saved primary profile, raw moments
and periodic half-contrast region/peak descriptors are independently reproduced;
scalar descriptors supply the intervening five-unit site series.
Maximum retained raw-moment discrepancy is
2.56e-13; maximum conserved-pool
relative drift is 8.52e-12.
All 243 complete statuses and the prescribed 12 independent seeds per group
are required; a failed or missing trajectory cannot be silently censored.
The signed integrals, full raw curves, independent-seed intervals and weak
vector/pair statuses are retained in [the analysis audit](supporting/metadata/vision_noise_analysis_audit.json)
and [the seed summaries](supporting/data/vision_noise_seed_summary.csv).

A [separate numerical amendment](../../../../data/runs/vision_noise_main_20261017/results/analysis_numerical_amendment/AMENDMENT.md)
preserves the original frozen estimator while replacing cancellation-prone tail
subtraction and weak FFT products by positive suffix sums and direct endpoint
products. It leaves the producer, all trajectories and scientific definitions
unchanged; no normalized correlation is clipped. All 1601 lags in all three
windows of all 27 actual mean/zero controls pass independent direct-product
regression, with maximum absolute correlation error below $1.281\times10^{-12}$.
Extreme tiny-tail/large-range fixtures also pass, with genuinely unresolved
endpoint energies left undefined. A separate
[center-observability adapter amendment](../../../../data/runs/vision_noise_main_20261017/results/statistics_centroid_amendment_20261017/AMENDMENT.md)
retains all field/count checks but compares positions only when a current
component exists. A flat control otherwise carries its previous scalar- versus
profile-cadence position, which is not an observed center. Such carried labels are
flagged as unobserved current centers in the validation audit and cannot meet
the unchanged direction/site gates.
All 243 cases and 194643 saved profiles pass the actual strict adapter validation;
no center tolerance or mathematical readout is changed. Seed/group tables retain numerical fallback
and sensitivity fields alongside weak-state statuses.

The frozen pilot numerical checks distinguish fixed-forcing ODE accuracy from
sensitivity to the interpolated input path. Tightening the ODE solver changes
P by at most $7.71\times10^{-8}$ over the three checked cases. Coupling the
5-to-2.5-unit interpolation changes maximum P by 0.00546/0.00176/0.01807 in
Wave-pinning/Goryachev/Holmes, with raw-profile differences of
2.31%/0.219%/12.74% of basal contrast. Those are finite-pilot path sensitivity
checks, not a theorem of long-run stochastic convergence or a claim that
local histories are numerically identical. The [reviewed pilot](../vision_noise_pilot_20261017/RESULTS.md)
retains every control, check and interpretation limit.

Allocated storage measured after generation: native main group
2,802,688,000 bytes and this report group 641,187,840 bytes before
this final prose/manifest. The root data-manifest owner receives exact measured
sizes; scientific importance is the complete replicated stochastic comparison.

A fresh-context actual-result critique will assess the curves, local histories,
uncertainty and scientific claims before this report is folded into the final
three-figure checkpoint. Output generation alone does not establish an accepted
manuscript figure or a general enlargement effect.
