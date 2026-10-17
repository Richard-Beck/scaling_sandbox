# Fresh actual-result critique: replicated noisy rings

Reviewed 17 October 2026 after all **243 native trajectories**, the independent
technical verification, complete statistical exports and the revised report
rendering finished. Scope includes the manuscript vision and instructions,
prospective design, reviewed pilot, both separately preserved analysis
amendments, all seven actual PNGs, the revised [report](../../RESULTS.md),
seedwise tables, independent-seed intervals, realized-input audits and exact
representative-history selection. The revised report reviewed here was
generated at 21:00:53 UTC. This review ran bounded analysis-only checks and
launched no integration.

**Disposition: no essential scientific, statistical or display correction
remains for this noisy-ring component of the review checkpoint.** It supports
a finite-window, within-model comparison of measurable boundary-activity
memory under the declared input. The result is model- and regime-dependent;
it is not a general calibrated migration-persistence law or acceptance of a
finished manuscript figure. All three models, both strengths, all 216 noisy
realizations and all 27 deterministic controls remain retained.

## Actual findings and their limits

**Wave-pinning and Holmes show a resolved higher-input angular-memory
contrast while broadly retaining activity-vector magnitude.** At lag 8000,
the independent checks reproduce these conditional direction cosines and
pointwise 95% whole-seed bootstrap intervals:

| Model, $f=0.01$ | C20 | C60 |
| --- | ---: | ---: |
| Wave-pinning | 0.593 [0.481, 0.693] | 0.972 [0.967, 0.978] |
| Holmes | 0.315 [0.147, 0.463] | 0.738 [0.614, 0.836] |

Direction-pair coverage exceeds 99% for these four comparisons. Final mean
relative magnitude is approximately 1.000/1.000 in Wave-pinning and
1.020/0.978 in Holmes at C20/C60, whereas final initial-axis memory is
0.287/0.912 and 0.045/0.334 respectively. Thus loss of initial-axis memory
is not equivalent to wholesale loss of the resultant. Rare weak moments and
changes in local organization still occur; ensemble means do not establish
that every realization stays strongly polarized.

**Angular memory and conditional physical motion give different comparisons.**
Wave-pinning's eligible physical RMS endpoint displacement at lag 8000 falls
from 3.249 [2.71, 3.74] at C20 to 2.254 [2.03, 2.48] at C60, with coverage
99.2%/94.6%. The earlier one-seed pilot's opposite physical comparison does
not carry over to this ensemble. Holmes gives 4.568 [3.88, 5.17] and
5.539 [4.25, 6.70], with overlapping intervals and coverage falling from
92.1% to 37.3%. Its stronger angular memory cannot be relabelled slower
physical site wandering. These comparisons condition only on eligible
endpoints, not continued site ancestry; they are neither diffusion
coefficients nor cell displacements.

**Goryachev C20 at the higher mean is a different forced-state regime.**
The uniform mean alone ends with $R=3.24\times10^{-17}$ and a numerically
flat profile, whereas C40/C60 retain $R=0.743/0.836$. The noisy C20 condition
has final ensemble mean $R=0.231$ and initial-axis memory $-0.02096$.
Its short restricted correlation area cannot be interpreted as noise-driven
turning of the same intact pole found in larger rings. The globally lowest
noisy moment is worker 0093 at $t=20635$: $R=0.001075$, while contrast is
4.700 and $q_2=0.2826$. Its displayed native profile demonstrates substantial
local activity despite cancellation of the first moment. The nearest profile
time, 20640, is explicitly distinguished from the scalar event.

Goryachev's frozen single-region/prominence criterion supplies essentially no
high-input physical-displacement comparison: at lag 8000 C20 has **one
eligible pair out of 57612 possible pairs**, while C40/C60 have none. A
collapsed bootstrap interval from that single observed pair conveys no
population drift precision. Replacing its fragmented physical plot by an
explicit availability note is appropriate; the unchanged raw values, statuses
and coverage remain in the CSV. No new statistical gate was introduced.

**Most decorrelation remains unresolved within the declared lag horizon.**
Sixteen of eighteen full-window noisy curves do not cross $1/e$ by lag 8000.
The two crossings are higher-input Holmes C20 (descriptive crossing about
6921) and the distinct Goryachev C20 regime (about 236). Higher-input
Wave-pinning restricted areas are 6256.5, 7627.2 and 7884.6 at C20/40/60;
Holmes gives 4303.9, 6331.5 and 7102.8. They remain signed restricted areas,
not fitted or extrapolated persistence times. Holmes C20's first/last-half
areas are 3975.8/4481.0, exposing window dependence. Uncentered correlations
also retain the prepared mean and continuing adaptation. No chosen burn-in,
exponential fit or stationarity claim is warranted. The lower-strength weak
responses remain visible on honest common scales.

## Independent statistical, selection and provenance checks

The [reproducible audit](review_actual_results.py) and
[complete checks](FRESH_RESULT_CHECKS.json) pass:

- All 243 cases reach the prescribed 32000 horizon on five-unit scalar grids;
  saved profiles use cadence 40. There are 45 complete groups, with exactly
  12 independent noisy realizations in each of 18 groups. No model, strength,
  failed case or outcome has been silently removed.
- All ten authoritative statistical CSV hashes match the closed analysis
  audit. All 1,167,129 seed-lag rows are accounted for. Explicit endpoint sums
  at lags 0, 5, 1000, 4000 and 8000 in all three windows of every seed agree
  with exported seed correlations within $1.36\times10^{-13}$.
- All 135 signed restricted areas and all 54 noisy integral bootstrap
  intervals reproduce exactly from retained seed-level endpoint products and
  energies. Independently resampled final trace intervals and lag-8000
  direction/site intervals, values and pair coverage agree to roundoff.
  Uncertainty resamples whole independent seeds, never time origins.
  Intervals are pointwise and conditional on the fixed models/protocol;
  they are not simultaneous bands or uncertainty over biological calibration.
- All 18 medoids reproduce from raw scalar descriptors and independently
  counted current half-contrast regions in the saved profiles. The globally
  lowest-resultant diagnostic is also exactly reproduced. All 18 medoids are
  shown, with one native activity colour scale within each model; examples
  are not selected for a favourable size effect.
- All 22 native frozen executable sources and 11 protocol copies retain their
  declared hashes. The original prospective analysis remains unchanged.
  The technical verifier passes all 243 archives, including exact regeneration
  of all 216 forcing realizations; the closed analysis checks 1944 worker-output
  hashes. Maximum native conserved-pool relative drift is
  $8.52\times10^{-12}$.

The final operative analysis source is
`61622d8d21102397d38e9dc4050cc83852a60169c0e08fdedff494e1ab415dc7`.
Its [numerical correction](FRESH_NUMERICAL_AMENDMENT_REVIEW.md) preserves the
declared estimator using positive tail energies and direct products where FFT
roundoff is inadequate; no correlation is clipped. The separate
[flat-center adapter review](FRESH_FLAT_CENTER_AMENDMENT_REVIEW.md) preserves
the original strict tolerance for observable centers. The **completed
all-243** export now confirms zero observable-center discrepancy and no
eligible flat frames. Only the homogeneous mean control has flat cached
labels: 5572 scalar samples and 697 saved profiles. They do not enter physical
endpoint measurements. The original failed analysis and both correction
records remain preserved.

## Realized input and final display assessment

The [input-audit check](FRESH_INPUT_AUDIT_CHECKS.json) independently reconstructs
all 216 realized local mean/second-moment/variance and integrated zero/first-mode
summaries from native knot-mode traces, verifies all 432 source hashes and
checks exported Gaussian/lognormal covariance values against their native
records. The largest scaled reconstruction discrepancy is
$8.39\times10^{-16}$. The report directly links every audit and states that
finite-realization deviations are observations, not normalization targets.

Fixed local input statistics do not fix total activation or whole-ring angular
forcing across circumference. Physical covariance width 1.5 and normalized
OU-knot interpolation with knot spacing 5/time scale 50 belong to the declared
law; the interpolation is not an exact stationary continuous-time OU process.
Native input and time units remain uncalibrated across models. These are
limits on generalization, not reasons to tune the completed experiment.

All seven actual PNGs were visually inspected. The corrected correlation
labels no longer overlap; the revised physical panel includes Wave-pinning
and Holmes bootstrap bands and explicitly exposes Goryachev's insufficient
support. Common scales preserve the weak lower-strength responses. Magnitude,
direction/coverage, full/half windows, all condition medoids and the local
cancellation diagnostic make the scientific distinctions assessable.
The [render record](FRESH_RENDER_CHECKS.json) verifies all 23 report artifacts
and 457 source references after the derived-only polish; five unchanged PNGs
remain byte-identical to their predecessor versions. All statistical CSVs and
native evidence remain unchanged by that polish.

The experimental bridge is appropriately a whole-boundary activity measurement
with vector magnitude, local contrast, site organization and physical positions.
It supplies an observable motivated by the WGD enlargement/persistence
hypothesis. Observed migration summaries, conditional OU velocity persistence,
these ring activity readouts and the still-pending new WGD motility experiments
remain distinct. No cell movement law or further campaign is needed to make
this bounded result reviewable.
