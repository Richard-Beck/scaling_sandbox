# Fresh actual-result critique: bounded noisy-ring pilot

Reviewed 17 October 2026 after all 30 pilot workers and all six numerical
audits completed. Scope: repository/manuscript instructions and vision, the
prospective design and design critique, frozen implementation and archives,
[actual report](../../RESULTS.md), all five rendered PNGs, analysis source and
tables, and the [numerical comparison](../../../../../../data/runs/vision_noise_numerical_audit_20261017/NUMERICAL_COMPARISON.json).
This review performed bounded read-only checks; it changed no simulation
source, inputs or trajectories and launched no integration.

**Disposition: the proposed main ensemble is scientifically interpretable as
an experiment on boundary activity organization under the declared local
fluctuations. No remaining essential display or numerical correction blocks
its freeze.** Retain all three models, both strengths and every deterministic
control. The high-input Goryachev response includes a change of forced-state
regime; a common directional-persistence estimate across all conditions is
therefore unavailable. This pilot supplies no ensemble size inference or
calibrated migration prediction.

## Actual evidence and essential interpretation constraints

1. **The mean controls determine the Goryachev interpretation.** At C20 and
   the higher input, the uniform mean alone ends with relative first-moment
   magnitude $R=3.01\times10^{-13}$ and activity contrast
   $2.95\times10^{-12}$. C60 retains $R=0.836$ under its matched mean. The
   noisy C20 trajectory reaches $R=0.00721$ at $t=3770$, while contrast is
   3.281, $q_2=0.1970$, and two half-contrast regions with three prominent
   peaks remain. Weak noisy first moment and homogeneous activity are distinct
   states. Preserve this regime distinction rather than attributing the full
   short/long contrast to noise-driven turning of otherwise comparable poles.

2. **Angular memory requires the physical-coordinate evidence.** At lag 1000,
   high-input Wave-pinning has C20/C60 direction cosines 0.9699/0.9949, while
   eligible physical RMS displacements are 0.7884/0.9853 native length units
   with pair coverage 98.7%/90.7%. Thus greater angular memory coexists here
   with a larger conditional physical displacement. Goryachev has zero
   eligible high-input single-region pairs at that lag; its low-input coverage
   is only 6.0%/2.6%. The embedded coverage/geometry figure now makes these
   limitations visible. Endpoint displacements do not establish site ancestry,
   incumbent replacement or a migration law.

3. **Preserve preparation and finite-window limitations.** Most curves remain
   close to one through lag 1000. Even the high-input Holmes C20 restricted
   integral depends materially on the window: 900.7 over the full interval,
   849.0 in the first half and 974.4 in the last half. These descriptors combine
   the prepared initial state, magnitude dynamics and axis dynamics. Keep raw
   signed curves, full/half windows, direction gates with pair coverage and
   unresolved-decay statuses; do not select burn-in or force exponential fits.
   Independent seeds, rather than time origins, supply ensemble uncertainty.

4. **The comparative presentation is adequate after correction.** The shared
   memory range contains the actual extrema, $-0.419$ to 1.132, and the
   magnitude/correlation ranges contain their data. Every noisy history is
   shown, with a shared native activity colour range within each model. The
   one-seed medoids are correctly described as individual histories rather
   than typical population examples. Mean and zero controls remain visible;
   no strong-effect history or retained model has been selected for its size
   trend. Fixed local statistics still permit circumference-dependent total
   input and whole-ring input modes.

## Execution and numerical assessment

Independent checks passed 240 pilot output hashes, 30 copied-input checks and
22 frozen-source hashes. All 30 cases reached 4000; the minimum concentration
was positive and the worst pool-relative error was $2.61\times10^{-12}$.
Direct cross-products and both endpoint energies at four lags in all three
windows of every case matched exported lag parts within $2.69\times10^{-14}$.
All six numerical workers and their output hashes also passed independent
checks.

| Audit at C60, higher input | Wave-pinning | Goryachev | Holmes 3 |
| --- | ---: | ---: | ---: |
| Maximum $|\Delta P|$, fixed forcing, tighter ODE solver | $4.53\times10^{-8}$ | $1.40\times10^{-9}$ | $7.70\times10^{-8}$ |
| Maximum $|\Delta P|$, coupled finer interpolant | 0.00546 | 0.00176 | 0.01807 |
| Maximum full-window correlation change, finer interpolant | 0.000110 | 0.0000176 | 0.000563 |

The second check changes the driving interpolation and quantifies path
sensitivity, rather than ODE error. Local maxima are more sensitive than the
global vector: finer interpolation changes Goryachev/Holmes prominent-peak
counts at 13.6%/5.9% of scalar samples, while half-contrast region counts agree.
Keep those thresholded peak counts descriptive. These three audited histories
support the declared five-unit interpolant; they establish no blanket
convergence guarantee for every future transition.

## Main freeze and optional improvements

The single bounded extension to horizon 32000 and maximum lag 8000 is justified
by unresolved pilot decay. Freeze that choice before execution; retain 12
independent seeds at C20/40/60 for both strengths, plus matched mean and zero
controls through the same horizon (243 cases total). Declare the main scalar
cadence 5, saved-profile cadence 40 and physical diagnostics from scalar
descriptors independently checked against saved profiles; the pilot physical
diagnostic retains its original 20-unit cadence. Longer observation may still
leave decay unresolved. No new amplitude search or model exclusion is needed.

Optional: a labelled absolute-scale close-up of the Goryachev C20 profiles
could make its remaining local activity easier to see beside the shared-range
kymographs. A centered vector covariance could later supplement adequately
resolved windows. Neither is a prerequisite for this bounded ensemble.

**Anticipated user critique:** “Does the larger ring preserve its local activity
site, or does its angular vector change less because the same physical motion
covers a smaller angle?” The actual Wave-pinning example demonstrates why the
distinction matters; the main experiment's vector, magnitude, local-pattern,
physical-displacement and coverage readouts can assess it without assuming a
cell velocity law.
