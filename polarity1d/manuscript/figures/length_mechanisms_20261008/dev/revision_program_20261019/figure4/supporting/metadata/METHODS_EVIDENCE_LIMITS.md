# Figure 4: methods, evidence and limits

This is a read-only reanalysis of the original 243 completed runs: 216 noisy
histories (3 models × 3 circumferences × 2 strengths × 12 independent seeds),
18 deterministic mean-input controls and nine no-input controls. All runs span
0–32000 native time units. No new simulation, burn-in selection, fitted
exponential time, event classifier or model parameter change is introduced.
Every native trace hash is checked against both producer status and the
previous audited seed table. Scalar cadence is 5; profile cadence is 40.

For circumference $C$, the native field and normalized activity moment are

$$
M_k(t)=\int_0^C u(s,t)e^{2\pi i k s/C}\,\mathrm ds,
\qquad P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|.\tag{1}
$$

The original direction table supplies B, retaining its frozen pooled-endpoint
estimator and 1000 whole-seed bootstrap draws. The full-vector raw correlation,
first/last-half windows and signed restricted areas remain in the unchanged
source report. New metrics use 2000 complete-seed bootstrap resamples, RNG seed
2026101904, and percentile 2.5–97.5 bounds. Time origins within a seed are never
independent resampling units; no confidence bands are supplied for one-run controls.

For the original resolved-site gate, both endpoints must have one current
half-contrast region, contrast ≥0.05 native units, at most one ≥0.05-prominence
peak and R≥0.1. A distinct sensitivity table removes only the peak-count veto.
Its less restrictive estimand does not repair or replace the primary gate.
The original full-profile independent descriptor audit remains authoritative;
all profiles displayed in A additionally reproduce their stored raw moments
within the numerical tolerance recorded in validation.json.

For shortest physical endpoint displacement, native site-angle memory and the
common-reference projection, use

$$
d_C(t,\tau)=\left[(s(t+\tau)-s(t)+C/2)\bmod C\right]-C/2,\tag{2}
$$

$$
K_C(\tau)=\left\langle\cos\left(\frac{2\pi d_C}{C}\right)\right\rangle_E,
\qquad
K_{60}(\tau)=\left\langle\cos\left(\frac{2\pi d_C}{60}\right)\right\rangle_E.\tag{3}
$$

The expectation pools the eligible endpoint pairs E across complete seeds.
Bootstrap resamples recompute the ratio of pooled seed numerator sums and
eligible-pair counts. Every native and remapped value uses identical E. The
reference circumference 60 is the largest tested, so every observed shortest
displacement maps to an angle at most π. This compares coordinates using
observed displacement distributions; it changes neither forcing nor dynamics.
It is not an additive causal attribution, a counterfactual evolved trajectory,
a diffusion estimate or evidence of uninterrupted site lineage. Original
shortest endpoint displacements already fold any travel beyond half the native
ring; remapping cannot recover winding or path length. Lags 5–8000, including
500, are exported so short-lag sensitivity to that folding is inspectable. Conditioning
on sparse single-region intervals can alter comparisons across sizes. CSVs
retain all sparse/undefined values, eligible-pair fractions, independent seeds
with pairs and bootstrap-defined fractions. Supplementary sensitivity curves
are shown only where at least six seeds contribute and coverage is ≥1%; this
display rule does not censor the tables. Goryachev's main physical curve remains
unavailable because the original peak gate is too restrictive.

Signed activity-direction increments use

$$
\Delta_\tau\phi(t)=\arg\left[P(t+\tau)P(t)^*\right],\qquad
\min(R(t),R(t+\tau))\ge r_0.\tag{4}
$$

Tables retain lags 5,40,100,500,2000,8000 and gates r₀=0.05,0.1,0.2. The
supplementary 5-unit survival curves average the pooled observed increment
distribution, with whole-seed uncertainty. Quantiles, maxima and the descriptive
count above π/4 are retained per seed. That threshold is a tail diagnostic,
not an abrupt-reorientation event definition. Small sampled increments do not
rule out subcadence fast changes; large wrapped increments alone would not
identify a physical site switch. The full native histories and R controls
must accompany interpretation.

Multiple-region and weak-R fractions use all 6401 sample times, with no chosen
adaptation exclusion. Half-contrast regions are defined relative to each
current profile's own minimum and range. Multiple regions need not be equal,
opposed, persistent, independent or functional competitors. Weak global
asymmetry can coexist with local activity. R is initial-moment normalized and
may exceed one; it is not a concentration, coherence ratio or migration speed.

The noise law is a declared model input; biological noise and migration are not
calibrated. These bounded histories do not establish stationary persistence,
universal enlargement effects or a force-to-motion relationship. The native
source campaign's fixed-forcing ODE and interpolation sensitivity checks remain
linked in its report; they are not long-horizon continuum convergence proofs.
