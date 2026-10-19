# Spring active amount and spatial first harmonic

This adds no integration or forcing change. All 81 original cases are retained:72 noisy histories (three sizes × two strengths × 12 seeds), six deterministic mean-input controls and three no-added-input controls. The scalar trace stores active-only M₀ separately from the conserved active-plus-inactive `amount`.

$$
M_0(t)=C(t)\,\operatorname{mean}_j u_j(t),\qquad
M_1(t)=C(t)\,\operatorname{mean}_j[u_j(t)e^{2\pi i q_j}],
\qquad q_1(t)=\frac{|M_1(t)|}{M_0(t)}.
\tag{1}
$$

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad
R(t)=|P(t)|=\frac{M_0(t)}{M_0(0)}\frac{q_1(t)}{q_1(0)}.
\tag{2}
$$

Every scalar record is retained at its actual cadence 0.5 over 0–3200 native units (6401 samples), without an R gate or adaptation exclusion. Every saved active/inactive profile is replayed at its actual cadence 4 (801 profiles per case). Physical M₀ and M₁ are reconstructed with instantaneous circumference; no profile interpolation is introduced. These cadence 4 checks support the stored scalar readouts, rather than claim an independent full-field reconstruction at every cadence 0.5 frame. Scalar q₁/P definitions and the paired factorization are checked independently. The exported conserved total is separately labelled and is never substituted for active M₀.

`case_index.csv` identifies all 81 native trace/profile hashes, exact initial physical quantities and exported case NPZ hashes. Each case NPZ contains physical M₀/M₁, q₁, M₀/M₀(0), q₁/q₁(0), R, the paired product/residual, C(t) and separately labelled conserved total. `ensemble_scalar_traces.csv.gz` supplies every scalar time point for all 15 conditions. Controls have one deterministic trajectory and no confidence intervals. Every noisy condition uses all 12 independent seeds.

Intervals use 2000 whole-seed multinomial bootstrap bags (RNG seed 2026101911). The same bags are applied jointly to every component and time; time samples are not independent replicates. Full/first-half/last-half windows are closed 0–3200, 0–1600 and1600–3200, with 6401/3201/3201 samples. The midpoint belongs to both half-windows; they are paired sensitivities rather than independent replicates. Their summaries are arithmetic sample means, not stationarity estimates or new correlation fits. Case minima/maxima are explicitly sampled extrema.

Mean R equals the mean of paired products. It generally differs from the product of separate component means. Both exports retain that latter product, the appropriate component covariance and the factorization residual; bootstrap window covariance is recomputed from the same jointly resampled cases. Covariance here is the exact mean-product difference (population-moment convention), without a Bessel correction. Cross-case covariance at a fixed time is distinguished from covariance over the joint case/time population within a window. No additive causal shares are inferred.

q₁ is the normalized first harmonic of a nonnegative activity field, not a proof of profile identity, a site-lineage label or a general coherence measure. R may exceed 1. Amount and first-harmonic changes are measured here rather than inferred from mean R alone. Native Spring units/forcing differ from the chemical campaign; their common M₀/P/q₁ definitions permit a methods comparison, not a pooled exposure, timescale or model ranking. The chemical exports use cadence 5 and a 32000-unit horizon; Spring uses cadence 0.5 and 3200. The supporting sheet changes no current main layout.

## Measured full-window means

All entries use the same 6401 scalar samples per case and all 12 noisy seeds. Brackets are pointwise 95% whole-seed intervals for the full-window sample mean. R and the paired product agree to roundoff.

| f | Resting C₀ | M₀/M₀(0) | q₁/q₁(0) | R |
|---|---:|---:|---:|---:|
| 0.001 | 1 | 1.000306 [1.000293, 1.000319] | 0.997416 [0.997167, 0.997710] | 0.997720 [0.997480, 0.998005] |
| 0.01 | 1 | 1.002852 [1.002759, 1.002948] | 0.997086 [0.995360, 0.998975] | 0.999917 [0.998226, 1.001750] |
| 0.001 | 2 | 1.000254 [1.000238, 1.000271] | 0.997595 [0.997463, 0.997714] | 0.997849 [0.997722, 0.997960] |
| 0.01 | 2 | 1.003225 [1.003047, 1.003381] | 0.961057 [0.956222, 0.965233] | 0.964002 [0.959225, 0.968133] |
| 0.001 | 3 | 1.000263 [1.000247, 1.000279] | 0.997760 [0.997682, 0.997851] | 0.998022 [0.997947, 0.998115] |
| 0.01 | 3 | 1.003498 [1.003310, 1.003695] | 0.904576 [0.890662, 0.918113] | 0.907436 [0.893544, 0.920956] |

For stronger input at C₀3, the measured full-window mean active-amount ratio is 1.00350, while the first-harmonic ratio is about 0.90458 and mean R is 0.90744. The lower mean relative first moment therefore accompanies lower normalized first harmonic rather than reduced mean active amount in this window. This does not assert constant amount at each frame, rigid translation, a stationary distribution or a causal decomposition of the size response. All deterministic controls and half-window values remain in the tables.
