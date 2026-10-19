# Figures 3–4: shared readouts and passive signatures

Discussion note, 2026-10-19. Existing sheets and audited statistics; no new layout, simulations or event definitions.

## Shared definition and suggested transition

Both experiments use active concentration, physical segment length and outward normal to form the same physical moments:

$$
M_0(t)=\sum_j u_j(t)\,\Delta s_j(t),\qquad
M_1(t)=\sum_j u_j(t)\,\Delta s_j(t)e^{i\theta_j},\qquad
\Delta s_j(t)=\frac{C(t)}{N}.
\tag{1}
$$

$$
P(t)=\frac{M_1(t)}{M_1(0)},\qquad
R(t)=|P(t)|,\qquad
q_1(t)=\frac{|M_1(t)|}{M_0(t)}.
\tag{2}
$$

Division by the initial complex moment rotates the reference axis once and scales by the initial magnitude. It does not normalize by current total activity. In particular,

$$
R(t)=\frac{M_0(t)}{M_0(0)}\,
      \frac{q_1(t)}{q_1(0)}.
\tag{3}
$$

M₀ is active amount, distinct from the conserved active-plus-inactive pool. Thus R combines amount and spatial balance; R≈1 does not prove a rigidly translated profile, and small R does not by itself prove activity loss. Spring uses instantaneous physical circumference in both moments. The [Figure 3 concept definitions](../figure3/standard_angles/vector_concept/README.md) and [Figure 4 Spring methods](SPRING_METHODS_EVIDENCE_LIMITS.md) agree.

Suggested short transition:

> Figure 3 used the normalized physical activity vector P and total active amount M₀ to distinguish cue-driven deflection, retained orientation and cancellation. Figure 4 applies the same definitions under ongoing spatially varying activation, asking how directional memory relates to resolved-site displacement and local organization. Native histories and matched mean-input controls qualify those comparisons.

These are different experiments. Figure 3 tests stationary localized cues and withdrawal; Figure 4 tests ongoing stochastic activation with mean-input and no-added-input controls. Chemical sizes are C20/30/40 in the current Figure 3 experiment and C20/40/60 in Figure 4. Common notation or overlapping sizes do not establish identical preparations, matched exposure or calibrated native clocks.

## Panel → measurement → interpretation

| Existing panels | Measurement | Supported reading and immediate limit |
|---|---|---|
| A1–A3 | Native active-field histories from fixed medoids | Spatial organization and change; one representative does not establish population frequencies or rigid translation. |
| B1–B3 | Global activity-direction cosine, conditional on both endpoint R≥0.1; both strengths | Resultant-direction retention among eligible pairs, distinct from site-centroid memory, full-vector correlation, amount and movement. |
| C1–C3; D1–D3 | Strictly gated shortest-endpoint RMS and native/remapped site-angle cosine on identical pairs | Physical versus angular coordinates of resolved shifts, with pair coverage. Neither uninterrupted site ancestry, traveled distance nor reconstruction/causal decomposition of B. Chemical lag 8,000; Spring C uses 800 and D uses 4. |
| E1–E3 | Per-seed fractions with multiple current half-contrast regions or R<0.1 | Local multiplicity and weak resultants; all 6,401 samples enter. Multiplicity does not establish functional competition. |
| Linked controls/windows/increments | Magnitude, raw-vector memory, mean/no-added-input controls, half-windows and sampled increment tails | Forced-state qualification and finite-window/cadence sensitivity; no stationary persistence or event classifier. M₀ remains separate. |

## Two useful observable pairs for the next discussion

**1. Resolved-site shift paired with angular remapping.** At stronger input and lag 8,000, Wave-pinning C20 site cosine changes from 0.593 to 0.944 when identical eligible shifts are expressed on C*=60; coverage is 99.2%. C60 has native/remapped cosine 0.972. Physical RMS also decreases 3.25→2.25: both shifts and coordinates change. Holmes' rising angular retention accompanies nonmonotone physical RMS 4.57/6.32/5.54 across C20/40/60, with coverage 92.1%/66.4%/37.3%; only 10 seeds contribute at C60. These are conditional populations.

Spring's lag 4 small/middle basal-arc RMS 0.0607/0.0657 accompanies native site cosines 0.940/0.980 and common-reference cosines 0.9923/0.9910. At lag 800, RMS 0.298/0.592/0.863 approaches uniform wrapped-angle references 0.301/0.590/0.879: folded-endpoint mixing, not travel or recovered winding. Site-centroid and global first-moment direction remain different even under a one-region gate.

**2. Local-region occupancy paired with relative-resultant magnitude.** Stronger-input Wave-pinning is almost always one current region (multiple-region times 0.001% at C20; 0.056% at C60), with mean R near 1 and no sampled R<0.1 times. Holmes' multiple-region fraction rises 1.15% [0.27%, 2.27%]→17.15% [4.38%, 32.93%], while weak-R time reaches only 0.53% [0%, 1.31%]. Substantial resultants can coexist with multiple local regions. Figure 3's translated final-profile proof cannot be carried over to noisy Figure 4 histories.

Spring reinforces the separation: stronger-input lag 800 global-direction cosines are approximately 0.019/−0.017/0.028, while mean R remains 1.000/0.964/0.907. Direction can decorrelate while the relative first moment remains substantial. At C₀3, multiple-region time is 14.12% [11.90%, 16.39%], whereas weak-R time is 0.40% [0.25%, 0.57%]. Weaker input instead gives direction cosines 0.584/0.832/0.958 with mean R≈.998 and no multiple-region samples. This is forcing-dependent behavior, not a universal enlargement law. The quoted R values do not establish unchanged M₀.

Goryachev C20 has 50.51% [49.23%, 51.85%] multiple-region time and 11.05% [10.26%, 11.83%] weak-R time at stronger input, but its mean-input control also loses polarity. This is not clean wandering around an intact original state. Strict lag 8,000 site pairs are 1/0/0 across C20/40/60, precluding a population physical-site curve. Both strengths, controls and native examples remain accessible.

## Limits and source check

Fractions and intervals above are the existing whole-seed estimates, with 12 independent noisy seeds per condition; temporal samples are not independent replicates. Region multiplicity does not establish opposed/equal/persistent competitors. Stronger Spring has rare >π/4 sampled increments at cadence 0.5, whereas the chemical diagnostic has none at cadence 5, conditional on both endpoint R≥0.1. This is supporting native context, not a cross-model event-rate comparison or evidence for/against subcadence events.

Chemical M₀/M₀(0) is in the [original ensemble table](../../vision_noise_20261017/supporting/data/vision_noise_ensemble_traces.csv.gz). Spring's active-only scalar M₀, separate from conserved `amount`, is now [exported for all 81 cases](supporting/data/spring_amount_spatial_balance/case_index.csv) at cadence 0.5 and independently checked against every saved profile at cadence 4. The [amount/spatial-balance supplement](figure4_spring_amount_spatial_balance.png) measures the previously implicit distinction: stronger-input C₀3 has full-window mean M₀/M₀(0)=1.00350 [1.00331,1.00369], q₁/q₁(0)=0.90458 [0.89066,0.91811] and R=0.90744. Lower mean R here accompanies lower normalized first harmonic, not lower mean active amount. The paired product is exact per case/frame; separately averaged factors require their covariance term. This does not assert constant M₀, rigid translation or causal shares. Current main artwork is unchanged.

- [Existing two-sheet caption](supporting/metadata/CAPTION_WITH_SPRING.md) and [current quantitative evidence](RESULTS_WITH_SPRING.md).
- [Chemical observation fractions](supporting/data/observation_statistics.csv) and [site/remapping estimates](supporting/data/endpoint_statistics.csv).
- [Spring observation fractions](supporting/data/spring_noise_observation_statistics.csv) and [short/long-lag endpoints](supporting/data/spring_noise_short_long_endpoint_summary.csv).
- [Chemical-only interim methods](supporting/metadata/METHODS_EVIDENCE_LIMITS.md), read together with the [Spring extension](SPRING_METHODS_EVIDENCE_LIMITS.md).
- [Spring amount/balance definitions, measured means and limits](SPRING_AMOUNT_SPATIAL_BALANCE.md), [joint condition/window summaries](supporting/data/spring_amount_spatial_balance/condition_window_summaries.csv) and [scalar/profile audit](supporting/metadata/spring_amount_spatial_balance_validation.json).
- [Independent scientific review](../supporting/metadata/INDEPENDENT_REVIEW.md).

