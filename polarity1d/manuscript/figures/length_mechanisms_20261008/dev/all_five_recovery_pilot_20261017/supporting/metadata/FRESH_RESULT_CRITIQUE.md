# Fresh actual-result critique: all-five linear recovery pilot

Reviewed 17 October 2026 from the [current report](../../RESULTS.md), the actual [PNG](../../all_five_recovery_pilot_curves.png), native sources, frozen producer, and saved full-field archives. This review ran no integrations and changed no source data.

**Verdict.** Accept this bounded pilot as evidence that the declared incumbent-targeted activity-loss assay supplies a resolved, amplitude-consistent early **local-core recovery** readout in all five retained linear models. Accept the recorded clocks as pilot decisions for subsequent review. This is not yet the all-length Figure 1 C result, a cross-model rate ranking, or evidence of general stability or challenger resistance. No remaining essential correction blocks the pilot checkpoint.

The useful scientific result is heterogeneous size dependence, rather than a large common enhancement. All 15 selected core deficits shrink, and actual probe core amounts rise. The displayed middle-to-long changes are:

| Model | Middle → long | Core deficit restored at the selected clock | Increase in the finite-window coefficient |
|---|---|---|---:|
| Wave-pinning | 19 → 39 | 13.39% → 13.73% | 2.56% |
| Goryachev | 19 → 39 | 5.40% → 18.59% | 243.94% |
| deBelly | 27 → 39 | 18.80% → 19.27% | 2.50% |
| Holmes 3 | 28 → 50 | 19.05% → 19.11% | 0.35% |
| Spring | L₀ 1.6 → 2.6 | 12.57% → 14.83% | 17.98% |

These are the actual pilot pairs, not interchangeable versions of Figure 1 A's highlighted pairs. The shortest Goryachev case restores only 0.0603% at its selected clock; its visually near-flat trace is resolved rather than missing.

## Independent evidence checks

- All 31 frozen-source hashes, all 15 copied and original checkpoint hashes, and all 10 report artifact hashes verify. Updated report source hashes also verify. Frozen native core sources match current native sources. The canonical harness differs from the producer only in metadata serialization/collection handling; the frozen integration implementation remains preserved.
- All **150** archived trajectories contain finite fields at every prescribed time. Initial kicks match the specified restart, conservative inactive-pool compensation and untouched other fields. Spring's maximum initial-output discrepancy is one floating-point ulp, 2.22 × 10⁻¹⁶. All 144 fieldwise comparisons of independently executed sham and zero-kick histories are identical. Independently recomputed core/whole amounts agree with archives within 3.56 × 10⁻¹⁵; Spring native RHS rates reconstructed at restart, selected and final times across its 30 trajectories agree exactly with saved rates.
- All **315 case × candidate** gates independently reproduce from the NPZ arrays, including numerical sensitivity, all three amplitude-pair comparisons and every sampled prefix. The recorded clock is the largest common passing candidate in each model. The next candidate fails saturation in at least one size; no candidate exclusion here comes from numerical or amplitude-linearity failure. At selected clocks, the maximum pairwise amplitude response spread is approximately 1.29%, below the prescribed 5% allowance. Numerical sensitivity is empirical, not a certified discretization bound.
- Spring completion is supported by raw evidence, not merely a green collector flag. Its original JSON files still say `producer=running`, `complete=false`; the fatal error logs show failure at final JSON serialization after the archive save. Separately hashed recovered statuses retain those originals and logs, and identify unavailable solver counters. The [recovery script](../../../../../../code/mechanisms/all_five_recovery_metadata_repair_20261017.py) reconstructs existing arrays and rates without integration. Archive, original-status and error-log hashes verify. The completion claim is scientifically defensible; the original jobs did not exit successfully.

## Interpretation and actual figure

The perturbation is a common question with declared native implementation differences. WP, Goryachev, Holmes and Spring transfer removed active material to their corresponding explicit inactive pool. deBelly has no such Rac pool: Rac alone is lowered while Rho, MCA, mechanics and four endpoint-controller values start unchanged. Its native 0.5-s staggered update and 0.98 endpoint relaxation remain active. The dense reference checks arithmetic of that update, not temporal convergence. Spring retains endpoint motion, dilution and each trajectory's own physical-length integration. The source implementations support these distinctions; the assay does not isolate a reaction pathway or attribute recovery causally to width.

The coefficient measures the deficit **inside the frozen original core**, including reaction and transport. The exported outside-core difference is scientifically material in some cases. At Spring L₀ 0.6, core recovery is 4.65%, whereas whole-domain active-deficit recovery, normalized to its own initial deficit, is 9.36%; the outside-core probe excess is 4.72% of initial core damage. At deBelly L39, the corresponding recovery fractions are 19.27% and 20.94%. Thus local-core recovery must not be relabelled whole-field restoration. The report already makes the distinction and exports the evidence; adding one quantitative example would make it easier to understand.

The shared 0–0.2 restored-fraction axis is honest. It correctly leaves the WP, deBelly and Holmes middle/long traces nearly coincident instead of enlarging tiny differences. Native clocks and all three sizes are visible, and the full short histories justify the early-window restriction. The ten-panel figure works as a supporting pilot diagnostic. At manuscript-page width, its five columns, repeated labels and long footer will need simplification. Retain all-five coverage in the eventual main panel; move full recovery histories to support rather than selecting arbitrary successful model examples.

Goryachev L39 has four sample centers in the interpolated core, versus eight at L19. Its strong increase is secure for the retained native discretization and numerically resolved under tighter solver tolerances; those checks do not establish spatial convergence. Likewise, maximum 0.5% attenuation of tapered basal excess is not equal integrated fractional core damage across models: primary damage ranges from approximately 0.052% to 0.285% of basal core active amount.

## Essential revisions and scope gates

**Resolved during this review:** the earlier phrase “faster later recovery” was unsupported by log-time steepening. Archived normalized paired core rates decrease from restart to the selected clock in all 12 continuous cases, for example Goryachev L39: 40.04 → 27.19 per native time, and Spring L₀ 0.6: 6.99 → 1.75. The current report correctly states continued deficit shrinkage and the later violation of the 20% range, explicitly excluding an acceleration inference. Regenerated report/generator hashes verify; the scientific curves and gates are unchanged.

**Before a main-text claim:** preserve the local-core definition, native clocks, deBelly intervention distinction and spatial-resolution limitation. Obtain and audit the authorized all-length result before describing an all-length trend; three sampled sizes cannot determine the intervening curve. A ring extension requires its declared size-range/challenged-state clock check and the same perturbation definition. Existing inactive-pool-return maps cannot silently fill the new recovery panel. These are limits on the next claim, not reasons to rerun or enlarge this completed pilot now.

## Optional improvements

- Put one of the quantitative core/whole comparisons above into the report and define the plotted restored fraction directly in the eventual caption.
- Label the nearly zero shortest Goryachev endpoint numerically if that trace remains in a main figure. Keep the common vertical scale.
- Use fewer panels for the main figure and retain these complete all-five histories as linked supporting evidence. Treat source/pathway attribution or mesh sensitivity as separately declared work if stronger mechanistic or continuum claims are later sought.

## Anticipated user critique

This is a prediction from the recorded user preferences, not a quotation. The likely questions are whether the figure states a clear comparative result, retains all five models, uses honest scales, and ties observations to displayed evidence. The defensible answer here is that local early recovery is measurable in every model, with a much larger middle-to-long change in Goryachev, a moderate change in Spring, and small changes in the other three. The user may reasonably want the full size curves before accepting Figure 1 C, rather than another example-based summary. Figure 1 should keep its linear organization/input-gain/recovery sequence; subsequent supported-three-model ring size/separation results should use the common assay and retain their separate geometry and evidence. None of this pilot establishes a migration-persistence mechanism.
