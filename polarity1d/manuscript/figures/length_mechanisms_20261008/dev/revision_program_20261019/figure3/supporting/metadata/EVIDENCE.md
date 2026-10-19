# Figure 3: claims and retained evidence

Generated: 2026-10-19T15:28:39.640809+00:00. Reanalysis finished; no simulations were commissioned.

The geometry claim is exact: circular shorter-arc separation d maps to outward-normal angle Δθ=2πd/C. In the competition grid d is the prescribed distance from the basal incumbent peak to the challenger cue centre, and the angle is between those prescribed positions. Neither is an observed separation between evolved sites. This adds a continuously valued two-dimensional directional consequence to a spatial polarity pattern; the established interval experiments do not supply these circular normal directions. This is a coordinate consequence, not an extra biochemical mechanism or a cell-movement prediction.

The controlled competition claim is narrower: in Holmes C40, changing only challenger cue-centre position from d12 to d16 changes the finite endpoint from one spanning region to two associated peaks. Original basal, local parameters, input peak/support, input duration and readout time agree. The C30,d12 comparison additionally changes basal organization and cannot isolate geometry from chemistry. Same-distance and same-angle contrasts show that neither one scalar coordinate completely describes the grid.

The circular-vector consequence is supported by the full native C40,d20 profile and independently replayed moments: two antipodal activity peaks have normalized first-moment magnitude 0.03993 while q₂ is 0.40065. This is cancellation of a spatial moment. The direction is unavailable by the prespecified R=0.1 gate, and the local field remains active.

## Descriptor analysis, including useful negative results

All 652 complete conditions are retained from 867 requested placements: 198 overlap exclusions and 17 unqualified preparations. Exclusions are not missing outcomes or zero activity. Within any fixed circumference d and Δθ are exact positive multiples, hence their ranks and one-dimensional ordering are identical. No within-size regression can distinguish their independent causal influence. The grid contains 17 circumferences in Wave-pinning/Goryachev and 16 with qualified, nonoverlapping Holmes outcomes.

For a descriptive transfer check, fit the same fixed 7-neighbour classifier separately within each model while withholding every placement at one entire circumference. Distances use only the training set's coordinate standard deviations. All tied neighbours at the kth distance enter; vote ties use the training class frequency and then stable class order. No tuning or random split is used. Accuracy is averaged equally over held-out sizes; condition-weighted and balanced-class accuracies, fold predictions, and k3/k15 sensitivity are exported. Deterministic placements are not independent biological samples and these accuracies carry no sampling confidence interval or significance test.

- Wave-pinning: placement changes the finite class in 9/17 qualified sizes; held-out-size accuracy is 0.493 for angle, 0.786 for physical separation and 0.923 for physical separation plus basal width.
- Goryachev: placement changes the finite class in 3/17 qualified sizes; held-out-size accuracy is 0.664 for angle, 0.792 for physical separation and 0.952 for physical separation plus basal width.
- Holmes 3: placement changes the finite class in 11/16 qualified sizes; held-out-size accuracy is 0.546 for angle, 0.752 for physical separation and 0.934 for physical separation plus basal width.

Physical separation transfers better than angle in this fixed analysis, but adding incumbent width improves transfer further. This is not causal superiority of physical distance: circumference covaries with available pool, incumbent organization and the unmeasured interactions of route lengths. The equally accessible physical-separation-plus-C control performs similarly; width is not isolated by prediction. Finite class rarity and excluded placements affect transfer, and the classifier does not establish asymptotic stability. No ring-vs-interval matched intervention identifies which competition mechanisms are uniquely circular. Therefore the classifier remains supplemental and the main figure focuses on exact geometry, controlled position contrasts and the activity-moment consequence.

## Definitions

For a uniform native ring with activity a(s), define its angular moments by

$$ M_n(t)=\int_0^C a(s,t)\exp\!\left(\mathrm{i}n\frac{2\pi s}{C}\right)\,\mathrm{d}s. \tag{1} $$

Use P(t)=M₁(t)/M₁(0), R(t)=|P(t)| and qₙ(t)=|Mₙ(t)|/M₀(t). The angle is reported only where R≥0.1, following the reviewed existing analysis. The native finite-volume quadrature is retained, with the original peak cell as angular origin. Uniform background contributes zero first moment. Neither a global angle nor a local maximum identifies cell movement.

## Reproduction and retained sources

Run `OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/revision_figure3_20261019.py` from repository root. This reads audited existing simulation and moment tables, replays four selected full native histories, validates exact phase continuity and four time-point moment summaries, exports complete placement and held-out-size tables and renders main/support sheets. Work is serial, below 1 GB resident memory and requires no new integrations.

The prior gain/recovery and standardized finite-cue ring sheets remain applicable supporting evidence, not replaced measurements. Sources include `dev/vision_aligned_checkpoint_20261017/vision_aligned_figure2_response.png`, `vision_aligned_figure2_competition.png`, and the common recovery tables under `data/runs/vision_ring_recovery_20261017/results`. These support local response measurements while the revised main figure concentrates on geometry.

See `validation.json` for numerical checks, measured allocation and artifact hashes, and `source_sha256.json` for immutable input hashes.

## Subsequent four-model distance-duration revision

This geometry sheet is the first of two main sheets. The [four-model continuation](../../distance_response/revision_figure3_distance_duration.png) adds controlled finite-duration responses and prominent native Spring profile evidence. Its [evidence and limits](../../distance_response/supporting/metadata/DURATION_EVIDENCE.md), [150 requested conditions with actual archive provenance](../../distance_response/supporting/data/distance_duration_all_models.csv) and [endpoint sensitivity table](../../distance_response/supporting/data/distance_duration_shortest_tested_sensitivity.csv) distinguish retained angle, vector magnitude and fixed-region activity loss. The [complete Spring size/placement sheet](../../spring_distance/spring_distance_duration.png) and [Spring field/geometry audit](../../spring_distance/supporting/metadata/validation.json) retain the full three-size extension. The first sheet's native competition illustration is explicitly Holmes-specific; no all-model comparison is inferred from those selected histories.

## Scope of historical chemical support and the current Spring counterpart

The retained 652-condition placement grid, held-out-size descriptor fits and exact interval/ring reflection maps cover Wave-pinning, Goryachev and Holmes only. They are explicitly chemical-only historical support. The [Spring standard-cue coverage record](../../spring_distance/supporting/metadata/descriptor_coverage.json) and [14-row table](../../spring_distance/supporting/data/spring_standard_long_cue_descriptor_coverage.csv) evaluate the same finite-duration question at all three qualified Spring sizes. All 14 eligible T=1000 geometries yield one cue-associated peak; continuous angle and old-core amount variation remain exported, but categorical descriptor ranking is uninformative for a constant class. No generic classifier accuracy is reported.

Exact nonuniform interval reflection is not applicable to the declared Spring mean-force closure: the native interval samples endpoint Rac forces, whereas the ring averages the nonlinear force around the perimeter. The validated uniform-state derivative and native stretching/dilution/diffusion carryover do not equate these nonuniform mechanical rules. The current Spring counterparts are its three-size basal organization/function measurements and its actual distance-duration trajectories, in explicitly native units.
