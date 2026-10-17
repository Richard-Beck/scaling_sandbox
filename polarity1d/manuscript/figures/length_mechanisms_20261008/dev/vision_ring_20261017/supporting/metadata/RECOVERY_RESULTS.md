# Common active-deactivation recovery across ring size and cue separation

The new common assay is complete for **50 qualified basal states and 652 exact challenged states** in Wave-pinning, Goryachev and Holmes. All 702 state assays pass the frozen numerical gates, at all three amplitudes. The old inactive-withdrawal/return assay is neither substituted nor renamed.

Incumbent recovery after a small targeted active-to-inactive transfer increases with circumference in all three models, with a substantially larger change in Goryachev over the displayed near-doubling pair. During challenge, this conditional response is lower than the same size's basal response in every tested state. It stays positive while the actual fixed incumbent region loses activity in all 652 unperturbed challenged trajectories during the assay window.

The coefficient therefore contributes a measurable local response to the size/separation framework; its sign alone does not predict whether a challenger wins, two sites remain, or activity merges. Those later results require the linked native histories and complete outcome matrix.

## New basal size response

The comparison below is entirely within native rings. C19 and C39 are circumferences, not the same-sized linear domains. Exact C = 2L counterparts are checked separately. Rates across models use distinct native clocks and are not a biochemical calibration.

| Model | K at C19 | K at C39 | C39/C19 | Qualified basal sizes | Challenged conditions |
| --- | ---: | ---: | ---: | ---: | ---: |
| wave_pinning | 0.65679738 | 0.75454239 | 1.14882 | 17 | 210 |
| goryachev | 2.3031468 | 10.06498 | 4.3701 | 17 | 280 |
| holmes_model3 | 1.0281111 | 1.0640923 | 1.035 | 16 | 162 |

## Same-assay native interval/ring correspondence

Thirteen available native basal pairs have circumference C = 2L and exactly the same active-deactivation definition and native observation clock. Maximum absolute relative K difference is **0.0477% for WP, 1.0164% for Goryachev and 0.000225% for Holmes**. These are independently prepared states at their retained native meshes: linear 201 cells versus ring spacing 0.05. Their centers and fractional core/taper geometry are independently measured. The largest difference is Goryachev L35/C70. This supports close baseline correspondence without asserting exact state identity or attributing the mismatch solely to integration error. The smaller empirical E checks integration sensitivity at each fixed discretization.

All 13 native linear source archives and checkpoint hashes were independently verified, and their F/K values were recomputed at the exact shared clock. The [pair table](../../../../../../data/runs/vision_ring_recovery_20261017/results/ring_linear_common_recovery_pairs.csv) and [13-pair audit](../../../../../../data/runs/vision_ring_recovery_20261017/results/ring_linear_common_recovery_pairs_audit.json) retain the exact sources. The separate matched-mesh operator fixture remains an algebraic identity check.

## Challenge-conditioned response

The challenged assay begins at exactly 10 native units after cue onset, with the same continuing cue applied to probe and sham. The original basal geometry and removal direction remain frozen. The scalar depends on both circumference and the complete state reached at that separation; basal recovery is measured only once per size.

| Model | Smallest K10/Kbasal − 1 | Largest K10/Kbasal − 1 | Positive conditional F | Actual sham core losses |
| --- | ---: | ---: | ---: | ---: |
| wave_pinning | -0.531107 | -0.000164174 | 210 | 210 |
| goryachev | -0.202047 | -4.32491e-05 | 280 | 280 |
| holmes_model3 | -0.0785095 | -0.00136672 | 162 | 162 |

All primary conditional repair fractions are positive and resolved above their empirical numerical sensitivity. Model-specific native windows are fixed across interval/ring geometry, size, cue location and basal/challenged stage: WP/Holmes 0.1778279410038923; Goryachev 0.005623413251903491. Native K, dimensionless F, complete deficits, actual sham/probe changes and RHS rates are retained.

## Assay and verification

The frozen original basal excess-active profile is multiplied by a compact raised-cosine taper within its half-contrast core, using independently measured left/right radii. The resulting direction h0 is removed from active species0 and added pointwise to inactive species1 at ε = 0.001, 0.005, 0.01. Other native fields are initially unchanged, and every field evolves during recovery. There is no later return, clipping, remeshing, retuning or state-dependent direction. The core integration uses fractional finite-volume edge-cell weights.

$$
D(t)=\int_{\Omega_0}\left[a_\mathrm{s}(x,t)-a_\mathrm{p}(x,t)\right]\,\mathrm{d}x,\qquad F_\tau=1-\frac{D(\tau)}{D(0)},\qquad K_\tau=\frac{F_\tau}{\tau}.
\tag{1}
$$

The 15-state preflight was reviewed before the full grid. All sampled positive-time prefixes and amplitudes pass the same tight/tighter numerical-sensitivity, saturation and amplitude-linearity gates. Full-grid maximum E = 5.5452026e-07; maximum |F|+E = 0.19121824. No state required an unavailable kick or another clock. Independent sham/zero-kick integrations agree bitwise. The empirical E is not a certified error bound, and K is neither an instantaneous eigenvalue nor generic stability.

The read-only data audit checked **7,020 native trajectory trials** and all 2,106 positive-amplitude summary rows. It independently reconstructed kicks, endpoint amounts and native RHS rates, conservation checks, F/K/E and prefix gates. Maximum conserved-pool drift was 6.661e-16; maximum native endpoint rate-allocation discrepancy was 1.066e-14. This is an equation/data audit rather than an independent re-integration. The fresh-context design reviewer separately verified actual preflight archives and gates.

Every challenged starting array agrees bitwise with its complete retained native cue history at time10; all 702 copied inputs have exact source hashes. Matched 201-cell linear/402-cell ring operator, core-weight, taper and cue-width reflection fixtures pass. Those fixtures do not imply pointwise equality between independently prepared native basals.

Full geometry bookkeeping remains **867 planned cells: 652 challenges, 198 pre-existing cue/core overlap exclusions, 17 unqualified Holmes C80 cells**. Source and executable inputs were frozen before SLURM integration; arrays 21126192 (preflight) and 21126207 (full grouped grid) completed. No new jobs remain.

## Data and provenance

- [basal.csv](../../../../../../data/runs/vision_ring_recovery_20261017/results/basal.csv) and [challenged.csv](../../../../../../data/runs/vision_ring_recovery_20261017/results/challenged.csv): all three positive amplitudes; primary ε = 0.005.
- [challenged_full_matrix.csv](../../../../../../data/runs/vision_ring_recovery_20261017/results/challenged_full_matrix.csv): all 867 planned geometry cells, with unavailable statuses retained.
- [recovery_audit.json](../../../../../../data/runs/vision_ring_recovery_20261017/results/recovery_audit.json), [independent_data_audit.json](../../../../../../data/runs/vision_ring_recovery_20261017/results/independent_data_audit.json) and [scientific_summary.json](../../../../../../data/runs/vision_ring_recovery_20261017/results/scientific_summary.json).
- [Complete source-state audit](RECOVERY_INPUT_AUDIT.json), [algebraic mirror checks](RECOVERY_ALGEBRAIC_PREFLIGHT.json) and [bounded preflight result](RECOVERY_PILOT.md).

The full state endpoints, five declared primary/sham snapshots and dense scalar deficit/amount/rate traces for every trial remain under `data/runs/vision_ring_recovery_20261017/workers/`. Allocated storage measured after collection: **794,918,912 bytes (0.7403 GiB)**, before this Markdown record. Root is responsible for updating the related output-group row in DATA_MANIFEST.md.
