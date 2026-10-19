# Figure 4: final four-model evidence and bounded claims

The two main sheets compare Wave-pinning, Holmes 3 and Spring in their separate native units. Complete Goryachev results and its forced-state controls remain linked supporting evidence. All 243 original chemical cases and all 81 new Spring cases are retained, including both declared strengths, every original replicate and matched mean/no-input controls.

## Two claims

1. Angular retention and physical endpoint displacement can respond differently to size. Wave-pinning's stronger-input angular retention rises while its site RMS decreases; Holmes' physical endpoint RMS is nonmonotone despite rising angular retention. Spring's weaker-input angular retention rises with size even though basal-reference arc RMS is nonmonotone. Stronger Spring input removes long-lag direction retention at all three sizes, so these results do not support a universal enlargement advantage.
2. Pairwise remapping and passive local-pattern/magnitude measurements add distinct information. Remapping identical eligible site endpoints changes their angular interpretation without changing their physical shifts. Multiple local regions, weak global resultants and sampled increment tails describe observed organization; they do not establish functional competition, an abrupt-event mechanism, site ancestry or cell movement.

## Spring quantitative comparison

All entries below use 12 independent seeds and lag 800 native units. Parentheses are pointwise 95% whole-seed bootstrap intervals. Site coverage is the fraction of candidate endpoint pairs passing the strict gate, not the fraction of time with proven site ancestry.

| f | Resting C₀ | Direction cosine | Basal-reference endpoint RMS | Eligible pairs | Mean R | Multiple-region time | Weak-R time |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.001 | 1 | 0.584 (0.473, 0.692) | 0.170 (0.142, 0.196) | 100.0% | 0.998 | 0.000% | 0.000% |
| 0.001 | 2 | 0.832 (0.768, 0.885) | 0.197 (0.160, 0.235) | 100.0% | 0.998 | 0.000% | 0.000% |
| 0.001 | 3 | 0.958 (0.946, 0.969) | 0.142 (0.121, 0.160) | 100.0% | 0.998 | 0.000% | 0.000% |
| 0.01 | 1 | 0.019 (-0.019, 0.057) | 0.298 (0.291, 0.304) | 99.8% | 1.000 | 0.048% | 0.003% |
| 0.01 | 2 | -0.017 (-0.074, 0.050) | 0.592 (0.565, 0.615) | 87.3% | 0.964 | 3.470% | 0.077% |
| 0.01 | 3 | 0.028 (-0.043, 0.099) | 0.863 (0.818, 0.908) | 57.7% | 0.907 | 14.119% | 0.405% |

At stronger input, lag-800 RMS values 0.298/0.592/0.863 approach the uniform wrapped-angle references 0.301/0.590/0.879. These folded endpoint measurements discard winding and cannot measure travel. In the main remapping panel Spring instead uses early lag 4: the small and middle rings have similar basal-arc RMS (0.0607/0.0657), native site cosines 0.940/0.980 and common-reference site cosines 0.9923/0.9910. The full curves and lag 4/50/800 summaries remain accessible. The CDF distance from uniform is descriptive only, with correlated temporal origins and no p-value.

The largest observed relative circumference excursion across all noisy trajectories is 0.119%. Basal-reference arc is the primary site metric; midpoint-current-circumference arc and circumference increments are retained as metric sensitivities. A fixed material site changes neither primary coordinate under uniform stretch.

## Sampled increment evidence

The original chemical diagnostic at cadence 5 contains no >π/4 increment conditional on both R≥.1; its maximum is 0.4074 rad. This does not exclude faster changes. Spring has rare large sampled global-angle increments at cadence .5. They are displayed with actual native field and scalar-vector context, selected as diagnostics rather than medoids or event classes.

| f | Resting C₀ | Eligible .5-unit pairs | Samples >π/4 | Largest sampled increment (rad) |
|---|---:|---:|---:|---:|
| 0.001 | 1 | 76800 | 0 | 0.1294 |
| 0.001 | 2 | 76800 | 0 | 0.0643 |
| 0.001 | 3 | 76800 | 0 | 0.0372 |
| 0.01 | 1 | 76796 | 17 | 2.5981 |
| 0.01 | 2 | 76729 | 6 | 1.5502 |
| 0.01 | 3 | 76450 | 3 | 0.9287 |

The largest stronger-input examples have endpoint R values 0.474→0.218 (C₀1), 0.140→0.107 (C₀2) and 0.148→0.145 (C₀3); their current-region counts are 1→2, 2→2 and 2→2. Thus a large global-angle increment alone does not identify uninterrupted relocation of a single activity site. Saved primary fields are spaced by 4 units and cannot resolve the detailed subinterval field evolution at scalar cadence .5.

## Goryachev forced-state qualification

At stronger input, the C20 deterministic mean control loses polarity. Its noisy trajectories therefore cannot be presented as clean spontaneous wandering of the original polarized state. Complete native histories, both strengths, vector/magnitude controls, eligibility and local-pattern statistics are preserved in the chemical interim sheet and original report. Strict lag-8000 site endpoints are too sparse for an inferential physical-site curve (1/0/0 eligible pairs for C20/40/60).

## Validation and source access

The main audit verifies all 81 Spring cases and 64,881 saved profiles against physical moments, counts, source hashes and exact time grids. The supplemental replay independently checks all saved dominant-region centroids (maximum material-coordinate error 2.23×10⁻¹⁶), regenerates all 72 exact noisy forcing paths, and retains first/last-half site and whole-seed increment-CDF analyses. Independent scientific review separately reproduces the core estimates and whole-seed intervals. Frozen pilot failures and later prospective estimator sensitivity gates remain explicit.

- [Two-sheet caption and integration metadata](supporting/metadata/content.json)
- [Spring methods and limits](SPRING_METHODS_EVIDENCE_LIMITS.md)
- [Core source and statistics audit](supporting/metadata/spring_noise_analysis_audit.json)
- [Centroid, half-window, CDF and forcing audit](supporting/metadata/spring_noise_supplement_audit.json)
- [Early/late endpoint measurements](supporting/data/spring_noise_short_long_endpoint_summary.csv)
- [Selected scalar increment records and original native archive paths](supporting/data/spring_noise_large_increment_native_context.csv)
- [Preserved original chemical interim](figure4_stochastic_revision.png)
