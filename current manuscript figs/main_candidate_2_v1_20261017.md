# Candidate 2: separation changes the fate of two activated regions

This is a provisional main-text candidate from the bounded circular-domain distance pilot. Fresh critiques and the subsequent mesh and imbalance checks are linked below.

![Holmes distance candidate](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/holmes_distance_candidate.png)

## Proposed message

Two challenged regions can activate similarly but interact differently after cue withdrawal. In this Holmes example, moving the cue from 12 to 16 µm changes the observed response from filling the valley into one broad activity region to two separated, narrower caps retained through the recorded cue-free window. This is a contrast in spatial reorganization, not an assumed law that a farther challenger activates more slowly.

The ring circumference is fixed at 40 µm. Both cases restart the exact same qualified N = 800 basal array and receive the same raised-cosine cue: total physical support 3 µm, native peak 1 µM/s, on for 1000 s, then off for 13,000 s. Native reactions, diffusion and conserved Rac/Rho pools are unchanged. The shorter and complementary routes are 12/28 µm and 16/24 µm. All panels share the same fixed display frame with the original incumbent at 0; only cue position changes. Holmes concentrations, time and length use the adapter's native µM, s and µm units.

## Why this design conveys the result

The two spatial histories resolve activation and the post-withdrawal filling of the intervening valley. The raw concentration profiles show the common starting state and its later organization, so activity remaining at the old position is not mistaken for an unchanged incumbent. The intermediate profile is 450 s after withdrawal; the final profile follows 13,000 cue-free seconds. Full 14,000-s histories remain in the supporting survey.

Early cue-region activity is almost matched: at 1 s the d = 12/16 means are 0.6594/0.6577 µM, and at 10 s they are 1.6938/1.6888 µM. Their later difference therefore occurs after substantial local activation in both cases. This descriptive comparison does not isolate the chemical cause of the difference.

The final d = 12 profile is a shifted single cap with the original width and peak concentration. At d = 16, the two peaks remain close to the original concentration, while each cap is narrower:

| Profile/site | Peak Rac (µM) | Half-local-contrast width (µm) | Active Rac above profile minimum in its valley-bounded region (µM·µm) |
|---|---:|---:|---:|
| Common basal cap | 1.484351 | 16.5379 | 18.5266 |
| d = 12 final single cap | 1.484351 | 16.5379 | 18.5266 |
| d = 16 final original-side cap | 1.484687 | 8.5226 | 9.5628 |
| d = 16 final cue-side cap | 1.484464 | 7.9209 | 8.8884 |

These are measurements of site organization. Active Rac alone is not conserved; the active-plus-inactive Rac pool is conserved. Whole-ring active amounts are 33.1293 µM·µm basally and 33.0451 µM·µm at the d = 16 endpoint. The table is not a causal allocation of a conserved active pool.

## What remains open

The d = 16 endpoint is not exactly symmetric: peak positions are approximately −0.10 and 16.00 µm, and both heights and amounts differ slightly. Two full-field late-window stationarity diagnostics pass, but the normalized cue-minus-original fixed-basin amount imbalance still changes by −4.78 × 10⁻⁶ over the last 1000 s. Two sites persisting over this finite interval do not establish a stable branch or exclude slower competition. The completed targeted checks below support this finite-window contrast.

This fixed-circumference intervention does not establish a monotone distance law, minimum successful cue duration, a front-establishment gate, isolated diffusion travel time, migration persistence or a biological WGD calibration. It also does not establish that the same two-site response occurs in every model.

## Targeted follow-up checks

Fresh N = 1600 preparations and challenges reproduce one cap at d = 12 and two caps at d = 16. Terminal full-field differences from N = 800 are below 8.97 × 10⁻⁵ on the frozen field scales. The largest transient difference is larger (0.01685), in active Rho during the d = 12 merger; the check supports the spatial outcome without establishing identical merger timing. [Mesh result and full profiles](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/holmes_distance_robustness_20261017/supporting/metadata/MESH_RESULTS.md).

From the N = 800 d = 16 endpoint, a concurrent unperturbed control and both signs of an active-Rac transfer equal to 1% of the smaller fixed core's above-background active amount each retain two caps through another 13000 cue-free seconds. Global Rac is conserved, and the other fields are initially unchanged. About 2.2% of the imposed imbalance difference remains at the endpoint relative to control; most of the reduction occurs early, followed by slight later regrowth. This supports robustness to this specified disturbance over that window. [Imbalance result, exact perturbation and full fields](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/holmes_distance_robustness_20261017/supporting/metadata/IMBALANCE_RESULTS.md).

## Supporting results and reproducible sources

- [Full scientific report with all eleven distance cases and three reflected-cue numerical controls](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/metadata/RESULTS.md). The controls are not biological replicates.
- [Early spatial survey](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/distance_spatial_early.png) and [full spatial survey](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/distance_spatial_full.png). Far Wave-pinning cases show loss near the original location and a remaining displaced cap; they should not all be labelled merger. Goryachev forms a separate cue-local maximum before original-site activity declines.
- [Native site organization](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/data/holmes_site_organization.csv), [early activation](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/data/holmes_early_activation.csv), [two-site heights, widths and valleys](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/data/holmes_two_site_metrics.csv), and [fixed-basin tail drift](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/data/holmes_two_site_tail_drift.csv).
- [Frozen run manifest](../polarity1d/manuscript/data/runs/chemical_ring_distance_pilot_20261017/manifest.json), [d = 12 full fields](../polarity1d/manuscript/data/runs/chemical_ring_distance_pilot_20261017/workers/10), [d = 16 full fields](../polarity1d/manuscript/data/runs/chemical_ring_distance_pilot_20261017/workers/11), and [verification](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/metadata/verification.json). SLURM array 21083222 completed all 14 workers. The renderer checked 72 relevant hashes, exact starts/phase continuity, native integrals, peak/valley values and conserved pools.
- [Figure generator](../polarity1d/manuscript/code/mechanisms/chemical_ring_distance_report_20261017.py) and [vector PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_distance_20261017/supporting/vectors/holmes_distance_candidate.pdf). The figure is generated from retained arrays without new integration.

[The first fresh critique](reviews/chemical_distance_r1.md) accepts the compact comparison as a candidate about post-activation spatial reorganization and prompted the targeted checks. [A different fresh critic](reviews/holmes_robustness_r1.md) independently checked those completed results and accepts their support for the finite-window claim. The sampled merger diagnostic differs by 5 s between grids; the perturbation residual relaxes rapidly and then grows slightly. Neither check establishes asymptotic stability or resolves the migration interpretation.
