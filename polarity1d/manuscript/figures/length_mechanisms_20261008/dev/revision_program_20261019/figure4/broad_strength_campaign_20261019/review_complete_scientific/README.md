# Complete scientific pilot: organization and angular memory

[Two-page result review](broad_strength_review.pdf) ·
[Overview PNG](broad_strength_exploratory.png) ·
[Four-model native profile context](broad_strength_transition_profiles.png) ·
[Raw/qualified tables](index.html) ·
[Exact export receipt](supporting/metadata/RENDER_RECEIPT.json)

All **576 scientific cases** are collected: 324 retained cases plus 252 new
cases. The overview is the snapshot observed at **19:36:02 UTC, 19 October
2026**, with three of four numerical companions complete. WP, Goryachev and
Holmes numerical screens have independently passed; Spring qualification is
still pending in this snapshot. No replication or range extension has been
launched. The [accepted earlier 229-output snapshot](../review/README.md) and
the preceding main Figure 4 sheets remain unchanged. This is a discussion
artifact, not an agreed final Figure 4 composition.

## What the broader strengths show

WP, Holmes and Spring reach frequent multiple-region organization with retained
contrast at stronger tested inputs. Large Goryachev instead often has low
contrast relative to its own basal profile; raw half-contrast region counts can
remain high there. This is a difference within the sampled forcing family and
range, not evidence that a model cannot share another organization elsewhere.
Small Goryachev also exhibits substantial contrast-qualified multiple-region
occupancy, so there is no model-wide binary separation.

The largest sizes make the change across the common **0.01 → 0.03 → 0.1** grid
interval particularly clear. Values below are full-window means across seeds.
The original 0.01 level has twelve seeds; 0.03 and 0.1 currently have two each.

| Model, largest original size | Qualified multiple-region occupancy | Low-contrast occupancy |
|---|---:|---:|
| WP, C60 | 0.06% → 21.6% → 95.8% | 0% → 0% → 0% |
| Goryachev, C60 | 0% → 6.3% → 3.5% | 0% → 85.5% → 92.6% |
| Holmes 3, C60 | 17.1% → 74.3% → 89.9% | 0% → 0% → 0% |
| Spring, C₀3 | 14.1% → 48.4% → 79.9% | 0% → 0% → 0% |

“Qualified” means more than one connected current half-contrast region and
contrast at least 10% of that case's initial basal contrast. Low contrast means
below that same 10% level. Both fractions use all saved times in the denominator.
The linked HTML shows raw and qualified fractions **adjacent**, preventing a
raw count in a nearly flat field from being mistaken for substantial peaks.
For example Goryachev C60 at f=0.3 has raw multiple-region occupancy 96.43%,
qualified occupancy 1.16%, and low-contrast occupancy 96.44%.

These patterns depend on size. At f=0.1, qualified occupancy across small,
middle and large sizes is **13.5/73.6/95.8% for WP**, **21.5/73.7/89.9% for
Holmes**, and **4.7/47.1/79.9% for Spring**. Goryachev instead gives
**58.6/25.4/3.5%**, with low-contrast occupancy **7.5/66.5/92.6%**. These
two-seed contrasts locate questions for replication; they do not establish
precise population frequencies or a general monotone size law.

The positive sampled runs are more than isolated threshold crossings in the
largest-size f=0.1 WP/Holmes/Spring paths. Their per-seed 90th-percentile sampled
spans are respectively **3900–5811**, **1048.5–1167**, and **21.5–24.6** native
units, observed at cadences 5, 5 and 0.5. A span is last minus first consecutive
positive sample time. Edge-censored runs remain marked in the export; these
are not lifetimes of identified sites, uninterrupted lineage or event rates.
The models' native clocks are not calibrated to one another.

Active amount and coherence distinguish additional responses. At f=1 the
smallest WP/Holmes cells have mean active-amount ratios 1.626/2.147 and mean
coherence ratios 0.544/0.587; their mean R is 0.884/1.278. Spring C₀1 instead
has mean amount 1.194, coherence 1.765 and R=2.144, each relative to its own
basal. R can exceed one. A higher or lower R alone does not identify an amount
change, a translated coherent profile or the number of active regions.

## Do organization changes align with angular-memory changes?

They overlap in several conditions but are not equivalent. The same
0.01 → 0.03 → 0.1 largest-size sequence gives the following global-direction
cosines at lag 8000 for chemicals and lag 800 for Spring, with eligible-pair
coverage alongside:

| Model | Direction cosine | Direction-pair coverage |
|---|---:|---:|
| WP | 0.972 → 0.570 → −0.045 | 100% → 98.1% → 92.9% |
| Goryachev | 1.000 → 0.152 → −0.149 | 100% → 22.2% → 7.9% |
| Holmes 3 | 0.738 → 0.048 → −0.024 | 99.2% → 98.9% → 96.2% |
| Spring | 0.028 → 0.008 → −0.030 | 99.3% → 97.9% → 97.9% |

Direction pairs require R≥0.1 at both endpoints. The estimator pools the
cosine sums and eligible counts across whole seeds; coverage divides eligible
counts by all available lag pairs. No seed or pair is dropped because its
correlation is unfavorable. Sparse estimates are labelled unavailable for
interpretation when fewer than 100 pooled eligible pairs or less than 1%
coverage remain; contributing/supported seed counts are exported. In particular
Goryachev C60 f=0.3/1 has insufficient support. Its apparent direction loss is
therefore qualified by the weak-resultant/mean-input regime, unlike the high
coverage in the WP/Holmes/Spring comparisons. Mean-input controls themselves
lose the first moment in relevant high-input conditions; this is not clean
spontaneous wandering around the original basal state.

Spring supplies a useful dissociation: already at f=0.003, C₀3 has lag-800
direction cosine 0.543 with 100% coverage, mean R=0.992, and only 0.66%
qualified multiple-region occupancy. Thus appreciable angular-memory loss can
precede frequent multiple regions. Conversely, rising R at strong forcing
does not restore directional memory: C₀3 R rises from 0.703 at f=0.1 to 1.280
at f=1 while the corresponding long-lag cosines remain near zero. These are
global first-moment directions, not site travel, cell movement or a causal
allocation of memory loss to changing organization.

## What the native context page shows

The separate four-panel page overlays exact basal and f=0.01/0.03/0.1 native
profiles at the largest size. It takes the smallest saved seed number in each
condition and the exact window midpoint: t=16000 for chemicals, t=1600 for
Spring. Selection is independent of the field, region count or R; some selected
intermediate-strength snapshots therefore still contain one region. Legends
state the actual sampled region count and R. These individual snapshots are
illustrations alongside the full-window occupancy, not representative marginal
states or paired forcing paths. They are unshifted and preserve native
concentration scales. Goryachev includes a native-scale enlargement; Spring
uses material position and records actual physical C(t) in the exported data.

The two-page core also retains its C20 f=1 WP/Gory diagnostic comparison with
the separately declared conditional selection at t=16040. The two selections
serve different purposes and are never pooled for prevalence estimates.

## Bounded replication proposal, not submitted

Bring three **common** new levels, **f=0.003, 0.03, 0.1**, from two to twelve
seeds across all four models and all three sizes. This is 36 condition groups
and **360 additional noisy workers**, with existing mean/basal controls reused
on exact protocol proof. The 0.003 level addresses early memory loss before
frequent multiple regions; 0.03 resolves the broad observed change; 0.1 checks
the higher-occupancy/low-contrast contrast across models and sizes. Keeping
common levels avoids selecting different favorable doses for each model.
The retained 0.001/0.01 levels already have twelve seeds. The remaining broad
endpoints retain exploratory status pending author discussion.

The current two-seed spreads and window sensitivity justify replication rather
than finer strength spacing first: WP C60 f=0.03 qualified occupancy ranges
4.2–39.0% across seeds, with a 31.1-percentage-point first/last-half difference
and a 0.341 directional-cosine difference. Holmes C60 f=0.03 has a 25.7-point
half-window occupancy difference. Spring C₀3 f=0.03 has only a 0.66-point
occupancy difference but a 0.129 direction difference. More seeds can estimate
finite-window variability; they do not by themselves establish stationarity.
Matched full/first-half/last-half estimates must therefore remain visible.

Measured ordinary-worker timings at these exact conditions predict **63.1 CPU
hours** for 360 additional cases. The longest of the 72 measured input runs is
2139 s. A proposed small-QoS array would retain one CPU, 3 GiB and a two-hour
cap per worker, without a concurrency throttle. Measured compressed/allocated
case storage predicts **3.83 GB** additional native output; a 30% reserve is
**4.98 GB**. These are estimates from two paths per condition, not guaranteed
upper bounds. Proposed seed identities and exact input hashes are exported.
No job or native protocol has been created or submitted for this proposal.
Spring numerical/observation acceptance, independent source preflight and root
review precede launch. Author feedback can steer allocation; it is not an
additional authorization gate.

## Data, definitions and scientific limits

- [All condition/window measurements and observed seed ranges](supporting/data/condition_window_snapshot.csv)
- [Memory: every model/size/strength/control/window/lag](supporting/data/memory_all_conditions_windows_lags.csv)
- [Joint organization/memory and coverage](supporting/data/organization_memory_joint_windows.csv)
- [Seed/window sensitivity](supporting/data/organization_memory_seed_window_sensitivity.csv)
- [Sampled positive-run spans and censoring](supporting/data/sampled_qualified_run_spans.csv)
- [Native change-context selection and source hashes](supporting/data/transition_profile_selection.csv)
- [Every change-context native profile cell](supporting/data/transition_profile_cells.csv)
- [Unsubmitted replication conditions and measured resource inputs](supporting/data/replication_proposal_conditions.csv)
- [Unsubmitted planned seed/input identities](supporting/data/replication_proposal_unsubmitted_cases.csv)
- [Replication proposal receipt](supporting/metadata/REPLICATION_PROPOSAL.json)
- [Preceding methods: physical M₀/M₁/P, amount/coherence factorization](../review/README.md)
- [Frozen native scientific plan](../PLAN.md)

All levels retain the original positive lognormal family and native spatial/time
correlations. Scaling f changes mean and fluctuation amplitude together; matched
mean controls use the declared ensemble mean, not a realized-path mean. P is
the normalized physical-activity vector M₁(t)/M₁(0), q₁=|M₁|/M₀ and R=AQ at
each frame. Mean R is the mean paired product, not the product of mean amount
and coherence. Actual Spring circumference enters physical quadrature.

Full chemical windows are 0–32000 at scalar cadence 5; Spring uses 0–3200 at
0.5. Native profile cadence is 40/4 respectively. Full/first-half/last-half
windows retain the same lag checkpoints; the midpoint is included in both
halves. Native amplitudes, time units and effective forcing are not calibrated
across models. The threshold and prominence sensitivities are descriptors,
not unique biological boundaries. Appearance, coexistence, low first moment,
low contrast and amount/coherence changes may overlap. This finite grid cannot
establish discrete phases, functional competition, site lineage, event rates,
causal pathways or absence of another organization outside the explored range.
