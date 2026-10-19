# Existing weak/strong noisy-ring regimes: author discussion inventory

19 October 2026. This bounded reanalysis responds to the question of whether
Figure 4 shows only one input regime and one wandering peak. It changes no
simulation, forcing, main artwork or main-panel selection. The stored histories
contain distinct local organizations, but the current main histories expose
that diversity unevenly. These are observed patterns, not functional competition
or mechanistic event classes.

## What the current main sheets actually show

| Current panels | Input and selection |
|---|---|
| A1–A3, native histories | Stronger input, f=0.01 only; existing condition medoids: WP C20 seed202610180018, Holmes C60 seed202610180205, Spring C0=2 seed202610190446. |
| B1–B3, directional memory | Both f=0.001 dashed and f=0.01 solid, all three sizes and all 12 seeds per group. Shading is the stronger-input whole-seed interval. |
| C1–C3, site endpoints | f=0.01 only; chemical lag8000, Spring lag800, with the stated site gates and pair coverage. |
| D1–D3, paired remapping | f=0.01 only; chemical lag8000, Spring early lag4. These are eligible site endpoints, not the global direction in B. |
| E1–E3, pattern/magnitude | f=0.01 only; all6401 scalar samples per seed, with no selected adaptation exclusion. |

Goryachev is absent from these two main sheets, but its complete histories,
both strengths and controls remain in the chemical supporting evidence.
The existing [18 chemical medoids](../../../vision_noise_20261017/vision_noise_selected_histories.png)
and [six Spring medoids](../figure4_spring_all_histories.png) already provide
matched-size weak/strong comparisons. The histories at different strengths use
independent seeds, not a paired common noise realization.

The displayed Holmes medoid has multiple half-contrast regions at **1.78%**
of samples; its group mean is **17.15%**. The medoid combines several frozen
descriptors and is not selected to reproduce every marginal group mean. The
displayed Spring C0=2 medoid has **3.08%** multiregion occupancy; its largest-size
group has substantially more. WP remains predominantly one coherent region.
Thus the current examples should not be read as an inventory of all behaviors.

## Strength increases both mean activation and fluctuation amplitude

Both campaigns use positive extra activation with a unit-variance Gaussian
field, stationary OU knots and the declared normalized interpolation:

\[
I(s,t)=fA_m\exp\!\left[Z(s,t)-\frac{1}{2}\right],
\qquad \mathbb{E}[I]=fA_m,
\qquad \operatorname{Var}(I)=(fA_m)^2(e-1).
\tag{1}
\]

Increasing f from0.001 to0.01 raises both the ensemble mean and standard
deviation tenfold, and variance100-fold; the local coefficient of variation
remains approximately1.311. This is **not a zero-mean variance sweep**.
The reference amplitudes A_m are WP0.1, Goryachev6.4, Holmes1 and Spring4.
Each strength has its own spatially uniform mean-input control; the basal
control has no extra input. Chemical correlation scales are ell1.5/tau50,
whereas Spring uses its declared native ell0.15/tau5 material law. These are
not calibrated equal biological exposures, and equal local covariance does
not fix integrated or first-angular-mode forcing across sizes.

## Prevalence and observed spans across every size

At **f=0.001**, every model/size group has zero sampled multiregion and zero
weak-resultant occupancy. This concerns the declared thresholds and observation
window, not an exclusion of unsampled brief changes or small secondary peaks.

At **f=0.01**, the following fractions use all6401 scalar samples per seed.
Brackets are the unchanged accepted95% whole-seed intervals. The last column
reports the median/maximum, across12 seeds, of each seed's longest consecutive
multiregion sample-run span, in that model's native time.

| Model, size | Multiregion time % [95% interval] | Seeds with any /12 | Weak R<0.1 time % | Longest sampled span: median / maximum |
|---|---:|---:|---:|---:|
| WP C20 |0.0013 [0,0.0039]|1|0|0 /0|
| WP C40 |0.0547 [0,0.1589]|3|0|0 /110|
| WP C60 |0.0560 [0.0130,0.1159]|6|0|2.5 /95|
| Goryachev C20 |50.5103 [49.23,51.85]|12|11.0530|600 /750|
| Goryachev C40 |0 [0,0]|0|0|0 /0|
| Goryachev C60 |0 [0,0]|0|0|0 /0|
| Holmes C20 |1.1535 [0.2721,2.268]|8|0.0625|97.5 /1935|
| Holmes C40 |2.8069 [0.7368,5.437]|11|0.1237|207.5 /2960|
| Holmes C60 |17.1458 [4.382,32.93]|12|0.5273|317.5 /21460|
| Spring C0=1 |0.0482 [0.0247,0.0729]|9|0.0026|0.25 /1.5|
| Spring C0=2 |3.4695 [2.937,4.085]|12|0.0768|24.5 /44.5|
| Spring C0=3 |14.1189 [11.90,16.39]|12|0.4049|101.75 /186.5|

A run is consecutive positive **observations**, with span equal to last minus
first positive sample time. Chemical scalar cadence is5 and saved-field
cadence40 over32000; Spring cadences are0.5 and4 over3200. Adjacent samples
do not establish uninterrupted persistence between them, and numerical spans
cannot be ranked across uncalibrated native clocks. Zero span can indicate
either no positive sample or an isolated positive sample; the seed counts
distinguish these cases. Censoring and neighboring negative times remain in
the [run table](supporting/data/sampled_runs.csv).

The Holmes C60 distribution is heterogeneous: the largest run spans21460 in
seed202610180212, while the median seed's longest run is317.5. The largest
run is descriptive supporting context, not a representative selection.

## Prominent peaks add a different kind of diversity

Multiple half-contrast regions are disconnected sets above each current
profile's minimum plus half its range, with periodic connectivity. Multiple
peaks use the separate native prominence threshold0.05. A secondary peak can
sit below half-height and leave only one half-contrast region.

| Model | Multiple-prominent-peak time %, weak input: small/middle/large | Strong input: small/middle/large |
|---|---:|---:|
| WP |0 /0 /0|0.379 /1.695 /2.806|
| Goryachev |72.13 /90.50 /94.77|99.80 /99.98 /99.98|
| Holmes |0 /0.0013 /0.0182|4.036 /17.99 /41.10|
| Spring |0 /0 /0.0195|0.0924 /6.595 /23.90|

In particular, weak-input Goryachev is not literally a smooth isolated peak:
many samples have additional resolved local maxima, although their half-height
organization remains one region. Its narrow dominant peak and common native
color scale can conceal these smaller amplitudes. Prominence is not a measure
of a second functional polarity site. New whole-seed intervals for this
descriptor are available in the [condition inventory](supporting/data/condition_inventory.csv).

## Concrete native examples, without changing representatives

These examples lie inside the **existing condition medoids**. Within each
medoid, the selected interval is its longest consecutive sampled run having
multiple regions, contrast at least0.05 and multiple0.05-prominence peaks;
the listed native field is the saved profile nearest the interval midpoint.
This explicit qualification is a descriptive sensitivity gate. Primary
occupancy fractions above remain ungated.

| Existing medoid | Positive sampled interval; span | Native profile | Measured organization at that profile |
|---|---|---|---|
| Holmes C60, f0.01, seed202610180205 |9385–9740;355|t9560|2 half-contrast regions,3 prominent peaks, contrast1.146, R0.795;9 saved fields in interval.|
| Spring C0=3, f0.01, seed202610190471 |130–280.5;150.5|t204|2 regions,2 peaks, contrast0.537, R0.159;38 saved fields in interval.|
| Goryachev C20, f0.01, seed202610180086 |16670–17300;630|t17000|3 regions,3 peaks, contrast2.224, R0.103;16 saved fields in interval. Mean-control qualification below is essential.|

The actual current A3 Spring C0=2 medoid also has a qualifying interval
2752–2765.5, with a saved field at2760 showing2 regions/2 peaks/R0.313.
WP's displayed C20 medoid has no multiregion sample; that is faithful to its
predominantly single-region regime. Native archives and their hashes, exact
field times/scales and the same selections at every size are retained in
[native medoid examples](supporting/data/native_medoid_examples.csv).

## Mean controls define the interpretation

Every no-extra-input control remains one region/one prominent peak with R
approximately1. All WP, Holmes and Spring mean-input controls also remain
one region/one peak and never enter R<0.1. Their R values need not remain1:
for example, stronger Spring mean controls settle near0.973/0.979/0.980.

**Stronger-input Goryachev C20 is a different forced-state regime.** Its uniform
mean control is weak-resultant at98.30% of sampled times, has contrast<0.05
at96.98%, and ends with effectively zero first moment. Its noisy multiple
regions cannot be presented as spontaneous breakup of an otherwise intact
mean-driven pole. Goryachev C40/C60 retain one region in their mean controls
but adapt to final R approximately0.743/0.836. Even weaker Goryachev means
shift the magnitude to approximately0.956/0.980/0.985. Mean controls separate
adaptation from fluctuations; their differences do not yield additive causal
fractions. R remains a relative first-moment magnitude, not active amount or
proof of coherent translation. The existing
[amount/coherence analysis](../SPRING_AMOUNT_SPATIAL_BALANCE.md) provides that
separate distinction for Spring.

## Scope, evidence and a possible next discussion view

The defensible inventory is: predominantly coherent-region wandering in WP;
additional local maxima and a small-ring forced-state qualification in
Goryachev; intermittent or sometimes long-lived sampled multiple-region
organization in strong Holmes; and coherent weak-input versus appreciably
more multiple-region organization in strong larger Spring. Counts alone
identify neither functional competition nor an abrupt switching mechanism,
and global activity angle is not cell movement or site lineage.

If the author wants a diagnostic next, a compact matched-size weak/strong
comparison can reuse the existing medoids at WP C20, Goryachev C20, Holmes C60
and Spring C0=3, with one native color scale per model, occupancy/run context
and the visible Goryachev mean-control caveat. A declared native-time zoom
would expose multiple-region examples more clearly than a full-history strip.
The author subsequently requested that separate review artifact: see the
[weak/strong native histories and snapshots](../noise_strength_review/README.md).
It is discussion support, not a main-Figure4 replacement.

The [read-only source](../../../../../../code/mechanisms/revision_figure4_pattern_regime_review_20261019.py)
checks324 trace hashes and324 primary-field hashes, independently replays
region counts in all259524 stored profiles, and reproduces the prior
per-seed fractions exactly. The five compact CSVs cover all324 cases,
60 condition/control groups and all24 unchanged medoids. New descriptor
intervals use2000 whole-seed resamples; deterministic controls receive no
bootstrap replicas. [Validation and hashes](supporting/metadata/validation.json)
also verify the six current main PNG/PDF/SVG files remain unchanged.
Runtime was approximately9seconds on one CPU with80MiB peak memory; no
SLURM-heavy computation or native integration was needed.
