> **deBelly correction (2026-10-08):** the old in-loop centered endpoint extension is flawed and deprecated. Do not use its simulations for current results. The [independent-endpoint replacement](DEBELLY_BOUNDARIES.md) regenerates all deBelly preparations and three-test responses; other four models retain unchanged audited evidence. The replacement is complete and independently audited; [current results](data/runs/debelly_corrected_three_test_20261008/README.md) supersede the historical counts and centered deBelly results below.

# Current five-model, three-test comparison

The main polarity heatmap and length-versus-cue-duration supplement now focus on
exactly three response experiments: fixed physical width 3, relative width 20%,
and a finite gradient with constant physical slope. Spring's fixed support is
0.2 in its original native units. No micron conversion is implied for Spring.
Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring remain the five models.

## One exact initial cell for all three challenges

At every model/length, all three tests restart the exact nominal-amplitude,
relative-20% corrected basal-release checkpoint. No state is reflected,
reconditioned, or replaced by a separately prepared fixed-width state. The
initial 110 model/length preparations include 97 eligible and 13 ineligible
states; each appears in all three tests. Preparation-ineligible states do not
count as reversal failures.

This supersedes the old comparison's separate preparation condition for each
response width. Those preparations were very close, but not generally bitwise
identical. Changed starts and changed stimuli are integrated afresh. Audited
relative-patch and gradient trials can be reused only where the full checkpoint
digest and cue are exactly identical. Imported immutable native archives and
metadata undergo the same complete gate/selection audit as new trials.

## Goryachev: stronger, length-independent patch amplitude

The bounded pilot tests fixed multipliers K = 4, 16, 64 at the shortest, middle,
and longest eligible lengths, under both patch geometries. Each combination
starts from the original full checkpoint and independently tests durations
0, 5.12, 51.2, 512 and 5,120, with the original recovery and two scalar gates.
The lowest factor achieving both gates at all six geometry/length pairs is
selected before examining the full-length campaign.

The independently audited 90-trial pilot selected K = 64: factors 4, 16 and 64
achieved both gates at 2, 4 and 6 pairs respectively. Both Goryachev patches
therefore have peak 6.4 at every length. This multiplier does not affect the
initial preparation, basal dynamics, or the existing gradient slope. All other
models retain their nominal patch peaks. This is a response-permissive assay
choice, not a calibrated biological stimulus or a parameter-robustness result.

The relative patch support follows 20% of physical length; its peak remains
constant. The fixed patch support remains 3 (Spring 0.2), with no amplitude or
support rescaling. The existing fixed-slope gradient retains one g per model;
its peak increases with length by definition. Neither patch amplitude scales
with length.

## Unchanged finite-pulse outcome

Each positive pulse withdraws to basal input and recovers for 10,000 native
units. The original final-500 scalar orientation and 1%-strength replacement
gates, the separate 0.1% native-field flag, base duration grids, independent
zero controls and local nonmonotonic refinement remain in effect. All sampled
fail-to-pass boundaries must resolve within 1%; the endpoint is the shortest
tested passing duration, not a guarantee against arbitrarily narrow intervals.

The continuously applied gradient report is a separate experiment and does
not supply this figure's third column. Cue duration here excludes preparation
and recovery. It need not equal the time needed to completely eradicate the
old pole or arrive at a new full-field steady state.

## Current main-figure statistic (2026-10-10)

Panel C reports Spearman rank correlation between length and shortest passing
tested finite cue duration. Include every audited scalar success under its
original gate, regardless of the separate native-field tail flag. Require at
least four distinct successful lengths and nonconstant durations, uniformly
across all models, geometries and gates. Ties use average ranks. No significance
or effect-size threshold is applied. Constant durations give undefined
correlation; insufficient coverage remains unavailable. Censored and ineligible
lengths never receive invented durations. All 30 current cells meet the rule.

In-cell n counts included scalar successes. Native-field pass/fail/unknown
counts, excluded statuses and included length ranges are exported in
`data/derived/length_correlations.csv`; open markers in the duration supplement
retain failing field flags. `data/derived/correlation_validation.json` checks
all coefficients independently with pairwise average ranks. Correlations
describe monotonic association conditional on observed success; they do not
establish full-field convergence, capture every nonmonotonic feature, or imply
a scaling law. The previous field-stable-only OLS table is retained for
provenance and no longer supplies panel C. See the
[SciPy definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html).

## Historical campaign and mechanism summary

The completed duration suite contains 74,187 independently audited trials:
43,134 exact unchanged imported trials/controls and 31,053 newly integrated
trials. Eight local refinement rounds resolved every sampled boundary. Of 660
scalar readouts, 536 succeeded, 46 were censored and 78 were preparation-
ineligible; no execution errors remain. Goryachev succeeds at all 21 eligible
lengths in both patch tests and both gates. All 25 available descriptive fits
were independently checked by closed-form regression.

`code/mechanisms/three_test_mechanisms.py` reconstructs short, middle and long
eligible lengths in each cell, adds sampled interior duration extrema, and
retains neighboring passing/failing native pulses for both gates. Measurements
separate initial coherent width, pole peak, old-pole excess amount, new-pole
activation, old-pole removal, edge velocity, pool gradients and scalar gates.

Fourfold mobility, matched-gradient-peak, cue-only MCA freezing, cue-only Spring
feedback removal and longest-length 401-point mesh checks are explicitly
diagnostic interventions. They do not modify native comparison inputs or fill
in native response curves. Fixed-dose interventions identify trajectory effects;
they are not fully refined counterfactual cue-duration thresholds.

The suggested constant-speed removal mechanism has three independent
requirements: old-domain width grows with length, physical edge speed is nearly
length independent, and the measured removal clock explains the relevant
response clock. Growing peak height in a narrow pole and commitment before
old-pole disappearance are distinct possibilities. Positive descriptive slopes
alone establish none of these mechanisms.

The completed spatial campaign contains 128 native replays and 107 diagnostic
controls across all 15 cells. The evidence distinguishes five useful mechanism
classes: growing coherent domains versus cue coverage; weaker initial countercue
activation and stronger immediate incumbent restoring response in Goryachev;
input-strength
acceleration opposed by transport/domain competition; spatial MCA/mechanical
memory in deBelly; and extension-dependent global gating in Spring. Each cell
has its own interpretation, evidence, limits and testable prediction.

The current [Goryachev length-mechanism report](figures/Goryachev_length_mechanism/index.html)
consolidates its starting-state and small incumbent-suppression evidence. It
supersedes the earlier Goryachev-specific exploratory interpretation: equal
local countercues initially activate less material in longer cells, and their
incumbents show stronger immediate restoring reaction per unit damage. Late
recovery is slower; these coupled effects are not exhaustive causal fractions.

Initial right-pole source calculations use each native model's reaction term,
not cue amplitude alone. At equal patch peak, the long/short source ratio is
about 0.10 for Goryachev and near 1 for Wave-pinning, Holmes and Spring. Thus the
same length-independent amplitude does not imply the same biochemical forcing.
Pole contrast is measured above the contemporaneous median in the central 5%
of the cell; the edge threshold adds half the original pole contrast to that
background. Absolute half-height arrival is retained separately, so broad
gradient activation is not conflated with a coherent new pole.

Cue-only MCA freezing defeats all 12 sampled native-successful deBelly doses
across the three tests. Spring's cue-only feedback removal defeats all four
sampled gradient doses, while longer sampled patch doses can still reverse.
Fourfold mobility rescues the native censored L=39 fixed-width Wave-pinning
case at the same 5,120-unit pulse. These are fixed-dose causal probes, not
exhaustive counterfactual threshold campaigns.

Longest-length 401-point remeshing preserves 14 of 15 orientation classifications
and all 15 strength classifications at native selected orientation doses or
censor caps. deBelly fixed-width orientation changes; its exact longest-length
threshold is resolution sensitive. These checks do not establish mesh
convergence and do not test every independently selected strength endpoint.

## Execution and artifacts

Use the existing `scaling_sandbox` environment. Compute runs on SLURM QoS small,
one CPU and 3 GB per duration worker (4 GB per spatial worker), without an array
concurrency throttle. Executable scientific sources remain hash-frozen while
producers run. HTML has no timed self-refresh; status exposes generation time
and producer state.

The canonical outputs are `figures/main_polarity.png/pdf/svg`,
`figures/supplement_width_response.png/pdf/svg`, and
`figures/three_test_mechanisms/index.html` with its complete PDF. The supplement's
existing filename is retained for manuscript links, while its columns now
represent the three tests, not three preparation widths.

Campaign inputs, status, endpoints and audits live in
`data/runs/three_test_reversal_20261006/`; the pilot is
`data/runs/goryachev_fixed_amplitude_pilot_20261006/`. Spatial histories, the full
15-cell evidence table, measured quantities, controls and audit are in
`data/runs/three_test_mechanisms_20261006/`.

After completed execution and audits:

```
python polarity1d/manuscript/code/three_test_figures.py
python polarity1d/manuscript/code/mechanisms/three_test_analysis.py
python polarity1d/manuscript/code/build.py verify
```

The two rendering commands use compact comparison summaries and local native
spatial histories respectively. Historical campaign folders remain immutable.
