# Manuscript foundation audit — 17 October 2026

This audit supports the requested vision-aligned outline. It inspects the current
files, figure renderers, native-result audit summaries and direct user messages;
it does not infer approval from an agent-authored status note. The canonical
five-model comparison and the integrated literature figure should be retained
as the first publication checkpoint. Their scientific roles are complementary:
length changes finite-cue polarity replacement, while published migration
summaries motivate a possible persistence connection. Neither result alone
establishes that enlargement causes persistent migration.

## Sources and direct author preferences

The governing sources are [MANUSCRIPT_VISION.md](../MANUSCRIPT_VISION.md),
[root agent notes](../AGENTS.md),
[manuscript agent notes](../polarity1d/manuscript/AGENTS.md),
[current figure collection](../current%20manuscript%20figs/README.md),
[current report status](CURRENT_STATUS.md),
[THREE_TESTS.md](../polarity1d/manuscript/THREE_TESTS.md),
[PROTOCOL.md](../polarity1d/manuscript/PROTOCOL.md), and
[ASSUMPTIONS.md](../polarity1d/manuscript/ASSUMPTIONS.md).

Direct user records were inspected in the local Codex sessions. These paths are
evidence locations, not new repository dependencies; entire conversations and
the collaboration email have not been copied into this audit.

| Date and direct record | Brief author wording | Consequence for the outline |
| --- | --- | --- |
| 3 Oct, `rollout-2026-10-03T13-01-00-01a102b6-4080-7f63-acb9-d9ae521c059f.jsonl`, line 1067 | “Very good”; then asks to make MSD the new B and enlarge tracks in C | The six-panel 3-Oct layout implements a positively received foundation and requested layout revision. This is stronger evidence than an agent claim that the whole final figure was approved. |
| 3 Oct, `rollout-2026-10-03T15-43-13-01a1034a-c066-7422-9912-c124934aef38.jsonl`, line 365 | “only 5 models”; excludes Otsuji, authorizes clean rerun | Retain Wave-pinning, Goryachev, corrected deBelly, Holmes 3 and Spring. |
| 6 Oct, `rollout-2026-10-06T11-20-10-01a111cd-03a2-7a32-a5fd-2c88eb57227e.jsonl`, line 150 | “application duration required to achieve persistent reversal” | The third comparison assay is a finite gradient followed by basal recovery, not continuous-gradient first passage. |
| 8 Oct, `rollout-2026-10-08T11-28-06-01a11c20-fb4f-7560-b10f-e82159617f2d.jsonl`, line 897 | “Seems good”; mark old extension flawed and repeat all three tests | Use only the corrected independent-endpoint deBelly replacement. |
| 9 Oct, `rollout-2026-10-09T12-46-24-01a1218f-0847-7cc3-ae14-586bb26ef0b1.jsonl`, line 8 | “prefer proper legends”; wants descriptions, hashes and reproduction instructions | Put detailed interpretation in legends and text, keep panel text economical, retain provenance. The same message requests attention to decreasing portions of complete response curves. |
| 10 Oct, `rollout-2026-10-10T15-55-51-01a12762-d4a0-7dc2-9cbd-9beedf4ffafc.jsonl`, line 9 | “Use Spearman”; requests fair inclusion and fewer plot annotations | The 10-Oct correlation figure is the intended latest main comparison, superseding the older OLS version. |
| 13 Oct, `rollout-2026-10-13T13-27-58-01a1364e-870d-7dd3-b0b8-0021c62eb57a.jsonl`, lines 148 and 214 | “reserve lines only for the model”; subsequently questions implied measured-lag count | A later panel-B redesign remains open. Distinguish samples from digitized lines from actual author-calculated lag points. No later approval of the proposal was found in that session. |
| 17 Oct, `rollout-2026-10-17T14-45-31-01a14b2e-f887-7790-86d6-fb01de6abc4f.jsonl`, lines 124, 215, 307, 337 | “theoretical basis for polarity modification”; “existing data supports”; “much closer”; then bundles mechanisms with measurement proposals and reserves the experimental gap | The outline needs a coherent theoretical/empirical argument, an integrated mechanism-to-measurement component, and a prominent planned-assay results figure. A mechanism checkpoint does not define the whole paper. |

All session filenames above are under `/home/4473331/.codex/sessions/2026/10/`
and the corresponding day directory. Role was checked as `user`, rather than
reading quoted agent summaries as user instructions. These are direction and
selective positive feedback, not blanket acceptance of every scientific claim.

## Exact assets to integrate

SHA-256 values were measured on 17 October. Canonical PDF sources are preferable
for the outline where the PDF assembler can preserve vector content.

| Asset | SHA-256 |
| --- | --- |
| `polarity1d/manuscript/figures/main_polarity.png` | `b9546b972968fc759d00fdea48b1a52bc53f595a9adc80667f97ee40e47684db` |
| `polarity1d/manuscript/figures/main_polarity.pdf` | `1348072eaa20a95081bd705c42cf9090d66d2a6adfd76dc162442ad26d9bb536` |
| `polarity1d/manuscript/figures/supplement_width_response.png` | `ee4802daaa8f257db99e8aceb2fbed359f7b4ecd68420d9dcdb6fb2120819aeb` |
| `polarity1d/manuscript/figures/supplement_width_response.pdf` | `567f69e42ab02a07d75eddc6636c34137c8240099ba1f41a93751f813d7c1257` |
| `wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU.png` | `91723a8c10acd88e3ce1ee545ae116332676a32e0e9c34812e540500d97ffa6a` |
| `wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU.pdf` | `31de9b1be5fd27c495dfd1d8116313e1f8667141b839e24dc04755e2979f03fb` |
| `polarity1d/manuscript/figures/CAPTIONS.md` | `bd6844fbc4b5a3d61ae01d1f4ea4d115be6f46fb456eacd3b49e91ac0afcb20c` |
| `polarity1d/manuscript/data/derived/length_correlations.csv` | `014ab5b97222cf5bdc00f07bb71af4e204facda54f612efd2ac2f1acb8172791` |
| `polarity1d/manuscript/data/derived/figure_generation_provenance.json` | `534ee6ad1a6110320ed2c822505cf7f933593a30fb85714b977d7b38934d3353` |
| `wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU_provenance.json` | `4b4e1fcca0a81aedc54c38db6ba13a3aa4d06573d79aa29ddbd23e753a67b046` |

The working `current manuscript figs/main_polarity.png` was stale at audit
start: SHA-256 `e647322f16a812c8bc60141962ae4e8757074061573b8cfd8c354d015343c644`,
the 8-Oct OLS-era render. It has now been replaced by a byte-identical copy of
the canonical 10-Oct Spearman figure, on the parent agent's request. The working
width supplement was already byte-identical to canonical and remains so. The
working migration PNG is also byte-identical to the canonical source. No
canonical scientific source, trial or figure was modified by this audit.

## Five models and three finite-cue assays

### Detailed proposed main-figure legend

**Cell length changes the finite stimulation required for lasting polarity
replacement in a model- and cue-dependent manner.** A, recorded Wave-pinning
establishment at nominal length 15: a finite left-end raised-cosine patch is
withdrawn, after which native active profiles and the scalar orientation
approach a verified basal polarized endpoint. The relative-20% support equals
3 at this length. B, the three rightward challenges illustrated at lengths 15
and 30: a fixed-support cosine patch, a cosine patch covering 20% of the domain,
and a finite ramp with constant physical slope. Fixed support is 3 in the
non-Spring models and native 0.2 for Spring; no micron conversion is assigned
to Spring. Patch peaks are independent of length. Goryachev uses a single
audited patch multiplier of 64, peak 6.4, at all lengths in both patch assays;
the other models use their native nominal peaks. The gradient peak increases
with length by construction while its slope stays fixed within each model.
At each model/length all three challenges restart the exact same qualifying
relative-20% nominal basal checkpoint, without reflection or reconditioning.
Every tested pulse starts independently from that checkpoint and is followed
by 10,000 native time units of basal recovery. C, Spearman rank correlation
between length and shortest tested passing cue duration for orientation
reversal and scalar strength replacement, shown separately for all five
models and three assays. Orientation requires all 201 samples of the final
500-unit observation window to have orientation above 0.05, contrast at least
0.05 and orientation range below 0.005. Strength additionally requires the
reversed scalar magnitude to agree with the initial magnitude within 1%.
Each correlation includes every audited scalar success, regardless of its
separate native-field tail flag, and requires at least four distinct successful
lengths with nonconstant durations. Ties receive average ranks. In-cell
counts give the included successful lengths; censored and preparation-ineligible
lengths receive no invented duration. All 30 cells satisfy the rule. The
uniform color range is −1 to +1, without a significance or effect-size cutoff.
The separate pointwise native-field range/drift diagnostic remains available
in the complete duration supplement. Scalar success does not establish
full-field steady state. Pulse duration excludes preparation and recovery;
success-conditional correlations summarize monotonic association, can omit
difficult lengths, and do not establish a scaling law. deBelly uses corrected
independent local endpoint commands, not the deprecated centered extension.

The source legend is [CAPTIONS.md](../polarity1d/manuscript/figures/CAPTIONS.md).
The [generation provenance](../polarity1d/manuscript/data/derived/figure_generation_provenance.json)
records generation on 10 October from the corrected 8-Oct campaign.

### Detailed proposed full-curve supplement legend

**Complete length dependence of the shortest passing finite cue.** Five model
rows are crossed with the three challenge geometries. Blue circles/solid lines
show orientation reversal and orange squares/dashed lines show scalar strength
replacement. Filled success markers pass the independent native-field tail
flag; open success markers fail it. Open triangles mark no successful pulse
within the tested duration range, at the explicit censor cap. Crosses below
the axis denote preparation-ineligible lengths. Cue duration is logarithmic
and excludes preparation and 10,000-unit recovery. Every sampled fail-to-pass
duration boundary is refined within 1%; endpoints remain the shortest tested
successes, so unobserved narrow success windows cannot be excluded. All three
assays use the exact shared preparation at each length. The curves expose
nonmonotonic structure and censoring that cannot be reconstructed from the
main correlation heatmap. Units and size conventions remain model-native;
curves do not calibrate comparable biochemical doses or clocks across models.

### Verified status and narrative consequences

The corrected campaign's [status](../polarity1d/manuscript/data/runs/debelly_corrected_three_test_20261008/results/status.json)
and [independent audit](../polarity1d/manuscript/data/runs/debelly_corrected_three_test_20261008/results/independent_audit.json)
agree: 330 model/length/assay cases, 660 scalar readouts, 74,183 audited trials,
536 successes, 46 censored and 78 preparation-ineligible readouts; no execution
errors or missing chunks, and an empty audit error list. The 30 correlations
were independently checked with pairwise average ranks. The correlation table
contains 431 field-tail passes among 536 scalar successes. The separate
field-stable-only OLS table is historical provenance for this main panel.

| Model | Fixed-width orientation / strength correlation | Relative-width orientation / strength | Constant-slope orientation / strength |
| --- | ---: | ---: | ---: |
| Wave-pinning | 1.000 / 1.000 | 1.000 / 1.000 | 0.608 / 0.632 |
| Goryachev | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| Corrected deBelly | 0.406 / 0.291 | 0.550 / 0.476 | 0.427 / 0.385 |
| Holmes 3 | 1.000 / 1.000 | 1.000 / 1.000 | 0.047 / 0.406 |
| Spring | 1.000 / 1.000 | 1.000 / 1.000 | −0.874 / −0.874 |

Four models show perfectly increasing rank order for both local-patch assays
among observed successes. Larger domains can therefore require longer finite
stimulation for persistent replacement in several established models, which
supports studying polarity as a contributor to persistence. The complete
response is nevertheless stimulus-dependent: fixed-slope Spring becomes
faster with length, Holmes gradient orientation is nearly uncorrelated, and
deBelly is nonmonotonic. Constant slope also increases total/peak input as the
domain grows. A universal rule that enlargement always protects polarity
would contradict the retained evidence.

Some correlations use restricted successful ranges: Holmes fixed-width strength
uses seven lengths from 10–20, while its orientation uses ten from 10–25;
Wave-pinning fixed-width strength uses 16 from 5–31. All 72 eligible deBelly
readouts succeed, but just 13 pass the separate field-tail diagnostic. Its
shortest durations are legitimate scalar results, not verified steady-state
replacement times. These distinctions belong in the legend and supporting
methods rather than being silently filtered out.

## Literature migration and conditional OU evidence

### Detailed proposed integrated-figure legend

**Published migration summaries support distinguishing directional persistence
from migration speed.** A, digitized observed and primary-fit OU-predicted
enlarged/polyploid-group-to-reference folds for nine summaries across Xuan 2018,
Mallin 2023, Jemaà 2017, Jemaà 2023 and the Hosny hypoxic baseline. The group
label WGD+/WGD− is a shorthand across heterogeneous studies, not an assertion
of matched isogenic WGD in every experiment. Xuan's MSD fold is interpolated
at the last lag shared by both digitized population curves, approximately
10.95 h; Hosny values are group medians, and the signed migration index is
distinct from unsigned directness. Directness is net displacement divided by
path length. B, full digitized Xuan population MSD traces and analytical fitted
OU ensemble MSD on log–log axes. The observed trace crossover near 2.32 h is a
descriptive interpolation of published raster curves; the fitted OU crossover
is about 1.11 h and does not reproduce the same crossing. Digitized trace
samples are correlated samples of published lines, not independent cells or
recovered original experimental lags. C, selected primary-fit model-generated
Xuan trajectories, sampled over the indicated time window; these paths are
simulations, not observed cells. D, relative data misfit and regularized total
objective when the enlarged-group persistence/noise parameters are varied with
reference-group parameters held fixed. The heuristic objective is not a
calibrated likelihood; its scale has no biological unit. E, primary fitted
velocity-memory times and velocity-noise amplitudes for each comparison.
These times are conditional OU parameters, not directly measured velocity
autocorrelation times. F, persistence folds under five prior/pooling choices,
including a pooled geometric mean. Scenarios are conditional sensitivity
analyses, not uncertainty intervals. The Jemaà observations constrain speed
alone, leaving persistence/noise tradeoffs supplied by priors and hierarchy.
Full-curve residuals, study-specific acquisition and grouping limitations, and
Hosny drift/leave-one-study-out sensitivity are retained in supporting
material. The figure motivates the persistence hypothesis; it does not
identify a unique polarity mechanism or isolate size as its cause.

Source [caption](../wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU_caption.md),
[observed-fold table](../wgd_langevin/refined_ou_20261003/figures/main_03_observation_folds.csv),
[persistence table](../wgd_langevin/refined_ou_20261003/figures/figure3_persistence.csv),
and [fit results](../wgd_langevin/refined_ou_20261003/results/RESULTS.md).

### Key results, strength and limits

Xuan's enlarged group has about 0.637 times the sampled speed but 2.076 times
the population MSD at the last shared lag. Mallin's PACC group has about 0.519
times the path length and 3.123 times the directness. Hosny's PACC baseline has
about 1.771 times the median path length and 1.673 times the median directness.
Jemaà's comparisons instead report higher speeds. These observations support
an argument centered on sustained displacement and directional organization
rather than a claim that WGD uniformly increases or decreases speed.

The primary inferred OU persistence folds are 9.247 (Xuan), 41.001 (Mallin),
13.836 (Jemaà 2017), 15.153 (Jemaà 2023) and 12.381 (Hosny); the pooled geometric
mean is 15.798. All five conditional sensitivity scenarios retain folds above
one. Leaving Hosny out and refitting yields a pooled fold of 19.71. These
numbers describe the fitted hierarchy. Data-only Jacobian rank is 18/22, while
regularization raises the primary rank to 22/22. Speed alone cannot determine
the Jemaà persistence times.

The Xuan reference fit systematically overpredicts the early digitized curve
and underpredicts the late curve; mean log prediction/observation is +0.267
early and −0.257 late. This exceeds the raster envelope and survives completed
optimization. Thus the fitted OU representation is an interpretable
conditional summary with visible model discrepancy. Raster envelopes are not
population confidence intervals. Stable tetraploid clones, giant polyploids
and treatment-induced PACC states also differ in biology and assay conditions;
their combination does not isolate enlargement as a causal intervention.

The source map is [studies.csv](../literature_review/migration_data/studies.csv).
Retained primary articles include Xuan 2018 (`10.1038/s41598-018-29817-5`),
Mallin 2023 (`10.1007/s10585-023-10216-8`), Jemaà 2017 (Oncotarget article 14291),
Jemaà 2023 (PMCID PMC10531202), and the retained Hosny preprint
`10.1101/2025.09.03.673887`. Hosny's September-2026 journal metadata was identified,
but its full text was not the source used in this extracted baseline; do not
silently cite uninspected final-version measurements as the retained inputs.

## Actionable second-checkpoint migration repair

The 13-Oct [exponent design preview](../wgd_langevin/refined_ou_20261003/figures/design_proposals/figure3B_20261013/figure3B_exponent_proposal.pdf)
is a later proposal, not a replacement of the integrated canonical figure and
not a subsequently user-approved asset. Its PDF SHA-256 is
`a8e5f7b84d34d67e26ebabc120d0bb7daf0c91fbfae9737db8c7c36f730c2615`.
It corrects a previous overly dense dot display: twelve points per group are
explicitly described as chosen display samples of raster curves, while six
lower-panel slopes are newly calculated window estimates. These are not
recovered author-measured lag locations or author-reported local exponents.

The best bounded revision preserves fitted data and the other five panels:

1. Write a new outline-owned migration renderer, leaving the canonical 3-Oct
   scientific campaign and figures immutable. Use the original integrated
   builder's axes/layout, with a two-row subgrid in its B slot.
2. Read the selected `fit_01.json` whose digest is recorded in
   `main_03_migration_OU_provenance.json`. Reuse the existing `figure_tracks.msd_expected`
   expectation. The exact B drawing block is in
   `scripts/build_manuscript_figures.py`, lines 82–113. Replace population
   connecting lines with small, truthfully labelled raster-curve display points;
   reserve continuous lines for the fitted OU expectation, as requested.
3. Reuse `build_exponent_preview.py`'s estimator and identical windows for the
   digitized curve and analytical OU prediction: six equal log-lag windows,
   width approximately 0.304 decades, and matching-window moving model slopes.
   A constant power-law exponent appears horizontal; use reference levels 1
   and 2. This shows curvature and model discrepancy without relabelling the
   fitted OU time as an observed exponent.
4. Retain an observed-crossover indication, explicitly interpolated from the
   full digitized traces, and the 2.32-h number. Do not move it onto the model
   crossover or imply the two coincide. Keep sparse display points separate
   from the full trace used for fitting/crossover calculations.
5. State once in the visible legend: “Points sample digitized published
   population curves; original experimental lag locations are unknown.” Put
   the nonindependence, chosen display count and slope-window details in the
   detailed caption. No sample-size or biological error bars can be computed
   from the selected display points.
6. Verify that selected fit and all original panel numerical inputs retain
   their hashes, that expected MSD at every original fitting lag is unchanged,
   and that known synthetic power laws return their true slopes. Compare
   original/revised panels visually at the final printed size. Record source
   and new output hashes and describe the revision as editorial reassembly,
   without a new fit.

There is **no current basis for choosing dots at actual author-calculated lag
locations**. The extraction [methods](../literature_review/migration_data/xuan_2018/extraction/METHODS.md)
show that the 510 samples came from native raster columns of two lines. The
movie's 10-minute/12-hour versus 16-hour timing discrepancy cannot recover
the published population-MSD lag grid. Inferring a regular grid from presumed
cadence would invent measurement locations. If true author-lag points are
required, they need an independently recovered source dataset first. Until
then, explicit curve-sample labels are the defensible repair.

The original paper does contain a separate mean logarithmic-slope summary in
Fig. 4E, digitized as approximately 1.065 (non-PGCC) and 1.473 (PGCC), with
reported SEM. It could supply an optional independent-looking *reported
summary* inset, carefully labelled as that paper's global statistic; it is not
six local exponents and cannot establish where individual MSD lags were
sampled. Adding it is unnecessary for the first outline checkpoint.

## Outline integration guardrails

Keep the main comparison, the integrated mechanism/measurement component,
existing migration evidence and a reserved substantial new-assay result figure.
The planned SUM-159 chemically induced WGD+/− motility pilot has a collaboration
and starting cell pair, but no results in this evidence set. Its placeholder
must say results pending, with proposed readouts separated from completed
observations. Do not replace that major experimental contribution by the
mechanism-linked measurement proposals.

The conceptual bridge has three different observables: a finite opposing cue
needed for 1D polarity replacement, spontaneous/noisy-ring polarity-vector
persistence if completed, and measured migration directness/MSD or conditional
OU velocity memory. State the hypothesized connection without equating their
numerical values or claiming a calibrated cell-velocity law. Mechanism work
can motivate activity/protrusion, site-fate, lifetime and separation measurements
without already predicting whole-cell movement.

This audit created no large outputs. It changed only this memo, a short README
note and the authorized stale working PNG copy; canonical data and figures,
Git index and existing scientific result directories were left untouched.
