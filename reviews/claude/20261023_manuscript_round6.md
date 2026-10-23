# Review: manuscript_round6 — 2026-10-23

Deliverables reviewed (all at commit `209a0c3`; none of the reviewed files differed in the working tree):
- `manuscript_outline/revision_round6_20261023/review_snapshot_initial/outline_v12.pdf`, pp. 1–30 (summary, Figs 1–6 with legends and text, supplement index, pilot observables table)
- `polarity1d/.../revision_round6_20261023/figure1/figure1_round6.png` and `legend.md`
- `figure2/figure2_round6.png` and `legend.md`
- `figure3/figure3_round6.png` and `legend.md`; `figure3_round6_supplement.pdf`, pp. 1–4
- `figure4/supporting/review_snapshot_initial/figure4_main.png` and `legend.md`; `figure4_supplement.pdf`, pp. 1–3
- Data checked: `figure3/supporting/data/alpha_bridge.csv`
- Baseline for comparison: `manuscript_outline/outline_v11.pdf`, pp. 1–18. It is unchanged between `b347520` and `209a0c3`.

Baseline: `b347520`, given in the request. The prior Claude review `20261023_figures3_4_revision.md` covered Figs 3–4 at `49754ab`, which predates this baseline. Its concern numbers 10–18 carry forward below, re-checked against the round-6 deliverables. New concerns start at 19.

Assumptions:
- **Request content.** The request gave only facts (target, commits, paths). It contained no rationale.
- **Path list vs diff.** The request lists `figure1_supplement.pdf`, `figure2_supplement.pdf` and `figure3_supplement.pdf`, but those files do not exist. The committed files are `figureN_round6_supplement.pdf`, and I assumed those are meant. The diff also adds deliverables that were not listed: the top-level `manuscript_outline/outline_v12.pdf`, `figures1_4_with_legends.pdf`, the per-figure `*_with_legend.pdf` and `*_legend.pdf`, `figure3_round6_new_supplement.pdf`, and the live (non-snapshot) Figure 4 exports. I reviewed the listed snapshot versions and assumed the others are assemblies or live copies of them.
- **Figure 4 is an in-progress snapshot.** Its running header, legend and readme state that the targeted modes campaign was still active at capture, and the conditional panels J–K are absent. I reviewed A–I as delivered. I did not count the missing J–K against the revision.
- **Revision brief.** `request/revision_brief_round6_20261023.md` is signed "@Richard" but sits outside `author/`. As in earlier reviews, I read it only in the rationale pass and cite it under "Codex's position" where it prescribed what I flag. The user adjudicates between the brief and the vision.

Revision scope (from the deliverables and the diff):
- **Figs 1–2:** re-lettered so that each lettered panel holds one plot; no new data.
- **Fig 3:** C–D rebuilt around a spatial susceptibility map g(x), which is new. New panels:
  - pulse tests away from the edge;
  - per-position drift contributions;
  - a Wave-pinning edge time course;
  - size-normalized predicted vs measured drift (shared translation only);
  - an implied-α bridge to Fig 4.
- **Fig 4:** re-laid out as A–I.
  - Closure is now shown at every input, with partial-identification bounds.
  - The heat map is replaced by remainder-share lines.
  - Takeover is shown per model at a single input.
  - Heading-correlation curves and birth rates moved to the supplement.
- **Outline v12:** rebuilt around these changes. Fig 5, the reserved Fig 6 and the pilot table are carried over.

## Summary
Verdict: **more vision-aligned than baseline: mixed. Revision executed in a vision-aligned way: partly.**

Fig 3 makes real progress:
- The reader can now see the susceptibility map that drives the drift prediction.
- The map is validated away from the edge (maximum error 1.69%).
- The drift is shown to come from the two domain edges.
- The main comparison is now size-normalized. This resolves prior Concerns 10, 11 and 15.

The rest of the round loses ground in four places, each a major concern:
1. **Fig 3H drops the key test (Concern 19).** The α bridge is empty for Wave-pinning, the one model where it is not trivial. The CSV predicts α ≈ 2.48–2.87, and Fig 4D shows observed α ≈ 2.5–3.4 in the drift regime. A cutoff of f ≤ 0.001 excludes every Wave-pinning observation.
2. **Fig 4 shows no model heading persistence (Concern 20).** No main panel shows a model ring holding its heading longer. The C(ℓ) curves were moved to the supplement and labelled "superseded", although the figure's title claim rests on them.
3. **Text regressions in Figs 3–4 (Concern 21).** The narrative and "Proposed measurements" text went backwards. Fig 4 now proposes no cellular measurement. Fig 3 lost its WGD+/− signatures and the half-turn photoactivation assay.
4. **Page 1 reads as a status log (Concern 22).** The summary and running header read like a run report, and the plain Fig 3–4 claims from baseline are gone.

Prior Concern 17 (revision-log wording in legends) is worse.

## Concerns
| # | Severity | Criterion | Concern | Status | Codex's position |
|---|----------|-----------|---------|--------|------------------|
| 10 | major | Claim strength / role | The main Fig 3 panel plotted the refuted independent-edge prediction. | **addressed** | Fig 3G shows only the shared-translation prediction; the independent-edge candidate is in the supplement. |
| 11 | major | Presentation | The five-decade 1:1 span did not test size. | **addressed** | 3G plots predicted vs measured D_s/D_s(C_min) only. Wave-pinning sits on 1:1 down to about 0.47. |
| 12 | major | Role (novel) / clarity | The size claim was missing from the Fig 3 title and narrative. | **partly reopened** | The baseline title carried both results ("Enlargement slows both turning routes and selectively suppresses boundary susceptibility"). The round-6 title, "Size-dependent susceptibility links deterministic perturbations to weak-input domain drift", drops the turning-route result of A–B and gives no direction. The narrative bullet says only "changes with size in a model-dependent manner". The per-model folds (0.701 / 0.914 / 0.989 / 0.088) survive only in legend C and in a text box inside panel F. No position recorded. |
| 13 | major | Claim strength | The E(iii) normalization artefact. | **addressed** (at baseline; kept) | 4H uses the native-amplitude convention. R ≈ 0 in the drift regime and rises with f and size. |
| 14 | minor | Claim strength | Takeover within 20τ conflates "slower" with "rarer". | partly addressed | Legend I: "Certification and finite follow-up limit inference about all births." The pilot table cites an 80τ cohort ("persists through 80τ"). The main panel still uses 20τ only. |
| 15 | minor | Presentation | The Fig 3C insets used free y-axes. | **addressed** | The insets are gone. 3C plots C20 and C60 maps on one axis per model. |
| 16 | minor | Mechanism → measurement | The pilot table was not keyed to Fig 2 routes, and the separated photoactivation assay had been dropped. | partly addressed; assay regressed | The table's last column and footnote now relate rows to Fig 2 routes. The baseline Fig 3 text proposed "vary photoactivation offset from straddling to separated placement; measure time to half-turn". Round 6 removed it, so Figs 3A–B again have no assay. |
| 17 | minor → **major** | Density / clarity | Legends carry methods detail and revision-log wording. | **open, worse** | New examples: "Exact preceding…", "All A data and timings are unchanged", "the existing edge pulse", "the historical pulse-normalized variant", "exactly as in the preceding edition", "unchanged round 5 transition inputs", "inherited fixed three-size neighbourhood", "Candidate J–K", and "the strong regime is currently unexplained" (twice in one paragraph). Legend 3E lists 32 four-decimal fractions. The running header on every page reads "REVIEW CANDIDATE; TARGETED MODE CAMPAIGN ACTIVE". |
| 18 | note | Consistency | Fig 3B and C–H use different size grids. | partly addressed | Supplement S1/S4 add B-matched sizes (37.8 / 47.6; Spring 3.15 / 3.97). The main panels still use different grids. |
| 19 | **major** | Narrative progress (bridge) / role (novel) | Fig 3H, implied α vs observed α, is empty for Wave-pinning and Goryachev. Only Holmes (2 points) and Spring (about 10) appear, and both models predict α ≈ 2 because their drift is size-independent. That agreement is the drift null, not a test. The CSV predicts α = 2.48–2.87 for Wave-pinning, size-dependent and essentially f-independent. Fig 4D shows observed Wave-pinning α ≈ 2.5–3.4 at f ≈ 0.003–0.02. The f ≤ 0.001 cutoff removes every Wave-pinning observation, although the same figure's 3G and the pilot table treat f ≤ 0.01 as the drift regime. The result that would explain why Wave-pinning gains more than drift geometry predicts (narrative §4: α 2.5–3) is computed but not shown. | open | Readme: "Wave-pinning and Goryachev have none under the existing shared-lag rule. Their deterministic alpha predictions remain in the CSV; they are not plotted against invented observations." Brief §6 asked for the comparison against "Fig 4B weak-input α" and did not set the f cutoff. |
| 20 | **major** | Role (novel) / claim strength | No main Fig 4 panel shows a model ring holding its heading longer. The baseline C(ℓ) row under the kymographs is now supplement S4.1b, labelled "Superseded". The nulls keep their C(ℓ) curves (B), so the pedagogical nulls get the direct display and the novel model result does not. The title claim ("Enlargement preserves heading") now rests on derived summaries (α in D, tolerance fold in E) and on kymographs whose x-axis is arc position, which cannot be read as heading. | open | Brief §3 layout lists A–K without a model C(ℓ) panel, and "Superseded main panels move to the supplement". |
| 21 | **major** | Narrative progress / mechanism → measurement | **Fig 4 text.** "Narrative points" is cut from three result bullets to one status bullet ("…The independent modes-only campaign is active; this review candidate claims no completed mode-collapse result"). "Proposed measurements" lists no cellular measurement ("Heading persistence and local size exponent, not migration persistence"; "Takeover probabilities are per-model…"; "Existing exemplar Fourier spectra…"). The baseline listed lateral-edge fluctuations, births per perimeter-time, finite-horizon takeover and straddling photoactivation. **Fig 3 text.** The proposals ("Map terminal domain shift per dose…"; "Measure both boundaries…") give no expected WGD+/− direction and no model discrimination. **Effect.** This reopens the narrative gap "Figs 3–4 lack experimental signatures", which baseline had largely closed. Only the pilot appendix still carries the signatures. | open | No position recorded. The brief prescribes figure content, not narrative or measurement text. |
| 22 | **major** | Presentation / vision alignment | The page-1 "Scientific argument" no longer states the Fig 3–4 claims plainly. Fig 3's paragraph is mostly linear-range brackets ("first >10% departures 0.02048 at C=20 … 0.65536"). Fig 4's paragraph says "this review candidate claims no completed mode-collapse result". Baseline page 1 gave readable claims for each figure: enlargement slows both routes; susceptibility folds by model; drift predicted without fitting; heading memory extended; takeover declines with size. The verb mismatch ("supply empirical motivation" vs the vision's "support") is unchanged and is listed under Outstanding. | open | Campaign status is declared deliberately (Fig 4 readme: "Immutable review snapshot…"). |
| 23 | minor | Consistency | The pilot table's cross-references are stale after re-lettering. They read "Figs 3D, 4D", "Fig. 4E(i)", "Fig. 4E(ii)", "Fig. 4B–D" and "Fig. 3C"; round-6 3D is the mean-state profile, 4E is the tolerance fold, and birth rates are now S4.1. Its "Model drift subset: f ≤ 0.01 and multiple-domain fraction ≤ 0.05" no longer matches main 4G, which covers all inputs. | open | Round-6 README: the "prospective pilot appendix … preserved". Brief: "Update legends and panel callouts to match." |
| 24 | minor | Clarity | **Fig 2.** Legend G says histories are "labelled kept / lost / lost then commits", but the round-6 artwork has no such labels (baseline E did). H–I name classes ("Initiation-limited", "Saturating activation", "Slow non-activity state") but not models, so the class-to-model mapping now exists only in the removed G labels. Panel I has no WGD−/WGD+ key. There is a large blank band between the E line-style key and the E plots. | open | Brief: "plain black facet labels", "re-letter and re-lay out only". |
| 25 | minor | Role / density | **Fig 3E.** Four facets of near-identical Gaussians, plus 32 numbers in the legend, convey one fact: about 99.6–100% of predicted drift comes from within two correlation lengths of the two edges. That fact is the mechanism → measurement point (measure lateral-edge fluctuations, not interior ones), but no text states it as a result. One sentence plus one facet would serve. **Fig 3F.** Carries a prose text box inside the panel area. | open | Brief §3 asked for E as its own panel with fractions in the legend. |
| 26 | minor | Presentation | **Fig 1.** Re-lettering turned the preparation exemplar into three lettered panels (A–C) across the top third. Panel A is mostly empty, and all three show *establishing* polarity, not the challenge the title is about. The shortest-pulse determination and the magnitude of the increase are still not shown (narrative §1). | content outstanding (out of scope); prominence new | Brief: Figs 1–2 style only. |
| 27 | minor | Clarity | **Fig 3 legibility.** The A key is cropped at the top edge, A's "50%" and "90%" labels overlap the peaks, the 3E legends sit on the curves, and 3B's light size labels are faint. **Fig 4 legibility.** In 4G, "observed rate unavailable" overlaps the crosses. **Axis ranges.** Goryachev's different 4F ordinate (10⁻⁶–1) is no longer declared in the legend; baseline legend D declared it. **Typography.** The Fig 4 title has a trailing period. | open | — |
| 28 | note | Consistency | The Fig 4 snapshot readme reports that exact own-mask moments reduce \|R\| "in only 7/325" comparisons and, later, "in 321/325". These are different denominators or conventions, but a reader of the supplement will see a contradiction. | open | The readme distinguishes "matched-denominator" from the "round 5 amplitude-scaled point convention". |

## Outstanding
These are outside this revision's scope or untouched by it, and are not counted against it.
- **Magnitude gap.** At α ≈ 2, one WGD gives about 1.6-fold, against a published OU fold of about 16. Not shown. The Wave-pinning predicted α of about 2.5–2.9 would raise the figure only modestly, to about 1.8–1.9-fold per WGD, which is worth stating once Concern 19 is fixed.
- **Fig 1.** Duration magnitude and the shortest-pulse curves are not shown in the main figure (scope: style only).
- **Fig 2.** Three items remain open: the mechanistic basis of the two scaling modes (candidate net-activation maps), deBelly's retained mechanical state, and the size contrast in the 2H schematics (2G is still L39 only).
- **Fig 3B.** Fold change per volume doubling is not plotted.
- **Strong-input regime.** The legend says plainly that it is unexplained, pending the mode campaign (J–K). This is honest and appropriate for now.
- **Verb.** Page 1 says "supply empirical motivation"; the vision says "support". Unchanged.
- **Supplement curation.** The Fig 3 supplement appends the whole round-5 supplement, and Fig 4 S4.1b is labelled "Superseded". A curation pass is needed before submission.

## Narrative targets
| Narrative item | Section | Status | Evidence in deliverable |
|---|---|---|---|
| Shortest pulse per length (not shown) | 1 | not addressed (out of scope) | Fig 1 E is still rank correlations only |
| Magnitude of the duration increase | 1 | not addressed (out of scope) | — |
| Mechanistic basis of the scaling modes | 2 | not addressed (out of scope) | — |
| deBelly retained state | 2 | not addressed (out of scope) | — |
| Size contrast in 2F (now 2H–I) | 2 | not addressed | G still L39 only |
| 3C candidate: edge response to a weak pulse, at each size | 3 | **addressed, extended** | 3C full spatial map at C20/C60, validated off-edge (1.69%) |
| Fold change per volume doubling | 3 | not addressed | — |
| Two regimes ↔ two ways of losing heading (speculative) | 3 | not addressed | — |
| Resolution of coexisting domains after withdrawal | 3 | partly, unchanged | 3A "cue off" row |
| Predicted vs measured drift (candidate panel) | 4 | **addressed** | 3G, size-normalized, shared translation, no fit |
| Drift as small-signal response; linear range | 4 | addressed (text) | Linear-range brackets in legend and page 1 (Concern 22 on placement) |
| Holmes/Spring "surprising" vs Fig 3 sliding | 4 | resolved | Susceptibility, not strong-cue sliding, sets drift |
| Wave-pinning α > 2 explained by size-suppressed drift | 4 | **computed, not shown** | Concern 19 |
| Births per perimeter; takeover by distance and size | 4 | partly | 4I per model, two distance bins; births moved to S4.1 |
| Exact decomposition of heading change | 4 | addressed | 4H remainder; S4.2 shares |
| Heading loss in the multiple-domain regime "not identified" | 4 | stated as unexplained | 4G bounds at all f; legend: "strong-input regime is currently unexplained" |
| Detector-threshold sweep | 4 | not addressed | — |
| Experimental signatures for Figs 3–4 | gaps | **regressed in figure text**; kept in pilot table | Concern 21 |
| Imposed → noise bridge | gaps | **advanced** for drift (3C → 3G), but the α step is missing (Concern 19) | 3C, 3E, 3G |
| Magnitude vs published fold | gaps | not addressed | — |

## Narrative claims
| Claim | Figure | Still supported? | Note (regression vs stale narrative) |
|---|---|---|---|
| Required duration rises with length (fixed and relative widths); ramp model-specific | 1E | supported | Stale narrative: panel now 1E, not 1C |
| Two scaling modes (broaden vs heighten) | 2A–C | supported | Stale: split into A–C |
| Goryachev initiation slows about 2.5-fold; others length-independent | 2D | supported | Stale lettering |
| Replacement 5–12-fold slower L19 → L39 | 2E | supported | — |
| Holmes narrow-cue 45-fold penalty | 2E | **qualified** | The legend says censoring "does not establish a growing penalty"; the narrative should soften this (unchanged since baseline) |
| Pulse-train response classes | 2F–G | supported | 2G labels now missing (Concern 24) |
| Two turning regimes; larger rings slower at every cue position | 3A–B | supported | Data unchanged; dropped from title (Concern 12) |
| Distance matters little beyond the regime boundary | 3B | supported | — |
| α ≈ 2 (Holmes, Spring), 2.5–3 (Wave-pinning), declining toward 1 | 4D | supported | Stale lettering (was 4C) |
| Tolerance folds 2.2 / 2.4 / 3.1 / 4.1 | 4E | supported | — |
| Larger ring holds heading longer near the transition | S4.1b only | **no longer shown in main** | Regression (Concern 20) |
| Size reduces drift in Wave-pinning (about 2×), strongly in Goryachev, not in Holmes/Spring | 4F; 3C | **strengthened** | Independently predicted by the 3C susceptibility folds |
| "If edges respond independently… edge pulse predicts drift" | 3G | premise refuted; shared-translation prediction succeeds | Suggested narrative update: drift is common-mode translation generated at the two edges (3E) |
| Drift + relocations within 15% in the drift regime | 4G | supported | Shown at f ≲ 0.01 within the all-input panel |
| 0/33 multiple-domain conditions within 15% | 4G | withdrawn | Bounds span unity at strong input; the legend calls the regime unexplained |
| Secondary domains spawn readily; incumbents rarely lose | 4I; S4.1 | partly supported | Takeover declines with size in all four models at one input; the finite-horizon caveat remains (Concern 14) |
| Heading lost via balance among coexisting domains (speculative) | 4H; S4.2 | untested | R rises with f and size; membership share rises (S4.2) |

## Result roles
| Result | Role(s) | Presented to maximize role? |
|---|---|---|
| Fig 1 preparation (A–C) | pedagogical | No: three lettered panels for setup (Concern 26) |
| Fig 1E duration correlations | novel | As at baseline |
| Fig 2A–G | novel; hypothesis-generating | Mostly; 2G lost its class labels (Concern 24) |
| Fig 2H–I schematics | hypothesis-generating | Partly: classes are not tied to models in the artwork |
| Fig 3A–B turning routes | novel; hypothesis-generating | Less than at baseline: not in the title, and the assay was dropped (Concerns 12, 16) |
| Fig 3C susceptibility map | **novel; hypothesis-generating** | Yes visually. Direction and folds are not stated in the narrative (Concern 12) |
| Fig 3E edge localization of drift | novel; hypothesis-generating | No: the result is not stated, and it is over-plotted (Concern 25) |
| Fig 3G size-normalized drift prediction | novel (bridge) | Yes |
| Fig 3H implied α | novel (bridge) | No: the informative model is empty (Concern 19) |
| Fig 4A–B nulls | pedagogical | Over-weighted: two rows, with the C(ℓ) display the models lack (Concern 20) |
| Fig 4D–F α, tolerance, drift | novel | Yes |
| Fig 4G closure bounds | diagnostic (where the domain picture works) | Acceptable; dense |
| Fig 4H remainder | diagnostic | Yes; artefact fixed |
| Fig 4I takeover by size | novel; hypothesis-generating | Partly: no takeaway sentence; 20τ caveat |
| Pilot observables table | hypothesis-generating | Mostly; stale references (Concern 23) |

## Detail

**Concern 19: Fig 3H.** From `alpha_bridge.csv`:
- **Wave-pinning.** Predicted α rises from 2.48 (C20 neighbourhood) to 2.87 (C60) and is identical to three decimals at f = 10⁻⁴, 3×10⁻⁴ and 10⁻³, as linear response requires. Every observed row reads `unavailable_no_shared_lag…`.
- **Fig 4D.** Wave-pinning has qualified α at f ≈ 0.003–0.02, with values about 2.5–3.4, where 3G and the pilot table say the drift regime holds.

Because the prediction does not depend on f, it can be compared directly with those observations. Two ways to do so:
- extend H to f ≤ 0.01 (or to the drift-regime subset already used in the pilot table);
- overlay the predicted α band on the Wave-pinning facet of 4D.

Either would show the round's most distinctive result: in Wave-pinning, enlargement suppresses edge susceptibility, which predicts α above the drift null, which is what is observed. Holmes and Spring α ≈ 2 is the drift-null expectation for size-independent susceptibility and cannot by itself test the bridge. The choice of f window is Codex's; the brief did not set one.

**Concern 20: Fig 4 heading display.** The null facets in B show the C20 vs C60 decorrelation a reader needs to see. The models get no such facet in the main figure. S4.1b (Wave-pinning f = 0.03, Holmes 0.02, Spring 0.005, Goryachev 0.01) is exactly that display, and it is clear: at C60 the heading correlation stays at 0.45–1.0 over the full lag window, while at C20 it is lost within about 3×10³. Suggested fix: restore those four curves as their own lettered panel, after C or merged as a facet row. This satisfies the brief's one-plot-per-letter rule, and the title claim then has a direct display. The brief's layout omitted the panel; whether to add it is the user's call.

**Concern 21: figure text.** Baseline (v11, p. 16) Fig 4 "Proposed measurements" read: "Measure lateral-edge fluctuations, secondary-protrusion births per actual perimeter-time, finite-horizon takeover probability and leading-edge displacement after boundary-straddling photoactivation." Round 6 (v12, p. 18) replaces this with three analysis statements. Suggested Fig 3 text, built from round-6 data:
- **Edge displacement per photoactivation dose:** smaller in WGD+ cells under Wave-pinning-like (×0.70 over a 3-fold perimeter) or Goryachev-like dynamics; unchanged under Holmes-like or Spring-like dynamics.
- **Where to measure fluctuations:** at the lateral edges of the leading region, which generate essentially all of the drift (3E). Interior fluctuations should not predict heading drift.
- **Restore the straddling vs separated half-turn assay** for 3A–B.

**Concern 22: page 1.** Suggested direction: one or two plain sentences per figure, as in baseline. The linear-range brackets belong in the Fig 3 legend or the supplement. Campaign status belongs in the Fig 4 scope section, not in the argument or the running header.

**Fig 3 progress (Concerns 10, 11, 15).** Verified in the deliverable:
- 3C shows C20 and C60 maps on one axis per model. The existing edge pulses (open squares) and the new off-edge tests (circles) land on the lines.
- 3F shows both Wave-pinning edges settling to the same permanent offset.
- 3G Wave-pinning points lie on 1:1 from 1.0 down to about 0.47. Holmes clusters at 0.8–1.0. Spring's prediction is flat at about 1, while measured values spread from 0.85 to 1.5. Goryachev fails on its own axis and is labelled "model fails".
