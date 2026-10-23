# Review: manuscript_round6 (follow-up) — 2026-10-23

Deliverables reviewed (all at commit `028d9f5`; the working tree was clean at launch, so the working-tree files I viewed are the committed versions):
- `manuscript_outline/outline_v12.pdf`, pp. 1–32. I read pp. 1–15, 17–32; p. 16 is the Fig 4 artwork, which I reviewed from the PNG.
- `manuscript_outline/revision_round6_20261023/pilot_observables.pdf`
- `polarity1d/.../revision_round6_20261023/figure1/figure1_round6.png` and `legend.md`
- `figure2/figure2_round6.png` and `legend.md`
- `figure3/figure3_round6.png` and `legend.md`
- `figure4/figure4_main.png` and `legend.md`; `figure4_supplement.pdf`, pp. 1–16

Baseline: `209a0c3`, given in the request. This is also the commit of the prior Claude review of this target (`20261023_manuscript_round6.md`). Concern numbers 10–28 carry forward from that review. New concerns start at 29.

Assumptions:
- **Request content.** The request gave only facts (target, commits, paths). It contained no rationale.
- **Path list vs diff.** Every listed path appears in `git diff --stat 209a0c3 028d9f5`. The diff also changes deliverables that were not listed:
  - `figureN_legend.pdf`, `figureN_with_legend.pdf` and `figures1_4_with_legends.pdf`;
  - `figure3_round6_new_supplement.pdf`;
  - `pilot_observables.png`;
  - `figure4/supporting/vectors/figure4_A_J.pdf`.

  I assumed these are assemblies or alternative renderings of the listed files, and did not review them separately.
- **Supplements.** I did not reopen the Fig 1–3 supplements. Their diffs are small, and none of this round's concerns rests on them.

Revision scope: a response to the prior review, with no new simulations in Figs 1–3.
- **Fig 2:** presentation fixes only. Class labels are restored on G, H and I name the models, and I gains a WGD key.
- **Fig 3:**
  - the title is restored to the size claim;
  - H is rebuilt with matched-input observations at f ≤ 0.01, so Wave-pinning is now included;
  - the legend and text are rewritten;
  - the cellular measurements are restored.
- **Fig 4:**
  - a new panel D shows model heading autocorrelation;
  - the later panels are re-lettered E–J;
  - the legend and text are rewritten;
  - the mode campaign is completed, and its Fourier tests (labelled K–L) stay in the supplement.
- **Outline v12:** a new page-1 argument, the running header fixed, and a new pilot table with updated references.

## Summary
Verdict: **more vision-aligned than baseline: yes. Revision executed in a vision-aligned way: yes, with minor residue.**

All four major concerns from the prior review are resolved in the deliverables, not just in the response note:
- **Concern 19 (Fig 3H).** It now shows the non-trivial test. In Wave-pinning, enlargement suppresses edge susceptibility, which predicts α of 2.48–2.88, above the drift null of 2. The 15 observed Wave-pinning α values (about 2.4–3.4) scatter around the 1:1 line.
- **Concern 20 (Fig 4D).** The figure now shows directly that the larger model ring keeps its heading longer.
- **Concern 21 (measurements).** Figs 3–4 again propose cellular measurements. They give WGD+/− directions and say which model classes each measurement would and would not separate.
- **Concern 22 (page 1).** Page 1 now states scientific claims, and uses the vision's verb "support" for the migration data.

What remains is minor:
- K–L lettering on supplement-only panels, and an inconclusive Fourier screen that is repeated in the Fig 4 conclusions (29);
- one assay that drops Goryachev from the models it cannot separate (30);
- page-1 gaps for Figs 1 and 5 (31);
- some remaining methods density in the legends (17, 32).

No major concerns are open.

## Concerns
| # | Severity | Criterion | Concern | Status | Codex's position |
|---|----------|-----------|---------|--------|------------------|
| 10 | major | Claim strength / role | Refuted independent-edge prediction in main Fig 3 | **addressed** (kept) | — |
| 11 | major | Presentation | Absolute drift span did not test size | **addressed** (kept) | — |
| 12 | major | Role (novel) / clarity | Size claim missing from Fig 3 title and narrative | **addressed** | Title: "Enlargement slows both turning routes and selectively suppresses edge susceptibility". The narrative bullet gives per-model folds and direction. Page 1 repeats it. |
| 13 | major | Claim strength | E(iii) normalization artefact | addressed (kept) | — |
| 14 | minor | Claim strength | 20τ takeover conflates "slower" with "rarer" | partly addressed (unchanged) | Legend J: "Lower 20-time-scale incidence can reflect slower or rarer takeover; lifetime takeover probability is unidentified." The 80τ cohort appears only in the pilot table. Codex: no new follow-up simulations. |
| 15 | minor | Presentation | Free-axis insets in 3C | addressed (kept) | — |
| 16 | minor | Mechanism → measurement | Pilot table not keyed to routes; separated photoactivation assay dropped | **addressed** | Fig 3 text (p. 14) and Fig 4 text (p. 19) restore the straddling vs separated half-turn assay, with the prediction that WGD+ turns more slowly by both routes. |
| 17 | major → **minor** | Density / clarity | Methods detail and revision-log wording in legends | **partly addressed** | **Fixed:** revision-log phrases are gone from the Fig 2–4 legends, and the running header is now scientific. **Remaining:** the Fig 1 legend ("deprecated in-loop centered extension", "not silently converted"); Fig 4 legend J (selection rule, 2000 resamples, 97.5% rule); a trailing amplitude-bracket paragraph in Fig 3 (see 32); page 2 ("retain their numerical evidence… with updated panel references"); the Fig 5 sources line ("direct-user feedback"); the supplement index ("earlier main artwork is preserved"). |
| 18 | note | Consistency | 3B and 3C–H use different size grids | partly addressed (unchanged) | Codex: the B-matched sizes are supplementary, and no claim of identical grids is made. |
| 19 | major | Narrative progress / role (novel) | 3H was empty for Wave-pinning | **addressed** | 3H shows 15 Wave-pinning points (f = 0.003–0.01), predicted 2.48–2.88 against observed about 2.4–3.4, around 1:1. The legend states "above the size-independent-drift expectation of 2". The Goryachev facet is empty and labelled as such. |
| 20 | major | Role (novel) / claim strength | No main panel showed model heading persistence | **addressed** | New 4D: smallest vs largest ring at the C inputs. In all four models the larger ring decorrelates more slowly. S4.1b is retitled by content, and "superseded" is gone. |
| 21 | major | Narrative progress / mechanism → measurement | Figs 3–4 text lost cellular measurements and WGD signatures | **addressed** | **Fig 3 (p. 14):** three assays with WGD+/− directions by model class. They include an edge-vs-interior fluctuation test and a caveat that the folds are for 3× circumference, not one WGD. **Fig 4 (p. 19):** five assays, covering heading autocorrelation, lateral-edge tracking, births per perimeter-time and takeover by distance, photoactivation offset, and Fourier spectra. See 30 for one discrimination gap. |
| 22 | major | Presentation / vision alignment | Page 1 read as a status log | **addressed** | Plain claims for Figs 2–4. "Existing migration observations support this hypothesis". The campaign status is removed. See 31 for residue. |
| 23 | minor | Consistency | Stale pilot-table references | **addressed** | The new table cites 3F–G, 4G, 4J, 4D–F, 3C/F and the certified-birth supplement. The footnote says correctly that 4H shows all inputs. |
| 24 | minor | Clarity | Fig 2 labels, model mapping and WGD key | **addressed** | G has kept / lost / lost then commits. H titles name the models. I has its own WGD−/WGD+ key. A blank band still separates the E key from the plots (cosmetic). |
| 25 | minor | Role / density | 3E over-plotted, result unstated; prose box in 3F | partly addressed | Legend E now leads with the result ("99.6–100%… suggesting lateral-edge fluctuations as the cellular readout"). The shares are rounded, and the 3F box is removed. The four near-identical facets remain. Codex: the brief asked for the facets. |
| 26 | minor | Presentation | Fig 1 setup prominence; duration magnitude not shown | prominence partly addressed; content outstanding | A–C still occupy the top third, and A is mostly white space. Codex: the preparation row is smaller and the duration magnitude is out of scope. |
| 27 | minor | Clarity | Cropping, overlaps, faint labels, undeclared axes | **addressed** | Verified: the 3A key and threshold labels are clear; the 3E legends are off the curves; the 4H unavailable-rate crosses are at the top; legend G declares the Goryachev ordinate; the title has no trailing period. |
| 28 | note | Consistency | Exact-mask remainder counts (7/325 vs 321/325) | addressed in text; residual note | Codex: the two counts belong to different conventions; full masks reconstruct activity algebraically. S4.5 now shows the spread. Across conventions the remainder ranges from about −2 to +0.7 and changes sign; under the exact native-total masks it is about 0 by construction. Legend I says only that the decomposition is "descriptive". Add that the magnitude and sign of R depend on convention, and why point/native-amplitude is the main one. |
| 29 | minor | Presentation / role | **Fourier panels.** The supplement labels them "K" and "L", using main-figure letters, and the main legend refers to "Fourier-mode tests (K–L)", although the main figure ends at J. The inconclusive strong-regime screen ("partial agreement… does not establish a robust strong-regime explanation") appears twice: as a legend paragraph and as a fourth "Narrative points and conclusions" bullet. A diagnostic with a null result is occupying conclusion space. | open | Codex: "Conditional Fourier panels are K–L"; "The Fourier result remains supplementary under the declared descriptive promotion rule." |
| 30 | minor | Mechanism → measurement | **Fig 4 lateral-edge assay.** The text (p. 19) says "Size-suppressed physical drift would favour a Wave-pinning-like susceptibility change". Goryachev also has strongly size-suppressed drift (4G; edge-gain fold 0.088), and the Fig 3 text and pilot table both say "Wave-pinning/Goryachev". The vision asks each measurement to say which mechanisms it fails to separate. Here that is Wave-pinning-like vs Goryachev-like, which the Fig 2 first-rise assay (Goryachev initiation slows with size) would then separate. | open | No position recorded. |
| 31 | minor | Presentation / narrative order | **Page 1 argument.** It has no Fig 1 sentence and no Fig 5 sentence; Fig 5 appears only through the "support" clause in the hypothesis paragraph. The Fig 2 paragraph omits its headline size result (replacement 5–12-fold slower from L19 to L39) and spends its second sentence on the Holmes censoring caveat. The Fig 2 "Narrative points" have the same gap. **Layout.** Page 6 is blank apart from one link. | open | Codex: Figs 1–2 are style-only in this round. |
| 32 | note | Clarity | **Fig 3 legend placement.** The final paragraph ("At the right edge, last tested linear peak amplitudes…") has no clear referent and belongs in methods. The all-model folds (0.701 / 0.914 / 0.989 / 0.088) sit in legend F, which describes a Wave-pinning-only panel. They belong with C, where the maps are shown. | open | Codex (17): amplitude brackets moved to methods. This paragraph remains in the legend. |

## Outstanding
Outside this revision's scope or untouched by it.
- **Magnitude gap.** A single WGD (1.26-fold in circumference) still gives about 1.6-fold at α = 2, or about 1.8–1.9-fold at the Wave-pinning α of about 2.5–2.9, against a pooled OU fold of about 16. Not shown. The Fig 3 text now notes that its folds are for 3× circumference. A per-WGD conversion of the predicted edge-displacement change would help judge feasibility. For Wave-pinning it is roughly 0.70^(ln1.26/ln3) ≈ 0.93 per WGD, a small effect for a cellular assay.
- **Strong-input regime.** Stated plainly as unexplained. This is appropriate.
- **Fig 1:** shortest-pulse curves and duration magnitude are not in the main figure.
- **Fig 2:** three items remain open:
  - the mechanistic basis of the two scaling modes;
  - deBelly's retained mechanical state;
  - the size contrast in the G histories (L39 only).
- **Fig 3B:** fold change per volume doubling is not plotted.
- **Fig 5:** two titles disagree. The image page reads "Migration observations and conditional OU interpretation"; the text page reads "Published migration summaries support…".
- **Detector-threshold sweep:** not done.
- **Supplement curation** before submission. The Fig 4 supplement mixes S4.x and K–L labels, and the Fig 3 supplement still appends the round-5 material.

## Narrative targets
| Narrative item | Section | Status | Evidence in deliverable |
|---|---|---|---|
| Shortest pulse per length | 1 | not addressed (out of scope) | — |
| Magnitude of the duration increase | 1 | not addressed (out of scope) | — |
| Mechanistic basis of the scaling modes | 2 | not addressed (out of scope) | — |
| deBelly retained state | 2 | not addressed (out of scope) | — |
| Size contrast in 2F (now 2H–I) | 2 | partly | Schematics now carry the WGD− / WGD+ key and model names; G is still L39 only |
| 3C candidate: edge response to a weak pulse, at each size | 3 | addressed (kept) | 3C maps at C20 and C60 |
| Fold change per volume doubling | 3 | not addressed | — |
| Two turning regimes ↔ two ways of losing heading (speculative) | 3 | partly | Assays now pair straddling vs separated photoactivation (3A–B) with edge fluctuations and births (Fig 4). The mapping is not stated as a hypothesis. |
| Predicted vs measured drift | 4 | addressed (kept) | 3G |
| Wave-pinning α > 2 explained by size-suppressed drift | 4 | **addressed** | 3H, 15 points around 1:1 |
| Larger ring holds heading longer | 4 | **addressed** | 4D |
| Births per perimeter; takeover by distance and size | 4 | partly (unchanged) | 4J; S4.1, S4.3, S4.4 |
| Exact decomposition of heading change | 4 | addressed | 4I; S4.2; S4.5 (see 28) |
| Heading loss in the multiple-domain regime "not identified" | 4 | stated as unexplained | Legend and text |
| Detector-threshold sweep | 4 | not addressed | — |
| Experimental signatures for Figs 3–4 | gaps | **addressed** | Fig 3 and Fig 4 measurement sections, pilot table |
| Imposed → noise bridge | gaps | **advanced** | 3C → 3E → 3G → 3H (drift and α, in the drift regime only) |
| Magnitude vs published fold | gaps | not addressed | — |

## Narrative claims
| Claim | Figure | Still supported? | Note (regression vs stale narrative) |
|---|---|---|---|
| Required duration rises with length; ramp model-specific | 1E | supported | Stale lettering (1C → 1E) |
| Two scaling modes | 2A–C | supported | Stale lettering |
| Goryachev initiation slows about 2.5-fold; others length-independent | 2D | supported | — |
| Replacement 5–12-fold slower from L19 to L39 | 2E | supported | Not stated in the outline text (31) |
| Holmes 45-fold narrow-cue penalty | 2E | qualified | The legend says censoring "does not establish a growing penalty". Suggested narrative update. |
| Pulse-train classes | 2F–G | supported | Labels restored |
| Two turning regimes; larger rings slower everywhere | 3A–B | supported | Now in the title |
| Distance matters little beyond the regime boundary | 3B | supported | — |
| α ≈ 2 (Holmes, Spring), 2.5–3 (Wave-pinning), declining | 4E | supported | Stale lettering (4C → 4E) |
| Tolerance folds 2.2 / 2.4 / 3.1 / 4.1 | 4F | supported | Text gives 2.16 / 2.37 / 3.09 / 4.12 |
| Larger ring holds heading longer near the transition | 4D | **supported** (restored) | — |
| Drift reduced in Wave-pinning (about 2×), strongly in Goryachev, not in Holmes or Spring | 4G; 3C | strengthened | Predicted by the susceptibility folds |
| "If edges respond independently… edge pulse predicts drift" | 3G | premise refuted; shared translation succeeds | Suggested update: drift is common-mode translation generated at the two edges (3E–F) |
| Wave-pinning α > 2 means "more than drift geometry predicts" | 3H | **now explained** | Suggested update: α above 2 is predicted by size-suppressed edge susceptibility |
| Drift plus relocations within 15% in the drift regime | 4H | supported | — |
| 0/33 multiple-domain conditions within 15% | 4H | withdrawn | Replaced by identification bounds |
| Secondary domains spawn readily; incumbents rarely lose | 4J; S4.1, S4.3 | partly supported | 20τ caveat (14) |

## Result roles
| Result | Role(s) | Presented to maximize role? |
|---|---|---|
| Fig 1 preparation (A–C) | pedagogical | Partly; still prominent (26) |
| Fig 1E | novel | As before |
| Fig 2A–G | novel; hypothesis-generating | Yes, labels restored |
| Fig 2H–I | hypothesis-generating | Yes: models named, WGD key present |
| Fig 3A–B | novel; hypothesis-generating | Yes: in the title, with the assay restored |
| Fig 3C, F susceptibility | novel; hypothesis-generating | Yes; the folds sit in legend F rather than C (32) |
| Fig 3E edge localization | novel; hypothesis-generating | Mostly: result stated, measurement proposed; the facets are redundant (25) |
| Fig 3G–H drift and α bridge | novel (bridge) | **Yes** |
| Fig 4A–B nulls | pedagogical | Acceptable; the models now have the matching display (4D) |
| Fig 4D–G | novel | Yes |
| Fig 4H–I closure, remainder | diagnostic | Acceptable; I needs a note on convention dependence (28) |
| Fig 4J takeover | novel; hypothesis-generating | Partly (14) |
| Fourier screen (supplement K–L) | diagnostic, inconclusive | Over-promoted into the conclusions (29) |
| Pilot table | hypothesis-generating | Yes |

## Detail

**Concern 19, verified.** 3H Wave-pinning: x spans 2.48–2.88 and y spans about 2.4–3.4, with intervals crossing 1:1. Holmes clusters at about 2.1–2.2 on x and Spring at about 1.9–2.1, which is the drift-null expectation. Within Wave-pinning the predicted range is narrow, so the panel tests the level of α (above 2) more than its trend. A one-line caveat would help, for example "the comparison tests the excess over 2, not the within-model slope".

**Concern 20, verified.** In 4D, the C60 dotted curves stay above 0.45 (Wave-pinning), 0.1 (Holmes) and about 0 at the end of the Spring window. The C20 curves fall below 0.1 within about 5×10³. The Goryachev C60 curve is flat at 1 (pinned). The legend states the conclusion: "These direct curves show longer heading memory in larger rings."

**Concern 29.** Suggested fix: relabel the supplement panels S4.7–S4.8 and remove "(K–L)" from the main legend. Keep one sentence on the strong-regime status in Scope/limitations, and delete either the legend paragraph or the narrative bullet. The Fourier-spectrum assay on p. 19 tests a detector-free description; it does not separate size mechanisms. Saying so would match the vision's mechanism → measurement standard.

**Concern 30.** Suggested wording: "Size-suppressed physical drift would favour Wave-pinning-like or Goryachev-like susceptibility. The two are separated by the Fig 2 first-rise assay, since initiation slows with size only in Goryachev."

**Concern 31.** Suggested additions: one sentence for Fig 1, for example "Under fixed- and proportional-width cues, the duration needed for lasting replacement rises with length in all five models." For Fig 2, lead with "replacement slows 5–12-fold from L19 to L39 in every model", before the class signatures.
