# Review: manuscript_round7 — 2026-10-26

Deliverables reviewed (all at commit `f815a87`; the working tree was clean at launch with HEAD = `f815a87`, so the files I viewed are the committed versions):
- `manuscript_outline/outline_v13.pdf`, pp. 1–29 (all pages read; artwork pages viewed as rendered in the PDF)
- `polarity1d/.../dev/revision_round7_20261026/figure3/main_legend.md`; Figs 3–6 as they are embedded in outline v13
- `figure3/pilot_conduction_penalty.png`, `pilots/holmes_mode_rates_native_fast_seed_CI.png`
- I did not open `revision_round7_20261026/index.html`. It is a 4-line viewer change.

Baseline: `122d306`, given in the request. It is newer than the narrative's baseline and than the latest prior review (`20261023_manuscript_round6_2.md`, which reviewed `028d9f5`). Concern numbers 14–32 carry forward from that review. New concerns start at 33.

Assumptions:
- **Request content.** The request gave only facts (target, commits, paths). It contained no rationale.
- **Path list vs diff.** Every listed path appears in `git diff --stat 122d306 f815a87`. The diff also changes deliverables that were not listed:
  - `figureN_legend.pdf`, `figureN_with_legend.pdf` and `figures1_6_with_legends.pdf`;
  - `figure6/figure6B_detail_round7.{png,pdf}`;
  - `figure5/noise_legend.md`, `noise_readme.md` and `supporting/viewers/noise_panels_EF*.png`;
  - the `request/` mockups.

  I assumed these are assemblies, alternative renderings or inputs, and reviewed the figures through outline v13.
- **Round-7 request documents.** `manuscript_outline/revision_round7_20261026/request/brief_round7.md` and `target_narrative_round7.md` carry the user's byline ("@Richard") but sit outside `author/`. I read them only in the rationale pass. I treated them as the user's statement of this round's scope, and the target narrative as the narrative the outline text was asked to follow. `author/vision.md` stays authoritative, and `author/narrative_claims.md` (round 3) is the authorised baseline narrative for the claim tables below.
- **No supplement.** The brief says "There is no supplement. Drop dropped panels". Panels removed since round 6 (noise tolerance, remainder, takeover by distance, edge share) are therefore treated as user-directed removals, not regressions.

Revision scope: Figs 1–2 unchanged. Figs 3–5 rebuilt around two proposed persistence mechanisms, cytosolic conduction and Goryachev jump suppression:
- **Fig 3:** turning, plus a susceptibility map, a linearity check, gain against size with a conduction law, and a D_v sweep.
- **Fig 4:** null baselines, model heading correlation at weak and strong input, kymographs and local α.
- **Fig 5 (new mechanisms figure):** drift and α bridge, closure, raw heading-change distributions, the strong-input patch field, and native vs fast-cytosol mode rates.

Old Fig 5 becomes Fig 6. The outline text should follow the target narrative.

## Summary
Verdict: **more vision-aligned than baseline: mixed. Revision executed in a vision-aligned way: partly.**

**What the figures gain.** They now carry a physical, perturbation-tested mechanism, which is exactly what the vision asks for ("suggest a physical basis"):
- **3F:** fast cytosol removes Wave-pinning's size-dependent boundary gain.
- **3E:** an isolated-front conduction law predicts the Wave-pinning ring gains without fitting.
- **5F:** fast cytosol removes most of Wave-pinning's strong-input size advantage.

The split between persistence (Fig 4) and mechanism (Fig 5) is also cleaner than round 6.

**What the text loses.** Three problems:
- The outline text largely does not state these results, and does not follow the target narrative's three-route synthesis (33).
- The cellular measurements that round 6 tied to WGD+/− predictions and model classes have collapsed to generic boilerplate (34). The new cytosolic-mobility measurement gets no predicted signature.
- The Goryachev jump-suppression half of the round's thesis is not shown. Its key panel (5D) cannot be read for its intended signature (35).

Three major concerns are open.

## Concerns
| # | Severity | Criterion | Concern | Status | Codex's position |
|---|----------|-----------|---------|--------|------------------|
| 14 | minor | Claim strength | 20τ takeover incidence conflates slower with rarer | moot (panel 4J dropped by the brief) | — |
| 17 | minor | Density / clarity | Methods detail and revision-log wording in the legends and front matter | **open; partly regressed** | See 37. The running header reads "ROUND 7: CYTOSOLIC CONDUCTION AND HEADING RELOCATIONS" again. Page 1 says "Figures 1–2 retain the preceding evidence unchanged" and "Figure 6 is the preceding migration Figure 5, renumbered". The Fig 6 scope text says "The requested observed-point/model-line … display is implemented". Legend 3E carries five- to six-figure constants. |
| 18 | note | Consistency | Different size grids for strong turning and weak response | addressed in text | The Fig 3 scope now says: "Strong turning and weak response use different amplitudes and size grids". |
| 25 | minor | Role / density | 3E facets over-plotted | moot (panel dropped) | — |
| 28 | note | Consistency | Remainder convention | moot (4I dropped) | — |
| 29 | minor | Presentation / role | Fourier screen promoted into the conclusions | **addressed** | The Fourier analysis is now a mechanistic test (5F, native vs fast cytosol), not an inconclusive screen. |
| 30 | minor | Mechanism → measurement | Lateral-edge assay omitted Goryachev | superseded by 34 | The assay text containing the fix is gone. |
| 31 | minor | Presentation / narrative order | Page-1 gaps | partly | Fig 1 and Fig 6 sentences are present. The Fig 2 headline (5–12-fold) has been dropped from page 1. Figs 3–5 are now described as procedures (33). |
| 32 | note | Clarity | Fig 3 legend placement | addressed | The legend is rewritten; the linearity brackets sit in D. |
| 33 | **major** | Role (novel) / narrative progress | **The outline text does not state the round's novel results or follow the target narrative.** Page 1 and the Fig 3–5 conclusions describe what was done, not what was found. For example: "Figure 5 … test[s] whether long-range transfer contributes to persistence"; "Strong-input patch-field and mode-rate experiments test inactive-pool transport". Missing or buried: (a) the three-route synthesis (drift suppression in Wave-pinning, jump suppression in Goryachev, geometry in Holmes and Spring); (b) that fast cytosol removes Wave-pinning's strong-input size advantage (heading α 0.41 → −0.03 at f = 0.3, 0.72 → 0.21 at f = 0.1; readme only — legend F gives only native-minus-fast differences); (c) that Holmes' strong-input advantage is *not* cytosolic (difference 0.08); (d) that the isolated-front law predicts absolute Wave-pinning ring gains within 6.9–8.0% (readme only). | open | The brief required "outline text following the target narrative". The readmes record the results in full; the outline was not written from them. |
| 34 | **major** | Mechanism → measurement | **Regression of the cellular measurements.** Fig 3 measurements fall from three assays with WGD+ directions and model classes (v12) to two terse lines with neither. The edge-vs-plateau fluctuation test is gone. Figs 4 and 5 share one identical generic sentence ("measure full-perimeter polarity activity, coherence, boundary displacement and cytosolic-pool mobility"). It gives no predicted signature and no discrimination. Cytosolic-pool mobility (e.g. FRAP), the measurement this round's mechanism most directly implies, is named but given no prediction. The target narrative gave one: "measured mobility and perimeter predict whether enlargement increases persistence by more than geometry". Concern 21, addressed in round 6, has reopened. | open | No position recorded. The target narrative (§3 last bullet, §5 observables) supplies the text. |
| 35 | **major** | Narrative progress / presentation | **The Goryachev jump-suppression route is not shown, and 5D cannot be read for it.** The brief frames Figs 3–5 around two mechanisms, one of which is Goryachev jump suppression. No outline text states what Goryachev's persistence mechanism is. The relocation onset shift (f ≈ 0.004 at C20 vs 0.025 at C60, per the target narrative) is not plotted. 5D has 16 facets on a symmetric-log x axis. Their flat-topped densities do not let a reader see a collapsing core vs a non-collapsing tail. The legend does not open with a result, and it says "Collapse is descriptive, not a test passed by the data". | open | Readme "Unexpected": Goryachev C60 strong input has median \|P\| 0.047, so the phase "cannot be identified uniquely as relocation… The new panel does not by itself prove that Goryachev heading loss is exclusively jump-mediated." The brief said "show the panel whatever it shows". The honest calibration is correct. The problem is that the reader is not told what it shows. |
| 36 | minor | Claim strength / presentation | **3E mixes predictions and fits without marking them.** The Wave-pinning curve is an unfitted prediction. The Holmes C\* is a regression on the plotted native gains. The Spring curve combines a sweep-fitted C\* with a measured empirical rise. All are drawn as the same kind of line. The legend's opening claim, "Goryachev loses gain much more steeply than the conduction law predicts", has no Goryachev conduction-law curve to compare against. The pilot `pilot_conduction_penalty.png` (penalty against C/D) is the more direct cross-model test: Holmes collapses cleanly, Wave-pinning approximately, and Spring does not collapse. | open | Readme: "Goryachev's native susceptibility curve is retained without fitting a conduction law". The Spring apparent crossover (142) "is not a transport length". |
| 37 | minor | Claim strength / presentation | **Page 1.** "Published migration *motivates* this hypothesis" is weaker than the vision's "support". The Fig 6 title keeps "support", so the two are inconsistent, and v12 page 1 used "support". Page 1 also describes procedures and revision state rather than claims (see 17 and 33). | open | — |
| 38 | minor | Claim strength / role | **What "heading" means at strong input.** In 4C and 5F the strong-input "heading" is the phase of the first Fourier mode of a patch field with no polarity domain (Wave-pinning C60 f = 0.3 median \|P\| 0.16; Goryachev C60 0.05). The outline should say plainly that the strong-input persistence advantage concerns the orientation of whole-ring activity imbalance, not persistence of an established polarity. That changes how it bears on the vision's polarity → migration link. The 4C conclusion, "decay more slowly for larger rings", also overstates for Spring at f = 0.3, where the curves nearly coincide, and for Goryachev at f = 0.005, where C40 and C60 are pinned at 1. | partly | The Fig 4 scope notes that "Strong-input patch fields can have small first-mode coherence". The readme notes that f = 0.005 sits in different regimes across models. |
| 39 | note | Role / clarity | **5E has no stated implication.** The legend reports that patch density falls about 20% while width rises 28% and speed 23% from C20 to C60 (f = 0.1). This contradicts the target narrative's "local statistics … are the same", and the deliverable reports it honestly. Neither the legend nor the text says what these trends mean for persistence. | open | Readme: the differences are reported, not corrected. |

## Outstanding
Outside this revision's scope or untouched by it.
- **Magnitude gap.** At α ≈ 2, about 1.6-fold per WGD, against a pooled OU fold of 15.8. Neither the outline nor the figures state it. The target narrative's "squared ratio of boundary responses predicts drift ratio" conversion is also absent.
- **Fig 3A–B role.** The target narrative itself flags the organisation of A–B as a known issue: strong-cue turning is not connected to weak-input gain.
- **Bridges not shown as panels.** Pool depletion (Fig 2D) → jump suppression. Wave-pinning conduction (Fig 3) → weak drift and strong low-q rates (Fig 5): this one is now shown in parts (3F, 5A–B, 5F) but not connected in the text.
- **Figs 3–5 have no predicted-readout panel** comparable to Fig 2H–I. The target narrative records this gap.
- **Fig 6** is unchanged. The artwork title ("Migration summaries and conditional OU dynamics") still differs from the outline title.
- **Figs 1–2** carry over the earlier outstanding items: duration magnitude, the scaling-mode mechanism, deBelly's retained state, and the size contrast in G.
- **Holmes strong regime** is unexplained: fast cytosol changes α by only 0.08.

## Narrative targets
From `author/narrative_claims.md` (round 3), with the round-7 target narrative where it supersedes.

| Narrative item | Section | Status | Evidence in deliverable |
|---|---|---|---|
| Shortest pulse per length; duration magnitude | 1 | not addressed (out of scope) | — |
| Mechanistic basis of the scaling modes; deBelly retained state | 2 | not addressed (out of scope) | — |
| Fold change per volume doubling (3B) | 3 | not addressed | — |
| 3C edge-pulse candidate | 3 | addressed (kept) | 3C–D |
| Predicted vs measured drift | 4 | addressed (moved) | 5A–B |
| Imposed → noise bridge | gaps | **advanced** for Wave-pinning (mechanistic) | 3E–F conduction → 5A–B drift → 5F strong-input rates; text does not connect them (33) |
| Strong regime "not identified" | 4 | **partly addressed** | Wave-pinning: cytosolic (5F). Holmes: not cytosolic. Goryachev: open (35). |
| Births and takeovers by distance; exact decomposition | 4 | removed by the brief | — |
| Experimental signatures for Figs 3–4 | gaps | **regressed** | 34 |
| Magnitude vs published fold | gaps | not addressed | — |
| Target narrative: three routes | 5′ | not stated | 33, 35 |
| Target narrative: Holmes and Spring obey the conduction law | 3′ | partly | Holmes fitted (3E) and collapses on C/D (pilot). Spring does not collapse; it has an intrinsic rise (36). |

## Narrative claims
| Claim | Figure | Still supported? | Note (regression vs stale narrative) |
|---|---|---|---|
| Fig 1–2 claims | 1–2 | supported (unchanged) | Stale lettering, as in the prior review |
| Two turning regimes; larger rings slower everywhere; distance matters little beyond the regime | 3A–B | supported | — |
| α ≈ 2 in Holmes and Spring, 2.5–3 in Wave-pinning, declining with input | 4E | supported | Stale lettering (4C → 4E) |
| Tolerance folds 2.2 / 2.4 / 3.1 / 4.1 | — | no longer shown | Removed by the brief. Suggested narrative update. |
| Larger ring holds heading longer | 4C | supported at weak input | Strong-input qualification (38) |
| Drift reduced in Wave-pinning, not in Holmes or Spring | 5A | supported | Spring drift ratio rises to about 1.35 at mid sizes; suggested narrative note |
| Drift plus relocations within 15% in the drift regime | 5C | supported | — |
| Target: strong-input local statistics size-independent | 5E | **not supported** | Density −20%, width +28%, speed +23%. Suggested narrative update. |
| Target: Holmes strong regime "consistent with the same [cytosolic] mechanism" | pilot | **not supported** | Fast cytosol changes Holmes α by about 0.08. Suggested narrative update. |
| Target: fast cytosol drops Wave-pinning α 0.9 → 0.4 (f = 0.1), 0.45 → 0.1 (f = 0.3) | 5F | supported, stronger at f = 0.3 | Mode-1: 0.90 → 0.36 and 0.46 → 0.02. Pilot comparisons are disclosed. |
| Target: isolated front predicts Wave-pinning gain within 6–9% | 3E | supported (6.9–8.0%) | In the readme, not the outline |

## Result roles
| Result | Role(s) | Presented to maximize role? |
|---|---|---|
| 3A–B turning | novel; hypothesis-generating | Figure yes; assay text lost its WGD prediction (34) |
| 3C map, 3D linearity | pedagogical / diagnostic | Acceptable. D could be smaller. |
| 3E conduction law | novel | Partly: prediction and fits are indistinguishable (36) |
| 3F D_v sweep | novel (mechanism by perturbation) | Figure yes. The text does not lead with "fast cytosol removes the size effect" in the conclusions; the legend does say it. |
| 4A–B nulls | pedagogical | Yes, labelled as such in the legend |
| 4C–E | novel | Yes at weak input; strong input needs the reframing in 38 |
| 5A–C drift, α, closure | novel bridge; diagnostic | Yes |
| 5D raw heading changes | intended as a mechanism test | No (35) |
| 5E patch field | descriptive | Implication unstated (39) |
| 5F fast-cytosol mode rates | novel; hypothesis-generating | Under-stated (33). Measurement implication missing (34). |
| Fig 6 | supporting evidence | Unchanged; verb mismatch with page 1 (37) |

## Detail

**33. Suggested fix.** Write the page-1 Fig 3–5 paragraphs and the conclusion bullets as findings, as in the target narrative:
- "Fast cytosol removes Wave-pinning's size-dependent edge gain; an isolated-front conduction law predicts the ring gains within 8% without fitting."
- "At strong input, fast cytosol abolishes Wave-pinning's size advantage (heading α 0.41 → −0.03 at f = 0.3) but changes Holmes by only 0.08."

State the route table, marking the Goryachev route as not established (35). In legend 5F, report the native and fast α values, not only their difference.

**34. Suggested fix.** Restore the v12 Fig 3 assay wording, with its WGD+ direction and the Wave-pinning/Goryachev vs Holmes/Spring discrimination, and the edge-vs-plateau test. Give Figs 4 and 5 separate assay lists. Add the cytosolic-mobility prediction: if the inactive pool's diffusivity and the cell perimeter put C near C\*, WGD+ should gain more persistence than geometry predicts; if C ≪ C\*, persistence follows geometry. Say which mechanisms this assay cannot separate. Goryachev-like jump suppression would need the Fig 2 first-rise assay.

**35. What the reader needs.** A one-line result for each model at the weak input: core collapses / tail present / unresolved. The Goryachev C60 weak condition has zero counted switches and \|P\| about 0.93; say that this is consistent with suppression, not evidence of it. A simpler display, such as tail weight beyond a fixed physical displacement against lag, would make the intended signature legible. Plot the relocation onset against f at C20 and C60 if it is to carry the Goryachev claim.

**36.** Draw the Wave-pinning curve solid and labelled "predicted (isolated front)", and the Holmes and Spring curves as "fitted". Either drop "than the conduction law predicts" for Goryachev or show the bound being exceeded. Consider replacing or supplementing 3E with the C/D collapse from the pilot. It tests the law across models directly, and it shows Spring's failure honestly.

**Verification of earlier fixes.** 29 and 32 are resolved by the restructuring, verified in the v13 artwork and legends. 30's fix was lost when the assay text was cut (now 34).
