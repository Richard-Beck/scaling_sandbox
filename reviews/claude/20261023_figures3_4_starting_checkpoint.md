# Review: figures3_4_starting_checkpoint — 2026-10-23

Deliverables reviewed (commit `214db68` unless noted):
- `manuscript_outline/outline_v10.pdf` (all 20 pages; Figs 3–4 compared page by page with `outline_v9.pdf` at the baseline)
- `manuscript_outline/revision_figure4_scoring_20261023/figure4_with_legend.pdf`
- `manuscript_outline/revision_round3_figure3_single_C_20261022/figure3_with_legend.pdf`
- `polarity1d/.../dev/revision_round4_20261023/figure4/legend.json` and `figure4_main.png`: **working-tree copies; these files are not tracked in `214db68`** (see Assumptions)
- Supporting CSVs used to check claims (working tree, round-4 `figure4/supporting/data/`): `closure_by_regime.csv`, `plotted_E.csv`, `jump_scoring_coverage.csv`, `drift_decomposition.csv`

Baseline: `87cb5f5`, given in the request. It is also the narrative's stated baseline. There are no prior Claude reviews of this target.

Assumptions:
- **Untracked deliverable.** The requester listed `revision_round4_20261023/figure4/legend.json` as a changed deliverable, but nothing under `revision_round4_20261023/` is tracked at `214db68`. The working tree was clean at launch, so the directory is ignored. I reviewed the working-tree file. Its legend text is identical to the committed `figure4_with_legend.pdf` and to outline_v10, so the review does not depend on which copy I used.
- **Figure 3 has not changed.** `figure3_with_legend.pdf` (single-C) is new to git, but its artwork and legend are identical to Figure 3 in `outline_v9.pdf`, which was committed at the baseline. Outline_v10 also reproduces Figure 3 unchanged. I therefore treat Figure 3 as having no content change in this revision.
- **The narrative is current for Figure 4.** `author/narrative_claims.md` names outline_v9 as its baseline, but its Section 4 describes features found only in the v10 / scoring-revision Figure 4: normalized drift in row 1, the dashed mixture fit, 0 of 33 and Goryachev ratios. In v9, row 1 shows absolute D_s and Goryachev is "Not estimable". The user evidently wrote Section 4 while looking at the then-uncommitted scoring revision. I treat the narrative as current for Figure 4 and stale only in its header.
- **"Next revision Figures 3 and 4.md".** This file is new at `214db68` and signed "@Richard", but it sits outside `author/`. Following CLAUDE.md, I did not treat it as an author document. I read it only during the rationale pass, and used it only to read the target name: this commit is the *starting point* for that plan, not its execution. Where the plan already covers a concern, the table says so.
- **Request content.** The request contained only facts (target, commits, paths). It included no rationale.
- **Out of scope.** The diff also deletes many unrelated repository files (`R/`, `julia/`, `docs/`, `MANUSCRIPT_VISION.md`, …) and adds review tooling. These are not manuscript deliverables and were not reviewed.

Revision scope: a checkpoint commit. The only change visible to a reader is Figure 4 (v9 → v10):
- Panel E row 1 now shows drift normalized to the smallest ring at the same input, D_s/D_s(C_min), instead of absolute D_s.
- Rows 2–3 now use counted replacements with jump angles taken from the "scorable" events, which makes 154 conditions estimable (v9: 42). The detector-free mixture fit is overlaid as dashed curves.
- The legend, the narrative points, the equations (adds 4.4) and the outline's page-1 summary of Figure 4 were rewritten to match.

Figures 1, 2, 3 and 5 and the reserved Figure 6 are unchanged.

## Summary
Verdict: **more vision-aligned than baseline: yes (modestly). Revision executed in a vision-aligned way: partly.**

Two changes are genuine gains:
- The normalized-drift row directly shows the size effect against the α = 2 drift reference.
- The closure test no longer reports agreement in "all 42 estimable conditions" (v9). It now reports the failures of the multiple-front regime openly.

There are three main problems:
1. **The headline closure failure is stated more strongly than the scoring supports.** The figure reports under-prediction "0/33" and does "not establish closure", and the narrative concludes that heading is lost through an unidentified process. The supplement's own missing-angle bounds span closure in every condition I checked. In those conditions only 33–54% of events are scored, and the scored events turn the heading by only a few degrees (mean cos Δθ 0.98–0.997).
2. **A new Goryachev result is shown but never stated.** Goryachev's event-based prediction falls short of the observed decay by 3- to 170-fold at every estimable condition, and the shortfall grows with size. The 57/154 pooled count hides this.
3. **The novel result is pushed further out of view.** Revision bookkeeping ("915,757 counted replacements include 321,423 scored angles", "No new simulations were run") replaced v9's statement of what Figure 4 shows, both on page 1 and in Figure 4's narrative points.

The next-revision plan builds its Figure 4E around the multiple-front heading loss. Concern 1 should therefore be settled before that work starts.

## Concerns
| # | Severity | Criterion | Concern | Status | Codex's position |
|---|----------|-----------|---------|--------|------------------|
| 1 | major | Claim strength | Under-prediction in the multiple-front regime ("0/33 event-based … within 15%"; "does not establish closure") is presented as a result. But only 33–54% of events are scored, the scored jumps are near-zero-angle (mean cos 0.98–0.997), and the supplement's missing-angle upper bounds bring k_pred to or above k_obs in every checked condition. The data cannot distinguish "a process the tracked measures do not capture" from "unscored large-angle jumps". | open | The legend's limitations state that the scorable-event assumption "is not a validation of missingness at random", and the bounds are in the supplement. The README: "The decomposition therefore retains its failures in the regime it seeks to explain." The main figure and narrative points still state the failure without the bound. |
| 2 | major | Narrative claims / clarity | Goryachev event-based k_pred/k_obs is 0.006–0.34 at all 37 estimable conditions. It is about 0.14 at C20 and about 0.007 at C60, so drift plus counted replacements explains at most a few percent of Goryachev's heading loss in large rings. This is newly visible in E row 3 (v9: "Not estimable"), but the legend, the narrative points and page 1 do not mention it, and the pooled "57/154" hides it. It undercuts the narrative's account that suppressed drift explains Goryachev's tolerance fold of 4.1. | open | No position recorded. |
| 3 | major | Role (novel) / clarity | The page-1 argument for Figure 4 dropped v9's description of what it shows (size exponents; tolerance folds 2.2–4.1). In their place it reports event-count bookkeeping and revision-log language ("now", "No new simulations were run"). All three Figure 4 narrative points are closure accounting. None says that larger rings hold their heading longer, that α ≈ 2 at weak input, or that the advantage shrinks as input strengthens. The figure title is still the method ("Null models…"), whereas Figures 1, 2, 3 and 5 have claim titles. | open | No position recorded. The README and legend.json repeat the same bookkeeping paragraph. |
| 4 | minor | Clarity | The Figure 4 subtitle in outline_v10 (p. 13) reads "main", a role field leaking into the text. In v9 it was a descriptive sentence. This is a regression. | open | — |
| 5 | minor | Information density | The detector-free mixture fit (dashed lines in E rows 2–3) disagrees with the counted estimates by up to roughly four orders of magnitude at weak input. In the main figure it serves only as a diagnostic, and it doubles the line count in two rows. | open | The plan in "Next revision Figures 3 and 4.md" (user-signed, outside `author/`) moves the λ_eff row and the detector-free fit to the supplement. |
| 6 | minor | Clarity | The running header on every page ("FIGURE 4: CORRECTED EVENT SCORING AND NORMALIZED DRIFT"), the subtitle ("…with corrected drift/jump comparisons") and the dateline are revision labels. They are not part of the manuscript argument. | open | — |
| 7 | minor | Terminology | The same event class is called relocation (A–B), jump (B, E, eq. 4.3), replacement (D, E), switch (D key, narrative), and scored, scorable or counted event (E, eq. 4.4). "Replacement" is also Figure 1–2's term for imposed-cue takeover. The revision added two more variants. | open | — |
| 8 | minor | Clarity | E row 1 uses a log axis from 10⁻⁴ to 1 for Goryachev, but a ¼–8 axis for the other three models. The legend does not say so. | open | — |
| 9 | note | Notation | Eq. 4.3 now uses C, while eq. 4.1 and v9's eq. 4.3 use C̄ (actual mean circumference). | open | — |

## Outstanding
These items were not changed in this revision and are not counted against it.
- **Figure 3 is unchanged.** It still has no plotted fold-change per volume doubling. Panel C repeats B (speed = half-angle/t₅₀). There is no edge-pulse panel and no systematic readout of how coexisting domains resolve. The next-revision plan targets C/D.
- **Figures 3–4 still offer no discriminating experimental signature** comparable to Fig. 2F. Figure 4's "Proposed measurements" text is unchanged from v9 and names no expected direction or mechanism. The next-revision plan's pilot-observables table addresses this.
- **Figure 4D Goryachev** kymographs are almost blank at the shared activity scale, and α₂₀,₆₀ is unavailable. This is unchanged from v9; the plan's per-size scale addresses it.
- **No figure bridges** imposed replacement (Figs 1–3) and noise-driven persistence (Fig. 4). Births vs takeovers are still not plotted.
- **Magnitude gap:** about 1.6-fold per WGD at α ≈ 2, against an OU fold of about 16. Not shown.
- **Verb mismatch on page 1:** published migration data "supply empirical motivation", but the vision and Figure 5's title say "support". Unchanged from v9.
- **Figures 1, 2 and 5** carry no change. Their narrative gaps (Figure 1 duration magnitude; Figure 2 scaling-mode mechanism, deBelly's retained state, the size contrast in 2F) remain open.

## Narrative targets
| Narrative item | Section | Status | Evidence in deliverable |
|---|---|---|---|
| Shortest-pulse duration curves; magnitude of the increase | 1 | not addressed (out of scope) | Figure 1 unchanged |
| Mechanistic basis of the two scaling modes | 2 | not addressed (out of scope) | Figure 2 unchanged |
| deBelly state retained across gaps | 2 | not addressed (out of scope) | — |
| Size differences in 2F | 2 | not addressed (out of scope) | — |
| Fold change per volume doubling | 3 | not addressed | Figure 3 identical to v9 |
| Regimes ↔ ways heading is lost (speculative) | 3 | not addressed | — |
| 3C candidate: edge-pulse response | 3 | not addressed (planned next) | 3C still shows mean angular speed |
| Resolution of coexisting domains after cue withdrawal | 3 | partly, unchanged | Single example in the 3A "cue off" row only |
| Directly measured drift, normalized to smallest ring | 4 | addressed in this revision | E row 1, D_s/D_s(C_min) with unity reference |
| Drift f-scaling as a small-signal response | 4 | not addressed | Not plotted |
| Detector-threshold sweep (candidate supplement) | 4 | not addressed | — |
| Heading loss in the multiple-domain regime "not identified" | 4 | made more prominent, but premise questionable | E row 3 under-prediction; see Concern 1 |
| Births vs takeovers; heading-change decomposition; birth-rate panel | 4 | not addressed (planned next) | Multiple-front fraction row only |
| Experimental signatures for Figures 3–4 | gaps | not addressed | Measurement text unchanged |
| Imposed → noise bridge | gaps | not addressed | — |
| Magnitude vs published fold | gaps | not addressed | — |

## Narrative claims
| Claim | Figure | Still supported? | Note (regression vs stale narrative) |
|---|---|---|---|
| Two turning regimes; larger rings slower at every cue position | 3A–B | supported | No change |
| Distance matters mainly via regime (WP/Holmes/Spring < ~10%) | 3B | supported | No change |
| α ≈ 2 at weak input (Holmes, Spring), 2.5–3 (WP); declines toward 1 | 4C | supported | No change |
| Tolerance folds 2.2 / 2.4 / 3.1 / 4.1 | 4C | supported | Values match the legend |
| Larger ring holds heading longer near the transition | 4D | supported (WP, Holmes, Spring) | Goryachev panel uninformative (pre-existing) |
| Size reduces drift in WP (~2×) and strongly in Goryachev, not in Holmes/Spring | 4E row 1 | supported | Consistent with the figure |
| Drift + counted relocations within 15% in the drift regime | 4E row 3 | supported for WP/Holmes/Spring | Not for Goryachev at any input (Concern 2) |
| 0 of 33 multiple-domain conditions within 15%; heading lost through an uncaptured process | 4E row 3 | weakly supported | Missing-angle bounds span closure (Concern 1). The narrative should carry this qualifier. |
| Goryachev suppressed drift ↔ tolerance fold above drift-only expectation | 4E/4C | weakened | Drift + events explain ≤ ~15% of Goryachev decay at C ≥ 30 (point estimate) |
| Detector-free fit over-predicts by orders of magnitude at weak input | 4E dashed | supported | — |
| Header: "outline_v9, commit 87cb5f5" | — | stale narrative | Section 4 matches v10, not v9 |

## Result roles
| Result | Role(s) | Presented to maximize role? |
|---|---|---|
| Fig 3 regimes and slowing with size | novel; hypothesis-generating | Partly. The photoactivation readout is named, but no expected signature is given (unchanged). |
| Fig 4A–B null exponents and mixtures | pedagogical | Mostly. Brief and illustrative, but takes two full rows; collapsing B (as planned) would help. |
| Fig 4C–D larger rings retain heading; α near 2 | novel | No. Not stated in the narrative points or on page 1 (Concern 3). |
| Fig 4E row 1 normalized drift | novel (mechanistic) | Yes. Clear unity reference tied to α = 2. This is the revision's best change. |
| Fig 4E rows 2–3 closure | novel / diagnostic | No. The headline failure rests on an untested scoring assumption, the Goryachev result is unstated, and the detector-free overlay adds density. |
| Fig 4E row 4 multiple-front fraction | hypothesis-generating | Partly. Descriptive only; births and takeovers are not separated. |

## Detail

**Concern 1: closure failure versus scoring coverage.** Legend E states that unscored events "are represented by the scorable-event mean". It also states that "this assumption is not a validation of missingness at random". The round-4 supporting tables suggest the assumption is unlikely to hold and that the bound is wide enough to reverse the conclusion. Representative primary rows from `drift_decomposition.csv` and `jump_scoring_coverage.csv`:

| Condition | multi-front frac. | scored frac. | mean cos Δθ (scored) | k_obs | k_pred (point) | ratio | λ_eff missing-angle upper | k_pred using upper |
|---|---|---|---|---|---|---|---|---|
| Holmes C60, f = 0.03 | 0.69 | 0.54 | 0.984 | 3.9e-4 | 1.1e-4 | 0.29 | 5.2e-4 | ≈ 6.3e-4 (> k_obs) |
| Holmes C60, f = 0.1 | 0.91 | 0.51 | 0.960 | 2.0e-3 | 8.9e-4 | 0.46 | 2.4e-3 | ≈ 3.2e-3 (> k_obs) |
| Goryachev C60, f = 0.03 | 0.08 | 0.42 | 0.993 | 2.1e-3 | 1.4e-5 | 0.007 | 1.8e-3 | ≈ 1.8e-3 (≈ 0.9 k_obs) |
| Goryachev C60, f = 0.1 | 0.03 | 0.33 | 0.997 | 1.2e-3 | 6.7e-6 | 0.006 | 1.5e-3 | ≈ 1.5e-3 (> k_obs) |

The upper k_pred values are my calculation: drift term (2π/C̄)²D_s plus the bound λ_eff. The scored events turn the heading by only about 5–16°, so they are not relocations in the null-model sense. Unscored events are the obvious candidates for large-angle takeovers. Until a scoring method covers them (for example, an amount-weighted heading at a fixed offset before and after each counted event), the honest wording is that closure in this regime is *undetermined*, not that it fails. Two things depend on this:
- the narrative's "(not identified)" conclusion and its "secondary domains shift balance" hypothesis;
- the next-revision plan's E(iii) decomposition, which is meant to locate an excess that may partly be a scoring artefact.

At minimum, the main figure should draw the bound band in row 3, or the legend should say the bounds span unity.

**Concern 2: Goryachev.** All 37 estimable Goryachev conditions lie between 0.006 and 0.34. They fall monotonically with size at fixed input; at f = 0.03 the ratio is 0.14 / 0.057 / 0.024 / 0.013 / 0.007 for C20–60. Goryachev's valid (dominant-domain) fraction is also only 0.24–0.50 at C60, so drift is measured over a minority of the time. These facts matter in three places:
- They bear on the narrative's reading that strongly suppressed drift accounts for Goryachev exceeding every null.
- They bear on the next-revision plan, which includes Goryachev in the Fig. 3D drift prediction and expects heading persistence to scale "roughly by size squared".
- They may share Concern 1's cause, because Goryachev's scored fraction is also below one half.

Either way, the figure shows this result, so the legend should state it.

**Concern 3: page 1 and narrative points.** v9's page 1 said Figure 4 "compares local size exponents and noise tolerance with drifting-peak and random-relocation nulls" and gave the fold range. v10 replaces this with the README paragraph. A reader of the outline no longer learns from the summary that enlarged rings keep their heading longer. Suggested direction: lead with the claim from the narrative's Section 4 title, keep one sentence on closure (qualified per Concern 1), and move event counts to the legend or methods.

**Concern 4.** `outline_v10.pdf` p. 13, directly under the Figure 4 title: "main".

**Concerns 5–9.** These are self-explanatory from the table. Concern 5 is already in the plan. Concerns 7 and 8 are cheap legend fixes. They are worth doing in the planned slimmed Figure 4, because new terms (birth, takeover, split) are about to be added to the same figure.
