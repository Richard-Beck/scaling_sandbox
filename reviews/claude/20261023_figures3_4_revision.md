# Review: figures3_4_revision — 2026-10-23

Deliverables reviewed (all at commit `49754ab`; the working tree matched the commit for every file viewed):
- `manuscript_outline/outline_v11.pdf`: pp. 1–3 (summary, figure list) and pp. 9–16 (Figures 3–4 artwork, legends and narrative points)
- `manuscript_outline/revision_figures3_4_20261023/figure3_with_legend.pdf`, `figure4_with_legend.pdf` and `pilot_observables.pdf`
- `polarity1d/.../revision_round5_20261023/figure3/figure3_round5.png`; `figure3_round5_supplement.pdf`, pp. 1–14 of 25
- `polarity1d/.../revision_round5_20261023/figure4/figure4_main.png`; `figure4_supplement.pdf`, pp. 1–8
- `manuscript_outline/outline_v11.json`, Figure 3–4 entries (the source of the legend and narrative text)
- CSVs used to check claims: `figure3/supporting/data/edge_prediction_comparison.csv` and `edge_response_size_folds.csv`; `figure4/supporting/data/birth_heading_condition_metrics.csv`

Baseline: `214db68`, given in the request. The prior Claude review of these figures (`20261023_figures3_4_starting_checkpoint.md`) reviewed the same commit, so its concern numbers 1–9 carry forward below. New concerns start at 10.

Assumptions:
- **Revision brief.** `request/revision_brief_20261023.md` (formerly "Next revision Figures 3 and 4.md") is signed "@Richard" but sits outside `author/`. As in the checkpoint review, I treated it as rationale: I read it only in the rationale pass and cite it under "Codex's position". Where a concern follows from something the brief prescribed, the table says so. The user adjudicates between the brief and the vision.
- **Request content.** The headless request contained only facts (target, commits, paths). It included no rationale.
- **Requester's path list vs diff.** All eight listed paths are additions in the diff. The diff also adds several deliverables that were not listed:
  - `figure3_legend.pdf`, `figure4_legend.pdf`, `figures3_4_with_legends_and_pilots.pdf`, `figures_2_4.pdf` and `pilot_observables.png`;
  - `figure3_round5_new_supplement.pdf` and `figure4_round5_new_supplement.pdf`;
  - the round-4 Figure 4 exports that were ignored at the checkpoint and are now tracked.

  I did not review these separately. I assume they are subsets or assemblies of the reviewed files.
- **Out of scope.** Figures 1, 2 and 5 and the reserved Figure 6 were not changed, per the brief. I checked only the outline's summary text about them.

Revision scope: Figure 3 replaces the old panel C (mean angular speed, a restatement of B) with two new panels:
- C: deterministic weak edge-pulse responses at the Figure 4 sizes;
- D: a no-fit prediction of Figure 4's tracked drift from those responses.

Figure 4 is slimmed:
- the three null-mixture α panels, the λ_eff row, the detector-free fit and the multiple-front fraction move to the supplement;
- kymographs are re-centred, and Goryachev gets per-size activity scales;
- the closure comparison is restricted to the drift regime;
- a new panel E adds certified secondary-domain birth rates (i), takeover probability by size and distance (ii), and an exact signed heading-loss allocation (iii).

A new prospective pilot-observables table is added, and the outline is rebuilt as v11 around these changes.

## Summary
Verdict: **more vision-aligned than baseline: yes. Revision executed in a vision-aligned way: partly.**

The revision makes real progress on the narrative's main gaps:
- Figures 3–4 now carry experimental signatures (the pilot table, and the Fig. 3C photoactivation counterpart).
- Fig. 4E(i–ii) measures births and takeovers directly.
- Fig. 4 now leads with its size claim.
- Most notably, Fig. 3C–D links an imposed-cue response to noise-driven drift. This is the first figure-level bridge between Figs 1–3 and Fig. 4.
- The checkpoint review's three major concerns are addressed or mostly addressed.

There are four new major concerns. Each one buries or blurs this round's main new results:
1. **Fig. 3D, wrong prediction in the main panel.** It shows a prediction the figure itself refutes (independent edges; the measured edge correlation is 0.98 or higher). The prediction that keeps the shared translation matches measured drift about 1:1, including its size dependence, but appears only in the supplement.
2. **Fig. 3D does not test size.** Its five-decade 1:1 span comes mostly from the f² scaling built into linear response, so the panel does not display the size test that matters.
3. **Fig. 3's title and narrative drop the size claim.** Neither states the size-dependent edge susceptibility (the mechanism → measurement link). The baseline title said "enlargement slows both"; the new title describes the method.
4. **Fig. 4E(iii), normalization artefact.** Its dominant colours are an artefact: a unit-modulus point-region moment set against |P| < 1 gives "motion" fractions of about 2 (Wave-pinning), 2.5 (Holmes) and 6–8 (Spring), with matching negative residuals. This happens exactly in the drift regime, where Fig. 4D shows that drift alone closes the decay.

## Concerns
| # | Severity | Criterion | Concern | Status | Codex's position |
|---|----------|-----------|---------|--------|------------------|
| 1 | major | Claim strength | Checkpoint: the multiple-front closure failure ("0/33") was stated more strongly than the scoring supports. | **addressed** | The main figure now shows closure only where f ≤ 0.01 and multiple-region frames are ≤ 5% (45/47 within 15%). Strong-input estimates and bounds moved to supplement S4.2. Codex's audit: the missing-angle bounds span unity in 18/33 conditions and allow agreement within 15% in 24/33; 9 remain outside even at the bounds. |
| 2 | major | Narrative claims | Checkpoint: Goryachev's event-based shortfall was unstated. | **mostly addressed** | The Fig. 3 narrative says "Goryachev does not give a general drift-based explanation of heading persistence". The Fig. 4 legend gives "Goryachev 0/1". Fig. 3D gives Goryachev's spread (0.0036–4.0). The 37-condition shortfall itself is visible only in S4.2. |
| 3 | major | Role (novel) | Checkpoint: Fig. 4's size result was buried under bookkeeping. | **addressed for Fig. 4** | The title is now "Enlargement extends heading memory…"; legend B and the first narrative point lead with the claim. The same problem has reappeared in Fig. 3 (Concern 12). |
| 4 | minor | Clarity | "main" leaked into the subtitle. | addressed | The role line is now descriptive. |
| 5 | minor | Density | Detector-free overlay in the main figure. | addressed | Moved to S4.2. |
| 6 | minor | Clarity | Revision labels in the header and subtitle. | addressed for the header; partly for legends | The running header is now scientific. Legends still contain revision-log language (Concern 17). |
| 7 | minor | Terminology | Event-class terms varied. | mostly addressed | Legend A now separates null relocations from measured replacements and birth-to-takeover events. Panel E adds "membership", "certified birth" and "takeover", each defined. |
| 8 | minor | Clarity | The Goryachev axis differed without explanation. | addressed | Declared in legend D. |
| 9 | note | Notation | C vs C̄. | addressed | Eq. 4.4 and legend D use C̄. |
| 10 | **major** | Claim strength / role (novel) | Fig. 3D's main panel plots the independent-edge prediction D_s = (D_L+D_R)/4. The legend itself reports edge-increment correlations of 0.980 / 0.983 / 0.997 (Wave-pinning / Holmes / Spring), so the plotted prediction is known to be wrong by a factor of 2, and the panel reads as a uniform two-fold under-prediction (medians 0.51–0.54). The prediction that keeps the shared translation (supplement p. 7; `full_covariance_over_observation`) is about 1.0, e.g. 0.82–1.13 across Wave-pinning conditions, with no fitting. The positive, novel result is therefore in the supplement: domain drift is rigid translation, and its rate is predicted from the deterministic edge susceptibility. The narrative's own speculative premise (independent edges) is refuted, and that refutation should be stated as a result. | open | The brief prescribed the independence formula and said "Disagreement is a result. Don't adjust anything to close it." METHODS_RESULTS: "we retain the failing independence prediction in the main panel"; the covariance version is "a separate, untuned sensitivity". |
| 11 | **major** | Presentation / claim strength | Fig. 3D spans about five decades in D_s. That spread comes almost entirely from input strength f, and linear response builds in D ∝ f². Close 1:1 agreement over that span is therefore expected and does not test the size question. The size dependence is a small within-cluster spread that cannot be read. The data do support a size test: in Wave-pinning at f = 0.001, measured D_s(C60)/D_s(C20) is 0.47 and the covariance prediction gives 0.49, consistent with terminal edge gain 0.70² ≈ 0.49. No panel shows this. | open | No position recorded. |
| 12 | **major** | Role (novel, hypothesis) / clarity | The Fig. 3 title changed from a size claim ("…enlargement slows both") to a method description ("Local edge pulses expose coupled boundary motion and test an untuned prediction…"). The new panel C contains a novel size result: terminal edge displacement per pulse dose falls with size in Wave-pinning (×0.70, C20→60) and strongly in Goryachev (×0.088), but barely in Holmes (×0.91) or Spring (×0.99). This matches which models show size-suppressed drift in Fig. 4D, so it is the mechanism → measurement link the vision asks for. It appears only as a numeric list at the end of legend C. The narrative points and the page-1 summary report only the 0.51–0.54 independence ratios. | open | No position recorded. The role line reads "Cue placement separates sliding from takeover; local boundary response tests a no-retuning prediction…". |
| 13 | **major** | Claim strength / presentation | Fig. 4E(iii)'s dominant pattern is a normalization artefact. Q = ΣA e^{iθ}/ΣA has unit modulus for one region, while the measured P = M₁/M₀ has \|P\| < 1. In single-domain conditions the "motion" fraction is therefore about 1/\|P\| and the residual about 1 − 1/\|P\|. Values from the source CSV: Wave-pinning C20 at f = 0.001 and 0.01 gives M = 1.98 / 2.00 and R = −0.98 / −1.00; Spring gives M = 6.2–7.8 and R = −5.2 to −6.8. This happens in the very regime where Fig. 4D row 2 shows drift alone closing the decay within 15%. \|P\| also changes with size (Spring M 7.8 → 6.2, C1 → C3), so comparisons between the size rows are confounded. The potentially informative trend (membership B rising with f and size; amount A small) is hard to read against this. The legend says only "convention-dependent… cancellations are retained". | open | METHODS l. 42: "In a broad single domain, \|Q\|=1 while \|P\|<1; substantial cancelling motion/residual contributions are consequently expected under this mandated normalization." RESULTS: values "are not physical percentages". The brief specified P = ΣAᵢe^{iθᵢ}/ΣAᵢ "plus a background residual". |
| 14 | minor | Claim strength | Fig. 4E(ii) measures takeover *within 20τ* (1000 units). In Holmes, the continuous-cue half-turn time at C47.6 is about 3×10³ (Fig. 3B), longer than the window. The decline of takeover with size may therefore partly reflect slower takeover, not rarer takeover. The pilot table turns this into "Fraction of secondary protrusions that take over: Lower". Certification coverage is also low for Spring at moderate f (about 7–14%; S4.3), and Holmes sits around 20–30% over much of the f range. | open | The legend says the decline "applies to this finite-horizon, certified subset; unequal certification prevents a population-level birth/takeover claim". The pilot table labels the row "Conditional model result". |
| 15 | minor | Presentation | The Fig. 3C insets use free, non-zero y-axes. Holmes (455 → 438, −4%) and Spring (29.6–30.1, ±1%) look as size-dependent as Wave-pinning (−24%). Spring even looks non-monotone. A shared normalized axis (A/A(C_min)) would show the real contrast, which is the point of Concern 12. | open | — |
| 16 | minor | Mechanism → measurement | The pilot table is organized by model, not by the Fig. 2 mechanisms (depletion-limited initiation; saturation and expansion; mechanical memory), and it does not say which mechanism each observable supports or fails to separate. In particular, the edge-response row separates Wave-pinning from Holmes, although both are expansion-constrained in Fig. 2. That this observable sorts the models differently from Fig. 2's routes should be stated. Two further points: given the near-perfect edge co-motion, row 1 would be better phrased as common-mode boundary translation; and the baseline's straddling-vs-separated photoactivation half-turn assay (baseline Fig. 3 legend; narrative §4) no longer appears in Fig. 3 or in the table. | open | The brief supplied the table's rows and allowed directions to be updated. Codex updated two directions from data (births "model dependent"; secondary-regime edge fluctuation "no common direction"). |
| 17 | minor | Density / clarity | The legends are methods-length: Fig. 3 runs about 1.5 pages plus eqs 3.1–3.4, and Fig. 4E is one paragraph of about 400 words plus eqs 4.1–4.3. They also retain revision-log wording, e.g. "Curves, gaps, brackets and three Figure 3 sizes are unchanged", "the existing periodic watershed detector", "the original 10⁻⁴ scaled-RHS tolerance", "the existing selected input… raw evidence are unchanged", "prescribed", "declared". Typo: "boundary).." in legend B. Eqs 3.3–3.4 and the relaxation-row protocol belong in methods. | open | — |
| 18 | note | Consistency | Within Fig. 3, panel B uses the Figure 3 sizes (C 30/37.8/47.6; Spring 2.5–3.97), while C–D use the Figure 4 sizes (20–60; Spring 1–3). A reader cannot compare B with C at a common size. | open | METHODS_RESULTS: "C uses Figure 4's size grid, not B's three sizes." The brief allowed both. |

## Outstanding
These items fall outside this revision's scope or were left unchanged by it; they are not counted against it.
- **Magnitude gap** (out of scope per the brief). At α ≈ 2, one WGD gives a heading-persistence fold of about 1.6, against a published OU fold of about 16. Not shown.
- **Fig. 3B fold change per volume doubling.** Still not plotted (B is unchanged).
- **Detector-threshold sweep.** Not done.
- **Resolution of coexisting domains after cue withdrawal.** Still a single exemplar in 3A; the round-2 resolution analysis survives only in the supplement.
- **Verb mismatch on page 1.** "Existing migration observations supply empirical motivation" against the vision's "support". Unchanged.
- **Figs 1, 2 and 5.** Their narrative gaps (duration magnitude; the scaling-mode mechanism; deBelly's retained state; the size contrast in 2F) remain open, by design.
- **Supplement curation.** The Figure 3 supplement is 25 pages and still carries "superseded" round-2 pages. A curation pass is needed before submission.

## Narrative targets
| Narrative item | Section | Status | Evidence in deliverable |
|---|---|---|---|
| 3C candidate: edge displacement after a weak, brief pulse at the incumbent's edge, at each size | 3 | **addressed** | Fig. 3C at five sizes in all four models (see Concerns 12 and 15 on presentation) |
| Predicted vs measured drift (candidate panel) | 4 | **addressed, presentation compromised** | Fig. 3D; Concerns 10–11 |
| Drift as a small-signal response (f-scaling) | 4 | addressed (legend only) | Legend D: exponents 1.87–1.99. MSD exponent not reported. |
| Fold change per volume doubling | 3 | not addressed | B unchanged |
| Two regimes ↔ two ways heading is lost (speculative) | 3 | partly | E(ii) distance classes (d ≤ 0, near, far) overlap, so takeover is distance-insensitive. Not linked to Fig. 3's regimes in the text. |
| Resolution of coexisting domains after withdrawal | 3 | partly, unchanged | 3A "cue off" row |
| Births per perimeter; takeover by distance and size | 4 | **addressed** (with identification limits) | E(i) certified lower bounds; E(ii) certified subset, finite horizon (Concern 14) |
| Exact decomposition of heading change | 4 | addressed in form; compromised | E(iii); Concern 13 |
| Heading loss in the multiple-domain regime "not identified" | 4 | partly | The 0/33 claim is withdrawn from the main figure. E(iii) hints that membership changes dominate at high f, subject to Concern 13. |
| Detector-threshold sweep | 4 | not addressed | — |
| Experimental signatures for Figs 3–4 | gaps | **largely addressed** | Pilot table; Fig. 3 "cellular counterpart" (Concern 16) |
| Imposed → noise bridge | gaps | **partly addressed** | 3C→3D links an imposed deterministic response to noise drift. E(ii) shows takeover declining with size, and the legend cites the Fig. 2B–C asymmetry. Neither is stated as the bridge. |
| Photoactivation straddling vs separated (schematic predictions) | 4 | partly | Straddling only; the separated assay was dropped (Concern 16) |
| Magnitude vs published fold | gaps | not addressed (out of scope) | — |

## Narrative claims
| Claim | Figure | Still supported? | Note (regression vs stale narrative) |
|---|---|---|---|
| Two turning regimes; larger rings slower at every cue position | 3A–B | supported | The data are unchanged, but the claim is no longer the figure title (Concern 12; regression in presentation). |
| Beyond the regime boundary, distance matters little (WP/Holmes/Spring) | 3B | supported | Consistent with E(ii)'s overlapping distance classes. |
| α ≈ 2 at weak input (Holmes, Spring), 2.5–3 (WP), declining toward 1 | 4B (was 4C) | supported | Stale narrative: panels renumbered. Goryachev now has 11 resolved local neighbourhoods (open = forced), so "too abrupt for local estimates" needs updating. |
| Tolerance folds 2.2 / 2.4 / 3.1 / 4.1 | 4B | supported | Legend: 2.155 / 2.373 / 3.088 / 4.122. |
| Larger ring holds heading longer near the transition | 4C (was 4D) | supported | Goryachev is now readable: the C60 exemplar does not decorrelate at all. |
| Size reduces drift in WP (~2×) and strongly in Goryachev, not in Holmes/Spring | 4D (was 4E row 1) | **strengthened** | Independently corroborated by the Fig. 3C edge gains (Concern 12). |
| Holmes/Spring "surprising" given Fig. 3 sliding | 4D vs 3B | resolved | The weak-edge susceptibility, not strong-cue sliding, predicts drift. The narrative's reconciliation is now shown (supplement p. 7 rather than main). |
| Drift f-scaling f^1.94–1.99 in WP, Holmes, Spring | 4D legend | supported | Stale numbers: Spring is 1.87–1.93. |
| "If edges respond independently… edge pulse predicts drift" | 3D | **premise refuted; prediction succeeds with shared translation** | Suggested narrative update: drift is common-mode translation set by edge susceptibility. |
| Drift + relocations within 15% in the drift regime | 4D row 2 | supported | 45/47 |
| 0/33 multiple-domain conditions within 15%; process "not identified" | — | **withdrawn from main** | Codex's audit: bounds span unity in 18/33. The narrative should drop or qualify this. |
| Larger rings spawn secondary domains more readily than they lose incumbents | 4E(i–ii) | partly supported | Per-perimeter certified birth rates are roughly size-independent in WP/Holmes, rise with size in Spring and fall in Goryachev. Takeover within 20τ falls with size in all four. Finite-horizon caveat (Concern 14). |
| Heading lost via "balance between coexisting domains" (speculative) | 4E(iii) | weakened / untested | The amount (A) fractions are small. Membership (B) rises with f and size. Interpretation is limited by Concern 13. |
| Goryachev suppressed drift ↔ tolerance fold above drift-only expectation | 4D / 3D | weakened | The revision itself says Goryachev "does not give a general drift-based explanation". |
| Detector-free fit over-predicts at weak input | S4.2 | supported (now supplement) | — |

## Result roles
| Result | Role(s) | Presented to maximize role? |
|---|---|---|
| Fig. 3A–B regimes and slowing with size | novel; hypothesis-generating | Less than at baseline. The title no longer carries it, and the half-turn photoactivation assay was dropped. |
| Fig. 3C size-dependent edge susceptibility | **novel; hypothesis-generating** | No. It is the clearest mechanism → measurement result in the set, but it is stated only as a number list, and the insets mislead (Concerns 12, 15). The pilot table does carry the direction. |
| Fig. 3D drift prediction from deterministic response | novel (bridge) | No. The main panel shows the refuted variant and does not display the size test (Concerns 10–11). |
| Fig. 4A nulls | pedagogical | Yes. Brief, with α-mixture diagnostics moved to the supplement. |
| Fig. 4B–C heading retention, α, folds | novel | Yes. Now leads the title, legend and narrative. |
| Fig. 4D normalized drift + drift-regime closure | novel (mechanistic) | Yes |
| Fig. 4E(i–ii) births and takeovers | novel; hypothesis-generating | Partly. Honest identification limits, but no takeaway sentence linking it to Fig. 2B–C / Fig. 3 (the legend only says the asymmetry "motivates measuring"). |
| Fig. 4E(iii) heading-loss allocation | intended novel | No. The artefact dominates the display (Concern 13). Better in the supplement until renormalized. |
| Pilot observables table | hypothesis-generating | Mostly. Pre-specified directions with conditions. Not keyed to mechanisms (Concern 16). |

## Detail

**Concerns 10–11: Figure 3D.** From `edge_prediction_comparison.csv`:
- Every Wave-pinning row has `prediction_over_observation` ≈ 0.41–0.57 and `full_covariance_over_observation` ≈ 0.82–1.13. Holmes is similar (about 0.48–0.71 and 0.96–1.42).
- The supplement (p. 4) shows both edges settling to the same persistent displacement after a pulse at either edge. The independence page (p. 5) shows increment correlations near 1 at quarter-window lag.

The figure therefore establishes three things:
- the incumbent's two boundaries move as one (neutral translation);
- a deterministic pulse at one edge fixes the gain of that translation;
- that gain predicts noise-driven drift with no free parameter, including its size dependence. In Wave-pinning at f = 0.001, measured D_s falls from 6.9×10⁻⁶ at C20 to 3.3×10⁻⁶ at C60 (×0.47). The covariance prediction falls from 7.1×10⁻⁶ to 3.5×10⁻⁶ (×0.49).

This is the novel bridge the narrative wanted. Plotting only the refuted variant turns a success into an apparent factor-2 failure.

Suggested fix: show the shared-translation prediction as the primary marker, with the independent-edge prediction as open markers labelled "tested assumption, rejected (edge correlation ≥ 0.98)". Add a size-normalized comparison: predicted vs measured D_s(C)/D_s(C_min) at each f, which could be overlaid on Fig. 4D row 1. This keeps the brief's "disagreement is a result" while showing which disagreement was learned from. The choice is the user's, because the brief prescribed the independence formula.

**Concern 12: Figure 3 role.** Suggested direction: a title that carries both results. For example: "Enlargement slows reorientation by either route, and a weaker boundary response to local input in some models predicts their reduced drift." The narrative point should give the C20→C60 terminal-gain folds, say which models they match in Fig. 4D, and name the cellular assay. Measure displacement per photoactivation dose at the leading-edge boundary in WGD+ vs WGD− cells: it is expected to be smaller in WGD+ cells under Wave-pinning-like or Goryachev-like dynamics, and similar under Holmes-like or Spring-like dynamics.

**Concern 13: Figure 4E(iii).** From `birth_heading_condition_metrics.csv`:

| Condition | motion | amount | membership | residual |
|---|---|---|---|---|
| WP C20, f = 0.001 | 1.977 | 0.000 | 0.000 | −0.977 |
| WP C20, f = 0.01 | 1.998 | 0.000 | 0.000 | −0.998 |
| WP C60, f = 0.1 | 0.284 | 0.372 | 1.169 | −0.826 |
| Spring C1, f = 0.001 | 7.780 | 0.000 | 0.000 | −6.780 |
| Spring C3, f = 0.001 | 6.154 | 0.000 | 0.000 | −5.154 |

In single-domain conditions, motion + residual = 1 exactly, with motion ≈ 1/|P|. A reader of the main heatmap sees dark red M rows and blue R rows everywhere. The natural reading is "motion over-explains heading loss and something unexplained cancels it", which is false: this is the regime where drift alone closes the decay.

Two fixes would make the panel usable:
- define the point-region moment with each region's own first-moment modulus, or rescale Q by |P| so the drift regime reads M ≈ 1, R ≈ 0;
- or report the components as heading-angle contributions rather than moment-projection fractions.

Until then, E(iii) belongs in the supplement (S4.5 already holds the clearer line plots). The main figure should keep a one-line statement of what survives, which is that membership changes, not amount shifts, grow with f and size in Wave-pinning and Holmes, and that this remains descriptive.

**Concern 14.** Takeover could instead be reported as a cumulative incidence curve against time since birth, by size, or the window could be scaled to each size's half-turn time. Either would separate "rarer" from "slower". The narrative's bridge ("incumbents become harder to remove") holds under both readings, but the pilot row "fraction that take over" assumes the first.

**Prior concerns 1–9.** Verified against the deliverables, not against Codex's response table:
- Main-figure closure is now restricted to the drift regime, with the bounds in S4.2 (Concern 1).
- The Goryachev statements are present in the Fig. 3 narrative and the Fig. 4 legend (Concern 2).
- The Fig. 4 title and legend B lead with the claim (Concern 3).
- The header and subtitle are scientific (Concerns 4 and 6).
- The detector-free fit is in S4.2 (Concern 5).
- Legend A separates the event classes (Concern 7).
- Axis ranges are declared in legend D (Concern 8).
- C̄ is used in eq. 4.4 (Concern 9).
