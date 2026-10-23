# Response to Claude: round6 manuscript revision

Claude reviewed commit `209a0c3867444077969cdd17c3789d89a66c3a53` against
the initial checkpoint `b34752047d53bc9d2a9e6cf6d98b99fe2a42ee8b`.
The [review](../../reviews/claude/20261023_manuscript_round6.md) is preserved
verbatim at its original path, with SHA256
`6dd620c7cd406625a473c248f80da028c6de15f2143c916fadf7e8325c17d041`.
[The integrity record](review_preservation.json) is checked by the manuscript
verification script. The frozen reviewed outline and Figure 4 snapshot remain
available; the current [outline_v12](../outline_v12.pdf) contains the responses.
The review request had three abbreviated supplement filenames; Claude located
the actual `figureN_round6_supplement.pdf` files. Current figure legends and the
reader use those complete paths, which are included in the link checks.

The principal changes restore the matched-input Wave-pinning exponent test,
direct model heading correlations, concrete cellular assays, and a scientific
opening argument. Figure 4 adds a separate model-correlation panel D; its
remaining main panels become E (alpha), F (tolerance), G (drift), H (closure),
I (remainder), and J (takeover). The Fourier tests are supplementary S4.7–S4.8. This small
extension of the brief's layout addresses Concern 20 while keeping one quantity
per letter and preserving the specified order of the other plots.

| Concern | Response and disposition |
| --- | --- |
| 10. Independent-edge prediction in the main comparison | Addressed in the reviewed candidate and retained. Figure 3G uses shared translation only. The independent-edge candidate remains supplementary with the boundary co-motion evidence against it. |
| 11. Absolute drift span obscured the size test | Addressed in the reviewed candidate and retained. Figure 3G compares size-normalized predictions and observations; the absolute comparison remains supplementary. |
| 12. Direction and turning-route result lost from the title/narrative | Accepted. Figure 3 restores enlargement slowing turning and selectively suppressing boundary susceptibility, states the model-specific directions, and connects edge-localized fluctuations to the drift prediction. |
| 13. Point-region normalization artefact | The corrected native-amplitude main comparison is retained. The original convention and exact-mask alternatives remain explicitly separated in the supplementary sensitivity analysis. |
| 14. Finite follow-up conflates slower with rarer takeover | Qualified throughout. The main panel reports incidence within 20 correlation times, not eventual success. The existing pooled-input common-cohort 20–80τ sensitivity is retained; its direction is separate from the single-input distance-resolved comparison, which is not monotonic in every bin. Certification selection, censoring and lifetime fate remain unresolved. No new stochastic follow-up was introduced. |
| 15. Free-axis insets | Addressed in the reviewed candidate and retained. Figure 3 uses separate map and time-course panels without the former displacement insets. |
| 16. Missing assay for the two turning routes | Accepted. Figure 3 restores the finite-cue photoactivation-offset assay: vary placement from straddling the incumbent boundary to separated placement and measure half-turn time. It is distinguished from the weak local susceptibility assay. |
| 17. Revision-log prose and excessive method density | Accepted. Figures 3–4 captions remove preceding-edition and campaign-log wording. The opening argument and running header state scientific claims. Detailed gates, amplitude brackets, source history and campaign status are in methods/provenance and dated status records. Requested spatial fractions remain in the Figure 3 legend at readable precision, with full values in CSV. |
| 18. Different size grids in 3B and 3C–H | Qualified and retained visibly. The requested B-matched deterministic sizes and three-size fold comparison are supplementary. The main endpoint maps retain the brief's C20/C60 comparison; no claim of identical grids is made. |
| 19. Wave-pinning missing from the exponent bridge | Accepted. Figure 3H now uses the existing matched-input native deterministic predictions and qualified observations for f ≤ 0.01, requiring multiple-domain fraction ≤ 0.05 across each fixed fit neighborhood and the native eligibility/linearity gates. There are 53/130 qualified comparisons, including 15 Wave-pinning points at f = 0.003, 0.005 and 0.01. Wave-pinning predicts alpha 2.484–2.878 versus observed 2.362–3.435. This is a comparison, not exact agreement. No extrapolation from the f = 0.001 map and no new noise histories are used. Unavailable Goryachev observations remain unavailable. |
| 20. No direct model heading-memory display in main Figure 4 | Accepted. New Figure 4D restores the model heading correlations at the declared transition inputs, immediately after the model kymographs. The nulls no longer have the only direct correlation display. The supplemental curves are described by their scientific content rather than as a superseded result. |
| 21. Cellular measurements and WGD signatures regressed | Accepted. Figures 3–4 again propose paired boundary tracking, weak photoactivation at multiple positions, straddling-versus-separated turning, births per actual perimeter-time, fixed-horizon takeover, and activity-heading correlations. Expected size directions are conditional on model/regime, not numerical predictions for one WGD or unique mechanism identification. Model activity heading, observed migration and conditional OU memory remain distinct. |
| 22. Opening argument read as a status log | Accepted. Page 1 now gives short scientific claims: both turning routes slow with enlargement; spatial susceptibility predicts drift without fitting; enlargement extends heading memory and tolerance; closure becomes incomplete beyond the single-domain regime. Linear-range brackets and campaign status are removed from that argument. Migration observations are described as supporting the hypothesis. |
| 23. Stale pilot references and conditioning | Accepted. A new round6 pilot table updates its references to 3F–G, 4G, the certified-birth supplement, 4J, 4D–F, and 3C/F. The original pilot directions are preserved. The drift gate still defines the conditional drift interpretation, while closure is explicitly displayed at all inputs. The main takeover conditioning is updated to one input per model; the pooled and 80τ common-cohort sensitivities are identified separately. The earlier pilot files remain intact. |
| 24. Missing Figure 2 labels, model mapping and WGD key | Accepted as presentation-only changes. Figure 2G restores kept/lost/lost-then-commits labels, H identifies the models associated with each schematic class, and I has its own WGD−/WGD+ line key. The width-series key is moved nearer its plots. All numerical source tables and schematic coordinates are unchanged. |
| 25. Edge-contribution density and prose in 3F | Partly accepted within the explicit brief. The four model facets and per-model/size fractions remain as requested, but the caption now states the common result plainly: 99.6–100% of predicted drift comes from within two input correlation lengths of the edges. Fractions are rounded for reading; signed bins and exact values remain in CSV. The prose box in F is removed, and cellular boundary-fluctuation measurements are stated in the text. |
| 26. Figure 1 setup prominence and missing duration magnitude | Presentation concern addressed: the preparation row is smaller and the correlation evidence larger. Establishment data, sampled durations and scientific analyses are unchanged. The missing main-figure duration magnitude/shortest-pulse curves remain an outstanding content request outside the brief's style-only scope. |
| 27. Cropping, overlaps, faint labels and undeclared axes | Accepted. Figure 3 moves its A key and threshold labels, increases B-label contrast and moves E labels off the curves. Figure 4 moves missing-rate annotation away from crossed values, declares Goryachev's distinct drift ordinate and removes the title's trailing period. Figure 2's added labels/keys are checked for overlap. Native vector PDFs are preserved during outline integration. |
| 28. Apparently contradictory exact-mask counts | Accepted. The readme presents both numbers together with their controlled conventions: exact own-mask substitution reduces absolute remainder in 7/325 cases under the original amplitude-scaled convention, and 321/325 with a native-total denominator. Full covering masks reconstruct represented activity algebraically under the latter denominator; that reduction does not independently validate a causal domain mechanism. |

The Fourier result remains supplementary under the declared descriptive
promotion rule.
The eight-point slope screen gives partial agreement in Wave-pinning, Holmes
and Spring, while Goryachev fails. Wave-pinning and Spring fail the narrower
three-point slope screen, so that agreement does not establish a robust
explanation of the strong-input regime. The [stencil sensitivity][stencil]
reports all three-, five- and eight-point results. The
[complex/phase comparison][phase] also retains its sparse-pair flags: the large
Spring point-rate ratio uses only four eligible terminal-lag pairs out of 24,
and is not evidence for a physical 330-fold amplitude effect. Phase-rate point
estimates have no bootstrap precision claim. Reaching the 96-seed cap completes
the authorized allocation rule even when a rate remains unresolved or its
log-rate SE remains above 0.1; those outcomes remain explicit.

The authorized mode-only campaign contains 200 conditions and 9,984 completed
histories: 128 conditions at 24 seeds and 72 at the 96-seed cap. No simulation
failures were recorded. At the cap, 37 conditions retain log-rate SE above 0.1
and 18 have unresolved precision. These disjoint outcomes are exported with the
final Figure 4 receipt; completed allocation is distinct from meeting the
precision target for every condition.

The [runtime provenance record][runtime] documents one limitation: the initial
full wrapper's launch-time hash was not recorded. Its first committed version
is recovered history, not contemporaneous proof of the complete initial file.
The native kernel has a prelaunch hash; the worker numerical body matches the
recovered committed version, and the [native-default comparison][kernel]
checks equivalence with the original solver for the actual noise/σ defaults.
The final precision additions use a runtime frozen and hashed before execution.
All original source/basal hashes are preserved.

The review's outstanding questions remain explicit rather than being claimed as
solved: the quantitative gap to published conditional OU folds, Figure 1's
duration magnitude, the mechanistic basis of the Figure 2 classes and retained
mechanical state, detector-threshold sensitivity, and final submission-level
supplement curation. The complete earlier supplementary evidence is preserved
as required by the brief. Native-model size contrasts are not calibrated
cellular folds.

Rebuild instructions, compact plotted values and numerical qualifications are
in the four figure folders. [Integration verification](verification.json)
checks source preservation, native PDF embedding, local links and the review
hash; [storage measurement](../../reports/round6_storage_20261023.json) records
allocated storage without deleting data.

## Follow-up review of the completed revision

Claude reviewed the completed, pushed commit
`028d9f5f7e3ae0c56d1d07a99bc5231f1913a472` against `209a0c3`.
The [follow-up review](../../reviews/claude/20261023_manuscript_round6_2.md)
is preserved verbatim, SHA256
`69e8c13963c31d912c1645c3fb85efa18c35d09db7aae96a43a9c05db7b4f5e2`.
It finds all prior major concerns resolved, with minor presentation and
interpretation edits remaining. Those edits use the existing numerical results.

| Concern | Response and disposition |
| --- | --- |
| 17. Remaining editorial wording and legend density | Accepted. Figure 1 removes the endpoint revision history and states native units directly. Figure 4J retains the scientific conditioning and finite-horizon interpretation, with bootstrap implementation details in methods. The introduction and supplement index remove revision-log wording. Figure 5 changes only an editorial source-link label; its scientific content and assets remain unchanged. |
| 19. Level versus within-model trend | The Figure 3H caption now states that the narrow Wave-pinning prediction range primarily tests the excess of alpha above 2, rather than the within-model trend. |
| 28. Convention dependence of the remainder | Accepted. Figure 4I now states that the magnitude and sign of the remainder depend on the decomposition convention and explains the displayed point-region/native-amplitude convention. The exact-mask identities and controlled 7/325 versus 321/325 comparison remain supplementary. |
| 29. Supplement letters and repeated Fourier conclusion | Accepted. Supplementary rate and slope plots are S4.7–S4.8, with subsequent Fourier diagnostics numbered consistently. Main text no longer calls them K–L. The inconclusive screen is removed from the main legend and conclusion bullets; the limitations retain the required statement that the strong-input regime remains unexplained. The spectrum assay tests a detector-free description and does not by itself distinguish size mechanisms. |
| 30. Drift assay cannot distinguish Wave-pinning from Goryachev | Accepted. The assay now names both model classes and points to the first-rise assay in Figure 2, where initiation slows with size only in Goryachev. |
| 31. Missing opening results and isolated source-link page | Accepted. The introduction adds explicit Figure 1 and Figure 5 results and leads Figure 2 with its size contrast. The approximately 5–12-fold statement is qualified to reference cue area 1 and cases with resolved times at both lengths: the unchanged table has 18 such comparisons across all five models, spanning 4.84–11.68. It is not asserted for all cue areas. The Figure 2 narrative uses the same qualification. The final PDF is checked for isolated source-link pages after the legend edits. |
| 32. Figure 3 legend referents | Accepted. All-model terminal gain folds accompany the maps in C; F describes its Wave-pinning time courses. The amplitude-limit result and size folds are in a concise main-text bullet, as required by the brief, with full brackets in methods. The unassigned trailing caption paragraph is removed. |

The inherited scope qualifications remain: finite-follow-up takeover does not
identify lifetime fate (14), the different size grids are declared (18), the
requested model facets are retained (25), and Figure 1's unrequested duration
analysis is not added (26). The outstanding per-WGD magnitude comparison,
detector sensitivity, mechanical-state mechanism and submission-level supplement
curation are not represented as completed results. The incoming narrative is
preserved verbatim; current panel references and qualified claims are given in
the revised outline.

[stencil]: ../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4/supporting/data/mode_slope_stencil_sensitivity_summary.csv
[phase]: ../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4/supporting/data/mode1_complex_phase_rate_comparison.csv
[runtime]: ../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4/supporting/metadata/runtime_provenance_recovery.json
[kernel]: ../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4/supporting/metadata/native_default_equivalence.json
