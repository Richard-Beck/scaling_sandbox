# Round-2 review interim (22 October 2026)

[**outline_v5**](outline_v5.pdf) · [**Round-2 reader**](revision_round2_20261022/index.html) · [**Figure 4 interim and data**](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4/index.html).

The new rate-based persistence analysis and old/new S0 decision tables are integrated as review-interim supplementary sheets. The brief requires author review here, before the Figure 4 main and the Figure 3/2 revisions. The preceding main artwork is retained and labelled; this is a partial work-package handoff. Figures 1, 5 and 6 remain identical to outline_v4.

Rebuild after the analysis and S0 producer finish, using the existing `scaling_sandbox` environment:

```sh
python -B manuscript_outline/prepare_revision_round2_20261022.py
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v5.json --output manuscript_outline/outline_v5.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v5.json --report manuscript_outline/outline_v5.validation.json
```

[Interim integration and protected-section receipt](revision_round2_20261022/integration.json). The 20 October edition below remains the preceding completed composition.

---

# Manuscript outline and revised figures

[**Integrated manuscript outline**](../manuscript_outline/outline_v4.pdf) · [**Revised Figures 2–4 and legends**](../manuscript_outline/revision_20261020/revised_figures_2_4.pdf) · [**Report reader**](../manuscript_outline/revision_20261020/index.html).

**20 October revisions:** Figures 2–4 have completed their native campaigns and independent scientific, source and display reviews. Figures 1, 5 and reserved Figure 6 retain their preceding content. The linked report is a working manuscript outline; it contains no invented experimental observations.

Figure 2 uses equal-area cue widths, initiation and pulse trains to distinguish response signatures. Holmes shares the narrow-cue penalty with Wave-pinning, so that signature identifies a model class. Calibration-unavailable Holmes trains and deBelly numerical sensitivities remain explicit in the [Figure 2 reader and complete supplement](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/signature_assays_20261020/index.html).

Figure 3 measures cue-driven handover and cue-free competition versus physical separation at fixed circumference. Complete active-region amounts include challenger activity beyond the cue. Distance effects vary by model, and handover differs from a complete turn of the measured activity heading. All bounds, pulse brackets and native histories are linked in the [Figure 3 reader](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance/index.html).

Figure 4 shows that a size advantage in heading persistence can coexist with multiple active regions. At input fraction 0.05, Wave-pinning, Holmes and Spring retain a persistence advantage across seven declared views, with pointwise whole-seed 95% uncertainty from 12 independent histories per size. These selected conditions reuse the same common-horizon cohort in both selectors. H1 and H3 remain inconclusive and H2 fails its overall criterion; local multiplicity is descriptive rather than a validated predictor. The [Figure 4 reader](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/signature_figure4_20261020/figure4_complete_reader.html) links the statistical core and all 13 native-atlas appendix volumes. Simulated polarity-heading persistence remains distinct from observed migration and conditional OU persistence estimates.

### Earlier 19 October consolidation

[**Consolidated Figures 2–4: main PDF and legends**](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/consolidated_figure_set_20261019/figures_2_4_main.pdf) · [**64-page comprehensive supplement**](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/consolidated_figure_set_20261019/figures_2_4_supplement.pdf) · [**HTML reader and native-profile appendices**](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/consolidated_figure_set_20261019/index.html).

[Figure selection and interpretation](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/BEST_EFFORT_FIGURES_20261019.md) · [Historical scientific record](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/RESULTS.md).

**19 October consolidation, superseded for Figures 2–4:** all selected experiments, main plates, legends and comprehensive supplements have passed independent scientific/source/display checks. Figure 2 covers five models and 97 qualified model–length states, three widths, five middle-width strength levels and five native withdrawal times: 3,492 core protocols plus two targeted deBelly withdrawals. It retains incumbent amount/extent/peak and links the all-length broadening response to a native deBelly example. All 53 failed finite field tails remain explicit.

Figure 3 reports all 645 stationary-cue protocols across a narrower circumference range, with initial internal/straddling/external footprints distinguished. The measured complex activity vector connects early reorientation rates to native field histories. All outcomes, time/duration results and four 15-page native-profile appendices are retained. Finite-tail status is 620 pass, eight Holmes fail and 17 unavailable.

Figure 4 reports the complete 936-case, nine-strength core campaign plus 26 separately analysed Goryachev transition cases. Its main combines short post-burn-in kymographs, heading memory with usable-pair coverage and spatial organization. All strength/size/window/control results remain in comprehensive supplements. Spring sampled-duration failures and the unresolved Goryachev 0.02–0.03 interval remain explicit. Portable redraws reproduce current scientific artwork without native archives.

Figures 1, 5 and reserved Figure 6 retain their preceding content. Earlier outline and review PDFs below are preserved development history; they do not contain this final consolidation. These are working manuscript resources, not a completed manuscript.


This package develops the argument in [MANUSCRIPT_VISION.md](../MANUSCRIPT_VISION.md)
using existing comparison and migration results, the delivered mechanism checkpoint,
and a prominent reserved section for new WGD motility experiments.

The six proposed main figures cover the five-model comparison, cue-shape and
timing signatures, fixed-size separation and competition, replicated noisy-ring
activity and heading persistence, literature migration evidence, and new WGD
motility results. Each figure has a detailed
legend, narrative conclusions, proposed measurements and inference limits. Main
figure numbering is proposed; the 20 October revisions replace Figures 2–4
with the completed, independently reviewed artwork. Mechanism-linked measurement proposals remain beside
their motivating theoretical results. The new motility experiments have their
own major results placeholder and supply no invented observations.

## Deliverables

- [Current integrated outline PDF](outline_v4.pdf) and [editable document](outline_v4.json).
- [Revised Figures 2–4 PDF](revision_20261020/revised_figures_2_4.pdf) and [reader with every supplement volume](revision_20261020/index.html).
- [Exact reviewed-source integration record](revision_20261020/integration.json) and [outline validation](outline_v4.validation.json).

- [Preceding integrated outline PDF](outline_v3.pdf) and [editable document](outline_v3.json):
  the earlier 19 October composition, superseded for Figures 2–4 by the consolidation above.
- [First integrated outline PDF](outline_v1.pdf) and [editable document](outline_v1.json):
  the initial integration preserves the delivered checkpoint figures.
- [Revised outline PDF](outline_v2.pdf) and [editable document](outline_v2.json):
  the critically assessed second edition, with expanded mechanistic explanations,
  physical-motion/coverage diagnostics and the requested migration-panel repair.
- [Revision decisions and rebuild instructions](REVISION_NOTES.md) and
  [independent coherence review](../reports/manuscript_coherence_review_v1_20261017.md).
- [Foundation and direct-user evidence audit](../reports/manuscript_foundation_audit_20261017.md)
  and [mechanism intake](../reports/manuscript_mechanism_intake_20261017.md) identify
  exact scientific sources, qualifications and author preferences.

The planned-experiment graphic is an explicit placeholder, not a simulated assay.
No completed manuscript or experimental-validation result is implied.

## Rebuild the current report

Use the existing `scaling_sandbox` Conda environment. These commands assemble
completed, reviewed portable exports; they launch no native integrations.

```sh
python -B manuscript_outline/prepare_revision_20261020.py
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v4.json --output manuscript_outline/outline_v4.pdf
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/revision_20261020/revised_figures_2_4.json --output manuscript_outline/revision_20261020/revised_figures_2_4.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v4.json --report manuscript_outline/outline_v4.validation.json --require-tracked
```

The intake checks exact final content, independent reviews, artwork, readers
and every native-atlas volume. Figures 1, 5 and 6 are preserved as identical
section content. Figure 4 keeps its minimum 12-inch printed width; the HTML
reader permits horizontal scrolling on smaller screens. Full native-field
archives remain on HPC with measured storage and exact public inventories.

### Earlier draft rebuilds

Use the existing `scaling_sandbox` Conda environment; Python, NumPy, Matplotlib,
Pillow and pypdf are sufficient. From the repository root:

```sh
python manuscript_outline/draw_experimental_placeholder.py
python manuscript_outline/prepare_outline.py --mechanism-sections manuscript_outline/mechanism_sections_v1.json
python scripts/build_manuscript_outline.py --input manuscript_outline/outline_v1.json --output manuscript_outline/outline_v1.pdf
python -B manuscript_outline/prepare_revision.py
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v2.json --output manuscript_outline/outline_v2.pdf
python -B manuscript_outline/prepare_revision_20261019.py
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v3.json --output manuscript_outline/outline_v3.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v3.json --report manuscript_outline/outline_v3.validation.json
python -B scripts/build_manuscript_outline_v2.py --input polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/supporting/metadata/revision_review.json --output polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/supporting/vectors/revised_figures_2_4.pdf
```

`manuscript_outline/build_distance_explorer_20261019.py` rebuilds the offline
distance-response explorer from retained analysis tables. Full scientific
reanalysis and new native integrations require the separately retained campaign
histories; they are distinct from document assembly.

The renderer uses only retained assets and prose; it launches no simulations.
The adjacent build manifest records source and image SHA-256 identities,
page/asset mapping, PDF identity, links and layout checks. Rendering timestamps
can differ on rebuild, so exact PDF-byte identity is not required. The scientific
figure files themselves remain unchanged unless a revision is explicitly named.

Bulk native histories and raw images are outside the compact checkout. Source
code, compact summaries, protocol/audit records and selected scientific assets
are retained for review and outline regeneration. Re-running integrations is a
separate SLURM workflow, not a prerequisite for rebuilding this document.
