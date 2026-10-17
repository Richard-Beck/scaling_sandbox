# Vision-aligned manuscript outline

This package develops the argument in [MANUSCRIPT_VISION.md](../MANUSCRIPT_VISION.md)
using existing comparison and migration results, the delivered mechanism checkpoint,
and a prominent reserved section for new WGD motility experiments.

The six proposed main figures cover the five-model comparison, all-five linear
mechanisms, supported rings, replicated noisy-ring activity memory, literature
migration evidence, and new WGD motility results. Each figure has a detailed
legend, narrative conclusions, proposed measurements and inference limits. Main
figure numbering is proposed; revised mechanism artwork changes only title
numerals. Mechanism-linked measurement proposals remain beside
their motivating theoretical results. The new motility experiments have their
own major results placeholder and supply no invented observations.

## Deliverables

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

## Rebuild

Use the existing `scaling_sandbox` Conda environment; Python, NumPy, Matplotlib,
Pillow and pypdf are sufficient. From the repository root:

```sh
python manuscript_outline/draw_experimental_placeholder.py
python manuscript_outline/prepare_outline.py --mechanism-sections manuscript_outline/mechanism_sections_v1.json
python scripts/build_manuscript_outline.py --input manuscript_outline/outline_v1.json --output manuscript_outline/outline_v1.pdf
python -B manuscript_outline/prepare_revision.py
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v2.json --output manuscript_outline/outline_v2.pdf
```

The renderer uses only retained assets and prose; it launches no simulations.
The adjacent build manifest records source and image SHA-256 identities,
page/asset mapping, PDF identity, links and layout checks. Rendering timestamps
can differ on rebuild, so exact PDF-byte identity is not required. The scientific
figure files themselves remain unchanged unless a revision is explicitly named.

Bulk native histories and raw images are outside the compact checkout. Source
code, compact summaries, protocol/audit records and selected scientific assets
are retained for review and outline regeneration. Re-running integrations is a
separate SLURM workflow, not a prerequisite for rebuilding this document.
