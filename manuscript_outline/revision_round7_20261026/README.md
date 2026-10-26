# Round 7 manuscript revision

[outline_v13](../outline_v13.pdf) · [Figures 1–6 with legends](figures1_6_with_legends.pdf) · [Reader](index.html).

Figures 3–5 follow the new conduction and relocation narrative. Figures 1–2
retain their exact outline sections, artwork and numerical evidence. The previous
migration Figure 5 is Figure 6; Figure 7 is reserved. Dropped panels and earlier
supplements are absent from this edition. Historical files are preserved.

Figure methods open with Unexpected. Compact CSVs, provenance and execution
records accompany the figures under `revision_round7_20261026`. The pilot PNGs
are checks, not manuscript panels. No output was tuned to pilot values.

Rebuild the individual figures first, then use the existing scaling_sandbox
environment to run `integrate_round7_20261026.py integrate`,
`build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v13.json --output manuscript_outline/outline_v13.pdf`,
and `integrate_round7_20261026.py validate` (paths are repository relative;
the integration script is in `polarity1d/manuscript/code/mechanisms/`).
