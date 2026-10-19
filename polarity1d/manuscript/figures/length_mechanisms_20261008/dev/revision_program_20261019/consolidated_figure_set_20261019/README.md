# Consolidated Figures 2–4, 19 October 2026

[Main figures and standalone legends](figures_2_4_main.pdf) · [Comprehensive supplementary figures](figures_2_4_supplement.pdf) · [HTML reader](index.html) · [Editable legends](LEGENDS.md) · [Three-page legend PDF](legends.pdf).

The main PDF contains three unchanged qualified vector plates, each followed by its standalone legend. The 64-page supplement contains Figure 2 on pages 1–7, Figure 3 on pages 8–29 and Figure 4 on pages 30–64. Four separately linked 15-page Figure 3 native-profile atlases remain data appendices. The HTML reader links every source, detailed method and compact reconstruction guide. Figures 1, 5 and reserved Figure 6 are outside this assembly.

The [selection manifest](supporting/metadata/selection.json) explicitly records inputs and their independent qualification receipts. The [assembly record](supporting/metadata/assembly.json) records input/output SHA-256 hashes, exact page ranges and layout checks. Rebuilding launches no native integrations:

```sh
python manuscript_outline/assemble_final_figures_20261019.py --manifest polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/consolidated_figure_set_20261019/supporting/metadata/selection.json --output polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/consolidated_figure_set_20261019
```

Run from the repository root in the existing `scaling_sandbox` environment. Input changes are rejected until explicitly requalified and selected. The dated experimental [selection and interpretation](../BEST_EFFORT_FIGURES_20261019.md) distinguishes complete results, qualification limits and future research. This is a working manuscript figure set, not a completed manuscript or author acceptance.
