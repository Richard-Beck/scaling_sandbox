<!-- completed Figure 4 integration 20261023 -->
**Current integrated edition (23 October):** [outline_v9](../outline_v9.pdf), [completed Figure 4 reader](../revision_round3_figure4_20261022/index.html), and [current Figures 2–4 plates](../revision_round3_figure4_20261022/figures_2_4.pdf). The Figure 2 integration record and its preceding-edition files are preserved below.

# Figure 2 — round 3 integration

[Revised outline_v6 PDF](../outline_v6.pdf) · [Editable outline](../outline_v6.json) · [Figure 2 and legend](figure2_with_legend.pdf) · [Reader](index.html) · [Main plates 2–4](figures_2_4.pdf)

The full revised Figure 2 replaces its preceding outline section and legend. A–C show size responses; D–E distinguish temporal classes; F gives per-class prospective biosensor predictions. The plate is 7.2 × 10 inches. All 450 width observations and 180 temporal coordinates are included. Broad Holmes cues already resolve every length, so no new simulations were needed. The deBelly L19/gap16 point was checked against native fields before plotting.

[Complete supplement](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure2/figure2_round3_supplement.pdf) preserves the superseded main, fate maps, ratios, tables and every earlier supplement page. [Source data and audit](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure2/README.md) include raw shares, visibility flags, exact endpoint checks and remaining numerical qualifications.

Figures 1, 3, 4, 5 and 6 are exact copies of their outline_v5 sections and use identical image files. The previous outline_v5 edition remains at its original path. [Integration receipt](integration.json).

Rebuild with the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure2_20261022.py
python -B polarity1d/manuscript/code/mechanisms/organize_length_dev.py --group revision_round3_20261022/figure2
python -B manuscript_outline/prepare_figure2_round3_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v6.json --output manuscript_outline/outline_v6.pdf
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/revision_round3_20261022/figure2_with_legend.json --output manuscript_outline/revision_round3_20261022/figure2_with_legend.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v6.json --report manuscript_outline/outline_v6.validation.json
python -B manuscript_outline/verify_figure2_round3_20261022.py
```

Generated UTC: 2026-10-22T18:55:52.599915+00:00. No timed browser refresh. No simulations or retuning are part of report assembly.
