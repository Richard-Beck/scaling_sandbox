# Figure 3 round 3 manuscript integration

[Revised outline_v8](../outline_v8.pdf) · [Figure 3 and legend](figure3_with_legend.pdf) · [reader](index.html) · [editable outline](../outline_v8.json) · [main plates 2–4](figures_2_4.pdf).

Figure 3 places eight native Wave-pinning profiles at C = 30 on the left, with time to 50% heading (B) above mean angular speed (C) on the right. B and C retain all four models and three sizes; i–ii identify the displayed profile placements. This revision reuses the four previously completed relaxation protocols and all existing readouts; no new simulations were run. The C = 47.6 profiles and their 20,000-horizon qualification remain in the previous main and legend appended to the supplement. [Complete supplement](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3_single_C/figure3_round3_supplement.pdf) retains every round-2 main and supplement page, plus 90/99% heading-time curves. [Native source data and sampling qualifications](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3_single_C/README.md).

Figures 1, 2, 4, 5 and 6 preserve their exact outline_v7 sections and artwork; Figure 2 retains its completed round-3 edition. outline_v7 and all preceding files remain intact. No timed refresh. Data generation: 2026-10-22T20:05:54.597143+00:00; producer complete.

Rebuild with the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_single_C_20261022.py
python -B manuscript_outline/prepare_figure3_round3_single_C_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v8.json --output manuscript_outline/outline_v8.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v8.json --report manuscript_outline/outline_v8.validation.json
python -B manuscript_outline/verify_figure3_round3_single_C_20261022.py
```
