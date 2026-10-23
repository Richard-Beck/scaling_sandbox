# Figure 3 round 3 manuscript integration

[Revised outline_v7](../outline_v7.pdf) · [Figure 3 and legend](figure3_with_legend.pdf) · [reader](index.html) · [editable outline](../outline_v7.json) · [main plates 2–4](figures_2_4.pdf).

Figure 3 replaces the round-2 main with sixteen native Wave-pinning profiles and continuous-cue half-heading time and mean angular speed across four models and three sizes. Only four fixed-parameter relaxation protocols were run. Example iii withdrew at 20,000; all four cue-free phases reached the stationarity stop. [Complete supplement](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3/figure3_round3_supplement.pdf) retains every round-2 main and supplement page, plus 90/99% heading-time curves. [Native source data and sampling qualifications](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3/README.md).

Figures 1, 2, 4, 5 and 6 preserve their exact outline_v6 sections and artwork; Figure 2 retains its completed round-3 edition. outline_v6 and all preceding files remain intact. No timed refresh. Data generation: 2026-10-22T19:48:39.746127+00:00; producer complete.

Rebuild with the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py collect
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py render
python -B manuscript_outline/prepare_figure3_round3_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v7.json --output manuscript_outline/outline_v7.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v7.json --report manuscript_outline/outline_v7.validation.json
python -B manuscript_outline/verify_figure3_round3_20261022.py
```
