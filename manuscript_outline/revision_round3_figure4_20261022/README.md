# Figure 4 round 3 — current manuscript outline

[**Revised outline_v9**](../outline_v9.pdf) · [**Figure 4 and legend**](figure4_with_legend.pdf) · [**Reader**](index.html) · [**Supplement**](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4/figure4_supplement.pdf).

The compact campaign is complete: 260 conditions and 45,696 independent histories (96–384 per condition). Local size exponents resolve in 125/260 neighborhoods. Of the 260 single-condition rates, 128 meet the precision target, 36 exceed it at the seed cap, 8 have insufficient valid bootstrap draws and 88 have no eligible rate. The drift/relocation closure is estimable in 42 conditions; all are within 15% of observed decay. Goryachev has no estimable closure ratio. Missing estimates and cross-method disagreement remain explicit.

Figure 4 compares circumference-dependent heading decay with peak drift and random-relocation nulls. Model hues encode size, and the direct drift/jump decomposition distinguishes measured ingredients from null expectations. Figures 1, 2, 3, 5 and 6 retain their exact outline_v8 sections and artwork. Length represents WGD-associated enlargement; simulated polarity heading, observed migration summaries and conditional OU velocity-memory estimates remain distinct readouts. No hypothesis gates are regenerated.

[Scientific methods and data status](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4/supporting/metadata/scientific_readme.md) · [Integration receipt](integration.json).

Data generation: 2026-10-23T10:03:39.273853+00:00; manuscript producer finished. No timed refresh.

Rebuild the figure with its single entry point, then this integration:

```sh
python -B polarity1d/manuscript/code/mechanisms/rebuild_figure4_round3_20261022.py
python -B manuscript_outline/prepare_figure4_round3_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v9.json --output manuscript_outline/outline_v9.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v9.json --report manuscript_outline/outline_v9.validation.json
python -B manuscript_outline/verify_figure4_round3_20261022.py
```
