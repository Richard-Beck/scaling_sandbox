# Figure 4: corrected event scoring and normalized drift

[**Revised outline_v10**](../outline_v10.pdf) · [**Figure 4 and legend**](figure4_with_legend.pdf) · [**Reader**](index.html) · [**Supplement**](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4/figure4_supplement.pdf)

Figure 4 now normalizes physical drift by the smallest size at the same input and estimates jump angles from scorable events. The existing 915,757 counted replacements include 321,423 scored angles. Event-based closure is estimable in 154 conditions; 57 are within 15% of observed decay. Among conditions with multiple fronts for most frames, 0/33 event-based and 3/33 detector-free comparisons are within 15%. The decomposition therefore retains its failures in the regime it seeks to explain. No new simulations were run.

Figures 1, 2, 3, 5 and 6 preserve their exact outline_v9 sections and artwork. Length represents WGD-associated enlargement; activity heading, migration summaries and conditional OU persistence estimates remain distinct evidence.

[Scientific methods and results](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4/supporting/metadata/scientific_readme.md). Generated 2026-10-23T13:03:14.134842+00:00; producer finished.

```sh
python -B polarity1d/manuscript/code/mechanisms/figure4_round4_reanalysis_20261023.py
python -B manuscript_outline/prepare_figure4_scoring_20261023.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v10.json --output manuscript_outline/outline_v10.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v10.json --report manuscript_outline/outline_v10.validation.json
python -B manuscript_outline/verify_figure4_scoring_20261023.py
```
