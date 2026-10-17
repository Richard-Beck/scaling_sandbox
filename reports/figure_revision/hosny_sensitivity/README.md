# Hosny sensitivity figure

This standalone analysis asks how the fitted WGD persistence/noise contrast
changes when the Hosny baseline drift contribution is fixed or profiled. It is
not integrated into the main Figure 3.

`build_figure.py` creates `hosny_sensitivity.png` from the analysis outputs in
`../figures/figure3_review/hosny_sensitivity/`. The figure has six panels:

- fixed-drift and refitted-drift relative persistence/noise surfaces;
- a profiled drift-contrast curve and selected near-optimal profile fits;
- persistence-fold and noise-fold comparisons for the five-study fit, the
  same fit scored on the remaining four studies, and a refit excluding Hosny.

The Δloss 1, 4, and 9 guide lines and selected profile points are descriptive
fit-sensitivity markers. The objective uses heuristic residual scales and
regularization penalties; it is not a likelihood. The displayed profile points
are sampled constrained profiles, not a joint confidence region or full
parameter envelope. The figure therefore reports sensitivity, not confidence
intervals or probabilities.

Source analysis and results are separate from the figure renderer. Run the
renderer in the existing environment after the analysis CSVs and JSON are
available:

```sh
conda run -n scaling_sandbox python reports/figure_revision/hosny_sensitivity/build_figure.py
```

The output directory contains the PNG, a short caption, and renderer
provenance. The original surface and other Figure 3 outputs are not modified.
