# Figure 3 relative persistence–noise loss surface

`figure3_relative_persistence_noise_surface.png` compares the WGD+ and WGD−
geometric means for persistence and noise. WGD− is held fixed. At each grid
point, all five WGD+ study-specific log persistence and noise parameters shift
by common offsets, preserving their fitted study deviations. Other parameters
and Hosny signed drifts remain fixed. Thus (1, 1), labeled **Equal group
means**, equates the two geometric group means; individual study parameters
can still differ.

The star marks the selected `fit_05.json` WGD+ point: persistence ratio
15.1583 and noise ratio 0.138107 (approximately 15.16 and 0.1381). The
conditional surface is not a refit or an uncertainty region. The data score
uses the fitting residual mask, including Hosny medians only. The regularized
total loss adds the existing primary pooling and persistence penalties. Each
panel reports change from the selected fit, whose scores are data misfit
11.473550794223446 and total loss 33.06008722062528. Negative data-score
changes are retained.

Both ratio axes are logarithmic. Color applies signed asinh(change / 0.25)
for visibility, while the colorbar ticks report raw score changes. The CSV
contains all valid and masked grid points and their raw scores. The JSON
records input hashes, selected-fit provenance, bounds, markers, and validation.

Reproduce from the repository root in the existing environment:

```sh
conda run -n scaling_sandbox python reports/figure_revision/make_figure3_relative_surface.py
```

This writes the PNG, CSV, and JSON to this folder. The companion prior
absolute-surface producer is `reports/figure_revision/make_figure3_surface.py`;
it now also writes PNG and support data here without creating PDFs.
