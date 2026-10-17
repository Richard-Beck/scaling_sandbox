# Figure 3 persistence–noise surfaces

The paired surfaces show conditional sensitivity around the selected primary
Gaussian-hierarchy fit (`fit_05.json`). WGD− and WGD+ each have their own
panel. Within a panel, the plotted coordinates are the geometric-mean
persistence time and velocity-noise strength for that condition. Each grid
point shifts all five study-specific log parameters by the same offsets,
preserving their fitted between-study deviations. The other condition and
Hosny signed drifts remain fixed.

The data-misfit score is recomputed from the existing analytic predictions and
the same scaled residuals used during fitting. Hosny contributes median rows
only, matching the primary fit. The companion objective adds the original
primary pooling and persistence-prior penalties. Both panels show change from
the selected fit score; negative data-score changes are retained. The
geometric-mean markers are descriptive condition means, not separately fitted
optima. This is a conditional sensitivity surface, not a refit, likelihood
surface, confidence region, or uncertainty interval.

Each condition uses its own physical axes, centered on its fitted geometric
means, so both local valleys remain visible despite the groups occupying
different parameter ranges. Persistence values outside any study-level model
bound are masked. Noise strengths are also checked against the model bounds.
Surface height and color use signed asinh(change / 0.25) to
keep both the fitted valley and distant poor fits visible; the CSV stores raw
scores. The JSON records input hashes, selected-fit provenance, axis ranges,
and score/bound validation.

Reproduce from the repository root in the existing Conda environment:

```sh
conda run -n scaling_sandbox python reports/figure_revision/make_figure3_surface.py
```

Outputs are the 3D PNG/PDF, a paired heatmap PNG/PDF, the raw grid CSV, and
the metadata JSON under `reports/figure_revision/figures/`.
