# Figure 3 handoff

## Final outputs

The approved integrated figure is [figure3_updated.png](figure3_updated.png). The canonical export, [figures/figure3.png](../figure3.png), is an identical copy; both have SHA-256 `f631057b15b298cc72ecd17d5921a5e546711ac5a127b7bf6509be2effc0cc3d`. The full caption is [figure3_updated_caption.md](figure3_updated_caption.md), and source hashes, fit identity, observation mapping, track checks, and surface source are recorded in [figure3_updated_provenance.json](figure3_updated_provenance.json). The canonical observation table contains all nine plotted comparisons in [figure3_observations.csv](../figure3_observations.csv); persistence scenario folds, including the pooled geometric mean row used in panel E, are in [figure3_persistence.csv](../figure3_persistence.csv), with primary-fit noise values in [figure3_primary_noise.csv](../figure3_primary_noise.csv).

The standalone Hosny sensitivity supplement is [hosny_sensitivity.png](hosny_sensitivity/hosny_sensitivity.png), with [caption](hosny_sensitivity/CAPTION.md), [renderer provenance](hosny_sensitivity/figure_provenance.json), and [analysis results/provenance](hosny_sensitivity/hosny_analysis.json). Its principal data files are [relative-surface grid](hosny_sensitivity/hosny_surface.csv), [focused surface grid](hosny_sensitivity/hosny_surface_zoom.csv), [drift profiles](hosny_sensitivity/hosny_drift_profiles.csv), and [leave-Hosny-out comparisons](hosny_sensitivity/hosny_leaveout_comparison.csv). [Relative-surface methods](SURFACE_METHODS.md) documents the integrated Figure 3 relative surface; the standalone Hosny methods and renderer notes are in [hosny_sensitivity/README.md](../../hosny_sensitivity/README.md) and [hosny_sensitivity/METHODS.md](../../hosny_sensitivity/METHODS.md).

The user-selected Option 3 track figure remains available as [figure3_annotated_tracks.png](figure3_annotated_tracks.png), with [track methods](TRACK_METHODS.md) and [track provenance](figure3_tracks_provenance.json). The broader relative persistence/noise surface is also retained separately as [figure3_relative_persistence_noise_surface.png](figure3_relative_persistence_noise_surface.png), with its [data](figure3_relative_persistence_noise_surface.csv) and [provenance](figure3_relative_persistence_noise_surface.json). The original absolute persistence/noise surface is [figure3_persistence_noise_surface.png](figure3_persistence_noise_surface.png), with matching [heatmap](figure3_persistence_noise_surface_heatmap.png) and [data](figure3_persistence_noise_surface.csv); these are reference outputs, not substitutes for the integrated figure. Earlier alternatives ([comparison.png](comparison.png), [option1.png](option1.png), [option2.png](option2.png), [option3.png](option3.png)) and the prior canonical figure [figure3_previous.png](figure3_previous.png) are review/reference material, not the final primary Figure 3.

Only PNGs were produced for the new Figure 3 and Hosny outputs. The existing `figures/figure3.pdf` predates this update (September 28) and is stale; it was left untouched. No report renderer was run. Other report figures and supplements were outside this work and should not be considered newly integrated or validated by this handoff.

## What changed and what was checked

The integrated figure uses the exact selected five-study Gaussian OU primary fit, `fit_05.json`. Panel A combines the former observed-versus-predicted A–C comparisons into one shared nine-row panel. It includes both Jemaa 2017 and Jemaa 2023 speed observations, alongside Xuan speed and MSD, Mallin path length and directness, and Hosny path length, directness, and signed migration index. Hosny uses group medians in the fit residual mask; its quartiles are not converted into fold error bars. The selected fit source hash and all nine source-to-row mappings are recorded in provenance.

Panel B shows the same representative model-generated Xuan tracks selected from 32 existing simulations per WGD group. They are median-near 12-hour paths sampled every 10 minutes, not observed individual cells. Speed is path length divided by 12 hours; directness is net displacement divided by path length. The fitted MSD curve uses the campaign’s authoritative integrated-OU displacement helper. During review, a missing tau factor in a draft MSD calculation was corrected; the displayed 11.111-hour values now reproduce the selected fit, approximately 2,357 µm² for WGD− and 4,504 µm² for WGD+.

Panel C reuses the precomputed relative surface: WGD− parameters are held fixed while common shifts vary WGD+ persistence and noise; the data-misfit surface and regularized total-objective surface are shown separately. The equal-geometric-means point is (1, 1), and the selected fit is approximately (15.158, 0.1381). Panel D retains paired primary-fit persistence and velocity-noise values. Panel E shows the actual WGD+ / WGD− persistence fold for each study and its pooled geometric mean across primary, 5-minute prior, 20-minute prior, strong-pooling, and weak-pooling fits. The pooled row is present in the CSV export as well as the heatmap.

The Hosny supplement profiles both signed drifts and their contrast while refitting the other parameters, and compares the full five-study fit, the same fit scored on the remaining four studies, and a four-study refit. Removing Hosny also removes its data, pooling, and prior contributions in the refit; the same-four-study composition control separates this change from the composition effect. Stable scaled-Bessel Rice calculations resolved failures in the drift-tail region, with a normal asymptotic fallback only at extreme drift-to-step-noise ratios. The broad relative surface has 3,969 points; a focused 3,721-point grid improves local contour resolution. All 288 profile starts and all eight leave-Hosny-out starts converged. Contour/profile Δ values describe sensitivity for the heuristic weighted objective; they are not likelihood-based confidence intervals. Sampled one-dimensional profile fits do not form a joint uncertainty region.

The integrated builder was corrected during review to preserve both primary-fit parameters in panel D, restore actual scenario persistence folds and the pooled row in panel E, and prevent track annotations from covering paths. Panel colors, surface labels, color scaling, and plot spacing were also revised after visual inspection. Figure 3’s normal selector now invokes the integrated PNG builder, and the provenance writer merges existing figure entries. Its validator checks Figure 3 as PNG-only. The new outputs do not generate or update a PDF.

## Reproduction

Run from the repository root with the existing `scaling_sandbox` Conda environment; no new environment or dependency installation is required. To rebuild the integrated primary figure and its canonical PNG, observation/persistence/noise exports, caption, and provenance:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/build_figures.py 3
```

The selected-track candidates and their support outputs can be regenerated with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/figure3_tracks.py
```

The integrated relative surface is generated by `reports/figure_revision/make_figure3_relative_surface.py`; the companion absolute-surface PNG and data are generated by `reports/figure_revision/make_figure3_surface.py`. These are separate support producers, not required to redraw the integrated figure from saved data.

To recompute the standalone Hosny analysis, then its focused grid, and finally render the PNG from the resulting CSV/JSON files:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/hosny_sensitivity/analyze.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/hosny_sensitivity/refine_surface.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/hosny_sensitivity/build_figure.py
```

The Hosny renderer alone can be rerun from saved analysis outputs with the final command. Analysis and focused-grid steps recompute numerical results; the archived provenance records input hashes, evaluated grid sizes, optimizer convergence, and validation summaries.
