# Hosny preprint extraction

This package contains readings from bioRxiv 2025.09.03.673887v1, posted September 5, 2025, not verified journal-version measurements. It covers baseline hypoxic PACC/non-PACC and normoxic PC3 motility, parallel/perpendicular FMI, HIF-1alpha and RhoA knockdowns in hypoxia and normoxia, and supplementary size associations. No downstream fit was changed.

## Files and regeneration

- `digitized_points.csv`: 225 pixel readings, representing median, Q1 and Q3 for 75 group/metric summaries in Figures 2, 3, 5–10. `measurements.csv` repeats these calibrated values without pixel fields.
- `size_associations.csv`: ten printed Spearman correlations from Figures 2L–M, 3E–F, S3A–D and S4A–B. All are printed as nonsignificant; no exact correlation p-values are available.
- `size_metric_points.csv`: 191 visible isolated supplementary marker centers; incomplete paired cell-area/metric observations. Counts per panel and exclusions appear in `scatter_qc.csv`.
- `rayleigh_tests.csv`: fifteen printed group p-values. `reported_values.csv`: four caption-reported size means and SEMs.
- `protocol.json`: immutable source hashes, full-image crop bounds, linear calibration anchors and uncertainty settings. `overlays/` marks all measured medians in green and quartile endpoints in orange, or isolated supplementary dots with green circles. Red masks omit text from FMI scatters.
- `VERSION_COMPARISON.md`: access evidence, version differences and unresolved mappings.

From repository root:

```
/home/4473331/.conda/envs/scaling_sandbox/bin/python literature_review/migration_data/hosny_2026/extraction/extract.py
```

This existing Conda environment supplies Pillow, NumPy and SciPy. No environment was created and no compute-heavy job was required. Raw caption text can be regenerated with `pdftotext -layout literature_review/migration_data/hosny_2026/source/hosny_preprint_v1.pdf literature_review/migration_data/hosny_2026/derived/hosny_preprint_v1_layout.txt`. The retained `../derived/hosny_preprint_v1_layout.txt` is that derivative, not numerical source data.

## Procedure and visual review

Original downloaded JPEGs were read at native resolution. Main-panel median and interquartile endpoint pixels were read manually against the colored/gray summary bars, after reading the caption. Every main-panel extraction overlay and both supplementary overlays were inspected. These are one operator's reviewed readings, not an independent repeated digitization. Main values are linearly transformed between printed y-axis anchors; x pixels identify groups and have no physical calibration. Crop rectangles refer to full-image pixels and are metadata; source images are never resized for measurement.

Main reading sensitivity is ±3 pixels; each axis anchor has ±1 pixel uncertainty. Separate first-order worst-case calibration sensitivity is recorded. These bounds are extraction precision, not SD/SEM, confidence intervals or biological variability. Q1/Q3 are the plotted IQR endpoints, never error-on-median estimates. Native scatter-marker reading sensitivity is ±1.5 pixels and axis sensitivity ±1 pixel, propagated separately to cell area and the metric. Retained decimals support regeneration and do not imply that precision.

Supplementary dots were found by thresholding grayscale below 100 and accepting connected components of 8–42 pixels, width/height 3–8 pixels and aspect ratio 0.55–1.8. Panel borders, axes and FMI text regions were masked. Overlapping, oversized and ambiguous components were omitted; no splitting, interpolation or invented points was attempted. Some close markers can still be ambiguous at this resolution. The recovered counts (36, 33, 37, 38 in S3; 28, 19 in S4) are subsets of the plotted cells, with selection biased toward separated markers. Do not estimate correlations or population distributions from these subsets. The printed correlations are preferable. Neither point order nor plot proximity identifies the same cell across panels.

## Assay and interpretation

PC3-GFP cells were acquired every 0.5 hours during a 16-hour experiment, in self-generated radial oxygen gradients. Each migration group contains 60 cells. Baseline Figure 2 specifies 20 per independent experiment in three experiments; perturbation captions state 60 cells without an experiment count, so that field remains blank for perturbations. These cells are not 60 independent experimental replicates. The preprint rescales every trajectory's positions by 16/T for a track of duration T hours and excludes tracks with fewer than three frames (less than one hour). Thus reported distances can be extrapolated from shorter observations to a 16-hour equivalent; observation_hours=16 must not imply all cells were followed for the full interval. Dividing distance by 16 gives the paper's sampled mean path speed only under that rescaling convention, not a resolved short-time velocity.

Accumulated distance sums successive displacements; Euclidean distance is endpoint displacement; directness is their ratio. Parallel FMI aligns displacement with the local outward oxygen gradient; perpendicular FMI measures lateral motion. Directional bias in this gradient is distinct from spontaneous directional persistence. No OU time was inferred. Although FMI plots are dimensionless, the preprint's written FMI equations show projections without explicitly stating normalization by accumulated distance; the extraction follows the plotted values without silently repairing the methods.

Time-lapse PACC classification requires cell area at 16 hours ≥1500 um2 and growth ≥3 times initial area. A large non-growing cell remains a non-PACC. Fixed-cell nuclear classification uses ≥2.5 times parental nuclear area. Nuclear appearance or cell area alone cannot establish ploidy independently.

## Source inconsistencies and remaining limitations

Figure 5 and 8 captions refer to migration K–M while images label J–L; this package uses the actual image labels. Figure 9's caption misreferences FMI C/D although the plotted panels are B/C; actual B/C labels govern. S3 caption refers to A–F but embedded figure has A–D, and S4 caption refers to A–C but figure has A–B. Directness axes erroneously display um; the ratio is dimensionless. The main Figure 2 caption gives a parental nuclear reference of 116.55 um2 whereas Supplement S2 reports 165 um2 ±0.64 SEM. These reference values are inconsistent and must not be pooled or resolved by assumption.

Figure 8K is visually the same plotted Euclidean-distance distribution and summary as Figure 5K, despite differing perturbations; those readings were retained as version-specific plotted evidence, not treated as independently confirmed numerical data. Original trajectories, raw plot tables and individual track durations are unavailable. Author contact was not authorized and was not attempted. Journal-version figure changes remain unverified because full text and supplement were blocked.
