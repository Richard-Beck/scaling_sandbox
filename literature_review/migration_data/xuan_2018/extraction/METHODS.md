# Xuan 2018 extraction

## Regeneration and retained evidence

From the repository root:

```sh
conda run --no-capture-output -n cpose python literature_review/migration_data/xuan_2018/extraction/extract.py
```

The existing `cpose` environment provides OpenCV, NumPy, SciPy, Pillow and Matplotlib. The script uses one OpenCV thread and no segmentation model. `protocol.json` records versions, input SHA-256 hashes, axis anchors, scale-bar endpoints, tracking settings and timing scenarios. Originals remain immutable. `review/` contains annotated contact sheets covering every decoded frame; `tracking_overlay.mp4` is a review rendition, not the tracking input. The initial exploratory contact sheets were removed after the final overlays were generated.

## Movie calibration and timing

All 74 AVI frames decode at 1040 × 1392 pixels; no consecutive exact duplicate frames were found. Frame numbers are zero based. The visible 25 µm white bar occupies columns 910–987, giving edge endpoints 909.5 and 987.5 at row 1291: 78 pixels, or 0.320513 µm/pixel. A one-pixel allowance on each endpoint implies approximately 1.8% quadrature scale sensitivity. This is reading uncertainty, not biological variability. Image coordinates have x rightward and y downward, with the upper-left corner as origin. `scale_calibration.png` preserves marked endpoints.

Main motility methods specify 10-minute intervals for 12 hours; Fig. 4D also describes a 12-hour assay. The supplement labels this movie 16 hours. Its 73 intervals do not exactly reproduce a 12-hour endpoint span at 10-minute cadence. There is no recovered acquisition timestamp metadata resolving the discrepancy. Thus `experimental_time_h` is blank in `tracks.csv`. Conditional speeds use three explicitly named, equally spaced timing scenarios: methods cadence (10 minutes; 12 h 10 min total), 12-hour endpoint span, and 16-hour endpoint span. None is asserted as the movie's established timing. Playback fps is not used for experimental time.

## Tracking, visual review and exclusions

Gaussian smoothing (9 × 9), intensity threshold 45, and a 220-pixel minimum component area isolate fluorescent nuclei. The scale bar and label region are masked. Centroids are intensity weighted above background 30. Consecutive detections are linked by Hungarian assignment with a 60-pixel gate and area ratio between 0.5 and 2; nearly competing assignments are rejected. Missing or ambiguous links terminate a segment. No gap interpolation or silent reconnection is performed. Segment IDs are stable within a segment; different IDs may be the same biological cell after an identity break. The 78 segments are not 78 independent cells.

All 74 annotated frames were visually inspected in the 13 contact sheets. Eight full-window tracks were accepted: N010 (apparent giant nucleus), N015, N016, N017, N018, N021, N024 and N025. These supply `equal_window_comparison.csv`, `measurements.csv`, `msd.csv` and the displacement plot. N015 and N017 approach closely late in the movie but remain visibly separate; their centroid paths remain distinct. N010 remains an identifiable large, bright, irregular nuclear object throughout. Nuclear appearance alone does not establish ploidy, so the smaller nuclei are not labeled experimentally proven diploid/non-PGCC cells.

Other segments are retained as candidates for further review, not accepted measurements. Examples requiring caution include N001/N002 merging near the upper edge early in the movie; N008 merging/division-like changes in the upper-right cluster; N013 with splitting/nearby objects; N019 with an apparent division around frames 12–14; and boundary-censored N011/N014. Exclusion of these full tracks avoids treating merged centroids or daughter identities as one nuclear trajectory. `track_qc.csv` records segment coverage and status. No manual coordinate edits or drift subtraction were applied. There are no convincing fixed specimen features with which to estimate stage drift; the overlaid scale bar is fixed graphics, not a stage reference. Population displacement was not used as a drift proxy.

The 3-pixel localization allowance in the track table is a working sensitivity allowance, not an experimentally estimated confidence interval. Nuclear deformation changes the intensity centroid even at fixed position; segmentation threshold and photobleaching can also shift it. A 3-pixel position allowance is about 0.96 µm; differences of two positions have about 1.36 µm quadrature sensitivity. This limits interpretation of small single-frame steps. The movie does not identify finer cadence behavior or an OU persistence time.

## Metrics and scope

Path length sums successive centroid distances. Net displacement is the endpoint distance. Directness is net displacement divided by path length. Sampled step is path length divided by 73 intervals; conditional sampled speed divides that step by the scenario interval duration. These are nuclear centroid motion measurements. MSD at lag k averages squared displacement over every valid within-track pair separated by k frames; `n_pairs` is explicitly 74 minus k. Long-lag estimates have few pairs and correlated errors. MSD lags and successive positions are not biological replicates.

Only equal frame windows 0–73 are compared. The apparent giant nucleus is one cell in one field; the seven smaller nuclei provide an illustrative within-field comparison, not a population WGD effect or a replacement for the publication's replicated results. No persistence model or downstream fit was run. Results should not be interpreted as final acquisition-timed speeds until the cadence discrepancy is resolved.

## Figure 4 digitization

The 2018 version-of-record PNG is retained in `figures/`. Panel B is mean speed (µm/h), C is mean fraction in the slow mode, D is the published MSD curve, and E is the logarithmic slope of MSD. E's dimensionless “persistence” is not an OU time constant. Caption and methods specify mean ± SEM; Fig. 4D reports n = 5 and 50+ cells. General statistics describe at least three experiments; precise B/C/E denominators are not independently assigned.

Bar tops and visible SEM caps were read on the full native image. `digitized_points.csv` retains each pixel position and transformed value. B calibration is 0 at row 469 and 30 at row 52; C shares these rows for 0 and 0.8; E uses 0 at row 1002 and 1.6 at row 586. Axis reading allowance is 1 pixel and glyph reading allowance 1.5 pixels, propagated separately from SEM.

Both complete visible MSD traces were extracted by saturated red/blue pixel masks, taking the median trace ordinate at each native x column. The legend is excluded. Logarithmic x anchors are column 206 at 1000 seconds and column 347 at 10000 seconds; logarithmic y anchors are row 1002 at 1 µm² and row 586 at 10000 µm². Panel crop and anchors are in the protocol. The 510 curve samples are samples of two plotted population curves, not independent observations or author trajectory data. The overlay was visually checked against the original: it follows the full blue and red paths. X reading uncertainty is approximately 3% for a 1.8-pixel combined allowance on the log axis; y uncertainty is recorded per point. No biological error envelope is inferred where none is displayed.

Speed bar readings imply about 57% greater non-PGCC speed, the slow-mode bars imply about 56% greater PGCC slow-mode fraction, and slope bars imply about 38% greater PGCC slope (paper text says approximately 35%). These comparisons are consistency checks on raster readings, not exact reported numerical datasets. Full-figure source hashes link all points back to the immutable asset.
