# Literature migration data extraction plan

The objective is to replace eyeball estimates with reproducible measurements and assess whether accessible movies provide additional migration information. The user designated Xuan 2018 retracking as high priority. This plan does not launch tracking, digitization or refitting, and does not change the polarity-reversal analysis being developed by another agent.

## Priority and expected outputs

| Priority | Study and task | Output |
| --- | --- | --- |
| 1 | Retrack the original Xuan movie, identifying the apparent single PACC/PGCC and multiple non-PACCs. | Audited pixel tracks, calibrated coordinates, annotations and an individual-cell migration comparison. |
| 2 | Digitize Xuan Figure 4, especially complete MSD curves and speed summaries. | Traceable curve points and summary measurements with extraction uncertainty. |
| 3 | Digitize Mallin Figure 2 and assess its two motility movies for retracking. | Migration distributions, size–motility point pairs where visible, and a movie feasibility assessment. |
| 4 | Recover Jemaà 2017 Figure 4 and Jemaà 2023 Figure 3 migration measurements. | Clone-specific speed summaries and resolvable individual values. |
| 5 | Resolve Hosny final-publication access, compare versions, and extract the relevant preprint main/supplementary plots. | Version comparison and better migration, FMI and size-association inputs. |

These priorities reflect the present assets and user direction. A newly recovered numerical dataset should take precedence over digitizing the corresponding plot.

## Xuan retracking

1. **Establish the input and calibration.** Use `xuan_2018/source/xuan2018_movie.avi`, verify its manifest hash and decode all 74 frames without resampling. Record frame indices, dimensions and any duplicated frames. Calibrate pixels from the visible 25 micrometer scale bar and save the marked endpoints. Resolve the supplement's 16-hour label versus the main methods' 12-hour duration; do not assign experimental cadence from 20 fps playback. If time remains unresolved, retain frame-index tracks and explicitly labeled timing scenarios rather than one assumed physical-time result.
2. **Identify cells and follow them.** Annotate the apparent large PACC/PGCC and all usable smaller nuclei with stable IDs. Save the rationale and uncertainty of each classification; nuclear appearance does not independently establish ploidy. Track nuclear centroids on original frames, excluding the scale bar. Use manual correction or a suitable existing tracking workflow, with visual review throughout. Record gaps, crossings, divisions, entry/exit, and exclusions. Do not silently interpolate through ambiguous identities.
3. **Review the tracks.** Produce overlays or a movie showing IDs and trajectories, plus a QC table and displacement-versus-frame traces. Review drift against observable image features; avoid subtracting population movement as if it were stage drift. Record any correction and its evidence. Preserve uncorrected coordinates.
4. **Measure motion.** Once timing is defensible, compute path length, net displacement, directness, sampled speed and MSD using explicit definitions. Compare equal observation windows and state the number of contributing displacement pairs at each MSD lag. Assess acquisition cadence and localization error before interpreting short-lag speed or directional memory. Match the paper's conventions where possible and document departures.
5. **Limit the claim to the data.** The apparent single PACC is one cell in one field, not a replicated population estimate. Successive frames or MSD lags are not independent biological replicates. Report this as an illustrative within-movie comparison; do not use it alone to replace published population summaries or infer a WGD effect. A persistence estimate should only be attempted if track duration and sampling support identification, with model dependence made explicit.

The first deliverable is a reviewable tracking package, not an updated OU fit. It should contain pixel coordinates even if the acquisition timing is still unresolved.

## Figure digitization

For each target, preserve the original downloaded figure and identify the publication version, figure/panel, condition, clone, units, sample size and reported summary statistic. Read the caption before choosing extraction targets; plotted means, medians, quartiles and fitted curves are different evidence.

Save the crop bounds, axis calibration points, linear/log scaling and annotated extraction overlay. Retain both pixel measurements and transformed values. Validate against printed values or independently repeated readings where possible. Record pixel-resolution and calibration uncertainty separately from biological variability. Error bars must retain their stated meaning; do not turn SEM, IQR or reading uncertainty into interchangeable intervals.

- **Xuan:** digitize both full MSD traces on their logarithmic axes, plus Figure 4B speed summaries and error bars. Figure 4E's “persistence” is the paper's MSD-derived measure; do not assume it is an OU correlation time. Points sampled from a plotted curve do not create independent observations.
- **Mallin:** extract accumulated distance, Euclidean distance and directness distributions, with plotted group summaries. Extract the paired size–metric points from scatterplots where separable. Separate marginal dot plots cannot recover per-cell pairing between different metrics. Mark merged or unresolved dots rather than inventing values.
- **Jemaà 2017:** use Figure 4C–D for the clone tracking comparison and preserve clone identities. Movies in the supplement concern other assays/populations. Their existence does not establish access to these clone trajectories.
- **Jemaà 2023:** use Figure 3F for clone-specific speeds; document overlapping dots and distinguish their recoverable distribution from exact source data. The supplementary PDF supplies characterization, not migration trajectories.
- **Hosny:** first compare the 2025 preprint with the September 2026 journal version when accessible. Preserve version-specific figure mappings. Extract the relevant baseline and perturbation migration/FMI summaries and supplementary size associations. Keep drift/gradient alignment distinct from spontaneous directional persistence. An author data request is an option, but contacting authors requires a separate user instruction.

## Mallin movie feasibility

The saved videos are publisher streaming renditions, not the advertised WMV originals. Establish spatial calibration, unique acquisition frames, cadence and experimental conditions before quantitative tracking. Inspect clustering and shape changes for centroid ambiguity, and determine which tracks have usable coverage. Do not assume the movies contain all cells or replicates underlying Figure 2. If streaming compression prevents reliable measurement, retain a feasibility report and prioritize figure digitization or obtaining the originals.

## Output records and reproducibility

Use each study's `extraction/` directory. Keep bulk working files under `extraction/work/`; put compact reviewed outputs and methods outside it. Suggested deliverables are:

| File | Required content |
| --- | --- |
| `protocol.json` | Input hashes, publication version, software/settings, calibration, timing evidence and assumptions. |
| `tracks.csv` | Study/movie, cell ID, frame index, experimental time if known, original pixel coordinates, calibrated coordinates and QC status. |
| `track_qc.csv` | Classification evidence, gaps, exclusions, corrections and usable intervals per cell. |
| `digitized_points.csv` | Study/version, figure/panel, condition/clone, series, pixel position, physical values, evidence type and uncertainty. |
| `measurements.csv` | Metric, units, observation window, cell/series identity, value, summary definition and source linkage. |
| `METHODS.md` | Extraction procedure, review findings, limitations and exact regeneration commands. |

Keep evidence types explicit: author-reported numerical value, digitized summary, digitized individual point, retracked representative movie, or model-derived estimate. Preserve source hashes and derivation links in the manifest as outputs are added.

## Integration into the manuscript analysis

Before any refit, map each proposed replacement to the current observation table and compare values, protocols and summary conventions. Record changes in extraction precision, scientific information and sample coverage separately. Freeze the reviewed extraction package, then make any downstream fit update as a separate analysis step, retaining the old inputs and results. No fit update is authorized merely by completing this storage task.

Check existing Conda environments before adding dependencies. Use SLURM for work meeting the repository's compute thresholds; do not start a full segmentation campaign on a login node. Current organization and checksum checks require no new environment.
