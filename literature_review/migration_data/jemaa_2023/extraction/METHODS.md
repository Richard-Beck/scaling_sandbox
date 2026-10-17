# Jemaà 2023 migration digitization

Regenerate from the repository root using the existing cpose environment:

```sh
conda run --no-capture-output -n cpose python literature_review/migration_data/jemaa_2023/extraction/digitize.py
```

Read `protocol.json` for immutable source hashes, publication identity, native crop coordinates, software versions and calibration. `readings.json` stores the audited clone identities, mean-line rows, visible error-cap rows and dot-selection settings. `digitized_points.csv` retains native full-image pixel positions alongside calibrated speeds. `measurements.csv` preserves clone summary lines and error-bar conventions. `component_qc.csv` records accepted isolated components and unresolved fragments/clusters. `panel_crop.png` is the native review crop; `digitization_overlay.png` enlarges it fourfold and marks accepted dots (blue circles), means (blue lines), SEM endpoints when visible (purple), and excluded clusters (gray boxes).

The source figure and caption were visually inspected before extraction; the resulting overlay was reviewed against the native raster. Color segmentation identifies red/orange or green glyphs in clone-specific columns. Components qualify only inside explicit area/width/height bounds, with shape and boundary checks. Components adjoining the mean line are excluded. Merged regions are not split, assigned multiplicities or converted into invented individuals. The partial isolated-dot subset is biased toward uncrowded regions: do not estimate a population mean, variance or exact cell count from it. Horizontal jitter carries no physical meaning. Separate clone dots are not paired trajectories.

Speeds use a linear vertical axis; the original-image origin is at upper left with y increasing downward. Axis anchors and crop bounds are saved in the protocol. Reading allowance is one native pixel; axis-anchor allowance is 0.5 pixel, propagated as a separate worst-case calibration sensitivity. Neither allowance is biological SEM. The source hashes, clone mapping, algebraic calibration and regenerated outputs were checked; these checks do not certify full distribution recovery.

Figure 3E provides RKO diploid/tetraploid clone context; 3F provides the speed comparison for D1/D2/T1/T2. Axis anchors: 0 µm/h at row 948.5 and 40 µm/h at row 746.5 (about 0.198 µm/h per pixel). Mean-line readings are 2.5743, 2.9703, 6.3366 and 6.7327 µm/h, respectively. The caption states SEM and n = 5; individual SEM cap endpoints cannot be resolved from the raster and remain blank, not zero. The green diploid point clusters cannot be separated into exact individual values. Six isolated tetraploid dots are retained; along with four means, the table has ten rows. Gray cluster boxes identify excluded components without inferred dot counts.

Methods describe a 24-hour observation with hourly imaging. The downloaded two-page supplementary PDF supplies characterization, not clone migration trajectories; no raw migration table was recovered. Summary means and partial dots remain distinct evidence.

These packages do not update manuscript observations, fit OU models, or infer persistence from speed. Published biological replication remains distinct from plotted individual-cell observations and extraction precision.
