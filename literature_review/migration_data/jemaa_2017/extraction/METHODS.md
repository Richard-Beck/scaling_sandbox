# Jemaà 2017 migration digitization

Regenerate from the repository root using the existing cpose environment:

```sh
conda run --no-capture-output -n cpose python literature_review/migration_data/jemaa_2017/extraction/digitize.py
```

Read `protocol.json` for immutable source hashes, publication identity, native crop coordinates, software versions and calibration. `readings.json` stores the audited clone identities, mean-line rows, visible error-cap rows and dot-selection settings. `digitized_points.csv` retains native full-image pixel positions alongside calibrated speeds. `measurements.csv` preserves clone summary lines and error-bar conventions. `component_qc.csv` records accepted isolated components and unresolved fragments/clusters. `panel_crop.png` is the native review crop; `digitization_overlay.png` enlarges it fourfold and marks accepted dots (blue circles), means (blue lines), SEM endpoints when visible (purple), and excluded clusters (gray boxes).

The source figure and caption were visually inspected before extraction; the resulting overlay was reviewed against the native raster. Color segmentation identifies red/orange or green glyphs in clone-specific columns. Components qualify only inside explicit area/width/height bounds, with shape and boundary checks. Components adjoining the mean line are excluded. Merged regions are not split, assigned multiplicities or converted into invented individuals. The partial isolated-dot subset is biased toward uncrowded regions: do not estimate a population mean, variance or exact cell count from it. Horizontal jitter carries no physical meaning. Separate clone dots are not paired trajectories.

Speeds use a linear vertical axis; the original-image origin is at upper left with y increasing downward. Axis anchors and crop bounds are saved in the protocol. Reading allowance is one native pixel; axis-anchor allowance is 0.5 pixel, propagated as a separate worst-case calibration sensitivity. Neither allowance is biological SEM. The source hashes, clone mapping, algebraic calibration and regenerated outputs were checked; these checks do not certify full distribution recovery.

Figure 4C establishes MFH152 clone identities: diploid Cl1/Cl2 and tetraploid Cl7/Cl8. Figure 4D is the quantitative speed panel. Axis anchors: 0 µm/h at row 522 and 80 µm/h at row 322 (0.4 µm/h per pixel). Mean-line readings are 16.2, 16.0, 22.6 and 28.2 µm/h, respectively; visible SEM endpoints are retained separately and appear slightly asymmetric because of raster resolution. Caption specifies mean ± SEM and n = 3; it is not the number of visible single-cell dots. Four, three, three and twelve isolated dots are resolvable by clone (22 total); 12 summary/cap positions bring the digitized table to 34 rows.

The caption describes 24 hours, whereas the main results describe 20 hours. This conflict remains explicitly unresolved; no trajectory cadence is assigned from it. Supplement movies concern other populations/assays and were not retracked as if they were these clone trajectories. No public clone-level numerical migration dataset was recovered.

These packages do not update manuscript observations, fit OU models, or infer persistence from speed. Published biological replication remains distinct from plotted individual-cell observations and extraction precision.
