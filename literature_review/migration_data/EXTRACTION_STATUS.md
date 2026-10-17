# Extraction products and evidence limits

Generated October 3, 2026. Each study's `extraction/` directory holds its final review package. Source images and movies remain immutable in the study directories. The accepted scope is extraction and documentation; current manuscript observation tables and OU fits were not modified.

| Study | Products | Remaining limits |
| --- | --- | --- |
| [Xuan 2018](xuan_2018/extraction/METHODS.md) | Original AVI decoded without resampling; 1,904 pixel/calibrated detections in 78 conservative segments; eight visually reviewed full-window tracks (one apparent giant and seven smaller nuclei); QC, all-frame overlays, scale calibration, equal-window metrics and MSD pair counts. Both full Fig. 4D MSD traces plus B/C/E means and SEM caps digitized. | 12-hour methods versus 16-hour movie label unresolved. Physical times blank; speeds conditional on named timing scenarios. Segments are not independent cell counts. Only eight reviewed tracks supply comparisons. Ploidy inferred from appearance remains unconfirmed; one giant is illustrative, not replicated population evidence. |
| [Mallin 2023](mallin_2023/extraction/METHODS.md) | Fig. 2 accumulated distance, net distance and directness glyphs and central bars; visible paired size/metric points; overlays, coverage and unresolved regions; streaming movie frame audit and feasibility assessment. | Overlapping dots prevent exact full distributions; central-bar statistic unspecified. Streaming movies lack established spatial calibration and acquisition coverage; quantitative retracking is deferred. |
| [Jemaà 2017](jemaa_2017/extraction/METHODS.md) | Fig. 4C–D MFH152 Cl1/Cl2/Cl7/Cl8 identity mapping, clone means/SEM and isolated individual dots, calibrated overlays and component QC. | 20-hour text versus 24-hour caption discrepancy. Published replicate n does not equal visible dot count. Supplement movies do not establish access to these clone trajectories. |
| [Jemaà 2023](jemaa_2023/extraction/METHODS.md) | Fig. 3F RKO D1/D2/T1/T2 means, resolvable individual dots and exclusion records, calibrated overlay and protocol. | Dense overlaps prevent complete distributions; SEM endpoints unresolved. Characterization supplement provides no migration trajectories. |
| [Hosny](hosny_2026/extraction/METHODS.md) | Recovered preprint PDF/main figures; 75 baseline/perturbation migration/FMI median-and-quartile summaries; ten printed size correlations; 191 partial supplementary paired points; printed Rayleigh tests and size means; overlays and [version audit](hosny_2026/extraction/VERSION_COMPARISON.md). | Journal full text/supplement blocked despite OA metadata. Numerical version equivalence unverified. Distances rescaled to 16-hour equivalents; written FMI normalization and nuclear-area reference require caution. |

## Reading the outputs

Start with each `METHODS.md` and `protocol.json`; they identify source versions, software, axis/crop calibration, hashes, exact regeneration commands and assumptions. `digitized_points.csv` retains pixel readings with calibrated values and evidence types. `measurements.csv` gives named summary conventions; other study tables retain individual-point or movie QC. File schemas differ by assay to preserve the actual evidence. Do not concatenate unlike summaries or uncertainty intervals without an explicit mapping.

Reading/calibration uncertainty is separate from SEM, IQR or biological spread. Isolated dots are incomplete, visibility-selected subsets; their average is not a replacement for a printed group mean or median. Dots from separate marginal plots are not paired measurements. Curve samples, MSD lags and movie frames are not independent biological replicates. No numerical public trajectory source was recovered for these comparisons; digitization does not create exact source data.

## Integrity, regeneration and storage

`manifest.csv` includes extraction hashes and derivation links. Original relocated records remain intact. Run study commands in their methods, then:

```sh
python literature_review/migration_data/audit/validate_extractions.py
```

This checks required products, CSV row counts, JSON protocols, Xuan track key uniqueness and equal-window/MSD consistency, verifies retained source hashes, and refreshes the output inventory and `audit/extraction_validation.json`. Study-specific validation records document additional reading/calibration checks. Visual review is documented separately and cannot be replaced by these file checks.

Contact sheets and extraction overlays are retained review evidence, not temporary scratch files. Exploratory images, caches and transient work directories created during this extraction were removed; bulk originals remain in `source/`, `figures/` or `derived/`. No new Conda environment or compute-heavy campaign was needed.

Before downstream replacement, compare each proposed value with the existing observation table's assay, population, units, timing and statistic. Distinguish improved precision from improved biological coverage. Freeze the chosen inputs and retain previous observations/results before a separately authorized fit update. Hosny final-version comparison and unresolved movie acquisition metadata remain open evidence tasks, not successful extractions.
