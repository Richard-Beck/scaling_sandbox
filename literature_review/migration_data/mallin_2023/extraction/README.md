# Mallin extraction package

Read [METHODS.md](METHODS.md) for interpretation and exact regeneration commands.

- [digitized_points.csv](digitized_points.csv): original pixels and calibrated values for resolved Figure 2 dots and paired size–motility scatter points; partial distributions.
- [measurements.csv](measurements.csv): six central bars with reading/calibration uncertainty; their mean/median definition is unspecified in the source.
- [coverage.csv](coverage.csv), [unresolved_regions.csv](unresolved_regions.csv): recovery and exclusions, with no inferred merged-dot multiplicity.
- [calibration.json](calibration.json), [protocol.json](protocol.json): axes, crop bounds, input hashes, video decoding audit and assumptions.
- [overlays/figure2_annotated.png](overlays/figure2_annotated.png): full native-size visual audit. Green = resolved dots, orange = unresolved regions, cyan = crops, blue = central bars.
- [overlays/movie1_contact_sheet.jpg](overlays/movie1_contact_sheet.jpg), [overlays/movie2_contact_sheet.jpg](overlays/movie2_contact_sheet.jpg): representative movie review images; playback times are not experimental times.
- [movie_frame_diagnostics.csv](movie_frame_diagnostics.csv), [movie_candidate_states.csv](movie_candidate_states.csv): candidate frame-state mapping, not calibrated tracks.
- [validation.json](validation.json), [threshold_sensitivity.csv](threshold_sensitivity.csv): internal consistency and raster sensitivity checks.

The sources are 2023 publication assets. Movies are publisher MP4 renditions rather than original WMVs. Original acquisition mapping and spatial calibration remain unresolved. No downstream fit or literature-observation table was changed.
