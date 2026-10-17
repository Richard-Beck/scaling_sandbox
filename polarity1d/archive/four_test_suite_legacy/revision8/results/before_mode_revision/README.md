# Steady-state suite, revision 8

Generated: 2026-10-01T15:48:50.330131+00:00. Status: {"ineligible": 35, "success": 517, "censored": 142, "unreliable": 48, "unsettled": 7, "paused": 21}. Preparations available: 441/441.

These figures show revision-8 results only. Pending or unusable measurements are grayscale, with no distinct failure symbols. Figures are preliminary until all expected records and finite-duration chunks have completed. Execution failures are reported separately in data and do not count as demonstrated failure to settle.

[Full test suite](four_test_all_models.png) · [PDF](four_test_all_models.pdf) · [Settling companion](time_to_settle.png) · [PDF](time_to_settle.pdf).

Suite: seven rows, competition / gradient reversal / finite countercue / periodic amplitude / periodic phase. Periodic preparations and scoring retain revision 3; sampling is P/80 for every model. Countercue recovery is 10000 units for every model. Spring competition remains paused. Finite thresholds are the shortest passing tested duration in the full revision-5 dense grid; incomplete sweeps do not produce thresholds.

Companion: detected settling time for each cue phase, sampled at 100-unit checkpoints after a shared 1000-unit confirmation window. Cap 10000 per phase. All evolving chemistry and mechanical state must have <=0.1% spatial pointwise change and projected final drift over the window. RD/spring drift uses the full RHS; deBelly uses the final finite difference. Native time scales are not converted to seconds for uncalibrated models. Gray denotes unavailable/unsettled/inapplicable phases, including paused spring competition. A settled phase is plotted even if the assay subsequently fails polarity/peak eligibility. The no-cue column shows finite autonomous preparation; spring additionally shows its retained spontaneous gradient preparation. Spring finite conditioning now uses a left Gaussian.

[Endpoints](endpoints.csv) · [Settling times](settling_times.csv). Full per-field diagnostics are in ../preparations; complete state checkpoints are in ../checkpoints. Figures have editable SVG exports and source build_report.py.
