# Steady-state suite, revision 6

Generated: 2026-09-29T21:36:56.601346+00:00. Status: {"ineligible": 35, "success": 518, "censored": 142, "unreliable": 48, "unsettled": 6, "paused": 21}. Preparations available: 441/441.

These figures show revision-6 results only. Pending or unusable measurements are grayscale, with no distinct failure symbols. Figures are preliminary until all expected records and finite-duration chunks have completed. Execution failures are reported separately in data and do not count as demonstrated failure to settle.

[Full test suite](four_test_all_models.png) · [PDF](four_test_all_models.pdf) · [Settling companion](time_to_settle.png) · [PDF](time_to_settle.pdf).

Suite: seven rows, competition / gradient reversal / finite countercue / periodic amplitude / periodic phase. Periodic preparation and harmonic reliability scoring retain revision 3. Spring competition remains paused. Finite thresholds are the shortest passing tested duration in the full revision-5 dense grid; incomplete sweeps do not produce thresholds.

Companion: detected settling time for each cue phase, sampled at 100-unit checkpoints after a 500-unit confirmation window (deBelly finite no-cue: 1000). Cap 10000 per phase. All evolving chemistry and mechanical state must have <=0.1% spatial pointwise change and projected final drift over the window. RD/spring drift uses the full RHS; deBelly uses the final finite difference. Native time scales are not converted to seconds for uncalibrated models. Gray denotes unavailable/unsettled/inapplicable phases, including paused spring competition and spring finite conditioning (spontaneous no-cue preparation). A settled phase is plotted even if the assay subsequently fails polarity/peak eligibility. The no-cue column includes spring spontaneous preparation for both reversal tests; curves can overlap.

[Endpoints](endpoints.csv) · [Settling times](settling_times.csv). Full per-field diagnostics are in ../preparations; complete state checkpoints are in ../checkpoints. Figures have editable SVG exports and source build_report.py.
