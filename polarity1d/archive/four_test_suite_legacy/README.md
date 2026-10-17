# Four-test model comparison

**Status 2026-09-29:** the suite is provisional and under ongoing length-response diagnosis. Read [current project status](../CURRENT_STATUS.md) and [the preparation steady-state audit](preparation_audit/README.md); many experiments are expected to rerun following audit review. Revision 5 is the current numerical snapshot, not the final consolidated scientific result.

For output locations and a short account of the revisions/debugging, see the [task handoff](TASK_HANDOFF.md).

The latest comparison is [revision 5](revision5/RESULTS_README.md), which refines finite-countercue duration sampling across all seven models without changing the length grid, cues, coefficients, endpoints, or periodic scoring. It retains all old finite trials, adds 31,802 new duration tests across 133 eligible cases, and records 21 ineligible cases with explicit skip reasons. The other 616 case rows are reused byte-identically. The finite-duration table includes all 33,367 evaluated old and new trials.

- [Seven-model comparison (PNG)](revision5/results/four_test_all_models.png) · [PDF](revision5/results/four_test_all_models.pdf)
- [Finite-countercue endpoints (PNG)](revision5/results/finite_countercue.png) · [PDF](revision5/results/finite_countercue.pdf)
- [Periodic harmonics (PNG)](revision5/results/periodic_harmonics.png) · [PDF](revision5/results/periodic_harmonics.pdf)
- [Endpoint table](revision5/results/endpoints.csv), [all finite-duration trials](revision5/results/finite_duration_tests.csv), and [periodic fit metrics](revision5/results/periodic_harmonics.csv)

Revision 4 is preserved in [`revision4/`](revision4/RESULTS_README.md) as the dense-length comparison. Revision 3 is preserved in [`revision3/`](revision3/RESULTS_README.md) as the original-length harmonic reanalysis. Revision 2 is preserved in [`revision2/`](revision2/RESULTS_README.md) as the measurement using per-switch latency. The original default-parameter first-pass report remains under [`results/`](results/) as a historical record; revisions 2–5 supersede its spring, Holmes, and deBelly endpoint assays.
