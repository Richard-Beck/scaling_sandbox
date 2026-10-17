# Migration and fresh reproduction — 2026-10-01

This document records the earlier v8 migration and exact reproduction. Its
original `results/` artifacts now live in
`../archive/four_assay_v8_pre_normalization_20261001/results/`. The active
`results/` contains the normalized two-assay revision-9 rollout; see README.md
and its paired readout-preservation validation for the current protocol.

The accepted revised-v8 reference is preserved in `baseline/`. A fresh complete
campaign using the new pipeline produced `results/`. Neither older trajectories
nor prepared checkpoints were reused for that campaign.

| Verification | Result |
|---|---|
| Scientific model files | 7 |
| Shared simulation/orchestration scripts | 1 (`runner.py`) |
| Plotting scripts | 1 (`plot.py`) |
| Fresh preparations | 441 |
| Fresh response jobs | 4096 |
| Cases compared | 770 |
| Finite-duration tests compared | 34,092 |
| Full-state settling checks compared | 7,273 |
| Eligible OLS series | 66 |
| Endpoint, settling, flag and regression table differences | **0** |
| Finite tests, periodic fits, preparation checks/eligibility | **Exact equality** |
| SLURM preparation/response/controller/finalization records | 4539 completed, all exit 0 |

`results/comparison_to_baseline.json` compares all four tables plus every finite
test, periodic fit and preparation phase/eligibility diagnostic. It uses numeric
tolerance rtol=1e-9, atol=1e-10; all table numeric differences are exactly zero.
`results/exact_scientific_equality.json` additionally confirms exact decoded
record equality for all finite tests, periodic fit diagnostics, settling checks
and eligibility diagnostics. Source hashes and baseline snapshot hashes were
also rechecked after relocation. A redraw and comparison after removing the
legacy reports tree verifies the new pipeline's independence from it.

The revised outcome counts remain 498 success, 161 censored, 48 unreliable,
35 ineligible, 7 unsettled and 21 paused. The 161 censored cases include the
19 distinct Goryachev resolved-at/near-withdrawal classifications. Spring
competition was never scheduled. No amplitude, cue, time horizon, initialization,
solver tolerance, mesh or endpoint criterion changed during migration.

## Source consolidation

| Active file | Retained implementation |
|---|---|
| `models/wave_pinning.py`, `otsuji.py`, `goryachev.py`, `legi.py` | Canonical model classes extracted from `src/polarity1d/jilkine.py`; generic metadata/grid/integrator reused from the package |
| `models/holmes_model3.py` | Selected corrected Eq. 5/6 class and homogeneous root initialization; required retained periodic training helpers |
| `models/debelly.py` | Exact v8 full-state checkpoint solver, required cue/measurement helpers; retains raw moving-boundary state |
| `models/spring.py` | Author Fig4c activation/force functions, checked moving-boundary RHS and integration; necessary initial-state helpers |
| `runner.py` | V8 preparation/restart engine, responses, exact retained Fourier scoring and periodic preparation, consolidated campaign orchestration |
| `plot.py` | Accepted revised-v8 mode coloring, regressions and plots; compact-record loading replaces revision-directory lookup |

The spring functions originate from Zmurchok et al. (2020), author repository
`zmurchok/neutrophil-tension`, commit f78b05. The author's GPLv3 license is
preserved as `models/SPRING_LICENSE.txt`; the original source and its source.json
remain in reports/trial_report_assets/zmurchok2020. No runtime access to those
report files is needed.

## Archive locations and recorded paths

- Old `reports/four_test_suite/` moved intact to
  `polarity1d/archive/four_test_suite_legacy/`. This includes original revision6,
  rejected revision7, accepted revision8, old adapters, diagnostics, frozen
  provenance and approved assay-explanation-v2 source/assets. No legacy outputs
  were deleted. Former report status is preserved as REPORT_STATUS_BEFORE_MIGRATION.md.
- Fresh worker scratch moved to
  `polarity1d/archive/four_assay_reproduction_20261001/`. It retains complete
  restart checkpoints, scalar worker JSON, diagnostic NPZs, manifests, logs and
  frozen source provenance. `four_assay/work/` contains only the final status and
  archive pointer.
- Original provenance/path strings are preserved. To locate an old legacy
  record's referenced suite asset, replace its `reports/four_test_suite/` prefix
  with `polarity1d/archive/four_test_suite_legacy/`. For fresh worker references,
  replace `polarity1d/four_assay/work/` with
  `polarity1d/archive/four_assay_reproduction_20261001/`.

Historical scripts retain their original relative paths for evidence; the
compact runner is the supported reproduction entry point. The original WP1
Figure1_WP1_cross_model and other polarity1d campaigns remain untouched.

Remaining model-specific protocol differences are documented in README.md.
Numerical equality validates the migration; it does not settle the biological
interpretation of length effects or supply values for censored observations.
