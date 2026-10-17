# Seven models, two assays: normalized gradient cues

The active suite runs gradient reversal and finite countercue at every existing
sampled length. User-authorized revision 9 replaces the four-assay v8 campaign.
Competition and periodic tracking are deferred and absent from the active task
list, runner, figures and slope tables. The directory name remains `four_assay`
so established paths to the active runner and figures continue to work.

- `models/`: the seven unchanged model implementations.
- `runner.py`: linear cues, native preparation, endpoint gates, SLURM execution,
  collection and paired baseline validation.
- `plot.py`: cumulative time facets and normalized component/cumulative slopes;
  PNG, PDF and SVG exports.
- `protocol.json`: 308 tasks (154 model–length combinations, two assays), original
  dense duration grids, normalization amplitudes and retained timing rules.
- `baseline/`: immutable accepted v8 records and reference artifacts.
- `results/`: current normalized two-assay outputs and paired validation.
- `work/`: current SLURM workers, checkpoints, diagnostic traces and logs.

The prior active results, work metadata and source snapshot are preserved in
`../archive/four_assay_v8_pre_normalization_20261001/`. The successful 58-case
normalization pilot is retained in `pilots/linear_shared_20261001/`; its frozen
source references describe the preceding v8 code and results, now archived.
`MIGRATION.md` records the earlier migration/reproduction, not this new protocol.

The completed full revision-9 sweep passed paired validation: 308 preparations,
308 assay cases and 34092 finite-duration tests. All 213 baseline-successful
readouts were preserved, with no losses or gains and unchanged classifications
for 70 censored, 19 ineligible and six unsettled cases. The unchanged controls
have zero mismatches, including 23280 pulse tests. All 34 eligible component
fits and 34 cumulative fits were independently checked against raw records.
Holmes pulse thresholds are unchanged; deBelly thresholds are 7–51% shorter and
Spring thresholds 8–42% shorter across their successful sampled lengths.

## Normalized cues and retained assay differences

Both assays use the same linear spatial cue within each model. With `q=x/L`,
left cue is `amplitude*(1-q)` and right cue is `amplitude*q`. Training and
countercue are aliases for these same left and right profiles. Left and right
are mirrors with equal amplitude; amplitudes are not equalized across models.

| Model group | Shared cue amplitude |
|---|---:|
| Wave-pinning, Otsuji, Goryachev, LEGI | 0.1 |
| Holmes 3 | 1.0 |
| deBelly | 3.6 |
| Spring | 4.0 |

Gradient reversal prepares a left-cued state, then uses the existing 20-unit
left-to-right crossfade and 20000-unit observation window. Spring retains its
preceding spontaneous preparation and reflection. Finite countercue prepares
under the same left cue, withdraws it until a qualified autonomous left state
settles, applies the matching right cue for each tested duration, then observes
10000 units of no-cue recovery. Its endpoint remains the shortest passing
*tested* pulse duration, with final 500-unit confirmation of autonomous right
polarity. It excludes recovery observation time.

Native initial states, solvers, meshes, units, response/settling gates and all
model-specific checks remain unchanged. Full-state preparation uses the original
1000-unit confirmation window, 10000-unit per-phase cap, 100-unit checks and
relative tolerance 0.001. deBelly retains its additional Rac and bound-MCA x50
checks. Holmes retains corrected Eq. 5/6 and basal Rac activation 3.5. Spring
retains nonlinear cue activation and dynamic physical domain length.

## Run, render and validate

Use the existing environment; no installation is required:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/4473331/.conda/envs/scaling_sandbox/bin/python
$PY polarity1d/four_assay/runner.py launch
$PY polarity1d/four_assay/runner.py status
```

The pipeline submits all 308 preparations, plans qualified response jobs in
10-duration chunks, then collects and renders complete outputs. All jobs use
SLURM QoS `small`, one CPU and 3 GB, with no array concurrency throttle. Every
worker checks frozen sources and checkpoint hashes. Missing output or incomplete
duration coverage is an execution error, not biological censoring. Preserve an
existing `work/` and `results/` before another full campaign; `baseline/` stays
intact. Custom `--work` and `--results` destinations propagate to SLURM workers.
The results landing page and status file show the current production stage and
generation time, without browser self-refresh.

To redraw complete results or rerun paired validation:

```bash
$PY polarity1d/four_assay/plot.py
$PY polarity1d/four_assay/runner.py compare
```

Validation joins records by model, assay and length, preserving the v8 task index
as `baseline_task_index`. Every sampled baseline-successful readout must remain
successful; new censoring, ineligibility or unsettled preparation is a loss,
without compensation by gains. All gradient-reversal cases and the four models
whose finite cues already matched are exact controls, including their full pulse
tests, settling checks and eligibility diagnostics. Intended finite-cue changes
in Holmes 3, deBelly and Spring may change numerical preparation/response times.
`comparison_to_baseline.csv` and `.json` distinguish those changes from failures.

## Figures and fits

`length_response_and_settling` has one column per retained assay and one row per
model. Each curve is cumulative: settled preparation stages add in chronological
order, and the final curve adds the successful response endpoint. Mode bands
run bottom to top in preparation order, with the response band at the top.

`length_regressions.csv` retains six component metrics;
`cumulative_length_regressions.csv` contains six cumulative metrics. OLS uses an
intercept and more than three valid points. Cumulative fits require every prior
preparation stage to settle and total-time fits require successful responses.
They fit their own valid subsets, rather than adding separate component slopes.

The heatmap shows component and cumulative fits on one dimensionless symmetric
color scale: `normalized_slope = slope * characteristic_length / characteristic_response`.
Characteristic length is the midpoint of each fit's valid length range;
characteristic response is its OLS prediction there. CSVs retain raw slopes and
both reference values. Unavailable fits or reference responses of absolute value
<=1e-12 appear gray. Fits do not estimate censored values or establish a scaling
law. Runtime imports have no dependencies on reports or historical archives.
