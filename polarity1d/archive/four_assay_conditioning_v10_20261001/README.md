# Seven models, two assays: normalized gradient cues

The active suite runs gradient reversal and finite countercue at every existing
sampled length. User-authorized revision 10 aligns preparation across all seven models, following
the normalized-cue revision-9 rollout.
Competition and periodic tracking are deferred and absent from the active task
list, runner, figures and slope tables. The directory name remains `four_assay`
so established paths to the active runner and figures continue to work.

- `models/`: the seven unchanged model implementations.
- `runner.py`: linear cues, native preparation, endpoint gates, SLURM execution,
  collection and paired baseline validation.
- `plot.py`: separate conditioning/release/response facets and normalized slopes;
  PNG, PDF and SVG exports.
- `protocol.json`: 308 cases (154 model–length combinations, two assays), original
  dense duration grids, shared cue amplitudes and common preparation order.
  Revision 10 reruns the 42 Spring cases and reuses the other 266 verified cases.
- `reference/normalized_v9/`: compact immutable input for the Spring amendment;
  prior source/results/work are in `../archive/four_assay_normalized_v9_20261001/`.
- `baseline/`: immutable accepted v8 records and reference artifacts.
- `results/`: current normalized two-assay outputs and paired validation.
- `work/`: current SLURM workers, checkpoints, diagnostic traces and logs.
- `reference/continuation_20261001/`: verified continuation profiles for both
  settling events at all 154 model–length combinations, full preparation
  histories, reference-arrival times and independent holdout checks. See its
  [assessment](reference/continuation_20261001/RESULTS.md). This analysis leaves
  the active assay inputs, clocks, figures and fits unchanged.

The prior active results, work metadata and source snapshot are preserved in
`../archive/four_assay_v8_pre_normalization_20261001/`. The successful 58-case
normalization pilot is retained in `pilots/linear_shared_20261001/`; its frozen
source references describe the preceding v8 code and results, now archived.
`MIGRATION.md` records the earlier migration/reproduction, not this new protocol.

The preceding full revision-9 sweep passed paired validation: 308 preparations,
308 assay cases and 34092 finite-duration tests. All 213 baseline-successful
readouts were preserved, with no losses or gains and unchanged classifications
for 70 censored, 19 ineligible and six unsettled cases. The unchanged controls
have zero mismatches, including 23280 pulse tests. Its component and cumulative fits were independently checked against raw
records. The current figure revision removes cumulative fits and deduplicates
shared conditioning.
Holmes pulse thresholds are unchanged; deBelly thresholds are 7–51% shorter and
Spring thresholds 8–42% shorter across their successful sampled lengths.

The revision-10 Spring amendment is complete: 42 fresh Spring cases and 266
unchanged reference cases. All 42 Spring readouts remain successful and their
measured response endpoints are unchanged. Both Spring assays condition for
1400 units at every length; the removed pre-settling formerly made gradient
preparation total 2400–2600 units. Across all 154 model–length combinations,
conditioning times and full-state settling checks now match between assays.
All 213 successful suite readouts remain, with zero control mismatches. All
26 eligible individual-phase fits were independently verified.

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

Every model starts from its native initial state and conditions directly under
the left cue until settled. There is no spontaneous pre-settling or reflection,
including for Spring. Gradient reversal then uses the existing 20-unit
left-to-right crossfade and 20000-unit observation window. Finite countercue prepares
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

For this Spring-only revision, the pipeline submits 42 preparations, plans
qualified response jobs in 10-duration chunks, and combines them with the 266
unchanged cases from the frozen revision-9 reference before rendering. Reference
reuse is specific to this amendment; changes to the other models require their
records to be rerun. The full case grid remains 308. All jobs use
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
without compensation by gains. Against v8, the six non-Spring gradient experiments and the four models whose
finite cues already matched are exact controls. Against normalized v9, every
case except Spring gradient reversal is an exact control, including full pulse
tests, settling checks and eligibility diagnostics. Shared conditioning times
are checked between assays at every model and length. The authorized removal of Spring gradient pre-settling/reflection may change
its numerical conditioning/response times. Earlier finite-cue normalization
already changed preparation/response times in Holmes 3, deBelly and Spring.
`comparison_to_baseline.csv` and `.json` distinguish those changes from failures.

## Figures and fits

The revised [assay explanation PNG](results/assay_explanation.png)
([PDF](results/assay_explanation.pdf), [SVG](results/assay_explanation.svg))
uses the approved keeper’s visual style for the four current preparation/response
phases. [Caption and reconstruction notes](results/ASSAY_EXPLANATION_NOTES.md)
explain the verified wave-pinning examples and rebuilding. This new revision
awaits user review; the archived approved keeper is preserved.

`length_response_and_settling` has four columns and one row per model: shared
conditioning cue time, countercue release/no-cue settling, gradient reversal
response time, and countercue pulse duration. Conditioning is shown once after
checking that both assays have identical conditioning times. Response curves
exclude conditioning and release settling; countercue duration also excludes
post-pulse recovery observation. All plotted times are individual phase times.

`length_regressions.csv` contains four metrics, with OLS intercepts and more than
three valid points. There are no cumulative fits or cumulative regression CSV
in active results. The heatmap shows only these four normalized slopes on one
dimensionless symmetric color scale:
`normalized_slope = slope * characteristic_length / characteristic_response`.
Characteristic length is the midpoint of each fit's valid length range;
characteristic response is its OLS prediction there. CSVs retain raw slopes and
both reference values. Unavailable fits or reference responses of absolute value
<=1e-12 appear gray. Fits do not estimate censored values or establish a scaling
law. Runtime imports have no dependencies on reports or historical archives.
