# Seven models, four assays: accepted v8 baseline

This is the active pipeline for the cross-model length comparison. It preserves
the user's accepted revised revision-8 protocols and analysis. The suite has
seven models, so there are seven model files rather than six.

- `models/{wave_pinning,otsuji,goryachev,legi,holmes_model3,debelly,spring}.py`:
  equations, selected parameters, native initialization and mechanical solvers.
- `runner.py`: shared cues, full-state preparation, eligibility, all four assays,
  Fourier reliability scoring, SLURM submission, collection and comparison.
- `plot.py`: mode shading, OLS lines, combined settling/response figure and slope
  heatmap. Each has PNG, PDF and SVG exports.
- `protocol.json`: frozen task list, model-specific duration grids and shared
  timing rules. Cue definitions are explicit in `runner.Cue`; model constants
  are in the corresponding model file.
- `baseline/`: compact accepted v8 records, numeric reference tables and figures.
- `results/`: freshly computed records, figures, tables and verification.
- `work/`: temporary SLURM workers, full-state checkpoints, diagnostic traces and
  logs. Completed work is moved to the archive after validation.

The combined figure shows cumulative time curves in one facet for each
competition, gradient-reversal and finite-countercue assay: each settled
preparation stage adds to the preceding stages, and the final curve adds the
successful response endpoint. Finite countercue includes cue-on then no-cue
settling; Spring gradient reversal includes spontaneous then left-gradient
settling. For finite countercue, the response endpoint is the shortest passing
pulse duration, excluding the subsequent recovery observation window. Mode
shading bands run from bottom to top in preparation order, with the response
band at the top. Periodic outputs remain separate.

The numerical grid, generic BDF integrator and stimulus primitives reuse
`../src/polarity1d`. Runtime execution has no imports or file dependencies in
`reports`, earlier revision directories, or the historical archive.

## Run and render

Use the existing environment recorded in `environment.json`; no new environment
or installation is required. From the repository root:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/4473331/.conda/envs/scaling_sandbox/bin/python
$PY polarity1d/four_assay/runner.py launch
$PY polarity1d/four_assay/runner.py status
```

`launch` submits 441 preparation tasks on SLURM, then plans 4096 response jobs,
then collects and plots 770 assay cases. Dependencies require successful worker
completion. Jobs use QoS `small`, one CPU and 3 GB; arrays have no concurrency
throttle. It never simulates spring competition. Source hashes and checkpoint
hashes are checked by each worker. Missing execution outputs and incomplete
finite-duration sweeps prevent finalization; they cannot become biological
censoring. `launch` refuses an already initialized work directory and finalization
refuses to overwrite existing results. Preserve/move an earlier `work/` and
`results/` before a new full campaign; `baseline/` stays intact.

To redraw the current figures without running simulations:

```bash
$PY polarity1d/four_assay/plot.py
$PY polarity1d/four_assay/runner.py compare
```

The plotter accepts `--results PATH` for another complete compact dataset. Both
rendering and comparison work after worker scratch is archived. Full records
are in one `records.json.gz`: 770 raw assay cases, 441 preparations, every tested
finite duration, every full-state settling check, periodic fit diagnostics and
frozen provenance. Raw Goryachev successes remain in these records; analysis
classification is exported in `endpoints.csv` with original status and endpoint.

## Protocol retained

All models use a 1000-native-time-unit full-state confirmation window, a 10000
cap per preparation phase, checks every 100 units, and relative tolerance 0.001
for pointwise field range and final drift projected over the window. Each field
uses its window maximum with floor 1e-8; mechanical length/displacement use the
reference length. Preparation saves every 5 units. deBelly drift uses the last
recorded finite difference; RD and spring use their full RHS. Autonomous
countercue recovery lasts 10000 units with a final 500-unit response confirmation.
Competition and reversal observations last 20000 units. Response sampling is
2.5 units, pulse sampling min(2.5, duration/20), periodic sampling P/80. Periodic
traces cover six periods, with separate P=100/P=400 scores and plots. Finite
thresholds are the shortest passing *tested* duration, using the complete original
model-specific dense grids.

| Model group | Competition | Reversal gradient | Finite training / countercue |
|---|---|---|---|
| Wave-pinning, Otsuji, Goryachev, LEGI | 0.2 half-cosine patches, outer 10%, 2% imbalance; 20-unit withdrawal taper | 0.1 linear, 20-unit crossfade | Left 0.1 linear then no cue; right 0.1 linear pulse |
| Holmes 3 | 0.2 half-cosine patches, outer 20%, 2% imbalance; abrupt withdrawal | 1.0 linear, 20-unit crossfade | Left 0.2 Gaussian, center 0.15L, width 0.1L, then no cue; right 1.0 linear pulse |
| deBelly | 3.6 half-cosine patches, outer 20%, 2% imbalance; abrupt withdrawal | 3.6 linear, 20-unit crossfade | Matching left/right 3.6 half-cosine patches, outer 20%; intervening no cue |
| Spring | **Paused** | Spontaneous left preparation, then 4.0 linear and 20-unit crossfade | Matching left/right 4.0 Gaussians, width parameter 0.2L; intervening no cue |

Cue application before nonperiodic assays continues until the full-state
criterion passes. Finite countercue also requires a qualified autonomous left
state after cue withdrawal. deBelly additionally checks Rac and bound-MCA x50
ranges <=0.002 over the last 1000 units. Holmes uses the corrected original
Eq. 5/6 factor of two with selected basal Rac activation 3.5. Spring retains
the author's nonlinear cue activation; physical domain length is dynamic.

Remaining differences are intentional baseline choices, rather than new
standardizations: amplitudes, cue shapes, withdrawal taper, native initial
states/units, extra mechanical/species checks, solver families, and periodic
preparation. Periodic starts remain homogeneous in the four Jilkine models,
symmetric algebraic equilibrium in deBelly, 1500 units of spontaneous preparation
and reflection in spring, and a 40-unit local cue followed by 2000 no-cue units
in Holmes. Spring gradient preparation also retains spontaneous polarity and
reflection. Stimulus strengths have **not** been equalized.

## Analysis and interpretation

Nineteen Goryachev competition cases are classified as censored with the distinct
`resolved_at_withdrawal` mode: 18 already meet the resolution rule at withdrawal;
one resolves within the first 2.5-unit sample. This is not ordinary right
censoring or failure to resolve. Colors also distinguish unobserved endpoints,
ineligible states, unsettled preparations, failed periodic reliability and pause.
Settled preparation times remain valid even if later assay eligibility fails.

OLS fits include an intercept and require more than three valid points. Each
period and periodic output is fitted separately. `length_regressions.csv`
retains the original component fits and adds normalization fields;
`cumulative_length_regressions.csv` contains direct fits to each cumulative
curve. Cumulative fits require every preceding preparation stage to settle and,
for total-time fits, a successful response endpoint. They use their own valid
subsets rather than sums of separately fitted component slopes.

The heatmap shows the 12 component metrics and eight cumulative metrics in two
groups, using one symmetric color scale for dimensionless normalized slopes:
`normalized_slope = slope * characteristic_length / characteristic_response`.
The characteristic length is the midpoint of each fit's valid length range;
the characteristic response is its OLS prediction at that length. Both values
and the raw slope are retained in the CSVs. The response denominator retains
its sign for signed periodic phase outputs. Unavailable fits and effectively
zero reference responses (absolute value <=1e-12) appear gray. The raw component
fits remain comparable to the accepted baseline; normalization and cumulative
fits are new analysis. Fits describe valid subsets and do not estimate censored
values or establish a mechanistic scaling law.

Verification details and the completed migration inventory are in
`results/comparison_to_baseline.json` and `MIGRATION.md`. The legacy revision tree
is retained under `../archive/four_test_suite_legacy/` solely for historical
provenance. The approved assay explanation v2 remains preserved there.
