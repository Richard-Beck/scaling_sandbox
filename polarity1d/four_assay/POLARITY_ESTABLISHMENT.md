# Polarity establishment from native rest

**Chosen basis for continuation:** the corrected basal-withdrawal pilot is the
selected basis for upcoming reversal experiments. See its
[decision note](pilots/polarity_establishment_ref15_relative20_basal_20261003/README.md)
for the preparation rationale and planned exclusion of pure LEGI, which cannot
provide persistent polarity after withdrawal. The other six models remain in
scope, using qualifying polarized preparations.

This companion campaign is separate from the accepted baseline, revision-11
preparation chain and 2026-10-02 finite-countercue width pilot. It tests how a
finite localized cue establishes persistent polarity, with the same spatial cue
rule available for subsequent opposite-end countercues in every model.

## Inputs

For each existing model/length, start independently at rest. Apply a left cosine
patch for duration D, then withdraw only the added patch, restoring basal
conditions. D=0 supplies a no-cue control. The six non-LEGI models retain their
original native starts and zero *extra* input after withdrawal.

LEGI starts at its exact homogeneous equilibrium under uniform positive S=0.1:
activator=0.1, inhibitor=0.1, response=0.5. During the pulse, S=0.1+patch; after
withdrawal, S=0.1. This implements the withdrawal rule in section 1.1 of the
[Mori/Jilkine/Edelstein-Keshet supplement](https://personal.math.ubc.ca/~keshet/pubs/Supplement.pdf).
That section specifies positive c0 without its numeric value; 0.1 is an explicit
campaign convention matching nominal LEGI cue amplitude, not a claimed published
value. It does not vary with L, A, width or duration. Native reaction parameters
and solvers are retained. LEGI therefore loses spatial polarity after withdrawal.
Setting total S=0 instead lets activator/inhibitor vanish and freezes the response;
the previous LEGI persistent endpoints were artifacts of that degenerate input.

The corrected report is
[pilots/polarity_establishment_ref15_relative20_basal_20261003/results/index.html](pilots/polarity_establishment_ref15_relative20_basal_20261003/results/index.html).
Its six unaffected models reuse identical, already-audited data; LEGI's 462
combinations are rerun. Source/case compatibility, source record hashes and final
checkpoint/arrival/selection checks preserve traceability. The preceding report
is retained as a superseded zero-total-input campaign. No persistent polarized
endpoint means both establishment timing cells are unavailable, rather than zero.

The seven models retain their native parameters, grids and solvers. The full
current campaign contains 3,234 model × length × amplitude × condition combinations:

- Two fixed widths: 10% and 20% of a **standardized reference length of 15** for
  all models except Spring, giving widths 1.5 and 3. Spring uses native source
  default reference 1, giving widths 0.1 and 0.2. This is a comparison convention;
  Spring's native length is not calibrated in microns. Widths are held fixed
  in native physical cue coordinates across lengths.
- A third condition uses **20% of cell length**, with no reference-length scaling.
  For Spring it follows the instantaneous physical domain during the pulse;
  for the other models it follows the native cue domain L.
- Seven peak amplitudes: 0.1, 0.25, 0.5, 1, 2, 4 and 10 times each model's current
  amplitude. Peak amplitude is independent of length in every condition.
- A fixed patch wider than a short domain is evaluated only on the available
  domain; neither its width nor peak is changed. This produces a less localized
  stimulus on short cells and is retained as part of the physical-width test.

The right patch is the exact spatial mirror of the left, with the same width and
peak. No reversal experiment is run by this establishment campaign. Spring maps
the patch onto its instantaneous physical domain; deBelly retains the native
input coordinates and sampling of the established solver.

## Duration scan and bounded refinement

Each combination starts with a no-cue control and 17 logarithmically spaced
positive durations. Initial ranges are 0.05–5,120 for the chemical models and
Holmes 3, 0.5–2,000 for deBelly, and 0.002–50 for Spring, in native time units.

The objective is the strongest **cue-aligned** final polarity, using the
existing orientation score with sign reversed for a left cue. Verified states
with contrast below 0.05 cannot set the polarized target: a high normalized
orientation of a vanishing signal should not defeat a stronger observable
response. A polarized candidate also needs cue-aligned strength at least 0.05.
Opposite-directed states are retained in records and endpoint clustering.

If the largest verified tested duration still improves strength by more than
1%, allow at most two fourfold upper extensions, up to 60,000. Two bounded
refinement rounds sample geometric midpoints around the three strongest
durations. Once the strongest observed endpoint is identified, refine each
sampled transition into its 1% strength band for up to eight steps, aiming for a
1% duration bracket. No global monotonicity assumption or bisection across the
entire duration range is made. Disconnected success intervals are permitted.

Select the shortest tested duration yielding at least 99% of the
strongest observable verified cue-aligned polarity. Upper-range increases,
unsettled releases and insufficiently refined duration boundaries remain
explicit flags. This is a bounded search for the strongest observed endpoint,
not proof of a global optimum. If the no-cue control already reaches that
cue-aligned endpoint, select duration zero and mark that no polarizing cue was needed.

## Release convergence and timing

Keep the existing common all-field convergence/holdout criteria: 1,000-unit
chunks, two consecutive passing windows at 0.0001 fixed-scale variation and
projected drift, followed by a 10,000-unit holdout allowing 0.001 excursion and
a passing final window. Maximum release continuation is 60,000. Sample every
5 units, with an initial 10-unit fine segment sampled every 0.25 units (0.05
for Spring). This resolves fast post-withdrawal changes without retaining an
entire long continuation at the fine interval.

Fixed field scales are the per-field maximum of the frozen revision-11
conditioning and release scales for each model/length. They remain fixed across
all amplitudes, widths and durations. deBelly displacement is centered and
physical length reconstructed, as in the existing settling rules.

After a reference is verified, polarity arrival is the first sampled entry into
1% of its final signed polarity that persists through the remaining history,
with at least 1,000 subsequent units. An absolute tolerance floor of 0.000001
handles near-zero final polarity. Also retain the separate all-field 1% arrival
for audit, but the heatmap uses the requested polarity clock. Total establishment
time is selected cue duration plus post-release polarity arrival, excluding
the extra confirmation and holdout time. A zero post-release arrival means the
release started within tolerance, not exact equilibrium.

## Final heatmaps grouped by model

Each model has one figure with four metric rows and three condition columns
(fixed width 10% of reference, fixed width 20% of reference, relative width 20% of L),
with L on the horizontal axis and actual native peak amplitude A on the vertical axis:

1. Selected cue duration.
2. Total establishment time.
3. Selected final cue-aligned polarity (or strongest observed strength when no
   qualifying polarized trial is available, explicitly marked).
4. Apparent number of endpoints.

The count is a diagnostic over **verified sampled complete native profiles**,
including the zero-duration control. Greedy complete-link clustering uses a
1% maximum pointwise distance on the same frozen field scales. It distinguishes
profile differences hidden by the scalar polarity score; continuous families
can generate many clusters. It is neither an exact attractor count nor
an additional search objective. No extra duration refinements chase the count.

Unpolarized, unsettled and no-qualifying-polarized cases have distinct symbols.
Flags expose unresolved upper strength, some unsettled releases, and duration
brackets wider than 1%. Colors share a scale across all three conditions within a model; native
time units differ between models. PNG, PDF and SVG are written **once at complete
finalization**, with one HTML collection organized into seven model sections.
Each model has PNG/PDF/SVG downloads, and a combined seven-page PDF contains the
whole collection. No intermediate heatmaps are produced. The preceding
model-specific-reference campaign and its two-column report remain preserved
under `pilots/polarity_establishment_20261003/`.

Every combination has the same symbols in all four heatmap sets. An open circle
means the gap between the selected passing cue duration and its nearest shorter
tested duration is greater than 1% of the selected duration. This measures how
finely the duration threshold was located; it is separate from the 1% polarity
strength/arrival criteria. The circle does not imply that the selected endpoint
failed to settle. Concurrent uncertainty flags are offset within a cell to avoid
hiding one another. Diamonds identify zero-cue selections throughout the collection.
Presentation updates retain the endpoint table and record rendering provenance
and figure hashes in `results/rendering_provenance.json`.

## Execution and retained evidence

Use the existing `scaling_sandbox` Conda environment. The launcher freezes and
hashes scientific sources and its manifest, submits one CPU/3 GB workers under
SLURM QoS `small`, splits arrays at 900 tasks, and uses no concurrency throttle.
Every trial retains a native restart checkpoint, compact sampled spatial
profiles, every release sample's polarity/contrast/full-field reference distance,
and all convergence/holdout diagnostics. Full spatial histories remain in memory
while verifying arrival; the stored distance history records the all-field audit.
The finalizer verifies checkpoints and recalculates every polarity arrival.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/4473331/.conda/envs/scaling_sandbox/bin/python
$PY polarity1d/four_assay/polarity_establishment.py launch
$PY polarity1d/four_assay/polarity_establishment.py status
```

Default destination: `pilots/polarity_establishment_ref15_relative20_20261003/`. Custom `--out`
destinations are supported. A SLURM monitor updates the status page and endpoint
table while workers run. The page has a clear generation timestamp and producer
state, with no browser self-refresh. `resume` submits only incomplete cases and
reuses completed per-duration trials. Missing execution is an error, not biological
censoring; final heatmaps require complete case coverage. `--reference-length`
(default 15) and `--spring-reference-length` (default 1) are configurable at
initialization/launch; workers use the frozen manifest values.

Validation includes algorithm checks in `tests/test_polarity_establishment.py`,
all-grid cue-mirror/cutoff checks, native release smoke tests for all seven models
and all three conditions, and a reduced SLURM integration campaign. Development smoke and
integration outputs are labeled separately from the full campaign.

## Corrected basal-withdrawal results

The basal-withdrawal campaign completed 3,234 combinations without execution
errors (2,476 polarized, 758 unpolarized). All 462 LEGI combinations yield one
homogeneous endpoint and no persistent polarity. All 8,316 LEGI trial profiles
return to A=I=0.1 and R=0.5; maximum final polarity score is 1.3e-16.
The 2,772 non-LEGI summary rows exactly match the preceding campaign.
Ten tests pass, and final integrity, arrival and selection audits cover all
105,947 combined trials. See `results/withdrawal_correction_validation.json`
for correction-specific comparisons and `results/validation.json` for the
full campaign audit. Downloaded primary references and their hashes are retained
under `reference_sources/`.
