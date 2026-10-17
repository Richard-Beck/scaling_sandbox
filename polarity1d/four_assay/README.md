> **Output retirement authorized (2026-10-17):** the user requested removal of large outputs from this deprecated tree. The [retirement record](OUTPUT_RETIREMENT_20261017.md) describes retained sources, scalar records and selected inputs, and loss of unselected native-array coverage. Historical work paths below must not be treated as complete local raw data. Current manuscript evidence remains under `../manuscript/`.

> **deBelly scientific warning (2026-10-08):** deBelly simulations in this historical suite using `center_boundary_motion=True` use a flawed, deprecated centered endpoint extension. Do not use them to infer length responses or populate current figures. They artificially couple endpoint commands and change internal strain/tension. Use the corrected independent-endpoint manuscript rerun instead; the other models are unaffected.

> **Deprecated as an active workflow (2026-10-03).** The user requested a clean five-model rerun from corrected finite-patch/basal establishment endpoints, excluding Otsuji and LEGI. Use [../manuscript/](../manuscript/README.md). This directory now supplies frozen provenance; the material below records the preceding workflow. See [DEPRECATED.md](DEPRECATED.md).

# Seven models, two assays: settled preparation chain

The [corrected basal-withdrawal pilot decision](pilots/polarity_establishment_ref15_relative20_basal_20261003/README.md)
marks that pilot as the chosen basis for subsequent reversal experiments.
We plan to exclude pure LEGI from that continuation because it loses polarity
under proper basal withdrawal and cannot supply a cue-free polarized starting
state. Its pilot results remain available as evidence for this decision.

The new [polarity-establishment companion](POLARITY_ESTABLISHMENT.md) scans
rest → finite patch → verified basal release, refines the strongest
observed persistent endpoint, and produces one final collection grouped by model
with four metric rows and three condition columns, including an approximate
endpoint-count diagnostic. The current rerun uses standardized reference 15
(Spring native reference 1), two fixed widths and a 20%-of-length condition.
LEGI restores positive uniform background after patch withdrawal and starts at its
homogeneous basal equilibrium; the previous zero-total-input LEGI endpoints are
superseded. Its separate output directory is `pilots/polarity_establishment_ref15_relative20_basal_20261003/`;
the preceding two-column establishment campaign is preserved under
`pilots/polarity_establishment_20261003/`. The existing results below remain the
settled preparation-chain comparison.

Revision 11 runs all 308 gradient-reversal and finite-countercue cases at the
154 existing model–length combinations. Competition and periodic tracking remain
deferred. The directory name stays `four_assay` for established figure paths.

## Preparation dependencies and clocks

| Timing test | Input | Reported time |
|---|---|---|
| Shared conditioning | Native resting initial state | First persistent 1% arrival at the conditioning reference |
| Cue-free release | Fully settled conditioning | First persistent 1% arrival at the release reference |
| Gradient reversal | Fully settled conditioning | Existing sustained right-polarity response from crossfade onset |
| Finite countercue | Fully settled release | Existing shortest passing tested pulse duration |

Conditioning is established once per model and length. Its exact saved full-state
checkpoint starts both release and gradient reversal. Release's exact saved
checkpoint starts every finite pulse. No pre-settling or reflection is introduced.
The first two reported clocks exclude confirmation/verification time, while
response experiments always start from the verified final profiles.

`settling.py` supplies one continuation and measurement ruleset for both settling
events and all models. Continue in 1000-unit chunks, sample every 5 units, require
two consecutive windows with each field's pointwise temporal range and final
sampled rate projected over 1000 units at most 0.0001 of a fixed scale. Then check
an additional 10000-unit holdout: excursions at most 0.001 and a passing final
window. The final holdout state is the reference. If a holdout fails, continuation
resumes; the total per-event cap is 60000. Unverified states remain unsettled.
Arrival is the earliest sampled entry after which all points of all monitored
fields stay within 0.01 of the reference on those same scales, through at least
1000 subsequent units. All histories are retained and audited independently.

Fixed scales are frozen in `reference/continuation_20261001/fixed_scales.json`,
calibrated from the previous native preparation transients/endpoints. Chemical,
velocity and tension fields use their maximum absolute transient amplitude with
floor 1e-8; mechanical displacement/length use nominal length. deBelly displacement
is mean-centered to remove rigid translation, with physical length reconstructed
from boundary separation; Spring retains dynamic physical length. This state
normalization is distinct from the length-slope normalization below.

Linear cues still match within each model: left = A*(1-x/L), right = A*x/L.
A is 0.1 for wave-pinning/Otsuji/Goryachev/LEGI, 1 for Holmes 3, 3.6 for deBelly,
and 4 for Spring. Native units, meshes, solvers and physical parameters are
unchanged. Eligibility remains stable left orientation below -0.05, sufficient
contrast, plus existing deBelly x50 checks. A converged unpolarized profile is
valid for measuring settling but cannot start an eligible reversal assay.
Gradient reversal retains the 20-unit crossfade, 20000-unit horizon and 20-unit
right confirmation. Countercue retains every original dense pulse duration,
10000-unit recovery and final 500-unit gate. No recovery time enters its readout.

## Files, execution and validation

- `runner.py`, `settling.py`: shared chain, response assays, SLURM, collection and audit.
- `protocol.json`: revision 11, original 308 tasks and duration grids, continuation rules.
- `plot.py`: four separate timing columns and normalized-slope heatmap.
- `work/checkpoints/`, `work/histories/`: full native states and complete phase histories.
- `reference/chained_v11_20261001/`: persistent final profile copies, preparation records,
  hashes and chain manifest, published at finalization.
- `results/`: current figures, complete records, clock/endpoint/fit tables and validation.
- `reference/normalized_v10/`: immutable compact previous comparison input.
- `reference/continuation_20261001/`: preceding exploratory continuation assessment;
  its release profiles originated from old conditioning, so revision 11 replaces
  them with the explicitly chained profiles.
- `baseline/`: immutable accepted v8 reference.

Use the existing environment:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/4473331/.conda/envs/scaling_sandbox/bin/python
$PY polarity1d/four_assay/runner.py launch
$PY polarity1d/four_assay/runner.py status
```

Preserve existing work/results before launching another campaign; custom `--work`
and `--results` are supported. The current campaign submits 154 paired preparation
workers, then response workers with ten pulse durations per chunk. All use SLURM
QoS `small`, one CPU/3 GB, without an array concurrency throttle. Sources and
inputs are frozen and checked throughout. Missing execution or duration coverage
is an error, not biological censoring. The landing page/status expose stage and
generation time, without timed browser refresh.

`comparison_to_previous.json/csv` compare against revision 10, report gains/losses
and endpoint changes, and audit full-history arrival calculations, eligibility,
checkpoint hashes and conditioning -> release dependencies. All response inputs
are intentionally refined, so numerical equality is not an exact-control criterion.
Losses of previously successful readouts are flagged explicitly. A descriptive v8
comparison remains in `comparison_to_baseline.json/csv`. Execution validity and
response-preservation outcomes are separate; neither establishes a scaling law.

## Figures and fits

`length_response_and_settling` shows shared conditioning arrival, cue-free release
arrival, gradient response and finite pulse duration. Settling fits include every
verified reference, including unpolarized states; response fits include successful
readouts only. `length_regressions.csv` retains OLS raw slopes/intercepts, valid
length ranges and R-squared, requiring more than three valid points. The heatmap
uses slope * characteristic_length / fitted_response, with characteristic length
at the midpoint of each fit's valid range. All four metrics are refitted; there
are no cumulative fits. Zero reference responses/unavailable fits remain gray.
Positive flat curves have zero normalized slope. Zero arrival means the starting
state was already inside the 1% neighborhood, not exact instantaneous equilibrium.

Revision 10's sources, work and figures are preserved in
`../archive/four_assay_conditioning_v10_20261001/`. Earlier v9 and v8 archives and
the approved `publication_figures/assay_explanation_v2` keeper remain intact.
The active explanatory figure is rebuilt from the current wave-pinning chain;
its notes and hashes identify the specific generation and records.

## Companion spatial profile figure

`results/length_response_profiles.png/pdf/svg` accompanies the four timing columns
with full native Rac/active profiles at initial, 50%, and 99% of each measured
time. Lengths and unavailable clocks are explicit in each panel; Holmes 3 and
Spring use labeled existing grid choices. See
[profile notes](results/LENGTH_RESPONSE_PROFILES_NOTES.md) for sampling, native
units, validation, and regeneration with `profile_companion.py`.
