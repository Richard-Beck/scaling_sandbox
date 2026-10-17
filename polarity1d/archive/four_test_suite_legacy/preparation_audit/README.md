# Initial-cue steady-state audit

**Completed 2026-09-29:** 308 cases across all 14 combinations; 117 steady,
173 not steady and flagged for later, 18 unverified under the spring pause.
SLURM array 20642205 completed all 287 authorized preparation reconstructions
with exit code zero. No endpoint campaigns or new spring competition runs occurred.

User-requested diagnostic of the original accepted competition and gradient-reversal
preparations across seven models and all existing lengths (308 cases / 14 model–test
combinations). The user confirmed that preparations **not at steady state** should
be flagged for later reruns. Flags are not authorization to launch those reruns.

## Scope and paused spring workflow

Only original preparation phases are reconstructed where full archived fields are
missing. No preparation is prolonged, no stimulus is corrected, and no competition
or reversal response campaign is rerun. Parameters, mesh, seed and initial states
match the accepted adapters. Original executable hashes are verified where recorded.
SLURM handles the 287 preparation reconstructions (QoS small, 1 CPU / 2 GB per task,
no concurrency throttle). Result tables include the original eligibility/status
separately from steady-state assessment.

**Spring competition remains PAUSED — NO CLEAR LENGTH CONTRAST.** That combination
is read-only: no additional competition simulations, steady-cue runs, tolerance
studies or rescue attempts are launched. Original 40-unit preparation fields exist
at L0=0.6, 1.8 and 2.6 and are assessed directly. The other 18 original preparation
states are explicitly not verified, rather than invented or reconstructed despite
the pause. The rejected steady-cue candidate is not substituted for the original
protocol. Any spring competition rerun flag remains blocked by the user's pause.

## Assessment time

| Model group | Competition assessment | Reversal assessment |
|---|---|---|
| Wave-pinning, Otsuji, Goryachev, LEGI | Time 20, immediately before the original 20–40 cue taper starts | Time 2000, before gradient cross-fade |
| Holmes Model 3 | Time 40, immediately before cue removal | Time 2000, before gradient cross-fade |
| deBelly | Time 80, immediately before cue removal | Time 2000, before gradient cross-fade |
| Spring (original protocol) | Archived time 40 states only | End of the original 2000-unit left-gradient phase, following spontaneous preparation |

Assessing the original four competition models at time 40 would confuse response
to withdrawal with steadiness under the initial cue. Their cue is already changing
from time 20 onward. The audit checks the end of the constant-cue plateau instead.

## Operational stationarity criteria

A state is labelled steady if **every** checked field passes both:

- Maximum pointwise range over the final constant-cue window is at most 0.1% of
  that field's maximum absolute value in the window.
- Its final instantaneous rate, projected over the same window, is at most 0.1%
  of that scale. deBelly uses the final stored finite-difference rate because its
  specialized staggered integrator does not expose one continuous combined RHS.

Competition windows are the final 10 native time units before withdrawal/taper.
Reversal windows are the final 500 units. Scale has a floor of 1e-8 in native field
units to avoid division by zero. All metrics and scales are retained for review;
the binary label is an operational tolerance check, not an assertion of exact
mathematical equilibrium. Different model-native units remain explicit.

- Reaction–diffusion models: **all state species**, using full reaction plus
  diffusion RHS for the final rate. Domain length is fixed by construction.
- Spring: active and inactive species **and dynamic physical length**, with full
  model RHS and exact cue transfer. Length changes are normalized by L0.
- deBelly: Rac, Rho, bound/unbound MCA, displacement shape, velocity, tension and
  **physical length = reference length + boundary separation**. Common translation
  is removed from displacement; deformation and length are retained. Displacement
  changes and length changes are normalized by reference length. Other fields use
  their own scales. Derived activation fields are not independent state variables.

A stationary but unpolarized/ineligible preparation may pass this stationarity
test. Original eligibility remains visible; it is not conflated with steady state.
Passing the earlier orientation/contrast gate alone does not prove all internal
fields and mechanics were steady.

## Key findings

- Every assessed competition preparation fails the stationarity criteria except
  Otsuji at L=1, which is stationary but was already ineligible for two-peak
  competition. The original four models are assessed before taper onset at time 20.
- Wave-pinning, LEGI, Holmes and spring reversal preparations pass at every length.
- Otsuji reversal is flagged at L=25, 27, 29, 31, 33, 35, 37 and 39.
- Goryachev reversal is flagged at L=9, 11, 12, 13, 15, 17, 18, 19, 21, 23, 25,
  27, 29, 31, 33, 35, 37 and 39.
- deBelly reversal is flagged at L=18, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37
  and 39. Bound/unbound MCA fail at all 12; tension also fails at nine. Rac and
  Rho pass these particular full-window checks. Orientation or Rac-only steadiness
  therefore misses the remaining internal preparation transients.
- The three archived original spring competition preparations fail; its other
  18 original lengths remain unverified. All spring competition flags are paused,
  not an invitation to retry the workflow.

## Outputs

- `results/REPORT.md`: 14-combination summary.
- `results/combo_summary.csv`: counts by model/test.
- `results/case_diagnostics.csv`: all lengths, assessments and rerun flags.
- `results/flagged_for_later_rerun.csv`: preparations that fail stationarity.
- `results/audit_summary.json`: completeness and scope checks.
- `data/task_*.json`: field-specific changes, rates, scales and provenance.
- `tails/task_*.npz`: retained preparation-window fields.
- `manifest.json`, `frozen_provenance.json`: cases and frozen execution sources.

Scientific flags do not change the accepted endpoint tables or launch any future
campaign. Review and authorization are separate from this audit. Preserve the
approved keeper `../publication_figures/assay_explanation_v2.*`.
