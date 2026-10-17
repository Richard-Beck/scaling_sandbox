# Revision 7: standardized modal assay protocols

Authorized 2026-10-01: remove undocumented preparation/cue/timing differences,
use modal choices and the longest existing observation/confirmation times.
Revision 6 is preserved. This campaign changes protocols, not model coefficients.
Spring competition remains paused under the standing explicit restriction.

## Shared rules

- Same 201-cell mesh and revision-6 lengths, periods 100 and 400 in native time.
- Cue geometry uses relative position q=x/L. All linear cues have amplitude 0.1:
  left 0.1*(1-q), right 0.1*q. No model-specific cue strength tuning.
- Competition: half-cosine cues in the outer 10% of each end, nominal amplitude
  0.2, left/right multiplier 1.01/0.99, held until steady. Shared 20-unit
  half-cosine withdrawal taper, then 20,000-unit autonomous observation.
- Gradient reversal: directly condition the model's initial state under the left
  linear gradient until steady; crossfade to right over 20 units; observe 20,000.
  Spring has no extra spontaneous preconditioning or orientation reflection.
- Finite countercue: same left-gradient conditioning, then no-cue settling.
  Challenge with the right linear gradient for each independently tested duration,
  withdraw abruptly, and recover without cue for 10,000 units for every model.
  Shared duration grid is the sorted union of all revision-5 grids (0.02–5120).
  Ineligible autonomous states are recorded, not rescued through special training.
- Periodic: directly settle under the left gradient, then six cycles of
  0.05*(1-cos(2*pi*t/P)*(2*q-1)). Score P through 6P, 80 samples per cycle,
  same edge-orientation harmonic reliability thresholds for every model.
- Each preparation phase has a 10,000-unit cap, 1,000-unit confirmation window,
  100-unit checks, <=0.1% maximum pointwise field range and projected drift.
  Every evolving species is checked. Physical length and displacement are checked
  for mechanical models (with rigid translation removed as a coordinate gauge).
  deBelly retains Rac/MCA half-mass stability checks over the same window.
- Shared polarity eligibility: left orientation below -0.05, orientation range
  below 0.005 over 1,000 units, absolute active contrast at least 0.05.
  Competition needs two peak excesses above 0.1 before/after taper.
- Shared response definitions retained: competition loser/winner <0.1 with winner
  >0.1 for 50 units; gradient right orientation/contrast >=0.05 for 20 units;
  finite recovery final 500-unit orientation range <0.005, contrast >=0.05 and
  right orientation >0.05. Periodic preparation must meet the same left eligibility.
- Ineligible starts, unsettled preparations, censored responses, unreliable fits,
  pauses and execution errors remain separate records. Gray in PNGs groups
  unusable outcomes; read CSV for their distinction. Missing finite chunks never
  produce a threshold. No endpoints are reused from previous revisions.

## Remaining differences and allowed state-specific checks

- Native time, concentration, and length units are not physically calibrated
  across models. Spring x-axis is resting L0; cue geometry follows its evolving
  physical length. The other models use their existing domains.
- Cue coupling follows each model's equations: spring alters basal activation in
  its nonlinear reaction; other models apply their existing activation terms.
  Equal numerical amplitude does not establish equal physical stimulus strength.
- Initial states and model coefficients remain established defaults/selected
  model regimes (including Holmes basal Rac activation 3.5, no PI layer).
  No special polarity training is added. Competition retains small seed-1 initial
  noise for RD/Holmes models; deBelly uses its solver's algebraic initial state.
- Settling checks include all species and mechanical variables rather than forcing
  a two-species state description. deBelly drift uses a saved finite difference;
  others use their full RHS. Spring cue drift now matches its nonlinear solver
  activation term; revision 6 used cue*v instead.
- Solvers remain model-specific: sparse BDF for RD/Holmes, LSODA for spring,
  specialized mechanics solver for deBelly. These are numerical implementations,
  not differences in stimulus histories. Gradient spring output is denser during
  the first 1,000 units; the same endpoint and observation horizon apply.
- Spring competition is the sole omitted assay, retained under the explicit pause.

## Execution

Existing scaling_sandbox Conda environment, one core and 3 GB per SLURM array slot,
QoS small, no array concurrency throttle. `launch.json` identifies jobs and stage.
Sources, model code, protocol notes, manifests and duration grids are frozen;
workers verify hashes. Do not launch twice or modify frozen files in flight.
The report producer regenerates two timestamped PNGs as records arrive and at
completion; HTML has no timed refresh. `results/status.json` records completeness.
Scientific interpretation remains subject to review even after numerical completion.
