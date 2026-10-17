# Full spring competition length sweep with steady cue preparation

**PAUSED BY USER (2026-09-29): NO CLEAR LENGTH CONTRAST.** Do not run, extend,
refine, retry, smooth or attempt to rescue this workflow without a new explicit
user request. Preserve this sweep as diagnostic evidence. The reproduction
instructions below document past work; they are not instructions to execute it.

Completed 2026-09-29 for all 21 accepted spring resting lengths: L0=0.6 to 2.6
in steps of 0.1. Three exact-protocol quick-test cases were reused; 18 additional
lengths were simulated using SLURM QoS small, one CPU and 2 GB per task, with no
array concurrency throttle. All original model and campaign files remain unchanged.

## Results

- `results/spring_competition_steady_full.png/pdf/svg`: original 40-unit cue curve
  beside the new raw endpoints, preparation times, and minimum post-release contrast.
- `results/endpoints.csv`: all 21 endpoints and key diagnostics.
- `results/audit.json`: completeness, convergence, and monotonicity checks.
- `results/all_results.json`: complete records with preparation histories.
- `data/preparation_L*.npz`, `data/released_L*.npz`: full simulated fields and readouts.

All 21 cases pass the chemical and physical-length steady-state checks, establish
both peaks at withdrawal, and reach a competition endpoint within the 500-unit
post-release observation window. Preparation takes 400–800 native time units.
Maximum relative mass error is 1.09e-10. Frozen original-source hashes match.

The full raw curve is **not nondecreasing**. Eight of 20 adjacent steps decrease.
Examples include 76 to 71 units from L0=0.9 to 1.0, 120 to 84.5 from 2.2 to 2.3,
and 155 to 110 from 2.4 to 2.5. No smoothing, interpolation, monotonic regression,
or endpoint replacement was applied. The deep minimum of the original 40-unit
preparation curve is absent, but a validated monotonic response is not established.

All cases undergo near-complete post-release homogenization: minimum active
spatial contrast ranges from 2.22e-16 to 3.61e-12, below the absolute integration
tolerance 2e-11. This makes numerical seeding of renewed polarity a substantive
concern. The earlier quick check already found that resolution time and winning
side can depend on tolerance and reflection when asymmetry is nearly erased.
The new irregularities should not be assigned biological mechanisms or used to
claim a length-scaling law without controls on numerical symmetry breaking.

## Protocol

The original 201-cell initial state, dual Gaussian cue (amplitude 4, width 0.2,
2% imbalance), reaction parameters, diffusion, mechanics and competition criterion
are retained. Integration uses rtol=2e-9 and atol=2e-11, as in the quick test.
Preparation runs in 100-unit chunks. The final 100-unit window must have maximum
range below 1e-6 at each spatial point for both active and inactive concentrations,
physical length range below 1e-7, instantaneous chemical rate below 1e-7, and
absolute length rate below 1e-9. Centre translation is a symmetry and is allowed.
Both chemistry and physical length are checked, not simply elapsed cue time.

After release the original LSODA procedure is used: 100 units sampled at 0.5,
then a restarted 400-unit segment sampled at 5. The endpoint is the first tested
time at which losing/winning peak excess remains below 0.1 for 50 units while
winning peak excess remains above 0.1. Late endpoints are quantized to the 5-unit
tail spacing. The original campaign curve in the plot uses its original, looser
tolerances; the earlier three-point report includes matching-tolerance controls.

## Execution and recovery

Array job 20642172 completed 17 new lengths. Its L0=1.3 task spent several minutes
in the second cue-preparation segment. A bounded retry with a smaller preparation
step (20642196) was also slow and was cancelled. A sparse-BDF preparation retry
(20642198) reached the same explicit steady-state checks and completed the
unchanged LSODA post-release procedure. Its endpoint is 79 units. The stalled
original task was cancelled after the replacement completed. The solver exception
is recorded in that result's `integration_variant`; the preparation equations,
cue, mesh, tolerances and steady-state criteria are unchanged. This case does not
determine the finding of nonmonotonicity: all eight decreasing adjacent pairs occur
elsewhere. Retry outputs and logs are retained for audit.

`manifest.json` identifies the three reused and 18 new cases.
`frozen_provenance.json` hashes the full-sweep launcher, quick-test implementation,
diagnostic imports and original sources. Individual new records include SLURM IDs
and source hashes. The BDF retry records its own source hash additionally.

Reproduce using the existing scaling_sandbox Python environment and
`../run_steady_cue_full.py` (`prepare`, `run --slot N`, `collect`), with the
`../steady_cue_full.slurm` array. The exceptional length is reconstructed by
`../retry_steady_cue_bdf.py`. No original endpoint CSV, figure, or frozen source
was overwritten. This full sweep is retained as a candidate protocol investigation.
