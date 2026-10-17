# Native noisy-ring implementation

This is an execution and archive contract for the experiment in
[DESIGN.md](DESIGN.md). It does not establish a size effect or a migration law.

The source is
[`vision_noise_20261017.py`](../../../../code/mechanisms/vision_noise_20261017.py).
The exact pilot source and inputs are frozen in
[`vision_noise_pilot_20261017`](../../../../data/runs/vision_noise_pilot_20261017/).
The [frozen manifest](../../../../data/runs/vision_noise_pilot_20261017/manifest.json)
declares all 30 jobs and their source/input hashes. SLURM array 21126258 contains
12 noisy cases, 12 deterministic mean-input cases and six zero-input cases;
each uses one CPU, 3 GB memory and small QoS. No main ensemble is part of this
pilot launch.

## Driving process

The periodic Gaussian convolution is discrete-L2 normalized and retains its
spatially uniform mode. The filter has standard deviation $\ell/\sqrt 2$, giving
a Gaussian spatial covariance width $\ell=1.5$. The exact finite-ring covariance
is exported. A stationary OU path is generated on a common 2.5-unit knot grid;
the primary 5-unit knots are every second point. This permits a coupled temporal
representation comparison without substituting an independent random path.

Variance-normalized linear interpolation has unit Gaussian point marginals and
therefore preserves the declared ensemble mean of the lognormal input at every
phase. Its temporal covariance between knots depends on interpolation phase.
It is not an exact stationary continuous-time OU process. The audit records
examples of that covariance, exact Gaussian and kernel hashes, and the generator
and parameter values. No random draw occurs inside the adaptive ODE solver.

Input audits retain both local statistics and raw integrated zeroth/first
angular input modes. Their variances need not be constant across circumference.
Neither local input fields nor their whole-ring modes are normalized using the
realized noise or cell circumference. Uniform mean controls use the same native
activation channel and can change patterned activity through patterned inactive
substrate.

## Worker archives

Workers are `workers/0000` through `workers/0029`. Each `status.json` names its
manifest job, terminal/completion status, runtime, accepted-step concentration
minimum and worst conserved-pool error, plus hashes of all other worker files.
Successful status is `status="complete"`, `complete=true`, `producer="terminal"`.
Failures preserve their partial output and exact last accepted full state.

- `trace.npz` has a uniform 5-unit `time` grid including zero and 4000. It retains
  raw `M0`, `M1_real`, `M1_imag`, `M2_real`, `M2_imag`, `P_real`, `P_imag`, `q1`,
  `q2`, contrast, local region/peak descriptors and sampled input statistics.
- `primary.npz` retains `time` every 20 units, `primary[frame,cell]`, cell centres
  `x` and `circumference`.
- `full_states.npz` retains complete native fields every 200 units, including
  zero and the prescribed end. `final_state.npz` records the exact final accepted
  solver state and its time.
- `scalars.csv` duplicates the scalar readouts. `local_patterns.jsonl` retains
  all circular region centres, widths and excess amounts plus all peak values,
  positions and prominences. Native prominence 0.05 is the frozen descriptive
  peak threshold.
- `noise_audit.json` records generator hashes and parameter/statistical audits.
  `noise_modes.csv` records integrated input modes and local input statistics at
  every forcing knot. Noise can be regenerated from the frozen source and seed.

Raw moments use $\theta=2\pi x/C$ and include the physical quadrature factor
$\Delta x$. $P=M_1/M_1(0)$ uses a complex division, so its initial value is one.
There is no framewise background subtraction or concentration clipping.
Undefined mass-normalized moments carry an explicit boolean flag and zero
sentinel; an undefined initial first-moment reference is a worker error.

Dominant-region position is the circular excess-weighted centroid of the
largest current half-contrast region. Its unwrapped displacement is descriptive:
the dominant region can change identity. Multi-region, multi-peak and weak-contrast
flags preserve that ambiguity; a dominant-site jump alone is not a replacement
claim. Full primary histories permit independent periodic alignment and site
tracking.

## Checks and numerical comparisons

The [frozen preflight](../../../../data/runs/vision_noise_pilot_20261017/results/frozen_preflight.json)
checks the actual frozen implementation/native imports, source/config/copied
input hashes, convolution covariance, discrete translation covariance, exact
coarse/fine shared knots, interpolation continuity/unit variance, all nine
retained basal references and native activation-channel/pool identities.
It performs no native trajectory integration. The read-only external basal lookup
is stated explicitly in that record because the six-case pilot does not copy
the three unused C40 basals.

ODE sensitivity compares the same forcing realization at tolerances tenfold
tighter and maximum step one against the primary result. Temporal representation
sensitivity compares the coupled 2.5-unit interpolant against the 5-unit
interpolant, both using that tighter solver. These are different checks: the
second changes the declared driving interpolation and is not an ODE error.
The audit cases are fixed at C60 and fraction 0.01 in each model. Their
configurations and sources are frozen separately from the pilot.

## Measured pilot execution

All 30 prescribed pilot trajectories completed at 4000 native time units.
The slowest worker used 119.59 seconds; the sum of worker wall times was
597.17 seconds. Worst accepted-step conserved-pool error was
$2.61\times10^{-12}$, and the minimum native concentration was positive
($2.30\times10^{-31}$). These are execution diagnostics, not evidence for a
particular size response.

The reproducible
[`vision_noise_execution_audit_20261017.py`](../../../../code/mechanisms/vision_noise_execution_audit_20261017.py)
verifier checks every output/source/input hash, exact time grids, full final
state agreement, independent raw moment/polarity recomputation, sampled native
budgets and exact regenerated Gaussian/coupling hashes. All 30 cases passed;
the [execution verification](../../../../data/runs/vision_noise_pilot_20261017/results/execution_verification.json)
records the casewise results and measured allocated storage (approximately
69.4 MiB at the verification time).

The [authorized main configuration](MAIN_CONFIGURATION.json) describes the
bounded 32000-unit experiment, with scalar cadence 5, primary cadence 40 and
full-state cadence 1000. It was frozen and launched as SLURM array **21126522**
after the final actual-pilot critique found no essential blocker and the
[root authorization](MAIN_ROOT_AUTHORIZATION.json) was recorded. The
[main manifest](../../../../data/runs/vision_noise_main_20261017/manifest.json)
contains exactly 216 independent noisy realizations, 18 matched mean controls
and nine zero-input controls. The prior 16000-unit proposal is preserved
separately as superseded. The
[resource projection](MAIN_PROJECTED_RESOURCES.json), based on measured pilot
file blocks and runtimes, estimates 3.76 GiB and about 23 summed CPU wall hours;
the slowest projected worker is approximately 16 minutes. Compression and
runtime can change with the longer dynamics; these estimates are not measured
main output or scheduler completion times.

The main
[plan freeze](../../../../data/runs/vision_noise_main_20261017/plan_freeze.json)
preserves 11 scientific records: the design, fresh design/pilot critiques,
final pilot report, numerical comparison and execution audits, analysis/verifier
sources and exact root authorization. Their copied hashes and the plan-record
digest are checked by every worker. The analysis contract is explicit in the
frozen configuration: maximum lag 8000, full/first/last-half windows, whole-seed
bootstrap rules, direction gates, scalar-cadence site diagnostics and descriptive
medoid selection. Main
[producer status](../../../../data/runs/vision_noise_main_20261017/status.json)
and exact worker records distinguish live work from completed output.

A separately recorded
[analysis numerical amendment](../../../../data/runs/vision_noise_main_20261017/results/analysis_numerical_amendment/AMENDMENT.md)
addresses floating-point cancellation in the original correlation implementation
for controls with extremely small late first moments. The frozen producer,
realizations, original analysis source and scientific estimator remain preserved.
The amendment contains a separate corrected source and exhaustive direct-product
regression, and requires its own fresh review before final statistics. Raw native
moment/noise verification does not use the affected lag-correlation arithmetic.
