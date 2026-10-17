# Poisson-kick velocity model: original restricted knockdown fit

The current exhaustive comparison is [here](exhaustive/README.md), with its
[Pareto PNG](exhaustive/results/hosny_poisson_pareto.png).

Original restricted-model figure: [hosny_poisson_joint.png](results/hosny_poisson_joint.png).
Only supplied `literature_review/hosny_2026.md` observations are used, transcribed in
`../data.py`. Previous Gaussian OU and constant-speed fits are retained separately.

## Model and interpretation

Position integrates velocity. Between events, dv/dt = -v/tau, with gamma = 1/tau.
Independent Poisson events arrive at rate lambda (/h). Every event adds a velocity
vector: magnitude Rayleigh(scale A, um/h), angle von Mises(concentration kappa)
centred on the oxygen-rich +y direction. No independently fitted directed velocity.
All cells in a condition have the same parameters. Normoxia has kappa=0.

Tau controls decay of the whole velocity vector, not an independently specified
turning resistance. Lambda controls physical kick frequency, independent of any
numerical or camera timestep. Kicks are uncorrelated events; lambda is not a
forcing correlation time. A is a velocity-jump scale, not stationary velocity SD.

Stationary mean forward velocity = lambda*tau*A*sqrt(pi/2)*I1(kappa)/I0(kappa).
The trace of the stationary velocity covariance = lambda*tau*A^2. Its lagged
covariance decays as exp(-lag/tau). These analytic quantities validate the simulator
and supply the MSD prediction, but quartile fitting uses simulated trajectories.

The simulator integrates exponentially decaying velocities exactly between events
and frame times. Initialization sums past kicks to 14 tau plus the first older
kick; neglected mean and variance tails are bounded relative to the stationary
values by exp(-14) and exp(-28). Von Mises inverse-CDF sampling uses a dense
numerically integrated table. No Euler timestep or Gaussian diffusion replacement.

## Joint parameter sharing

For each of non-PACC and PACC, fit six parameters simultaneously:

- Baseline tau, lambda, A, kappa.
- HIF KD kappa, constrained no greater than baseline; other parameters shared.
- RhoA KD tau, constrained no greater than baseline; other parameters shared.

Baseline data include untreated, HIF siNTC and RhoA siNTC as three separate sets of
constraints on the same parameters. HIF KD and RhoA KD observations also constrain
the parameters they share with baseline. No baseline-first/frozen-baseline stage.
Normoxia fits its own three parameters with uniform kick directions.

Expanded models are separate joint refits: add HIF tau; add RhoA kappa; add both;
release all KD sharing/inequality constraints. Added parameters in these expansions
are unrestricted within the search bounds. The primary model remains the main
plotted fit regardless of the expanded model scores. Scores are descriptive, not
likelihoods, significance tests or penalized model selection.

## Observation process and optimization

16-hour trajectories are sampled every 30 minutes, matching the notes. The fitted
quantities are accumulated path, directness (endpoint/path), and parallel fraction
(forward displacement/path), each at Q25, median and Q75. Squared residuals are
scaled by each observed IQR (floors 1 um or 0.1), then averaged over quantiles,
metrics and source datasets. The three baseline source datasets each have the same
weight as a KD dataset. Quartiles reflect population variability, not standard
errors. A shared amplitude is profiled analytically from all relevant path targets
at every optimizer evaluation; it is still jointly fitted, not fitted beforehand.

A 27 x 21 x 25 bank in tau, lambda and kappa uses 2,048 cells/point; interpolation
uses log(tau), log(lambda), log(1+kappa). Bounds: tau 0.005–16 h; lambda 0.12–24 /h;
kappa 0–8. Three differential-evolution seeds and local polishing are used. Nested
model candidate parameters are retained so releasing constraints cannot worsen the
bank objective. Independent 32,768-cell predictions check every compared fit.
Numerical interpolation/Monte Carlo discrepancies are saved, not treated as data
uncertainty. No claim of a proven global optimum or parameter confidence interval.

Baseline-rate profiles fix lambda and jointly refit all other primary parameters,
including amplitude and the KD modifications. They diagnose rate/amplitude tradeoffs
without establishing statistical confidence. Analogous baseline-memory profiles
fix tau and jointly reoptimize the remaining shared and KD parameters. Medians-only, omission of the parallel
fraction, and alternative parallel normalization (projection/endpoint) are separate
joint sensitivity fits. The notes do not resolve the normalization denominator.

## Checks and limits

Endpoint and perpendicular fractions are held out of fitting and displayed alongside
all source observations. Rayleigh p-value simulations (60 independent cells) and the
reported qualitative decrease of positive/high-endpoint PACCs are checks, not loss
targets. The notes' biological replicate clustering cannot be reproduced without
raw data. The high-motility threshold is the simulated shared-control PACC endpoint
Q75; positivity is the sign of the parallel displacement.

MSD curves and cue-removal simulations are explicitly labelled predictions, not
reported data. The model has no spatial or temporal oxygen response, cell area,
area-based retrospective PACC selection, or replicate effects. The supplied radial
FMI means and qualitative temporal response cannot constrain this homogeneous
stationary model without further assumptions. No numerical area/motility correlation
targets were supplied. No literature speed/persistence estimates are invented.

Nine checks test analytic mean displacement, MSD, stationary variance, velocity
covariance, identical event trajectories under changed sampling, amplitude scaling,
directional fraction identity, von Mises angular moment and unbiased mean motion.

## Reproduction

From the repository root, using the existing `scaling_sandbox` environment:

    sbatch wgd_langevin/hosny/poisson_velocity/run.slurm

The script requests 4 CPUs, 6 GB, red partition/small QoS and runs the optimization
then the PNG renderer. Source hashes, bounds, seeds and objectives are recorded in
`results/metadata.json`. CSV/JSON files retain predictions and optimization evidence;
PNG is the only figure export. Bank reuse checks the grids and population size.

The initial lower tau bound (0.04 h) was expanded to 0.005 h after a boundary
solution. Very short memory times trade off with large kick amplitudes; the
30-minute data do not establish the instantaneous speed spikes implied by any one
point estimate. See `results/memory_profiles.csv` and panel G. The final analysis
was run as SLURM job 20511938, reusing and extending the bank from job 20511931.
