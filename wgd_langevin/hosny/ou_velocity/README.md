# Current Hosny model: integrated OU velocity with shared controls

[Standalone results figure](results/hosny_ou_shared_baselines.png).
This supersedes the Hosny constant-speed active Brownian fit in the parent folder.
The original analysis remains available as provenance, not the current model.

## Model and interpretation

    dr = V dt
    dV = -(V - u e_y) dt / tau + sqrt(2 sigma² / tau) dW

V is the two-dimensional velocity. Both speed and direction fluctuate.
Sigma is the stationary standard deviation of each random velocity component;
sqrt(2)*sigma is the random-motion RMS speed shown in panel A. Tau is the
velocity-autocorrelation time. The directional bias u is the stationary mean
velocity toward the locally oxygen-rich +y direction, in µm/h. It is **not** the
angular steering parameter kappa from the preceding model. The relaxation models
motility memory, not physical inertial coasting of a cell.

Every cell in a condition has the same parameters; initial velocities are drawn
from their stationary Gaussian distribution. This is a homogeneous stochastic
population with time-varying speeds, not cell-specific permanent speed parameters.
Position and velocity are jointly sampled using the exact integrated OU transition,
including their correlated innovations. No Euler timestep approximation is used.
The observation cadence is 30 minutes over 16 hours.

## Shared controls and optimization

The only experimental source is `literature_review/hosny_2026.md`.
The existing audited transcription in `../data.py` is reused.

Seven parameter sets fit eleven source conditions:

- O: normoxia.
- N: untreated hypoxic non-PACC plus HIF-siNTC and RhoA-siNTC non-PACC controls.
- P: untreated hypoxic PACC plus both PACC siNTC controls.
- HN, HP: HIF knockdown of N and P.
- RN, RP: RhoA knockdown of N and P.

Untreated and siNTC medians/quartiles remain separate constraints on the shared
baseline. They are not pooled into an invented empirical distribution. No
separate siNTC parameters or experimental source offsets are fitted.

Sigma, tau and u are fitted jointly. Computationally, the bank uses sigma=1 and
beta=u/sigma; sigma is profiled analytically at every tau/beta candidate. This is
a joint optimization, not calibration of no-bias motion followed by a steering fit.

For each data source, Q25, median and Q75 of accumulated distance, directness and
parallel directional fraction enter the loss. Squared quantile residuals are
scaled by the square of that source's IQR, averaged over the three quantiles, and
summed across metrics. Shared-group loss averages the source contributions. IQR
is spread and an objective scale, not a confidence interval. No likelihood,
parameter confidence interval or significance-based model selection is claimed.

The directional fraction is assumed to be the net displacement projection divided
by total accumulated distance, the usual forward-migration-index interpretation.
The source notes leave the normalization unresolved. Refits with that metric
omitted, a different normalization, and medians alone are retained in sensitivity.csv.

Panel C consists of independent constrained refits, not post-hoc ablations. Its
columns specify which parameters may differ between the two compared conditions;
unlisted parameters are shared and optimized, not set to zero. The baseline
comparison contains six datasets (three per population); a knockdown comparison
contains the three shared-baseline datasets and its one treatment dataset. Scores
are averaged over all participating source datasets and the three metrics.
The full-model baseline used for prediction stays the same across comparisons;
constrained pair-specific refits in C are alternatives, not extra baseline states.

## Predictions and verification

Exact Gaussian transition tests verify integrated MSD, mean directed displacement,
velocity autocorrelation, the per-track directional-fraction/directness identity,
the zero-noise deterministic limit, and invariance under finer integration followed
by downsampling. All six pass. Independent predictions use 16,384 cells per state,
versus 4,096 cells per simulation-bank point. The bank contains 49 persistence
values spanning 0.02–32h and 41 bias/SD ratios spanning 0–2. Scalar amplitude is
profiled continuously; no parameter confidence bounds are implied by this domain.

Net displacement and perpendicular directional fraction are prediction checks.
Rayleigh diagnostics are retained in CSV using 60-cell simulated samples, but the
unexplained p-value panel has been removed. The qualitative high-motility subgroup
check is also retained numerically and omitted from the main figure; the notes do
not provide numerical proportions for validation. The figure shows new MSD
predictions (including a labeled >16h extrapolation), traces of fluctuating speed,
trajectory realizations, and the counterfactual removal of directed mean velocity.

The OU model produces a much wider path-length distribution than constant-speed
motion, but still underpredicts much of the reported spread. It also misses some
reserved endpoint checks. One stationary OU timescale couples speed fluctuations
and directional memory; it cannot independently adjust persistent speed episodes.
More elaborate OU/state models or between-cell heterogeneity would be additional
hypotheses, not silently added to this requested revision. The supplied radial
profile and qualitative temporal rise-and-fall in directional response cannot be
explained by this constant-coefficient local model. PACC identity is a retrospective
size/growth classification; these fits do not establish a causal WGD mechanism.

## Reproduce

From the repository root, using the existing Conda environment:

    sbatch wgd_langevin/hosny/ou_velocity/run.slurm
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/hosny/ou_velocity/plot.py

SLURM job 20511831 used four CPU cores, QoS small, and a 4GB memory request.
Metadata, source hashes, data-source mapping, numerical checks, sensitivity refits,
and prediction tables are in `results/`. The only figure format exported is PNG.
