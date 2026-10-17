# Current presentation: four shared-mechanism representatives

The current presentation enforces matching HIF/RhoA parameter subsets across PACC
and non-PACC, then compares the four baseline parameters across the shared frontier
and within-model ratio profiles. See [architecture](poisson_velocity/shared_mechanisms/README.md)
and the [four-representative review PNG](poisson_velocity/shared_mechanisms/presentation/results/hosny_four_representatives.png).
The [presentation notes](poisson_velocity/shared_mechanisms/presentation/README.md)
map the grouped bars, data comparisons, trajectories, MSD and condensed sensitivities.
The [previous compact parameter comparison](poisson_velocity/shared_mechanisms/results/hosny_shared_parameter_comparison.png) is retained.

The [population-specific exhaustive sweep](poisson_velocity/exhaustive/README.md)
and its [overview PNG](poisson_velocity/exhaustive/results/hosny_poisson_pareto.png) are retained.

The [original restricted Poisson comparison](poisson_velocity/README.md) is retained.

The [Gaussian OU revision](ou_velocity/README.md) is retained as a superseded analysis.

The following records the superseded constant-speed analysis.

# Hosny 2026: directed Langevin optimization

The requested standalone artifact is [results/hosny_optimization.png](results/hosny_optimization.png).
All experimental inputs are transcribed from `literature_review/hosny_2026.md`.
No external sources, trajectories or additional data were used.

The model has constant propulsion v, rotational diffusion 1/tau, and angular
restoring drift -kappa*sin(theta) toward the local oxygen-rich +y direction.
Initial headings follow the stationary von Mises distribution with concentration
kappa*tau. All cells in a condition share these parameters. The primary FMI
assumption is projected net displacement / accumulated distance; an alternative
normalization by net distance is analyzed as a sensitivity, not asserted as fact.

For each of 11 conditions, the Q25, median and Q75 of accumulated distance,
directness and parallel FMI are fitted. Each metric contributes its mean squared
quantile error divided by its observed IQR squared. Thus IQR is population spread
and a dimensionless loss scale, not an SEM/CI. Propulsion is profiled analytically;
persistence and steering are optimized on an interpolated simulation bank using
differential evolution with local polishing. Normoxic steering is fixed to zero.

The bank spans tau=0.04–12h and concentration=0–3, with 528 combinations and 1,536
cells each. Independent validation/predictions use 8,192 cells per condition.
Analytic unbiased MSD, stationary directed displacement and per-track FMI identities
pass. Baseline timestep and initial-orientation sensitivities are retained.
There is no likelihood, model-selection significance test, or parameter CI here.

Each knockdown uses its own figure-specific siNTC control. Joint comparisons share
steering *rate* kappa or persistence tau as specified, not the dimensionless product.
Endpoint distance and perpendicular FMI are reserved checks. Rayleigh tests are
simulated for independent 60-cell batches; biological replicate dependence is not
available and is not modeled. The radial profile's group mixture is unspecified
in the notes, so it is displayed as an out-of-scope observation rather than used
as a quantitative spatial fitting target.

The fitted homogeneous model is inadequate for the broad observed path-length
IQRs and misses some reserved checks. It must not be described as reproducing all
the observations or establishing a unique molecular mechanism. In particular,
knockdown steering estimates depend on fitting quartiles versus only medians.

Reproduce from the repository root using the existing environment:

    sbatch wgd_langevin/hosny/run.slurm
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/hosny/plot.py

The completed first optimization was SLURM job 20510306, QoS small, 4 CPUs, 4 GB
request. Runtime for analysis was about 71 seconds after process startup.
Source/input hashes, seeds and settings are in `results/metadata.json`.
CSV/JSON files retain the optimization and prediction evidence; the only figure
export is PNG. Full simulated trajectories are not retained, only 32 illustrative
tracks per condition.
