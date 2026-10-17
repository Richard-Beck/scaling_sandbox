Historical Poisson comparison: superseded for these studies by [Gaussian OU](../gaussian_literature/README.md). Hosny remains Poisson.

# Earlier literature reanalysis with the common Poisson model

Current figure: [literature_poisson_overview.png](results/literature_poisson_overview.png).
This replaces the constant-speed interpretation of the original uniform-medium
contrasts in `wgd_langevin/results`, whose files are retained as historical output.

## One simulator across the project

This pipeline imports `wgd_langevin/hosny/poisson_velocity/model.py` directly.
Between independent Poisson kicks, velocity decays with memory tau; positions
integrate velocity. Kick magnitudes are Rayleigh(scale A) and kick angles are
von Mises(kappa). Initial velocities are stationary to the same truncated-history
accuracy as the Hosny implementation. Every cell within a population shares
parameters; speed varies over time. No permanent speed heterogeneity is fitted.

The original assays are uniform-medium migration, including Mallin's Fig. 2
uniform 20% serum assay. They therefore use the isotropic special case kappa=0.
The fitted parameters are tau, kick rate lambda, and kick magnitude scale A.
Zero bias is a symmetry assumption for these assays, not a measurement of zero
cellular polarity. Higher track directness does not by itself imply chemotaxis.

We retain each study's original cadence and duration: Xuan 10 min/12 h; Mallin
30 min/24 h; Jemaà 2023 1 h/24 h; Jemaà 2017 30 min/nominal 24 h (the supplied
notes also mention 20 h; stationary mean-speed calibration is duration-independent).

## Data and observation model

The original observation table is reused. Xuan fits mean sampled speed and one
late-lag MSD, not an entire MSD curve. Mallin fits mean sampled path length and
mean per-track directness. Endpoint displacement is a held-out summary on those
same tracks. Jemaà only constrains sampled mean speed; its output consists of
viable families, not arbitrary reported persistence estimates.

Reading intervals are approximate extraction tolerances, not standard errors or
confidence intervals. The original loss is retained: sum of squared log residuals
scaled by squared half-log reading widths. This differs from Hosny's quartile loss
because these older inputs are means, not measured distribution quartiles. Model
consistency does not require inventing quartiles for studies that did not supply
them. Independent predictions use 65,536 cells/state; fitting banks use 2,048.

Tau spans 0.005–48 h, rate 0.05–64 /h. Scale is profiled analytically in log space
using the exact linear/quadratic amplitude powers. Grid compatibility checks ask
whether **some common amplitude** satisfies each population's fitting intervals;
they are not confidence regions. Monte Carlo bank error can matter for narrow
intervals, particularly Mallin parental directness, so independent checks are saved.

For Xuan and Mallin, eight joint population comparisons allow none, or every subset
of the three parameters, to change. Parameters not selected to change are shared
and jointly refitted. No WGD effect sign is imposed. These descriptive comparisons
are not likelihood model selection. Jemaà remains speed-only families, without
additional mechanistic fitting.

## Retained checks and ambiguities

Xuan's early MSD, approximate 2–3 h crossover, reported slopes and slow-time fraction
are checks, not silently added precise fitting targets. The slow fraction uses a
common unknown one-hour displacement threshold and separately considers net
displacement versus sampled accumulated distance. The notes do not establish
separate group-specific cell lengths; the old table's assertion of group dependence
is corrected in this pipeline's derived observation copy. Thresholds are displayed
as sensitivity curves, not independently known calibration or a switching model.

For the isotropic stationary model, mean squared displacement is

    MSD(t) = 2 * lambda * A^2 * tau^2 * [t - tau*(1-exp(-t/tau))].

Thus rate and amplitude only rescale the curve at fixed tau. They can alter the
relation between sampled mean speed and MSD through intermittency, but cannot
produce an MSD crossover if tau is shared. The reported Xuan crossover adds
information that fitting two amplitudes alone cannot supply. Mallin directness
also depends on kick frequency, so its previous memory-only inference must be
reassessed under this expanded model.

No Mallin gradient/acrylamide data, Xuan perturbation assays, collective migration,
mechanical measurements or new external sources are added to this scope. These
remain distinct experiments requiring their own observation models and constraints.

## Reproduce

    sbatch wgd_langevin/poisson_literature/run.slurm
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python wgd_langevin/poisson_literature/plot.py

The SLURM job uses red/small QoS, four CPUs and 6 GB. Existing compatible banks are
reused; do not reuse them after changing the simulator, grids or observation
process. Source hashes, bounds and scope are stored in results/metadata.json.
The same nine physical/numerical simulator checks used for Hosny run before fitting.
