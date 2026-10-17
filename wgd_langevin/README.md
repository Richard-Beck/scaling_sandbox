> **Current manuscript migration evidence:** the [3 October refined OU analysis](refined_ou_20261003/README.md) supplies [main_03_migration_OU](refined_ou_20261003/figures/main_03_migration_OU.png), also copied into [current manuscript figs](../current%20manuscript%20figs/README.md). This is the starting point for the evidence motivating increased migration persistence after WGD-associated enlargement. Read its reviewed observations, conditional persistence/noise fits and source limitations before interpreting older campaigns. No completed manuscript is implied.

The earlier [direct WGD fits with persistence-only differences](persistence_only/results/report.html)
remain a separate model comparison. OU-to-polarity fitting and parameter transfer
were abandoned; the current polarity-mechanism work seeks a physical explanation
for the migration phenotype rather than treating fitted OU memory as an already
identified molecular parameter.

# Retained historical source

The earlier active-Brownian, Gaussian/Poisson, feedback, persistence-only and
OU-to-polarity branches are preserved as reviewed scientific source and
documentation. Their fitting code, tests, launchers and derivation notes record
how the analysis developed; the refined OU workflow above supplies current
manuscript evidence. The exact source selection is recorded in
[`reports/publication_historical_sources_20261017.txt`](../reports/publication_historical_sources_20261017.txt).

Historical raw results, native trajectories, caches and report delivery copies
remain local and ignored. A source checkout does not contain every output
linked by older READMEs, and those links do not establish current scientific
acceptance or complete rerun coverage. Preserve source/input hashes, the
configuration and a separate verified data location when reusing an old result.
The existing root model invariant tests can run with the modelling environment;
older branch tests should run in their own workstream context.

# Earlier model comparisons

Xuan, Mallin and Jemaà now use **Gaussian velocity OU**, with two parameters per
population: velocity memory and Wiener noise strength. See the [methods](gaussian_literature/README.md),
[results](gaussian_literature/REPORT.md), and [overview PNG](gaussian_literature/results/literature_gaussian_overview.png).

Hosny retains its Poisson model and [shared-mechanism presentation](hosny/README.md).
The previous [Poisson literature analysis](poisson_literature/README.md) is preserved
as a historical comparison. The material below and root `results/` document the
superseded constant-speed analysis.

---

# WGD-associated migration as active Langevin motion

A standalone, two-dimensional point-particle model constrained by the repository's
literature observations. Stable tetraploid clones, giant polyploid cells, and
treatment-induced PACCs are analyzed separately. This is a phenomenological
comparison, not an inference that all differences were caused by WGD.

Start with [the results report](REPORT.md), [observation/model figure](results/observations_vs_model.png),
and [parameter profiles](results/parameter_profiles.png).

## Model

Each cell has a position r, heading angle theta, and constant intrinsic speed v.
The stochastic equations, in plain text, are:

    dr = v * (cos(theta), sin(theta)) * dt + sqrt(2 * Dt) * dW_position
    dtheta = sqrt(2 * Dr) * dW_angle
    persistence time tau = 1 / Dr

W increments are independent Gaussians with variance dt. Distances are micrometres
and times hours. Initial headings are uniform. There is no inertia, directed cue,
cell interaction, proliferation, or confining geometry. The primary fit sets Dt=0
and localization error=0 because these are not independently constrained here.
The simulator supports both as separate quantities for future trajectory fitting.

Optional lognormal cell-specific speeds have mean v and coefficient of variation
CV. Their log-standard-deviation is sqrt(log(1+CV²)). Speed and heading dynamics
are independent; there is no fitted correlation between speed and persistence.

For an unconfined population, the exact ensemble MSD is:

    MSD(t) = 4 * Dt * t + 2 * v² * (1 + CV²) * tau² * (t/tau - 1 + exp(-t/tau))
    long-time effective diffusion = Dt + v² * (1 + CV²) * tau / 2

Angular increments are sampled exactly and position is integrated using midpoint
headings. Internal timesteps are no greater than tau/8 and cadence/32 in the
fitting pipeline. Observations are sampled at each study's imaging cadence.
Measured speed is mean sampled step length divided by cadence; it is generally
smaller than v. Directness is calculated per track, then averaged (never a ratio
of population-average endpoint distance and population-average path length).

If v is interpreted as propulsion force / effective drag, migration alone does
not separate those two mechanical parameters. Neither cell size nor doubling of
DNA automatically fixes a parameter in this model. Nuclear microrheology is not
substituted for whole-cell motility noise.

## Evidence and fitting

[data/observations.csv](data/observations.csv) contains the actual input values,
figure-reading intervals, experimental cadence, duration, roles, and provenance.
`source` paths are relative to this directory. The table is a structured extraction
of pre-existing local review notes, not new digitization or raw experimental data.
Approximate central values are treated as means for this first analysis; exact
summary definitions should be verified against original figures before stronger
statistical inference. RKO D1/T1 and sarcoma Cl1/Cl8 are the primary clone contrasts;
other clones are retained as context rather than pooled as independent studies.

Primary sources and local records:

- [Jemaà 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10531202/),
  [local extraction](../literature_review/jemaa_2023.md): RKO selected clones;
  speed at 1h cadence over 24h. Most direct induced-tetraploid system here,
  but independent clone selection remains a confounder.
- [Jemaà 2017](https://www.oncotarget.com/article/14291/text/),
  [local extraction](../literature_review/jemaa_2017.md): selected sarcoma clones;
  speed at 30min cadence. Nominal 24h assay, with a 20h discrepancy in the text.
  Mean sampled speed calibration does not depend on this duration under stationarity;
  illustrative distribution widths do.
- [Xuan 2018](https://www.nature.com/articles/s41598-018-29817-5),
  [local extraction](../literature_review/xuan_2018.md): giant polyploid cells;
  speed and late-lag MSD fit; early MSD and crossover excluded from fitting.
  Slopes are diagnostics because their exact fit window is unavailable.
- [Mallin 2023 local extraction and primary-source links](../literature_review/mallin_2023.md):
  cisplatin-induced PACCs; path length and directness fit; endpoint displacement
  is an excluded-metric check on the same tracks, not independent validation.

Xuan 2020, collective closure, membrane crossing, metastatic incidence, shape,
and stiffness are outside the initial single-cell model. Xuan 2018's slow fraction
is retained but not fitted: the threshold is one cell length per hour and group
cell lengths were not reported. A switching model cannot be identified from those
occupancies alone. No unobserved switching rates are invented.

Fit score: sum of squared log prediction/observation residuals divided by squared
half-log-widths of the reading intervals. These widths are NOT sampling SEMs.
The score is a descriptive distance, NOT a likelihood, chi-squared statistic, or
basis for p-values/AIC. Compatible points satisfy every fitting interval separately.
Range endpoints are grid projections, not confidence intervals or independent
parameter boxes. Values within a marginal range need not be compatible in every
combination.

For each assay, 81 persistence values are simulated with 2,048 cells and a common
random seed, then interpolated on log tau to 321 values. Propulsion is profiled on
1,401 log-spaced values from 0.1 to 500 µm/h. Tau spans cadence/64 to 48h.
Grid limits are computational choices, not biological limits. Jemaà persistence
ranges reach both edges: numerical best-fit tau values there have no biological
interpretation. Mean speed alone leaves a continuous family of solutions.

Comparisons include shared parameters, speed only, persistence only, both speed
and persistence, and shared persistence with changes in mean speed plus polyploid
heterogeneity (reference CV=0; polyploid CV scanned 0 to 2.5). The last model is an
alternative explanation, not a new observation. Fixed CV sensitivity and doubled
reading widths are also recorded.

The Mallin sample-size sensitivity adds 1.96 times the model-predicted standard
error of mean directness (50 reference, 49 PACC tracks) to each reading interval.
It is explicitly conditional on the model, independent tracks, and interpreting
the plotted central value as a mean. It is NOT an empirical confidence interval
and does not account for biological replicate correlation. It is reported separately
from the narrower reading-only results.

## Reproduce

Use the existing Conda environment; no installation is needed:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/run_analysis.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/validate_analysis.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python -m unittest discover -s wgd_langevin/tests -v

Run from the repository root. Each command uses a single compute thread and small
memory; this pilot is below the repository's SLURM threshold. Any expansion beyond
those limits should use SLURM. Source hashes, package versions, seed and arguments
are in `results/summary.json`. Validation uses an independent seed and 8,192 cells,
with coarse/fine integration. Core tests check ballistic motion, active MSD,
diffusion/localization noise, and heterogeneous-speed second moments.

The compact CSV/JSON outputs and figures are retained; full simulated banks of
trajectories are regenerated rather than stored. `heterogeneity_sensitivity.png`
shows illustrative CV scenarios, not fitted tail distributions. Published maxima
without sample sizes cannot identify a tail likelihood or a lognormal CV.
