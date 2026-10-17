# Exhaustive joint knockdown mechanisms and Pareto frontier

Standalone figure: [results/hosny_poisson_pareto.png](results/hosny_poisson_pareto.png).
The Poisson-kick simulator, supplied observations, loss and simulation bank are
unchanged from the [Poisson model](../README.md). No external literature was used.

## Search architecture

Each knockdown can change any nonempty subset of tau (velocity memory), lambda
(kick rate), A (Rayleigh velocity-jump magnitude scale) and kappa (von Mises
concentration). There are 15 subsets per KD, hence 225 joint HIF/RhoA combinations
per PACC class. A joint no-change null adds one model per class: 452 fits in total.
The search does not add models in which only one KD has an empty subset. The
available complexity budgets are therefore 0 and 2 through 8.

Every model jointly refits the four baseline parameters and the KD-specific
parameters. Nonselected KD parameters equal the corresponding baseline parameter.
Selected parameters may increase or decrease within the established bounds. Each
untreated/siNTC observation remains a separate constraint on the shared baseline;
KD observations constrain all parameters they share with that baseline.

Complexity counts all KD-specific scalar parameters, including analytically
profiled amplitudes, in addition to the four baseline parameters. It does not
count source datasets or optimizers. Source weights and fitted quartiles are
unchanged. Normoxia is retained from the previous fit and contributes no constant
term to the separately computed PACC/non-PACC mechanism frontiers.

## Numerical procedure

The established 2,048-cell, 27 x 21 x 25 tau/rate/concentration simulation bank is
reused after source-provenance checks. A Numba implementation of its interpolated
loss is checked against the original loss for all 226 sharing structures. This
acceleration changes execution time, not the model, observation process or loss.

Each model uses three differential-evolution searches with polishing, plus local
refinement of embedded smaller-model fits and previous fitted starts. All
subset/superset bank-loss inequalities are checked. Independent 32,768-cell
simulations at each fitted state then evaluate every model on the same data
summaries. These are fresh Monte Carlo samples, NOT new experimental validation
data and NOT uncertainty estimates for the observations.

The displayed frontier is defined using these independent simulated scores. A
point is Pareto-dominated if an equally or less complex model has lower or equal
error, with at least one strict inequality. The best model at every budget is also
reported, even if it is dominated. Numerical noise can remove a more flexible
model from this empirical frontier; do not interpret tiny differences as evidence
that extra parameters physically worsen a perfectly optimized model.

All best-at-complexity models receive a second independent simulation. Differences
between these runs were at most about 0.0006 in the initial completed analysis.
Near-equal model rankings remain uncertain: the top four non-PACC two-parameter
models differ by only about 0.00026. No likelihood, significance test, confidence
interval, complexity penalty or unique mechanism claim is introduced.

## Figure interpretation

- A: all 225 combinations plus the null; independent-score Pareto frontier;
  original signed HIF-kappa/RhoA-tau hypothesis marked with a black cross.
- B: best HIF and RhoA subsets at every parameter budget; filled dots distinguish
  points actually on the frontier. Subsets need not be nested along the frontier.
- C: direct data comparison at the original two-parameter budget (one parameter
  per KD). The best subset is a coloured band; the original hypothesis is a
  dashed IQR/median. Both are regenerated using 65,536 cells/state. This equal-
  complexity comparison is not a claim that two parameters is the global optimum.
- D: fitted parameters for those same two-parameter models, with KD-specific
  entries highlighted; shared entries are literal shared baseline values.
- E/F: homogeneous baseline trajectories and analytic MSD predictions from those
  same illustrated models. Short tracks are prefixes of the long tracks.

Endpoint and perpendicular fractions remain held-out checks. Population IQRs are
not standard errors. Additional Rayleigh and qualitative subpopulation checks are
retained in CSV. Spatial/temporal oxygen responses and cell area are not modeled.
The directional-index normalization and sub-frame memory/amplitude tradeoff remain
limitations of the supplied data, not resolved by exhaustive subset enumeration.

## Output and reproduction

    sbatch wgd_langevin/hosny/poisson_velocity/exhaustive/run.slurm

Uses the existing `scaling_sandbox` environment, four CPUs, 6 GB, red partition,
small QoS. `fits.json` checkpoints completed fits and independent validations;
rerunning reuses these and regenerates checks, detailed predictions and the PNG.
Do not reuse checkpoints after changing model, loss, optimizer conventions or
source data; use a new results directory for such a revision.

`all_models.csv` retains every score; `frontier.csv` and `best_at_complexity.csv`
provide the compact summaries. `fits.json` includes every parameter set and search
seed outcome. `original_reference.json` preserves the original constrained fit.
`metadata.json` records source/bank hashes and the architecture. PNG is the only
figure export. Analytical/model checks and accelerated-loss equivalence are in
`tests.csv`.
