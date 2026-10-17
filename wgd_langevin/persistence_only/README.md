# Direct WGD fits with persistence-only group differences

This is the active four-variant OU campaign. It supersedes all OU-to-1D-polarity
calibration and transfer campaigns. Previous outputs are preserved as historical
records and are never loaded as data, starts, constraints or priors here.

Open [the report](results/report.html), [findings](results/FINDINGS.md),
[comparison CSV](results/comparison.csv), and [prediction CSV](results/predictions.csv), and [current best parameters](results/current_best_parameters.csv).
The parameter export includes both hypotheses; visible parameter tables emphasize
persistence-only alternatives. Pending fits have no invented parameter estimates.
The report is regenerated atomically whenever a fit completes, records its
creation time, and has no automatic browser refresh.

## Question and scope

Can each WGD-associated migration contrast be described by changing only persistence
times? Each experiment has its own shared baseline parameters. Within its WGD−/WGD+
pair, tp may differ for Gaussian OU; tp and tm may differ for the other variants.
All remaining parameters, including noise, feedback, saturation, self-activation,
readout amplitude and directional forcing, are shared within that pair. A nested
null shares the persistence times too. No parameters are tied across unrelated
experiments. WGD labels abbreviate the reported reference/polyploid or PACC groups;
they do not imply all studies are clean causal WGD experiments.

Nine separate paired experiments: Xuan 2018; Mallin 2023; Jemaà 2017; Jemaà 2023;
Hosny baseline; Hosny HIF control; Hosny HIF knockdown; Hosny RhoA control; Hosny
RhoA knockdown. Hosny treatments are fitted independently, not predicted by a
knockdown mechanism. Unpaired normoxia is excluded. Jemaà uses the designated
Cl1/Cl8 and D1/T1 contrasts; context clones are not pooled as independent replicates.

Inputs are the repository's existing transcriptions in `../data/observations.csv`
and `../hosny/data.py`, with their recorded provenance and approximate reading
intervals. No new digitization or raw-track inference is claimed. Imaging cadence
and duration: Xuan 1/6 h and 12 h; Mallin 0.5 h and 24 h; Jemaà 2017 0.5 h and
24 h; Jemaà 2023 1 h and 24 h; Hosny 0.5 h and 16 h. Jemaà 2017 has a documented
20/24 h source discrepancy; stationary mean speed does not depend on track duration.

## Four variants and stochastic definitions

The existing four dynamical variants are retained; numerical functions have been
copied into this standalone campaign so it imports no polarity calibration code.
Time is hours; spatial readout is micrometres. p and q are 2D latent vectors.
Velocity is a shared readout amplitude times p, and position integrates velocity.

```
Gaussian: dp = (-p/tp + force*ey) dt + dW

Memory mean dynamics:
  dp/dt = (-p + r*q)/tp + force*ey
  dq/dt = (p-q)/tm

Coupled:
  dp = [(-p + r*q/(1+norm(q)/vs))/tp + force*ey] dt + eta*dW
  dq/dt = (p/(1+norm(p)) - q)/tm

Self-activation adds to dp/dt:
  self_a*p*norm(p)/(self_K^2 + norm(p)^2)
```

Linear memory retains the existing stationary substitute-covariance realization:
for each spatial component A is the stated 2x2 drift matrix, S is diagonal with
entries tp/2 and tp^2/(2*r*tm), and Q = -(A*S + S*A-transpose). Gaussian forcing
has covariance Q. This specifies a correlated-noise closure; its derived covariance
changes with the allowed time constants. There is no additional fitted group-specific
noise parameter. It is not an independently driven q-noise model. The instantaneous
p-noise convention is one for the linear variants; a fitted shared readout amplitude
restores a shared physical noise scale. Nonlinear models use independent isotropic
white noise in p only, with fitted shared eta. Fixing ps=1 removes an exact latent
scale redundancy; saturation, eta and readout amplitude remain shared.

Force is fitted only for Hosny, enters along the recorded parallel direction,
and is identical for both groups. Cue gain and cue strength are not separately
identifiable. tp and tm are model relaxation times; feedback can make the observed
autocorrelation time much longer. All four models are phenomenological candidates.

## Fitting, numerical checks and identifiability

72 fits: 36 shared nulls followed by 36 persistence-only alternatives. Eight fresh
starts per fit, including broadly randomized parameter values, with 120 evaluations
per start. The two best are polished with up to 80 evaluations each. Each alternative
also polishes its null solution, ensuring the nested baseline is tried explicitly.
No previous fit supplies a start. Common random numbers apply across groups and
optimizer evaluations; groups with identical parameters produce exactly identical
predictions. Each saved result independently audits that only tp/tm differ.

Search bounds (hours for times): tp, tm 0.05–12; memory r 0.005–0.99; nonlinear
r 0.005–8; vs 0.03–10; eta 0.03–5; self_a 0.0001–2; self_K 0.02–2; force
0.0001–10. These are finite computational search bounds, not biological priors or
proofs of exclusion. Results flag parameters within 3% of either end of their log
range. Parameters optimized in log coordinates; shared readout amplitude is profiled
analytically. Linear limits approached at small positive r/self_a are not exact
zero-boundary submodels.

Linear displacement covariances, unforced mean speeds, MSDs and mean endpoints are
analytic; directness and forced quartiles use sampled Gaussian increments. Linear
search/polish/validation use 512/2048/8192 independent cells. Nonlinear simulations
use stochastic Heun at a fixed 0.01 h maximum step and 48 h burn-in, with
96/512/2048 independent cells. Validation fixes all fitted parameters and readout
amplitude and uses an independent seed. Nonlinear checks additionally use 0.005 h
steps and 96 h burn-in. Changing discretization also changes random trajectories,
so differences include sampling variation. These checks diagnose numerical
sensitivity; they do not establish stationarity for every possible long-memory fit.

Xuan fits speed and late MSD (11.111 h); early MSD is held out at 1/6 h. Mallin fits
path and directness; endpoint is held out. Jemaà fits speeds only. Hosny fits
path/directness/parallel quartiles; endpoint and perpendicular summaries are held
out. Held-out metrics come from the same experiments and are correlated checks,
not independent biological validation. Missing Hosny quartiles remain missing.

Non-Hosny objectives sum squared log residuals scaled by half log reading-interval
width. Hosny uses squared residuals scaled by observed IQR (floor 0.04; missing IQR
scale 0.1). All are descriptive losses, not likelihoods: no AIC, p-values or sampling
confidence intervals. Scores compare variants only within the same experiment.
Adequacy for Xuan/Mallin requires every fitted metric inside its reading interval
both at fit and all validation settings. Held-out agreement is shown separately.
Hosny has no statistical adequacy threshold.

The report exposes local Jacobian rank after profiling amplitude, near-bound
parameters, and ratio ranges among search starts within one objective unit of the
best search start. Those are identifiability diagnostics, not confidence intervals.
Speed-only data cannot identify persistence or distinguish dynamic mechanisms.
Sparse summaries leave many nonlinear parameters underdetermined even when fitted
measurements agree. Finite multistart misses do not establish impossibility.

## Reproduce

Use the existing `scaling_sandbox` environment; no dependencies were installed.
Run from the repository root. Jobs use SLURM red/small, one CPU, 3 GB, two hours;
there is no array concurrency throttle. The null array must finish successfully
before the alternative array runs.

```
python -m unittest discover -s wgd_langevin/persistence_only -p 'test_*.py' -v
python wgd_langevin/persistence_only/campaign.py --manifest
sbatch --array=0-35 wgd_langevin/persistence_only/run.slurm
sbatch --dependency=afterok:NULL_JOB_ID --array=36-71 wgd_langevin/persistence_only/run.slurm
python wgd_langevin/persistence_only/figures.py
python wgd_langevin/persistence_only/publish.py
```

Result files are resumable only with identical code and input hashes. Preserve
previous results and use a new output directory after changing inputs or model code.
The manifest records source hashes and numerical package versions; each fit records
its SLURM job ID, starts, convergence, constraints, bounds and numerical checks.

## Xuan memory example tracks

[Track report](results/xuan_memory_tracks.html) shows 32 independent stationary
tracks per condition at the original 10-minute cadence over 12 hours, plus
unselected individual examples with 1-minute resolution. It uses the saved
Xuan memory fit without refitting. [Reproduction script](xuan_memory_tracks.py)
exports PNG/PDF figures, coordinates, and checks against analytical speed/MSD.
Run with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python wgd_langevin/persistence_only/xuan_memory_tracks.py`.

## Ordered-time Xuan memory background fit

[Live constrained fit](results/xuan_memory_ordered/report.html) fits Xuan alone
with 0 < tp− < tp+ < 12 h and 0 < tm− < tm+ < 12 h. It uses 32 starts,
500 evaluations per start, shared feedback/readout, and analytic speed/MSD.
A 1e-6 h numerical floor implements the positive open interval; it does not
retain the original 0.05 h floor. The best completed start and report are
checkpointed after every start. The original fitted tracks are unchanged.
