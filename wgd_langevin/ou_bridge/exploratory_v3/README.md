> Active campaign (2026-09-26): [direct WGD fits with persistence-only differences](../../persistence_only/results/report.html). OU-to-polarity fitting and transfer have been abandoned. The material below is historical.

# Exploratory length trends and reference-only WGD predictions

This is a new descriptive pass. It replaces held-out-first presentation with freely fitted spatial parameter trends, and replaces joint WGD calibration with a fixed length intervention after reference-only calibration.

Open `results/spatial.html` and `results/msd.html`. They are regenerated as completed jobs arrive; manual refresh loads the latest files. Timestamps and completion counts are shown. No browser refresh timer.

## Spatial fits

Nine lengths for each source: 5, 7.5, 10, 12.5, 15, 20, 25, 32, 39 µm; DeBelly uses 21, 23, 25, 27, 29, 31, 33, 36, 39 µm. Existing trajectories are reused at the previous five lengths; four intermediate lengths per source are newly simulated. No spatial mechanism constants change.

Every OU family is fitted separately at every length using all seven protocols and the weak-plus minus weak-minus contrast. Each of those eight observables has equal total weight; trajectories retain early/late phase balancing. All active parameters are free within the existing physical/numerical bounds; none is shared across lengths. Eight starts per fit, including nearest prior spatial fits as initialization, with up to 240 least-squares evaluations each. There are 180 independent fits. Near-best alternatives (within 5% objective) are displayed rather than interpreted as confidence intervals. This pass does not impose new power laws.

## Reference-only transfer

The 20 existing source × OU candidate power laws are retained for this visualization. Their exponents were learned in the earlier spatial campaign. All identifiable reference baseline parameters and a shared physical amplitude are reoptimized using **only reference observations**, with six independently initialized starts and higher-sample polishing. No previous WGD fit is used as an initialization. The five source laws in each family share the same reference-calibration problem, so 20 distinct baseline fits (four families × five studies) yield 100 predictions (20 candidates × five studies).

The user selected exactly sqrt(2) for the WGD/reference length ratio. Each affected parameter is multiplied by sqrt(2) raised to the selected exponent; all other parameters, noise and physical readout scale stay fixed. No WGD statistic is fitted or used to select curves. Unforced gain changes have no observable effect; Hosny's fixed directional cue identifies force = gain × cue.

Gaussian and memory MSDs use exact stationary formulas. Nonlinear trajectories use the previous isotropic two-vector closure, independent cells, maximum integration step 0.01 hours (reduced for fast transferred times), and burn-in max(96 hours, 16 times the largest tp/tm). MSD is averaged across cells and all origins on a 1/6-hour sampling grid out to 24-hour lag. Main curves use 1024 cells and alternatives 512. All curves at different lengths share random seeds to reduce Monte Carlo differences. Up to three alternative near-optimal reference fits produce descriptive envelopes, not confidence bounds. Numerical validation of reference summaries uses a fresh seed and finer integration.

MSD measurements are overlaid only where they actually exist in the transcribed data (Xuan). Other studies supply reference speed or path/directness/quantile observations, not measured MSD curves. Jemaà speed-only calibration cannot identify persistence. Finite bounds and local multistart search limit reference fit quality; inspect the reference-observation table, especially Mallin. A plotted WGD curve is a conditional prediction, not a claimed successful fit.

The initial cancelled reference submission used an old joint fit as one initialization. Its one completed Mallin/Gaussian result was moved to `results/archive_initialization/` and rerun with reference-independent seeds; it is excluded from reports.

## Jobs

- Spatial generation: 20596704 (45 lengths, 25 reused).
- Free fitting: 20596749 (180 fits after spatial completion).
- Reference-only fits: 20596770 plus replacement task 20596771 (20 fits); 20596750 cancelled.
- Report producer: 20596969.

Existing `scaling_sandbox` Conda environment; SLURM red/small, one CPU, 3 GB per task, no concurrency throttle. Previous reports and results are retained.

The few nonconverged best spatial starts received up to 960 continuation evaluations (array 20596971). All 180 best spatial solutions subsequently converged. Hosny/self-activation reference polishing received a separate continuation (20596979); the report explicitly flags a remaining iteration limit if that continuation does not converge. Alternative-start scatter and envelopes remain descriptive, not parameter uncertainty estimates.
