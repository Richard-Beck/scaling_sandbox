# Representative Xuan movie checks

These diagnostics compare **nuclear centroid paths**, not population migration fitting targets. Eight visually accepted tracks cover the same 74 frames (73 intervals): N010 is the apparent giant nucleus; N015, N016, N017, N018, N021, N024 and N025 are the seven smaller nuclei. Ploidy is not measured. Their comparison with the fitted WGD+ and WGD− population models is conditional on appearance, not evidence that their experimental groups are known. The 78 extracted segments are not 78 independent cells. No short, partial, merged, boundary-censored, ambiguous or unreviewed tracks enter these checks; no missing positions are interpolated.

## Outputs

- `xuan_movie_trackplots.png` and `.pdf`: observed equal-window paths and the first unselected simulated paths, with identical spatial limits. Paths begin at a common origin, without endpoint alignment or drift subtraction.
- `xuan_movie_diagnostics.png` and `.pdf`: displacement CDFs, raw step-vector correlation, whitened prediction residual quantiles, step sequences and longest slow/fast episodes.
- `movie_track_checks.csv`: one row per accepted track, timing scenario and localization sensitivity. Includes conditional speeds, centroid-scale sensitivity, residual statistics and episode lengths.
- `movie_correlations.csv`: observed and model vector correlations with finite-trajectory simulation ranges.
- `movie_quantiles.csv`: displacement and standardized residual quantiles.
- `movie_simulated_examples.csv`: the first eight unselected simulated paths per model group for the displayed scenario.
- `movie_provenance.json`: exact fit/input/script hashes, generation time, excluded-segment counts, simulation verification and limitations.
- `FINDINGS.md`: interpretation of the final regenerated diagnostics.

## Sampling and interpretation

The stationary integrated Gaussian OU displacement covariance is evaluated exactly at each scenario's frame interval. Its one-coordinate displacement variance at lag t is sigma² tau³ [t/tau + exp(-t/tau) − 1]. Adjacent increment covariance follows second differences of this variance. A Cholesky factor generates 1,024 independent 73-interval paths per model group. The script checks simulated step second moments against analytic predictions with the finite-ensemble Monte Carlo standard error; it never changes fitted parameters. The seeds are fixed, and example paths are the first draws, without selecting visually similar tracks.

The displayed physical-time scenario is the methods' 10-minute cadence (73 intervals span 12 h 10 min). CSVs also evaluate equally spaced 12-hour and 16-hour endpoint spans. None establishes the true acquisition timing; playback frame rate does not supply physical time. Spatial plots and measured per-frame displacements remain identifiable without resolving the timing.

Two centroid sensitivity calculations use 0 and 3 pixels. The reviewed extraction's 3-pixel working position allowance (0.962 µm) is **assumed to be an independent per-coordinate Gaussian standard deviation only for this sensitivity calculation**; neither independence nor a Gaussian distribution nor that standard deviation was measured. It adds diagonal increment covariance 2a² and adjacent covariance −a². No centroid noise is subtracted from the observed trajectories. Deformation, threshold changes and photobleaching could produce structured centroid errors beyond this assumption. The approximately 1.8% spatial scale reading sensitivity is reported around each observed mean step in the CSV; it is not a biological uncertainty interval.

Raw zero-model-mean step-vector correlation is the average dot product at lag k divided by the mean squared step. This avoids imposing a short-track sample-mean subtraction. The model expectation is the corresponding increment-covariance ratio; finite-trajectory ratios need not be unbiased. For the seven smaller nuclei, the displayed envelope uses independent sets of seven simulated tracks, matching the observed number; for N010 it uses individual simulated tracks. The shaded 2.5–97.5% ranges describe fixed-parameter model realizations, **not parameter confidence intervals or biological population confidence intervals**. Lags and successive steps are correlated observations, not replicates.

Whitened residuals solve L r = observed increments for the full OU-plus-assumed-localization increment covariance C = L Lᵀ. Each component is the model's forward conditional prediction error divided by its conditional standard deviation. Under the fitted model and sensitivity assumptions, these residuals should be independent standard normals. Their quantile plot, RMS, coordinate means and lag-one correlation provide diagnostic checks, without refitting or assigning a biological population p-value.

Slow and fast episodes are consecutive steps below the model radial-step 25th percentile or above its 75th percentile. These thresholds are defined by the fitted OU plus localization sensitivity and remain distinct from the paper's experimental slow-mode classifier. Longest-run lengths are compared with fixed-parameter simulated trajectories of the same 73 intervals. The representative movie cannot identify sub-frame velocity noise, establish an OU persistence time, resolve timing, or demonstrate a population WGD effect.

## Reproduction

From the repository root, after final selection of the primary fit:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/refined_ou_20261003/scripts/movie_diagnostics.py --fit wgd_langevin/refined_ou_20261003/results/selected_fits.json
```

This is a small exact-covariance computation (one CPU thread, under a minute, well below 8 GB); no SLURM allocation is needed. All source extraction tables remain unchanged.
