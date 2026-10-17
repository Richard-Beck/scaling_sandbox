> Active campaign (2026-09-26): [direct WGD fits with persistence-only differences](../../persistence_only/results/report.html). OU-to-polarity fitting and transfer have been abandoned. The material below is historical.

# Revised joint WGD refits

This replaces the frozen-parameter literature screen, not the spatial simulations or their calibrations. Open `results/report.html` for completed WGD fits and `results/spatial.html` for experiment-level heatmaps and grouped trajectory inspection.

## Transferred information

For each spatial source × OU candidate, select the smallest calibration objective among the tested **nonconstant** laws. The constant law remains a null. Record the calibration improvement, including negative improvements caused by imperfect optimization; do not treat forced selection as evidence for scaling. Only Otsuji's nonlinear candidates passed the previous 10% calibration-improvement threshold.

The selected law determines which parameters may differ between reference and polyploid/PACC groups. WGD fitting reoptimizes all identifiable parameters, independently for each study. No numerical prior or direction constraint is imposed. Fits with identical family and eligible parameter sets are the same inference problem and are computed once. Length-law signs are compared afterward, assuming larger effective length in the polyploid/PACC group. Ratios alone are exploratory: parameter identifiability and uncertainty have not been established.

All other parameters, including noise amplitude and readout scale, are shared. Gaussian and linear-memory instantaneous noise is set to one as a units convention; fitting a shared readout scale exactly restores a free shared physical noise amplitude. Nonlinear ps=1 removes an exact latent-state scale redundancy; noise in latent units and physical readout scale are both fitted. Linear-memory noise uses the published-style substitute covariance closure already declared by the original bridge, rather than independent noise in both states. Nonlinear models retain the declared isotropic two-vector extension, noise in p only. These are phenomenological choices, not measured spatial residual noise.

The fixed cue in Hosny identifies only gain × cue amplitude. In unforced observations gain is unidentifiable and is excluded. No knockdown mechanism is introduced. Source power-law magnitudes, absolute lengths, time conversion and parameter constants are not frozen across systems.

## Fitting and adequacy

Four bounded least-squares starts plus higher-sample polishing. Null fits run first and initialize one differential start. All optimized coordinates are log-positive. Parameter bounds are stored per result. Common random numbers stabilize stochastic fitting; both WGD groups use exactly the same random draws, ensuring a shared null produces exactly equal population predictions.

Xuan/Mallin/Jemaà use log residuals divided by half the log figure-reading interval width; a shared physical amplitude is profiled analytically. Hosny fits available control quartiles for path length, endpoint, directness and the two directional ratios, dividing residuals by each observed IQR (missing IQR uses 0.1). These are descriptive objectives, not likelihoods. Hosny summaries from the same tracks are correlated; repetitions across control blocks are not independent evidence.

Report fit values and independent-seed validation values. Nonlinear validation also halves dt from 0.02 h to 0.01 h and doubles burn-in from 48 h to 96 h. Search/polish/validation ensembles are 64/256/1024 independent cells (linear 256/1024/4096, with exact speed and MSD where unforced). The simulation was checked against the exact OU limit: MSD relative error 4.45% with 2048 independent cells; identical-group/common-random-number test passes exactly.

All reading intervals must be met by the fitted and validation predictions before the report calls a fit adequate. Separately, the JSON records whether an amplitude exists satisfying all intervals, by intersecting the amplitude intervals analytically. This avoids confusing an optimized-amplitude miss with infeasibility. A failed bounded search does not prove impossibility. Jemaà speed-only results are explicitly underdetermined. Hosny has no binary statistical acceptance criterion; inspect its quartile residuals.

No confidence intervals, bootstrap or sign-identifiability profile is claimed in this first pass. Finite sample/step/burn checks expose simulation sensitivity, not optimization completeness. Existing source calibration nonidentifiability remains relevant.

## Spatial reporting

The report now leads with one condensed overview: spatial sources × OU candidates, comparing fitted experiments at fitted lengths, new experiments at fitted lengths, and all experiments at unseen lengths. Independent, constant and power-law fits are separate columns. Every heatmap uses post-event normalized RMSE, aggregated by RMS over the included experiments and lengths. Two collapsed detail figures preserve experiment-level fitting and unseen-length errors. The redundant optimizer-objective plot is removed; its data remain in CSV.

Trajectory inspection shows ten physical contact sheets: for each of five spatial sources, a 20-facet calibration image and a 15-facet held-out-experiment image. Experiments are rows; five lengths are columns. All four OU candidates and the spatial trajectory appear together. Three global fit buttons switch every sheet between independent, constant and power-law fits; there are no trajectory dropdowns. Independent predictions at unseen lengths are explicitly unavailable. Row scales are shared across lengths; the plots focus on the cue-change response after conditioning, while pulse and competition retain their complete trajectories. Each sheet links to a full-resolution PNG.

## Execution

SLURM arrays 20596643 (20 nulls) and 20596644 (36 differential fits; depends on successful null completion). One CPU and 3 GB per task, red/small, one-hour wall limit, no concurrency throttle. Each completed result atomically republishes the current report under a filesystem lock. No browser refresh timer. Manifest and result files are retained under `results/`.

Run with the existing `scaling_sandbox` Conda environment:

    python campaign.py --manifest
    python campaign.py INDEX
    python spatial_report.py
    python publish.py
