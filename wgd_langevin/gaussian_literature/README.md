# Gaussian velocity OU: Xuan, Mallin and Jemaà

[Overview PNG](results/literature_gaussian_overview.png) · [Results and limitations](REPORT.md)

This is the current model for the original uniform-medium literature contrasts.
Hosny's Poisson analysis remains unchanged. The preceding constant-speed and
Poisson literature outputs are historical comparisons, retained in their folders.

Model: dx = v dt; dv = -v/tau dt + eta dW.
Each axis receives independent standard Wiener noise. All cells in a condition
share tau and eta; initial velocities are stationary Gaussian. There is no drift,
permanent cell heterogeneity, constant-speed constraint, or kick-rate parameter.

Tau is in hours; eta is in micrometres/hour^(3/2). Stationary per-axis velocity SD
is eta*sqrt(tau/2), and RMS speed is eta*sqrt(tau). The latter is a derived quantity,
not a third fit parameter. The existing exact integrated Gaussian simulator in
`../hosny/ou_velocity/model.py` is imported read-only; its `sigma` argument denotes
stationary per-axis SD, NOT the Wiener coefficient eta. No Hosny files are modified.

## Fitting and observations

The original local observations and their reading intervals are reused. Xuan fits
mean sampled speed and one late-lag MSD; Mallin fits mean sampled path and mean
per-track directness. No entire MSD curve is fitted. Early MSD, crossover, endpoint,
slow fractions and window-dependent MSD slopes remain checks or diagnostics.
Jemaà's selected clone speeds define viable curves, not fitted unique parameters.
Additional clone observations remain context, as in the previous analysis.

The objective is summed squared log residual divided by half-log reading width
squared. These widths are reading tolerances, not SEM, confidence intervals or a
likelihood. At each tau, eta is profiled analytically in the joint objective.
Mean speed, path, endpoint and MSD use exact Gaussian moments at the reported
camera cadence. Directness is the mean of individual endpoint/path ratios, fitted
using a 161-point log-tau simulation bank with 16,384 cells per point, interpolated
with PCHIP. It is not approximated by a ratio of means. Final simulated checks use
65,536 fresh cells. Memory bounds are 0.0005–48 h. The bank is only for Mallin's
0.5 h cadence and 24 h duration.

Panel B jointly refits both populations under four hypotheses: shared parameters,
only tau differs, only eta differs, both differ. A statement about 'noise only'
means the Wiener coefficient eta, not stationary velocity SD or RMS speed. Such
restricted comparisons depend on parameterization.

Xuan slow fractions are evaluated using overlapping one-hour windows, with both
net displacement and sampled-path interpretations of an unknown common cell-length
threshold. No threshold is fitted and no group-specific cell lengths are asserted.
Source definitions remain ambiguous. Slopes are evaluated over three candidate
lag windows because the reported regression window is unavailable.

## Reproduction

From the repository root, in the existing scaling_sandbox Conda environment:

    python wgd_langevin/gaussian_literature/run.py
    python wgd_langevin/gaussian_literature/plot.py

A fresh run takes approximately 15 seconds on one CPU in this environment and does
not require SLURM. Only PNG is exported as a figure. CSV/JSON/NPZ files retain
parameters, predictions, joint refits, validation and source hashes. Delete or
rename the local bank before rerunning after changes to the simulator or bank
construction; its cache only checks the tau grid. Nine numerical checks test
stationary variance, integrated MSD and camera-sampled speed across three memory
scales. Independent validation additionally checks the analytical means and the
Mallin directness interpolation.
