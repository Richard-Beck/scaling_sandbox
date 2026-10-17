# Holmes Model 3 provenance and validation

This record documents the two-GTPase Holmes adapter and its focused numerical
checks. It is not a substitute for the complete seven-model screen.

## Source model and defaults

The adapter implements the explicit no-phosphoinositide Rac–Rho subsystem in
Holmes et al. (2012), “Modelling Cell Polarization Driven by Synthetic
Spatially Graded Rac Activation,” PLoS Computational Biology 8(6):e1002366,
DOI [10.1371/journal.pcbi.1002366](https://doi.org/10.1371/journal.pcbi.1002366).
The source is the Rac–Rho mutual-inhibition Model 3 section: Eqs. (1), (5),
and (6), with relevant values from Table 1. The model is implemented in
[`polarity1d/src/polarity1d/holmes.py`](../../polarity1d/src/polarity1d/holmes.py)
and registered as `holmes_model3`.

The active Rac and Rho fields each exchange material with a composite inactive
pool. Rac activation is suppressed by Rho through a Hill function; Rho
activation is suppressed by Rac. Each active/inactive pair has an equal and
opposite reaction transfer, and the no-flux diffusion solver conserves each
pair over the domain.

| Parameter | Default | Source meaning |
| --- | ---: | --- |
| Total Rac, total Rho | 7.5, 3.1 μM | Table 1 total levels |
| Rac/Rho basal activation | 0.2, 6.6 μM/s | Î_R1, Î_r |
| Rho inhibition level for Rac / Rac inhibition level for Rho | 1.25, 1 μM | a1, a2 |
| Hill exponent | 3 | n |
| Active-state decay | 1 s⁻¹ for each | d_R, d_r |
| Active / composite-inactive diffusion | 0.1 / 50 μm²/s | D_m and D_mc(L0) |
| Reference length | 20 μm | L0 |

This comparison uses the two-GTPase Model 3 subsystem without PI feedback.
Cdc42 and PI-related parameters in Table 1 do not enter this model. The
published table omits membrane-cycling parameters (kon, Dc, V, and koff),
stating that their effects are represented by the composite length factor and
effective inactive diffusivity. Since those omitted values are needed to
calculate their length dependence, this first pass holds c-tilde(L)=1 and
D_mc(L)=50 μm²/s at every length. That explicit approximation means these
Holmes results screen biochemical length effects only through geometry and
cue profiles; they do not reproduce the source's potential length-dependent
membrane cycling.

The initial condition is the low-activity homogeneous fixed point of the
no-cue two-species equations: active Rac 0.034298 μM and active Rho 2.109251
μM, with inactive concentrations set by the tabulated totals. Time is seconds,
length is micrometers, concentrations are micromolar, and the external
stimulus is an additive Rac activation rate in μM/s.

## Numerical checks

- The homogeneous default state has a no-cue reaction residual below 9e-16.
  Both active/inactive reaction pairs sum to zero to floating-point precision.
- A seeded 0.001 μM spatial perturbation on 41 cells, evolved for 0.2 s with
  BDF, remained finite and nonnegative. Relative integrated-pool errors were
  2.22e-16 for Rac and 1.44e-15 for Rho.
- Initial periodic smoke runs at 41 and 81 cells gave matching latencies under
  an earlier scorer. The full suite now records cue-following latency separately
  from whether the previous orientation was established, and applies the
  0.05 active-contrast floor at the response endpoint. Those preliminary
  refinement values are superseded and should not be used as suite results.
- The 201-cell SLURM check at L=20 μm completed the gradient-reversal task in
  7 s: eligible, with a 15 s response latency from cross-fade onset. The
  full-screen period-following results at L=20 μm have median cue-following
  latency 10 s at P=100 s and 37.5 s at P=400 s, with 12/12 switches reaching
  orientation above 0.05 and active contrast above 0.05. Periodic checks
  continue through 6.25P so the last scheduled switch has its complete
  response window.

Across the full length sweep, Holmes passes all 12 periodic switch windows at
every tested size. At P=100 s, latency is 12.5 s at L=10 μm and 10 s from
L=15–50 μm; at P=400 s, it is 37.5 s at all sizes. The direct gradient
reversal time remains 17.5 s at all seven lengths, which is the sampling
resolution set by the 2.5-s output interval. Competition and finite-countercue
cases are ineligible at every tested length under the suite's stated
eligibility checks: two competing peaks do not both form in the former, and a
stable no-cue prior polarity is not retained in the latter.

The periodic score retains per-switch pre-direction, signal-strength, and
alignment diagnostics. At the period-400 center-length check, orientation
immediately before each zero-crossing was approximately -0.003 and active-field
contrast was about 0.0015. This cue-following result therefore starts from a
nearly neutral pattern, rather than a strongly retained front. The endpoint is
the measurable response to the new cue direction; it is distinct from the
suite's direct and finite reversal tests, which require established prior
polarity.

## LEGI persistence check

The finite-countercue result for LEGI at L=9 in the suite is consistent with
its native equations. The LEGI response field has no direct decay term; after
activator and inhibitor clear, the response profile can retain spatial memory.
For the 201-cell screen, orientation was -0.1540 after conditioning and
-0.1544 after the 2000-unit no-cue interval. A 20-unit opposing cue changed
orientation to +0.0732; after the 5000-unit recovery it was +0.0755. The final
response ranged from 0.311 to 0.559 (contrast 0.247), so this outcome is not
caused by normalizing a vanishing response field. In the final periodic screen,
LEGI follows at L=9 μm for P=100 s (20 s latency) and at L=9–21 μm for P=400 s
(47.5–95 s latency); the remaining sizes are censored. This is consistent
with its native response-memory term and the length-dependent adaptation time.

The final full-array reference outputs are
[`task_0171.json`](data/task_0171.json) (gradient reversal),
[`task_0173.json`](data/task_0173.json) (period 100 s), and
[`task_0174.json`](data/task_0174.json) (period 400 s). The dedicated
201-cell gradient-reversal check ran as job 20618939_171 on one CPU.
