Current Xuan/Mallin/Jemaà results: [Gaussian OU](gaussian_literature/REPORT.md) and [PNG](gaussian_literature/results/literature_gaussian_overview.png). Hosny is unchanged. Material below is historical.

# Historical report: constant-speed model

This report is superseded by [the Poisson reanalysis](poisson_literature/REPORT.md),
which uses the same simulator as Hosny. The numerical conclusions below refer to
the old model and should not be transferred to the new parameters.

---

# Results: a Langevin explanation of WGD-associated migration

The implemented active-particle model captures two distinct routes to altered
migration: greater sampled speed in stable tetraploid clones, and slower but more
directional motion in giant polyploid populations. The current literature supports
these routes to different degrees. It does not identify a universal WGD parameter
transformation or a unique molecular cause.

All numerical fits below are conditional on the measurement definitions and
approximate figure values in [the observation table](data/observations.csv).
[Methods and reproduction commands](README.md) describe assumptions and provenance.
No new experiments or raw cell trajectories were analyzed.

## Main parameter interpretation

| System | What the model can explain | Parameter interpretation |
|---|---|---|
| RKO, Jemaà 2023, D1/T1 | About 2.6-fold higher sampled speed | A 2.6-fold propulsion increase is sufficient if persistence is fixed. A persistence-only change also fits; speed alone cannot choose. |
| Sarcoma, Jemaà 2017, Cl1/Cl8 | About 1.7-fold higher sampled speed | A 1.7-fold propulsion increase is sufficient at fixed persistence. Again, a persistence-only explanation is compatible. |
| Giant polyploid cells, Xuan 2018 | Lower sampled speed but greater late MSD | Representative homogeneous fit: propulsion 27.6 to 15.1 µm/h; persistence 0.150 to 0.974h; angular diffusion 6.68 to 1.03 per hour. |
| PACCs, Mallin 2023 | Lower accumulated travel but greater directness and net displacement | PACC propulsion about 7.1 µm/h and persistence about 2.6h; parental motion lies near the limit where turns are unresolved between 30min frames. Its intrinsic speed and persistence are weakly constrained. |

The Xuan primary fit corresponds to a propulsion ratio of 0.55 and a persistence
ratio of 6.5. The long-time dispersal coefficient is about 1.96 times larger despite
lower speed. These are parameters of a homogeneous effective population, not a
claim that every giant cell behaves alike. Compatible marginal ranges from the
reading intervals are 0.056–0.277h for reference persistence and 0.653–1.451h for
PGCC persistence. Allowing heterogeneity changes these ranges; they are not
confidence intervals.

Jemaà persistence profiles span the full search grid. The numerical minimum in a
CSV is one arbitrary member of a degenerate family, not an estimate of biological
persistence. Sampling explains the degeneracy: rapid turns within an hour lower
recorded frame-to-frame speed even when intrinsic propulsion stays the same.
Other clones in the table preserve the reported ranges (RKO approximately
2.3–2.8-fold and sarcoma approximately 1.5–1.8-fold); they were not pooled into the
primary clone contrasts.

![Observation and model comparison](results/observations_vs_model.png)

## What predicts observations excluded from fitting?

**Xuan MSD crossover.** The primary mean-speed/late-MSD fit predicts a crossover
at 1.33h, earlier than the reported approximate 2.5h (tested range 2–3h). However,
17,300 pairs on the grid satisfy the original fitting intervals and also cross
between 2 and 3h. Thus the point optimum misses this check, but the compatible
model family contains the observed crossover. This count is an existence check,
not a probability or a posterior fraction. See
[crossover sensitivity](results/crossover_sensitivity.json).

**Xuan early MSD.** At 600s the primary fit predicts approximately 15.1 µm² for
reference cells and 6.0 µm² for PGCCs. The recorded early ranges are 15–25 and
9–15 µm², respectively. The PGCC check misses at 600s. The underlying note places
these readings only approximately within 600–1000s, so this is not strong evidence
of model failure without precise curve digitization. At 1000s the predicted PGCC
MSD is approximately 16.1 µm², illustrating why the uncertain horizontal coordinate
matters.

**Xuan reported MSD exponents.** A log-log fit over 1–6h gives model slopes 1.07
and 1.45, matching the reported approximate 1.07 and 1.46. A fit over 10min–6h gives
1.23 and 1.67, which does not. The original fitting window is not in our extraction;
we therefore report window sensitivity rather than declaring a validated exponent
fit. The 1–6h agreement is a post-fit diagnostic, not an independently specified test.

**Mallin endpoint displacement.** Fitting accumulated path length and directness
predicts approximately 43 and 68 µm endpoint displacement, versus the extracted
ranges 35–45 and 65–70 µm. These are checks on other summaries from the same tracks,
not independent experimental validation.

The Mallin parental directness estimate, 0.13–0.14, falls below the model's fitted
value of about 0.142; independent larger simulations give about 0.146. Consequently
there is no primary reading-only parameter set satisfying every parental interval.
This discrepancy is small relative to finite-track variability: conditional on
50 independent tracks, the model-based sampling sensitivity admits parental tau
from the grid floor up to about 0.165h. For 49 PACC tracks, the analogous persistence
range is approximately 1.87–3.60h. These sensitivity ranges are model-based and do
not include biological-replicate dependence. The parental lower limit remains
unresolved; its point-fit intrinsic speed should not be interpreted literally.

## Does population heterogeneity change the explanation?

Yes for MSD amplitude, but not for all the observed time dependence.

An alternative with a common persistence time of 0.150h, mean propulsion changing
from 27.6 to 17.6 µm/h, and PGCC speed CV about 1.85 also fits Xuan's mean speed and
late MSD. The reference CV is fixed to zero in this comparison. Thus those two
summaries alone do not establish longer persistence.

However, when persistence is shared and Dt=0, heterogeneous constant speeds only
multiply the MSD curve by a constant. That alternative cannot produce an MSD
crossover or different MSD slopes over the same time window. The reported crossover
therefore favors a persistence change over this speed/heterogeneity-only explanation
within the tested model class. More complex correlations or switching models have
not been ruled out.

For the Jemaà fast tails, speed heterogeneity provides a plausible explanation:
illustrative CV=0.7 at fixed 2h persistence produces an RKO 99th-percentile sampled
speed around 25 µm/h while keeping the mean near 7 µm/h. Sarcoma CV=0.35 produces
an approximately 56 µm/h 99th percentile near a mean of 27 µm/h. These are scenarios,
not fitted CV estimates. The source notes do not provide the individual values or
complete sample sizes required for a defensible distributional fit. Higher CV is
not established merely because an unnormalized maximum is larger.

![Parameter identifiability](results/parameter_profiles.png)

## Scope and verification

- Five physical/numerical tests pass: ballistic limit, active-motion MSD,
  diffusion plus localization noise, heterogeneous-speed second moment, and
  invalid time-grid handling.
- Independent simulations use 8,192 cells and a different seed. Sampled speed
  predictions agree with the fitting bank to within roughly 0.4%; endpoint and
  directness differences are a few percent, consistent with the limited precision
  of a Monte Carlo bank. These errors matter for narrow reading intervals.
- The independent verification repeats integration at twice the internal temporal
  resolution. [Numerical results](results/numerical_validation.csv) retain both
  values and Monte Carlo standard errors; they do not conflate stochastic sampling
  variation with a deterministic timestep error.
- [Example simulated trajectories](results/example_trajectories.png) use the same
  spatial axes within each study. They are illustrations, not experimental tracks.
- No interactions, wound-edge effects, chemotactic gradients, cell division,
  deformability, or explicit size scaling are fitted. The existing Potts and
  polarity-model campaigns are unaffected.

No slow/moving switching extension is fitted in this first version. The reported
slow-state fractions use an unknown cell-length threshold, so they do not provide
a defensible test of the homogeneous model or identifiable switching rates. The
remaining discrepancies justify better extraction before adding unconstrained
parameters. The most informative next data would be exact MSD curves and fit
windows, framewise tracks, and per-cell speed/directness pairs. Those would test
whether persistence is genuinely altered, whether a broad fast-cell tail suffices,
and whether cells actually switch motility states.

## Files

- [Simulator](model.py) and [analysis pipeline](run_analysis.py)
- [Independent verification and alternative-model analysis](validate_analysis.py)
- [Parameter fits](results/parameter_fits.csv), [profiles](results/parameter_profiles.csv),
  [model comparisons](results/model_comparison.csv), and
  [heterogeneity comparison](results/heterogeneity_model_comparison.csv)
- [Predictive checks](results/predictive_checks.csv), [sensitivity](results/sensitivity.csv),
  [slope windows](results/slope_window_sensitivity.csv), and
  [heterogeneity scenarios](results/heterogeneity_scenarios.csv)
- [Run metadata and parameter ranges](results/summary.json)
