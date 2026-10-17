# Revised older-literature results: the common Poisson model

The earlier Xuan, Mallin and Jemaà work now uses the **same simulator and physical
parameter definitions as Hosny**, with isotropic kicks (kappa=0) for the original
uniform-medium assays. The old constant-speed results remain explicitly archived.
[Overview PNG](results/literature_poisson_overview.png).

The important change is interpretive: intermittency can change sampled mean speed
and track directness without a permanent heterogeneous population. The original
small sets of fitting summaries no longer identify a unique memory change.

## What the comparisons establish

For Xuan, allowing either tau+rate, tau+amplitude, or rate+amplitude to differ can
match mean speed and the single late-MSD target essentially exactly. A shared-memory
rate+amplitude explanation therefore exists for **those fitting targets** even
with homogeneous cells. Single-parameter changes alone are inadequate.

For Mallin, all three two-parameter alternatives also match the fitting targets as
well as the fully flexible model, including a shared-memory rate+amplitude model.
However, every tested version retains the parental-directness discrepancy described
below. “As well as full” must not be read as “within every reported interval.”

Xuan's reported MSD crossover remains informative: at a common tau, changing rate
and amplitude only rescales the entire stationary MSD curve, so it cannot cross.
Thus the qualitative crossover supports a memory difference within this model,
while the fitted speed/late-MSD amplitudes alone do not. This is narrower than the
previous blanket interpretation of longer persistence in both studies.

Jemaà remains speed-only compatible families. At fixed tau and rate, the sampled
speed ratios are explained by matching kick-amplitude ratios (about 2.64 for RKO
D1/T1, 1.69 for sarcoma Cl1/Cl8). Neither memory nor rate is identified. The plotted
contours are one rate slice, clearly labelled, not estimates of that rate.

## Illustrative numerical minima, not identified constants

| Study/population | tau (h) | rate (/h) | kick scale A (um/h) |
|---|---:|---:|---:|
| Xuan reference | 0.0525 | 6.13 | 81.8 |
| Xuan PGCC | 0.724 | 25.0 | 4.06 |
| Mallin parental | 0.0101 | 38.4 | 116 |
| Mallin PACC | 0.217 | 0.325 | 83.5 |

These are one numerical optimum per population in an underdetermined family. The
three unknowns are typically constrained by only two fitted means. They must not
be used to assert identified fold changes or compared literally to old intrinsic
propulsion speeds. Compatible tau/rate regions and profiled rate fits are saved.
Large kick scales and short memory also reflect sub-frame compensation.

## What the independent predictions miss

Independent simulations use 65,536 cells per fitted state, separate from the
2,048-cell fitting banks. The MSD values are computed from the exact stationary
moment formula. Independent predictions give:

- Xuan mean speeds about 22.72 and 14.67 um/h, within the original reading ranges.
  The single late MSD targets are exactly matched by the illustrative parameter
  minima, but this does not represent a full-curve MSD fit.
- Xuan early MSD at 600 s is about 26.3 and 7.7 um², outside the nominal 15–25 and
  9–15 ranges. The source only places these readings around 600–1000 s.
- The illustrative Xuan curves cross at about 1.15 h, rather than the reported
  approximate 2–3 h interval. This is a failed check of the chosen minima, not a
  proof that every compatible Poisson family fails that observation.
- Mallin parental path is about 304 um and mean directness 0.1466; the directness
  remains above the extracted 0.13–0.14 interval. No parental grid point meets all
  original fitting intervals. The interval reflects reading precision, not an
  empirical uncertainty interval for a sample of cells.
- Mallin PACC path is about 162.6 um and directness 0.4058, within fitting intervals.
  Its held-out mean endpoint is about 62.1 um, below the 65–70 interval. Parental
  endpoint is about 44.6 um, within 35–45.

The endpoint failure may distinguish members of the compatible family, but it was
not silently promoted into the fitting loss. Narrow reading intervals and a
Monte Carlo interpolation bank also limit numerical precision: sampled mean speed
or path shifts by up to about 1.4% between bank and independent simulation; the
Mallin parental directness shifts by about 0.0025. Bank scores are descriptive
optimization distances, not likelihoods or uncertainty estimates.

## Xuan slow/fast observation

The supplied notes specify one average cell length in one hour, with no numerical
length and unresolved exact implementation. The old extraction's assertion that
thresholds were group-dependent was unsupported; the derived input copy corrects
it. The revision checks a **common unknown threshold** under two interpretations.

For the illustrative fit, calibrating a net-displacement threshold to reference
slow fraction 0.45 gives approximately 10.5 um and predicts PGCC slow fraction
0.433, below the reported roughly 0.695. For one-hour sampled path, the corresponding
threshold is about 21 um and PGCC slow fraction about 0.861, above the observation.
Neither interpretation has a common threshold meeting both reported fraction
intervals on the displayed 0–60 um grid for these point fits. These checks do not
rule out all compatible parameter families or identify switching rates. Occupancy
is computed from overlapping one-hour windows, matching a fraction of time under
a stationary homogeneous process rather than inventing permanent fast/slow classes.

Reported MSD exponents remain window-sensitive diagnostics; numerical slopes for
several windows are retained in `slope_windows.csv`, not treated as measured tau.

## Scope and verification

The same nine physical/numerical simulator checks pass. The simulator is imported
directly from the Hosny Poisson implementation, not copied or approximated by an
active-Brownian surrogate. The observation process matches each original cadence
and uses mean per-track ratios for directness. Means keep the original reading-
interval objective; Hosny's quartile loss is not imposed on unavailable quartiles.

All eight joint parameter-sharing alternatives are retained for Xuan and Mallin.
Jemaà only contributes speed contours as before. No new literature, gradient/KD
assays, cell interactions, localization noise or permanent heterogeneity are added.
Old results and their methods remain clearly marked as historical, with current
entry points redirected here and to Hosny.

The initial run completed on SLURM job 20518799, red/small QoS, four CPUs, in about
241 seconds. Source/input hashes, model path, bounds and independent numerical
checks are retained alongside the PNG and CSV/JSON outputs.
