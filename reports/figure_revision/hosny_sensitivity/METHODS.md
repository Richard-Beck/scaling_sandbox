# Hosny sensitivity analysis

`analyze.py` writes the computational support files to
`reports/figure_revision/figures/figure3_review/hosny_sensitivity/`. Run it
from the repository root in the existing environment:

```sh
conda run -n scaling_sandbox python reports/figure_revision/hosny_sensitivity/analyze.py
```

The full relative surface is deliberately broad for context. The dense local
contour grid used for the primary panels is generated separately, without
re-running the profiles or leaveout fit:

```sh
conda run -n scaling_sandbox python reports/figure_revision/hosny_sensitivity/refine_surface.py
```

It evaluates a 61 by 61 logarithmic grid spanning tau ratios 10–21 and noise
ratios 0.10–0.19, inserting the exact selected-fit coordinate. Its CSV repeats
the same fixed-drift and two-drift-refit columns as the broad surface. The JSON
records input hashes, exact dimensions, convergence/bounds checks, and
sampled Δ1/4/9 coordinate extents for both the broad and refined grids.

The selected primary fit is `fit_05.json`. Its data score and regularized
objective are recomputed with the exact fitting residual mask (Hosny medians
only) and checked against the stored values before any sensitivity analysis.
The fit objective separates data misfit, pooling penalty, and persistence
prior penalty.

The relative tau/noise surface fixes WGD− and preserves the five fitted WGD+
study deviations while sweeping group geometric-mean ratios over tau 0.5–30
and noise 0.05–2, clipped to study-level model bounds. It compares the stored
Hosny drifts with a bounded reoptimization of both drifts at every point. The
same data residual mask and primary regularization terms are used. Because
drifts are unpenalized, reoptimization should not increase the total objective.

The three drift profiles hold in turn the reference drift, WGD drift, or their
contrast fixed. At each sampled value, all 20 log parameters and the remaining
drift coordinates are reoptimized. A coarse pilot identifies the interval
within the descriptive Δ objective 9 threshold; a dense 31-point grid covers
that interval, with zero and the unrestricted fitted value included exactly.
Each point uses three starts: fit_05, the preceding profile solution, and a
deterministic perturbation. Exported ranges are sampled profile sensitivity,
not joint confidence intervals or full uncertainty limits. Objective changes
of 1, 4, and 9 are descriptive thresholds selected before profiling.

The leave-Hosny-out fit removes Hosny from data residuals, pooling centers and
penalties, and persistence penalties, then refits the remaining 16 study
parameters using the same primary shrinkage scales. The comparison CSV also
reports the full five-study pooled folds and the same four studies evaluated
at their full-fit parameters and pooled only over those four.

The analysis locally replaces SciPy's overflow-prone Rice mean calculation
with the equivalent scaled-Bessel expression. For noncentrality above 10,000,
it uses the large-noncentrality normal asymptotic for the Rice mean and
quantiles. This avoids non-finite optimizer trials while agreeing with SciPy
at ordinary noncentralities; the exact fit score and selected multistart
fits are validated. Per-start parameter vectors and loss decompositions are
retained in the CSVs, with input hashes and convergence checks in
`hosny_analysis.json`.
