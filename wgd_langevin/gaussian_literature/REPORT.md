# Gaussian OU results for the older literature

[Standalone overview PNG](results/literature_gaussian_overview.png)

All fits are homogeneous, stationary and unbiased. Hosny is unchanged.
The original two fitting summaries per population are retained.

| Study / population | Memory tau (h) | Wiener coefficient eta (µm/h^1.5) | Derived RMS speed (µm/h) |
|---|---:|---:|---:|
| Xuan reference | 0.1067 | 99.91 | 32.63 |
| Xuan polyploid | 0.7255 | 20.29 | 17.28 |
| Mallin reference | 0.0005, lower bound | 14347 | 320.81 |
| Mallin polyploid | 2.120 | 5.538 | 8.062 |

**The Mallin reference row is an unresolved boundary fit.** Its instantaneous
noise and speed values are not biological estimates. As tau falls below the
camera interval, eta can increase to retain the observed step sizes, approaching
position diffusion. Lowering the bound cannot meaningfully solve the directness
mismatch; it mainly changes the reported instantaneous parameters. Two parameters
for two targets does not guarantee that those targets are jointly attainable.

## Fit and checks

Xuan's two fitted central values are reproduced for both populations. The fitted
polyploid memory is about 6.8 times reference, and RMS speed about 0.53 times
reference, conditional on this model and the supplied figure-read summaries.
No confidence intervals or biological significance are claimed.

Mallin path means are reproduced. Reference directness is approximately 0.1443
in independent simulation, versus the reported 0.130–0.140 interval. Polyploid
directness is approximately 0.4098, within 0.400–0.420. The two endpoint checks
are now within their supplied intervals: 44.04 versus 35–45 µm, and 68.82 versus
65–70 µm. Endpoint and fitted summaries come from the same tracks, so this is a
held-out summary check rather than independent experimental validation.

Xuan's predicted early MSDs at the nominal 600-second lag are 18.71 and 7.69 µm².
The reference meets its reading interval; the polyploid misses 9–15 µm². The source
abscissa is uncertain (approximately 600–1000 seconds). The predicted crossover
is 1.05 h, earlier than the reported 2–3 h. Those observations were not optimized
against. Their mismatch remains visible rather than being absorbed by extra fit
parameters. Slow-fraction threshold curves and candidate-window slopes are retained;
unknown numerical cell length and observation definitions limit interpretation.

In the joint refits, allowing both tau and eta to differ gives the best score for
both studies. Mallin retains residual error even with both free. This comparison
is specifically in tau/eta coordinates; 'tau only' holds the Wiener coefficient
fixed, not the stationary RMS speed. Scores are descriptive discrepancies, not
likelihood ratios. Jemaà speed-only observations still cannot distinguish memory
from noise intensity, and are shown as constraint curves.

## Figure guide

A: optimal tau/RMS points and regions meeting both fitting intervals. No reference
Mallin compatible region exists on the grid. B: joint restricted refits. C: tau,
eta and derived RMS values, with the Mallin reference boundary flagged. D: fitting
targets and Mallin endpoint checks. E: Xuan MSD and crossover check. F: slow-fraction
sensitivity to a common unknown threshold (solid net displacement, dashed sampled
path). G: Jemaà speed-only curves, with camera cadence included. H: simulated tracks,
with short realizations taken as prefixes of long realizations.

Nine analytical simulator checks pass. Independent simulations agree with exact
mean predictions; Mallin bank-to-validation directness differences are about
0.00029 and -0.00021, below one validation standard error. Reducing parameter count
has removed the unconstrained Poisson rate but has not removed observational or
model-specification limitations.
