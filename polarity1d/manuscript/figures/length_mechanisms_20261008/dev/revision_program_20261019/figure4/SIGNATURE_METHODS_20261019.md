# Figure 4 signature methods and scope

These methods describe the broad-strength experiment and its derived signature
package. Completion, generation times and source hashes are recorded in each
delivered package; this methods note does not replace those receipts. The
preceding accepted reports and the original failed Spring sampling receipt
remain preserved.

## Native experiment and input

Each condition begins from its exact qualified polarized basal. Chemical
circumferences are 20, 40 and 60, with cell width 0.05. Spring resting
circumferences are 1, 2 and 3 with 100, 200 and 300 cells; physical basal
circumferences are approximately 1.04140310, 2.04418073 and 3.04516825.
Spring's physical circumference evolves during the experiment. Its mean-force,
global-extension closure does not represent local mechanical-force propagation.

The added positive input at each position is

$$
I_f(q,t)=f I_{\mathrm{ref}}
\exp\!\left[\sigma Z(q,t)-\frac{\sigma^2}{2}\right],
\qquad \sigma=1 . \tag{1}
$$

The reference input is 0.1 for Wave-pinning, 6.4 for Goryachev, 1 for Holmes 3
and 4 for Spring. Increasing the common factor scales the entire signal,
including its ensemble mean and standard deviation. It is not a variance-only
experiment. Uniform controls receive exactly the declared ensemble mean
$fI_{\mathrm{ref}}$, without normalizing individual realized paths; no-added-input
controls are retained.

All models use the common factors 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03,
0.1, 0.3 and 1. The five middle factors have twelve independent noise seeds in
the replicated package; the four outer factors have two exploratory seeds.
The complete design comprises 816 noisy cases, 108 deterministic mean controls
and twelve no-added-input controls. Exact original cases are reused only with
verified source, basal, seed and forcing-protocol identity.

The Gaussian field is spatially homogeneous, with stationary OU knots and
analytic variance-normalized linear interpolation between knots. The spatial
kernel is the periodized Gaussian, discretely normalized in L2 and applied by
circular convolution to independent spatial Gaussian innovations; its exact
covariance is the kernel autocorrelation. The spatial zero mode is retained.
Normalized interpolation preserves unit marginal Gaussian variance, but its
two-time covariance can depend on phase within a knot interval. This is a
declared finite-knot law, rather than an assertion of exact continuous-time OU
stationarity.

| Quantity | Chemical rings | Spring ring |
|---|---:|---:|
| Native horizon | 32000 | 3200 |
| Spatial kernel length | 1.5 physical units | 0.15 basal physical units |
| OU knot relaxation time | 50 | 5 |
| Applied knot interval | 5 | 0.125 |
| Latent coupling interval | 2.5 | 0.125 |
| Scalar output cadence | 5 | 0.5 |
| Native profile cadence | 40 | 4 |
| Full native-state cadence | 1000 | 100 |

Spring forcing is fixed in material coordinates: its physical correlation
length scales as $0.15C(t)/C_{\mathrm{basal}}$. These native units and forcing
references are not a biological calibration across models. Solvers, transport,
pool constraints and source hashes are frozen in the
[pilot manifest](../../../../../data/runs/revision_noise_strength_broad_20261019/manifest.json)
and [scoped replication manifest](../../../../../data/runs/revision_noise_strength_replication_scoped_20261019/manifest.json).

## Physical activity vector and its separate readouts

For active concentration $u_j$, material fraction $q_j$ and physical cell width
$\Delta s_j(t)=C(t)/N$, the same moments used in Figure 3 are

$$
M_n(t)=\sum_j u_j(t)\,\Delta s_j(t)\,
e^{2\pi i n q_j},\qquad
P(t)=\frac{M_1(t)}{M_1(0)} . \tag{2}
$$

$P$ is the normalized physical-activity vector: division rotates its initial
axis to zero and scales by its initial first-moment magnitude. It is not
normalized independently by each frame's active amount. The physical active
amount $M_0$ is distinct from the conserved total of active and inactive
species. Actual Spring $C(t)$ enters every moment. The magnitude obeys the
paired-frame identity

$$
R(t)=|P(t)|=A(t)Q(t),\qquad
A(t)=\frac{M_0(t)}{M_0(0)},\qquad
Q(t)=\frac{q_1(t)}{q_1(0)},\quad
q_1(t)=\frac{|M_1(t)|}{M_0(t)} . \tag{3}
$$

The mean of the product is retained; separate means of $A$ and $Q$ must not be
multiplied to reconstruct mean $R$. Raw basal $M_0$, $|M_1|$, $q_1$ and $q_2$
remain recoverable in the context table. A weak first moment can arise from
cancellation, low coherence, amount depletion or other spatial changes.

Directional memory is the pooled cosine of endpoint heading differences,
conditional on $R\geq0.1$ at both endpoints. Its usable-pair coverage divides
eligible pairs by all available pairs. Groups below 100 eligible pairs or 1%
coverage remain unavailable. The main fixed lag is 8000 chemical or 800 Spring
units; the supplements retain all six native lags and full/half windows.
The uncentered full-vector correlation is a separate descriptive quantity
that includes magnitude and forced-state adaptation; it is not a centered
stationary-process fit or migration autocorrelation.

## Spatial organization, weighting and observation windows

A threshold region is a periodic connected arc of cells at or above the current
minimum plus half the current contrast. Flat contrast at most $10^{-12}$ gives
zero regions. Raw multiple-region occupancy retains all scalar frames. Its
amplitude-qualified counterpart requires current contrast at least 10% of the
same case's initial basal contrast. Low-contrast occupancy is reported
separately. Absolute prominence 0.05 is a native concentration threshold, not
5% of each field's range; a separate 5%-of-basal-contrast sensitivity is derived
from saved profiles. These counts do not assign competition or exclusive phases.

Physical widths join seam-crossing regions before counting cells and multiply
the cell fraction by actual $C(t)$. The largest width contributes once per
amplitude-qualified profile frame; the second-largest contributes once per
qualified multiple-region frame. Condition means weight available seed-level
descriptors equally. Frame coverage, region counts, contributing seeds and
material-fraction widths accompany conditional physical widths. Widths are
limited by the native cell mesh and actual profile cadence; they do not measure
continuous lifetimes or site lineage.

Main summaries discard the first half, using chemical times 16000–32000 and
Spring times 1600–3200. Full, first-half, second-half and quarter windows remain
available. Inclusive adjacent windows share their boundary sample. Late-quarter
differences compare paired seed descriptors; their shorter heading lag is
2000 chemical or 200 Spring units. Window sensitivity is retained rather than
assuming equilibrium or certifying stationarity.

Noise-condition intervals use 2000 whole-seed percentile resamples, with the
same bags paired across readouts, windows and lags. Direction is recomputed
from resampled sufficient-statistic sums, with coverage and support checked
within each bag. Intervals remain unavailable when fewer than 95% of bags have
support. Two-seed outer strengths retain descriptive seed ranges; deterministic
controls are not assigned bootstrap replicates. The intervals condition on
the exact basals, model parameters, forcing family and observed windows.

## Numerical qualification and biological interpretation

Four largest-size, factor-1 paths were compared using identical finite forcing
under standard/tighter ODE solves and denser saved observations. The original
Spring check failed twelve sampled multiple-region run-q90 observation gates;
the corresponding ODE and selected occupancy, vector and memory checks passed.
The [original passed:false receipt](../supporting/metadata/FIGURE4_BROAD_NUMERICAL_SPRING_RING_INDEPENDENT_CHECK.json)
is unchanged. A separately reviewed and
[explicitly accepted ten-metric scope](../../../../../data/runs/revision_noise_strength_replication_scoped_20261019/results/STABLE_ESTIMAND_ACCEPTANCE.json)
excludes both run-q90 descriptors. This does not relabel the full original gate
as passed. All run spans remain cadence-dependent sampled descriptors; no event
rate, establishment or continuous lifetime is inferred. Instantaneous widths
and endpoint-concurrency geometry remain exploratory derived descriptions
outside that numerical acceptance.

Goryachev's C20, factor-0.01 mean control already loses its first moment, and
larger high-input mean controls often become low contrast. Such conditions are
forced-state changes rather than spontaneous wandering about an unchanged
incumbent. Across models, a finite strength grid locates observed changes over
sampled intervals, not common calibrated doses or phase boundaries.

Activity-heading memory provides a measurable activity-persistence analogue to
the empirical WGD migration motivation. It is distinct from cell movement.
Migration predictions require an independently measured activity-to-motion
coupling; multiple activity regions alone do not establish functional
competition. Folded physical endpoint displacement and angular-remapping
evidence from the original strengths remains accessible, with its existing
shortest-arc, coverage and winding limitations.

## Separate targeted Goryachev refinement

The broad grid left the C60 interval between factors 0.01 and 0.03 unresolved.
A separately approved refinement added factors 0.015 and 0.02, each with twelve
independent noise realizations and one matched-mean control: 26 fresh cases.
The exact basal, mesh, native law, correlations, solver, horizon and cadences
were unchanged. Its [frozen manifest](../../../../../data/runs/revision_gory_strength_transition_20261019/manifest.json)
and [two-page supplement](gory_transition_refinement_20261019/gory_transition_review.pdf)
remain separate from the immutable 936-case broad-grid closure.

Both added levels retained qualified contrast and a coherent heading in the
second half; no contrast-qualified multiple-region scalar sample was observed
there. Mean contrast relative to basal was 0.7765 and 0.6459, respectively, and
lag-8000 heading cosines were 0.9991 and 0.9966 with full usable-pair coverage.
At the retained factor-0.03 anchor, low-contrast occupancy was 0.8560. These
sampled points narrow the unresolved interval to (0.02, 0.03); they do not
exclude an intervening organization or establish a phase boundary. No further
strength ladder was pursued. Exact compact fields and CI tables support a
native-archive-free redraw of this separate supplement.
