# Spring ring noise: definitions and bounded interpretation

Spring participates in every Figure 4 analysis family: native field histories,
initial-axis memory, physical first-moment magnitude, raw-vector and
magnitude-gated direction correlations, full/first/last-half restricted memory,
endpoint site shifts, pair coverage, geometric remapping, local-region and
weak-resultant occupancy, sampled increment distributions and deterministic
mean/no-input controls. Its native units and declared forcing differ from the
chemical models; the analyses provide within-model size comparisons.

## Exact preparations and the declared extension

Resting circumferences 1, 2 and 3 use 100, 200 and 300 uniform material cells.
Their exact prepared physical circumferences are 1.0414031003, 2.0441807295 and
3.0451682524. The middle state is the retained independently audited preparation.
The new small/large states start from uniform native concentrations, receive
the same peak-4, full-physical-support-0.2 preparation for 10 native time units,
then qualify on a 13000-unit cue-free continuation. The independent preparation
audit checks two consecutive 1000-unit windows and a 10000-unit holdout.
The original full states and hashes are preserved.

The approved carryover retains native reaction, dilution and diffusion under
uniform material stretching. The replacement of two endpoint forces by the
mean of the nonlinear local Rac-force law is the declared ring closure.
It couples periodic chemistry to one global circumference; it supplies no
local force-propagation distance, elasticity field or velocity law for a cell.

## Imposed stochastic field and numerical evidence

The positive additive native activation is

$$
I(q,t)=4f\exp\!\left[Z(q,t)-\frac12\right],
\qquad f\in\{0.001,0.01\}. \tag{1}
$$

The unit-variance Gaussian field is constructed by circularly convolving
independent Gaussian spatial innovations with a periodized Gaussian filter.
The filter is normalized by its discrete Euclidean norm, so the exact discrete
covariance is its circular autocorrelation. The spatial zero mode is retained.
The covariance-width parameter is 0.15 in each preparation's basal physical
coordinates, resolving approximately 15 native cells. No realized spatial
mean removal, amplitude rescaling or chemical clipping is introduced.

The final input uses exact stationary temporal OU knots with correlation time
5 and spacing 0.125 native units. Between knots, linear interpolation is
divided by its analytic marginal standard deviation:

$$
Z(q,t)=
\frac{(1-a)Z_j(q)+aZ_{j+1}(q)}
{\sqrt{(1-a)^2+a^2+2a(1-a)e^{-h/5}}},
\qquad a=\frac{t-jh}{h},\quad h=0.125. \tag{2}
$$

This is a specified finite interpolated process, rather than an exact
continuous-time OU path. Spatial homogeneity and stationary OU knots do not
imply exact continuous-time stationarity: two-time covariance between knots
can depend on within-knot phase despite the normalized marginal variance.
Time 5 is one tenth of the coefficient-only
circumference relaxation time 50 at fixed force; this is a native-model scale choice, not a biological
calibration or equal exposure across models.

The imposed field remains attached to material coordinates. Its instantaneous
physical covariance width therefore changes as

$$
\ell(t)=0.15\,\frac{C(t)}{C_{\mathrm{basal}}}. \tag{3}
$$

Fixed local basal statistics do not equate whole-ring input amounts or angular
forcing modes across sizes. The realized-input table separately retains exact
forcing-regeneration hashes, local moments, physical-width ranges and basal
physical integral/angular-mode moments.

The original 0.5/0.25 pilot is preserved, including its failed pointwise
vector-difference gate. A subsequent four-pair check uses two seeds at each
strength and exact conditional OU bridges that preserve every original
0.25-unit knot while inserting 0.125-unit midpoints. Before those runs,
the estimator sensitivity tolerances were frozen: at lags 0–100, vector and
direction correlations differ by at most 0.05; weak-R and multiregion
occupancies by 0.02; strict endpoint coverage by 0.05; and reference-arc RMS
by the larger of 0.01 native length or 10% of baseline RMS. Every pair passes
separately. Opposite changes are not averaged to obtain a pass.

The largest paired correlation difference is about 0.0152, whereas the largest
pointwise vector difference is about 0.193 and the largest local field
difference is about 21.9% of initial contrast. The maximum relative full-field
mean absolute difference in that pair is 3.50%. These remain recorded
finite-forcing-law path sensitivities. The separate same-law tight-ODE check
has maximum vector difference about 1.17×10⁻⁶. Neither kind of finite pilot
establishes long-horizon or cross-size stochastic convergence. The final
ensemble uses the declared finer 0.125 law without changing strengths,
correlation time, spatial scale or horizon to seek a preferred size trend.

## Physical moments and memory

For material positions \(q_j\), native active concentration \(u_j\), and current
physical circumference \(C(t)\), the physical moments are

$$
M_k(t)=\frac{C(t)}{n}
\sum_{j=1}^{n}u_j(t)e^{2\pi i k q_j},
\qquad
P(t)=\frac{M_1(t)}{M_1(0)},\quad R(t)=|P(t)|,\quad
q_k(t)=\frac{|M_k(t)|}{M_0(t)}. \tag{4}
$$

The Spring convenience helper's shape-normalized vector is not substituted
for this physical first-moment normalization. Its direction agrees, but its
magnitude generally differs when active amount or circumference changes.
Every saved primary profile independently reproduces physical moments,
current half-contrast-region counts and prominence-based peak counts.

The full-vector correlation is the uncentered ratio of pooled within-seed
endpoint products to the square root of pooled endpoint energies. It includes
weak resultants, with weak endpoint RMS flagged separately. Direction
correlations use the unit activity vector only when both endpoint magnitudes
meet the declared gate 0.05, 0.1 or 0.2. Neither quantity is a cell-velocity
correlation or a fitted migration-persistence time.

The experiment spans 3200 native time units; scalar cadence is 0.5 and primary
field cadence is 4. All full, first-half and last-half analyses use the same
maximum lag 800. Restricted integrals retain the signed area through lag 800
without an exponential fit, positive-only truncation or extrapolation.
All 12 complete independent seeds per size/strength enter each noisy group.
Intervals resample whole seeds and recompute the original ratio. Time origins
within a trajectory are never counted as independent replicates.

## Material site change, reference arcs and changing circumference

The dominant current half-contrast region is selected by its excess active
material, and its circular centroid is computed from that same excess.
An independent saved-profile replay checks these centroids, rather than
validating only region counts or moment values.

The wrapped endpoint material shift and two arc metrics are

$$
\Delta q(t,\tau)
=\operatorname{wrap}_{[-1/2,1/2)}
  [q_{\mathrm{site}}(t+\tau)-q_{\mathrm{site}}(t)], \tag{5}
$$

$$
d_{\mathrm{basal}}=C_{\mathrm{basal}}\Delta q,\qquad
d_{\mathrm{mid}}=\frac{C(t)+C(t+\tau)}{2}\Delta q,\qquad
\Delta C=C(t+\tau)-C(t). \tag{6}
$$

The first is the primary basal-reference arc shift; the second is an explicit
instantaneous-metric sensitivity. Circumference change is reported separately.
Subtracting \(C(t)q(t)\) coordinates across different circumferences would mix
uniform stretch with site change and depend on the choice of ring cut.
Neither arc metric is molecular transport, a whole-cell displacement or the
unwrapped distance travelled by an activity structure.

The original strict endpoint gate requires one current half-contrast region,
contrast at least 0.05, at most one peak with prominence at least 0.05, and
relative first-moment magnitude at least 0.1 at both endpoints. A separate
one-region sensitivity removes the peak-count veto. Both tables retain pair
coverage and contributing seeds; a sparse or unavailable gate is not evidence
for stability or instability.
Population curves are displayed only when at least six independent seeds
contribute and at least 1% of endpoint pairs qualify; every computed value,
including sparse or undefined results, remains in the tables.

The native site-angle cosine and common-reference projection are

$$
K_{\mathrm{native}}(\tau)
=\left\langle\cos(2\pi\Delta q)\right\rangle_E,\qquad
K_*(\tau)
=\left\langle\cos\!\left(\frac{2\pi d_{\mathrm{basal}}}{C_*}\right)\right\rangle_E,
\quad C_*=3.0451682524. \tag{7}
$$

Both use exactly the same eligible endpoint set \(E\) and pooled seed
weights. The largest basal circumference provides the reference, so remapped
shortest displacements remain at most a half-turn. Original native shortest
endpoints already fold any travel beyond a half-turn; the remapping cannot
recover winding or path length. It is a geometric coordinate comparison, not
an evolved counterfactual, additive causal partition or reconstruction of
global activity-vector memory. Even in a single-region profile, the region
centroid and global first moment can differ because background and tails
contribute to the latter.

For a uniform wrapped endpoint angle, the corresponding basal-arc RMS is
$C_{\mathrm{basal}}/\sqrt{12}$. This is a reference for circular mixing, not an
upper bound or an inferred diffusion coefficient. The main Spring remapping
panel uses early lag 4; supporting curves retain the full lag range and explicit
4/50/800 comparisons. Late shortest-arc differences near the uniform reference
cannot recover the physical travel or winding discarded by wrapping.

## Passive signatures and limits

Multiple current half-contrast regions and weak-resultant occupancy use all
6401 scalar sample times per seed. Multiple regions are not automatically
equal, opposed, persistent, independent or functional competitors. Small
global asymmetry can coexist with local activity. The sampled angle-increment
distributions retain magnitude gates and whole-seed intervals. They do not
define an abrupt-event class, and small increments at cadence 0.5 cannot
exclude faster unobserved changes.

The two main claims concern differences between physical and angular readouts,
and the added information in geometric remapping plus passive local-pattern
and magnitude measurements. They do not establish a universal size law,
stationary persistence, biological noise calibration or a mechanism for
whole-cell migration.

The [full audit](supporting/metadata/spring_noise_analysis_audit.json),
[centroid/half-window/increment/forcing supplement](supporting/metadata/spring_noise_supplement_audit.json),
[prospective sensitivity gates](supporting/metadata/spring_noise_refinement_plan.json)
and [actual paired audit](supporting/metadata/spring_noise_bridge_audit.json)
identify source hashes and complete coverage. Native worker records and
production state remain in the [manual-refresh status reader](spring_noise_status.html).
