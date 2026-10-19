# Broad-strength pilot: exploratory result review

[Two-page review PDF](broad_strength_review.pdf) ·
[Organization overview](broad_strength_exploratory.png) ·
[Native profile comparison](broad_strength_native_profiles.png) ·
[Exact generation receipt](supporting/metadata/RENDER_RECEIPT.json) ·
[Saved producer observation](supporting/metadata/generation_snapshot.json)

This separate discussion artifact displays existing collected native results.
It leaves both preceding Figure 4 sheets and the accepted weak/strong medoid
review unchanged. Native production and numerical qualification are still in
progress: the overview states the exact producer observation time and completed
output count. This static export does not claim to update itself. The
[manual-refresh collector](../../../../../../../data/runs/revision_noise_strength_broad_20261019/results/partial/index.html)
contains newer observations while its read-only collector runs.

The common strength axis is f = 0.0001, 0.0003, 0.001, 0.003, 0.01,
0.03, 0.1, 0.3, 1. Each model retains its exact original basals, three sizes,
positive lognormal forcing family, native correlation scales and time units.
Increasing f scales the positive input mean and fluctuation amplitude together.
The two original strengths retain twelve independent seeds; each new strength
initially targets two. Mean-input controls are deterministic single cases at
the declared ensemble mean, without normalization to a realized noisy path.
These common f values do not calibrate effective noise strength across models.

## Reading the overview

The left column shows the fraction of all saved scalar times with more than
one connected current half-contrast region. Solid curves retain every sampled
profile; dashed curves additionally require contrast at least 10% of that
case's initial basal contrast. A threshold region is a connected arc above
the current minimum plus half the current max–min contrast. The difference
between these curves exposes counts in nearly flat fields. It does not define
exclusive regimes, establish multiple comparable peaks or demonstrate functional
competition. Native and basal-relative prominence sensitivities, low-contrast
fractions, q2, sampled positive-run spans and half-window checks remain in the
linked source tables.

The middle column shows the mean magnitude of the normalized physical-activity
vector R, with the matched deterministic mean-input result dashed. Constant-mean
loss is an outcome in its own right; for example the previously audited Goryachev
C20 mean control loses its first moment at f = 0.01. Mean input can also remove
the first moment at broader strengths in other models. The noisy and mean-input
histories need not undergo the same response.

The right column separates physical active amount A from normalized shape
coherence Q. All quantities retain the same physical quadrature used in Figure 3,
including the actual Spring circumference at each time:

$$
M_0(t)=\sum_j u_j(t)\,\Delta s_j(t),\qquad
M_1(t)=\sum_j u_j(t)e^{i\theta_j}\,\Delta s_j(t). \tag{1}
$$

$$
P(t)=\frac{M_1(t)}{M_1(0)},\quad
q_1(t)=\frac{|M_1(t)|}{M_0(t)},\quad
A(t)=\frac{M_0(t)}{M_0(0)},\quad
Q(t)=\frac{q_1(t)}{q_1(0)},\quad R(t)=A(t)Q(t). \tag{2}
$$

R can exceed one. Its time/seed mean is the mean of the paired framewise product,
not the product of the separate means of A and Q. An increased active amount can
coexist with reduced coherence and a substantial first moment. No movement law
or angle-to-cell-velocity mapping is assumed.

Noisy symbols encode available replication: filled circles are complete two-seed
groups, squares are retained twelve-seed groups, and open data symbols are
one-path partial groups. Open markers below the left panel mark groups with no
available noisy path; they are availability marks, not zero occupancy.
Whiskers are observed seed ranges, not confidence intervals. Lines connect only
available measured values and guide the eye across a finite grid. No discrete
phase boundaries, final rare-pattern frequencies, monotone strength ordering or
cross-model regime equivalence are inferred.

## Native profile example

The companion uses the first prespecified pilot seed at f = 1 and C = 20 for
Wave-pinning and Goryachev, not a medoid and not a selected maximum. For the WP
path, among saved profiles with more than one current half-contrast region,
contrast at least 10% of basal and R at least 0.75, it selects the frame nearest
the full-window midpoint. This gives the exact native time 16040. Goryachev uses
that same native time on its own independent first path. The time equality
does not calibrate model clocks. This explicitly conditioned diagnostic is not
a typical marginal state or an estimate of its prevalence.

At that frame WP has two threshold regions with R about 0.885, A about 1.922
and Q about 0.461. Goryachev also has two threshold regions, but with R about
0.072, A about 1.091 and Q about 0.066. Their native profiles make the difference
visible. Exact values, native cell samples, seed identities and source hashes
are exported below. Both panels preserve actual concentration amplitudes and
include the exact basal and the matched mean-input profile at the selected
time. Goryachev includes a labelled native-scale enlargement alongside the
full zero-based axis; no independent peak normalization is used.

## Sources and limits

- [Exact collected condition snapshot: all windows, cases and seed ranges](supporting/data/condition_window_snapshot.csv)
- [Profile selection, scalar context and native source hashes](supporting/data/profile_selection.csv)
- [Every plotted profile cell](supporting/data/native_profile_cells.csv)
- [Frozen campaign plan](../PLAN.md)
- [Read-only renderer](../../../../../../../code/mechanisms/revision_figure4_broad_strength_review_20261019.py)
- [Canonical case/window descriptors and sampled-run spans](../../../../../../../data/runs/revision_noise_strength_broad_20261019/results/partial/case_window_descriptors.csv)
- [Actual profile-cadence prominence sensitivity](../../../../../../../data/runs/revision_noise_strength_broad_20261019/results/partial/case_profile_readouts.csv)
- [Prospective numerical/observation checks](../../../../../../../data/runs/revision_noise_strength_broad_20261019/results/partial/upper_strength_numerical_checks.csv)

Overview means cover each complete case's full native window: 0–32000 at scalar
cadence 5 for chemical models, or 0–3200 at cadence 0.5 for Spring. Every case
contributes equally within a condition; these means include the initial forced
adaptation. Full/first-half/last-half summaries are available, with the midpoint
included in both halves. Native profiles have cadence 40 or 4 respectively;
the prominence sensitivity reconstructed from those fields has that actual
cadence, distinct from scalar readouts. A sampled consecutive-positive span
is the last minus first saved positive time, with window-edge censoring retained;
it is not a resolved event lifetime or event rate. Frozen numerical checks
must be reviewed before broader replication. No new native run is launched by
this renderer.
