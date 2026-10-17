# Spring ring preparation: reviewed development closure

Root reviewed and approved equation (3) for a bounded first feasibility pilot on
17 October 2026. The fresh frozen preparation run is
[spring_ring_preparation_pilot_20261017](data/runs/spring_ring_preparation_pilot_20261017/),
submitted as SLURM **21083215**. This is a development model. The first milestone
is one self-maintained polarized ring state at native-like local parameters.
No challenger sweep or final figure claim is selected here.

**First pilot outcome:** SLURM 21083215 completed successfully in 41 seconds.
The first 10-unit cue produced one active region and one resolved peak. Both
chemical fields and circumference passed two successive 1000-unit settling
windows, followed by a 10000-unit holdout; total cue-free observation was 13000
native time units. The final circumference is 2.04418, and active Rac spans
0.70806–1.20630. Maximum sampled physical-amount error is
\(3.80\times10^{-12}\); maximum scaled holdout excursion is
\(3.70\times10^{-14}\). No fallback cue or challenger was run.

![Spring circular global-extension preparation: cue-selected active profiles, evolving extension, cue-free field history and physical-amount conservation.](figures/length_mechanisms_20261008/dev/spring_ring_preparation_pilot_20261017/spring_ring_preparation_feasibility.png)

[Single-page vector PDF](figures/length_mechanisms_20261008/dev/spring_ring_preparation_pilot_20261017/supporting/vectors/spring_ring_preparation_feasibility.pdf)
· [Caption](figures/length_mechanisms_20261008/dev/spring_ring_preparation_pilot_20261017/supporting/metadata/caption.txt)
· [Saved-state audit](data/runs/spring_ring_preparation_pilot_20261017/results/audit.json)

The figure establishes this preparation feasibility point. It displays the
active region together with the global extension that feeds back on signalling,
so successful polarity is not attributed to a silently fixed mechanical state.
The kinetic, force and dilution coefficients remain native; the circular
force-averaging closure remains an explicit adaptation. The final four-panel
layout is a development check, not a proposed final manuscript layout.

The SciPy BDF implementation emitted one internal-difference `RuntimeWarning`;
its exact text is retained in the audit. Every saved field and the final basal
right-hand side are finite, the final derivative is below the settling
tolerance, and all saved-state checks pass. The internal cause of the warning
is not established here. No independent replay or spatial-convergence claim is
made from this single discretization.

The native source is the authors’ immutable
[Figure 4C script](../../reports/trial_report_assets/zmurchok2020/authors/Fig4/Fig4c.py),
recorded at commit `f78b05a16a99d13d3ef86fc919d49ebdee479cee` in
[source.json](../../reports/trial_report_assets/zmurchok2020/authors/source.json).
The active interval adapter is [spring.py](code/core/models/spring.py).

## What the original mechanics specify

Adding the original two endpoint equations gives

\[
\dot L=F(u_-)+F(u_+)-2k(L-L_0), \qquad k=0.01.
\tag{1}
\]

The local force switch, global deactivation and dilution are

\[
F(u)=\frac{10^{-3}}{1+\exp[-20(u-1)]},\qquad
\delta=3+80(L-L_0),\qquad
\text{dilution}=-\frac{\dot L}{L}(u,v).
\tag{2}
\]

Viscosity is one. The two endpoint forces have a defined mechanical role;
chemical periodicity removes those endpoints and therefore does not, by itself,
specify the moving-domain Spring model.

The [original paper’s 2D formulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC7642449/)
applies Rac-dependent traction normal to the boundary, replaces length excess
with area excess for global deactivation, and retains deformation-dependent
chemical transport and dilution. Its mechanics use a deforming 2D continuum,
not a one-variable circumference. The authors’ local README states that the 2D
code is available upon request. This supplies a source-based motivation for
distributed boundary forcing, but does not uniquely prescribe the ring below.

## Recommended minimal closure

Let the ring have circumference \(C(t)\), resting circumference \(C_0\), and
uniform material coordinate \(q\in[0,1)\). Use periodic chemical transport and
the same pointwise Rac-to-force switch. Replace the mean of the two native
endpoint forces by the mean around the ring:

\[
\dot C=2\left[\langle F(u)\rangle_q-k(C-C_0)\right],
\qquad \langle F(u)\rangle_q=\int_0^1F(u(q,t))\,dq.
\tag{3}
\]

Equation (3) is the **new mechanical closure**. It retains the uniform-field
native length equation exactly and assigns every ring position the same
mechanical status. It averages the local force after the nonlinear switch;
\(F(\langle u\rangle)\) is a different rule and is not used. The factor two
preserves the native scalar extension relaxation time rather than inserting a
new circumference-dependent gain. This is an explicit global-extension
analogue, not the authors’ native continuum mechanics or a derived reduction of
their 2D elasticity.

The chemical equations are

\[
\begin{aligned}
\partial_tu\big|_q&=\frac{0.01}{C^2}\partial_q^2u
 +\left[4+s(q,t)+\frac{5u^6}{1+u^6}\right]v
 -[3+80(C-C_0)]u-\frac{\dot C}{C}u,\\
\partial_tv\big|_q&=\frac{10}{C^2}\partial_q^2v
 -\left[4+s(q,t)+\frac{5u^6}{1+u^6}\right]v
 +[3+80(C-C_0)]u-\frac{\dot C}{C}v.
\end{aligned}
\tag{4}
\]

These preserve the source reaction coefficients, diffusivities, baseline
activation, force amplitude and threshold, extension-dependent deactivation,
and dilution. With periodic fluxes, the physical amount obeys

\[
\frac{d}{dt}\left[C(t)\int_0^1(u+v)\,dq\right]=0.
\tag{5}
\]

No concentrations are clipped or renormalized. Starting at \(C=C_0\), the force
bounds imply \(0\le C-C_0\le0.1\); feedback therefore remains present rather
than being frozen at a chosen extension. This is a bound, not an observed
trajectory.

## Why a single protrusive gap is not the first Spring construction

Unlike the deBelly implementation, native Spring does not omit the protrusive
region from its chemical domain. Opening a gap would require a new chemical
transport or reservoir rule. Retaining only two actuator flanks at an incumbent
gap would privilege that position: a challenger elsewhere would have no direct
native actuator. Opening a second gap would additionally require rules for
combining endpoint forces, elastic coupling, material removal and dilution
across segments. Such an arc model can be constructed, but it adds more choices
than equation (3) and is unnecessary for the initial periodic polarization test.

## Bounded first preparation pilot

Use \(C_0=2\) in native Spring length units, within the existing interval range,
and 200 periodic cells. This is a feasibility point, not a matched conversion
from interval length to cell perimeter. Start from the source concentrations
\(u=0.05645\), \(v=2-u\), \(C=C_0\), without importing an interval checkpoint
or adding the nonperiodic cosine perturbation. Initial physical amount is
\(2C_0\).

Apply one symmetric cue with native extra activation peak four and total support
0.2, then remove it exactly. A raised cosine changes the cue geometry to respect
the periodic domain; it is not dose-matched to an endpoint half-patch. An initial
10-unit pulse is a bounded setup choice, followed by native basal evolution.
Retain full active/inactive fields, circumference, extension, global deactivation,
mean force, and physical amount throughout. Use the existing ring pilot’s
single-site, complete-field settling and holdout criteria, additionally requiring
settling of circumference. A failed preparation remains a failed first attempt;
it is not evidence that Spring cannot polarize on a ring.

All trajectory integration runs via SLURM in the existing `scaling_sandbox`
environment, one CPU, 3 GB, QoS `small`, a one-hour allocation and fresh frozen
sources. The BDF solver uses a full dense numerical Jacobian: the mean-force
and dilution terms couple all active cells globally, so a purely local periodic
sparsity mask would omit real dependencies. Root approved this bounded pilot
after reviewing equation (3) and the algebraic implementation.

For a later distance experiment, all challenger positions would share one exact
prepared state and starting circumference. Circumference could subsequently
evolve differently with challenger location through the retained feedback. Such
changes would be measured as responses of the model, rather than silently
constraining circumference and removing extension feedback.

## Implemented checks and current limit

The [algebraic prototype](code/mechanisms/spring_ring_closure_20261017.py) imports
the native local force/reaction functions and has no integration routine. Its
[check record](data/runs/spring_ring_closure_review_20261017/operator_checks.json)
tests physical-amount conservation at the derivative level, cyclic-rotation and
reflection covariance, exact homogeneous-field agreement with the complete
native derivative, and averaging after the force switch. These checks establish
implementation identities only. The fresh trajectory producer has its own
saved-state audit and figure exports. Spatial convergence and distance effects
remain outside this first pilot.
