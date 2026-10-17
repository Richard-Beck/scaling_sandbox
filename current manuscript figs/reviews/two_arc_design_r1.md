# Independent review: asymmetric deBelly two-arc distance design

17 October 2026. Fresh-context scientific review; no challenger integrations were
run and no implementation sources were changed. I did not design this adapter.

## Decision

**The construction is acceptable for a bounded development pilot, subject to
the numerical checks below. It is not yet a main-text explanation of greater
persistence in enlarged cells.** Tangential mechanics on a fixed reference loop
is a legitimate declared reduction; a full normal-shape or axisymmetric solver
is not a prerequisite for this question.

Its strongest prospective contribution is testing how a localized challenge
interacts with an incumbent-organized chemical and mechanical landscape while
total reference circumference and total MCA remain fixed. A result consisting
only of longer response times at larger separation would leave the physical
explanation unresolved.

## What the intervention physically means

The proposed start is one qualified, cue-free incumbent on the C = 62 µm,
402-cell, one-cut reference loop. A new point cut at a selected original grid
face creates two segments of lengths d and C − d. No material is removed.
Each segment has an incumbent flank and a challenger flank, with its local
coordinate pointing from incumbent to challenger. The existing local Rac and
tension laws independently drive all four endpoints. Rac/Rho and MCA have no
diffusive passage through either unmodelled protrusive region. Mechanical fields
are solved separately on the segments, and one uniform MCA correction spans
both segments after their synchronized transport/exchange stages.

This is a specific two-segment mechanochemical model. It does not contain a
resolved cap, force balance between a cap's two flanks, shared cap Rac dynamics,
or finite-speed transfer through a cap. Independence of the flanks is an
acceptable explicit extension of the corrected endpoint model. Calling their
two signals a single biological site's activity requires displaying their
agreement, or acknowledging their disagreement. Averaging them before applying
the nonlinear command law would change the model and should remain excluded.

“Fixed circumference” must mean fixed **reference** circumference. The native
membrane displacement and endpoint relaxation remain dynamical. The drawing
does not establish a freely moving, physically circular contour, and these
experiments do not measure migration velocity or observed migration persistence.

## Decisive interpretive limits

1. **Separation and boundary geometry change together.** Moving the cut changes
   both arc lengths, the mechanical boundary-value problem, transport boundary
   conditions and the local incumbent-state values sampled by the new
   controllers. The latter are part of the incumbent's spatial landscape and
   may be the interesting mechanism; they must be measured rather than hidden
   inside a distance label. Both routes, d and C − d, matter. Neither the
   shortest route alone nor a straight-line chord is automatically the relevant
   length.

2. **A cut-only control does not subtract the geometry away.** At each d it
   estimates what introducing the gap and its native controllers does without
   the added Rac cue. The cue-versus-cut contrast is the cue's effect within
   that geometry; its variation across d includes the interaction of cue and
   geometry. That is a valid controlled comparison, but not an isolated
   molecular travel-distance effect. If a cut-only trajectory recruits a new
   front, weakens the incumbent, or selects a different branch, that is a
   substantive result and a limit on the cue-driven interpretation.

3. **The model already has instantaneous long-range operations.** Quasistatic
   elliptic mechanics responds throughout each segment on a mechanical solve.
   Uniform MCA correction changes both segments during the same staggered step.
   These operations preclude an unsupported narrative that an inhibitory
   mechanical signal or replenished MCA simply takes longer to travel farther.
   A transport explanation must identify the relevant evolving pool, local
   depletion/accumulation or profile change, and its measured dynamics. MCA
   centroid displacement is not molecular travel.

4. **The newly imposed relaxation can create a distance-dependent transient.**
   Initializing new controller displacements to the precut face displacement
   preserves the state at the instant of opening. It does not preserve the
   previous material velocity: native relaxation immediately acts on that
   displacement, even when the new Rac command is negligible. At a nonsymmetric
   cut this can differ from the diametric case. Record the new endpoint
   membrane velocities and commands relative to the inherited precut motion,
   and show the corresponding early cut-only response. This is an intended
   boundary intervention, not a restart identity failure.

5. **This single-dose pilot is not a threshold or size experiment.** Three
   separations, one 1000-s cue and 13000-s recovery can reveal a useful contrast,
   coexistence, switchback, or no resolved contrast at that dose. They cannot
   locate a minimum switching duration, establish a scaling law, or by
   themselves show what enlarging the cell does. A later main-text argument
   must connect a measured separation mechanism to the existing length
   experiments, or test its predicted consequence when circumference changes.

## Numerical checks that matter

These checks are requirements for trusting the new adapter, not claims that
they have already passed. The reviewed draft included exact seven-field
reconstruction, physical global-mass checks and diametric native one-step
comparisons; additional checks were still being added.

| Check | What a convincing pass establishes |
|---|---|
| Exact start and old controllers | Reassembling both arcs reproduces every original bulk sample exactly; both original raw and relaxed incumbent-controller states are preserved with their orientation signs. The new controller initialization is audited separately against the declared face interpolation. |
| Transparent interface | Splitting and reconnecting an interior interface with continuous mechanics and the original chemical diffusion/advection face flux, without new endpoint controllers, reproduces the original one-cut operator. A closed gap is not a transparent interface. |
| Equal-arc native equivalence | At d = C/2, independently commanded native endpoints with `center_boundary_motion=False` reproduce all eleven state fields for both cue and cut-only steps. This checks the declared extension, not its biological validity. |
| Unequal-arc orientation | Explicitly reverse/reindex a nonsymmetric test state and cue: displacement and velocity change sign; chemistry and scalar tension do not; endpoint labels and flux signs swap consistently. Test all stages. Equal-arc equality alone can conceal orientation errors. |
| Shared MCA accounting | Sum physical quadrature over both arcs once, then apply the same concentration increment everywhere. Retain per-arc pre/post-correction masses and all four signed bound/unbound endpoint fluxes. A conserved total alone can conceal erroneous transfer between arcs. Separate the sequential-exchange residual from physical boundary loss. |
| Cue identity and event handling | Every d receives the same sampled peak, profile, support and integrated dose on the unchanged grid; withdrawal occurs exactly at 1000 s. No unintended overlap with the incumbent cue/core or duplicate endpoint dose occurs. |
| Paired control identity | Cue and cut-only trials at each d share exact bulk/controller initial state, geometry and all noncue parameters. Preserve finite, nonnegative fields and native timing throughout. |

Native equivalence is most informative at the individual staggered stages and
over a short multistep trace, including withdrawal; one endpoint comparison at
one time is weaker. Spatial/time refinement can wait until a consequential pilot
contrast exists. It then becomes necessary before presenting a small quantitative
separation effect as a result.

## What a persuasive figure would actually show

Keep the main result to one readable figure if possible. A geometry panel
should label both routes, the one shared MCA accounting domain and four local
flanks. Beside it, show the unchanged starting bound-MCA and tension profiles
with the three challenger positions. This makes the candidate physical landscape
visible before presenting an outcome.

The result panels should compare cue and matched cut-only traces at every d
using common physical windows. At minimum display challenger and incumbent Rac
for both flanks, challenger local bound MCA, and the mechanical variable that
supports the proposed mechanism. Local regions must have fixed physical width
across d. If the proposed explanation involves transport, show its measured
local bound-MCA budget, retaining advection, diffusion, net binding and uniform
correction until their relative sizes justify simplification. A full budget may
belong in the linked methods, but a causal sentence needs a displayed observable.

Define outcomes before inspecting the traces. Challenger activation alone is
not incumbent replacement. Keep incumbent survival, challenger-only selection,
coexistence, unpolarized outcomes, switchback and unresolved recovery distinct.
A circular first moment can rotate because a second site grows while the first
persists; it cannot be the only switching criterion. At nonantipodal positions,
“site replacement” is more precise than “reversal.” Claim a stable final field
only after the required field-level stationarity assessment; a fixed recovery
horizon alone does not provide it.

An acceptable prospective claim would be: **at fixed reference size and total
MCA, the incumbent's spatial MCA/mechanical state changes how a second site
responds to the same localized cue.** The actual sign, relevant field, time
window and post-cue outcome must come from the pilot. A causal allocation to one
pathway would need an evolved intervention; a source budget on coupled native
trajectories only describes their contributions.

## Anticipated user critique

**Anticipated user critique — hypothetical, not an actual user quotation:**
“What exactly has to move farther, and why should that make a larger cell harder
to repolarize, when your mechanics acts everywhere immediately and your MCA
replacement is uniform? Have you measured persistence, or only the consequence
of cutting the cortex at a different place?”

The present design can address the second part by quantifying cut-only effects
and preserving flank-level incumbent/challenger outcomes. It can answer the first
only if the data identify a concrete local landscape or transport bottleneck
that bridges to the existing length result. Declaring the geometry clearly and
passing numerical equivalence checks are necessary; neither supplies that
mechanistic bridge by itself.

## Sources and review boundary

Read the repository and manuscript `AGENTS.md` instructions,
[current report status](../../reports/CURRENT_STATUS.md), the
[four-assay workflow warning](../../polarity1d/four_assay/README.md),
[current figure context](../README.md),
[circular preparation result](../circular_polarization_result_20261017.md),
[corrected endpoint interpretation](../../polarity1d/manuscript/DEBELLY_BOUNDARIES.md),
and the [single-incumbent protocol](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_single_incumbent_prototype.md).

Source inspected: [two-arc implementation draft](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_two_arc_distance.py),
SHA-256 `26b080935dcd4254a8d2fa4c409f4ee93deda654ad26fb55ca5c6cfa2093b1b3`,
and the current native deBelly timestep. The implementation was still changing;
this review does not certify later source versions, an executed preflight, or
any challenger result. No challenge evidence was available to this review.
