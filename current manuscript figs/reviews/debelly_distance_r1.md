# Independent review: executed deBelly two-arc distance pilot

17 October 2026. Fresh-context scientific critique of the six completed trials.
I inspected both rendered figures, the full-field archives, the native-interval
accounting summaries, the frozen adapter and the preparation/design records.
No simulations, tuning or changes to the authored draft were performed.

## Scientific decision

**Keep this as a useful construction result and discuss the coupling assumptions
before selecting it for a main-text separation mechanism.** The implementation
checks support the declared experiment. The unequal-cut result is more specific
than disagreement between flanks: by 14000 s the short arc is spatially uniform
in Rac and the long arc is nearly uniform at a different activity level, in both
cue and gap-only trials. The construction has therefore changed the nature of
the activity pattern. It does not presently give a clean assay of competition
between two localized sites on the same cortex.

The legitimate finding is that an unequal partition of this qualified
single-incumbent state, with independent local endpoint laws and shared MCA
correction, produces markedly different arc-wide chemical/mechanical states
without an added Rac cue. That is scientifically informative about the proposed
extension. It is not yet an explanation of increased migration persistence in
WGD-associated enlarged cells. Full normal-shape mechanics or a new 2D campaign
is not required to reach this decision.

## What the actual fields establish

The [main result figure](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_two_arc_result/debelly_two_arc_distance.png)
correctly shows both flanks and separate spatial arcs. Its lower row is decisive:
the near and intermediate cases are essentially two plateaus, separated by the
gaps. They are not profiles with a localized old peak and a localized new peak.
I independently calculated the following ranges over **every saved spatial cell
at 14000 s**, using the unchanged archives. The initial full-loop Rac contrast
is 0.27291.

| Short arc d (µm) | Trial | Short-arc Rac range | Long-arc Rac range | Fraction of total MCA in short arc |
|---:|---|---:|---:|---:|
| 18.6617 | Cue + gap | 2.89 × 10⁻¹⁵ | 0.001750 | 94.18% |
| 18.6617 | Gap only | 2.89 × 10⁻¹⁵ | 0.001813 | 93.99% |
| 24.8308 | Cue + gap | 1.23 × 10⁻¹⁴ | 0.000657 | 97.42% |
| 24.8308 | Gap only | 1.33 × 10⁻¹⁴ | 0.000922 | 96.39% |

Short-arc Rac is approximately 0.664847 throughout; long-arc Rac lies around
0.892–0.893. The long-arc spatial ranges are only 0.24–0.66% of the initial
whole-loop contrast. Calling this only “flank asymmetry” misses the larger
change: each nominal physical site borders a low-activity arc on one side and
a high-activity arc on the other, while the distinction between incumbent and
challenger largely disappears within either arc. A site average or circular
first moment would conceal this failure of a common site identity.

This is a Rac-pattern observation, not a claim that all fields are homogeneous
or settled. At 14000 s the high-Rac long arcs retain cortical flow of both signs
and nonzero tension. The [native-field companion](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_two_arc_result/debelly_two_arc_native_fields.png)
also shows continuing MCA redistribution. The full-field audit reports maximum
scaled pointwise temporal ranges of 0.008616–0.01929 over 13000–14000 s for the
unequal cuts. Endpoint Rac has stopped being an adequate settling diagnostic.
No eventual fate or newly qualified steady state should be assigned.

The cue has a clear transient chemical effect at all three locations; the
maintained challenged-flank Rac exceeds 4.57 by withdrawal. In the unequal cuts,
however, cue and gap-only trials approach very similar Rac patterns by the
displayed endpoint: their maximum final spatial Rac differences are 5.81 × 10⁻⁵
and 2.47 × 10⁻⁴. They are not identical complete states: bound MCA differs by up
to 0.000504 and 0.002010, respectively, and relaxation remains unresolved.
Thus “the cue has no effect” would also be too strong.

At d = 31 µm, both arcs retain Rac contrast of approximately 0.27291. The cue
trial has a challenger-centered profile, whereas the control retains an
incumbent-centered profile. Their displayed final windows are quiet across all
eleven fields, but no independent recovery holdout was prescribed. This is
finite-time evidence of cue-dependent site selection in the equal-arc geometry.
It should not be promoted to a qualified asymptotic fate or a three-point
distance law.

The independent [read-only field check](debelly_distance_r1_field_check.json)
records the final min/max/range of all eleven fields on both arcs, full-field
tail ranges, MCA fractions and cue/control differences. All six archive hashes
match the executed independent audit.

## The shared MCA operation is a central modeling choice

There is no diffusive passage across either gap. Nevertheless, the arcs are
coupled through the uniform correction after their synchronized native
transport/exchange stages. Each receives the same concentration increment,
so the returned amount is proportional to arc length, regardless of where MCA
left the modeled cortex. This coupling must be considered alongside the
independent flank controllers, rather than treating the arcs as wholly
independent chemical systems.

The [integrated per-arc budgets](../../polarity1d/manuscript/data/runs/debelly_two_arc_distance_20261017/integrated_arc_pool_budgets.json)
make the magnitude clear. In the **24.83-µm gap-only control**, the short arc's
total MCA rises from 2.31922 to 5.97627, while the long arc falls from 3.88078
to 0.223727. The short arc receives 4.11691 by uniform correction, loses 0.711457
at its endpoints, and has a positive 0.251598 native sequential-exchange
residual. These terms reproduce its net gain of 3.65705. No cue is needed for
this major partition of the conserved pool.

The final short-arc bound MCA is approximately 0.188–0.244; the long-arc value
is approximately 0.00336–0.00675. Relative to the unchanged native MCA feedback
threshold of 0.068, these are substantially different local activation regimes.
High bound MCA reduces the native mechanical increment to Rac activation.
The resulting low-Rac, inactive short arc and low-MCA, active long arc are
consistent with coupled retention/loss and feedback. The recorded budgets
establish the accounting and identify a plausible mechanism to discuss. They
do not independently prove which feedback initiated the partition or how much
of the later pattern each pathway caused.

Uniform correction is an implementation operation; it should not be renamed
biological recycling or a finite-speed shared cytoplasmic reservoir. Likewise,
the native sequential-exchange residual is not an additional biological source.
It is substantial in the integrated budgets and is correctly kept separate.
If this partition later becomes a quantitative mechanistic claim, its sensitivity
to the native time discretization will matter. A new convergence campaign is
not needed merely to report this bounded construction finding.

## Coding checks, startup surgery and the remaining uncertainty

The passed checks materially strengthen trust in the result: exact original
bulk fields and incumbent controllers, native equal-arc equivalence, unequal-arc
orientation invariance, the transparent pre-cut operator, matched starts,
event timing, finite/nonnegative chemistry and explicit pool balances. The
observed arc contrast should not be dismissed as an unverified sign error,
mass leak or omitted controller restart.

Opening a new boundary still changes the model immediately. At the unequal
cuts, the inherited new-flank membrane velocity is −0.006880 or −0.003211 µm/s
in arc 1, whereas its first gap-only step gives +0.019325 or +0.009899 µm/s.
The corresponding jumps are +0.026205 and +0.013110 µm/s; the diametric jump is
zero. Initial new cortical commands are zero in all cases. This is a genuine,
geometry-dependent response to the newly imposed displacement-relaxation
controller, even though position is transferred exactly.

It is reasonable to call this boundary surgery or a startup impulse. It is
not established that this impulse alone causes the later arc-wide activity.
Opening also permanently removes chemical face couplings, changes the elliptic
mechanical boundary problem and introduces independent endpoint dynamics. The
existing cue/control experiment separates added cue from gap installation; it
does not separate these components of the gap intervention. Extending recovery
or repeating the same cue would not resolve that causal ambiguity.

## Decision to discuss before more experiments

The most useful next step is to agree what one unmodelled protrusive site
means when its two flanks cease to agree. Three related physical choices need
to be explicit:

1. **Site identity.** Are the two flanks autonomous boundary patches, or do they
   belong to one protrusive region with shared chemical dynamics, mechanics,
   or force balance? A shared endpoint command alone need not restore a
   localized Rac site. Neither averaging before the nonlinear law nor enforcing
   equality afterwards is a neutral repair.
2. **Chemical connection.** Which fields can exchange through the omitted
   protrusion or cell interior, and what represents that exchange? No
   cross-gap diffusion together with globally uniform MCA replacement is a
   particular hybrid coupling, whose consequences are now visible.
3. **Meaning of creating a challenger.** Is opening a permanent new boundary,
   including its controller relaxation, the intended biological intervention?
   If yes, the claim should concern challenger creation under that declared
   boundary geometry. If the intended question is separation between sites
   on a connected cortex, the interface assumptions require discussion first.

No choice needs to be tuned until it yields a favorable monotone result. The
present output has already supplied the information needed for that modeling
discussion. It also explains why an additional flank coupling cannot yet be
advertised as the uniquely required fix.

## Figure use and connection to the manuscript question

The current figures are useful evidence sheets. For a compact discussion with
the user, show the same starting incumbent and two routes, the near-cut cue and
gap-only **full Rac profiles**, the equal-arc contrast, and a small MCA-per-arc
summary. The arc-wide plateaus and shared MCA operation deserve the visual
emphasis; a page dominated by cue-on flank Rac can make this look like an
ordinary two-site competition assay. Keep both unequal separations in support.

Do not make this result carry either main-text candidate on its own. The
enlargement–incumbent-architecture question and the independent challenger
separation question remain useful distinct hypotheses. This fixed-reference
circumference experiment examines a proposed geometry for the latter. It has
not measured a minimum successful cue, spontaneous polarity switching,
cellular migration or directional persistence, and it does not vary cell size.
The observed migration summaries, conditional OU persistence estimates and
these deterministic polarity patterns remain different readouts.

## Anticipated user critique

**Anticipated user critique — hypothetical, not an actual user quotation:**
“I asked whether two protrusions farther apart compete differently. In the
unequal cases you have given me two almost uniform arcs and moved most of the
MCA into one of them, even without the cue. Are those really two protrusions
on one cell, and which physical assumption made that happen?”

That critique is supported by the saved fields. The appropriate response is
to show the arc profiles and pool accounting, explain the declared coupling,
and agree the intended site/interface physics before further experiments.
Neither treating this as a biological distance effect nor discarding it as a
generic failed simulation would answer the useful question it has exposed.

## Review boundary and sources

Read the repository and manuscript `AGENTS.md` instructions,
[current report status](../../reports/CURRENT_STATUS.md),
[four-assay deprecation warning](../../polarity1d/four_assay/README.md),
[biological figure context](../README.md),
[migration figure caption](../../wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU_caption.md),
[circular preparation result](../circular_polarization_result_20261017.md),
[prospective two-arc review](two_arc_design_r1.md),
[executed result and protocol](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_two_arc_distance_result.md),
[startup audit](../../polarity1d/manuscript/data/runs/debelly_two_arc_distance_20261017/startup_audit.json),
and [full-archive audit](../../polarity1d/manuscript/data/runs/debelly_two_arc_distance_20261017/independent_audit.json).
The frozen `worker.py`, native feedback laws and unchanged six full-field
archives were inspected directly. This is a scientific interpretation review,
not an additional integration, asymptotic qualification or mesh-convergence
certification.
