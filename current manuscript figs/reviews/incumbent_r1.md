# Fresh-context critique: incumbent candidate, round 1

Reviewed the [candidate Markdown](../incumbent_candidate_review_20261017.md), its actual rendered four-panel PNG, technical caption, measurement CSV, methods and verification metadata, then the central discussion and current report's relevant text. Read root and manuscript AGENTS.md, the current-status guidance and deprecated four-assay README before interpreting. This is an editorial/scientific review of retained evidence; no simulations were run.

## Anticipated user critique

> This is more specific, but I still do not think we have established what the main figure is supposed to teach the reader. Goryachev gives a substantial effect and Wave-pinning gives a useful contrast. That does not yet make “increasing incumbent scale” an explanation: the enlarged incumbent, depleted inactive pool and changed reaction response occur together. If the message is that enlargement can alter incumbent maintenance and challenger responsiveness, show those quantities precisely. If the message is greater stability, this figure has not measured it. Also, the perturbation behind the most striking panel is invisible. I should not have to reconstruct it from a technical caption before I can tell what 15.7-fold means.

## Scientific judgment

**Retain this as a serious chemical-mechanism candidate; do not yet select it as the one main-text figure.** Its strongest result is specific: in the displayed Goryachev regime, longer cells have less inactive substrate at the challenged end and a larger immediate incumbent restoring reaction per unit induced damage. Wave-pinning shows that a broader incumbent can accompany nearly unchanged initial local activation and a much smaller change in that restoring diagnostic. That contrast is informative: increased active amount alone does not guarantee these two response trends.

The defensible message is:

> Enlargement accompanies different incumbent organizations. In Goryachev, concentration into a narrow incumbent coexists with depleted challenger substrate and enhanced immediate restoring reaction; broadening in Wave-pinning has different local response trends.

“Accompanies” and “coexists” matter. Neither panel independently changes architecture while holding the rest of the polarized state fixed. The comparison supports a mechanistic interpretation through known kinetics, but does not identify incumbent width or height as the causal control variable.

| Claim | Judgment from the displayed evidence |
| --- | --- |
| The incumbent becomes taller in Goryachev and broader in Wave-pinning | Supported, although Goryachev's width is poorly resolved in this rendering. |
| The same local added cue initially drives less activation in longer Goryachev cells | Supported: the added source is local cue strength times local inactive concentration; the L39/L9 ratio is 0.333. This is the source contribution, not total net growth. |
| Longer Goryachev cells give a larger immediate restoring reaction per damage in this perturbation | Supported: the L39/L9 ratio is 15.7; Wave-pinning is 1.28. |
| Enlargement causes greater incumbent stability or a stronger restoring force | Unsupported here. The diagnostic is reaction-only, direction-specific and immediate. It excludes diffusion, supplies no global stability measure and is not a mechanical force. |
| Competitor formation or takeover is suppressed | Unresolved. Initial source strength does not establish formation, survival after withdrawal or takeover. |
| These changes explain WGD migration persistence | Unresolved. This is a candidate polarity mechanism, with no calibrated WGD fold change or migration readout. |

Panel C is an exact re-expression of the inactive concentration already shown in A–B. That is legitimate when used to explain the functional meaning of substrate depletion, but it should not be counted as independent evidence. Panel D is the genuinely additional probe.

The whole-inactive-profile intervention is identical in fractional removal and hold duration across lengths; the resulting active-core damage is not matched. In Goryachev its fractional size is approximately 0.0010% at L9 and 0.00049% at L39. Dividing by damage yields a useful immediate response coefficient, but does not establish resistance to finite competitor-induced damage. These small denominators do not by themselves invalidate the measured trend: the raw excess reaction also increases strongly.

## Most important concrete revision

**Make D intelligible as an experiment, and align the claim with that experiment.** Show a small protocol inset: remove 0.5% of the whole inactive profile for 0.1 native time units, then return that original profile pointwise; measure excess basal reaction within the original active core at return. Label the outcome “Immediate restoring reaction per damage,” including “reaction excludes diffusion” in the figure or adjacent caption. Use “relative to each model's own L9” explicitly. Do not label it stability or restoration speed.

If a small retained core-deficit trace can fit, it would make damage and return visible, but distinguish its net dynamics from the reaction-only numerator. More numerical precision or additional length samples would not resolve this communication problem.

Secondary rendering revisions:

- Give Goryachev a clearly labelled incumbent-end enlargement/inset. Its roughly 0.7–0.8-unit core is compressed into a tiny part of a 0–40 axis. The purported organization contrast should be visible without inference from the caption.
- Resolve the color-role switch: purple denotes L21 above and Goryachev below. Direct model labels on C–D and a clear separate model key would reduce the reassignment burden.
- Mark incumbent and challenged ends and state the fixed mean pool concentration/local kinetics assumption. Total available material grows with length by construction; “more incumbent material” needs that context.
- Preserve Wave-pinning's near-null C result. It is an informative comparison, not a failed demonstration. Avoid enlarging the small deviation from unity into an apparent mechanism.

## Priority and additional work

This figure is suitable if the chosen main-text question is **how enlargement changes local chemical responsiveness in contrasting polarity systems**. It is less suitable if the question is **why a prompt challenged-end response can coexist with a persistent incumbent**: it contains neither a challenged-end trajectory nor incumbent retention during a challenge. The retained Goryachev/deBelly time-course contrast is more direct for that question. Combining every diagnostic or all five models would obscure the decision about which insight matters most.

**No additional simulation is necessary to make this candidate scientifically acceptable under the narrow claim.** Clarify the intervention, improve profile readability and use the retained sources. Causal attribution to incumbent scale, global stability or successful challenger suppression would require a separately defined test; those are stronger questions, not prerequisites for presenting these local observations honestly. The circular distance experiment addresses a separate hypothesis and should not be used to patch this candidate's interpretation.
