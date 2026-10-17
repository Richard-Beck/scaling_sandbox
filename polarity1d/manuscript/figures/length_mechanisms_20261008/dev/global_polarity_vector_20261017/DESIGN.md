# Global activity polarity on the supported chemical rings

**Recommendation.** Describe the primary active field by its first angular moment, normalized to its own initial incumbent. Show both its vector magnitude and orientation. This connects the local competition experiments to a measurable whole-cell activity descriptor without assigning velocity or force. Eight existing histories qualify the definition below; no integration or duration campaign was performed.

## Definition and background choice

Let \(u(s,t)\geq0\) denote active species in Wave-pinning/Goryachev, or active Rac in Holmes. Let \(C\) be circumference and \(s_0\) the initial incumbent coordinate. On the idealized circular boundary,

\[
\theta(s)=\frac{2\pi(s-s_0)}{C},\qquad
M_k(t)=\int_0^C u(s,t)\exp\!\left[ik\theta(s)\right]\,\mathrm{d}s.
\tag{1}
\]

Use the dimensionless complex vector

\[
P(t)=\frac{M_1(t)}{M_1(0)},\qquad
\mathbf{P}(t)=\bigl(\operatorname{Re}P(t),\operatorname{Im}P(t)\bigr),\qquad
R(t)=|P(t)|,\quad\phi(t)=\arg P(t).
\tag{2}
\]

Thus the initial incumbent points along \((1,0)\). The denominator both rotates the reference and fixes its amplitude. **Values above one are allowed:** \(R\) is the net oriented activity relative to the initial incumbent, not a fraction of active material. The initial state must have a numerically resolved nonzero moment; an unavailable basal polarity remains unavailable rather than being regularized into a direction.

A spatially uniform background contributes zero to \(M_1\), even if its level changes over time. No Rac threshold, per-frame minimum subtraction or positive clipping is needed. This invariance is an instantaneous descriptor property, not a claim that uniform activation has no effect on subsequent dynamics. Use the cell-width-weighted sum on the retained mesh. Define the supplementary mass-normalized coherence

\[
q(t)=\frac{|M_1(t)|}{M_0(t)},\qquad 0\leq q(t)\leq1,
\tag{3}
\]

and retain \(M_0\), \(|M_1(0)|\) and \(q(0)\) with all curves. Unlike \(P\), \(q\) falls if uniform background activity rises while the oriented component is unchanged. The baseline normalization also removes absolute differences between initial incumbents; those remain visible in the organization panels and these supplementary values. Subtracting a constant background would leave the numerator unchanged but make a total-excess denominator depend on that choice. It offers no advantage here.

## What to show

Three directly interpretable readouts are: **sampled maximum angular excursion** from the initial axis; **minimum relative magnitude and its recovery**; and **retained orientation plus magnitude at a declared post-cue time**. Use the existing onset clock, withdrawal at 1000 and endpoint at 14000, in each model's native units. These times do not calibrate model kinetics. Recovery to \((1,0)\) means recovery of this global descriptor, not identity of every field or survival of a particular site.

Show Cartesian vector trajectories or components alongside \(R(t)\). For an angle trace, withhold \(\phi\) wherever \(R<0.1\); keep those samples visible in Cartesian/magnitude plots and mark them as weak resultant. The 0.1 rule is a display convention relative to initial polarity, not a biological detection threshold. Do not connect angles across masked intervals or infer rotation direction through a zero. Sensitivities at 0.05 and 0.2 are retained. In particular, an apparent antipodal flip can depend strongly on the gate. Report only sampled extrema: recovery is sampled every 5 units for original pilot histories and every 50 for newer histories, so exact crossing times and continuous maxima are not established.

A weak resultant is not necessarily weak local activity. Opposed sites can cancel. A broad merged site can have the same first moment as a very different two-site pattern. Retain the matching profile/kymograph or local-site annotation; \(q_2=|M_2|/M_0\) is a useful cancellation diagnostic, not a unique site counter. No universal turning, migration-speed or persistence law follows from these vectors.

## Declared retained perturbations and illustration logic

The existing complete 652-condition matrix varies circumference and cue placement while holding each model's cue peak, width and duration fixed. The full primary-field histories are retained, so vector heatmaps can be computed without new simulations. Each fixed-size location comparison has the exact same complete basal state. Across size, material, architecture, both routes and angular cue width change together.

The following eight-history audit subset was agreed with the ring-campaign owner before computing these moments. It is not a new prespecified experiment: all underlying outcomes were already known. Do not display eight unrelated examples as a perturbation campaign.

| Comparison or role | Conditions | What is held fixed / what changes |
|---|---|---|
| Common three-model anchor | WP/Goryachev/Holmes C40 d16 | Same nominal ring geometry; each model has its own kinetics, state and native input scale. Includes cue-single and two-site outcomes. |
| Same-state location control | Holmes C40 d12 versus d16 | Identical full basal, circumference, cue strength, width and duration; only cue location changes. Merger versus two sites. |
| Explicit Holmes cue-single example | Holmes C30 d12 | Documents the missing fate; comparison with C40 d12 also changes size/material/architecture. It is not a pure separation control. |
| Same angular perturbation at two sizes | Goryachev C40 d20 versus C80 d40 | Antipodal cue, same physical cue width/strength/duration; size, material, routes and angular cue width change together. Replacement versus incumbent return. |
| Cancellation safeguard | Holmes C40 d20 | Same full basal as the C40 location controls, antipodal cue. Two active sites can produce almost no first moment. |

For a compact main illustration, the Goryachev size pair and the Holmes same-size location pair answer distinct questions and tie directly to existing heatmaps. WP C40 d16 or Goryachev C40 d16 supplies an unambiguous cue-site single peak in the same geometry as the Holmes examples. Holmes C30 and the antipodal cancellation case can remain linked diagnostics if space is limited. A new duration sweep is not necessary to establish this first checkpoint.

The ring-baseline comparison to linear Figure 1 must be treated separately: reflecting an end-polarized interval of length \(L\) gives a ring of circumference \(2L\), with doubled whole-incumbent width and amount. The vector should use the actual ring state, with its own basal normalization. A physical width-3 end cue reflects to width 6, whereas the retained ring cue has width 3; therefore baseline correspondence does not establish identical finite-challenge protocols.

## Retained-data check

The [analysis code](../../../../code/mechanisms/global_polarity_vector_20261017.py) reads immutable cue/recovery arrays, verifies exact primary-field continuity at withdrawal, and exports [all vector samples](supporting/data/global_polarity_traces.csv), [summary statistics](supporting/data/global_polarity_summary.csv), [source hashes](supporting/data/global_polarity_sources.csv) and an [audit](supporting/metadata/global_polarity_audit.json). It covers 20,108 samples across eight histories. Analytic nonnegative test profiles verify uniform-background cancellation, known rotation, scaling invariance, opposed-site cancellation and distinct peak structures with identical first moments; maximum error is \(7.5\times10^{-15}\). Existing stored mass-normalized harmonic magnitudes agree where available.

| Retained history | Minimum sampled \(R\) | Final \(R\) | Final orientation | Observed local organization |
|---|---:|---:|---:|---|
| WP C40 d16 | 0.337 | 1.000 | 142.9° | Cue-site single peak |
| Goryachev C40 d16 | 0.306 | 1.000 | 144.0° | Cue-site single peak |
| Holmes C40 d16 | 0.369 | 0.379 | 65.9° | Two retained sites |
| Holmes C40 d12 | 0.728 | 1.000 | 52.7° | Merged/bridged region |
| Holmes C30 d12 | 0.414 | 1.000 | 143.3° | Cue-site single peak |
| Goryachev C40 d20 | 0.00488 | 1.000 | 180° | Replacement through a near-zero resultant |
| Goryachev C80 d40 | 0.728 | 1.000 | 0° | Initial-region activity returns; vector error \(<6\times10^{-14}\) |
| Holmes C40 d20 | 0.000295 | 0.0399 | Withheld | Two sites; cancellation, not inactivity |

For the last case, total active Rac is 33.13 initially and 33.05 finally, while final \(q_1=0.0165\) and \(q_2=0.401\). Its maximum reported angle changes from 180° at the 0.05 gate to approximately 0° at 0.1/0.2. This is precisely why a heading-only summary is inadequate. No fresh scientific critique of this new descriptor has yet been performed.

## Experimental bridge

The analogous measurements are whole-boundary reporter asymmetry, its transient weakening and its retained axis after a localized challenge, alongside incumbent/challenger readouts. Reporter intensity requires suitable background correction and calibration before being interpreted as molecular amount. The model maps arclength to angle on a perfect circle; it does not predict cell shape, protrusion orientation or traction. In a deforming cell, the geometric weighting and coordinate registration must be declared; equal-arclength circle angles are not automatically the actual boundary normals. Association with measured velocity orientation can subsequently be tested experimentally, but the activity vector is not assigned the role of velocity in these figures.
