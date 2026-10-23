# Next revision: Figures 3 and 4

Oct 23, 2026 · @Richard

## Context and rationale

Figure 4 is the weakest figure. The nulls helped locate where the models depart from simple expectations, but they don't explain the size effect, and the figure doesn't connect to Figures 2–3 or suggest anything to measure in cells. Two findings from the round-3 data give it a clearer shape.

- **The weak-input regime behaves like a small-signal response of the incumbent's two edges.** For f ≤ 0.01, drift diffusivity scales as f^1.94–1.99 at every size in Wave-pinning, Holmes and Spring. The domain's MSD is close to diffusive (exponent about 1.1), and multiple domains are essentially absent (≤ 0.1% of frames in Wave-pinning; up to 7% in Holmes and 13% in Spring, only at f = 0.01). If the two edges respond independently to local input, drift is set by how far an edge moves per unit of input at that edge. That response can be measured directly, without noise, and belongs with Figure 3's other cue responses.
- **Beyond the drift regime, larger rings start many more secondary domains but rarely lose the incumbent.** In Holmes at f = 0.03, the multiple-front fraction rises from 7% at C20 to 69% at C60, while counted takeovers only double. This is the asymmetry of Figure 2 (initiation independent of length, replacement 5–12× slower) appearing under noise. Births, takeovers and the heading loss they cause need to be measured directly; the multiple-front fraction is not enough to see what is happening.

The aim is a Figure 3 that ends with a prediction for drift, and a slimmer Figure 4 that tests it and then shows what happens once secondary domains appear. Two of the new measurements have direct experimental counterparts, and these become pre-specified pilot observables.

## Figure 3: add the edge response and a drift prediction

Panels A (the two turning regimes) and B (half-turn time against cue gap) stay. Panel C (mean angular speed) is a restatement of B and goes. Two panels replace it.

| Panel | Content | Source |
| --- | --- | --- |
| A | Two turning regimes, Wave-pinning C = 30 (kymographs optional) | Existing |
| B | Half-turn time against cue gap, per model and size | Existing |
| C | Edge response to a weak, brief pulse at the incumbent's edge | New deterministic runs |
| D | Drift predicted from C against drift measured in the Figure 4 campaign | New analysis |

1. **Run the edge pulse.**
   - Models: Wave-pinning, Holmes, Spring, Goryachev. No retuning.
   - Sizes: the Figure 4 rings (C20, 30, 40, 50, 60; Spring 1–3), since D compares against Figure 4 drift. Also run the Figure 3 sizes (30, 37.8, 47.6; Spring 2.5, 3.15, 3.97) if C is to match B.
   - Start from each model's basal polarized state at the Figure 4 mean input for the chosen f, with no noise.
   - Pulse: centred on the incumbent's half-contrast boundary, with width and duration matched to the noise (ℓ = 1.5, τ = 50 chemical; ℓ = 0.15, τ = 5 Spring).
   - Amplitude small enough to be linear. Check by halving and doubling it: the response must scale proportionally. Report the amplitude used.
   - Track both edges (half-contrast crossings) through the response until they relax.
2. **Plot C.** For each model, plot pulsed-edge displacement against time, one line per size on the model's light-to-dark ramp. Summarize the response as integrated displacement per unit input against circumference, in an inset or adjacent subpanel.
3. **Compute D.**
   - Predict edge diffusivity from the measured response kernel and the declared noise covariance (σ, ℓ, τ, f) under linear response.
   - Take the domain's drift as D\_s = (D\_left + D\_right)/4, assuming independent edges.
   - Plot predicted against measured D\_s, using the Figure 4 tracked drift at f ≤ 0.01, with a 1:1 line, coloured by model and shaded by size.
   - Disagreement is a result. Don't adjust anything to close it.
4. **Check edge independence (supplement).** From the Figure 4 exemplar profiles, track both edges and report the cross-correlation of their displacements and the edge MSDs. A positive correlation means the shared pool couples the edges, and the factor of 4 above becomes smaller.
5. **Text link.** The legend for D says it predicts the drift measured in Figure 4. C's experimental counterpart is photoactivation straddling the leading-edge boundary.

## Figure 4: slimmer, with births and takeovers

| Panel | Content | Source |
| --- | --- | --- |
| A | Null kymograph pairs and their C(ℓ) | Current A |
| B | Model α against input, plus noise-tolerance fold | Current C |
| C | Model kymographs | Current D, fixed |
| D | Drift regime: normalized drift against input, and predicted/observed decay in that regime | Current E rows 1 and 3, trimmed |
| E | Beyond drift: births, takeovers and the source of heading loss | New analysis of compact records |

1. **A.** Keep the null kymographs and C(ℓ). Collapse current B's three α panels into one, or move them to the supplement. Keep the tolerance-fold panel next to B.
2. **C.** Centre every kymograph on its initial heading so a domain that wraps the seam doesn't read as two. Give Goryachev a per-size activity scale (C20 max 9.2, C60 max 83) and say so in the legend.
3. **D.**
   - Plot D\_s/D\_s(C\_min) against f, as now.
   - State the fitted f-exponent in the drift regime (about 2 for Wave-pinning, Holmes and Spring) in the legend, not on the plot.
   - Show predicted/observed decay only where multiple domains are rare (state the cutoff used).
4. **E(i), births.** Plot birth rate of secondary domains per unit perimeter per unit time against f, one line per size.
   - A birth is a tracked region that appears without overlapping any region in the previous frame and persists at least 5 frames.
   - Count splits separately and report them in the supplement.
5. **E(ii), takeovers.** Plot the probability that a birth becomes dominant, against size, binned by its distance from the incumbent at birth.
   - Use centroid-to-centroid arc distance minus the incumbent's half-width, with half-width estimated from amount over plateau height.
   - Use the same 5-frame hold as the switch counts.
6. **E(iii), source of heading loss.** Write P = Σᵢ Aᵢ e^{iθᵢ} / Σᵢ Aᵢ over tracked regions, plus a background residual.
   - Over the condition's analysis lag, split the heading change into domain motion at fixed amounts, shifts in amount between coexisting domains at fixed positions, births and losses, and the residual.
   - Plot the fractions against f, by size.
7. **E text link.** The legend points back to Figure 2B–C: challengers start readily at any size, but incumbents become harder to remove. No separate panel for this.
8. **Move to the supplement:** the λ\_eff row, the detector-free fit and its disagreement, and the multiple-front fraction row (superseded by E(i)).

The usual style rules apply: no titles on the figure, black panel labels, plain facet labels, symbols defined beside their panels.

## Pilot observables, fixed before the SUM-159 data arrive

Each maps onto a model quantity in Figures 3–4. Writing the expected directions down now stops the "strongly suggestive" criterion being judged after the fact.

| Observable in cells | Model counterpart | Drift regime: WGD+ vs WGD− | Secondary-domain regime: WGD+ vs WGD− |
| --- | --- | --- | --- |
| Lateral-edge fluctuation of the leading protrusive region | Edge diffusivity (Fig. 3D, 4D) | Similar | Larger |
| Secondary-protrusion birth rate per unit perimeter | Birth rate (Fig. 4E(i)) | Rare in both | Higher |
| Fraction of secondary protrusions that take over | Takeover probability (Fig. 4E(ii)) | Not applicable | Lower |
| Heading autocorrelation time | Heading decay (Fig. 4B, C) | Longer, roughly by size squared | Advantage shrinks |
| Edge displacement after photoactivation straddling the leading-edge boundary | Edge response (Fig. 3C) | Model-dependent: similar (Holmes, Spring) or smaller (Wave-pinning) | Not applicable |

These directions come from the models as they stand; update them if the new panels change them.

## Scope and data

- **In scope:** Figure 3 C–D, Figure 4 as above, the pilot observables table.
- **Out of scope this round:** Figure 1, Figure 2 (including 2F), the magnitude comparison with the published folds, and the literature check on geometric persistence.

| Work | Needs |
| --- | --- |
| Fig. 4 E(i)–(iii), Fig. 4D trims | Round-3 compact records (per-region centroid, amount, track identity per frame), on HPC |
| Edge-independence check | Round-3 exemplar profiles (seeds 0 and 1), on HPC |
| Fig. 3C edge pulse | New deterministic runs; cheap, no noise |
| Fig. 3D | Fig. 3C output plus the existing Fig. 4 drift table |

No new noise campaign. If an analysis breaks or a measure is undefined for some conditions, build the figure anyway, plot what exists, and list the cases in the readme.
