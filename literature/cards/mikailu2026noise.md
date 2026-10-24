---
key: mikailu2026noise
title: "Stochastic models for cell polarization and crawling: journey from external fluctuations in reaction-diffusion system to internal noise effects"
authors: Mikailu BB; Montañola-Sales C; Ramírez-Piscina L; Alonso S
year: 2026
venue: Eur Biophys J
doi: 10.1007/s00249-026-01858-z
url: https://link.springer.com/article/10.1007/s00249-026-01858-z
type: model
readout: [simulated_polarity, observed_migration]
size_varied: none
kinetics_fixed: na
noise: both
role: [background, method]
themes: [noise_driven_domain_motion]
sections: [fig4, methods]
access: full_text
status: read
checked: 2026-10-23
---

# Mikailu et al. 2026: external and internal noise are not interchangeable

**Why it matters to us.** Shows that the choice of noise model changes polarity dynamics qualitatively. Our Fig 4 uses one extrinsic noise model, so this bears on how far our results generalize.

## What they did
A bistable (Schlögl-type) membrane reaction with a global feedback that holds total coverage near a target, coupled to a phase field for cell shape, previously calibrated to Dictyostelium. Noise is either external (Ornstein–Uhlenbeck, correlation time about 10 s, at the boundary or in the bulk) or internal, from a chemical Langevin equation for finite molecule numbers. Also simulated without the phase field on a periodic square.

## Results
- [R] With external noise, a larger reaction coefficient k strengthens the deterministic mechanism against fixed noise, giving faster, more persistent crawling.
- [R] With internal noise, noise terms also scale with the reaction rate (as √k), so raising k makes patterns noisier and cells slower. The dependence on k is inverted.
- [R] Temporally correlated external noise generates small transient patches, which the authors treat as pseudopod precursors. White noise of the same variance, and internal noise, do not produce them.
- [R] At high noise the initial circular domain dissolves to a homogeneous state. The transition is sharper and less variable without temporal correlation.

## What transfers
- Temporal correlation of the input, not only its variance, controls whether transient secondary patches appear.
- The authors conclude that intrinsic noise can drive polarity and motility but probably lacks the correlation needed for multiple transient pseudopodia.

## Relation to our claims
- (I) Fig 4 sub-threshold satellites and secondary domains near input hot spots may depend on the correlation time of our input. State the input's spatial and temporal correlation in the Fig 4 methods, and consider a correlation-time check in the supplement.
- (I) Background for the intrinsic-versus-extrinsic distinction raised in `syntheses/noise_driven_domain_motion.md`.

## Caveats
2D, single bistable species, no explicit competition analysis, no size variation.
