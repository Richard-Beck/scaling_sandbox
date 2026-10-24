---
key: kryvoruchko2026intrinsic
title: Intrinsic stochasticity in cell polarity and contact inhibition of locomotion
authors: Kryvoruchko M; Camley BA
year: 2026
venue: arXiv:2604.18851 (preprint)
doi: 10.48550/arXiv.2604.18851
url: https://arxiv.org/abs/2604.18851
type: model
readout: [simulated_polarity]
size_varied: copy_number
kinetics_fixed: yes
noise: intrinsic
role: [precedent]
themes: [noise_driven_domain_motion, repolarization_and_memory, persistence_and_migration_theory]
sections: [fig3, fig4, discussion]
access: full_text
status: read
checked: 2026-10-23
---

# Kryvoruchko & Camley 2026: copy-number noise makes a wave-pinning peak diffuse

**Why it matters to us.** The closest published work to our Fig 4 drift analysis. It measures noise-driven rotational diffusion of a single wave-pinning peak on a ring, but at one fixed size.

## What they did
Particle-based stochastic version of the Mori wave-pinning model on a 1D periodic membrane of length L = 20 µm. Individual Rho GTPases diffuse (membrane 0.1, cytosol 10 µm²/s, cytosol simulated explicitly) and bind or unbind stochastically. Total copy number N varied from 10³ to 10⁵ at fixed geometry. The Hill constant is rescaled with N so the large-N limit recovers deterministic wave-pinning. Heading is the vector sum of membrane-bound positions on the ring, the same construction as our P = M1/M0. Contact inhibition of locomotion is modelled as a local, transient rise in the off rate.

## Results
- [R] Polarity depends strongly on N. Near N ≈ 2×10³ no stable peak persists and multiple transient peaks appear. Near 5×10³ a peak forms and wanders visibly. Above about 10⁵ the peak is nearly stationary.
- [R] The peak's centre of mass performs an unbiased random walk. Angular diffusion D_P, fitted from angular MSD over 0–100 s, falls with N: about 3.8×10⁻⁴ to 1.4×10⁻⁴ rad²/s for N between 1.5×10⁴ and 4×10⁴.
- [R] Power-law fit over N = 5×10³–2.5×10⁴ gives D_P ∝ N^−1.14±0.02. Restricted to N > 10⁴, the exponent is 1.03 ± 0.03, consistent with 1/N.
- [R] An independent-molecule argument (their Appendix B) gives D_P = D_φ/N_m, with N_m the number of membrane-bound molecules in the peak.
- [R] Persistence time taken as 1/D_P is about 0.7 h at N = 1.5×10⁴. Near N = 2.6×10⁴ it matches a fibroblast estimate of about 1.2 h.
- [R] Reorientation after a contact falls smoothly as the contact moves from the peak towards the rear. It becomes negligible once the inhibited region no longer overlaps the polarized domain. Wider contacts increase reorientation up to a plateau when the contact spans about half the membrane.
- [R] The mean stochastic reorientation is slightly larger than the deterministic one and converges to it as N → ∞.
- [R] An inhibitor that spreads by diffusion beyond the contact (their Appendix D) gives reorientation even from contacts at the rear, and drives near-180° reversals.
- [R] The authors estimate physiological Rac1 and Cdc42 copy numbers near 10⁶, above their simulated range.

## What transfers
- The heading definition and the fitting approach match ours.
- Size is never varied. They hold L fixed and change N, which changes concentration. They note that Walther et al. 2012 instead varied N together with cell width to keep concentration fixed.
- Noise is intrinsic only. Our Fig 4 input noise is extrinsic.

## Relation to our claims
- Precedent for Fig 4: domain drift as the source of heading loss, and τ ≈ 1/D as the bridge to migration persistence. Cite as the nearest neighbour.
- (I) Their 1/N_m scaling suggests a third null for Fig 4. If a broadening domain's content grows in proportion to circumference, content-averaged noise gives D_s ∝ 1/C and α = 3. Edge-limited drift gives α = 2. See `syntheses/noise_driven_domain_motion.md`.
- (I) Their overlap dependence is a single-size precedent for the two Fig 3 regimes: input that overlaps the incumbent displaces it; input away from it does little unless it nucleates a competitor.

## Caveats
- Preprint. One model. D_P fitted over a 100 s window only.
- Their persistence numbers rest on copy numbers below their own physiological estimate.
