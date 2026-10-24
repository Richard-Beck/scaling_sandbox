---
key: guan2023relocation
title: Particle-based simulations reveal two positive feedback loops allow relocation and stabilization of the polarity site during yeast mating
authors: Guan K; Curtis ER; Lew DJ; Elston TC
year: 2023
venue: PLoS Comput Biol 19:e1011523
doi: 10.1371/journal.pcbi.1011523
url: https://doi.org/10.1371/journal.pcbi.1011523
type: model
readout: [simulated_polarity]
size_varied: none
kinetics_fixed: na
noise: intrinsic
role: [counterpoint]
themes: [noise_driven_domain_motion]
sections: [fig4, discussion]
access: abstract
status: stub
checked: 2026-10-23
---

# Guan et al. 2023: molecular noise alone may not relocate a strongly self-reinforcing site

**Why it matters to us.** A counterpoint. In their core circuit, intrinsic noise is too weak to move the site; relocation needs an extra feedback loop.

## What they did
Particle-based simulations of the core yeast polarity circuit, with and without a second pathway that recruits pheromone receptors to the membrane.

## Results
- [R] Molecular-level fluctuations in the core circuit are unlikely to overcome the strong positive feedback needed for polarization, so they do not produce relocating polarity sites.
- [R] Adding a second positive feedback loop, receptor recruitment coupled to gradient sensing, produced relocating sites resembling those observed. The same loop lets cells stabilize the site once it aligns with the gradient.

## What transfers
Whether noise drives relocation depends on circuit structure, not only on noise strength.

## Relation to our claims
- (I) Our Fig 4 relocation (jump) channel uses extrinsic input noise, which can be made strong. A reviewer may ask whether realistic intrinsic noise would ever reach the relocation regime. Our answer should be that f is a scan over input strength, not a claim about intrinsic noise levels.

## Caveats
Yeast mating; abstract only.
