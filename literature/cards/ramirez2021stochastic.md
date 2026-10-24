---
key: ramirez2021stochastic
title: A novel stochastic simulation approach enables exploration of mechanisms for regulating polarity site movement
authors: Ramirez SA; Pablo M; Burk S; Lew DJ; Elston TC
year: 2021
venue: PLoS Comput Biol 17:e1008525
doi: 10.1371/journal.pcbi.1008525
url: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8315557/
type: model
readout: [simulated_polarity]
size_varied: none
kinetics_fixed: na
noise: intrinsic
role: [precedent]
themes: [noise_driven_domain_motion]
sections: [fig4]
access: abstract
status: stub
checked: 2026-10-23
---

# Ramirez et al. 2021: when molecular noise can move a polarity site

**Why it matters to us.** States the conditions under which intrinsic noise relocates a polarity cluster, which bears on what sets drift in larger domains.

## What they did
A reaction-diffusion master equation method for mechanistic yeast Cdc42 polarity models, used to study how molecular fluctuations move the polarity site. Key results checked with full 3D particle-based simulations.

## Results
- [R] Molecular noise relocated clusters effectively when clusters held low numbers of limiting polarity factors and when Cdc42 dwell times at the site were short.
- [R] Clusters were stabilized by changes in abundance or binding rates that lengthened dwell times or increased the number of polarity molecules in the cluster.

## What transfers
The principle that more molecules per cluster means less noise-driven movement.

## Relation to our claims
- Background for Fig 4.
- (I) Consistent with the content-averaging null in `syntheses/noise_driven_domain_motion.md`: models that broaden their domain with size put more material in it, which should damp intrinsic-noise drift.

## Caveats
Yeast models; intrinsic noise; abstract only.
