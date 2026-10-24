---
key: dyer2013wandering
title: Tracking shallow chemical gradients by actin-driven wandering of the polarization site
authors: Dyer JM; Savage NS; Jin M; Zyla TR; Elston TC; Lew DJ
year: 2013
venue: Curr Biol 23:32-41
doi: 10.1016/j.cub.2012.11.014
url: https://pmc.ncbi.nlm.nih.gov/articles/PMC3543483
type: experiment
readout: [observed_polarity]
size_varied: none
kinetics_fixed: na
noise: na
role: [background, method]
themes: [noise_driven_domain_motion, measurement_methods]
sections: [fig4, fig6]
access: abstract
status: stub
checked: 2026-10-23
---

# Dyer et al. 2013: the yeast polarity patch wanders, and the wandering is useful

**Why it matters to us.** Direct observation that a polarity domain performs a biased random walk in living cells, with a mechanism and a functional role.

## What they did
Budding yeast responding to sub-saturating or graded mating pheromone. Imaging of the Cdc42 polarity patch, combined with modelling and genetics.

## Results
- [R] In sub-saturating pheromone the polarity patch wanders around the cortex. Cells tracking a gradient in microfluidics also show a wandering patch, mostly on the up-gradient side.
- [R] Modelling and genetics suggest actin-directed membrane traffic drives wandering by diluting local polarity factors.
- [R] Receptor-activated Gβγ interacting with polarity regulators appears to bias the wandering. Artificially blocking wandering impairs gradient tracking.
- [R] The authors propose that the patch undergoes an intracellular biased random walk, allowing noise filtering by time averaging.

## What transfers
- A cellular precedent for measuring domain-centroid wandering, which our Fig 4 drift readout proposes for SUM-159.
- The driver is vesicle traffic, not chemical noise. Our models have no traffic.

## Relation to our claims
- Background for Fig 4 and for the proposed edge-diffusion observable.
- (I) A non-chemical driver of wandering could scale with size differently from our input noise. Worth one sentence in the discussion.

## Caveats
Yeast, non-motile. Abstract only.
