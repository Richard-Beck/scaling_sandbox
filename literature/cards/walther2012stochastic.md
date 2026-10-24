---
key: walther2012stochastic
title: Deterministic versus stochastic cell polarisation through wave-pinning
authors: Walther GR; Marée AFM; Edelstein-Keshet L; Grieneisen VA
year: 2012
venue: Bull Math Biol 74:2570-2599
doi: 10.1007/s11538-012-9766-5
url: https://pmc.ncbi.nlm.nih.gov/articles/PMC3480592
type: model
readout: [simulated_polarity]
size_varied: copy_number
kinetics_fixed: yes
noise: intrinsic
role: [precedent]
themes: [noise_driven_domain_motion, size_and_polarity]
sections: [fig4]
access: abstract
status: stub
checked: 2026-10-23
---

# Walther et al. 2012: stochastic wave-pinning at low molecule number

**Why it matters to us.** The first stochastic treatment of wave-pinning, and the origin of the "pinning position fluctuates" picture behind our drift analysis.

## What they did
Gillespie-type stochastic version of the single-GTPase wave-pinning model in 1D. They use local perturbation analysis to predict where polarity initiates and where wave-pinning holds, then compare deterministic and stochastic simulations.

## Results
- [R] They determine a threshold number of molecules needed for robust polarization in a given reaction volume.
- [R] As molecule number falls, wave-pinning is lost. The transition zone widens and the pinning position fluctuates more, until the wave broadens beyond what can be sustained and collapses. The high and low plateau levels are much less affected.

## What transfers
- The qualitative picture: noise moves the domain boundaries, and the boundary position is the soft variable.
- Per Kryvoruchko & Camley [@kryvoruchko2026intrinsic], this study varied molecule number together with cell width to hold concentration fixed. If so, it varies size at fixed concentration, which is close to our design. Not yet confirmed from the full text.

## Relation to our claims
- Background for Fig 4 drift in Wave-pinning.
- (I) If the width variation is confirmed, check whether they report pinning-position variance against width. That would be direct precedent for our size-dependent drift.

## Caveats
Abstract only. Full text needed for the size design and any quantitative fluctuation results.
