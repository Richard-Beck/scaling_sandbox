---
key: rocco2001front
title: Diffusion coefficient of propagating fronts with multiplicative noise
authors: Rocco A; Casademunt J; Ebert U; van Saarloos W
year: 2001
venue: Phys Rev E 65:012102
doi: 10.1103/PhysRevE.65.012102
url: https://arxiv.org/abs/cond-mat/0106081
type: theory
readout: [none]
size_varied: none
kinetics_fixed: na
noise: extrinsic
role: [background, method]
themes: [noise_driven_domain_motion]
sections: [fig3, fig4, methods]
access: abstract
status: stub
checked: 2026-10-23
---

# Rocco et al. 2001: how noise makes a front diffuse

**Why it matters to us.** The theory behind treating domain drift as a small-signal response. Our edge-pulse → drift prediction is an application of this framework.

## What they did
Small-noise expansion for reaction-diffusion fronts under multiplicative noise. They justify an earlier multiple-scale result (Armero et al.) for the front's effective diffusion coefficient.

## Results
- [R] Bistable fronts (between a stable and a metastable state) and pushed fronts wander diffusively about their mean position under noise. Pulled fronts are subdiffusive.
- [R] The derivation splits front motion into diffusion of a suitably defined front position plus fluctuations about the mean profile. The position is defined so that those fluctuations are orthogonal to the left translation mode.
- [R] The authors present the method as a general strategy for fluctuations of fronts and other coherent structures.

## What transfers
- Front diffusion is a projection of noise onto the translation (Goldstone) mode. At small noise it scales with noise variance.
- Our domains are pinned, not propagating, and in mass-conserved models the two edges couple through the cytosol. The framework needs adapting for that.

## Relation to our claims
- (I) Background for Fig 4: drift diffusivity ∝ f^1.94–1.99 and a near-diffusive MSD in Wave-pinning, Holmes and Spring are what this theory predicts. Goryachev's MSD exponent near 0.75 is not.
- (I) Method for Fig 3C–D: the edge-pulse response is the response function onto the translation mode. Present the drift prediction as applying this known approach to size scaling, not as a new idea.

## Caveats
- Abstract and introduction only. Armero et al. (the original result) is in `queue.md`; its citation is not yet confirmed.
- Year: issue dated January 2002, posted 2001. Confirm before citing.
