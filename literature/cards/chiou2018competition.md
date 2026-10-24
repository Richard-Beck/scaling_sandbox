---
key: chiou2018competition
title: Principles that govern competition or co-existence in Rho-GTPase driven polarization
authors: Chiou JG; Ramirez SA; Elston TC; Witelski TP; Schaeffer DG; Lew DJ
year: 2018
venue: PLoS Comput Biol 14:e1006095
doi: 10.1371/journal.pcbi.1006095
url: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5916526/
type: theory
readout: [simulated_polarity]
size_varied: none
kinetics_fixed: na
noise: none
role: [background, precedent]
themes: [competition_coexistence]
sections: [fig2, fig4]
access: abstract
status: stub
checked: 2026-10-23
---

# Chiou et al. 2018: saturation sets how fast polarity peaks compete

**Why it matters to us.** The saturation rule is the established basis for our Fig 2 expansion route and for coexisting domains in Fig 4.

## What they did
Analysis of mass-conserved activator–substrate (two-component reaction-diffusion) models of Rho-GTPase polarity, asking what sets the timescale of competition between clusters.

## Results
- [R] Some such models show rapid winner-takes-all competition; others allow long coexistence. The timescale of competition varies enormously with parameters.
- [R] One factor explains most of that variation: how close the peak active-GTPase concentration comes to a saturation point set by the model parameters.
- [R] Peaks well below saturation compete effectively and resolve quickly to one. Peaks near saturation compete ineffectively and become metastable.
- [R] The authors argue that saturation and its effect on competition are general properties of Rho-GTPase polarity machinery, whatever the feedback mechanism.

## What transfers
A model-independent rule linking peak shape to competition speed.

## Relation to our claims
- Fig 2: the expansion route in Wave-pinning and Holmes, where active concentration sits at a fixed plateau and a challenger must widen rather than grow taller, is the saturated case. Present as extending this rule to size.
- (I) Fig 4: the rule predicts the direction of the multiple-domain trend per model. Models that broaden at a fixed plateau (saturated) should keep secondary domains longer in larger rings; Goryachev, which raises its peak instead, should resolve them quickly. That is what Fig 4E shows.

## Caveats
Abstract only. Check whether the paper discusses domain size or cell size explicitly.
