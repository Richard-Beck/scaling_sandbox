---
key: hubatsch2019threshold
title: A cell size threshold limits cell polarity and asymmetric division potential
authors: Hubatsch L; Peglion F; Reich JD; Rodrigues NTL; Hirani N; Illukkumbura R; Goehring NW
year: 2019
venue: Nat Phys 15:1078-1085
doi: 10.1038/s41567-019-0601-x
url: https://discovery-pp.ucl.ac.uk/id/eprint/10080678
type: experiment
readout: [observed_polarity, simulated_polarity]
size_varied: threshold
kinetics_fixed: yes
noise: none
role: [background, precedent]
themes: [size_and_polarity]
sections: [fig1, discussion, methods]
access: abstract
status: stub
checked: 2026-10-23
---

# Hubatsch et al. 2019: polarity kinetics do not scale with size, so small cells lose polarity

**Why it matters to us.** The best-known demonstration that cell size changes polarity at fixed kinetics. It also gives empirical support for our "no retuning" design.

## What they did
PAR polarity in the C. elegans germ lineage, where cells shrink with each division. Kinetic parameters and pattern length scales measured in cells of different sizes. Several reaction-diffusion polarity models analysed across system size and total pool, including wave-pinning and Goryachev models.

## Results
- [R] PAR protein behaviour fails to scale with cell size.
- [R] In theory, this lack of scaling gives a size threshold below which polarity is destabilized. A threshold existed in every model examined, and it varied linearly with membrane diffusion.
- [R] In empirically constrained models the threshold lies near the size at which germ-lineage cells switch from asymmetric to symmetric division.
- [R] Reducing germ-lineage cell size genetically or physically triggered premature loss of polarity at the predicted thresholds.
- [R] In the associated thesis, measured kinetics showed moderate adaptation to size at most.

## What transfers
- Justification for holding local kinetics fixed when size changes: measured polarity kinetics need not adapt to size.
- The boundary length scale is set by diffusion and reaction rates, not by cell size.

## Relation to our claims
- Background for the whole model design (Methods) and for Fig 1.
- (I) Their result concerns the small-size end, where polarity fails. Ours concerns the large-size end, where established polarity becomes harder to replace. Together they frame size as a control parameter of polarity at fixed kinetics.
- (I) Check that our smallest lengths (L19, C20) sit well above any polarity threshold in each model, or a reviewer may attribute the size effect to proximity to the threshold.

## Caveats
Embryonic PAR system, not migration. Abstract and thesis abstract only.
