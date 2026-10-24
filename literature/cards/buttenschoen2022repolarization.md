---
key: buttenschoen2022repolarization
title: "Cell repolarization: a bifurcation study of spatio-temporal perturbations of polar cells"
authors: Buttenschön A; Edelstein-Keshet L
year: 2022
venue: Bull Math Biol 84:114
doi: 10.1007/s11538-022-01053-z
url: https://link.springer.com/article/10.1007/s11538-022-01053-z
type: model
readout: [simulated_polarity]
size_varied: none
kinetics_fixed: na
noise: none
role: [precedent]
themes: [repolarization_and_memory]
sections: [fig1, fig2, fig3]
access: abstract
status: stub
checked: 2026-10-23
---

# Buttenschön & Edelstein-Keshet 2022: routes to reversing an established polarity

**Why it matters to us.** The closest precedent for our Figs 1–2: how an already polarized cell responds to a reversing stimulus, analysed across model classes.

## What they did
PDE bifurcation analysis of two polarity models under spatio-temporal perturbation: the single-GTPase wave-pinning model and a mutually antagonistic Rac–Rho model. Compared with full PDE simulations, and with simulations in deforming 1D model cells.

## Results
- [R] The sensitivity of a polar cell to a reversal stimulus was not well understood before this work.
- [R] The two models reverse by distinct routes.
- [R] Full simulations agree with the bifurcation predictions, and deforming-cell simulations are consistent with experiments.

## What transfers
The bifurcation approach and the finding that reversal mechanism is model-class dependent.

## Relation to our claims
- Precedent for Figs 1–2. Our three routes (initiation, expansion, mechanical memory) extend the idea that reversal mechanism depends on model class. What we add is how each route scales with size.
- (I) Their routes may correspond to our two Fig 3 regimes (displacing the incumbent versus nucleating a challenger). Check once read.

## Caveats
Paywalled; abstract only. Check whether domain length is varied. Their deforming-cell code builds on Buttenschön, Liu & Edelstein-Keshet 2020, which the repo has implemented (`literature_models/buttenschoen_2020_moving_domain/`) and which is top of `queue.md`.
