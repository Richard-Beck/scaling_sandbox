---
key: mitterwallner2020gle
title: Non-Markovian data-driven modeling of single-cell motility
authors: Mitterwallner BG; et al.
year: 2020
venue: Phys Rev E 101:032408
doi: 10.1103/PhysRevE.101.032408
url: https://doi.org/10.1103/PhysRevE.101.032408
type: method
readout: [observed_migration, ou_inference]
size_varied: none
kinetics_fixed: na
noise: na
role: [method, background]
themes: [persistence_and_migration_theory]
sections: [fig5, methods]
access: secondary
status: stub
checked: 2026-10-23
---

# Mitterwallner et al. 2020: long persistence from a short memory kernel

**Why it matters to us.** An alternative to the OU model we fit in Fig 5, and a warning that a fitted persistence time need not match any single memory timescale in the cell.

## What they did
Generalized Langevin equation for the 1D tangential velocity of MDA-MB-436 cells on a circular track, with a memory kernel and coloured random force. Because the two cannot be separated from trajectories alone, they fit a substitute kernel.

## Results
- [R] The data give a kernel made of instantaneous friction plus a delayed negative-friction term: Γ(t) = 2a·δ(t) + b·exp(−t/τ_m) with b < 0. The OU or persistent-random-walk model is the b = 0 limit.
- [R] Cell-averaged fit: a = 0.15 min⁻¹, b = −0.023 min⁻², τ_m = 6.02 min, ⟨v²⟩ = 0.14 µm²/min².
- [R] These give two emergent relaxation times, 3.23 min and 164 min, and a long-time diffusivity of 12.7 µm²/min.
- [R] Across individual cells the near-cancellation a ≈ |b|·τ_m holds closely. It is what makes the emergent persistence time long.
- [D] With the exponential kernel the model reduces to a two-variable Markov system, dv/dt = −av − bh + F_R and dh/dt = v − h/τ_m (from the repo note).

## What transfers
- A persistence time of hours can emerge from a kernel that remembers velocity for minutes. Persistence in migration data need not correspond to a slow molecular process.
- The near-cancellation makes the long timescale sensitive to small parameter changes. A modest change in friction balance with size could change persistence a lot.

## Relation to our claims
- (I) Fig 5: our OU velocity-memory fold (about 16) is conditional on the OU form. If the Xuan curves need a non-Markovian kernel, which their unreproduced crossover suggests, the fold may change. Consider a GLE fit as a supplement robustness check.
- (I) Keep OU and GLE estimates as the "conditional OU inference" readout, separate from simulated persistence.

## Caveats
- From the existing repo note, not the paper. Full author list not confirmed.
- The note gives the emergent persistence time as both about 2.7 h and 164 min; these agree.
