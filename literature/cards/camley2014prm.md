---
key: camley2014prm
title: Polarity mechanisms such as contact inhibition of locomotion regulate persistent rotational motion of mammalian cells on micropatterns
authors: Camley BA; Zhang Y; Zhao Y; Li B; Ben-Jacob E; Levine H; Rappel WJ
year: 2014
venue: Proc Natl Acad Sci USA 111:14770-14775
doi: 10.1073/pnas.1414498111
url: https://pmc.ncbi.nlm.nih.gov/articles/PMC4205601
type: model
readout: [simulated_polarity, observed_migration]
size_varied: none
kinetics_fixed: na
noise: extrinsic
role: [method]
themes: [persistence_and_migration_theory]
sections: [fig4, methods]
access: abstract
status: stub
checked: 2026-10-23
---

# Camley et al. 2014: persistence set by a fluctuating inhibitor field

**Why it matters to us.** A widely used way to give a reaction-diffusion polarity model a tunable persistence: drive it with a fluctuating inhibitory field. Close in spirit to our fluctuating input.

## What they did
Phase-field model of pairs of mammalian cells on adhesive micropatterns. Reaction-diffusion polarity inside each cell, coupled to adhesive and repulsive forces between cells. A fluctuating inhibitory field that deactivates the polarity species controls each cell's persistence. Four polarity mechanisms compared: contact inhibition of locomotion, neighbour alignment, velocity alignment and a default with none.

## Results
- [R] The fluctuating inhibitor controls persistence of motion.
- [R] Polarity mechanisms strongly regulate persistent rotational motion of cell pairs; small differences in mechanism give large differences in collective rotation.

## What transfers
Precedent for treating extrinsic polarity noise as a tunable input whose strength sets persistence, as our f does.

## Relation to our claims
- Method background for Fig 4. Kryvoruchko & Camley [@kryvoruchko2026intrinsic] cite it for the polarity-vector diffusion coefficient.

## Caveats
Abstract only. Cell size not varied.
