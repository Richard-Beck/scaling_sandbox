# Measurements that could test the model predictions in cells

Serves Fig 2F and the observables pre-specified for the SUM-159 WGD+/− pilot (Fig 6). Drafted 2026-10-23.

Each model mechanism should come with the cellular measurement it motivates. This file lists published precedents for those measurements.

## Precedents

| Model quantity | Cellular measurement | Precedent |
|---|---|---|
| Drift of the dominant domain (Fig 4E row 1) | Wander of the polarity patch or leading-edge centroid; lateral edge positions of the leading active region over time | Patch wandering in yeast [@dyer2013wandering] |
| Number of active regions; multiple-domain fraction (Fig 4E row 4) | Protrusive unicity index: 1 over the mean number of simultaneous protrusions | [@vaidziulyte2022persistent] |
| Heading versus domain alignment | Alignment of polarity axis (nucleus–Golgi) with direction of motion | [@vaidziulyte2022persistent] |
| Births, position and fate of secondary domains (proposed Fig 4 panels) | Protrusion origin (split from existing, or de novo), position, duration, maintenance or retraction, commitment | [@bosgraaf2009pseudopodia] |
| Overlapping versus separated turning regimes (Fig 3) | Response to a stimulus straddling the leading-edge boundary versus one separated from it | Overlap dependence in a model [@kryvoruchko2026intrinsic]; optogenetic tools (cards to be written: Wu 2009, de Seze 2025) |
| Mechanical memory across gaps (Fig 2D–F) | Direction kept after a stimulus switch in 1D channels; turnover of a rear structural marker | [@prentice2016memory] |
| Heading persistence | Direction autocorrelation decay time | [@vaidziulyte2022persistent], and standard practice |

## Our interpretation (I)

- **Adopt published indices where they exist.** Unicity and alignment [@vaidziulyte2022persistent] give the pilot values comparable with RPE1 and HeLa, and the same model-inference strategy (measure structure, read off persistence).
- **Classify protrusions by origin.** The split versus de novo distinction [@bosgraaf2009pseudopodia] separates the two Fig 4 channels in cells: edge displacement versus secondary-domain birth. Prediction for WGD+ cells, if our models are right: more de novo births per cell (larger perimeter, Fig 4E) but a lower fraction that commit, and smaller angular change per edge displacement.
- **Record stimulus position relative to the front.** Kryvoruchko & Camley note that comparing models with collision data requires knowing where contact occurred relative to the front [@kryvoruchko2026intrinsic]. The same holds for any photoactivation test of the Fig 3 regimes.
- **Size must be measured per cell.** WGD changes size and copy number together, and published size effects can have either sign (`size_and_polarity.md`). Per-cell area or volume lets the pilot regress persistence on size within each population as well as between them.

## To check

- Imaging cadence needed to resolve edge diffusion and protrusion births in SUM-159, given the timescales in [@vaidziulyte2022persistent] (persistence of hours; trafficking lag about 19 min).
- Whether the SUM-159 pair can carry a Rac or Cdc42 biosensor without changing migration.
- Cards for Wu 2009 and de Seze 2025, from the archived CPM notes.
