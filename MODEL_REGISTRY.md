# Model registry

This catalogue records the scientific lineage and implementation boundary for the public reports. “Extension” always means an experiment or assumption added in this repository, not a result claimed by the source paper.

| Report / model | Literature lineage | What is implemented | What is new here | Status |
|---|---|---|---|---|
| Eroumé–Marée Rac/Cdc42/Rho polarity | Marée et al. (2012); Eroumé et al. (2021) | Exact reaction system parsed from the public Eroumé VCML fixed-shape nine-PDE model: active/inactive Cdc42, Rac, Rho plus PIP/PIP2/PIP3 | Matched 0.5×/1×/2× implied-volume experiment in two teardrop aspect ratios | **Published VCell model ported + WGD-oriented size extension** |
| Holmes Rac/Cdc42 size scaling | Holmes et al. (2012) | Rac/Cdc42/Rho circuit and membrane-accessibility reduction | Geometrically similar enlargement and response decomposition | **Published model adapted + size extension** |
| Buttenschön cell size and Rac polarity | Buttenschön, Liu & Edelstein-Keshet (2020) | Published moving-domain variants with size-dependent diffusion and dilution | WGD-oriented size-versus-dilution comparison | **Partial reproduction + extension** |
| Wave-pinning polarity | Mori, Jilkine & Edelstein-Keshet | Minimal mass-conserved wave-pinning reaction–diffusion model | Domain-size sweep at fixed biochemical assumptions | **Published minimal model reimplemented + extension** |
| Single-cell oxygen | Sedlack et al. (2022) | 3D medium/cell oxygen transport and kinetic targets, with homogeneous cellular consumption | Uniform-sink simplification and 1×/2× geometry experiment | **Literature-derived simplification + extension** |
| Membrane-to-nucleus signalling | Feinstein et al. (2012) and later finite-element reproduction | PMVEC-inspired geometry, cAMP transport, membrane source, and PDE baseline | Abstract partially absorbing nucleus and size/sensitivity sweeps | **Literature-derived baseline + abstract extension** |
| Microtubule networks and cargo | Shariff, Murphy & Rohde (2010); Li et al. (2012), with junction motivation from Gao et al. and Chen et al. | Persistent constrained growth and rebound/regeneration components | Density pruning, 3D junction detection, illustrative cargo rules, and 1×/2× experiment | **Components adapted + custom transport assumptions** |
| Imposed-flux radial scenarios | Generic diffusion–reaction principles | External supply and nuclear-export radial scenarios | Alternative source/sink scaling assumptions | **Exploratory / generic; not a literature reproduction** |
| De Belly mechanochemical Rac/Rho polarity | De Belly et al. (2026) | Nothing currently implemented | Candidate alternative: local Rac/Rho inhibition plus long-range mutual activation through protrusion/tension and contraction/cortical flow | **Considered/evaluated; not implemented** |

## Key provenance records

Marée AFM, Grieneisen VA, Edelstein-Keshet L (2012). How Cells Integrate Complex Stimuli: The Effect of Feedback from Phosphoinositides and Cell Shape on Cell Polarization and Motility. *PLoS Computational Biology* 8(3): e1002402. [doi:10.1371/journal.pcbi.1002402](https://doi.org/10.1371/journal.pcbi.1002402). Foundational parent model; not separately reimplemented.

Eroumé K, Vasilevich A, Vermeulen S, de Boer J, Carlier A (2021). On the influence of cell shape on dynamic reaction-diffusion polarization patterns. *PLoS ONE* 16(3): e0248293. [doi:10.1371/journal.pone.0248293](https://doi.org/10.1371/journal.pone.0248293). The public VCell model, including `Kerbai_PLoSone_2021_teardrop_polarization_extended`, is the direct source for the port. Eroumé et al. omit the cytoskeletal module because shapes are fixed and retain the Marée parameter values. The size experiment in this repository is new.

De Belly H, Gallén AF, Strickland E, et al. (2026). Long-range mutual activation establishes Rho and Rac polarity during cell migration. *Nature Cell Biology* 28:1244–1257. [doi:10.1038/s41556-026-01965-1](https://doi.org/10.1038/s41556-026-01965-1). This experimental and minimal mechanochemical model is retained as an unimplemented alternative, not presented as part of the Eroumé port.
