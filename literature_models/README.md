# Literature models

Each model has the same top-level shape: `src/`, `scripts/`, `tests/`,
`params-or-config/`, and `results/`, plus a provenance README. This first pass
preserves known source and reproduction tooling without claiming that every
legacy workflow is already path-independent in its new location.

Model-specific migration1d adapters and campaigns live in `../migration1d/`;
the folders here preserve standalone/native implementations.

| Folder | Model lineage | Restored in this pass |
|---|---|---|
| `mori_2008_wave_pinning` | Mori, Jilkine & Edelstein-Keshet | report; migration implementation is in `migration1d` |
| `holmes_2012_rac_cdc42_rho` | Holmes et al. | report; migration implementation is in `migration1d` |
| `taniguchi_2013_excitable_waves` | Taniguchi et al. | Python source, scripts, SLURM files, tests |
| `chiou_2015_mcas` | Chiou et al. | R source, initial state, scripts, report |
| `giese_2015_bulk_surface_wave_pinning` | Giese et al. | bulk-surface/Giese source, scripts, tests |
| `buttenschoen_2020_moving_domain` | Buttenschön et al. | Python source, scripts, tests |
| `eroume_2021_maree_polarity` | Eroumé/Marée | Python source, scripts, tests |
| `debelly_2026_mechanochemical` | De Belly et al. 2026 mechanochemical polarity model (mechanics inherited from 2023) | Python and Julia source, scripts, SLURM files, tests |
| `M2072` | Algorta et al., eukaryotic motility driven by actin waves | Retired branch: compact historical record, original V1–V6 and selected leader parameters, and two self-contained Morpheus XML model examples; campaign infrastructure and bulk outputs removed under the 2026-10-17 cleanup. |

The shared, performance-selected WP1 adapters for the Jilkine repertoire and
De Belly 1D now live in `../polarity1d/`. The per-publication folders here
remain provenance/native-reproduction bundles rather than a screening API.
