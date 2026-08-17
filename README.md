# Cell-size scaling sandbox

This sandbox is a small, reproducible workspace for asking how cellular
processes change as cell size changes. It is intentionally organized around a
separation of concerns:

- `R/` contains reusable geometry, diffusion, scaling, and plotting functions.
- `scenarios/` defines named model instances: their geometry, boundary
  conditions, sources, sinks, and scaling assumptions.
- `reports/` contains R Markdown analyses that load one or more scenarios and
  present their results.

## Core question

For a process with source `S` and capacity `K`, we make its scaling assumption
explicit:

```
S \u221d V^alpha
K \u221d V^beta
```

The relative loading therefore scales as `V^(alpha - beta)`. Geometry can add
an independent transport constraint: for geometrically similar spherical
cells, surface area scales as `V^(2/3)` and radius as `V^(1/3)`.

## Quick start

Open `reports/size_scaling_overview.Rmd` in RStudio and knit it. The spatial
models use `ReacTran` for finite-volume diffusion and `deSolve` for adaptive
method-of-lines integration. Install them once if necessary:

```r
install.packages(c("deSolve", "ReacTran", "ggplot2", "rmarkdown"))
```

To work interactively, source the reusable functions and then a scenario:

```r
source("R/scaling.R")
source("R/diffusion.R")
source("R/plotting.R")
source("scenarios/external_supply.R")

scenario <- external_supply_scenario()
```

## Adding a scenario

Create a descriptive file in `scenarios/`. It should return a list with a
name, geometry, boundary condition, and a `parameters` list. Keep numerical
choices and biological assumptions there; keep numerical methods in `R/`.

Generated HTML, figures, and other rendered material should go in
`reports/output/`, which is excluded from version control.

## Reports

- `reports/size_scaling_overview.Rmd` documents the exploratory radial
  boundary-flux scenarios and their limitations.
- `reports/sedlack_single_cell_oxygen.Rmd` runs a Sedlack-inspired 3D Cartesian
  oxygen model with a flattened adherent cell, finite-permeability membrane,
  uniform cytoplasmic sink, one-time calibration, and 1x/2x volume comparison.
- `reports/membrane_to_nucleus_signaling.Rmd` tests whether distributed PDE
  hydrolysis redirects a membrane-generated cAMP-like signal away from a
  localized, partially absorbing nuclear target as cell volume increases.

## Publishing

GitHub Actions builds the Pages artifact from `docs/` on each push to `main`.
The workflow regenerates the report HTML and plot-level CSVs, but those derived
files under `docs/data/` and `docs/reports/` are deliberately untracked.
