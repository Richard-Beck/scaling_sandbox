# Cell-size scaling sandbox

**[Browse the rendered scientific report catalogue](https://richard-beck.github.io/scaling_sandbox/)**

This repository asks how cellular transport, signalling, polarity, and intracellular organization change as cell size changes. The landing page is the primary entry point: it identifies each source model, what was modified, the main result, and whether the analysis is a reproduction, port, extension, or exploratory mechanism.

## What is here

- Literature-derived polarity models from Mori, Holmes, Buttenschön, and the Eroumé implementation of the Marée Cdc42/Rac/Rho–phosphoinositide system.
- Transport studies of oxygen, membrane-to-nucleus signalling, and microtubule-mediated cargo delivery.
- Generic radial scenarios that isolate source, sink, geometry, and capacity scaling assumptions.
- Rendered reports with downloadable plot-level data and explicit limitations.

The compact literature and implementation catalogue is maintained in [`MODEL_REGISTRY.md`](MODEL_REGISTRY.md).

## Scaling question

For a process with source `S` and capacity `K`, the sandbox makes the assumed exponents explicit:

```text
S ∝ V^alpha
K ∝ V^beta
```

Relative loading scales as `V^(alpha - beta)`. For similar cells, surface area scales as `V^(2/3)` and length as `V^(1/3)`. Each report states which quantities are held fixed and which are allowed to scale.

## Repository organization

- `R/` contains reusable geometry, diffusion, model, and plotting functions.
- `scenarios/` defines geometry, boundaries, parameters, and scaling choices.
- `reports/` contains the scientific narratives and analyses.
- `reports/assets/` contains small, versioned report snapshots where a full simulation is too expensive for a page build.
- `docs/` is the GitHub Pages site root.
- `scripts/` contains rendering, data export, and expensive-generation tools.

## Build the public reports

Install R packages once:

```r
install.packages(c("deSolve", "ReacTran", "ggplot2", "gridExtra", "RANN", "rmarkdown", "Matrix", "rootSolve", "digest"))
```

From the repository root:

```text
Rscript scripts/check_reports_do_not_simulate.R
Rscript scripts/export_public_plot_data.R
Rscript scripts/render_public_reports.R
```

GitHub Actions runs the same public build on pushes to `main`. Generated HTML and plot CSVs under `docs/reports/` and `docs/data/` are intentionally untracked.

## WGD / Eroumé spatial snapshot

The WGD report ports the exact nine-PDE reaction system from the public Eroumé et al. VCell model, which descends from Marée et al. It then performs a new 0.5×/1×/2× geometry experiment; that size experiment is not attributed to Eroumé et al.

Public builds never run this expensive PDE solver. They validate and read the committed, checksummed files in `reports/assets/wgd_spatial/v1/`. Regenerate the full checkpoints and publish a new compact snapshot only when the model or scenario changes:

```text
Rscript scripts/generate_wgd_spatial_results.R
Rscript scripts/publish_wgd_snapshot.R
```

Pass `--reuse` to the generation script to retain completed worker checkpoints. Full raw results remain untracked under `reports/output/wgd_spatial_scale/raw/`; a missing or invalid public snapshot stops the build rather than silently rerunning or substituting a simulation.

## Add a scenario

Create a descriptive file in `scenarios/` that returns a list containing its name, geometry, boundary conditions, and parameters. Keep biological and scaling assumptions in the scenario and numerical methods in `R/`. Add a report that states the source model, modification, status, result, and limits, then register it in `MODEL_REGISTRY.md` and the landing page.
