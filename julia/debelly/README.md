# De Belly mechanochemical-polarization reproductions

This directory is a self-contained Julia/Gridap implementation adapted from
the authors' CC BY-SA 4.0 release. It does not read from or refer to `dev/`.
See `LICENSE` for attribution and redistribution terms.

It contains the published 1D model, its baseline and front-to-back
optogenetic protocols (mechanochemical and local-inhibition controls), and the
2D axisymmetric surface model with its representative spherical-cell example.
The locked Julia environment under `axisymmetric/` is shared by both models.

## Setup and runners

Julia 1.8.2 is supported and matches the HPC module used for verification:

```bash
scripts/run_debelly.sh instantiate
scripts/run_debelly.sh baseline smoke
scripts/run_debelly.sh mechanochemical smoke
scripts/run_debelly.sh local-inhibition smoke
scripts/run_debelly.sh axisymmetric smoke
```

`smoke` uses a short time horizon, coarse spatial discretization, two chemical
equilibration iterations in 2D, and no figures. It checks installation and
solver plumbing; it is not a scientific reproduction. Each run prints a
machine-readable `DEBELLY_SUMMARY` line with terminal Rac/Rho spatial metrics.

Use `published` for the released parameters (`100` 1D elements, or a `40×40`
2D background mesh; `T=350`; 100 2D equilibration iterations) and generated
plots/VTK data:

```bash
scripts/run_debelly.sh baseline published
scripts/run_debelly.sh mechanochemical published
scripts/run_debelly.sh local-inhibition published
sbatch workflow/run_debelly_axisymmetric.slurm
```

Outputs are deliberately gitignored. Set `DEBELLY_WRITE_OUTPUTS=0` to suppress
them during a full numerical run. Other useful overrides are `DEBELLY_T`,
`DEBELLY_DT`, `DEBELLY_PARTITION`, `DEBELLY_AXIS_N`, `DEBELLY_OPTO_TIME`, and
`DEBELLY_EQUILIBRATION_STEPS`.

## Expected qualitative behavior

The baseline settles to a polarized Rac/Rho state. Localized Rac activation
perturbs that polarity. With mechanochemical feedback enabled, the tension/MCA
coupling transmits the perturbation over the domain; the local-inhibition
control disables those mechanical terms and therefore remains much more local.
The axisymmetric example develops nonuniform Rac/Rho/MCA and cortical-flow
fields along the spherical meridian in response to localized activation.

## 2D runtime

The axisymmetric code prints its solver-only elapsed time and the runner prints
total wall time. Runtime depends strongly on mesh size, Julia compilation,
PETSc configuration, filesystem load, and whether PNG/VTK output is enabled.
Use this reproducible timing command on a compute node:

```bash
DEBELLY_WRITE_OUTPUTS=0 scripts/run_debelly.sh axisymmetric smoke
```

Measured timings and a rough full-run estimate for this cluster are recorded
below after verification. Treat extrapolation from the smoke mesh as an order-
of-magnitude planning number, not a benchmark guarantee.

<!-- RUNTIME_RESULTS -->

## Source layout

- `one_d/Mechanochemical_general_code.jl`: baseline-capable 1D model
- `one_d/protocols/`: published mechanochemical and local-inhibition optogenetics
- `axisymmetric/src/`: Gridap/GridapEmbedded/PETSc axisymmetric implementation
- `axisymmetric/examples/SurfaceViscousFlows.jl`: representative 2D reproduction
- `axisymmetric/{Project,Manifest}.toml`: locked dependency environment
