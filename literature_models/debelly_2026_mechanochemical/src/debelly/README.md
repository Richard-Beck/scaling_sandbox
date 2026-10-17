# De Belly mechanochemical-polarization reproductions

This directory is a self-contained Julia implementation of the De Belly
mechanochemical models. It does not read from or refer to `dev/`. See
`LICENSE` for attribution and redistribution terms.

It contains the published 1D model, its baseline and front-to-back
optogenetic protocols (mechanochemical and local-inhibition controls), and the
axisymmetric surface model with the supplied spherical-cell example. The 1D
sources retain the authors' Gridap finite-element implementation. The 2D model
is an independent body-fitted finite-element implementation of the authors'
full axisymmetric weak forms. It retains the two-component cortical velocity
and membrane displacement, surface viscosity, active stress, MCA-dependent
friction and transport, nonlinear membrane elasticity, projected tension,
Rac/Rho feedback, leading-edge protrusion, boundary conditions, and staggered
backward-Euler update order. Only the cut-cell embedding and PETSc machinery
are replaced. In the supplied run the reconstructed level-set contour changes
by less than 0.001 micrometres, so a fixed body-fitted meridian represents the
same surface without rebuilding a two-dimensional background mesh every step.

## Setup and runners

Julia 1.8.2 is supported and matches the HPC module used for verification:

```bash
scripts/run_debelly.sh instantiate
scripts/run_debelly.sh baseline smoke
scripts/run_debelly.sh mechanochemical smoke
scripts/run_debelly.sh local-inhibition smoke
scripts/run_debelly.sh axisymmetric smoke
```

`smoke` uses a short time horizon, coarse spatial discretization, and no output
files. The axisymmetric smoke case includes a shortened optogenetic pulse, so
it checks both solver plumbing and spatial response. Each run prints a
machine-readable `DEBELLY_SUMMARY` line with terminal Rac/Rho metrics.

Use `published` for the released 1D parameters or the complete 350-step
axisymmetric protocol at 80 meridional elements:

```bash
scripts/run_debelly.sh baseline published
scripts/run_debelly.sh mechanochemical published
scripts/run_debelly.sh local-inhibition published
scripts/run_debelly.sh axisymmetric published
```

Outputs are deliberately gitignored. The axisymmetric runner writes one tidy
CSV rather than VTK and PNG files. Set `DEBELLY_WRITE_OUTPUTS=0` to suppress
output. Useful overrides are `DEBELLY_T`, `DEBELLY_PARTITION`,
`DEBELLY_AXIS_N`, `DEBELLY_OPTO_TIME`, `DEBELLY_OPTO_START`, and
`DEBELLY_OPTO_STOP`.

The 1D mechanochemical source also accepts `DEBELLY_LENGTH`,
`DEBELLY_MCA_DENSITY`, `DEBELLY_MCA_NODAL_TOTAL`, and
`DEBELLY_OPTO_WIDTH`. Setting `DEBELLY_TRAJECTORY_OUTPUT` writes tidy chemistry,
velocity/displacement, and tension CSV files without enabling the original VTK
and plotting path. These controls are used by the WGD adapter; their defaults
leave the published reproduction unchanged.

## Expected qualitative behavior

The released 1D baseline starts exactly flat and remains symmetric. Localized
Rac activation perturbs that state. With mechanochemical feedback enabled, the
tension/MCA coupling transmits the perturbation over the domain; the local-
inhibition control disables those terms and remains more local. The
axisymmetric example develops front-biased Rac, reciprocal Rho, redistributed
MCA, cortical flow, tension, and membrane displacement during localized
activation, followed by post-stimulus evolution.

## 2D runtime

The axisymmetric code prints solver time and the runner prints total wall time.
Use this reproducible timing command:

```bash
DEBELLY_WRITE_OUTPUTS=0 scripts/run_debelly.sh axisymmetric smoke
```

On the verification compute node with Julia 1.8.2, the complete `T=350`,
80-element simulation takes about four seconds after package precompilation.
The 2D project uses only Julia standard libraries and precompiles in about one
second.

The implementation was compared directly with the authors' saved VTUs at
times 75, 100, 150, 200, 300, and 350. Across those times, full-profile
span-normalized RMSE is below 1.6% for Rac, Rho, and bound MCA; below 1.3% for
tension; below 4.1% for velocity magnitude; and below 1.6% for displacement
magnitude. At time 350, published versus reproduced extrema are:

- Rac: `0.2156--1.4471` versus `0.2156--1.4541`;
- Rho: `0.6335--1.9101` versus `0.6288--1.9101`;
- bound MCA: `0.02210--0.05145` versus `0.02154--0.05156`;
- velocity magnitude: `0--0.02050` versus `0--0.02127`;
- displacement magnitude: `0--0.08720` versus `0--0.08233`.

For a checkout that also has the upstream VTUs, the optional regression command
probes both trajectories on the same meridian and fails above 5% NRMSE:

```bash
conda run -n vtk-render python scripts/compare_debelly_axisymmetric.py \
  /path/to/upstream/output /path/to/axisymmetric.csv
```

## Source layout

- `one_d/Mechanochemical_general_code.jl`: baseline-capable 1D model
- `one_d/protocols/`: published mechanochemical and local-inhibition optogenetics
- `axisymmetric/src/DeBellyAxisymmetric.jl`: full body-fitted axisymmetric model
- `axisymmetric/MODEL_EQUIVALENCE.md`: equation mapping and quantitative validation
- `axisymmetric/examples/SurfaceViscousFlows.jl`: representative 2D reproduction
- `axisymmetric/Project.toml`: shared Julia dependency environment
