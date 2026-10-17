# Axisymmetric model equivalence

The implementation in `src/DeBellyAxisymmetric.jl` is a specialized,
body-fitted discretization of the model released in
`MechanochemicalAxisymmetricVector.jl` and `WeakForms.jl`. It is not a reduced
surrogate. This note records what is preserved and what is only a numerical
implementation detail.

## Preserved model

- Geometry: the published upper spherical meridian of radius 6, from
  `(-R,0)` to the truncated leading edge `x=0.9R`, with the axisymmetric
  integration weight `y`.
- Cortical flow: both Cartesian velocity components; covariant strain,
  meridional divergence, azimuthal strain/divergence terms, viscosity,
  MCA-dependent friction, Rho-dependent active stress, weak tangentiality,
  and the two endpoint Dirichlet conditions.
- Membrane mechanics: both displacement components; the complete `aM` and
  `aL` forms including their old-displacement nonlinear terms; MCA friction;
  velocity forcing; normal penalty; projected tensor tension; and the relaxed,
  Rac-triggered leading-edge boundary displacement.
- Chemistry: the released Rac and Rho reaction-diffusion forms, smooth MCA and
  tension thresholds, 100-step initial equilibration, localized optogenetic
  activation at steps 50--150, and the same physical constants.
- MCA: bound and unbound species, binding/unbinding, both diffusion constants,
  bound cortical advection and dilation, and the released mass-residual source.
- Coupling: the original staggered backward-Euler sequence—velocity,
  displacement, tension, boundary update, Rac, Rho, bound MCA, unbound MCA.

The authors' displacement bookkeeping assigns `x_old = x` before the MCA and
next velocity solves. Consequently the nominal displacement-rate terms in
those two forms are identically zero in the released run; the rewrite preserves
that executed behavior.

## Replaced numerical machinery

The released solver represents the one-dimensional meridian as a cut through
a two-dimensional Cartesian Gridap mesh, extends fields into a narrow band,
and solves with PETSc/MUMPS. The rewrite assembles the same surface weak forms
directly on linear meridian elements and solves the sparse systems with Julia's
SuiteSparse interface. Cut-cell extension stabilization is therefore neither
needed nor included; it is a conditioning device, not a physical term.

Although the released code recalculates a level set, direct contours of its
saved `_sur_1`, `_sur_75`, `_sur_150`, and `_sur_350` fields are identical to
within the output resolution. Their radial RMS deviation from the radius-six
sphere is `0.00090` micrometres at every checkpoint. The rewrite consequently
uses that fixed meridian. This is an observed property of the supplied model
run, not an assumption that the displacement field is zero: membrane
displacement remains a fully solved two-component field.

## Numerical validation

`scripts/compare_debelly_axisymmetric.py` probes the authors' VTUs and the new
CSV on the same 240 meridian points. At times 75, 100, 150, 200, 300, and 350,
the worst span-normalized profile RMSE from the default 80-element run is:

| Quantity | Worst NRMSE |
|---|---:|
| Rac | 0.45% |
| Rho | 0.23% |
| bound MCA | 1.56% |
| tension | 1.26% |
| velocity magnitude | 4.02% |
| displacement magnitude | 1.51% |

The larger velocity percentage occurs at time 75 while velocities are still
small. At time 350 the velocity NRMSE is 3.10%. Runs with 80, 160, and 320
elements give time-350 Rac minima of `0.215633`, `0.215633`, and `0.215633`,
and maxima of `1.4541`, `1.4527`, and `1.4524`, respectively.
