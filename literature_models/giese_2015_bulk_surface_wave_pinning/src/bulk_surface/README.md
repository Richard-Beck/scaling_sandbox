# Common fixed and moving bulk--surface framework

This package provides one 2D cell cross-section with a phase-field-weighted
Cartesian bulk and a separate ordered, periodic 1D membrane. It is intended as
the common numerical base for configurable polarity reaction networks.

Implemented components:

- implicit conservative diffusion in the 2D bulk;
- implicit tangential diffusion around the closed membrane;
- bilinear bulk-to-membrane traces;
- exactly paired trace/transpose-scatter exchange updates;
- Numba-compatible local reaction kernels and exchange stoichiometry;
- wrapped membrane stimuli;
- Strang-split midpoint reaction and implicit diffusion stepping;
- named conservation pools and common compressed trajectory output.

The membrane exchange rate is evaluated once at each reaction stage. If `T`
is the bulk-to-membrane trace, the paired amount rates are

```text
surface:  M_surface J
bulk:    -T.T M_surface J
```

so a declared bulk--surface pool is conserved independently of grid spacing or
membrane sampling.

Two small configurations are in `examples.py`: reversible exchange and a
membrane-active/bulk-inactive positive-feedback system.

The preferred moving backend retains the same compartment structure: a phase
field on a fixed laboratory grid supports the dynamic 2D bulk, while an
ordered material curve remains the explicit 1D membrane. Membrane nodes move
with and are projected onto the `phi=0.5` contour. Their material amounts are
preserved during stretching, their periodic 1D diffusion operator follows the
curve, and moving bilinear trace/transpose-scatter operators conserve exchange
with the bulk. See `moving_curve.py`.

`diffuse.py` also retains the earlier experimental diffuse-surface backend,
where membrane material occupies the `|grad(phi)|` band. It is kept for
comparison, but is not used by the Mori release workflow.

```python
from bulk_surface import DiffuseGeometry, SimulationConfig
from bulk_surface.diffuse_examples import released_mori_curve_case

geometry = DiffuseGeometry.circle(
    radius=1.5915494309, dx=0.1, interface_width=0.1, box_half_width=4.0
)
case, provenance = released_mori_curve_case(
    "dev/bulk_surface/mori_chiou_shared/v1/runs/mori/volume_1x", geometry
)
result = case.simulator.run(
    case.initial,
    SimulationConfig(end_time=2.0, dt=0.01, output_interval=0.05),
)
```

The same example can be run from the repository root with
`conda run -n scaling_sandbox python scripts/run_bulk_surface_mori_release.py`.

The moving backend uses conservative upwind transport for bulk material with a
smooth extension of the measured interface-normal velocity. The explicit
membrane markers move with that velocity and are projected back to the
`phi=0.5` contour. This is a reduced kinematic closure, not a cytoplasmic
Stokes-flow model; a hydrodynamic velocity extension can replace it without
changing the reaction or exchange interfaces.

```python
from bulk_surface import SimulationConfig
from bulk_surface.examples import wave_pinning_case

case = wave_pinning_case()
result = case.simulator.run(
    case.initial,
    SimulationConfig(end_time=10.0, dt=0.02, output_interval=0.5),
)
result.save("results/bulk_surface/wave_pinning")
```

Run the example driver from the repository root with the existing environment:

```bash
conda run -n scaling_sandbox python scripts/run_bulk_surface_examples.py
```
