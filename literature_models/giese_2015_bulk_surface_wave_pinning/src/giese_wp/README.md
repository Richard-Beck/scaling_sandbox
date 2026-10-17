# Giese bulk--surface wave pinning

This module implements steps 1--4 of `docs/jobs/01_giese.txt`: the native
wave-pinning branch from Giese et al. (2015), without a WGD adapter.

The cytoplasmic field `v` uses conforming P1 elements on a triangular disk.
The membrane field `u` uses P1 elements on the disk's boundary edges.  A
two-point edge quadrature assembles the nonlinear reaction load once; that
same vector enters the membrane and bulk weak equations with opposite signs.
Consequently the discrete total

```text
1' M_membrane u + 1' M_bulk v
```

is conserved to sparse-solver roundoff.

The published parameters and initial concentrations are the defaults in
`GieseParameters`.  The native design is

```text
diameter: 1.5, 3, 4.5, 6, 9, 12, 15 um
single site: 0.44 um/s at angle 0
two sites: 0.44 um/s at angle 0 and 0.40 um/s at angle pi
each site: 5% of the circumference for 10 s
```

`POL` follows the published peak/mean/perimeter formula.  Published `PF` is
one minus the relative length of the above-mean cluster around the global
maximum.  Numerically, that cluster is the connected above-mean component
containing the maximum; this is the paper's small-neighbourhood condition
without an arbitrary radius.  Profiles indistinguishable from homogeneous
at a relative range of `1e-5` receive `PF = 0`.

`half_mass_pf` implements the conflicting preprint/Results-section definition:
one minus twice the relative length of the smallest connected patch around the
global maximum containing half of the membrane mass.  Both PF definitions are
retained explicitly rather than silently substituting one for the other.

The default solver remains the exactly conservative reaction-explicit IMEX
scheme.  `scheme="giese_semi_implicit"` selects the paper's asymmetric split
(implicit membrane `u` with lagged `v`, then implicit bulk `v` with lagged
`u`) for comparison.  `make_gmsh_circular_system` generates an irregular disk
with explicit refinement toward its outer membrane.

Use the repository environment:

```bash
conda env create -f environment.yml
conda run -n scaling_sandbox pytest -q
conda run -n scaling_sandbox python scripts/run_giese_validation.py
```

The 14-condition native sweep is a compute workflow.  Submit it with
`bash workflow/submit_giese_native.sh`, then consolidate its checkpoints with:

```bash
conda run -n scaling_sandbox python scripts/run_giese_native.py --collect
conda run -n scaling_sandbox python scripts/plot_giese_native.py
```

Raw checkpoints are written under `reports/output/giese_wp/native/` and are
intentionally untracked.
