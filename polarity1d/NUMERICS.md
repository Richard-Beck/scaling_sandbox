# Numerical strategy and measured selection

## Decision

The production defaults are:

1. sparse adaptive BDF on a conservative finite-volume spatial discretization
   for wave-pinning, Otsuji, Goryachev, and LEGI;
2. a specialized conservative, banded staggered solver for De Belly;
3. job-level parallelism under SLURM for future screens, not threads inside a
   single small 1D solve.

The common framework is the model/protocol/result/phenotype contract, not a
requirement that every equation use the same time integrator.

## Strategies tested

### Sparse method of lines + BDF (selected for Jilkine)

The no-flux PDEs use a cell-centred finite-volume Laplacian. The semi-discrete
system is passed to adaptive BDF with a supplied block-local Jacobian sparsity
pattern. This handles the 100-fold active/inactive diffusion contrast and the
fast Otsuji reaction without a parabolic CFL restriction.

Default versus tight-tolerance BDF terminal NRMSE was `6.4e-6` or lower for all
four 300-second benchmark cases. Median time was 0.053–0.120 seconds.

### Exact diffusion + Strang/RK4 (retained oracle)

The finite-volume Neumann operator is diagonalized by a DCT-II. Diffusion is
therefore exact for the semi-discrete linear operator and unconditionally
stable; local reactions use RK4 with a model-specific step and positivity
retry. It agreed with tight BDF to 0.02–1.14% terminal NRMSE, but was 10–27×
slower on the benchmark. It remains valuable because its error mechanism and
code path are independent of BDF.

### Explicit adaptive RK45 (rejected)

At only 41 cells and 20 simulated seconds, RK45 needed 28,000–29,000 RHS calls
and 1.27–1.30 seconds. BDF needed 171–247 calls and 0.024–0.031 seconds. The
diffusive stiffness makes an explicit production backend unattractive.

### Released Julia/Gridap finite elements (retained reference)

The authors' staggered finite-element implementation is the strongest
provenance source and remains runnable. A fresh-process six-second, 20-element
smoke case took 59.7 seconds on the benchmark node, dominated by Julia/Gridap
compilation and assembly startup. That cost is inappropriate for one process
per small screening task.

### Conservative banded De Belly port (selected)

The optimized port retains the released staggered update order and physical
variables while replacing general FE assembly with 1D banded solves:

- membrane and cortical force balances: tridiagonal elliptic systems;
- bound/unbound MCA transport: backward-Euler diffusion plus conservative
  first-order upwinding (an M-matrix, hence positivity-preserving);
- MCA on/off exchange: exact local relaxation, preserving positivity and total;
- Rac/Rho: positive backward-Euler diffusion/linear decay with staggered mutual
  inhibition;
- total MCA: conserved to `4.4e-16` in the 350-second benchmark.

At 100 cells, `dt=0.5 s` took 0.25 seconds and differed by 0.97% terminal NRMSE
from `dt=0.25 s`; `dt=0.25 s` differed by 0.43% from `dt=0.125 s`. Spatial
NRMSE fell from 1.60% (50→100 cells) to 0.83% (100→200 cells). The released
0.98-per-second boundary relaxation is exponentiated by `dt`, so refinement
preserves its physical timescale. Use `dt=0.25 s` or smaller to confirm results
near classification boundaries.

The exact measurements, versions, hostname, and repetitions are in
[`results/benchmarks/backend_benchmark.json`](results/benchmarks/backend_benchmark.json).

## Stability gates

The implementation fails rather than silently accepting:

- non-finite states;
- materially negative concentrations;
- incompatible state/grid shapes;
- failed adaptive integration;
- invalid physical parameters.

Regression tests additionally check conserved material, baseline symmetry,
force-response ablations, mesh refinement, the independent solver comparison,
and common phenotype scoring.

## Size simulations

Kinetics and physical diffusion coefficients remain fixed when `length` is
changed. This is equivalent to the `D/L²` formulation on a unit coordinate.
Two meshes have distinct roles:

- fixed cell count: closest numerical reproduction of the published unit-grid
  experiment;
- fixed physical `dx`: required refinement control so enlargement is not
  confounded by spatial coarsening.

Long-time results must be paired with explicit resolution-time metrics and
diffusive-time rescaling; a state unresolved at an arbitrary endpoint is not
automatically a new stable state.
