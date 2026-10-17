# Marée et al. (2006) reconstruction

This folder contains a standalone, close reimplementation of the model described
in [`spec.md`](spec.md). It uses a custom C++ simulator because the local frozen
Morpheus 2.4.1 executable provides the standard CPM Hamiltonian but does not expose
XML hooks for the paper's orientation-selected pushing energy or the conservative
per-copy `F/B/P` transfers. The implementation does not require a new environment
or third-party dependency; field updates use OpenMP while the CPM sweep is serial.

## What is implemented

- 500×500 periodic axial hexagonal lattice; six nearest directions and all 36
  graph-distance ≤3 neighbors for the CPM interfacial energy.
- Seven biochemical fields, six `F`, six `B`, and six membrane `P` fields with the
  published parameters.
- No-flux hexagonal diffusion and first-order upwind actin transport. The default
  multirate solver uses explicit reaction/advection steps of 0.002 s and implicit
  backward-Euler diffusion at 0.02 s for active GTPases and 0.01 s for Arp2/3 and
  inactive pools. The original 0.0002-s explicit solver remains available as a
  reference.
- Random-order, one-attempt-per-site CPM sweeps, the published yield-offset
  Metropolis probability, area/contact energies, direction-selected protrusion,
  and positive-part Rho contraction.
- Immediate `B↔P` and oriented retraction transfers plus post-MCS conservative
  moving-domain normalization.
- The 2012-derived 15% `I_C` pulse for the first 10 s.

## Build, test, and run

```bash
cd literature_models/Maree_2006
make test

# A short local run (10 seconds of model time)
OMP_NUM_THREADS=4 ./maree2006 --steps 100 --output runs/short

# Slow numerical reference
OMP_NUM_THREADS=4 ./maree2006 --integrator explicit-reference \
  --steps 100 --output runs/explicit_short

# Full 20-minute run on the cluster
sbatch run_full.slurm runs/polarized_seed1

# The same job script can run a bounded 10-second pilot
sbatch --time=00:20:00 run_full.slurm runs/pilot_seed1 100
```

`metrics.csv` contains every MCS. `snapshots.csv` stores all intracellular sites
at the requested interval (10 s by default), and `manifest.txt` records the run
settings and largest pre-correction moving-domain errors. Use `--no-pulse` for the
resting control, `--frozen` for fixed-domain chemistry, and `--no-snapshots` for a
small diagnostic run. Run `./maree2006 --help` for the entry point; the source lists
all supported flags in `parse()`.

Summarize a completed run with:

```bash
python3 summarize_run.py runs/polarized_seed1
```

See [`BOUNDARY_COPY_DIAGNOSTIC.md`](BOUNDARY_COPY_DIAGNOSTIC.md) for the targeted
copy-energy tests and the diagnosed area-penalty scaling problem.

## Reconstruction choices and limitations

These choices correspond to the unresolved items in `spec.md`, rather than claims
about the historical source code:

1. Hex sites use axial coordinates. The Laplacian is
   `2/(3 dx²) sum(neighbor-center)`. A site has area `dx²` for the published CPM
   area calculation, even though the geometric Voronoi area of an axial lattice is
   `sqrt(3) dx²/2`. The initial cell is the nearest 30,000 sites to the center, so
   its effective diameter is about 18.2 µm rather than exactly 20 µm.
2. Barbed ends use conservative donor-cell upwinding. A `B` flux crossing the
   membrane becomes `P`; only orientations with positive dot product with an
   attempted displacement contribute to protrusion, without cosine weighting.
3. The dimensional conversion is `P=B*dx` and `B=P/dx`. Retraction transfers each
   orientation toward the opposite lattice direction. If that preferred receiver
   is unavailable, the first adjacent intracellular site is used.
4. Rho contractility and membrane production use positive-part geometry. A local
   six-neighbor connectivity check rejects potentially fragmenting retractions.
5. All seven scalar fields are normalized separately after geometry changes.
   `F_m` is conserved separately, while `B_m*dx² + P_m*dx` is conserved for each
   orientation. The pre-correction discrepancy is reported in `manifest.txt`.
6. The third-order contact stencil is all 36 sites of hex graph distance 1–3;
   copy proposals use the six nearest neighbors.
7. The accelerated default uses a conservative conjugate-gradient solve of
   `(I-dt*D*L)q_new=q_old` with tolerance `1e-9`. Diffusion is normalized only for
   roundoff-level mass error after each solve. `--integrator explicit-reference`
   retains the specification's inferred explicit scheme for comparisons.

The simulator has not been claimed as a pixel-for-pixel recovery. In particular,
the exact historical stencil, advection rule, connectivity policy, and boundary
conversion remain unavailable. A biologically long run is supplied as a SLURM job;
the explicit-reference mode remains substantially more expensive.

## Validation currently automated

`make test` checks active+inactive GTPase conservation with both integrators, the
published actin steady-state arithmetic, a two-MCS moving-domain smoke run, and
the paper's analytic thermal-ratchet calibration (maximum relative curve mismatch
below 9% on `b=0.05..0.25 nm^-1`). The longer resting, polarization, migration,
actin-architecture, contractility-sweep, and turning validations require scheduled
runs and remain empirical acceptance tests rather than unit tests.

## Native-size timing and agreement

The explicit-reference 100-MCS (10-s) pilot took 2:06 on four CPUs. The implicit
multirate version took 18.6 s, a 6.8-fold speedup. A coupled 200-MCS smoke test took
39.2 s, projecting to roughly 39 minutes for 20 model minutes. Four CPUs remain the
best tested allocation: a 40-MCS polarized benchmark took 26.1, 13.7, 7.3, and
15.5 s on 1, 2, 4, and 8 threads respectively.

At 10 s, every accelerated field had correlation above 0.999998 with the explicit
reference. Relative L2 differences were at most 0.055% (`C`); all other scalar and
total actin fields were below 0.046%. Use `compare_runs.py` to repeat this check.
Neither solver produced a boundary copy in the first 20 model seconds. The stopped
explicit long run also remained immobile through 252.1 s, so migration mechanics
still require diagnosis independently of the solver acceleration.
