> **Flawed deBelly extension deprecated (2026-10-08):** simulations using the old two-end extension with `center_boundary_motion=True` must not be used for scientific inference or current figures. It forces equal-and-opposite commands even under one-end Rac activation and changes strain/tension. The corrected reflection-symmetric extension uses independent local endpoint commands (`center_boundary_motion=False`); see [manuscript](manuscript/README.md). The published one-sided model is a separate boundary assay and is not implicated by this centering defect.

> **Manuscript workflow moved (2026-10-03):** use [manuscript/](manuscript/README.md) for the clean five-model corrected-preparation reversal comparison, its main figure, supplements and intermediate data. The seven-model workflow described below is deprecated for manuscript work; see [the retention inventory](manuscript/RETENTION.md).

# polarity1d — WP1 model-development framework

## Accepted seven-model four-assay baseline

The current four-assay length comparison is in [four_assay/](four_assay/README.md),
with seven model files, one runner and one plotter. Its fresh full rerun matches
the accepted revised-v8 baseline. [Combined figure](four_assay/results/length_response_and_settling.png) ·
[Slope heatmap](four_assay/results/length_response_gradients.png) ·
[Numeric verification](four_assay/results/comparison_to_baseline.json).
Historical report campaigns are preserved in `archive/four_test_suite_legacy/`.

> **Scientific scope decision:** the current Jilkine qualification is frozen as
> an equation-level, qualitative reconstruction. Exact figure matching is not
> pursued because the publication does not specify a reproducible
> symmetry-breaking mechanism and omits other numerically consequential
> details. Results from symmetric cues must not interpret solver-seeded winner
> selection as biology. The optimized De Belly port is accepted as the WP1
> screening model without an additional native-Gridap trajectory cross-check.

`polarity1d` is the replacement 1D framework for Work Package 1. It contains
the four canonical models compared by Jilkine and Edelstein-Keshet (2011) and
an optimized 1D port of the retained De Belly mechanochemical implementation.
The older `migration1d/` package is deprecated and is not imported here.

The model-building pass finished common protocols/phenotypes, numerical
qualification, and performance selection. Production campaigns build on that
frozen baseline rather than tuning the implementations to published figures.

## Included models

| Config name | State / observable | Literature behavior represented |
|---|---|---|
| `wave_pinning` | active + inactive GTPase / active | bistable mass-conserved wave pinning |
| `otsuji` | active + inactive GTPase / active | mass-conserved redistribution instability |
| `goryachev` | active + inactive Cdc42 / active | reduced positive-feedback Cdc42 polarity |
| `legi` | activator + inhibitor + response / response | local excitation and global inhibition |
| `debelly` | Rac, Rho, MCA, displacement, velocity, tension / Rac | reciprocal mechanochemical Rac–Rho communication |

All models return the same `Trajectory` object and accept the same physical
domain, output-time, and scalar-stimulus abstractions. Model-specific solvers
remain behind that interface: forcing unrelated equations through one
integrator was slower and less stable.

## Quick use

Use the existing repository environment:

```bash
conda run -n scaling_sandbox pytest -q polarity1d/tests

conda run -n scaling_sandbox python polarity1d/run.py \
  polarity1d/configs/wave_pinning_dual.json \
  /tmp/wp1_wave_pinning
```

Without an editable install, add the source directory for an interactive run:

```bash
PYTHONPATH=polarity1d/src conda run -n scaling_sandbox python -c '
from polarity1d import LocalizedPulse, WavePinning, phenotype_summary, simulate
run = simulate(
    WavePinning(),
    protocol=LocalizedPulse(sides=("left", "right"), asymmetry=0.02),
    length=10, cells=101, end_time=300, backend="bdf"
)
print(phenotype_summary(run, analysis_start=40))
'
```

The common protocol repertoire includes single/dual localized pulses, Gaussian
opto pulses, persistent gradients, gradient reversal, and alternating gradients
for response-bandwidth work. `noisy_initial_state()` supplies reproducible,
zero-mean perturbations and preserves declared molecular pools.
`phenotype_summary()` reports both physical and `time/L²` establishment and
resolution clocks so kinetic delay can be separated from a changed endpoint.

## Production campaign

The manifest-driven Phase 1 campaign adds explicit endpoint classes,
event-establishment checks, right-censoring, reversal conditioning, normalized
tracking metrics, and regional mechanochemical response metrics. Tasks are
digest-locked, resumable, independently executable by SLURM array index, and
collected with missing/failure validation. See [CAMPAIGN.md](CAMPAIGN.md) for
the complete assay definitions, parameter grids, replication, phenotype rules,
and reproduction commands. The [completed v1 report](results/campaigns/phase1_production_v1/README.md)
separates all seven assay families and pairs summary scans with representative
trajectory figures.

## Selected numerics

- Jilkine models: cell-centred conservative finite volumes plus SciPy's sparse
  adaptive BDF. A Strang/RK4 solver with exact discrete-Neumann diffusion is
  retained as an independent regression oracle.
- De Belly: cell-centred finite volumes, conservative implicit upwinding,
  exact positive MCA binding exchange, staggered backward Euler, and banded
  linear solves. The default `dt=0.5 s` follows the measured speed/convergence
  tradeoff; use `dt=0.25 s` as the routine refinement check near a phenotype
  boundary.
- Use approximately fixed physical `dx` across size comparisons (the compact
  qualification uses `dx=0.1 µm`). Holding `cells` fixed while increasing `L`
  intentionally reproduces Jilkine's unit-computational-domain setup, but it
  should not be the only mesh used for a scientific size result.

On the recorded HPC node, a 300-second, 101-cell Jilkine trajectory takes
0.05–0.12 seconds with BDF. A 350-second, 100-cell De Belly trajectory takes
about 0.25 seconds at `dt=0.5 s` and conserves MCA to machine precision. See
[the measured backend table](results/benchmarks/backend_benchmark.md) and
[NUMERICS.md](NUMERICS.md) for the complete selection rationale.

## Qualification assets

- [Canonical Jilkine panel](results/qualification/jilkine_qualification.png):
  raw profiles for all four models at selected times on a 10-µm domain, a
  common-stimulus `T=2000 s` length row, and a separate `T=2000 s`
  figure-guided length row. Both length rows use `L=1, 3, 17, 25, 39 µm`.
  See `MODEL_ASSUMPTIONS.md` and the machine-readable summary for the distinct
  protocol conventions.
- [De Belly panel](results/qualification/debelly_qualification.png): Rac/Rho
  propagation and independent reciprocal-arm ablations.
- [Machine-readable summary](results/qualification/qualification_summary.json).

These are compact deterministic qualification runs, not a production screen.
Regenerate them with:

```bash
conda run -n scaling_sandbox python polarity1d/scripts/qualify_models.py
conda run -n scaling_sandbox python polarity1d/scripts/benchmark_backends.py
```

If/when a broad screen is authorized, submit it through SLURM. The per-run
timings here are suitable for sizing jobs; no array concurrency throttle should
be added unless scheduler or filesystem constraints justify one.

## Scientific caveats

The optional [fixed-substrate clutch extension](FIXED_SUBSTRATE_CLUTCH.md)
adds adhesion turnover, elastic bond extension, and freely translating
L=21 µm cells to the symmetric de Belly assay. Its illustrative parameters,
numerical checks, and moving kymographs are documented separately from the
literature-model reproductions.

The Jilkine review does not specify every reproduction constant or a
reproducible symmetry-breaking mechanism. The framework therefore freezes and
discloses the unresolved Goryachev and wave-pinning assumptions rather than
tuning them. The 2026 de Belly port is equation- and parameter-derived. Its
native Julia/Gridap trajectory comparisons are explicit audit gates rather
than a waived equivalence claim. The authors' retained Julia source remains
under `literature_models/debelly_2026_mechanochemical/` as the provenance reference.
See [MODEL_ASSUMPTIONS.md](MODEL_ASSUMPTIONS.md).
