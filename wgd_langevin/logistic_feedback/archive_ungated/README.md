# Ungated polarity–velocity model

The logistic gate is removed from both feedback terms. k and z0 are no longer parameters.

    dp = [-p/tau_p + beta*v/(1+|v|/vs)] dt + sigma_p dW_p
    dv = [-v/tau_v + alpha*p/(1+|p|/ps)] dt + sigma_v dW_v
    dx = v dt

Norms are Euclidean. Independent Wiener noise is applied per coordinate.
Previous gated source snapshots are in `archive_gated/`; previous results remain
in `results/` and `results_dimensional/`.

## Requested 2D examples

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/simulate.py --dimensional --tau-p 1 --tau-v 10 --alpha 1 --beta 1 --vs 1 --ps 1 --duration 1000 --sample-interval 1 --noise 0.15 0.45 --out wgd_langevin/logistic_feedback/results_dimensional_ungated
```

[Overview PNG](results_dimensional_ungated/trajectories_overview.png).
Six trajectories per noise level, with the same seeds and initial conditions as
before: three resting, three aligned with magnitudes 2 at 0/120/240 degrees.
All positions start at zero. No burn-in or selection. Euler–Maruyama uses dt=0.005
and old-state drift and velocity; output is exactly t=0,1,...,1000. Noise
amplitudes are state units per square-root time. No physical unit names were supplied.

NPZ keys t,p,v,x store time and arrays shaped time × track × coordinate.
Individual PNG panels have equal X/Y scaling but independent limits. Diagnostics
show speed and polarity magnitude. Metadata records parameters, initial states,
seeds, versions, generation time and source hash. Python supports unequal noise
strengths; CLI uses equal p/v noise. Without --dimensional, gp/gv/r specify the
ungated nondimensional equations.

## Calibration and tests

[Ungated matching report](debelly_matching/results_ungated/REPORT.md).
The 2D examples use the user's parameters, not the De Belly fitted values.

```sh
OPENBLAS_NUM_THREADS=1 conda run -n scaling_sandbox python -m unittest discover -s wgd_langevin/logistic_feedback -p test_simulate.py -v
```
