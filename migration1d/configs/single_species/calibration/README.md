# Single-species calibration

Calibration uses the common membrane--cytoplasm strip model and proceeds in
three filters rather than a broad seven-dimensional sweep.

1. `local_bistability.json` scans roots of `J(u,v)=0` with `v` fixed and labels
   each root by local stability. A three-root interval is a necessary filter
   for classical wave-pinning.
2. `regime_map.json` applies one transient Gaussian activation pulse near the
   left end and maps the late-time outcome over `(mu, Lambda)`.
3. The regime map is repeated at `h=1,2,4` to expose the effect of cytoplasmic
   volume per membrane area.

Each plotted `h` slice holds `mu=mean(u+h*v)/K` fixed. A biological scaling
hypothesis that instead preserves membrane and cytoplasmic concentrations
should be evaluated as a path through the map, with `mu(h)=u/K+h*v/K`, rather
than by comparing fixed-`mu` pixels across slices.

The baseline holds `rho=100`, `kappa=0.067`, `Gamma=1`, and `n=2`. Production
maps run for 300 seconds. Every PDE row reports polarization time, terminal
state, amplitude, pole count, high-domain width, endpoint contrast, late-time
change, and conservation error.

Classification uses the late trajectory. The normalized amplitude is

```text
P = (max(u) - min(u)) / mean(u)
```

and high domains are connected regions above the midpoint between the current
minimum and maximum. Nonstationary late trajectories are separated into
coarsening and oscillatory/dynamic outcomes. Thresholds are recorded in
`../../_standards/calibration_v2.json`.

The production regime map is compute-heavy and is configured for SLURM. Small
grids can be run locally by overriding `--mu-count`, `--lambda-count`, and
`--h-values`.

Both calibration campaigns write figures and section metadata automatically.
The report intentionally emphasizes the categorical outcome map and direct
outcome counts by geometry. Near-uniform diagnostic maps (including the
stimulus-onset time and solver conservation error) remain available in the CSV
but are not promoted as scientific figures.
