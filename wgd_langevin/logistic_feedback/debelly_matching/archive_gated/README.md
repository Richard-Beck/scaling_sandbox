# De Belly half-mass polarity matching

Completed pilot: [results and interpretation](results/REPORT.md),
[calibration overlays](results/calibration.png), and
[held-out countercues](results/heldout_countercues.png). The establishment fit
is close, but the fitted parameters fail the held-out reversal test.

This deterministic pilot reduces the retained
`Figure2_standalone_v2_julia_activation_defaults_linear50.png` data to the
two-variable logistic feedback model. It does not rerun or modify De Belly.

The figure's interpolated half-mass coordinates are used without renormalization:

    p = Rac_x50/L - 0.5
    v = 0.5 - bound_MCA_x50/L

Both states are relative polarity coordinates, not calibrated velocity. Times
are seconds. A positive state indicates rightward Rac polarity / rightward MCA
depletion. Plotting reverses the transformation to retain the original x50/L
coordinates. Relative displacement is not included in the calibration.

    dp/dt = -p/tau_p + beta*S(p*v)*v/(1+abs(v)/vs) + cue_gain*u(t)
    dv/dt = -v/tau_v + alpha*S(p*v)*p/(1+abs(p)/ps)
    S(z) = 1/(1+exp[-k*(z-z0)])

The signed input u(t) is the analytical right-half minus left-half mean of the
manifest's continuous cosine cue. A single left cue of amplitude 3.6 and width
fraction 0.2 has u=-0.72; the two-sided 15% asymmetric cue has u=-0.108. The
original time envelope, including half-cosine taper, is retained. The fitted
gain is a projection assumption; symmetric activation is omitted by this
two-state reduction. No target trajectory is used as an external forcing.

## Calibration and validation

- Reference length: 29 um, chosen as the middle of the sustained polarized
  lengths in panel C. Selection uses polarity status, not fitting performance.
- Calibration: panel A/full (160-s tapered left cue) and panel C/beta_1x
  (80-s tapered, asymmetric two-sided cue), both at 29 um. Initial p=v=0.
- Loss: both fields, scaled by each reference trajectory's maximum absolute
  centered amplitude, with a floor of 0.02. Each assay receives equal weight;
  early data (through 500 s after withdrawal) and later data receive equal
  weight. Samples are approximately 5 s early and 50 s late. This balances
  transients against long plateaus. The loss is descriptive, not a likelihood.
- Eight reproducible bounded starts optimize nine parameters. Positive
  parameters are optimized in logarithms; k*z0 is optimized directly and can
  be negative. Bounds and the actual dimensional z0 are saved. Both noise
  strengths are fixed at zero. Optimizer convergence does not establish
  parameter identifiability or mechanistic equivalence.
- Held out: all panel E countercues (150/300/600 s), and all 3,000 duration-sweep
  outcomes. No countercue data enter parameter selection.
- Transfer: frozen reference parameters are tested across all 20 lengths and
  the three feedback settings. Since the reduced model has no native length
  or De Belly feedback parameter, these are explicit tests of unchanged
  parameters, not an assumed mapping of De Belly beta to reduced beta.
- Conditional calibrations at 17/21/25/33/39 um (1x), and at 29 um (0x/2x),
  each start from the reference fit. They use A/C only; the 2x setting has only
  C available. These are exploratory single-start mappings, not globally
  optimized estimates. Their reversal predictions remain held out.

The RK4 integrator aligns with all pulse events and sampling times and uses
a maximum 0.5-s step. Reference calibration and held-out trajectories are
checked at half that step. Coordinates are not clipped to their physical
half-mass bounds; excursions are explicitly recorded as model inadequacy.

## Reversal criteria

The final reference is the median over the final 1,000 s. A state must enter
and remain within 0.002 of that reference after cue withdrawal, with at least
1,000 s confirmation. Reversal requires the opposite sign, a final value
within 0.05 of the reflected conditioning reference, and both states at least
0.02 from the midpoint. These are the figure's criteria.

For a strict external comparison, primary reduced-model reversal labels use
the same De Belly conditioning reference as the retained audit. A separate
self-conditioned label uses the reduced model's own last pre-countercue
state, so poor conditioning-amplitude transfer is distinguishable from failure
to change sign. Missing thresholds mean no successful tested duration, not an
estimated infinite threshold; unpolarized cases are not labeled reversals.

## Run

The existing `scaling_sandbox` environment provides all dependencies. The
full multistart and transfer campaign uses a single CPU and 3 GB on SLURM,
with `small` QoS selected for the small resource footprint.

```sh
mkdir -p wgd_langevin/logistic_feedback/debelly_matching/results
sbatch wgd_langevin/logistic_feedback/debelly_matching/run.slurm
```

Run from the repository root. `--resume` reuses completed parameter fits;
it is suitable only when inputs and fitting code have not changed. JSON
checkpoints are saved after every optimization start. No arrays or concurrency
throttles are used. The script prints progress to the SLURM log.

```sh
OPENBLAS_NUM_THREADS=1 conda run -n scaling_sandbox python -m unittest discover -s wgd_langevin/logistic_feedback/debelly_matching -p 'test_*.py' -v
```

Outputs include the best parameter set, all multistart candidates, conditional
parameter mappings, full trajectory error metrics, all predicted duration-sweep
labels, threshold comparisons, PNG overlays, timestep checks, and input/source
hashes. The retained CSV data and manifests are the sources of truth; no new
half-mass estimator is substituted.

After the job completes, generate the report, representative trajectory archive,
gate diagnostic and extended timestep checks with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/debelly_matching/summarize.py
```
