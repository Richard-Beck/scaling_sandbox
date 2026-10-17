# Logistic polarity–velocity feedback

## Dimensional simulation

The `--dimensional` option integrates the original equations directly:

    dp = [-p/tau_p + beta S(p dot v) v/(1+|v|/vs)] dt + sigma_p dW_p
    dv = [-v/tau_v + alpha S(p dot v) p/(1+|p|/ps)] dt + sigma_v dW_v
    dx = v dt
    S(z) = 1/(1+exp[-k(z-z0)])

Reproduce the dimensional example using the existing Conda environment:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/simulate.py --dimensional --tau-p 1 --tau-v 10 --alpha 1 --beta 1 --vs 1 --ps 1 --gate-z0 0.5 --gate-k 1 --duration 1000 --sample-interval 1 --noise 0.15 0.45 --out wgd_langevin/logistic_feedback/results_dimensional
```

Here `--gate-k` and `--gate-z0` specify dimensional gate parameters. The noise
values are dimensional Wiener amplitudes (state units per square-root time),
independently applied to each coordinate of p and v. The internal step is 0.005;
plots and saved NPZ data use exactly 1,001 samples at t = 0, 1, ..., 1000.
NPZ keys are lowercase `t`, `p`, `v`, `x`. Six tracks per noise level use the
same resting/aligned initialization described below, now in dimensional units.
No physical unit names were supplied. The earlier nondimensional results are
preserved in `results/`; new figures and metadata are in `results_dimensional/`.

## Original nondimensional example

Implements the supplied model in 2D, with independent additive Wiener noise
in every polarity and velocity component. This is an illustrative simulation,
not a fitted biological model.

    dP = [-P + gp S(P dot V) V/(1+|V|)] dT + sigma_P dW_P
    dV = r [-V + gv S(P dot V) P/(1+|P|)] dT + sigma_V dW_V
    dX = V dT
    S(z) = 1/(1+exp[-K(z-Z0)])

The noise terms are SDE increments, not ordinary time derivatives. The
nondimensional noise amplitudes are sigma_P = sigma_p sqrt(tau_p)/p_s and
sigma_V = sigma_v sqrt(tau_p)/v_s. In particular, r multiplies the velocity
drift only. Position X = x/(v_s tau_p); no dimensional distances or times are
assumed.

Defaults: gp = gv = 3, r = 1, K = 8, Z0 = 0.25, hence G = 9.
Because noise and initial conditions were unspecified, examples use equal
polarity/velocity noise amplitudes 0.15, 0.45, and 0.9. Each scenario has six
independent tracks: three with P = V = 0 initially and three with aligned
P and V of magnitude 2 at angles 0, 120, and 240 degrees. All positions start
at the origin. The aligned start is an initialization choice, not a computed
equilibrium. All transients are shown; there is no burn-in or path selection.

Euler–Maruyama uses simultaneous old-state drift evaluation, Gaussian increments
with variance dt, dt = 0.005, and duration T = 100. Position uses the old velocity.
These are finite-step illustrative paths, not a convergence-qualified inference.

Run from the repository root using the existing environment (one CPU, small RAM):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/simulate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python -m unittest discover -s wgd_langevin/logistic_feedback -p 'test_*.py' -v
```

CLI options include `--gp`, `--gv`, `--r`, `--k`, `--z0`, `--noise`, `--dt`,
`--duration`, `--seed`, and `--out`. The Python `simulate` function additionally
supports unequal noise amplitudes and custom initial states.

Outputs in `results/`: three six-panel time-coloured trajectory PNGs, an
ensemble overview PNG, and a speed/gate diagnostic PNG. Individual panels use
equal X/Y scaling but independent limits to reveal their geometry. NPZ files
contain T, P, V, X (state arrays shaped time × track × coordinate), sampled
approximately every 0.05 time units including the final step. `metadata.json`
records all parameters, initial states, seeds, software versions, generation
time, and the simulator source hash.
