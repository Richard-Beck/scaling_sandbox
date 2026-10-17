> Active campaign (2026-09-26): [direct WGD fits with persistence-only differences](../persistence_only/results/report.html). OU-to-polarity fitting and transfer have been abandoned. The material below is historical.

# OU bridge: bounded polarity-to-motility campaign

This implements the agreed steps 1–4. Prospective experiment selection is deferred.
The accepted Gaussian WGD and Hosny Poisson fits are unchanged. Poisson forcing is
not a candidate in this campaign. Main question: which compact dynamical models
describe length-dependent spatial polarity responses and remain consistent with
the literature motility summaries under explicit transfer assumptions?

Start with the [interpretation](results/CONCLUSIONS.md), [results/REPORT.md](results/REPORT.md), the offline interactive
[trajectory report](results/report.html), and [performance matrix](results/performance_matrix.png).
Reports show their generation time; no browser refresh timers are used.

## Candidate definitions and the memory distinction

All deterministic calibration times use a numerical unit of 100 seconds for the
Jilkine implementations and 1000 seconds for De Belly. This is conditioning of
the optimization, not a biological equivalence. The observable p is relative
interpolated half-mass position minus 0.5. The cue u is the continuous right-half
minus left-half mean of the spatial activation input. It enters p directly.

```
Gaussian: dp/dt = -p/tp + gain*u

Memory:   dp/dt = (-p + r*q)/tp + gain*u
          dq/dt = (p-q)/tm,  0 <= r <= 0.995

Coupled:  dp/dt = (-p + r*q/(1+abs(q)/vs))/tp
                   + self_a*p*abs(p)/(self_K^2+p^2) + gain*u
          dq/dt = (p/(1+abs(p)/ps)-q)/tm
```

The coupled variant fixes self_a=0. The self-activation variant fits self_a>=0
with n=2; both remain development candidates. q is latent, not measured migration
velocity and not automatically bound-MCA polarity. Fixing its drive coefficient
to 1/tm removes a latent-state scaling degeneracy; it does not remove an original
feedback degree of freedom. The rescaling equivalence is tested. De Belly MCA
coordinates are retained in the spatial files for future joint-observable work;
the current cross-model campaign fits Rac polarity only.

Eliminating q from the memory model gives instantaneous damping a=1/tp and an
exponential kernel b=-r/(tp*tm), with decay time tm. This is a specified causal
response realization inspired by Mitterwallner et al. (2020), DOI
10.1103/PhysRevE.101.032408. We do not claim their substitute kernel uniquely
identifies a cue-response mechanism. The nonlinear coupled model has an analogous
linearized memory near a stable state, so these candidates need not be distinguishable
under weak perturbations alone.

## Spatial experiments and split

- Spatial models: wave-pinning, Otsuji, Goryachev, LEGI, accepted De Belly port.
- Jilkine calibration lengths: 5, 15, 39 um; held-out lengths: 10, 25 um.
- De Belly calibration lengths: 21, 29, 39 um; held-out lengths: 25, 33 um.
- Fixed physical mesh spacing <=0.2 um (minimum 25 cells); middle-length checks
  halve the spacing and refine temporal integration.
- Calibration: finite left pulse/withdrawal, +/-5% weak modulation around a
  cue-maintained state, and one alternating-gradient period.
- Held out: four-times-longer tracking period, strong gradient reversal, and
  2%-asymmetric two-front competition. The original retained De Belly Figure 2
  single/dual-cue and countercue trajectories are additional untuned transfer checks.
- Cue peak and relative geometry are held fixed across lengths, unlike a fixed
  physical gradient slope whose total amplitude increases with length. Consequently
  this is a controlled new experiment set, not a relabeling of all retained figures.
- Weak probes are around a tonic-gradient state, not necessarily a spontaneous
  polarized equilibrium. Conditioning is 4000 s for the Jilkine models except
  Otsuji (20000 s), and 10000 s for De Belly. The Otsuji extension followed a
  pre-fit conditioning-drift audit. It was not selected using test-fit performance.
- Front number and spatial contrast are retained. A signed coordinate cannot
  distinguish an unpolarized profile from symmetric competing fronts. Competition
  therefore remains a spatial-reduction diagnostic, not a full-profile prediction.
- LEGI zero-mass initial profiles map to p=0 and are explicitly flagged.

No implementation constants of the spatial models are fitted. Their existing
reconstruction caveats continue to apply. De Belly uses the retained Julia-default
activation values and reflection-symmetric extension.

## Calibration and length laws

Eight seeded bounded starts per model/family/calibration length. Bounds in numerical
time units: tp,tm 0.01–1000; r 0–30 (memory capped 0.995); ps,vs 0.003–2;
gain 1e-5–20; self_a 0–5; self_K 0.003–0.5. Maximum 240 least-squares evaluations
per start. Unconverged starts, bounds, singular values and near-best alternatives
are retained; no global optimum or confidence intervals are claimed.

The objective balances protocols and early/late phases, using a common calibration
polarity scale with a 0.02 floor. A fifth observable explicitly fits the difference
between weak-plus and weak-minus trajectories, with its own amplitude scale and
1e-4 floor. This prevents a large conditioning baseline from swamping the weak
response. RK4 is aligned to stimulus events; best fits are checked at half steps.
Coordinates are never clipped to improve agreement.

Step 3 jointly fits all three calibration lengths, comparing constant parameters,
a power law for tp, power laws for tp and tm, and power laws for tp and cue gain.
At most two exponents vary, each bounded -3 to 3; all other parameters are shared.
Reference parameters are reparameterized so the pilot's physical bounds hold at
every tested length, including the held-out lengths (their responses are not used).
The tm law is omitted for Gaussian OU. Three deterministic starts use the three
pilot solutions. A more elaborate law must improve calibration cost by >=10%
to replace the incumbent. This is a declared descriptive rule, not AIC or a
significance test. Held-out lengths and protocols do not select exponents.

## Literature consistency assumptions

Step 4 keeps the selected polarity length law fixed. It scans all ten increasing
pairs among the five tested lengths, a common time conversion per study, and
declared stochastic assumptions. These are effective-length scenarios, not fitted
cell-size measurements. They establish neither size causation nor WGD causation.
The conversion from area or nuclear size to effective length is not identified by
these data; conclusions requiring it remain conditional.

For the Gaussian model, the same p Wiener coefficient applies at both lengths.
For the memory model, the stationary substitute closure is used: B=tp/2 and
state covariance S=diag(B,B*tp/(r*tm)), with driving covariance
Q=-(A*S+S*A-transpose). Q is positive semidefinite for 0<r<1. The r=0 limit is
ordinary OU. This correlated-noise realization reproduces the published kernel
form's velocity covariance. It is distinct from assuming independent white noise
on both states. Camera-integrated increments are sampled with exact stationary
Gaussian covariance; analytical speed/MSD/endpoint means are used where possible.

The nonlinear extension uses isotropic vector p and q, Euclidean saturation,
and independent Gaussian noise on the two coordinates of p only. Noise coefficients
0.05, 0.2 and 0.6 are sensitivity scenarios shared across lengths, not measured
biological noise. Each bank records 2^23 integration steps after a 2^19-step
burn-in, with half-bank variance and residual-mean audits. Windows occupy at most
1/64 of the bank, and identical parameter sets reuse identical banks. Common random
numbers across lengths prevent simulation variation from masquerading as a size
effect. Sampling overlapping windows does
not create independent biological replicates. Failed stationarity diagnostics
prevent a compatible claim. Independent-seed/step-refinement checks are required
before treating candidate grid hits as robust.

Migration velocity is a scalar multiple of p; position integrates that velocity.
The primary hypothesis shares the amplitude and time conversion between populations.
A second hypothesis releases the population amplitude, retaining the common clock
and length law. For nonlinear models that amplitude is a migration readout scale,
not a change in internal noise. Length-only thus has a precise, restrictive meaning.

Xuan and Mallin use the accepted original targets and reading intervals, at each
study's cadence/duration. Jemaà is reported as insufficiently constrained because
speed alone does not identify dynamics. Scores are descriptive discrepancies.
Finite-grid misses are not proofs that no compatible parameters exist.

Hosny retains the actual path/directness/parallel quartiles and six control-source
summaries. The screen introduces a common directional velocity offset for the
two populations; it is a declared nuisance term, not an inferred oxygen mechanism.
It compares baseline shape errors, holds endpoint/perpendicular summaries out,
and predicts knockdowns with unchanged parameters to expose the need for treatment
effects. It does not claim to refit HIF/RhoA mechanisms. Unpaired normoxia is
outside this length-contrast screen. IQR-scaled losses have no supplied sampling
uncertainty threshold, so Hosny receives a descriptive score, not a binary
statistical compatibility decision.

## Reproduction and audit

Use the existing scaling_sandbox Conda environment. Spatial, fitting, scaling
and literature arrays have one CPU per task, small QoS, no concurrency throttle.
Run from the repository root:

```
sbatch --array=0-24 wgd_langevin/ou_bridge/spatial.slurm
# After spatial completion:
sbatch --array=0-59 wgd_langevin/ou_bridge/fit.slurm
# After fitting:
sbatch --array=0-19 wgd_langevin/ou_bridge/scaling.slurm
# After scaling:
sbatch --array=0-19 wgd_langevin/ou_bridge/literature.slurm
sbatch --array=2,7,12,17,22 wgd_langevin/ou_bridge/spatial.slurm --refine
```

Existing result files are resume checkpoints. Do not reuse them after changing
model code or objectives: preserve the old result directory and start fresh.
`report.py` independently recomputes saved pilot objectives. The `archive_initial/`
directory holds excluded development outputs: initial conditioning, a failed RK4
audit, and pre-weak-contrast objective versions. They are never read by the report.

`test_models.py` checks the analytical response, latent rescaling, nesting, cue
projection, half-mass reflection, stationary memory closure and sampled MSD.
The final provenance lists source hashes, package versions and SLURM job IDs.
