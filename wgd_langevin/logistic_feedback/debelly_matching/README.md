# De Belly matching with n=2 polarity self-activation

[Report](results_self_activation/REPORT.md), [calibration](results_self_activation/calibration.png),
[held-out countercues](results_self_activation/heldout_countercues.png).
Previous ungated results remain in results_ungated/ and source snapshots in
archive_ungated/. Gated results remain in results/.

    dp/dt = -p/tau_p + beta*v/(1+abs(v)/vs) + a*p*abs(p)/(K*K+p*p) + cue_gain*u(t)
    dv/dt = -v/tau_v + alpha*p/(1+abs(p)/ps)
    p = Rac_x50/L - 0.5
    v = 0.5 - bound_MCA_x50/L

The exact retained Figure 2 half-mass observables are used, without normalization.
Nine parameters are fitted, with n=2 fixed and both noises zero. The new a
is nonnegative (including zero) and K is positive; there is no logistic gate.
Time is seconds; v is relative mechanical polarity, not calibrated cell velocity.
No De Belly simulations are rerun or modified.

## Design

- Reference fit: A/full and C/beta_1x at 29 um, initially p=v=0.
- Input: continuous right-half minus left-half mean of the original cosine cues,
  preserving their timing and taper. A single left cue projects to -0.72;
  the asymmetric two-sided cue to -0.108. Symmetric activation is omitted by
  this reduction. One cue gain is fitted.
- Loss: equal weight per assay and early/late phase, split 500 s after withdrawal;
  samples approximately every 5 s early and 50 s late. Each field is scaled by
  its reference maximum centered amplitude with a 0.02 floor.
- Eight bounded starts optimize eight logarithmic positive parameters and linear
  a in [0,0.05]. K ranges from 0.003 to 0.5 in centered half-mass units.
  Starts are based on the prior ungated fit, including one with a=0. Bounds,
  convergence, and Jacobian singular values are saved; no confidence intervals.
- All countercues and duration-sweep labels remain held out.
- Frozen transfer across all 20 lengths and three feedback settings; conditional
  single-start fits at 17/21/25/33/39 um (1x), and 29 um (0x/2x). The 2x fit has
  only C available. De Belly beta is not equated to reduced beta.

## Numerics and criteria

Event-aligned RK4, max dt=0.5 s. Step-halving checks reference/conditional
trajectories and all reference-length sweep labels. Coordinates are not clipped;
physical-range validity is recorded.

Settling uses the final 1,000-s median, tolerance 0.002, and 1,000-s confirmation.
Reversal requires opposite sign, final state within 0.05 of reflected conditioning,
and polarity at least 0.02 before and after. Primary labels use the retained
De Belly conditioning reference; self-conditioned diagnostics use the reduced
pre-countercue state. Missing thresholds mean no successful tested duration,
including unpolarized or unsettled cases.

## Reproduce

```sh
mkdir -p wgd_langevin/logistic_feedback/debelly_matching/results_ungated
sbatch wgd_langevin/logistic_feedback/debelly_matching/run.slurm
# After completion:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/debelly_matching/summarize.py
OPENBLAS_NUM_THREADS=1 conda run -n scaling_sandbox python -m unittest discover -s wgd_langevin/logistic_feedback/debelly_matching -p test_campaign.py -v
```

One CPU, 3 GB, small QoS. Resume only with unchanged model and inputs; never
reuse checkpoints from a different model. Outputs include parameters, optimization records,
trajectory metrics, all sweep predictions, paired trajectories, plots and hashes.
