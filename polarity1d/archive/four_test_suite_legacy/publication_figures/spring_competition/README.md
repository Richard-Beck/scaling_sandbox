# Spring competition: initial mechanistic investigation

**PAUSED BY USER (2026-09-29): NO CLEAR LENGTH CONTRAST.** Do not resume, extend,
refine, retry or try to rescue this workflow without an explicit user request.
The observations below are retained historical diagnostics, not a validated
length effect. Read `AGENTS.md` and `WORKFLOW_STATUS.json` before doing any work.

This is a diagnostic investigation of the revision-5 U-shaped curve, not a
replacement campaign. Accepted results and executable sources remain unchanged.

## Figures

- `spring_competition_trajectories.pdf/png/svg`: accepted endpoint curve plus
  reconstructed native active-field trajectories at L0=0.6, 1.8, and 2.6.
  Cyan lines mark the accepted endpoint. The dashed blue polarity trace uses
  tighter solver tolerances at L0=0.6.
- `spring_competition_diagnostics.pdf/png/svg`: peak ratio, asymmetry,
  inactivation coefficient, extension, inactive pool, and withdrawal profiles.
- `spring_inactivation_cap.pdf/png/svg`: exploratory intervention, kept separate
  from native results. Prevents the reaction inactivation coefficient from
  exceeding its withdrawal value after the cue is removed.

## Findings established by the reconstructed trajectories

The three native endpoints exactly reproduce revision 5: 42.5, 12, and 25 native
time units after cue withdrawal, respectively. Mesh, cue, initial conditions,
solver tolerances, 100-unit early response segment and 50-unit persistence
criterion match the accepted spring protocol. No long 20,000-unit tail is needed
to confirm these endpoints. `summary.json` records original source hashes and
endpoint agreement. All simulations use the existing scaling_sandbox Conda
environment, one BLAS/OpenMP thread, and small bounded local runs.

| L0 | Accepted endpoint | Minimum postcue active contrast | Time of minimum |
|---|---:|---:|---:|
| 0.6 | 42.5 | 6.26e-9 | 8 |
| 1.8 | 12 | 0.0408 | 9 |
| 2.6 | 25 | 0.00395 | 15 |

Competition does not consist of one unchanged end peak gradually consuming the
other across these lengths. Cue withdrawal is followed by substantial loss of
contrast and then recovery of polarity. The short domain becomes almost uniform;
the intermediate domain preserves more asymmetry and recovers sooner. The long
domain retains its driven peaks longer, then passes through a deeper low-contrast
state than the intermediate domain before polarity recovers. At L0=2.6, an early
low peak ratio is not an accepted endpoint: the winner is too weak and the ratio
subsequently rises again. The peak-excess floor and persistence criterion matter.

At withdrawal, absolute extension is almost identical (approximately 0.0550) at
all three lengths, giving an inactivation coefficient near 7.402. Thus different
initial tension is not an explanation for the U-shape. Withdrawal active profiles
and inactive concentrations are similar but not identical. In the model equations,
diffusion in normalized coordinates scales as 1/(instantaneous length squared).
The rapid short-domain homogenization is consistent with stronger spatial mixing.
Subsequent tension feedback is not identical: inactivation falls immediately at
L0=0.6, whereas it initially increases at L0=1.8 and 2.6 before falling. Establishing
the relative causal roles of active diffusion, inactive redistribution, and
mechanical relaxation requires additional controlled interventions.

## Numerical qualification

`sensitivity.json` records the tolerance and symmetry checks. At L0=0.6 the native
solver uses rtol=2e-7 and atol=2e-9; the endpoint is 42.5 and the eventual winner
is the initially slightly weaker right side. Tighter tolerances change the winner
to the initially stronger left side and change the endpoint to 51.5 (rtol=2e-9,
atol=2e-11) or 48.5 (rtol=2e-11, atol=2e-13). A mirrored native run resolves at 47
rather than the symmetry-equivalent 42.5. These discrepancies arise after the
spatial contrast has almost vanished; good mass conservation alone does not
establish robustness of symmetry breaking. The 1.8 and 2.6 endpoints remain 12
and 25 at the tighter 2e-9/2e-11 tolerances. The tested checks preserve the
qualitative U-shaped ordering, but they do not establish converged short-domain
latency or winning-side selection. Do not promote the exact short-length branch
as a validated quantitative mechanistic result without further numerical work.

Extended sensitivity traces preserve the original 100-unit fast segment, then
restart with 5-unit output samples, as in the campaign. This allows the 50-unit
confirmation window to finish for endpoints later than 50, without replacing the
original fast segment with a different integration interval.

## Exploratory inactivation intervention

`check_inactivation_cap.py` starts from the unchanged driven withdrawal state and
caps the postcue reaction inactivation coefficient at its withdrawal value.
Mechanical motion, diffusion, dilution, basal activation and reaction feedback
otherwise remain active. The inactivation correction transfers equal amounts
between active and inactive species, preserving total mass. Original source files
are not edited. This diagnostic changes the feedback law after withdrawal, so
it is not an accepted model regime and does not enter the seven-model comparison.

The cap prevents the deep loss of contrast seen in the native 1.8 and 2.6 runs.
At 1.8 the endpoint changes only modestly, from 12 to 10.5. At 2.6 the two peaks
persist substantially longer before one wins: the endpoint increases from 25 to
82.5 native units. This supports a role for tension-dependent inactivation
in the native collapse-and-recovery trajectory; it does not support the simpler
claim that more inactivation necessarily lengthens competition. Direct competition
of maintained peaks and recovery after their collapse are different trajectories.

## Reproduction

Use `/home/4473331/.conda/envs/scaling_sandbox/bin/python`, with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for simulations. Run `investigate.py`,
`check_sensitivity.py`, and `check_inactivation_cap.py`, then the two
`make_*_figure.py` scripts. Source records are revision-5 `endpoints.csv` and the
unchanged revision-2 spring regime/mechanics adapter/author-code functions.
Stored NPZ traces contain actual independently simulated fields, not interpolation
between campaign endpoints. The figure uses fractional spatial position because
the physical domain moves; L0 denotes native resting length, not micrometres.
