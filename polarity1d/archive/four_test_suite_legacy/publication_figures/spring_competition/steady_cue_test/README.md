# Quick test: settle the spring model under the dual cue before release

Follow-up: the full 21-length sweep is in `../steady_cue_full/`. Its dense curve
is nonmonotonic despite the three-point ordering below; use that follow-up for
the current assessment of length dependence.

Only the cue preparation duration changes. The original 201-cell mesh, spring
equations, L0, basal activation, tension feedback, dual Gaussian cue amplitude 4,
width 0.2 and 2% imbalance, initial state, and competition endpoint criterion are
retained. Original campaign data and sources remain unchanged. This is a
three-length diagnostic, not a dense replacement curve.

Preparation runs in 100-unit chunks, with 5-unit output samples. Before release,
the final 100-unit window must have maximum range below 1e-6 at every spatial
point for both active and inactive concentration, and physical length range
below 1e-7. Instantaneous maximum chemical rate must be below 1e-7 and absolute
length rate below 1e-9. The original RHS plus the exact additive cue is used to
evaluate those rates. Concentrations are compared in normalized material
coordinates. Centre translation is a symmetry and is allowed; physical length
itself is checked explicitly. Integration uses rtol=2e-9 and atol=2e-11, matching
the earlier tighter-tolerance diagnostic.

| L0 | Preparation time | Physical length at release | Campaign endpoint, 40-unit prep | Tighter-solver endpoint, 40-unit prep | Endpoint after steady prep |
|---|---:|---:|---:|---:|---:|
| 0.6 | 400 | 0.6901788 | 42.5 | 51.5 | 73.5 |
| 1.8 | 800 | 1.8997831 | 12 | 12 | 80 |
| 2.6 | 800 | 2.6998830 | 25 | 25 | 91 |

All times and lengths are native model units. Both peaks remain established at
withdrawal in all three cases. The competition scorer is unchanged: losing /
winning end-peak excess must stay below 0.1 for 50 units, while the winning excess
stays above 0.1. Free evolution preserves the original first 100-unit segment
sampled every 0.5 units, then restarts for 400 units sampled every 5 units. All
reported crossings occur within the finely sampled first segment.

The U-shaped ordering is absent in this three-point check. The fast intermediate
case in the 40-unit preparation protocol is therefore sensitive to whether the
cue-held mechanical and chemical states have settled. At 40 units, physical
lengths were approximately L0+0.055; the steady preparations instead reach
extensions around 0.0902, 0.0998 and 0.0999. Thus release from the earlier protocol
occurs while the cue-held mechanics are still evolving. The new runs change both
chemical preparation and mechanical preparation; they do not isolate one of those
as the unique cause.

Important numerical qualification: after release, all three steady-prepared cases
pass through an almost uniform state, with minimum active spatial contrast around
1e-12. Exact re-emergence times and winner selection may depend on numerical
perturbations. The reported 73.5/80/91 ordering is an observed diagnostic outcome,
not yet a validated new length-scaling law. This quick test demonstrates that
steady cue preparation qualitatively changes the trajectories and the sampled
endpoint ordering. Further numerical controls would be required before replacing
the accepted campaign or treating the new timings as converged.

`summary.json` records all convergence metrics and source hashes. Individual
`result_L*.json` files retain the preparation convergence history. NPZ files
contain actual preparation and release trajectories. The generator is
`../check_steady_cue.py`, run using the existing scaling_sandbox Python environment
with one BLAS/OpenMP thread. `steady_cue_comparison.pdf/png/svg` shows preparation
length and end peaks alongside the post-release readouts and endpoint comparison.
