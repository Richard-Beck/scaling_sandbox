# Illustrative 2D tracks from fitted parameters

[Overview](trajectories_overview.png), [state amplitudes](state_amplitudes.png).

Uses the n=2 self-activation 29-um De Belly calibration without parameter changes. The
external calibration cue is absent. Scalar states become 2D vectors with
Euclidean saturation; this isotropic extension and its noise are assumptions,
not calibrated turning dynamics or experimentally measured migration.

The positive deterministic fixed point has |p|=0.0617752 and |v|=0.108669.
Each scenario starts six tracks there, aligned at evenly spaced angles; x=0.
No burn-in or path selection is applied. Integration is Euler–Maruyama,
dt=0.1 s (about tau_p/233), duration 10000 s, sampled every 1 s.

| Noise fraction | sigma_p | sigma_v |
| --- | --- | --- |
| 0.1 | 0.00128079 | 0.000333663 |
| 0.3 | 0.00384238 | 0.00100099 |
| 0.6 | 0.00768477 | 0.00200198 |
| 0.9 | 0.0115271 | 0.00300296 |
| 1.2 | 0.0153695 | 0.00400395 |
| 1.8 | 0.0230543 | 0.00600593 |

Independent isotropic noise uses sigma=f*state_amplitude/sqrt(relaxation_time).
For an uncoupled 2D OU process this gives RMS magnitude f times the reference
amplitude. Coupling changes the actual variance; these are noise conventions,
not measured percentages of variation. The p and v noise amplitudes differ
because their fitted amplitudes and relaxation times differ. No clipping is used.

Displacement is the integral of relative v over seconds, with no calibrated
distance unit. Equal coordinate scaling preserves path geometry; individual
panels have independent limits. NPZ arrays contain t,p,v,x with state shape
10001 x 6 x 2. Metadata records seeds, initial conditions, assumptions and hashes.

Reproduce from repository root using the existing environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/fitted_tracks.py --model self_activation --noise-fractions 0.1 0.3 0.6 0.9 1.2 1.8
