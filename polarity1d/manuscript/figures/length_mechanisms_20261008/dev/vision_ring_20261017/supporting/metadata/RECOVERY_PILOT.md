# Common active-deactivation ring port: bounded numerical preflight

All 15 predeclared states pass the same frozen observation clocks used by the linear assay. This author preflight permits review of the numerical port; it does not itself launch the full grid. The assay measures the perturbation-induced loss that is repaired relative to an exact paired sham, rather than generic incumbent stability.

The direction and half-contrast core come from the original complete basal state and remain fixed during challenge. The raised-cosine taper uses independently measured left and right radii. Primary active activity is transferred pointwise to its explicit inactive field at ε = 0.001, 0.005, 0.01. There is no later return, state reset, clipping, or parameter change. All native fields evolve. Challenged states are the exact complete fields at cue-on time 10 and retain the original cue during the short assay.

The finite-window coefficient is

$$
D(t)=\int_{\Omega_0}\bigl[a_\mathrm{s}(x,t)-a_\mathrm{p}(x,t)\bigr]\,\mathrm{d}x,\qquad F=1-\frac{D(\tau)}{D(0)},\qquad K_\tau=\frac{F}{\tau}.
\tag{1}
$$

One clock applies across all sizes, separations, basal/challenged stages and linear/ring geometry within each model: WP/Holmes τ = 0.1778279410038923; Goryachev τ = 0.005623413251903491 in native time. These clocks do not calibrate time across models.

## Primary-amplitude results

| Model | Circumference | Stage | Cue fraction | F | K | Actual sham core change |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| wave_pinning | 10 | basal | — | 0.070999557 | 0.39925985 | 0 |
| wave_pinning | 30 | basal | — | 0.13056927 | 0.73424496 | 5.3290705e-15 |
| wave_pinning | 80 | basal | — | 0.13738234 | 0.77255768 | 3.5527137e-15 |
| wave_pinning | 10 | challenged_t10 | 0.325 | 0.035652127 | 0.20048664 | -0.016139368 |
| wave_pinning | 80 | challenged_t10 | 0.5 | 0.13735979 | 0.77243085 | -0.00076903089 |
| goryachev | 10 | basal | — | 0.0028461041 | 0.50611684 | 3.5527137e-15 |
| goryachev | 30 | basal | — | 0.033822281 | 6.0145466 | 0 |
| goryachev | 80 | basal | — | 0.19029083 | 33.839026 | 2.8421709e-14 |
| goryachev | 10 | challenged_t10 | 0.225 | 0.0022710579 | 0.40385755 | -0.0013717003 |
| goryachev | 80 | challenged_t10 | 0.5 | 0.1902826 | 33.837562 | -1.2432876e-05 |
| holmes_model3 | 10 | basal | — | 0.17377085 | 0.97718532 | -2.6645353e-15 |
| holmes_model3 | 30 | basal | — | 0.18756616 | 1.054762 | 0 |
| holmes_model3 | 75 | basal | — | 0.19093508 | 1.0737069 | -4.2632564e-14 |
| holmes_model3 | 10 | challenged_t10 | 0.375 | 0.16257442 | 0.91422316 | -0.021020094 |
| holmes_model3 | 75 | challenged_t10 | 0.5 | 0.19067413 | 1.0722394 | -0.0092748509 |

The basal coefficients increase over the pilot size range. During challenge, positive conditional repair coexists with declining actual sham core amount in these six sampled challenged states. The paired coefficient therefore cannot substitute for actual incumbent loss or for the later spatial history.

## Numerical and input checks

Every sampled positive-time prefix passes the reviewed empirical sensitivity, saturation and amplitude-linearity gates over all three positive amplitudes. The largest empirical sensitivity E is 2.8596e-07; the largest |F| + E is 0.19121824. There are no unavailable kicks. Independently executed sham and zero-kick integrations agree bitwise in all native fields at both tight and tighter solver settings. Conserved Rac and Rho pools, as applicable, pass their relative-drift gate. These are empirical solver-sensitivity checks, not certified error bounds or stability eigenvalues.

An algebraic 201-cell interval → 402-cell ring mirror check verifies the full native RHS (maximum discrepancy 2.06 × 10⁻¹⁵), fractional-volume core weights, tapered kick and conserved transfer. Reflecting a linear end cue of support 3 yields total ring support 6; the retained challenge campaign uses total ring support 3 and is labelled separately. Constant and nonconstant discrete Fourier diffusion modes also pass. This exact operator fixture does not assert that independently prepared basals are pointwise identical.

Before integration, all 702 planned available state inputs were checked for exact source hashes and nonnegative conservative kicks. Every one of the 652 challenged snapshot arrays agrees bitwise with the complete native cue history at time 10. The full bookkeeping remains 867 planned geometry cells: 652 completed challenges, 198 pre-existing overlap exclusions, 17 unqualified Holmes C80 cells. Basal recovery is measured once for each of the 50 qualified model/size preparations and has no future-separation variable.

Evidence: [algebraic preflight](RECOVERY_ALGEBRAIC_PREFLIGHT.json), [complete source/input audit](RECOVERY_INPUT_AUDIT.json), [pilot numerical audit](../../../../../../data/runs/vision_ring_recovery_20261017/results/pilot_audit.json), [45 primary/sensitivity rows](../../../../../../data/runs/vision_ring_recovery_20261017/results/pilot.csv). Frozen producer and launch records are under `data/runs/vision_ring_recovery_20261017/`. Pilot SLURM array: **21126192**, all 15 tasks completed.
