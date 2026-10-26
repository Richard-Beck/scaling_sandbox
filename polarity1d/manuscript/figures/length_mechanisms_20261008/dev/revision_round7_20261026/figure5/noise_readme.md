# Unexpected

The new campaign uses 64 seeds per condition. Numerical values are reported from the frozen protocol without tuning. The forcing is the unchanged repository NoiseInput on the native grid, retaining its spatial zero mode and normalized OU-knot interpolation. The geometric mean of per-seed half-correlation rates is the declared estimator.

The small Holmes fast-cytosol effect is now resolved: the mode-1 size exponent falls by about 0.08 at both strong inputs, with paired seed intervals excluding zero. Wave-pinning shows a larger decrease, about 0.54 at input fraction 0.1 and 0.44 at 0.3. The fast-cytosol Wave-pinning heading exponent at 0.3 is slightly negative. Thus disappearance of the positive size effect is supported, while an exactly zero exponent is not imposed.

Patch count grows approximately with circumference, but the local statistics are not size independent. At threshold 0.65 and input fraction 0.1, C20 to C60 density falls 20.5%, mean width rises 28.2%, and RMS matched-patch speed rises 23.4%. At input fraction 0.3 the changes are −9.6%, +11.1%, and +18.7%. The same directions occur at the other thresholds. The measured speeds, about 0.04–0.05 at fraction 0.1 and 0.085–0.10 at fraction 0.3, exceed the pilot's approximate 0.03 and 0.05. Our declared centroid-overlap estimator and unchanged repository forcing differ from the pilot's unspecified overlap details and standalone forcing, so these differences are reported rather than corrected.

## Exponents and pilot comparisons

| Model | f | Quantity | Native alpha (95% CI) | Fast alpha (95% CI) | Native minus fast (95% CI) | Pilot difference |
|---|---:|---|---|---|---|---:|
| Wave-pinning | 0.1 | mode1 | 0.902 [0.848,0.957] | 0.358 [0.317,0.400] | 0.544 [0.500,0.592] | 0.49 |
| Wave-pinning | 0.1 | heading | 0.717 [0.657,0.777] | 0.211 [0.166,0.257] | 0.505 [0.455,0.555] | 0.65 |
| Wave-pinning | 0.3 | mode1 | 0.458 [0.434,0.482] | 0.020 [-0.000,0.039] | 0.438 [0.425,0.452] | 0.30 |
| Wave-pinning | 0.3 | heading | 0.406 [0.376,0.435] | -0.027 [-0.052,-0.004] | 0.433 [0.417,0.450] | 0.41 |
| Holmes | 0.1 | mode1 | 0.819 [0.759,0.881] | 0.738 [0.681,0.799] | 0.081 [0.046,0.117] | 0.05 |
| Holmes | 0.1 | heading | 0.676 [0.607,0.747] | 0.592 [0.526,0.664] | 0.084 [0.042,0.128] | 0.10 |
| Holmes | 0.3 | mode1 | 0.238 [0.209,0.267] | 0.158 [0.128,0.186] | 0.081 [0.076,0.086] | 0.10 |
| Holmes | 0.3 | heading | 0.129 [0.093,0.164] | 0.055 [0.019,0.089] | 0.074 [0.066,0.081] | 0.18 |

A confidence interval containing zero leaves that native/fast exponent difference unresolved. The 64-seed design is fixed; no seeds were selected by outcome. The local exponent here is the negative least-squares slope across three log sizes and geometric mean rates, rather than a differential exponent at one circumference.

## Threshold sensitivity

All patch statistics use native Wave-pinning fields, with active concentration above 0.65, and components spanning at least 10 cells (0.5 length units). Checks repeat the identical analysis at 0.55 and 0.75; no births or deaths are counted. Density is mean count divided by circumference; width pools all patch widths within each seed; speed is the RMS signed physical centroid velocity of mutual-best-overlap matches. Sampling is every 5 time units. Ambiguous splitting/merging and nonoverlapping regions supply no velocity. Values below are means across 64 seeds.

| f | Threshold | Density C20,C40,C60 | Width C20,C40,C60 | RMS speed C20,C40,C60 |
|---:|---:|---|---|---|
| 0.1 | 0.55 | 0.05809, 0.04822, 0.0465 | 6.598, 7.945, 8.244 | 0.04134, 0.04936, 0.05002 |
| 0.1 | 0.65 | 0.05742, 0.0474, 0.04563 | 6.095, 7.493, 7.812 | 0.04008, 0.04876, 0.04947 |
| 0.1 | 0.75 | 0.05662, 0.04661, 0.04478 | 5.565, 6.994, 7.337 | 0.04003, 0.05069, 0.05124 |
| 0.3 | 0.55 | 0.07334, 0.06711, 0.06634 | 6.324, 6.874, 6.957 | 0.09183, 0.1047, 0.1072 |
| 0.3 | 0.65 | 0.07289, 0.06667, 0.06588 | 5.679, 6.224, 6.31 | 0.0849, 0.09854, 0.1008 |
| 0.3 | 0.75 | 0.07202, 0.06613, 0.06549 | 5.016, 5.522, 5.596 | 0.08054, 0.09446, 0.09637 |

Threshold changes alter absolute widths and detected density. These finite-threshold summaries describe the patch field; they do not establish a threshold-independent point process or identify patch creation mechanisms.

## Numerical and statistical protocol

Wave-pinning runs last 32,000 time units and discard the first 16,000; Holmes runs last 20,000 and discard 4,000. Both use spatial step 0.05 and time step 0.05; all recorded frames are 5 time units apart. Native inactive diffusivities are 10 and 50 respectively; the fast control is 1000. Reaction terms, total-pool concentrations, initial frozen fields and all other kinetics are unchanged. Exact frozen basal inputs, worker source hashes, job IDs, and per-seed compact results are retained under ../noise/. The solver applies explicit native reaction and implicit periodic finite-difference diffusion, with a uniform roundoff correction restoring each species' diffusion-preserved zero mode.

The maximum relative conserved-pool deviation was 2.87e-12; minimum concentration was 0.000701. Native/fast forcing hashes agreed in all 768 pairs.

Four paired time-step-halving controls were run for each model and diffusivity at C40 and both input fractions 0.1 and 0.3, with the identical stochastic input. These controls assess discretization sensitivity of stochastic observables, not long-trajectory pointwise agreement. See EF_step_controls.csv. The native spatial mesh is retained; no stochastic mesh-refinement campaign was performed. Seed intervals therefore exclude spatial-discretization and estimator systematic uncertainty.

| Quantity | Median absolute relative change | Maximum absolute relative change |
|---|---:|---:|
| mode1_rate | 0.064% | 1.542% |
| heading_rate | 0.051% | 4.416% |
| P_abs_mean | 0.012% | 0.810% |

Confidence intervals resample 64 independent seed blocks 10,000 times, jointly across sizes and diffusivities. Rates are extracted per seed, so intervals summarize independent-run variability without treating correlated frames as independent. Plotted mode rates are geometric means. Mode autocorrelations are uncentered because nonzero spatial modes have zero ensemble mean by rotational symmetry. First-frame decay gives a rate lower bound of log(2)/5; nondecay by the longest lag gives an upper bound. The full censor status is retained per seed.

All plotted quantities derive from full-precision fields; float16 kymographs are illustrations only. The post-burn-in Fourier coefficients and coherence are float32/complex64 compact records. Neither heading nor modes are conditioned on coherence. Heading persistence is an activity-vector readout; it is distinct from simulated polarity replacement, observed cell migration, and conditional OU velocity-memory inference.

## Mathematical definitions

$$a_n(t)=\frac{1}{C}\int_0^C u(x,t)e^{-2\pi i n x/C}\,\mathrm{d}x,\qquad q_n=\frac{2\pi n}{C}.\tag{1}$$

$$R_n(s)=\frac{\operatorname{Re}\left\langle a_n(t+s)a_n(t)^*\right\rangle_t}{\left\langle |a_n(t)|^2\right\rangle_t},\qquad \gamma_n=-\frac{\log R_n(s_{1/2})}{s_{1/2}}.\tag{2}$$

The selected lag is the sampled lag closest to correlation 0.5 immediately before or at its first crossing. Mode-1 heading uses the unit complex coefficient in the same autocorrelation. Its coherence is the magnitude of mode 1 divided by the mean mode. Low coherence is retained; it makes heading a less sharply defined activity readout.

$$\alpha=-\frac{\sum_{j=1}^3(\log C_j-\overline{\log C})\,\overline{\log\gamma_j}}{\sum_{j=1}^3(\log C_j-\overline{\log C})^2},\qquad C_j\in\{20,40,60\}.\tag{3}$$

The overbar on log rate averages independent seeds. The native-minus-fast difference uses the same resampled seed blocks.

## Reproduction

Run the producer prepare and submit commands with the existing scaling_sandbox Python, then the separate strongest-input control prepare command. The immutable manifests record their SLURM commands. After every worker has stopped, run round7_noise_analysis_20261026.py analyze, then its validate command. Plot helpers are imported by the Figure 5 renderer. Reproduction writes into a fresh output location; the delivered frozen run must not be overwritten.

Generated 2026-10-26T14:58:37.682854+00:00. Producer COMPLETE; no timed browser refresh.
