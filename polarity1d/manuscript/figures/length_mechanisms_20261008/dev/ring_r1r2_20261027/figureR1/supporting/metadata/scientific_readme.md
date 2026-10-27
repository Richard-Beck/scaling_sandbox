# Unexpected

- The full seven-diffusivity law does not satisfy a universal 10% agreement check: at Dᵥ=2.5 the relative error is -23.5–-17.2%, and at Dᵥ=5 it is -14.4–-10.9%. The native Dᵥ=10 error is -8.17–-6.79%; the Dᵥ=1000 error is -0.111–0.184%. The per-row stationary θ and isolated-front constants were checked; no parameters were changed.
- The C60 native solver check has k(round 7)/k(round 3)=0.995187 (95% paired seed CI 0.993912–0.997484), a -0.481% difference whose interval excludes exact equality. The seed, Gaussian field, observation window and estimator are identical; the fixed IMEX time discretization differs from native BDF. This small solver difference is reported without adjustment.

## Key numbers per panel

A–B: Dₛ = 0.002, C20/C60, T = 24000, sample-path cadence 20. Analytic k(C20) = 0.00019739209, k(C60) = 2.1932454e-05; ratio = 9. No null-ensemble simulation was added.

C–D: f = 0.005, post-burn window 16000–32000, observation cadence 5, 384 native seeds per condition. The largest absolute C60 departure from the fast-cytosol fixed-drift guide is 0.00537 over the observed lag grid.

E: weak-dose slopes C20: G=2.2824821, C60: G=1.5976335. Higher-dose observations remain plotted.

F: qualitative transport-path schematic; no inferred measurement.

G: exact-ring gain C60/C20 at Dᵥ=1000 = 0.99776192; Dᵥ=2.5: law relative error [-0.23477, -0.17196], Dᵥ=5: law relative error [-0.14386, -0.10939], Dᵥ=10: law relative error [-0.081721, -0.067893], Dᵥ=20: law relative error [-0.044237, -0.040005], Dᵥ=40: law relative error [-0.02363, -0.021292], Dᵥ=100: law relative error [-0.0099058, -0.0070992], Dᵥ=1000: law relative error [-0.0011132, 0.0018406].

H: Native point estimates span 2.44–3. Fast-cytosol estimates span 1.78–2.22.

| D | f | C_label | status | alpha | alpha_ci_lower | alpha_ci_upper | lag | n_seeds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 0.003 | 30 | resolved | 2.4947 | 2.210816 | 2.781329 | 16000 | 384 |
| 10 | 0.003 | 40 | resolved | 2.764718 | 2.361702 | 3.161092 | 16000 | 384 |
| 10 | 0.003 | 50 | resolved | 2.683761 | 2.202432 | 3.180357 | 16000 | 384 |
| 10 | 0.005 | 30 | resolved | 2.612036 | 2.31671 | 2.926668 | 16000 | 384 |
| 10 | 0.005 | 40 | resolved | 2.441706 | 2.053862 | 2.834282 | 16000 | 384 |
| 10 | 0.005 | 50 | resolved | 2.999913 | 2.491836 | 3.5083 | 16000 | 384 |
| 10 | 0.01 | 30 | resolved | 2.861988 | 2.530154 | 3.210864 | 16000 | 384 |
| 10 | 0.01 | 40 | resolved | 2.952646 | 2.5507 | 3.353732 | 16000 | 384 |
| 10 | 0.01 | 50 | resolved | 2.636408 | 2.135994 | 3.112808 | 16000 | 384 |
| 1000 | 0.003 | 30 | resolved | 1.869646 | 1.58166 | 2.168169 | 16000 | 384 |
| 1000 | 0.003 | 40 | resolved | 2.144561 | 1.718953 | 2.547327 | 16000 | 384 |
| 1000 | 0.003 | 50 | resolved | 1.981545 | 1.456065 | 2.491978 | 16000 | 384 |
| 1000 | 0.005 | 30 | resolved | 2.057372 | 1.748339 | 2.346873 | 16000 | 384 |
| 1000 | 0.005 | 40 | resolved | 1.92531 | 1.504725 | 2.32545 | 16000 | 384 |
| 1000 | 0.005 | 50 | resolved | 2.219131 | 1.734532 | 2.702081 | 16000 | 384 |
| 1000 | 0.01 | 30 | resolved | 2.027732 | 1.751396 | 2.314383 | 11315 | 384 |
| 1000 | 0.01 | 40 | resolved | 2.102512 | 1.666497 | 2.544967 | 16000 | 384 |
| 1000 | 0.01 | 50 | resolved | 1.778175 | 1.248023 | 2.287218 | 16000 | 384 |

## Methods

No model retuning. Campaign 2 uses the unchanged round7_noise_20261026.integrate kernel (IMEX Euler dt 0.05, periodic native dx 0.05). Fast runs start from each exact native basal state and use the native per-condition seed, horizon, burn-in, NoiseInput parameters and Gaussian knots. Each completed worker verifies the exact latent SHA256 against its native noise audit before integration. Native histories exist only in RAM and are discarded after the unchanged round-3 region tracker, compact descriptor extraction and exemplar reduction. Replicates 0 and 1 alone retain ≤2000 frames and ≤400 bins, as float16 illustration.

All autocorrelations use the unchanged round-3 endpoint cosine estimator and |P|≥0.1 endpoint rule. H uses the unchanged lag grid, fixed adjacent-size neighborhoods and shared-lag selector nearest correlation 0.5, OLS of log rate against log actual circumference, and 2000 whole-seed bootstrap resamples. H plots only three-size interior neighborhoods (C20/30/40, C30/40/50 and C40/50/60); the two-size endpoint estimates remain in the extraction table but are excluded from the requested three-size plot. H keeps f separate, avoiding a pooling assumption. Bands in D are pointwise 95% whole-seed intervals. Dashed D guides and B/H predictions are closed form, with no simulated null ensembles. A reuses the original fixed drift sample paths.

The deterministic conduction calculation uses isolated-front κ=1.281489794, λ=5.304509755, Δu=0.960849846. Each plotted exact ring supplies its own θ; the prediction line interpolates θ between retained circumference states, without fitting gains. E fits only the declared weak-dose range through the origin, keeping every larger-dose point.

### Solver check

Native Dᵥ=10 is rerun with the same round-7 worker for original replicates 0–7 at C20/C60, f=0.005. k is compared at the round-3 cohort analysis lag, using paired whole-seed bootstrap weights for the k ratio.

| C | n_seeds | lag | native_k | round7_k | ratio | ratio_ci_lower | ratio_ci_upper |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 8 | 16000 | 9.955963e-06 | 9.92388e-06 | 0.9967776 | 0.9875984 | 1.000794 |
| 60 | 8 | 11315 | 1.624306e-06 | 1.616488e-06 | 0.9951868 | 0.9939115 | 0.9974836 |

### Rebuild and provenance

Data generation time: **2026-10-27T21:04:57.744064+00:00**. Producer complete: **True**. Portable rebuild: `polarity1d/manuscript/code/mechanisms/rebuild_figureR1_20261027.py` with no flags. To refresh compact data after campaign completion, run `ring_r1_campaign_20261027.py collect`, then `rebuild_figureR1_20261027.py --extract` under SLURM. Rebuild writes PNG, PDF and SVG at 180 mm width. Extraction hashes, exact source tables, frozen campaign design and submission receipts are beside this file.

Biological scope: size denotes WGD-associated enlargement. These simulated activity-heading and edge-motion readouts remain distinct from measured migration summaries and conditional OU velocity-persistence estimates. No migration-velocity law is assumed.
