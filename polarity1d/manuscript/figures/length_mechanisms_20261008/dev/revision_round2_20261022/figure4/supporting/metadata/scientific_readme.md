# Figure 4 round-2 interim: heading decay and geometry

The decay-rate estimand exposes size advantages at weak and intermediate fluctuating input that the censored 1/e persistence time could not resolve. With coexisting fronts, Wave-pinning, Holmes and Spring show size ratios below the single-front geometric expectation. These are simulated activity-heading readouts; they remain distinct from empirical migration summaries, conditional OU persistence estimates, and one-dimensional polarity-reversal measurements. Length models WGD-associated enlargement, motivated in part by greater observed migration persistence.

This is the mandatory review interim. **No new Figure 4 main has been composed.** The final case selection and full S0 criteria are in `interim_summary.json` and the decision tables. Existing simulations, calibrated parameters, original observation selectors, and the four-model Figure 4 campaign are retained. deBelly was absent from that campaign and has no invented Figure 4 arm.

## Changes since 20 October

The primary readout is the rotational heading decay rate, evaluated at six requested lags per view:

$$
k(\ell)=-\frac{\log C(\ell)}{\ell},\qquad
\rho_k(\ell)=\frac{k_{\mathrm{small}}(\ell)}{k_{\mathrm{large}}(\ell)}.
\tag{1}
$$

Both endpoints require polarity magnitude at least 0.1. Whole seeds are resampled with the original 2,000 multinomial draws, preserving the exact original condition bootstrap key and seed. The same weights apply across all lags and views within a condition. Sizes remain independent ensembles. Chemical lags are 250, 500, 1000, 2000, 4000 and 8000; Spring lags are tenfold shorter. All pairs inside the declared view are used; the previous estimator's half-window maximum is not imposed. Off-grid lags are explicitly unavailable rather than temporally interpolated.

Usable lags require endpoint coverage at least 0.5 in both sizes and a 95% bootstrap interval for 1−C excluding zero in at least one size. Primary lag selection uses the shortest usable observed lag with C_small≤0.98 and positive autocorrelations. Selection is held fixed within bootstrap draws. Nonpositive autocorrelation does not define a real logarithmic decay rate. A rate-ratio interval expands to [0,∞] if more than 2.5% of bootstrap draws have an undefined ratio; the defined-draw fraction remains in the table. Near-zero large-cell rates use the requested one-sided lower bound, where applicable. Rate ratios drifting by more than 1.5 across usable lags are flagged non-exponential. The secondary free-intercept ordinary least-squares slope ratio uses observed lags with 0.2≤C≤0.98 and coverage at least 0.5.

Geometric normalization uses actual mean circumference, with circumference uncertainty propagated using the same whole-seed draws:

$$
\rho_{\mathrm{geo}}=\frac{\rho_k}{(\bar{L}_{\mathrm{large}}/\bar{L}_{\mathrm{small}})^2}.
\tag{2}
$$

This ratio is **relative to the single-front geometric expectation**, not a site diffusivity. Chemical geometry is ninefold; Spring's actual geometry is approximately 8.55-fold. The review plot uses the model's actual geometric expectation, with ninefold as the nominal threefold-circumference reference.

Noise tolerance is the first downward crossing of absolute heading autocorrelation 0.5 at lag 2000 (chemical) or 200 (Spring), interpolated linearly in log input. Every designed input is retained, including 0.005, 0.015, 0.02 and 0.05 where designed. Missing, low-coverage and excluded forced inputs cannot be silently bridged: their crossing is bounded. Left/right censoring and bootstrap intervals are explicit. The forced-excluded analysis requires a confirmed false forced flag; missing qualification remains unavailable.

## Review exports

- Root PNGs and matching `supporting/vectors` PDF/SVG: rate ratio, geometric normalization, noise tolerance, and old-versus-new comparison.
- `figure4_reanalysis_supplement.pdf`: all seven views, both observation arms, both forcing choices for noise tolerance, and C/k/coverage for every model–size–input condition.
- `legacy_seven_view_robustness.pdf`: the superseded main robustness panel, preserved verbatim as a supplement.
- `legacy_S0_S12_decisions.pdf` and legacy CSVs: preceding estimand and criteria preserved.
- `supporting/data/condition_lag_metrics.csv`: C, k, confidence limits, coverage, valid/possible pair counts and unavailable statuses.
- `all_lag_size_ratios.csv`: every requested lag and lag-consistency flags. `size_ratios.csv`: selected primary and secondary estimates. `noise_tolerance.csv`: complete fixed-lag tolerance grid.
- `seed_lag_summaries.csv`: cosine sums and pair counts per complete seed; `old_tau_vs_new_rate.csv`: matched old/new comparisons.

Open triangles and arrows denote bounds. Finite interval bounds display both endpoints. Crosses show unavailable comparisons at a labelled plotting baseline, never as point estimates. Open rate-ratio circles/squares flag non-exponential conditions; marker shape gives regime, grey gives forced status. CIs are pointwise rather than simultaneous. Sparse historical two-seed ensembles remain identified in source tables. Lag selection and the C≤0.98 threshold can make the weakest inputs unavailable even when their heading is visibly stable.

Native scalar archives and the 164 MiB whole-seed bootstrap NPZ are retained on the HPC filesystem. The bootstrap NPZ is excluded from Git due to per-file hosting limits. The source tables and all review artwork are portable; the renderer reads the retained summaries. The analysis receipt identifies native scalar hashes, accepted stage manifests, original metric table, and owned analysis outputs. No new large-output group approaches 2 GiB; no simulation outputs were created, enlarged, moved or deleted.

The old τ-based H3 persistence-loss criterion and A6 direct-τ predictor have not been relabelled as new rate-based tests; their dependency is explicitly retained in the new S0 documentation. The broader original decision criteria can therefore remain inconclusive or fail even when particular rate ratios resolve a size effect. No local gain is substituted for the frozen overall H1 gate.

The final confidence-interval audit conservatively expanded C intervals to [−1,1] and k intervals to [0,∞] for 79 finite-C Goryachev lag rows with more than 2.5% undefined whole-seed bootstrap correlations. Only three of these rows had observed coverage at least 0.5. The valid-draw fraction is exported. This correction changed no usable-lag flags, primary rate ratios, noise tolerances or S0 calls. The affected Goryachev supplement pages were regenerated from the corrected CSV.
