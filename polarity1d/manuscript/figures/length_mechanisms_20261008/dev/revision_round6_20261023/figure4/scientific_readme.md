# Unexpected

The four-model mode-collapse hypothesis fails the declared transition-and-stronger rate/slope screens; model-specific agreement is reported separately and the strong-input regime remains unexplained. Exact region moments also do not generally eliminate the remainder under the original native-amplitude convention.

The same 325 eligible seed-condition pairs give two controlled comparisons. Under the original region-total normalization and native-amplitude convention, exact own-mask moments reduce |R| in 7/325 pairs. Under a shared native-total denominator with matched masks, exact moments reduce |R| in 321/325 pairs relative to point moments. These are different conventions applied to the same cohort, not conflicting estimates. The original convention assigns within-region shape changes with motion; the native-total comparison assigns them with amount/shape. Finite-width moments can be attenuated again by native-amplitude scaling built for unit-modulus point regions, so the first result does not establish failure of a domain mechanism. Conversely, complete own masks partitioning all excess reconstruct P algebraically; R then vanishes by construction. Its shrinkage is not causal validation. Contour-mask sensitivity and both point/exact pairs are exported.

The birth-count requirement shifts the takeover comparison well above the smallest ring’s half-memory input: Wave-pinning f=0.05, Holmes f=0.05, Spring f=0.05, Goryachev f=0.025. These are per-model common inputs; far-bin estimates generally decline with size, while near-bin estimates can be nonmonotone. A lower 20τ incidence cannot separate slower from rarer takeover; lifetime probability is unidentified.

## Existing exemplar evidence

Completed conditions: 260/260, two saved profiles each. Phase-only $a_1$ and same-endpoint native heading correlation differ by at most 0.0116053; finite float16/binning and $|P|$ eligibility differences remain recorded. Full 0–32 complex correlations and equal-time spectra are exported. Dominant nonzero mode $n^*$ and wavelength $C/n^*$ are reported in spectrum_dominant_mode_summary.csv and every condition's mode-rate rows.

Across all designed inputs, Wave-pinning: n*=1–3, C/n*=16.7–60 native length units; n*>1 in 4/60 conditions; Holmes: n*=1–1, C/n*=20–60 native length units; n*>1 in 0/60 conditions; Spring: n*=1–1, C/n*=1.04–3.05 native length units; n*>1 in 0/60 conditions; Goryachev: n*=1–6, C/n*=10–60 native length units; n*>1 in 27/80 conditions. A dominant first harmonic does not establish a single detected domain. The separate complex/phase rate check measures whether the amplitude-weighted complex first mode has the same rate as unit-weight heading correlation; their equality is tested, not assumed (mode1_complex_phase_rate_comparison.csv and supplement S4.7b). Among finite transition-and-stronger point-rate comparisons, complex/phase rate medians and full ranges are Wave-pinning: 0.99 (0.922–1.23; n = 24); Holmes: 1 (0.959–1.35; n = 30); Spring: 1.08 (0.962–330; n = 40); Goryachev: 1.22 (0.34–4.25; n = 37). These are diagnostic point ratios without phase-rate bootstrap precision. Sparse endpoint estimates can dominate the full range: Spring C = 3, f = 0.1 selects lag 1600 with only 4 eligible pairs out of 24 under the unchanged nearest-half rule, giving a ratio near 330. This is not evidence for a 330-fold physical amplitude effect. Selected-lag pair counts and coverage are exported; no alternative lag is selected to improve agreement.


## Current targeted campaign

Generated 2026-10-23T22:58:48.104839+00:00. Initial design 200 conditions × 24 seeds; up to 96 per condition where SE ln γ1 >0.1 (unresolved estimates also top up). Same generator, qualified basal, BDF numerics, half-horizon burn-in and full horizons as round 3. Only modes 0–32, |P|, time and circumference are durable. No new native profile history is written. Native integration runs through SLURM, using existing scaling_sandbox Conda.

Targeted precision-analyzed conditions: 200/200. The current precision receipt contains 9984 histories: 128 cohorts at 24 and 72 at 96. At the 96-history cap, 37 remain above SE = 0.1 and 18 have unresolved precision. Allocated-cohort precision audit complete: True; the SE target need not be attained at the cap. See supporting/metadata/modes_submission.json and modes_precision.csv. The independent precision controller updates mode evidence and rebuilds artifacts after the authorized additions. The source round 5 deliverables are unchanged.

The following point-estimate diagnostics use each model's predeclared transition input and all stronger inputs; all-input results are retained separately. These descriptive screens are not significance tests, and the table remains provisional while allocated cohorts are incomplete.

| Model | Rate N | α N | ln-rate RMSE | α–slope RMSE | Rate screen | Slope screen |
|---|---:|---:|---:|---:|---|---|
| Wave-pinning | 24 | 16 | 0.0873 | 0.474 | True | True |
| Holmes | 30 | 25 | 0.121 | 0.471 | True | True |
| Spring | 40 | 18 | 0.213 | 0.278 | True | True |
| Goryachev | 39 | 11 | 0.389 | 0.771 | True | False |

| Model | ln-rate RMSE (3 / 5 / 8) | α–slope RMSE (3 / 5 / 8) |
|---|---|---|
| Wave-pinning | 0.123 / 0.118 / 0.0873 | 0.762 / 0.493 / 0.474 |
| Holmes | 0.126 / 0.12 / 0.121 | 0.442 / 0.398 / 0.471 |
| Spring | 0.253 / 0.241 / 0.213 | 0.963 / 0.503 / 0.278 |
| Goryachev | 0.385 / 0.395 / 0.389 | 0.741 / 0.736 / 0.771 |

Model-specific rate-screen agreement: Wave-pinning, Holmes, Spring, Goryachev. Slope-screen agreement: Wave-pinning, Holmes, Spring. All allocated precision cohorts are terminal; unresolved rates and the 96-history cap remain explicit. The eight-point stencil is primary; three- and five-point sensitivities are exported without replacing it. Agreement of complex-mode turnover with phase-only heading turnover is separately tested, because amplitude weighting and coherence eligibility differ.

## Rebuild

`/home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4/rebuild_figure4.py`

This one entry point reads source data, collects completed analyses, exports every plotted point and rebuilds main PDF/PNG, supplement PDF and prose. Scheduler submission is separate and explicit: `/home/4473331/.conda/envs/scaling_sandbox/bin/python -B polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4/scripts/modes_campaign.py submit` (do not repeat while arrays are active).

Biological scope: length represents WGD-associated enlargement. Simulated one-dimensional activity-heading persistence is distinct from observed migration summaries and conditional OU velocity-persistence estimates.
