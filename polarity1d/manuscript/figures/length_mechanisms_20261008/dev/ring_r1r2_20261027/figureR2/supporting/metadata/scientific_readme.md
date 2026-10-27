# Figure R2 scientific readme

## Unexpected

46 conditions have incomplete execution; these are marked × pending, not treated as rate censoring. Existing-arm extension checks must finish before this is final inference.

4 α estimates are censored; open triangles/arrows preserve one-sided bounds, and paired opposing open triangles/arrows identify two-sided or unbounded slopes. No horizon caps or omitted bounds.

Pilot check Goryachev: α ≥ 3.645, status lower_bound; earlier check above 4. Bootstrap resolved fraction 0. No parameter or seed was changed.

Pilot check Wave-pinning: α = 1.64, status resolved; earlier check about 1.6. Bootstrap resolved fraction 1. No parameter or seed was changed.

Pilot check Holmes: α = 0.8958, status resolved; earlier check about 0.9. Bootstrap resolved fraction 1. No parameter or seed was changed.

Pilot check Spring: α = -0.1375, status resolved; earlier check about 0. Bootstrap resolved fraction 1. No parameter or seed was changed.

Wave-pinning per-cell versus per-perimeter α comparison remains pending

Spring per-cell versus per-perimeter α comparison remains pending

Existing Spring event times depend on model state through its changing physical circumference and cannot be regenerated exactly from seeds alone. Panel G therefore uses the new per-cell Spring outer sizes only, with the shared centre unavailable; chemical panels use regenerated existing-arm events.

Spring fixed per-cell intensity is ρ times the exact basal physical circumference at size 2. The shared centre keeps its original state-dependent per-perimeter hazard, as required by the no-rerun instruction; it is not an exactly constant-intensity centre trajectory.

The β = −1 sample paths retain the original round-7 literal Bernoulli arrival transition at interval 20. Panel B is the requested continuous-time analytic exponential, so this tiny finite-step discrepancy is explicit.

Cadence validation used the 46 completed target-model cohorts available before launch, including the shared centre and every rate/amplitude combination at C20. Eight high-rate target-model cohorts were still executing and were not used to choose cadence.

The finite Fourier cutoff produces small negative ringing in reconstructed Goryachev activity (minimum about −0.423). These band-limited reconstruction values are not negative native solver concentrations. Saved Fourier coefficients, including mode zero M₀, remain unchanged, and no reconstruction values are clipped. Its shared C20/C60 grayscale uses an asinh transform with linear width 1 native activity unit, explicitly keyed; neither ring is normalized separately.

## Key numbers by panel

Data generation time: 2026-10-27T21:10:40.493128+00:00. Current campaign and extension states are retained in provenance.json; each pending condition is explicitly listed in coverage.csv.

A: λ₀ = 1/7000; duration 24,000; original sample interval 20; physical peak SD 1.5. C20/C60 widths are 1:3.

B: Exact continuous-time λ(C) = λ₀(C/40)^β and C(ℓ) = exp[−λ(C)ℓ]. No null ensembles were run.

C: Existing-arm replicate 0 only, reconstructed from all stored active Fourier modes at 240 spatial positions; 321 displayed times. Goryachev low rate/low amplitude; Wave-pinning low rate/medium amplitude.

D: Whole-profile heading θ = arg P, P = M₁/M₀; endpoint |P| ≥ 0.1; 2000 whole-seed bootstrap draws. Guides use the C20 curve directly (dotted) and its 1/9 power (dashed).

E: 36 model/rate/amplitude combinations; 23 resolved, 9 pending, 4 bounded.

F: 18 model/rate/amplitude combinations; 0 resolved, 18 pending, 0 bounded.

G Goryachev: 18260 scored event-window pairs across displayed amplitude/size conditions.

G Wave-pinning: 18423 scored event-window pairs across displayed amplitude/size conditions.

G Holmes: 18459 scored event-window pairs across displayed amplitude/size conditions.

G Spring: 3157 scored event-window pairs across displayed amplitude/size conditions.

## Methods

All times, lengths, rates and activities use native units. The models retain their distinct native spatial and temporal scales. Biological migration persistence, conditional OU estimates and these 1D activity-heading readouts are distinct.

The existing Fourier rate estimator is imported unchanged from the frozen Poisson campaign: uncentered actual-pair-count mode autocorrelation, raw numerators and powers pooled over four seeds, first linearly interpolated half-height crossing. The record extends together for all four seeds when its half-height lag exceeds record/8, through the original fourfold maximum. Original upper/lower rate statuses are retained.

Mode 1 is sufficient for γ₁ by definition; active mode 0 is retained to compute P. The validation compares cadence strides 1,2,4,8,16,32,64,128 against unchanged full-cadence rates and exact four-seed bootstrap CIs. Selected strides are 64 (Wave-pinning, cadence32) and32 (Spring, cadence1.6). P remains at original cadence .5/.05 to score w=2,5,20 original-frame windows. No inactive modes or full profiles are added to the existing campaign.

For E/F, α is minus the OLS slope of log γ₁ against log actual mean physical circumference over the three sizes. 2000 whole-seed bootstrap draws independently resample the four seeds at each size. The unchanged mode estimator is applied to each resampled pooled numerator. A censored γ gives its mathematical log-rate interval; signed OLS coefficients propagate it to a one-sided or unbounded α interval. Bootstrap censoring is recorded rather than replaced by a finite rate.

G compares circular mean headings in [t−w,t) and (t,t+w] around each forcing event using w original stored frames, excluding any sample exactly at the event. Every finite nonzero circular mean enters, without a coherence gate. The reported value is pooled over scorable events; confidence intervals resample whole seeds. Overlapping event windows are retained. Chemical event regeneration uses the three frozen independent SeedSequence streams and checks total event counts against continuation checkpoints.

R2 kymographs use the first predeclared replicate with no seed search. Fourier reconstruction uses the retained half-cell phase convention directly: the nonnegative-index active coefficients and their conjugate partners reconstruct the spatial signal. Finite retained mode cutoff makes this a band-limited reconstruction, never a rerun. The stored coefficients and their mode-zero contribution are unchanged; small negative reconstructed values are Fourier ringing rather than native concentration negativity. Goryachev uses one shared asinh activity scale with linear width 1; all negative ringing is retained. Nulls and Wave-pinning use linear grayscale. D displays the common unrescaled lag range 0–12,000, with full curves retained in CSV. Its fractional-power guide is undefined at nonpositive measured C20 and those guide values remain missing, with no fitting or smoothing.

Rebuild from compact records: python polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render. Re-extract immutable HPC evidence under SLURM with the same entry point and command all. Campaign manifests, current producer state and code SHA are in provenance.json.

## Event-window sensitivity (readme only)

| Model | Arm | Size | Amplitude | w=2 | w=5 | w=20 |
|---|---|---:|---|---:|---:|---:|
| Wave-pinning | per_perimeter | 20.0 | low | 8.87775e-05 | 0.000426853 | 0.00484563 |
| Wave-pinning | per_perimeter | 20.0 | medium | 0.000901231 | 0.00476105 | 0.0582812 |
| Wave-pinning | per_perimeter | 20.0 | high | 0.00859907 | 0.0391358 | 0.210129 |
| Wave-pinning | per_perimeter | 40.0 | low | 0.000157552 | 0.000801799 | 0.0053133 |
| Wave-pinning | per_perimeter | 40.0 | medium | 0.00205921 | 0.0103284 | 0.0455074 |
| Wave-pinning | per_perimeter | 40.0 | high | 0.0122574 | 0.0372104 | 0.112868 |
| Wave-pinning | per_perimeter | 60.0 | low | 0.000621346 | 0.000894382 | 0.00293211 |
| Wave-pinning | per_perimeter | 60.0 | medium | 0.00144549 | 0.00612533 | 0.030799 |
| Wave-pinning | per_perimeter | 60.0 | high | 0.00955418 | 0.027663 | 0.0907786 |
| Goryachev | per_perimeter | 20.0 | low | 8.36945e-05 | 0.000329769 | 0.00245571 |
| Goryachev | per_perimeter | 20.0 | medium | 0.000219661 | 0.000736123 | 0.00633752 |
| Goryachev | per_perimeter | 20.0 | high | 0.000409129 | 0.000945024 | 0.00644367 |
| Goryachev | per_perimeter | 40.0 | low | 2.40657e-06 | 8.20203e-06 | 5.36352e-05 |
| Goryachev | per_perimeter | 40.0 | medium | 1.73346e-05 | 4.56259e-05 | 0.000245444 |
| Goryachev | per_perimeter | 40.0 | high | 0.000289078 | 0.000463877 | 0.00140404 |
| Goryachev | per_perimeter | 60.0 | low | 4.88012e-07 | 1.72146e-06 | 1.12351e-05 |
| Goryachev | per_perimeter | 60.0 | medium | 2.04157e-06 | 5.38683e-06 | 2.66976e-05 |
| Goryachev | per_perimeter | 60.0 | high | 5.57041e-06 | 1.25545e-05 | 5.45958e-05 |
| Holmes | per_perimeter | 20.0 | low | 0.00351156 | 0.00682558 | 0.0292989 |
| Holmes | per_perimeter | 20.0 | medium | 0.0131352 | 0.0425168 | 0.140143 |
| Holmes | per_perimeter | 20.0 | high | 0.0779848 | 0.19836 | 0.417141 |
| Holmes | per_perimeter | 40.0 | low | 0.00225396 | 0.00584778 | 0.0240176 |
| Holmes | per_perimeter | 40.0 | medium | 0.0148695 | 0.0390179 | 0.121432 |
| Holmes | per_perimeter | 40.0 | high | 0.0598004 | 0.137798 | 0.284145 |
| Holmes | per_perimeter | 60.0 | low | 0.00157357 | 0.00356458 | 0.012852 |
| Holmes | per_perimeter | 60.0 | medium | 0.0130569 | 0.0331344 | 0.0879718 |
| Holmes | per_perimeter | 60.0 | high | 0.0476032 | 0.119417 | 0.240993 |
| Spring | per_cell | 1.0 | low | 0.00527997 | 0.0334693 | 0.241731 |
| Spring | per_cell | 1.0 | medium | 0.0432715 | 0.169049 | 0.488225 |
| Spring | per_cell | 1.0 | high | 0.182096 | 0.386432 | 0.683897 |
