# Figure R2 scientific readme

## Unexpected

10 conditions have incomplete execution; pending estimates are × on a separate thin band below the α axes, with no assigned α value.

Campaign execution or extension checks remain unfinished in the extraction provenance; this snapshot is not final inference. The existing dependency finalizer rebuilds R2 again after both producers finish.

4 α estimates are censored. Open triangles and arrows preserve one-sided bounds; paired opposing triangles/arrows show two-sided or unbounded slopes.

Pilot check Goryachev: α ≥ 3.645 (lower_bound); earlier check above 4. Bootstrap resolved fraction 0. No parameter or seed was changed.

Pilot check Wave-pinning: α = 1.64 (resolved); earlier check about 1.6. Bootstrap resolved fraction 1. No parameter or seed was changed.

Pilot check Holmes: α = 0.8958 (resolved); earlier check about 0.9. Bootstrap resolved fraction 1. No parameter or seed was changed.

Pilot check Spring: α = -0.1375 (resolved); earlier check about 0. Bootstrap resolved fraction 1. No parameter or seed was changed.

Wave-pinning per-cell minus per-perimeter α ranges from 0.7211 to 1.042 across 7 matched resolved conditions (earlier check about +1).

Spring per-cell minus per-perimeter α ranges from 0.2705 to 1.108 across 9 matched resolved conditions (earlier check about +1).

Existing Spring event times depend on its changing physical circumference and cannot be regenerated exactly from seeds alone. G therefore retains the logged per-cell outer-size histories for Spring; chemical models use exact regenerated per-perimeter events. Supply arms differ across these model facets.

The shared Spring C40 (size 2) centre retains its original state-dependent per-perimeter hazard. The per-cell outer sizes instead have constant intensity equal to ρ times the exact size-2 basal physical circumference; the centre is not an exactly constant-intensity trajectory.

The β = −1 sample paths retain the original literal Bernoulli arrival transition at interval 20, while B shows the requested continuous-time analytic exponential. This finite-step mismatch is retained.

The prelaunch cadence selection used 46 completed target-model cohorts, including the shared centre and every rate/amplitude combination at C20. Eight high-rate target cohorts were still executing and were excluded from that selection.

Goryachev Fourier reconstruction reaches -0.422766 activity because of its finite mode cutoff. Display values are clipped at zero on the shared linear viridis scale only; signed compact reconstruction data, native Fourier coefficients and all heading calculations remain unchanged.

G Goryachev: large-minus-small event-triggered excess ranges from -0.014276 to -0.000363256 across w=2–64 original frames.

G Wave-pinning: large-minus-small event-triggered excess ranges from -0.338487 to 0.000837668 across w=2–64 original frames; its sign changes with window width.

G Holmes: large-minus-small event-triggered excess ranges from -0.279713 to -0.0306729 across w=2–64 original frames.

G Spring: large-minus-small event-triggered excess ranges from -0.254042 to -0.126403 across w=2–64 original frames.

Two original allocation records remain stuck in failed prolog on hpctpa3pc0067 after successful replacement runs completed. The existing initial controller dependency was repaired to retain every remaining live worker and exclude only those stale records; no simulation or scientific record was cancelled or changed. The exact receipt is scheduler_dependency_repair.json.

The additive cadence audit checked 49 complete target cohorts: 49 passed, 0 failed, and 5 were pending at that audit. Frozen prelaunch selection unchanged: True. Full comparisons remain in record_cadence_postlaunch_validation.csv.

## Key numbers by panel

Data generation time: **2026-10-28T00:32:22.239424+00:00**. Producers complete: **False**. Exact execution/extension status is in provenance.json and coverage.csv.

A–B: λ₀=1/7000, duration 24000, original sample interval 20, physical peak SD 1.5. Widths follow C20:C60=1:3. β=0 curves coincide.

C: Existing-arm replicate 0 only, all retained active Fourier modes, 240 spatial positions, 321 displayed times. Both models use low rate and low amplitude; no example search or rerun.

D: Whole-profile heading arg P, P=M₁/M₀, endpoint |P|≥0.1, 2000 whole-seed resamples. Dashed guide is C20 correlation to the power 1/9; shared display lag axis is 0–12000, and full retained-record curves remain in CSV.

E: 36 model/rate/amplitude combinations; 24 resolved, 8 pending, 4 bounded.

F: 18 model/rate/amplitude combinations; 16 resolved, 2 pending, 0 bounded.

G Goryachev at w=2: 4020 event windows, 402000 random control windows over both outer sizes; four seeds per size.

G Wave-pinning at w=2: 4039 event windows, 403900 random control windows over both outer sizes; four seeds per size.

G Holmes at w=2: 4167 event windows, 416700 random control windows over both outer sizes; four seeds per size.

G Spring at w=2: 2133 event windows, 213300 random control windows over both outer sizes; four seeds per size.

## Methods

All quantities use native model units. Size represents WGD-associated enlargement. Observed migration summaries, conditional OU velocity persistence and simulated activity-heading persistence remain distinct readouts.

The frozen Fourier estimator is unchanged: actual-pair-count uncentered mode autocorrelation with raw numerators and powers pooled over four seeds; γ₁ is ln 2 divided by the first linearly interpolated half-height lag. All four seeds extend together when the half-height lag exceeds record/8, through the existing fourfold cap. Original finite-record statuses are preserved.

E/F use the negative OLS slope of log γ₁ against log actual mean physical circumference over three sizes. The 2000 bootstrap draws independently resample whole seeds at each size and apply the unchanged rate estimator to pooled numerators. Signed OLS coefficients propagate log-rate intervals into α bounds; bootstrap censoring is retained. Pending symbols have their own band and no numerical α coordinate.

G retains the original event-window observable. Circular mean headings are computed from the w original stored frames immediately before and after each event, excluding a sample exactly at the centre. Every finite nonzero circular mean enters without a coherence gate. The widths are 2,3,4,6,8,12,16,24,32,48,64 frames. This corresponds to model-specific native times (chemical cadence 0.5; Spring cadence 0.05). Overlapping windows remain eligible.

For each seed and each width, 100 times as many control centres as scorable events are sampled uniformly from the same run’s fully observed window support. These random windows retain incidental real events and therefore estimate the same-run background, rather than an event-free trajectory. Control timing uses independent deterministic seeds recorded with seed-level event/control sufficient statistics in G_event_heading_change.csv. Neither event timing nor the model trajectory is changed.

The plotted event-triggered excess is the difference between the pooled event mean and pooled random-window mean:

$$\mathcal{E}(w)=\left\langle1-\cos\Delta\theta_w\right\rangle_{\mathrm{event}}-\left\langle1-\cos\Delta\theta_w\right\rangle_{\mathrm{random}}.\tag{1}$$

Confidence bands are pointwise 95% intervals from 2000 whole-seed bootstrap draws. The same seed multiplicities weight event and random sums/counts, preserving their covariance. They are conditional on the sampled random centres; 100 controls per event reduce timing Monte Carlo variability. Excess remains signed and is an event association, not a causal pathway attribution.

Chemical events are regenerated from the three frozen SeedSequence streams, with exact event counts, retained event tails and all RNG states checked against continuation checkpoints. Spring uses its saved per-cell event logs. No new simulation campaign is launched.

One common angular rotation per kymograph pair is determined from the circular mean of the two initial dominant-domain positions. The same angle is applied at every time to both sizes, rounded only to the display grid. Every kymograph uses linear viridis with a colourbar; the same model within the figure shares one activity scale. Negative Fourier ringing is clipped only for display, with signed plotted-source CSVs unchanged. D uses the same low-rate, low-amplitude conditions as C; the fractional-power guide is missing where the measured C20 correlation is nonpositive.

Rebuild portable data with `python polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render`. Read immutable HPC evidence under SLURM using command `all`. The existing dependency finalizer loads the revised renderer after both campaigns and their extensions finish. Frozen simulation sources, seeds, dynamics, amplitudes and rate estimators are unchanged.

## Event-triggered excess endpoints

| Model | Size | w (frames) | Event mean | Random mean | Excess | 95% seed CI |
|---|---:|---:|---:|---:|---:|---|
| Wave-pinning | 20 | 2 | 0.00859908 | 0.00020381 | 0.00839527 | [0.00695151, 0.00977782] |
| Wave-pinning | 20 | 64 | 0.518752 | 0.0377593 | 0.480993 | [0.445555, 0.526516] |
| Wave-pinning | 60 | 2 | 0.00955418 | 0.00032125 | 0.00923293 | [0.00660038, 0.0123401] |
| Wave-pinning | 60 | 64 | 0.176076 | 0.03357 | 0.142506 | [0.122976, 0.162726] |
| Goryachev | 20 | 2 | 0.000409135 | 4.05334e-05 | 0.000368602 | [0.00022396, 0.000486171] |
| Goryachev | 20 | 64 | 0.0318038 | 0.0174249 | 0.0143789 | [0.0125126, 0.0175969] |
| Goryachev | 60 | 2 | 5.577e-06 | 2.31088e-07 | 5.34591e-06 | [5.03271e-06, 5.68065e-06] |
| Goryachev | 60 | 64 | 0.000281614 | 0.00017874 | 0.000102874 | [8.94393e-05, 0.000118663] |
| Holmes | 20 | 2 | 0.0779848 | 0.00056147 | 0.0774234 | [0.0688534, 0.0853973] |
| Holmes | 20 | 64 | 0.58073 | 0.0350643 | 0.545666 | [0.539738, 0.551774] |
| Holmes | 60 | 2 | 0.0476032 | 0.000852774 | 0.0467505 | [0.038967, 0.0554462] |
| Holmes | 60 | 64 | 0.317167 | 0.0512146 | 0.265953 | [0.234243, 0.296606] |
| Spring | 1 | 2 | 0.182096 | 0.000793813 | 0.181302 | [0.167093, 0.196826] |
| Spring | 1 | 64 | 0.843508 | 0.0493789 | 0.794129 | [0.774217, 0.815087] |
| Spring | 3 | 2 | 0.055372 | 0.000472154 | 0.0548998 | [0.0455766, 0.061125] |
| Spring | 3 | 64 | 0.627837 | 0.0373019 | 0.590535 | [0.570328, 0.611886] |
