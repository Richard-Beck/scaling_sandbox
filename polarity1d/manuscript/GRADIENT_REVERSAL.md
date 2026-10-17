# Finite fixed-physical-slope countergradient

Added 2026-10-06 to the existing five-model manuscript comparison. This is a
second response geometry using the exact same corrected basal-release steady
states, not a new gradient-conditioned preparation. The stimulus is switched
on at pulse onset; there is no preceding left gradient or crossfade.

Every positive pulse applies added stimulus `g*x` over the complete native
domain, then removes it and recovers under native basal input for 10,000 units.
The added stimulus is zero at the left boundary and rises toward the right.
One positive slope g per model is fixed across all lengths and all three
preparation conditions. Peak equals g times length; mean equals half that peak.
Spring uses instantaneous physical position, so slope stays fixed during growth.
deBelly uses its existing nominal cue coordinates. No clipping, normalized
spatial coordinate, special half-domain rule, reconditioning or length-dependent
amplitude is introduced. Basal reactions and native parameters are unchanged.

The existing preparation conditions remain separate because they supplied
distinct selected checkpoints. Their width labels describe preparation only;
the countergradient always spans the complete domain. All 330 original
preparation conditions remain, including ineligible preparations.

## Bounded slope selection

The pilot uses the shortest, middle and longest eligible relative-20%-prepared
cell in each model. The common candidate factors are 1, 4 and 16. For each model,
candidate g equals factor times its nominal patch peak divided by its existing
reference length (15; Spring native reference 1). At the reference length,
the ramp peak therefore equals factor times the nominal patch peak.

Each pilot independently tests zero plus four positive durations: the existing
model's maximum pulse duration divided by 1000, 100, 10 and 1. Recovery and both
outcome gates are identical to the full comparison. Every pilot zero control
must retain its complete initial state; independent array audits recompute all
scalar gates and field flags. Select the lowest factor succeeding at all three
lengths for both gates. If none does, maximize the number of full-strength
successes, then orientation successes, then prefer the lower factor. The same
rule applies to every model. Keep failures explicit; no amplitude retuning
after observing the full-length campaign is allowed.

This is a pragmatic response-permissive calibration, not a fitted biological
gradient or a robustness claim. Pilot successes are based on a coarse duration
grid and do not establish the shortest pulse; the full campaign determines that.

The completed 225-trial pilot passed independent native-array audit. All five
selected slopes achieved both gates at all three sampled qualifying lengths:

| Model | Selected slope g | Ramp peak at reference length | Selected factor |
| --- | ---: | ---: | ---: |
| Wave-pinning | 0.0066666667 | 0.1 | 1 |
| Goryachev | 0.0266666667 | 0.4 | 4 |
| deBelly | 0.24 | 3.6 | 1 |
| Holmes 3 | 0.0666666667 | 1 | 1 |
| Spring | 4 | 4 | 1 |

Slopes use each solver's native stimulus and length units; they are not
cross-model calibrated physical strengths. No fallback selection was needed.
The full campaign completed at 2026-10-06 16:39 UTC: 74,295 trials, all 291
qualifying preparations successful under both gates (582 successful readouts),
78 preparation-ineligible readouts, no execution errors or unresolved sampled
boundaries. Five refinement rounds resolved all sampled fail-to-pass boundaries
to the existing tolerance. All no-cue controls and the independent complete
native-array audit passed. Closed-form OLS independently verified all 24
eligible fits among the 30 series.

The selected native-field diagnostic passes for 489 of 582 scalar successes.
deBelly passes it at only 3 of 36 selected endpoints per gate; its six series
therefore have insufficient stable-field points for descriptive fits. The
remaining failing flags occur in Wave-pinning and Holmes 3. These failures do
not change the retained scalar reversal criteria; they remain visible as open
markers and excluded fit members. Scalar persistence does not establish a
verified full-native-field steady state.

See the [completed supplement](figures/supplement_gradient_response/index.html),
[full readouts](data/runs/fixed_slope_reversal_20261006/results/endpoints.csv),
[trial audit](data/runs/fixed_slope_reversal_20261006/results/independent_audit.json),
and [fit validation](data/runs/fixed_slope_reversal_20261006/results/fit_validation.json).

## Shared machinery and outcomes

`code/gradient_reversal.py` changes only the cue constructor and campaign
manifest/entry point in `code/core/corrected_reversal.py`. Exact restarts,
checkpoint hashes, fresh no-cue controls, dense/log duration grids, nonmonotonic
local refinement, recovery, independent full-trial audit and both outcomes are
reused. The orientation and 1%-scalar-strength replacement gates, final-500
sample coverage and separate 0.1% native-field flag retain [their definitions](PROTOCOL.md).
Duration excludes both preparation and recovery. This is persistent reversal
after a finite countergradient, not time to a newly verified steady state.

The independent input audit additionally enforces a single positive model slope
and g-times-length nominal peaks at every preparation and length. Model units
and slope magnitudes differ; larger cells receive more total added input.

## Execution and artifacts

Use the existing `scaling_sandbox` Conda environment. Pilot and full campaign
run on SLURM QoS small, one CPU and 3 GB per worker, with no array throttle.
Do not edit frozen executable sources while either campaign runs.

```
python polarity1d/manuscript/code/gradient_reversal.py pilot-init
python polarity1d/manuscript/code/gradient_reversal.py pilot-launch
```

The pilot controller audits results, saves slope selection, then automatically
initializes and launches the full comparison. The full controller refines,
audits, and renders the separate gradient supplement without replacing the
existing patch figures. Files are under `data/runs/fixed_slope_pilot_20261006/`
and `data/runs/fixed_slope_reversal_20261006/`. Both retain frozen source copies
and generated-at status. Final figures appear in `figures/supplement_gradient_response/`.
There is no timed browser refresh. Opening or manually refreshing status loads
the latest available data.

## Continuous application: first recovery of original scalar magnitude

The additional experiment requested on 2026-10-06 keeps the identical selected
g*x gradient on continuously. It uses one exact relative-20% basal-release
checkpoint per model and length. There is no preceding gradient, crossfade,
withdrawal, basal recovery or reconditioning. The outcome is the first recorded
time when rightward orientation P is at least the magnitude of the original
leftward orientation. The finite-pulse 1% band, contrast gate, tail stability
gate and confirmation window are not part of this newly requested readout.
Contrast at the crossing remains available as a diagnostic.

All cases share a maximum application horizon of 20,000 native units. A case
without a crossing is marked not reached within that horizon, rather than given
an invented reversal time or altered slope. Preparations remain excluded before
integration. Full native histories and first-crossing brackets are retained;
the earliest passing sample is the reported duration, with an interpolation
estimate saved separately. Geometric sampling respects deBelly's native 0.5
step. Through the existing Spring wrapper, denser output also restricts adaptive
maximum step; native tolerances and physical parameters remain unchanged.

All 97 eligible integrations completed with no missing cases or errors; 13
preparations were ineligible. An independent audit reconstructs every orientation
from native activity profiles, verifies exact starting fields/checkpoint hashes,
first crossings, censor horizons, conservation and the continual cue.

| Model | Reached original magnitude | First passing duration range |
| --- | ---: | ---: |
| Wave-pinning | 0 / 20 | Not reached by 20,000 |
| Goryachev | 0 / 21 | Not reached by 20,000 |
| deBelly | 12 / 12 | 4–6.5 |
| Holmes 3 | 0 / 23 | Not reached by 20,000 |
| Spring | 21 / 21 | 0.401–0.494 |

The nonpassing models become right-biased, but their forced profiles give a
smaller first moment than their original basal polarized profiles. The
[single PDF](figures/continuous_fixed_slope/continuous_gradient_report.pdf)
shows the length curves, complete scalar histories and shape comparisons.
Forced endpoints are not independently certified steady states.

Fresh integration/audit code is `code/continuous_gradient.py` and
`code/continuous_gradient_audit.py`; the final rendering extension is
`code/continuous_gradient_report.py`. Original frozen numerical inputs and the
initial rendering are retained under `data/runs/continuous_fixed_slope_20261006/`.
All integrations ran on SLURM QoS small, one CPU/3 GB each, without an array
throttle. To reproduce, choose a fresh destination:

```
python polarity1d/manuscript/code/continuous_gradient.py init --out NEW_RUN
python polarity1d/manuscript/code/continuous_gradient.py launch --out NEW_RUN
python polarity1d/manuscript/code/continuous_gradient_report.py --out NEW_RUN
```

Run the final report command after the automatic execution/audit finalizer has
completed. Status and records expose generation time and producer state without
timed refresh.
