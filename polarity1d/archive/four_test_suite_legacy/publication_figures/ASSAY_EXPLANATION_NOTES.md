# Assay explanation figure — draft

`assay_explanation.pdf` and `assay_explanation.svg` are vector outputs. The SVG retains editable text. `assay_explanation.png` is a 300 dpi preview. The editable generator is `make_assay_explanation.py`.

Rebuild using the existing environment, without installing packages:

```
/home/4473331/.conda/envs/scaling_sandbox/bin/python reports/four_test_suite/publication_figures/make_assay_explanation.py
```

The figure is 7.2 × 8.8 inches, intended for full page / two-column width. It is an initial publication draft, not a final typeset figure. The left side shows schematic preparation → challenge → readout. The right side uses actual accepted results, all from wave-pinning so the measurement examples remain easy to compare. These examples illustrate assay definitions, not an across-length result or model ranking. The spatial profile in the readout key is illustrative.

## Data provenance

The four records are revision-5 `task_0015.json` (competition, L=7 µm), `task_0056.json` (gradient reversal, L=19 µm), `task_0062.json` (finite countercue, L=21 µm), and `task_0064.json` (periodic following, L=21 µm, P=400 native model time units). Their archived diagnostic paths are followed directly. Competition and reversal traces are from revision 4; the periodic trace is from revision 3. Those records are the retained non-countercue results used by revision 5. All source records/traces are hashed in `assay_explanation_provenance.json`.

Competition and reversal traces are shown only near the endpoint for clarity; no smoothing is applied. The finite-countercue panel shows the directly simulated final orientations for independent tested durations between 10 and 32 native units. The full set of durations is retained in the source JSON. Each point is one trial, not one time sample. No finite-countercue time course is invented: these wave-pinning finite trials do not archive full trajectories. The accepted shortest pass is 22 units. Periodic data are plotted through all six cycles alongside the fitted first harmonic. Phase lag is signed; positive lag is response delay relative to cue. The cue trace shows normalized direction rather than concentration.

## Protocol distinctions to retain in the caption

- Competition uses end-peak excess above the spatial minimum; peaks are measured within the outer thirds of the discretized domain. The loser/winner ratio must remain below 0.1 for 50 native time units while the winner exceeds 0.1. Cue shape, duration, and preparation differ across adapters (notably deBelly); the schematic does not impose one shared preparation time.
- Gradient reversal uses a 20-native-unit gradient cross-fade and sustained first-moment orientation in the new direction, with sufficient active spatial contrast. The plotted wave-pinning crossing is ≥0.05 for 20 units with contrast ≥0.05. Some other adapters use strict > inequalities. Time is measured from the start of the cross-fade, not its end.
- Finite countercue tests susceptibility to a temporary opposing cue. Its endpoint is a tested cue duration, not reversal latency. Eligibility requires established left polarity that persists without cue. Passing requires rightward orientation, sufficient contrast and a settled tail after free recovery. Wave-pinning recovery lasts 5,000 native units and is scored over its last 500 units. deBelly uses a localized countercue and 10,000 s recovery; its preparation additionally checks Rac and MCA x50 stability. Other model-specific protocols are in their accepted adapter records.
- Periodic preparation is model-specific: the original reaction-diffusion models start from their default state, deBelly from a symmetric algebraic fixed point, and spring / Holmes from prepared polarity. Periodic scoring uses the outer-20% edge observable, rather than the first moment used for reversal endpoints. It mean-centres that observable and performs the archived Fourier projection over inclusive P..6P. The first cycle is excluded. The example is a reliable fit, with amplitude 0.8215 and lag +0.07595 cycles. The figure does not display reliability rules or invent categories for missing results.
- Time units are native to each model; spring's length units remain uncalibrated. The figure's wave-pinning examples therefore use “model units.” All assays use each model's selected fixed parameter regime across tested lengths, not one numerically identical cue or initialization across seven models.

## Suggested caption

Four assays challenge polarity in one-dimensional domains. Schematic stages show preparation, imposed challenge and measured outcome. (A) After two end peaks are established, cues are withdrawn and the time to persistent dominance of one peak is recorded. (B) The conditioning gradient reverses over 20 native time units; reversal time begins at the start of this transition and ends at sustained opposite orientation. (C) Independent temporary countercues of different durations are followed by cue-free recovery. The endpoint is the shortest tested duration producing settled opposite polarity after withdrawal. Points show final orientation from individual trials. (D) A periodically alternating gradient produces an oscillating edge orientation. The first harmonic over cycles 2–6 yields amplitude and signed phase lag; cycle 1 is excluded. Blue shading in A and B shows the required persistence interval. Spatial first-moment and edge readouts are distinguished above. Measured examples are wave-pinning results retained in revision 5; model-specific preparation, cue strength and native units are preserved in the full comparison.

## Verification and limitations

The generator runs successfully in the existing `scaling_sandbox` Conda environment, without simulations or package installs. The raster preview was inspected for clipping, overlap, cue/response phase convention, and endpoint annotations. Selected competition/reversal endpoints, finite shortest pass and periodic harmonic values were checked against the archived scorer and records. The draft is readable at its intended full width; shrinking to a single column would require redesign. No existing result or figure was altered or retired.
