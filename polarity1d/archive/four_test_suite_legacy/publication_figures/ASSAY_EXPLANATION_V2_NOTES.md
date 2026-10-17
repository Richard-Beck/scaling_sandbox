# Rebuilt assay explanation figure

**KEEPER — user-approved on 2026-09-29.** The user described this as an excellent
figure and requested that it be marked as a keeper. Preserve the figure, sources,
fields and provenance; do not redesign or replace it without a user request.

The first design was rejected. Its files remain unchanged as `assay_explanation.*`. This replacement, `assay_explanation_v2.*`, follows the visual language of the old Figure 1 in `reports/figure_revision/figures/figure1.png` and its generators: spatial cue kymographs and overlaid active-species profiles coloured by time. There are no flowcharts. Protocol explanation is kept in this caption and notes.

## Files and rebuilding

- `assay_explanation_v2.pdf`: publication vector output, with rasterized cue maps.
- `assay_explanation_v2.png`: 300 dpi preview.
- `assay_explanation_v2.svg`: editable text, vector profiles/response traces, embedded cue maps.
- `make_assay_explanation_v2.py`: editable plotting source.
- `reconstruct_assay_fields.py`: exact-protocol reconstruction source.
- `assay_fields/*.npz`: reconstructed active fields used in the figure.
- `assay_fields/verification.json`: source hashes and endpoint/outcome checks.
- `assay_explanation_v2_provenance.json`: hashes of figure inputs, sources and measured quantities.

Run the plotting source with the existing environment:

```
/home/4473331/.conda/envs/scaling_sandbox/bin/python reports/four_test_suite/publication_figures/make_assay_explanation_v2.py
```

Reconstruct fields only if needed:

```
/home/4473331/.conda/envs/scaling_sandbox/bin/python reports/four_test_suite/publication_figures/reconstruct_assay_fields.py
```

Both use the already-installed `scaling_sandbox` environment. The reconstruction uses one core and 201 cells; the complete reconstruction ran in seconds, below the repository's SLURM thresholds. It writes exclusively into the new figure directory. No parameter search, spring work, campaign mutation, or package installation occurred.

## Caption

Four assays of one-dimensional polarity. Each row shows the applied spatial cue (left), active-species profiles at selected times (centre), and the measured response (right). Cue intensity is shown in grayscale, with white denoting zero and black the maximum applied intensity for that assay. Colours identify the same times in spatial profiles, vertical lines on cue kymographs, and points on response traces. Numerical colour keys give native model time, or time divided by the period in D. Profile scaling uses one common minimum and range across all displayed snapshots within each row; profiles are not individually rescaled. All examples use wave-pinning in its accepted fixed regime, with domain lengths indicated above the profiles.

(A) Two end cues establish competing peaks and taper to zero at time 40. The time to selection, tcomp, is measured from cue withdrawal until the loser/winner peak-excess ratio falls below 0.1 while the winning peak remains above 0.1; this condition must persist for 50 native time units (shaded interval). (B) After preparation under a leftward gradient, the gradient is cross-faded to the right over 20 native time units starting at time zero. The reversal time, trev, ends when first-moment orientation reaches +0.05 with sufficient active contrast for 20 units (shaded interval). The rightward gradient remains on. (C) Established left polarity persists without cue before time zero. A temporary rightward gradient is applied for duration d and then removed. The minimum passing tested duration, d=22, produces lasting reversed polarity during 5,000 native units of cue-free recovery; the neighbouring tested duration d=20 fails (dashed response). The cue and coloured profiles show the passing trial. The final 500-unit scoring interval is shaded. This endpoint measures required cue duration, rather than reversal latency. A symmetric logarithmic time axis is used in C, with a linear central region through |time|=30. (D) A sinusoidally alternating linear gradient drives periodic polarity. The first cycle is hatched and excluded; the first harmonic over P through 6P yields amplitude A and signed phase lag. Positive lag denotes delayed response. P=400 native model units. The dashed response-panel curve shows normalized cue direction; amplitude and lag annotations refer to the fitted response harmonic. Coloured snapshots are from the second cycle.

The endpoint assays use active-signal first-moment orientation, except competition, which uses end-peak excess above the spatial minimum. Periodic orientation compares active-signal means in the outer 20% of the left and right domains, using (right − left) / (|right| + |left|). Eligibility requires appropriate prepared states; orientation endpoints also require sufficient spatial contrast. Preparation, cue amplitudes, shapes, recovery durations and native units differ across the seven models; these examples illustrate measurements rather than prescribe one universal preparation. The full suite retains those model-specific choices.

## Provenance and verification

The accepted numerical records are revision-5 tasks 0015 (competition, L=7 µm), 0056 (gradient reversal, L=19 µm), 0062 (finite countercue, L=21 µm), and 0064 (periodic following, L=21 µm, P=400 native time units).

The existing competition/reversal diagnostics store scalar readouts but not spatial active fields; reaction-diffusion finite trials retain outcomes rather than full fields. Therefore A–C use bounded reconstructions of the exact accepted protocols, mesh, seeds, parameters, preparation schedules, output intervals and original integration horizons. The executable adapter and all polarity1d source hashes were verified against the accepted archived records before reconstruction and checked again afterward. Reconstructed fields are saved separately from campaign outputs.

- Competition's reconstructed endpoint is exactly 47.5 units after withdrawal, matching task 0015.
- Gradient reversal's reconstructed endpoint is exactly 32.5 units from cross-fade onset, matching task 0056.
- Finite d=20 fails and d=22 passes; final orientations match the corresponding task-0062 trial records within 1e-8 absolute / 1e-7 relative tolerance. The shortest passing duration in the full accepted trial list is 22.
- D uses the actual archived periodic spatial field identified by task 0064: `revision3/diagnostics/wave_pinning/task_0024.npz`. Reapplying the accepted harmonic scorer reproduces A=0.8215282370117998 and lag=+0.07595382735180142 cycles, with a reliable fit. No periodic field is reconstructed.
- Cue kymographs are evaluated using the actual `LocalizedPulse`, `_ImmediateGradientCrossfade`, `LinearGradient`, and `AlternatingGradient` functions used by these protocols. No stylized field approximation is used. Kymograph and response axes share time origins and plotted extents within each row.
- Spatial snapshots are sampled at the same times used for cue lines and response markers. Reconstruction snapshots A–C fall on stored output times; the final periodic snapshot at 1.9P is linearly interpolated between archived samples. No time-course smoothing is applied.

The PNG was visually checked after rebuilding for clipping, overlapping labels, profile/cue correspondence, finite-cue withdrawal and recovery, and periodic phase convention. The 9.3 × 8.2 inch design is intended for a broad scientific figure; final journal dimensions and font settings remain a typesetting choice. The examples do not compare all seven models, display length trends, or establish mechanisms. They explain what the assays impose and measure.
