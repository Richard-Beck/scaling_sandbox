# Revision 5: finer finite-countercue duration sampling

Revision 5 refines the tested duration grid for the finite-countercue assay across all seven models. It preserves the length grid, stimulation regimes, model parameters, endpoint criteria, and periodic scoring from revision 4. The other 616 case rows are reused byte-identically. Of the 154 finite-countercue cases, 133 receive additional duration simulations; 21 known-ineligible cases skip new trials with a recorded reason. Existing trial records are retained in both groups. New outcomes are simulated directly, not interpolated. The manifests and collection map record each source row, updated row, and provenance hash.

The finite-duration grids are defined in [`duration_grids.json`](duration_grids.json). Reaction-diffusion models test 291 candidate durations by subdividing each revision-4 interval into at least 10 segments (and at least one segment per 20 time units); Holmes tests 341 candidates using the same interval rule. Spring tests 111 candidates by subdividing each existing interval into 10 segments. deBelly tests the union of 50 evenly spaced durations from 10 to 1,000 seconds and its nine revision-4 candidates. Each finite-duration row records every tested candidate, outcome, endpoint bracket, and scorer diagnostics. Known-ineligible cases retain old trial records and list planned but unevaluated candidates separately; they have no fabricated outcome rows.

The report retains the three endpoint assays—competition, gradient reversal, and finite countercue—and periodic harmonic amplitude and phase columns. Periodic scoring is unchanged from revision 3: the shared fit over P..6P excludes the first cycle, applies the established reliability mask, and uses the same signed phase convention. Unreliable fits remain in the metrics table with their rejection reasons and are omitted from plotted metric scales. The finite-duration CSV records each tested duration and result. The endpoint is the shortest tested passing duration; no monotonicity is assumed, and any later failure after a pass is flagged in the row-level notes and report manifest.

The campaign audit found nonmonotonic wave-pinning outcomes near the transition: at L=35 the 48-s cue passed and 52 s failed; at L=37, 56 s passed and 60 s failed; at L=39, 60 s passed and 68 s failed. Holmes also has four adjacent pass/fail reversals: at L=26 µm, 6 s passed and 6.5 s failed; at L=28 µm, 6.5 s passed and 7 s failed; at L=30 µm, 7 s passed and 7.5 s failed; at L=36 µm, 9 s passed and 9.5 s failed. These later Holmes trials retained positive orientation and contrast but exceeded the existing 0.005 orientation-range settling limit. The selected endpoint remains the shortest passing duration, and all sampled outcomes are retained without smoothing or threshold changes.

## deBelly protocol note

This revision changes finite-countercue duration sampling only; it does not revise deBelly's selected stimulation regime or endpoint definitions. The historical Figure 2 protocol differs in several ways: it selected challenge onset by a settling criterion rather than using the fixed 10,600-s preparation, screened 50 cue durations from 10 to 1,000 s, and scored target x50 stability rather than the current first-moment orientation, contrast, and settling criteria. Revision 5 adds those 50 durations alongside the nine existing candidates while retaining the current protocol and scorer; it does not claim to reproduce Figure 2's other protocol choices.

For deBelly, the shortest passing duration changes from the revision-4 coarse grid's 250 s to 232.245 s at L=21 µm, and from 800 s to 636.327 s at L=39 µm. These agree with the earlier Figure 2 beta1× Rac x50 thresholds at those lengths; at L=39 µm the MCA x50 threshold is 656.53 s, so the agreement is specific to the Rac observable.

## Figures and data

- [Seven-model dense comparison (PNG)](results/four_test_all_models.png) · [PDF](results/four_test_all_models.pdf)
- [Competition](results/competition.png) · [Gradient reversal](results/gradient_reversal.png) · [Finite countercue](results/finite_countercue.png)
- [Periodic harmonics (PNG)](results/periodic_harmonics.png) · [PDF](results/periodic_harmonics.pdf)
- [Endpoint and provenance table](results/endpoints.csv) · [finite-duration trials](results/finite_duration_tests.csv) · [periodic fit metrics](results/periodic_harmonics.csv)
- [Campaign manifest](manifest.json) · [collected-row provenance](collection_manifest.json) · [plot/report counts](results/report_manifest.json)
