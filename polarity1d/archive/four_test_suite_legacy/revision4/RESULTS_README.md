# Revision 4: dense length sampling

Revision 4 samples the existing seven-model, four-assay suite more densely while keeping each model's stimulation regimes, coefficients, and scoring rules fixed. It reuses an existing result only when model, assay, length, period, seed, and cell count all match exactly; other task keys are simulated and are not interpolated. The manifest, collection map, and report manifest record task counts and per-row source provenance.

The length grids are defined in [`grid.json`](grid.json). The reaction-diffusion models and deBelly use the 1–39 μm odd-numbered 2-μm grid plus retained 12- and 18-μm cases. Holmes Model 3 uses 10–50 μm in 2-μm steps plus its retained 15- and 25-μm cases. Spring uses L0=0.6–2.6 in 0.1 increments. Sparse tick labels keep the denser x axes readable while every simulated point remains plotted; retained off-grid cases are called out here and in the figure's data tables.

The overview has three endpoint assays—competition, gradient reversal, and finite countercue—and periodic harmonic amplitude and phase columns. Periodic scoring uses the same revision-3 shared fit over P..6P, excluding the first cycle, with the established reliability mask and signed phase convention. Unreliable fits stay in the metrics table with their rejection reasons and appear in the labeled status strip, not the amplitude or phase scales. The finite-duration and harmonic CSVs retain each tested duration and fit statistic.

## deBelly protocol note

This revision changes length sampling only; it does not revise deBelly's selected stimulation regime or endpoint definitions. The historical Figure 2 protocol differs in several ways: it selected challenge onset by a settling criterion rather than using the fixed 10,600-s preparation, screened 50 cue durations from 10 to 1,000 s rather than the current nine-duration grid, and scored target x50 stability rather than the current shared edge-orientation endpoint. Those protocol differences are retained and documented; no new deBelly protocol changes were introduced for dense sampling.

## Figures and data

- [Seven-model dense comparison (PNG)](results/four_test_all_models.png) · [PDF](results/four_test_all_models.pdf)
- [Competition](results/competition.png) · [Gradient reversal](results/gradient_reversal.png) · [Finite countercue](results/finite_countercue.png)
- [Periodic harmonics (PNG)](results/periodic_harmonics.png) · [PDF](results/periodic_harmonics.pdf)
- [Endpoint and provenance table](results/endpoints.csv) · [finite-duration trials](results/finite_duration_tests.csv) · [periodic fit metrics](results/periodic_harmonics.csv)
- [Campaign manifest](manifest.json) · [collected-row provenance](collection_manifest.json) · [plot/report counts](results/report_manifest.json)
