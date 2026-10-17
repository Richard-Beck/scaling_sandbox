# Figure-only revision

Reference report: `docs/trial_report.pdf`. Feedback target: `reports/trial_report.pdf`.
Both PDFs are preserved byte for byte; this directory contains standalone figures,
new assay data and reproducible figure builders. No report renderer was invoked.

## Figures for review

- [Figure 1 — competition, gradient reversal and finite counter-cues](figures/figure1.pdf) ([PNG](figures/figure1.png)). The original A–F layout is retained. C and F now compare the native model with physical peak preservation. The former periodic-tracking G–I are replaced by a finite opposing-cue protocol, stimulus-free recovery traces and numerical cue-duration brackets.
- [Figure 2 — mechanochemical polarity](figures/figure2.pdf) ([PNG](figures/figure2.png)). A contains full-horizon Rac, bound MCA, Rho and tension kymographs at L=21 µm, plus the original full-model Rac/MCA polarity trajectories across all lengths. B retains the original settling-time summary for Rac. C shows full-model Rac at default Rho feedback. D retains the Rac settling summary across Rho-feedback strengths. E retains both Rac and MCA for 150-, 300- and 600-s opposing cues. F retains the Rac threshold and the “No eligible initial polarity” region. Feedback labels are “Rho feedback 0×/1×/2×”.
- [Figure 3 — migration summaries and OU dynamics](figures/figure3.png) (PNG only; [caption](figures/figure3_review/figure3_updated_caption.md), [provenance](figures/figure3_review/figure3_updated_provenance.json)). The integrated figure compares nine observed and primary-fit WGD folds, selected model-generated tracks with fitted MSD, WGD−-fixed relative persistence/noise surfaces, paired primary-fit persistence and velocity-noise estimates, and persistence folds across prior/pooling scenarios. [Hosny sensitivity supplement](figures/figure3_review/hosny_sensitivity/hosny_sensitivity.png) is standalone.
- [S1 — saturation and pool controls](figures/supplement1.pdf) ([PNG](figures/supplement1.png)). Adds comparison of three pool controls and explicit checks of active peak mass and inactive mean.
- [S2 — spring/tension model](figures/supplement2.pdf) ([PNG](figures/supplement2.png)). Keeps validation and the reversal-size sweep; adds peak competition at seven resting lengths and two resolutions, a complete example with stimulus/Rac/tension kymographs, protrusion and spring retraction, Rac polarity, and the losing-peak diagnostic.
- [S3 — demoted periodic tracking](figures/supplement3_tracking.pdf) ([PNG](figures/supplement3_tracking.png)). Retains the audited amplitude/phase heatmaps and adds stimulus, polarity and fitted-harmonic traces at L=15 and periods 100 and 400 s.

## What Chiou's control means

[Chiou et al. (2018), Figure 6](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006095)
add enough GTPase to fill the added cytoplasmic domain while preserving protein
content in the peaks. This isolates increased separation from increased peak
content. A literal match of mean inactive concentration alone does not guarantee
this: the same mean can support different peak masses, and a low-density branch
can instead become nearly homogeneous.

The main figure therefore uses a **Chiou-inspired physical peak control**. At a
reference length of 15 µm, take the native competition state at 200 s, or the
conditioned gradient state at 2,000 s. Extend the domain by separating the two
competition peaks or extending the background to the right of the conditioned
front. Copy the physical profiles and fill the added interval with the reference
background active/inactive concentrations. Correct the cell-average interpolation
by preserving integrated active excess above that background. Diffusion constants
and reaction parameters remain fixed. Larger control domains retain the reference
physical grid spacing, with endpoint refinements.

This is a restart experiment with prescribed initial peak geometry, not a claim
that the enlarged state is already a steady state, nor an exact reproduction of
Chiou's periodic-domain calculation. Native comparison states are restarted at the
same assay time. Competition is stimulus-free after restart, and its plotted
latency includes the 160 s between original withdrawal and restart. Gradient
reversal starts its 20-s cross-fade after 5 s on the restarted state; its plotted
latency starts at that cross-fade. Final reversal must be sustained for 20 s.

The control preserves active peak excess mass to floating-point precision and
peak heights to within 0.008%. The inactive background is preserved; its **domain
mean** can differ slightly because the same local trough now occupies a smaller
fraction of the enlarged domain. At L=39 this difference is about 0.10% for
wave-pinning competition, 3.19% for Otsuji competition, and 0.38% for Goryachev
competition. No pool is clamped during a simulation.

S1 separately includes a literal mean-pool match, obtained by adjusting initial
abundance until the inactive mean matches the L=15 value at 200 s for competition
or 1997.5 s for gradient reversal. S1G–I show why that condition alone is not a
peak-preserving control. An additional abundance-law bulk-fill screen is retained
in the source data, including failures to establish polarity; it is not used as
the principal isolated-distance test.

## Main quantitative result

All values below use the revised restart protocol and the 201-point reference
resolution. Competition values summarize three conserved-noise realizations by
the median; gradient reversal is deterministic. Censoring is not a finite event.

| Assay and model | Native L=15 | Native L=39 | Same physical peaks, L=39 |
|---|---:|---:|---:|
| Competition, wave-pinning | 2,645 s | Unresolved at 20,000 s absolute time | 4,740 s |
| Competition, Otsuji | 1,975 s | Unresolved at 20,000 s absolute time | 4,865 s |
| Competition, Goryachev | 415 s | 5,160 s | 560 s |
| Gradient reversal, wave-pinning | 32.5 s | 102.5 s | 47.5 s |
| Gradient reversal, Otsuji | No reversal within window | Initial state ineligible | Eligible initial state, no reversal within window |
| Gradient reversal, Goryachev | 752.5 s | No reversal within 20,000 s | 492.5 s |

Preserving peak content weakens competition slowing in all three systems. It also
reduces wave-pinning reversal slowing and removes the Goryachev large-domain
failure in this protocol. Otsuji remains resistant to gradient reversal. Thus the
1F phenotype changes with abundance/peak scaling and is not a universal distance
clock. The control does not by itself prove that saturation explains every
response: active peak content, geometry and cue exposure also change together
under the native scaling. The original saturation profiles remain essential.

Finite counter-cues in Figure 1 use a 0.1-peak opposing linear gradient. First
condition with a left gradient for 2,000 s, remove it for 2,000 s, and require
stable autonomous polarity before challenging. Apply durations 10–5,120 s, then
recover for 5,000 s with no stimulus. Require reversed first-moment polarity
above 0.05 and final-tail range below 0.005. Refine successful duration brackets
until they are at most 1 s or 1% wide. At this cue strength wave-pinning reverses
throughout the tested sizes, Goryachev reverses only at the smaller tested sizes,
and Otsuji does not exhibit lasting reversal. Ineligible states and failed
challenges remain visible. Brackets are numerical resolution, not biological
uncertainty. Figure 1H deliberately includes failed replacement traces.

Spring-model competition starts with two seeded Rac peaks with a 2% relative
imbalance, local total concentration 2, and a 40-unit dual activation cue. The
assay asks whether established peaks compete after withdrawal; it does not test
spontaneous peak formation. Resolution requires smaller/larger peak excess below
0.1 persistently to the end, with at least 100 time units of confirmation. All
14 cases resolve. Latency versus resting length is nonmonotonic (31, 20, 17, 14,
13, 23 and 39 model-time units), unlike the monotonic reversal-threshold increase.
Both resolutions give the same integer output-time events. The uniform spring
feedback shown in S2F is 80 times the extension, and spring retraction is 0.01
times the extension. Protrusion uses the original authors' boundary force
function. Model time and forces are not calibrated to experimental units.

## Reproduction and checks

Use the existing `scaling_sandbox` Conda environment. SLURM arrays use one core,
2 GB and QoS `small`, without an array concurrency throttle.

1. `sbatch reports/figure_revision/assays.slurm` — completed as 20602793 (122 tasks).
2. `sbatch reports/figure_revision/fixed_peaks.slurm` — completed as 20602925 (60 tasks), after step 1 finishes.
3. `OPENBLAS_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/build_figures.py` — figures only.
4. `OPENBLAS_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/validate.py` — coverage, conservation, peak preservation, matched-mean accuracy, endpoint grid refinement, sampled finite-cue monotonicity, spring resolution agreement and protected PDF hashes.

[Validation summary](validation.json), [figure provenance](figures/provenance.json),
[physical peak-control diagnostics](figures/fixed_peak_control_diagnostics.csv),
[mean-pool diagnostics](figures/pool_control_diagnostics.csv), and
[primary-fit noise values](figures/figure3_primary_noise.csv) accompany the figures.
Complete fields and per-case outcomes are under `data/`; the exact assay manifest
is `manifest.json`. These files represent model experiments and numerical checks,
not new biological measurements or fitted experimental uncertainty.
