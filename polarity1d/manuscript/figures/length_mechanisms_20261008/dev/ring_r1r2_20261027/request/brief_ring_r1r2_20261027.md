# Brief: ring figures R1 and R2

Oct 27, 2026 · @Richard

## Task

Build two new ring figures, R1 and R2. Run the campaigns they need, render them, and push them to the repo. No model retuning.

- **Base:** `codex/round7-mechanisms-20261026` at `802a0ce`. Branch from it as `codex/ring-r1r2-20261027`.
- **Outputs:** `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/figureR1/` and `…/figureR2/`.
- **Code:** `polarity1d/manuscript/code/mechanisms/`, suffix `_20261027`, with one rebuild entry point per figure that regenerates it from compact records.
- **Scope:** R1 and R2 only. Don't touch Figs 1–6 or the outline. Don't move anything into a supplement.
- **Reuse first:** most R1 inputs already exist in round 7 (paths given per panel). Only the runs listed under Campaigns are new.
- **Expected values** in this brief are checks, not targets. Where a result differs, build the figure anyway and list it under Unexpected. Never tune anything to close a gap.

## Figure R1

A drift null, then Wave-pinning only, at native and free cytosol. `R7` below = `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026`.

| Panel | Plot | Data |
| --- | --- | --- |
| A | Drift-null kymographs, C20 and C60 side by side, width ∝ circumference, time upward, greyscale | `R7/figure4/supporting/data/A_drift_C20.csv.gz`, `A_drift_C60.csv.gz` |
| B | Drift-null heading autocorrelation C(ℓ), C20 against C60, analytic: C(ℓ) = exp(−(2π/C)²·D\_s·ℓ), same D\_s as A | Analytic |
| C | Wave-pinning kymographs at f = 0.005, C20 and C60, faceted by D\_v = 10 (native) and D\_v = 1000. Same seed and noise realization in both facets; shared time axis and activity scale | Native: `R7/figure4/supporting/data/D_wave_pinning_f0.005_C*.csv.gz`. D\_v = 1000: new (Campaign 2) |
| D | Wave-pinning heading autocorrelation at f = 0.005, C20 against C60, faceted by D\_v. Dashed guide in each facet: the C60 curve expected if physical drift were size-independent, C20(ℓ)^((20/60)²) | Native: `R7/figure4/supporting/data/condition_lag_metrics.csv`. D\_v = 1000: new |
| E | Wave-pinning permanent domain displacement against edge-pulse dose, C20 and C60; fitted slope (edge gain) drawn as a line through each | `R7/figure3/supporting/data/D_linearity_plot.csv` |
| F | Schematic: ring with one domain, a pulse at one edge, and the two cytosolic paths between the edges (across the domain and around the back) drawn as resistors in series | Drawn |
| G | Wave-pinning edge gain against C, one series per D\_v (2.5 to 1000), points from the exact ring calculation; lines from G = (λ/2)/(1 + C/C\*), C\* = 2D\_v/(κ·Δu·θ(1−θ)), with κ, λ, Δu from the isolated front and θ from the ring states. No fitting | `R7/figure3/supporting/data/conduction_wave_pinning.csv`, `isolated_front.csv`, `R7/figure3/supporting/metadata/isolated_front.json` |
| H | Wave-pinning local α against C at weak input: observed at D\_v = 10 and D\_v = 1000 (points, 95% CI), and the prediction α = 2 + 2(C/C\*)/(1 + C/C\*) at native D\_v (line). Guide at 2 | Native: `R7/figure4/supporting/data/local_alpha.csv`. D\_v = 1000: new |

Panel-specific encoding:

- **G:** colour encodes D\_v (light = slow, dark = fast), not size, since size is on x. Say so in the key.
- **H:** use the weak inputs available at all five sizes (f = 0.003, 0.005, 0.01). One marker shape per D\_v; pool over f only if the three agree within CI, otherwise facet by f.

## Figure R2

Jump nulls, then the event-noise campaign (Campaign 3). Jump nulls relocate the peak to a uniformly random position at rate λ(C) = λ₀·(C/40)^β, labelled by supply × success:

- **β = 1:** supply per perimeter, constant success.
- **β = 0:** supply per perimeter with success ∝ 1/C, or supply per cell with constant success.
- **β = −1:** supply per perimeter, success ∝ 1/C².

| Panel | Plot | Data |
| --- | --- | --- |
| A | Jump-null kymographs, one facet per β, C20 and C60 side by side, rendered as R1A | β = 1, 0: `R7/figure4/supporting/data/A_perimeter_C*.csv.gz`, `A_same_rate_C*.csv.gz`. β = −1: new sample path, same generator and λ₀ |
| B | Jump-null heading autocorrelation, analytic C(ℓ) = exp(−λ(C)·ℓ), C20 against C60, one facet per β | Analytic |
| C | Event-noise kymographs, C20 and C60: Goryachev at low rate, low amplitude; Wave-pinning at low rate, medium amplitude | Campaign 3 exemplars |
| D | Heading autocorrelation for the same two conditions, C20 against C60. Two guides: dotted, no size dependence (C60 = C20); dashed, drift null (C20(ℓ)^((20/60)²)) | Campaign 3 |
| E | Local α against event amplitude, one facet per model (all four), one line per event rate, per-perimeter arm. Guides at 2, 0, −1 | Campaign 3 |
| F | Local α against event amplitude, per-perimeter against per-cell arm, Wave-pinning and Spring, faceted model × event rate | Campaign 3, both arms |
| G | Mean heading change per event, ⟨1 − cos Δθ⟩, against C, one facet per model, one line per amplitude, low rate, per-perimeter arm | Campaign 3 event logs |

Panel-specific encoding:

- **E–G:** colour encodes event rate (E), supply arm (F) or amplitude (G), as a light-to-dark ramp; say which in each key.
- **E–F:** α from the mode-1 decorrelation rate at C20, C40, C60 (Spring 1, 2, 3; definition under Estimators). Where any size's rate is only a bound, draw α as a bound: open triangle with arrow, never a point.
- **G:** Δθ is the change in mean heading between the windows \[t − w, t) and (t, t + w\] around each event, w = 5 frames. Report w = 2 and 20 frames in the readme only.

## Campaigns

### 1. Nulls

- Every null quantity (B panels, D guides, H prediction) is computed in closed form. No null simulation beyond the kymograph sample paths.
- New: one β = −1 sample path at C20 and C60, using the round-7 null generator, its λ₀ and a fixed seed. Don't search seeds.

### 2. Wave-pinning, fast cytosol at weak input (R1C, D, H)

- **Model:** Wave-pinning, D\_v = 1000; all other parameters native.
- **Grid:** C20, 30, 40, 50, 60 × f = 0.003, 0.005, 0.01.
- **Seeds:** the same seeds and noise realizations as the native runs in `polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022`, same count per condition.
- **Horizon and burn-in:** as those native runs.
- **Solver:** the round-7 worker (`round7_noise_20261026.py`), with D overridden. It already pairs native and fast runs on one noise realization.
- **Solver check:** rerun native D\_v = 10 with the same worker for 8 seeds at C20 and C60, f = 0.005. Compare k with the round-3 values, and report it in the readme.
- **Records:** the round-3 compact record per run (polarity vector, contrast, tracked regions, dominant identity). Full profiles for seeds 0 and 1 only, decimated as in round 3, for R1C.

### 3. Event-noise campaign (R2C–G)

**Existing arm (per perimeter).** Richard's running campaign: exponential inter-event times at a rate ∝ C, exponential amplitudes, a fixed physical footprint per event, four models, C20, 40, 60 (Spring 1, 2, 3), 3 rates × 3 amplitudes, 4 seeds. It is on the HPC and not yet in the repo.

- Let it finish, including its extension checks. It already stores the full Fourier decomposition of each profile at high cadence (about 40 GiB), so add nothing to it.
- R2C kymographs: reconstruct from the stored modes. No new exemplar runs.
- Event logs for R2G (time, position and amplitude per event): regenerate from the frozen seeds if the event sequence doesn't depend on model state. Don't rerun existing conditions to get them. If regeneration isn't possible, R2G uses the new arm only, and the legend says so.

**New arm (per cell).**

- **Models:** Wave-pinning and Spring.
- **Change:** the event rate is fixed at its C40 value (Spring size 2) at every size. Everything else is identical: amplitudes, footprint, seeds, horizons and extension checks.
- **Grid:** all 3 rates × 3 amplitudes, at C20 and C60 (Spring 1 and 3). The C40 runs are shared with the existing arm, so don't rerun them.
- **Records:** don't copy the existing arm's output. Store the event log, P(t), and only the low modes the γ estimator needs, at the coarsest cadence it resolves. Choose the mode cutoff and cadence on existing-arm runs: γ\_1 from the reduced record must match the full record within bootstrap CI. Estimate the footprint before launch and aim for under a tenth of the existing arm.

### Storage and submission

- Estimate the footprint first and stay inside quota with margin. Native profiles are deleted after extraction except the exemplars.
- Update `DATA_MANIFEST.md` and `RETENTION.md`.
- If SLURM submission fails from your session, write the submission scripts and a one-line command for Richard. Keep building everything else, and test the pipelines on existing data.

## Estimators

Use the existing repo estimators unchanged; don't write new versions.

- **Heading:** θ = arg P, where P = M1/M0 of the whole activity profile (round 3).
- **Heading autocorrelation:** C(ℓ) = mean cos(θ(t + ℓ) − θ(t)) over pairs with |P| ≥ 0.1 at both ends (round 3).
- **Local α (R1H):** the round-3 three-size estimator: shared lag where every member is nearest C(ℓ) = 0.5, α = −(OLS slope of ln k against ln C), whole-seed bootstrap, 2,000 resamples.
- **Mode decorrelation rate γ\_n (R2E–F):** the estimator behind the event campaign's preliminary spectra, with its upper and lower rate bounds. Use mode 1. α = −(OLS slope of ln γ\_1 against ln C) over the three sizes, using actual mean circumference. Bootstrap over seeds as above.
- **Censoring:** a bound is drawn as an open triangle with an arrow, never as a point, and never omitted. No horizon caps.
- **Units:** native. Say "native units" once in each legend.

## Figure style

- Double-column width, 180 mm; height as needed. Size each panel to its content.
- One plot per lettered panel. Facets of the same plot are fine.
- Facets of one plot share x scaling, and y wherever the ranges allow. No rescaled axes (no displacement/√lag and similar).
- Size encoding: light = small, dark = large, as a ramp of the model's fixed hue (Wave-pinning blue `#2878b5`, Holmes green `#36944c`, Spring gold `#c59a25`, Goryachev red `#c83b3b`). Nulls use a grey ramp. Where colour encodes something else (R1G, R2E–G), say so in that panel's key.
- Don't colour-code models when they already sit in separate facets.
- Panel labels: bold, black, same size and position throughout. Facet labels: plain black.
- No titles on the figure; the title goes in the text legend. No tables as panels. No caption text on the figure beyond axis labels, facet labels, keys and annotations.
- Keep symbols to a minimum. Define each one in a legend inset or beside its panel. Every guide line and band is labelled on the panel or in its key.
- Thin lines (about 1 pt), small markers, recessive grey guides, axes left and bottom only.

## Deliverables

In each of `figureR1/` and `figureR2/`:

- `figureR*_main.png`, `.pdf`, `.svg`.
- `supporting/metadata/main_legend.md`: a title line, then one paragraph per panel. Each states what is plotted and how it is encoded, then the observed result in one sentence. Define k, α, f, D\_v, C\* and λ(C) once, where first used.
- `supporting/metadata/scientific_readme.md`: Unexpected first, then key numbers per panel, then methods (including the R1 solver check and the R2G window sensitivity).
- `supporting/data/`: compact CSVs for every plotted point, curve and bound, with CIs and the lag used.
- Provenance: code commit and campaign manifests.

Then commit code, compact data, figures and metadata, and push `codex/ring-r1r2-20261027`. Large records stay on the HPC, logged in `DATA_MANIFEST.md`.

## Unexpected results

No stop points. If an estimator fails, a condition is missing, a rate is unresolved or a value is far from the checks below, build the full figure anyway. Plot what exists, mark what doesn't, and list each case under Unexpected: what happened, which conditions, and what you checked.

Checks (from earlier pilots and preliminary plots; record agreement or disagreement):

- R1D: at D\_v = 1000, the C60 curve lies close to the dashed guide.
- R1G: law against exact ring gain within about 10%; at D\_v = 1000, gain flat in C.
- R1H: native α at weak input about 2.4–3.4; at D\_v = 1000, about 2.
- R2E: low rate, low amplitude: Goryachev α above 4, Wave-pinning about 1.6, Holmes about 0.9, Spring about 0.
- R2F: per-cell α above per-perimeter α by about 1.