# Next round brief: Figures 3 and 4

Oct 23, 2026 · @Richard

## Scope

Running brief for the round after commit `769eac59` on `codex/round3-outline-compact-20261023` (Fig 3 round 5 review edition, Fig 4 round 5). Items get added as review goes on. Not sent yet.

- Fig 3: C–D don't show the thing that does the work. D compares drift to a prediction built from a susceptibility map the reader never sees. Fix that, and test the map away from the edge.
- Fig 4: restructure around a description of the strong-input regime that doesn't need domain detection. Keep the region-based diagnostics visible as the evidence for where the domain picture works.
- Figs 1 and 2: style only. One plot per lettered panel (see Style). No new analysis.
- Fig 3 work is deterministic and cheap. Fig 4 gets one targeted, low-storage noise rerun that stores Fourier modes only. No other new noise campaign.
- No model retuning. Existing Fig 3 A–B data stay as they are.

## Figure 3

The drift prediction in D is a linear-response calculation. A weak, brief input at ring position x shifts the whole domain permanently by g(x) per unit dose. Noise is many such inputs, so drift follows from g and the input covariance. The edge pulse in C is one point on g and contributes almost nothing (measured/native gain 0.997–1.004). The figure shows neither g nor where the drift comes from, so D can't be followed. Put g at the centre of C–D.

A and B stay as they are.

### 1. Susceptibility map (new main panel, replaces C time courses)

- Per model, at f = 0.001, C20 and C60 (Spring 1 and 3). All five sizes go in the supplement.
- Plot the map as the response to the standard pulse (Gaussian SD 1.5, τ = 50; Spring 0.15, 5) centred at each position: terminal displacement per dose, signed. Raw g as a thin line in the supplement.
- x axis: arc position in physical units, centred on the domain centre. Mean-state activity profile drawn underneath in light grey on its own axis.
- Light-to-dark model ramp for size, as everywhere else.
- The existing C edge pulses go on the map as points. They should land on the curve.
- Add a separate small panel, not an inset: Wave-pinning C20 and C60 time courses, pulsed and opposite edge on the same axes, showing the permanent shift and both edges moving together. Other time courses to the supplement.
- Drop the integrated-displacement inset. Report terminal gain fold only; that is what squares into D. Integrated response to the supplement.

### 2. Test the map away from the edge

- Same pulse, nominal amplitude, from the settled mean state, matched unpulsed control subtracted. Same numerics as round 5.
- Pulse centres: domain centre; halfway from centre to edge; edge ± 1 and ± 3 pulse SD; a quarter of the way into the basal arc; the antipode.
- Every model, C20 and C60 (Spring 1 and 3), f = 0.001. Track both edges, centroid and width.
- Plot measured terminal displacement per dose as points on the map from item 1. Disagreement is a result; don't tune.

### 3. Where the drift comes from

- Split the predicted Dₛ into per-position contributions, c(x) = g(x) · Σₖ gₖ ∫K(x,k) dv, so that Σ c = Dₛ. Signed; negative values stay.
- Draw c(x)/Dₛ as its own panel directly under the map, same x axis, C20 and C60.
- Report in the legend the fraction of Dₛ from within 2 input correlation lengths of each edge, from the plateau interior and from the basal arc, per model and size.

### 4. Where linear response stops

- At the right edge, every size C20–C60 (Spring 1–3), f = 0.001.
- Amplitude: nominal × 2ᵏ from k = −1 upward until terminal gain departs more than 10% from nominal, the domain is replaced, or the peak reaches the panel B continuous-cue peak. Whichever comes first.
- Duration at fixed dose aτ: τ = 12.5, 25, 50, 100, 200 (Spring ÷10).
- Supplement: gain against amplitude and against duration, per model, size ramp. Main text gets the linear-limit amplitude and its size fold.

### 5. Match sizes with B

- B uses 30/37.8/47.6 (Spring 2.5/3.15/3.97); C–D use 20–60. Run items 1 and 2 (edge and centre only) at 37.8 and 47.6 (Spring 3.15 and 3.97) as well. Mean states: settle from the nearest existing basal, same eligibility rules.
- Supplement: B half-turn-time fold (separated-cue plateau) against map-amplitude fold, at the three matched sizes.

### 6. Rework D

- Main D = the current second row (size-normalized drift, predicted against measured), shared-translation prediction only.
- State on the panel or in the legend that the predicted ratio equals the squared edge-gain fold from C (Wave-pinning 0.70² = 0.49).
- Add a separately lettered panel: implied local α from the predicted drift, α = 2 − d ln Dₛ/d ln C (same three-size local fit as Fig 4), against the Fig 4B weak-input α, per model. 1:1 line. This is why D matters; show it.
- Absolute predicted against measured Dₛ moves to the supplement. The independent-edge candidate goes with it; it is refuted by edge co-motion and doesn't need main-figure space.
- Goryachev stays, on its own scale. It fails; that is a result.

### 7. Legend

- Open D with the plain statement: a weak input anywhere shifts the domain by g(x) times dose; noise is many such inputs; drift follows from g and the input covariance; nothing is fitted.
- Say once that the pulse confirms g at the edge and that the map elsewhere comes from the linearized native model, tested by item 2.

## Figure 4

Fig 4 shows that larger rings hold heading longer. At strong input it reports each panel on its own terms and explains nothing. The domain picture breaks down there: the E iii remainder reaches 0.37 (Holmes C60, f = 0.03). Describe the strong regime without detecting domains. Keep the region diagnostics visible, because they show where the domain picture works and where it stops.

### The idea to test

Heading is the phase of the first Fourier mode a₁ of the activity profile (P = a₁/a₀). Its persistence is set by how fast that mode decorrelates. Let γ(q) be the decorrelation rate of the mode at physical wavenumber q = 2πn/C. A ring of circumference C reads its heading at q = 2π/C, so k(C) = γ(2π/C) and α = d ln γ/d ln q at that q.

- Drift regime: one domain fills the ring and drifts. γ ∝ Dₛq², so α = 2.
- Patchy regime: activity is a field of local patches with correlation length ξ. Below q ≈ 1/ξ, modes turn over at a size-independent rate. γ is flat, so α → 0.
- Larger rings read heading at smaller q, so they leave the drift regime at lower input.

Prediction: at fixed model and input, γ(q) from all five sizes and all modes falls on one curve once rings are patchy, and Fig 4 α equals its local log slope. In the drift regime, only n = 1 lines up. If it doesn't collapse, that is a result. Don't tune.

### 1. First test on existing exemplars

- Seeds 0 and 1 of every condition: the saved full profiles (≤ 400 bins × 2,000 frames, float16, on HPC).
- Per frame, Fourier modes aₙ, n = 0–32 (fewer if the bins limit it). Spring: material coordinate, q from mean circumference.
- Per mode: normalized complex autocorrelation, and γ from the lag where it crosses 0.5 (same lag rule as Fig 4). Also the phase-only autocorrelation of a₁; it should reproduce the Fig 4 heading C(ℓ). Report any difference.
- Equal-time spectrum S(q) = ⟨|aₙ|²⟩/⟨a₀²⟩ against q, all sizes overlaid. Report the dominant nonzero mode n\* and λ = C/n\*.
- From the same profiles, recompute the E iii remainder using exact region moments (integrate each region's own mask) instead of point regions. Report whether R shrinks. That says whether the region scheme fails or only the point reduction does.
- Supplement pages plus readme summary. Build everything else without waiting for item 2.

### 2. Targeted rerun storing modes only

- Same models, all five sizes, inputs f ≥ 0.003. Same input generator, burn-in and horizons as the round-3 campaign. No retuning.
- Store a₀–a₃₂ (complex float32) per frame at the analysis cadence, plus |P|. No profiles.
- 24 seeds per condition to start. Top up to 96 where the SE of ln γ(2π/C) is above 0.1.
- Estimate the footprint before launch. Cut stored cadence or mode count before seeds. If SLURM fails, write the scripts and a one-line command for Richard.

### 3. Main layout

One plot per lettered panel; facets across models are fine (see Style). In order:

- **A.** Null kymographs, C20 and C60, faceted by null.
- **B.** Null heading autocorrelation with α annotated, faceted by null.
- **C.** Model kymographs, C20 and C60, at the transition input (same input rule as round 5), faceted by model.
- **D.** Local α against f, faceted by model.
- **E.** Noise-tolerance fold against C. One plot.
- **F.** Dₛ/Dₛ(Cₘᵢₙ) against f, faceted by model.
- **G.** Closure: k\_pred/k\_obs against f at every input, not only f ≤ 0.01. k\_pred = drift + counted takeovers. Unscored takeover angles become an interval (lower and upper bound), never dropped or assumed. Grey line at 1. Faceted by model.
- **H.** Remainder share R against f, one line per size, faceted by model. Replaces the E iii heatmap. Motion and membership shares as lines in the supplement.
- **I.** Takeover probability against size at one transition input per model, not pooled. Use the designed input nearest the smallest ring's f½ with at least 200 certified births at every size. Drop the d ≤ 0 bin. All inputs in the supplement.
- **J.** Mode decorrelation rate γ against q, all sizes overlaid, at a weak and a transition input. Guide slopes 2 and 0. Faceted by model.
- **K.** Fig 4 α against the local slope of γ(q) at q = 2π/C, all sizes and inputs, 1:1 line. Faceted by model.

J and K go in the main figure only if items 1–2 support them. Otherwise they go to the supplement, and the legend says plainly that the strong regime is not explained. G and H stay in the main figure either way.

Supplement: certified birth rates (old E i), E iii shares as lines, pooled takeover, the exact-moment remainder check, spectra S(q).

## Style and deliverables

- **One plot per lettered panel.** Facets showing the same plot across models, sizes or nulls are fine. Plots with different quantities or axes get their own letters. No insets carrying a different quantity. Violations to fix in this round (not exhaustive; check every panel):
  - Fig 1A: cue kymograph, profiles and orientation trace under one letter.
  - Fig 2A: amount, extent and peak. Fig 2F: biosensor and replacement-time rows.
  - Fig 3C: time courses with insets. Fig 3D: two rows of different plots.
  - Fig 4: every panel, A–E.
- Figs 1 and 2: re-letter and re-lay out only. Data unchanged. Update legends and panel callouts to match.
- Standing figure rules apply: title in the text legend, black panel labels, plain black facet labels, no tables as panels, minimal symbols each defined beside its panel, panels sized to content, fixed model colours and light-to-dark size ramps.
- Censored values as open triangles with arrows. Superseded main panels move to the supplement.
- Outputs in new `revision_round6_<date>/figure1/` to `figure4/` next to round 5; keep round 5 intact. Compact CSVs for every plotted point, including the map, contribution and mode curves. One rebuild entry point per figure.
- `scientific_readme.md` opens with an "Unexpected" section. No stop points.
