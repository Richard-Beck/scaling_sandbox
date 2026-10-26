# Next round brief: Figures 3–6 (round 7)

Oct 26, 2026 · @Richard

Build the figure set proposed in the narrative. Figs 3–5 are rebuilt around two mechanisms, cytosolic conduction and Goryachev jump suppression; the current Fig 5 becomes Fig 6. Figs 1–2 are unchanged.

## Inputs

Both are in `manuscript_outline/revision_round7_20261026/request/`, next to this brief.

- `target_narrative_round7.md`. The target narrative. The panel tables in Sections 3–5 define Figs 3′, 4′ and 5′.
- `pilot_analyses_round6.md`. Methods and numbers for the pilot runs behind the new narrative. The pilot code is not in the repo. Reimplement from the methods given there, using the repo's model code and frozen states.

The primes (Fig. 4′, Fig. 3D′) are a convenience, used only to keep the proposed figures apart from the round-6 figures that the narrative also cites. Drop them in the outline: Fig. 4′ becomes Fig. 4, Fig. 3D′ becomes Fig. 3D, and so on.

Where your numbers differ from the pilots, report the difference. Don't tune toward them.

## Figures

- **Fig 3′ (A′–F′).** Build it as the panel table describes.
  - D′: use the linearity brackets in pilot Section 5, not the narrative's "3–25%".
  - E′: show gain against size for all four models, with conduction-law curves. Take Wave-pinning C\* from the isolated front, not a fit. Give the source of each model's C\* in the legend.
  - F′: plot Wave-pinning gain against circumference at the seven D\_v values of pilot Section 1.
- **Fig 4′ (A′–E′).** Build it as the panel table describes. The baselines open the figure and α is the last panel. Models sit in separate facets, so don't colour them; within a facet, size is a light-to-dark ramp.
- **Fig 5′ (A′–F′), new.** Build the mechanisms figure as the panel table describes.
  - B′: use one marker shape and pool over f.
  - D′: heading-change distributions for all four models, at a weak input near each model's C20 transition and at strong input (state which, ideally share the weak and strong input values across models), at C20 and C60 (Spring 1 and 3). Use raw heading changes with no coherence cut, binned by physical displacement C·Δθ/2π. Compute them at lags τ, 4τ and 16τ, with τ the model's input correlation time, and rescale displacement by √lag. Facet by model, input and size, with one line per lag, on a log y-axis. Drift collapses across lags; relocations leave a tail at fixed displacement that doesn't collapse and grows with lag. There is no pass criterion: show the panel whatever it shows. In the readme, give the coherence |P| distribution per condition, and say where 16τ is not well below the typical time between relocations.
  - E′: show the patch-field kymographs, plus patch density, width and speed against C, at the pilot threshold. Recheck them at two other thresholds and report that in the readme. Don't report births or deaths.
  - F′: plot Wave-pinning mode rates against q, native and at D\_v = 1000.
- **Fig 6′.** This is the current Fig 5, renumbered and otherwise unchanged. Fig 7′ stays reserved.
- **Figs 1–2.** Leave them as they are.

## Runs

- Wave-pinning: the D\_v sweep for 3F′, and the isolated front (κ, λ) for C\*.
- Holmes and Spring: the inactive-diffusivity sweeps of pilot Section 2, for C\* in 3E′.
- Strong-input noise runs in Wave-pinning and Holmes at C20, C40 and C60, f = 0.1 and 0.3, with native and fast cytosol. Use enough seeds for bootstrap CIs on α. The pilot's two seeds can't resolve differences of 0.1–0.2.
- Changing D\_v tests the mechanism; it is not retuning. Every main result stays at native kinetics.

Make two pilot pngs, which are checks and not manuscript panels:

- Holmes mode rates, native against fast cytosol, with seed CIs.
- The conduction penalty against C/D for Wave-pinning, Holmes and Spring.

## Rules

- Each lettered panel is one plot. Facets across models or sizes are fine.
- There is no supplement. Drop dropped panels; don't move them anywhere.
- Define every plotted quantity on the figure or in the legend.
- Legends open with the result. Titles go in the legend, not on the figure.
- Panel labels are black. Facet labels are plain black.
- Draw censored values as open triangles with arrows.
- `AGENTS.md` applies.

## Deliverables

- outline\_v13, with Figs 3′–6′ renumbered and the outline text following the target narrative.
- For each figure: png, pdf and svg, `main_legend.md`, and a `scientific_readme.md` that opens with an "Unexpected" section. Add compact CSVs and provenance.
- Code goes in `polarity1d/manuscript/code/mechanisms/`, with `_20261026` suffixes. Outputs go under `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026/`.
- Commit everything to a new branch off `122d306b`.
