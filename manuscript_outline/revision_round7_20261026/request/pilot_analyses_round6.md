# Round-6 review: pilot analyses, methods and results

Oct 26, 2026 · @Richard

## Scope

This records the pilot analyses run during the round-6 review (2026-10-24 to 2026-10-26), whose results now feed the bold-italic narrative in the [narrative and claims doc](https://claude.ai/artifact/NGu2ELqKdbRrPDSwEMT3BN). It gives methods and results at the level needed to judge and rerun them, not a full protocol. The questions came from the [review discussion doc](https://claude.ai/artifact/BgxkMc1oeTj5QFzVecedCC). None of these runs are in the repository; they are pilots, and every figure produced was a mockup or check, not a manuscript panel.

Code and data sources:

- **Repository:** `Richard-Beck/scaling_sandbox`, branch `codex/round3-outline-compact-20261023` at `122d306b` (round 6, outline\_v12), read-only. Model definitions, frozen basal states, the noise generator and the adjoint kernel (`phase_kernel`) were taken from it unchanged.
- **Local Wave-pinning ring:** a standalone 1D periodic implementation with the repository's parameters (D\_u = 0.1, D\_v = 10, k₀ = 0.067, γ = 1, K = 1, δ = 1, total mass per length 2.2683). Its linear-response gains match the repository's to four digits.
- **Holmes and Spring:** the repository's frozen native modules for Fig. 4 (`figure4_round3_campaign_20261022.native_module`), with inactive diffusivity changed where stated.

All gains are terminal domain displacement per unit dose (amplitude × duration) for a weak Gaussian pulse (SD 1.5; Spring 0.15) centred on one domain edge, computed from the adjoint of the translation mode of the linearized native model at a weak uniform mean input. This is the same quantity as the round-6 Fig. 3 edge gain.

## 1. Wave-pinning edge gain: cytosolic conduction

**Result.** Wave-pinning edge gain falls with ring size because the material a pushed edge consumes has to travel through the cytosol from the far edge. When the cytosol conducts freely (D\_v = 1000) the size effect disappears: gain is 2.640, 2.641 and 2.634 at C20, C40 and C60.

**Methods.**

- **Steady states:** Newton solves of the local ring from the repo's frozen basal states, with continuation in C for rings above C60 and in D\_v for the diffusivity sweep.
- **Gain:** adjoint of the translation mode at a weak uniform mean input, applied to a Gaussian edge pulse. Checked against direct pulse simulations: gain does not depend on pulse duration at fixed dose (durations 12.5 to 200).
- **Isolated front:** a single front with the cytosol clamped, with no ring and no fitting, gives the front speed coefficient κ = 1.2816 and the per-dose push λ = 5.2958 (pinning value v\* = 1.8011, front height Δu = 0.961). The domain covers a fraction θ ≈ 0.31 of the ring at every size.
- **Ring diagnostics:** time-integrated cytosolic perturbation after a right-edge pulse, and the material delivered to each edge against the material it converts.

**The law.** The two cytosolic paths between the edges (across the domain and around the back) act as resistors in series with the front:

```latex
G = \frac{\lambda/2}{1 + C/C^*}, \qquad C^* = \frac{2 D_v}{\kappa\,\Delta u\,\theta(1-\theta)}
```

With native D\_v = 10 this gives C\* ≈ 76; the supply length D\_v/(κΔu) is about 8 units. The Fig. 3–4 rings (C20–60) sit below the crossover.

**Findings.**

- The isolated-front prediction matches the exact ring gain to within 6–9% from C15 to C200. Predicted C60/C20 fold 0.708, exact 0.698.
- 1/G is linear in C, as resistors in series require. Fitted crossover 72, against 76 predicted.
- The time-integrated cytosolic perturbation is piecewise linear between the edges (within 1–3% of its range). Material delivered to each edge matches what that edge converts within 1–4%.
- The far edge keeps about half the post-pulse retraction (0.51 at C20 to 0.46 at C60). This refutes the earlier mixing-time explanation.
- Gain against D\_v collapses on C/D\_v, not C/√D\_v (table below). The remaining 6–9% shortfall from λ/2 = 2.648 depends on D\_v but not on size, consistent with finite front and pulse widths.
- Drift scales as gain squared, so the law predicts α = 2 + 2(C/C\*)/(1 + C/C\*): 2.36–2.91 over C20–60, against 2.48–2.88 in round-6 Fig. 3H. Only the level matches; the measured trend is flatter (Section 6).

| D\_v | Gain at C20 | Gain at C40 | Gain at C60 |
| --- | --- | --- | --- |
| 2.5 | 1.69 | 1.07 | domain lost |
| 5 | 2.03 | 1.48 | 1.16 |
| 10 (native) | 2.28 | 1.88 | 1.59 |
| 20 | 2.45 | 2.19 | 1.98 |
| 40 | 2.54 | 2.40 | 2.27 |
| 100 | 2.60 | 2.54 | 2.48 |
| 1000 | 2.640 | 2.641 | 2.634 |

This test is the basis of proposed panel 3F′ in the narrative.

## 2. Holmes and Spring follow the same conduction law

**Result.** All three drift models obey the same two-resistor law. They differ only in where their native ring sizes sit relative to the crossover C\*, which is why the Wave-pinning size effect is large and the Holmes and Spring effects are small.

**Methods.**

- **Holmes:** the repo's frozen native module, with only the inactive diffusivity D changed (native 50; tested 5 to 5000). Same adjoint edge gain as Section 1, mean input 1e-3, pulse SD 1.5, C20–60.
- **Spring:** the native module hard-codes v diffusivity 10, so a subclass adds (D\_v − 10) × Laplacian to the right-hand side and Jacobian, including the length-derivative term. Inputs from the round-5 edge-pulse states, mean input 4e-4, pulse SD 0.15, sizes 1–3. Tested D\_v 0.1 to 1000.
- **Penalty:** G(D→∞)/G − 1, with D = 5000 (Holmes) or 1000 (Wave-pinning, Spring) as the free-conduction reference, plotted against C/D.

**Findings.**

- **Holmes:** 1/G is linear in C at every D tested. C\*/D is 7.5–8.4 for D = 5–200, so C\* ≈ 400 at native D = 50, about five times the Wave-pinning crossover. Native gain falls only from 0.460 (C20) to 0.420 (C60); at D = 5000 it is flat at 0.4775.
- **Spring:** native gain is flat (0.302, 0.307, 0.305, 0.302, 0.299 at sizes 1–3). Slowing the cytosol brings back the size effect: at D\_v = 1 gain falls from 0.276 to 0.190. Penalty × D\_v collapses for D\_v ≥ 1, with a fitted C\* ≈ 43 (range 30–80) in Spring's units. At D\_v 0.3 and 0.1 the pattern changes, so the law covers D\_v ≥ 1 only.
- **Spring small-ring rise:** with a free cytosol (D\_v = 1000), gain rises about 5% from size 1 to 3 (0.305 to 0.320). This rise does not involve the cytosol and partly offsets the native conduction penalty.
- **Collapse:** the penalty is close to C/(8D) in all three models. Over the Fig. 3 size range, the native penalty is 26–79% for Wave-pinning, 4–13% for Holmes and 1–7% for Spring.

| Holmes D | Gain at C20 | Gain at C60 | C60/C20 |
| --- | --- | --- | --- |
| 5 | 0.351 | 0.207 | 0.59 |
| 12.5 | 0.416 | 0.310 | 0.75 |
| 25 | 0.444 | 0.376 | 0.85 |
| 50 (native) | 0.460 | 0.420 | 0.91 |
| 200 | 0.473 | 0.462 | 0.98 |
| 5000 | 0.477 | 0.478 | 1.00 |

| Spring D\_v | Gain at size 1 | Gain at size 3 | Size 3/size 1 |
| --- | --- | --- | --- |
| 1 | 0.276 | 0.190 | 0.69 |
| 3 | 0.294 | 0.259 | 0.88 |
| 10 (native) | 0.302 | 0.299 | 0.99 |
| 100 | 0.305 | 0.318 | 1.04 |
| 1000 | 0.305 | 0.320 | 1.05 |

The goal here was the mechanism, not a quantitative match for round-6 Fig. 5; that was left aside as agreed.

## 3. Wave-pinning at strong input: a patch field

**Result.** At strong input the Wave-pinning ring no longer holds one drifting domain. It holds a field of patches whose number grows with circumference while their local behaviour does not change with size. Larger rings still hold heading longer, because heading is the ring-scale (mode-1) component of that field and longer wavelengths decorrelate more slowly.

**Methods.**

- **Input:** the repo's Fig. 4 fluctuating input, reimplemented locally: lognormal (σ = 1) with Gaussian spatial kernel SD 1.5, Ornstein–Uhlenbeck in time (τ = 50) on 2.5-unit spatial knots sampled every 5 time units, mean = f × the reference input.
- **Solver:** semi-implicit spectral (implicit diffusion, explicit reaction), dt = 0.05, dx = 0.05. Run length 32000 with the first 16000 discarded; frames every 5. One control run at dt = 0.025 (C40, f = 0.1) checked the step size.
- **Design:** C20, C40, C60 at f = 0.1 and 0.3 (plus f = 0.03 at C40 as a weak-input reference), two seeds each, at native D\_v = 10 and at D\_v = 1000 (Section 4).
- **Patches:** contiguous regions with u above 0.65 (halfway between the basal level 0.175 and the plateau), at least 10 grid cells, tracked frame to frame by overlap.
- **Mode decorrelation:** for each Fourier mode n of u, the decay rate γ of its autocorrelation at the lag where it is nearest 0.5. Wavenumber q = 2πn/C. Heading is the phase of mode 1.

**Findings.**

- **The field is extensive.** Mean patch count rises roughly in proportion to C: 1.2, 1.9 and 2.8 at f = 0.1, and 1.5, 2.7 and 4.0 at f = 0.3 (C20, C40, C60). The fraction of frames with more than one patch is 17%, 72% and 97% at f = 0.1.
- **Local statistics do not depend on size.** Births per unit perimeter per unit time are 1.3–1.4e-4 at f = 0.1 and 3.5–4.1e-4 at f = 0.3 at every size. Patch velocity (rms about 0.03 and 0.05) and its correlation time (about 25 and 15) are also near-constant.
- **No repulsion between neighbours.** The pair spacing shows only a short-range attraction below about 4–5 units (merging), with no preferred spacing. The field is not a Turing-like array.
- **Mode rates collapse on wavenumber.** γ(q) from all three sizes falls on one curve. A larger ring's mode 1 sits at a smaller q, so it decorrelates more slowly. The local α is the slope of that curve at q1 = 2π/C, which is why α falls toward 1 at strong input rather than staying near 2.
- **Coherence falls with size** (mean |P| 0.42, 0.31, 0.22 at f = 0.1), so heading is defined less sharply in larger rings even while it persists longer.

Patch figures: `8_neighbour_interaction.png`, `9_local_statistics_and_mode_collapse.png` and the kymographs `10_kymographs_sizes_f0.1.png` and `f0.3.png` (the latter embedded as the proposed Fig. 4′ kymograph panel in the narrative).

## 4. Fast-cytosol test at strong input

**Result.** In Wave-pinning, a fast cytosol removes most of the extra persistence that larger rings show at strong input, so cytosolic conduction still sets it there. In Holmes the same change does little, as expected for a model whose native rings sit far below its crossover.

**Methods.** The noisy campaigns of Section 3, repeated with the inactive species made fast. Wave-pinning: D\_v 10 against 1000, run length 32000 (16000 discarded). Holmes: inactive diffusivity 50 against 1000, using the repo's native reaction terms in the same semi-implicit solver, run length 20000 (4000 discarded); total pools were conserved to machine precision. Both at C20, C40, C60 and f = 0.1 and 0.3, two seeds. Local α is the negative slope of log rate against log C over the three sizes, for the mode-1 rate and for the heading rate.

| Model | f | Mode-1 α, native | Mode-1 α, fast | Heading α, native | Heading α, fast |
| --- | --- | --- | --- | --- | --- |
| Wave-pinning | 0.1 | 0.93 | 0.44 | 1.11 | 0.46 |
| Wave-pinning | 0.3 | 0.45 | 0.15 | 0.46 | 0.05 |
| Holmes | 0.1 | 0.82 | 0.77 | 0.95 | 0.85 |
| Holmes | 0.3 | 0.34 | 0.24 | 0.60 | 0.42 |

**Findings.**

- **Wave-pinning:** with a fast cytosol, α roughly halves at f = 0.1 and nearly vanishes at f = 0.3. The patch field itself is unchanged (fraction of frames with several patches 16%, 64%, 92% against 17%, 72%, 97% at f = 0.1); only how fast its large-scale modes decorrelate changes. Small rings are barely affected; the rates rise mainly in C40 and C60.
- **Holmes:** α drops by 0.05–0.2, small next to Wave-pinning and possibly within seed noise (seed spread not yet estimated). Most of the Holmes size effect at strong input therefore does not run through conduction.
- **What remains in Wave-pinning** with a fast cytosol (α ≈ 0.45 at f = 0.1) is close to the Holmes level. One reading, untested: this is a shared part of the strong-input size effect that comes from the patch field's geometry alone.

The Holmes runs completed after the earlier wrap-up and are reported here for the first time; the narrative does not yet use them. Figure: `11_mode_rates_Dv10_vs_1000.png` (Wave-pinning only).

## 5. Reanalyses of round-6 repo data

**Result.** The repo's own Fig. 1–4 data already contain the views the narrative now asks for. No new simulations were needed for these; all come from the round-6 supporting CSVs.

**Turning routes against edge gain (Fig. 2 and 3 data).** Over the matched sizes C30, C37.8 and C47.6 (Spring 2.5, 3.15, 3.97), the two-region half-turn time grows far faster than edge gain falls. Turning by a second region is therefore not set by the linear edge gain.

| Model | Half-turn time fold | Edge-gain fold |
| --- | --- | --- |
| Wave-pinning | 2.18, 5.94 | 0.927, 0.850 |
| Holmes | 3.56, 23.9 | 0.981, 0.960 |
| Spring | 2.76, 9.35 | 0.987, 0.967 |
| Goryachev | 1.95, 4.20 | 0.595, 0.358 |

Each cell gives the middle and the largest size relative to the smallest (source: matched\_B\_susceptibility.csv).

**Edge gain against size (Fig. 3 map amplitudes).** Normalised to the smallest ring, gain falls by a factor of 0.70 in Wave-pinning, 0.91 in Holmes, 0.99 in Spring and 0.09 in Goryachev over a threefold size range. The conduction-law curves with the C\* values of Sections 1–2 pass through the three drift models' points. This is the mockup used for proposed panel 3E′.

**Linearity (Fig. 3 amplitude and duration sweeps).** Gain is independent of pulse duration at fixed dose to within 0.2%. The amplitude at which displacement departs from linear by 10% lies between 5% and 20% of the reference input for Wave-pinning, 2.6–5% for Holmes, 1.3–2.6% for Spring and 2.6–20% for Goryachev (brackets are neighbouring grid points). The mockup for 3D′ shows Wave-pinning displacement against pulse amplitude at C20 and C60.

**Heading memory (Fig. 4D model ACFs).** Heading autocorrelation at lag 2000 (Spring: 200), C20 against C60: Wave-pinning 0.32 against 0.90; Holmes 0.24 against 0.73; Spring 0.09 against 0.72; Goryachev 0.01 against 1.00. Larger rings hold heading longer in every model.

**Mode rates at strong input (Fig. 4 all-input mode rates).** Plotted against q = 2πn/C:

- **Holmes:** collapses well above q ≈ 0.5; local α about 0.8 at f = 0.1 and 0.2 at f = 0.3, in line with Section 4.
- **Spring:** the mode-1 rate is flat with size (α ≈ 0).
- **Goryachev:** collapses almost exactly (at q = 0.31 the rates are 5.5, 5.6 and 5.6e-3 across sizes).

**Goryachev jumps.** The figures used in the narrative (losses follow local input extremes, odds ratio 2.4; 89% of losses preceded by several regions; switch-on input rising with size) come from the coding agent's issue-5 report, not from this session. The jump-share check is still pending.

## 6. Caveats and open items

**Caveats.**

- **α from conduction matches in level only.** The law gives α = 2.36–2.91 over C20–60, close to the measured 2.48–2.88 range, but the measured trend with size is flatter than predicted.
- **Strong-input α rests on three sizes and two seeds.** No seed spread or error bars were computed, so differences of about 0.1–0.2 (for example the Holmes fast-cytosol change) are not resolved.
- **Patch counts depend on the threshold.** Birth and death counts in particular are inflated by patches flickering across u = 0.65. The size-independence of rates holds at the one threshold used; the narrative should not quote these numbers.
- **Spring's law holds only for D\_v ≥ 1.** The \~5% small-ring gain rise with a free cytosol is not explained.
- **Goryachev jump statistics** are the coding agent's (issue-5 report), not reproduced here.

**Open items.**

- [ ] Derive the shape of γ(q) at strong input; α is its slope at q1, but why that slope falls with input is not explained.
- [ ] Goryachev jump-share check (discussion doc, Section 5).
- [ ] Seed spread for the strong-input α values, in Wave-pinning and Holmes.
- [ ] Correct the narrative's 3D′ bullet, which gives 10% departures at 3–25% of the reference input; the brackets in Section 5 are the accurate ones.
- [ ] Decide whether the Holmes fast-cytosol result (Section 4) belongs in the narrative's Fig. 5′.

**Pilot files.** All code, data and figures are in this session's workspace, not in the repo, and that workspace is temporary. Code: `wp/` (Wave-pinning ring, front, noise, patch and mode analysis) and `hs/` (Holmes and Spring gains, Holmes noise). Figures embedded in the narrative:

| Narrative panel | Pilot figure |
| --- | --- |
| 3D′ linearity | `hs/fig_linearity_mock2.png` |
| 3E′ gain against size | `hs/fig_gain_size_mock_A.png` |
| 3F′ fast-cytosol test | `wp/fig3Y_mock.png` |
| Fig. 4′ kymographs | `wp/figs/10_kymographs_sizes_f0.3.png` |
| Fig. 4′ mode rates | `wp/figs/11_mode_rates_Dv10_vs_1000.png` |

The three-model conduction figure (`hs/fig_conduction.png`) was made but is not used in the narrative.
