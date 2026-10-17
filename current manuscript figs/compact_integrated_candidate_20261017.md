# Compact candidate: incumbent growth, conditional response and spatial competition

This new draft integrates the earlier framework and the systematic chemical-ring comparison into **two figures with different jobs**. Existing reviewed figures are unchanged. It uses retained data only; no new simulations, metrics or parameter tuning were performed.

## Figure 1: more incumbent material does not imply the same functional change

![Incumbent organization and honest shared-scale ring-native functional curves](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/compact_incumbent_and_function.png)

**Message.** Enlargement increases incumbent active material in all five retained preparations. The organization differs: Goryachev predominantly raises concentration, whereas the other models predominantly broaden the incumbent. In separately prepared native rings, Goryachev also shows substantially lower initial cue gain and stronger conditional repair with size; Wave-pinning and Holmes have nearly unchanged gain and much smaller repair changes. These are model-specific consequences, not a universal stability result.

**Why this presentation.** Panel A preserves the scientific content of the previously liked organization panel. B–C compare all three chemical models as size curves, using one common axis per observable and each model's own C19 reference. Initial gain varies by less than 0.09% over eligible placements within any size, so heatmaps add little. The repair axis explicitly shows logarithmic folds; there are no narrow model-specific axes that magnify Holmes's small change. Native values are retained in the [curve table](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/data/ring_native_size_curves.csv).

| Native-ring C19 → C39 | Initial input gain fold | Basal conditional repair fold |
|---|---:|---:|
| Wave-pinning | 0.99993 | 1.1586 |
| Goryachev | 0.52752 | 4.2672 |
| Holmes 3 | 0.99966 | 1.0147 |

**Caption.** A, excess-active amount, half-contrast core width and absolute peak concentration at the original near-doubling size pairs. These are separate observables, not an exact amount = width × peak identity. Spring uses resting length L₀ 1→2 and measured physical domains. B, immediate extra activation per unit raised-cosine input in the antipodal width-3 patch: exactly cue-weighted inactive substrate for Wave-pinning/Goryachev, or inactive Rac divided by 7.5 for Holmes. This is incremental input gain, not total activation propensity. C, immediate relative deficit damping after the specified inactive-pool hold/return assay, measured once per basal model/size. **Positive conditional repair does not imply net core gain or incumbent survival:** all 652 unperturbed challenged cores lose active material at t10. Open rings identify C19 and C39. B uses a common linear axis; C uses a common log₂ fold axis. Native rates are not calibrated across models. Holmes C80 is unavailable after failed single-incumbent preparation. A's interval checkpoints and B–C's ring checkpoints are different preparations; ring core amounts also rise approximately 2.1–2.2-fold over C19→39.

## Figure 2: separation modifies the later activity pattern

![Complete three-model endpoint maps linked to two exact same-state Holmes kymographs](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/compact_size_location_and_histories.png)

**Message.** Almost placement-independent initial gain does not imply placement-independent competition. The final activity pattern changes with both size and placement. Holmes provides a clear fixed-size location comparison: moving the challenger changes merger into two retained sites, despite the same complete starting state. The maps give the all-model context for choosing this illustration.

**Why this presentation.** The only main-figure heatmaps encode spatial outcomes. Circumference uses its actual numeric spacing, and categorical colors have the same meanings in all three models. Letters D/E mark the exact map cells illustrated below, rather than leaving the kymographs as unexplained examples. The two Holmes histories share one raw concentration scale and linear time. The coordinate cut is moved away from the incumbent: 0 marks its original center, and −10/30 are the same periodic point. The initially single incumbent is therefore not split across the plotted edges.

**Caption.** A–C, spatial association at 14,000 native units, after a width-3 raised-cosine cue lasting 1000 units and 13,000 cue-free units. Cue peaks remain model-specific: Wave-pinning 0.1, Goryachev 6.4, Holmes 1. All 867 planned cells are represented: 652 complete challenges, 198 cue/core overlap exclusions, and 17 unavailable Holmes C80 preparations. Colors describe endpoint geography; a sole old-associated peak does not establish tracked incumbent survival. D–E, exact native h = 0.05 Holmes C40 histories for the highlighted d12 and d16 cells, shown through 4000 units. Their full four-field basal states are identical. The cue position changes routes 12/28 into 16/24; both routes change together. White/cyan dotted lines mark original/cue centers; the dashed horizontal line marks withdrawal. Merger versus distinct-site retention agrees with the full 14,000-unit endpoints. These are finite-window outcomes, not stable-branch classifications.

## What the repair readout means, and what stays in support

The probe removes 0.5% of the whole current inactive-primary profile for 0.1 native unit, then returns that exact profile pointwise. All fields evolve with native kinetics. The matched **sham** has no inactive-pool withdrawal and otherwise experiences the same ongoing cue. Let D be sham minus probe active material inside the **fixed original basal half-contrast core**, immediately after return. The conditional coefficient is

$$
k_{\mathrm{cond}}
=\frac{\dot M_{\mathrm{probe}}-\dot M_{\mathrm{sham}}}{D}
=-\frac{\dot D}{D}.
\tag{1}
$$

A positive value means that this particular deficit is shrinking relative to the sham at that instant. It does not mean the challenged core is gaining material or that a distinct incumbent survives. The **sham net rate** is the actual signed active-material rate in that fixed region, including reaction and diffusion; it is not a second repair coefficient.

The previous t10 heatmaps apply the **same assay** to a state that has experienced ten cue-on units. Basal and t10 measurements differ in the state being probed, not in the definition of repair. They do not supply independent measures of generic stability. Their distance effects are also much less informative than the spatial outcomes here. We therefore retain these [challenged-state probes and actual regional-rate budgets](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/RESULTS.md) as diagnostics. Notably, all 652 unperturbed early cores lose material despite positive conditional coefficients. The region is open: diffusion can move active material out of it without destroying the site.

The additional [Goryachev C40/C80 histories](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/ring_functional_Goryachev_C40_C80_antipodal_history.png) remain linked support. They directly establish replacement versus incumbent retention with transient challenger activation in those particular runs; the endpoint map alone cannot establish that history. They are omitted from the main layout to keep one clear same-state separation comparison. No result here isolates growing incumbent size from the accompanying changes in material, local state and routes across circumference.

## Relation to the supplied PDF and experimental readouts

Figure 1A condenses supplied-PDF Figure 2. The original PDF already supports substrate-limited input and conditional recapture, and **its Figure 6 already shows Wave-pinning challenger activation/collapse, old-front return and retained reversal**. The new contribution is the systematic same-state, three-model ring size/position comparison, including merger and two-site outcomes. Ring-native numbers must not be presented as the original interval measurements.

Evidence for greater migration persistence in WGD-associated enlarged cells motivates this work; the figures aim to distinguish **experimentally measurable polarity behaviours in migrating cells**. Relevant readouts are incumbent extent, peak and integrated reporter signal; early challenger activation; local incumbent loss and rebound; times to extinction, replacement or merger; and two-site lifetimes, stratified by cell size and boundary separation. Integrated reporter signal requires appropriate calibration before it can represent molecular amount. Actual loss/rebound during competition must be measured independently of the conditional hold/return coefficient. For spontaneous events, compare rates per cell and per unit accessible boundary length and time using common detection criteria: different upstream fluctuations can change event rates without changing input gain. Cell turns can be an optional observational association; these figures neither predict trajectories nor require a force-to-movement mapping, and the proposed activity measurements would not by themselves explain WGD-associated persistence. The [experimental bridge note](migration_consequences_bridge_20261017.md) develops these distinctions.

## Sources and review

- [Earlier integrated draft](integrated_framework_candidate_20261017.md), [complete ring evidence and critique](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/RESULTS.md).
- [Original all-five values](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/data/all_five_original_interval_folds.csv), [complete endpoint table](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/data/complete_endpoint_map.csv), [exact selected native arrays](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/data/Holmes_C40_d12_exact_native.npz).
- [Verification and scope](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/metadata/verification.json), [input/output hashes](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/metadata/manifest.json), [render-only generator](../polarity1d/manuscript/code/mechanisms/compact_integrated_candidate_20261017.py). Vector PDF/SVG exports accompany each PNG under `supporting/vectors/`.

The [fresh-context critique](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017/supporting/metadata/FRESH_CRITIQUE.md) supports this two-figure structure as a compact discussion candidate. Its essential requests are resolved: the input-gain definition and positive-repair/net-maintenance distinction now travel with Figure 1, and the label/footer collision is removed. The reviewer independently checked the shared complete Holmes basal, highlighted cells and unclipped common concentration scale. The remaining limitation is scientific: these comparative observations do not establish an isolated causal chain from incumbent growth to migration persistence. This package is a reviewable candidate, not an accepted main-text figure.
