# Figure 2: evidence for discussion after author feedback

19 October 2026. **Planning note, not a replacement figure plan.** The reviewed Figure 2 is a rejected presentation prototype; its artwork and content are preserved. No new simulation, layout or rendering accompanies this note. The next draft must restore incumbent peak/size organization, remove the duplicate reversal-duration panel and activation/fate scatter, and reconsider the selected native histories. The author has not selected the challenge, optical actuator or replacement presentation.

## Exact retained organization evidence

The [current manuscript copy](<../../../../../../../../current manuscript figs/initial_polarity_length.png>) is the original four-row, five-model incumbent-size figure: active amount, width, mean/peak concentration and opposite-end inactive concentration. Its PNG is byte-identical to `figures/initial_polarity_length/initial_polarity_length.png` and the report copy; SHA-256 starts `b94d05c1caab`. The [later editorial arrangement](../../../../revision_20261010/editorial/assets/initial_polarity_length.png) contains the same 97-checkpoint measurements (PNG SHA-256 starts `873f1302020b`). Both are retained alternatives, not new artwork.

The [source table](../../../../../../data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv) contains all 97 eligible corrected basal states: WP 20, Goryachev 21, deBelly 12, Holmes 23 and Spring 21. Width is the boundary-connected half-contrast region; amount and mean use that region, with field-minimum background. Amount is a one-dimensional concentration integral, not a calibrated molecule count. Spring uses physical width/amount against resting length; the original deBelly organization figure uses nominal chemical-solver coordinates.

The strongest passive observation is **more active material accommodated by different spatial changes**. At L19→39, WP excess amount increases 2.124×, width 2.042× and peak 1.000×; Goryachev amount increases 2.152×, width 1.119× and peak 1.903×. The other three also broaden over their selected pairs. This motivates a challenge prediction, but morphology alone does not establish its cause.

This was already requested: the [17 October composition note](../../../vision_aligned_checkpoint_20261017/supporting/metadata/COMPOSITION_PLAN.md) says “Keep the all-five incumbent organization summary” and “do not substitute threshold-duration curves” for the long-cue observations. The [revision program](../../../../../../../../docs/tasks/29_10_19_revision_program.md) likewise lists incumbent width/concentration/amount among evidence to preserve. The [earlier integrated narrative](../../../../revision_20261013/chapters.json) explicitly used the editorial organization artwork as Figure 2.

## Three possible observable signatures

These are options for choosing an experiment-facing argument. The first two have the clearest existing cross-size quantitative contrasts. The third has the clearest retained spatial illustration. None is a stand-alone molecular classifier.

### 1. Early challenger activity after a characterized local activation pulse

**Perturbation/readout.** Apply the same physical activation footprint and characterized native input strength within a size series. Measure the stimulated region's activity increment over a declared early interval, alongside incumbent peak and width. An experimental activity reporter would supply the activity change; light intensity alone does not establish the biochemical dose.

**Existing evidence.** The saved native fixed-support trajectories at L19/39 use support 3 and native peak 0.1 for WP, 6.4 for Goryachev. During the first saved 0.1-unit interval, Goryachev challenger peak increments are 0.065977 and 0.043054 across sizes, a 0.65256 fold response. WP's baseline-subtracted increment is 0.017520→0.017486, a 0.99805 fold response. Exact snapshots and later reads are extracted in [observed_native_responses.csv](supporting/data/observed_native_responses.csv), with native archive identities. The extracted peak is measured in the fixed mirror of each original basal half-contrast core. These are finite-time field responses, not analytic source coefficients or an inferred infinitesimal slope.

**Mechanism distinction.** Both highlighted native channels add a local coefficient multiplying inactive substrate: extra active production is `stimulus × inactive`, with equal loss from inactive. Goryachev's taller incumbent coincides with less available inactive substrate; WP retains almost constant inactive concentration. The declining versus nearly invariant response is therefore a measurable consequence consistent with substrate-limited entry versus domain expansion under these particular channels.

**Ambiguity/limit.** Local feedback, deactivation and diffusion already contribute during 0.1 units. Substrate limitation is not uniquely established by a smaller activity increment; weaker actuator coupling or stronger inhibition could mimic it. Recruitment of a GEF, direct addition of active material and optogenetic inhibition need not implement the same channel. In particular, activity forcing independent of inactive substrate would not have this native instantaneous substrate factor. No matched experimental actuator has been established. Compare size folds within a model; the two native amplitudes and clocks are not calibrated across models. The narrow Goryachev39 incumbent has only four sampled centres in its initial core; these records do not establish spatial convergence.

### 2. Restoration after a small local incumbent activity loss

**Perturbation/readout.** Briefly deactivate a small measured fraction of incumbent activity, with an otherwise identical unperturbed control. Measure the induced deficit and its recovery in the original spatial core, retaining peak and edge profiles. This requires a perturbation of activity; bleaching the reporter alone would not implement it.

**Existing evidence.** The completed common assay removes a smooth 0.5% basal-excess direction inside the original core and transfers active material pointwise to inactive in WP/Goryachev/Holmes/Spring. deBelly changes Rac alone because it has no explicit inactive Rac pool. The [all 97 series](../supporting/data/all97_gain_recovery.csv) gives Goryachev restored fraction 0.05404→0.18588 at L19→39 (3.439×), versus WP 0.13391→0.13735 (1.026×). Their separately qualified native clocks are 0.0056234 and 0.177828, respectively. The [assay definition](../../../all_five_linear_function_20261017/DESIGN.md) preserves the frozen core/taper, native fields, paired sham and amplitude sensitivities.

**Mechanism distinction.** In combination with organization and the excitation response, stronger recovery with size supports a second obstacle to replacing the tall Goryachev incumbent: after a specified small loss, the coupled system restores more of that loss over its fixed early window. Nearly invariant WP restoration indicates that this particular increasing-restoration signature is absent there.

**Ambiguity/limit.** This assay combines native reaction and transport; it does not isolate an autocatalytic constant. Positive restoration can coexist with an incumbent declining in both sham and probe, and does not prove resistance to a finite challenger. Goryachev substrate depletion and stronger restoration covary, so their contributions to reversal time are not additively identified. The older whole-inactive-pool removal/return assay is a different disturbance and must remain separate. A physically matched optical deactivation needs its actual affected pool, spatial profile and amplitude characterized before these model results become experimental predictions.

### 3. Domain advance and withdrawal survival, rather than cue-end peak alone

**Perturbation/readout.** Use two specified local activation durations, then switch the cue off. Measure the full spatial activity profile: challenger height, expansion beyond the illuminated patch, incumbent edge retreat, and post-withdrawal survival. This directly links the baseline width/height distinction to a visible response rather than an endpoint scatter.

**Existing evidence.** The retained [WP L15 short/long pulse illustration](../../../../revision_20261010/editorial/assets/wave_pinning_withdrawal.png) uses the exact same basal checkpoint, support 3, peak 0.1 and durations 20/60, followed by 10000 native units without cue. Both create a high challenger peak (1.3295/1.4725 at withdrawal). The short pulse collapses and returns to the old state; the longer pulse leaves a wider new domain and ends at the reflected state. Rechecking both native chemical fields gives maximum fixed-scale endpoint errors below 4.6×10⁻¹³. The old artwork already uses identical axes and shows interface positions clearly.

**Mechanism distinction.** A fixed high-active plateau with advancing/retreating edges predicts that spatial domain growth, not cue-end peak height alone, separates transient excitation from retained polarity. This is a concrete observable expected for the retained WP domain-competition explanation.

**Ambiguity/limit.** The two pulses change dose as well as duration and are one size. Dose-dependent survival also occurs in other mechanisms; this pair is an illustration, not cross-model discrimination by itself. Corrected deBelly L39 also returns after support 3/T500 but retains challenged orientation after T1000/2000, despite strong cue-on activity in all three. That similarity is an explicit limit: return after withdrawal alone cannot diagnose slow mechanical rewriting. A mechanical claim would require a concurrent mechanical/attachment observable and an appropriate intervention, not relabelling activity fate as mechanics.

## Calculated gain and the visually weak Goryachev challenger

The existing gain panel calculates how the native input changes the **initial production term at the unchanged basal state**, averaged with the cue footprint. It is not a measured reporter slope, accumulated activation, challenger height or establishment probability. Option 1 instead uses actual saved activity after cue onset; the first interval includes all native processes. No optical recruitment or reporter kinetics have been simulated.

For the existing E case (Goryachev L39, support 3, peak 6.4, T1000), activity is positive: at withdrawal the challenger peak is 11.7598 and old peak 100.0520, against initial incumbent peak 117.3903. A full-field zero-based colour range near 117 makes the challenger hard to see. It is not absent.

The existing archive already contains exact spatial fields and dense early times, so no new simulation is needed to show this honestly. A possible display would pair the full profile with a clearly labelled stimulated-region enlargement or native-unit local activity curve, using the same local response scale across the compared sizes. Baseline-subtracted activity is another declared readout; a sham is needed for cue-attributable change if basal drift is material. Independent per-history peak normalization would conceal the amplitude contrast and should not imply equal responses. Experimental detectability remains unknown without reporter noise, dynamic range, time response and actuator calibration. The investigator must also choose whether the useful question is small early activity, late challenger height, spatial spread or post-withdrawal survival; selecting one is not authorized by this note.

## Inference boundary and provenance

Two inspected primary studies provide concrete experimental context for GEF
recruitment and activity-reporter characterization. Their source identities,
observations and limits are in
[experimental_precedents.json](supporting/metadata/experimental_precedents.json).
They inform the discussion of readouts; they do not select a tool or establish
a quantitative optical mapping for the planned cell experiments.

The existing evidence can connect different incumbent organizations to different conditional responses. It does **not** causally allocate the shared size-dependent reversal trend to uniquely identified mechanisms. The WP 4× inactive-diffusion control improves retained orientation at one fixed dose, acts during cue plus recovery, and still fails the stronger full-polarity-strength criterion in both conditions; it cannot close that gap by being called a general replacement rescue.

[evidence.json](supporting/metadata/evidence.json) records exact source paths/hashes, recomputed field samples, the retained artwork identities and unchanged prototype hashes. This note leaves challenge choice, optical-channel mapping, main-figure emphasis and any new measurements for discussion with the author.
