# Current manuscript figure copies

**Foundation copy audit, 17 October 2026:** `main_polarity.png` now matches the
latest [canonical five-model comparison](../polarity1d/manuscript/figures/main_polarity.png),
generated on 10 October. Panel C reports Spearman rank correlations using all
audited scalar successes; it supersedes the older OLS/no-fit working copy.
`supplement_width_response.png` is also byte-identical to its canonical source.
Both retain the corrected independent-endpoint deBelly replacement, the separate
native-field stability diagnostic, and all three shared-checkpoint finite-cue
assays. See the [canonical legends](../polarity1d/manuscript/figures/CAPTIONS.md)
for the inclusion rule and success-conditional interpretation.

This is the working figure collection for developing a manuscript; there is no
completed manuscript. The biological motivation is **increased migration
persistence in WGD-associated enlarged cells**. The 1D length experiments examine
size-dependent polarity organization and competition under unchanged local
kinetic parameters. Their experimental endpoint is measurable local activity
and protrusion behaviour; predicting cell movement is outside the present
figure scope.

Start with [main_03_migration_OU.png](main_03_migration_OU.png). This copy is
byte-identical to the
[3 October refined OU figure](../wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU.png).
Its [caption](../wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU_caption.md),
[observations](../wgd_langevin/refined_ou_20261003/figures/main_03_observation_folds.csv),
[persistence estimates](../wgd_langevin/refined_ou_20261003/figures/figure3_persistence.csv),
and [results and limitations](../wgd_langevin/refined_ou_20261003/results/RESULTS.md)
provide the empirical and modelling context:

- Xuan: the polyploid group has approximately 0.64 times the sampled migration
  speed, but 2.08 times the population MSD at the last shared lag (about 10.95 h).
  The full digitized population MSD curves cross at approximately 2.32 h.
- Mallin: the PACC group has approximately 0.52 times the path length and 3.12
  times the directness (net displacement divided by path length).
- Hosny baseline: the PACC group has approximately 1.77 times the median path
  length and 1.67 times the median directness. Signed migration index is a
  separate directional quantity.
- The OU fits infer longer velocity memory in all five fitted comparisons.
  Persistence folds are conditional fitted parameters, not directly measured
  autocorrelation times. The Jemaà studies constrain speed alone, so their
  persistence estimates depend on the hierarchy and priors. The displayed
  example tracks are simulated.

These observations motivate explaining sustained displacement and directional
persistence separately from speed. The experiments span stable tetraploids,
giant polyploids and treatment-induced PACCs; the collection does not establish
that size alone causes their migration differences.

For the present mechanism-figure discussion, use the
[current length-mechanism PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf)
and its [reader](../polarity1d/manuscript/figures/length_mechanisms_20261008/index.html).
The most recently completed review package is the
[updated checkpoint](updated_checkpoint_20261017.md). **Figure 1 uses all five
established linear 1D models; Figure 2 uses supported rings to evaluate
separation.** Figure 1 now includes every retained linear input-gain value and
a new local-core recovery pilot at three sizes per model. Figure 2 spans two
continued review sheets: linear/ring correspondence and shared-scale gain maps,
then complete outcome maps, three highlighted native histories and their global
activity vectors. The matching new ring-recovery measurements remain explicitly
pending. The [fresh assembled-figure critique](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/updated_checkpoint_20261017/supporting/metadata/FRESH_CHECKPOINT_CRITIQUE.md)
finds no essential handoff blocker. This is a review checkpoint, not a finished
publication layout.

The [revision record](figure_revision_plan_20261017.md) now defines the active
three-figure expansion: both linear and ring figures feature gain/recovery and
justified long-cue responses, followed by a separate replicated noisy-ring
polarity-persistence figure. The user has resumed this work; the
[completion requirements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/VISION_CHECKPOINT_REQUIREMENTS.md)
track the evidence needed for the next checkpoint. The
[complete all-five linear evidence](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_linear_20261017/RESULTS.md)
and [complete supported-ring evidence](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_ring_20261017/RESULTS.md)
have passed independent actual-result reviews; the replicated noise ensemble
and final package remain active. The
[preceding linear-then-ring checkpoint](linear_then_ring_checkpoint_20261017.md)
is preserved; its selected two-model functional coverage prompted this update.

The [preceding compact two-figure candidate](compact_integrated_candidate_20261017.md)
is retained as history; it placed ring response curves in Figure 1.
It retains the all-five incumbent-organization panel, uses shared-scale size
curves for all three chemical models, and links highlighted endpoint-map cells
to their exact merger/two-site kymographs. The figure explanations and fresh
critique distinguish conditional repair from actual incumbent maintenance.
The [earlier integrated framework](main_text_figure_candidates_20261017.md)
and [detailed panel discussion](integrated_framework_candidate_20261017.md)
retain the preceding six-panel draft and its raw evidence.
The [archived supporting-component overview](supporting_components_before_integration_20261017.md)
preserves the prior chemical/spatial comparison, supplied-PDF evidence table and
other-model notes. The [earlier discussion](length_mechanism_figure_drafts.md)
retains the preceding drafts and development of the circular-domain test.
The older mechanism PDF copies in this folder do not supersede that report.

The [systematic chemical-ring follow-up](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/RESULTS.md)
now compares Wave-pinning, Goryachev and Holmes using shared native-ring
preparations for their input response, conditional repair and spatial histories.
Its [fresh critique](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/supporting/metadata/RESULT_CRITIQUE.md)
checks the distinction between a restoring contribution and actual incumbent
retention, including a Goryachev example with transient challenger activation
and restoration of the original incumbent. Two-site retention also occurs in
larger Wave-pinning rings. These are subsequent development results beyond the
original PDF; their numerical values come from ring rather than interval states.

**deBelly scope (user direction, 17 October):** established corrected 1D
results may participate in main-text panels, including incumbent organization.
Its ring and two-arc extensions are a separate sidecar pending a convincing
ring demonstration. Main ring comparisons proceed with Wave-pinning,
Goryachev and Holmes; deBelly development does not hold up these figures.

The sidecar [mixed-boundary deBelly two-arc pilot](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_shared_chemical_arcs_20261017/RESULTS.md)
and its [fresh critique](reviews/debelly_mixed_arc_result_r1.md) document the
subsequent shared-chemistry baseline: localized chemical repolarization works
on unequal arcs, while mechanical flanks remain independent. All tested cue
positions select the challenger at this one size and dose, so the pilot supplies
a workable construction rather than a demonstrated distance-dependent fate.
The [same-state duration contrast](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_pulse_duration_contrast_20261017/RESULTS.md)
shows challenger activation followed by incumbent return after a 50-s cue,
versus repolarization after a 1000-s cue at the same location and strength.
Both use the prescribed permanent-cut baseline. In the subsequent
[Rac-triggered-gap sidecar](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_rac_continuity_20261017/RESULTS.md),
each connected above-threshold Rac region without a break creates one at its
maximum. Breaks stay at their birth
locations and close when their own local Rac falls below threshold. Merging
active regions retain existing breaks until that local closure condition holds.
The bounded preparation polarizes but remains time-dependent: none of 30
release windows qualifies. Its [fresh review](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_rac_continuity_20261017/FINAL_REPORT_REVIEW.md)
supports that diagnostic result; no new-rule duration challenge was run.
The subsequent user-authorized
[moving-break trial](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_moving_breaks_20261017/RESULTS.md)
also retains localized Rac but passes none of 30 actual or field-only
stationarity windows. Its report and fresh critique expose the mechanical
controller mismatch caused by relocating a break. This bounded diagnostic is
complete; it supplies no qualified ring checkpoint or subsequent challenge.
The [experimental-readout discussion](migration_consequences_bridge_20261017.md)
connects incumbent organization, challenger establishment and later competition
to possible activity measurements in migrating cells. Cell turns are an optional
observational association, rather than a required model-to-movement bridge.

**deBelly warning (2026-10-08):** the old two-end extension with `center_boundary_motion=True` is flawed and deprecated. Do not use its simulations for scientific inference. `model_debelly.png`, the old deBelly portions of the mechanism PDFs, and the deBelly row in the initial-profile report retain that historical method. They are not corrected-model results.

The corrected independent-endpoint replacement campaign is in ../polarity1d/manuscript/data/runs/debelly_corrected_three_test_20261008/. `main_polarity.png` and `supplement_width_response.png` now contain the completed independently audited corrected deBelly replacement; prior copies are archived in archive/centered_debelly_20261008/. Corrected mechanism figures are grouped in ../polarity1d/manuscript/figures/deBelly_mechanism/. Read ../polarity1d/manuscript/DEBELLY_BOUNDARIES.md for the boundary method and executable entry point.
