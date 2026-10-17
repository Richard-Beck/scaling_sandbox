# Manuscript polarity work

- Latest checkpoint scope (user clarification, 17 October 2026): develop three
  figures. Both the all-five established linear 1D figure and the supported
  three-model ring figure must feature gain/recovery measurements AND results
  from a clearly justified standardized long-cue protocol. A separate third
  figure presents replicated noisy-ring experiments measuring polarity-vector
  persistence across sizes, with declared noise statistics and magnitude/local
  pattern diagnostics. Transient assays alone must not dominate the argument.
  See the updated goal at the top of
  `../../current manuscript figs/figure_revision_plan_20261017.md`.

- Main-text figure direction (user clarification, 17 October 2026): Figure 1
  uses the established linear 1D framework and includes all five retained
  models. Figure 2 moves to supported ring models to evaluate separation.
  Keep interval and ring preparations, assays and numerical values distinct;
  ring results must not replace the linear evidence in Figure 1. deBelly's
  corrected established 1D results remain eligible, while its ring/two-arc
  development is a sidecar pending a convincing demonstration. Develop one
  reviewable checkpoint for detailed user feedback. Experimental implications
  concern measurable local activity/protrusion behaviour. Subsequent user
  feedback requests all-five coverage in Figure 1B-C where possible, ring
  counterparts of Figure 1A-C, response heatmaps versus separation and size,
  a cue-associated single-peak example, and a global polarity vector under
  relevant perturbations as an experimental bridge. Define the vector and
  examine magnitude as well as direction; this does not require assuming a
  cell velocity law. See `../../current manuscript figs/figure_revision_plan_20261017.md`.

- Equations in reports and documents require scientific typesetting: centred
  display equations with right-aligned numbers, proper mathematical fractions,
  derivatives, subscripts and superscripts. Keep parameter descriptions in the
  surrounding prose. Render LaTeX/Markdown math in the final document rather
  than exposing its source. The plain-equation preference applies only to CLI
  chat, not to PDF, HTML or manuscript documents.

- Scientific report guidance from the 13 October 2026 annotations:
  use direct, formal descriptions of measured changes. Every mechanistic claim
  must identify the observable, spatial region, cue and time/window, and point
  to displayed or linked evidence. Replace vague speed words with measured
  times. Do not introduce an undefined "memory variable", describe an internal
  Rho maximum as a pole reversal, or treat a spatial average as whole-field
  agreement. Explain exceptions and the full nonmonotone curve before focusing
  on a restricted length interval. A selected successful pulse is not evidence
  that withdrawal has no effect, and continuous state at withdrawal does not
  mean unchanged subsequent dynamics. Avoid unsupported strawman comparisons.
- For the deBelly illustrative sequence, use the relative-20% right patch from
  each length's exact corrected basal checkpoint; explicitly label other cue
  geometries in comparison figures. Use corrected independent endpoints only.
  Keep cortical endpoint commands distinct from membrane endpoint velocities.
- MCA budgets must name the pool. A bound-plus-unbound budget cannot establish
  which terms control bound MCA: physical exchange cancels in the continuous
  total-pool equation, while the native sequential updates leave a numerical
  exchange residual. Show that residual separately when checking total mass.
  When discussing local bound MCA, separate profile advection, flow-gradient
  dilution/concentration, bound diffusion, net binding and uniform replacement.
  Quantify their magnitudes before dropping terms. Endpoint outflux is distinct
  from local flux divergence. Preserve the native time staggering and document
  the exact finite-volume donor/face convention used to split advection.
- Mechanical source superposition allocates the flow at a fixed coupled state;
  it is not an evolved pathway ablation or a causal partition of reversal time.
  Centroid redistribution differs from molecular travel and local front
  depletion. Keep finite-amplitude polarized-branch existence distinct from
  instability of the homogeneous state, and do not infer a universal withdrawal
  threshold from selected successful trajectories.
- Prefer fewer readable panels to dense pages. Give MCA profiles adequate
  space, expose boundary-cell behavior instead of hiding it with axis limits,
  and keep diagnostic details in linked methods. Figure prose must state the
  observed result and its mechanistic implication, rather than narrating panel
  locations. Review new sections against the full annotation checklist under
  `figures/length_mechanisms_20261008/revision_20261013/editorial_review/`.

- Goryachev mechanism direction (2026-10-08): use the single report under
  `figures/Goryachev_length_mechanism/`. The two supported effects are lower
  inactive concentration reducing initial equal-countercue activation, and a
  stronger immediate incumbent restoring response per unit damage in longer
  cells. Preserve the distinction between concentration and total amount, and
  between immediate rebound and slower late recovery. These coupled effects
  are not an exhaustive or additive attribution of reversal duration. Previous
  exploratory Goryachev presentations are archived; compact report inputs and
  the cleanup inventory are in `data/runs/goryachev_length_mechanism_20261008/`.

- Alternative perturbation (2026-10-08): whole inactive-profile bolus removal
  and exact pointwise return, with native kinetics unchanged. Both Goryachev
  and Wave-pinning figures are under their `inactive_pool_perturbation/`
  subfolders. Primary hold time is 0.1; sensitivity is 1 native unit, with
  reductions 0.1%, 0.5%, 1%. Data and methods are in
  `data/runs/inactive_pool_figure_exports_20261008/`; all 288 archives audit.
  Keep these distinct from the earlier local activation-suppression figures
  and their different recovery tails. Regenerate using
  `code/mechanisms/inactive_pool_mechanism_figures.py render`.

- Current user direction (2026-10-06): read THREE_TESTS.md. All three response
  tests share the exact relative-20% nominal basal checkpoint at each length:
  fixed support 3 (Spring 0.2), relative 20%, finite constant-physical-slope ramp.
  Goryachev's audited K=64 patch multiplier is fixed across all lengths and both
  patch tests. This explicit direction supersedes the old matched-peak/matched-
  preparation-width rule below. Mechanism analysis covers all five × three cells.

- Read README.md, PROTOCOL.md, ASSUMPTIONS.md and the fresh campaign status before
  changing inputs or interpreting the comparison.
- Active reversal models are Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring.
  User excluded Otsuji and LEGI. Keep their corrected establishment evidence and
  old field-explorer diagnostics, but do not include them in reversal fits.
- Start every right challenge from its exact qualifying corrected basal-release
  checkpoint at nominal amplitude. Preserve preparation/countercue mirror,
  native background, width and peak; no reflection or substituted old preparation.
- Original retained evidence is immutable. Write fresh runs under data/runs;
  frozen executable inputs cannot change while their producer runs.
- Runtime dependencies belong under code/. Historical sources/records retained
  for diagnostics are provenance, not active manuscript pipeline inputs.
- Generate final comparison figures only after complete execution and independent
  checkpoint/gate/selection audits. Missing work is an error; it is not censoring.
- Scalar strength replacement and native field stability are separate measures.
  Preserve unavailable statuses and field-stability flags in figures and fits.
- Use the existing scaling_sandbox environment; full integration runs use SLURM,
  a suitable QoS and no unjustified array concurrency throttle.
- No timed HTML self-refresh. Show data-generation time and producer state.
