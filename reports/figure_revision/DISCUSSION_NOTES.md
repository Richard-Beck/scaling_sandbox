# Discussion notes on the current revised figures

Recorded 2026-09-29. These notes concern `figure_revision/figures` and exclude
the initial revision request in this discussion. They distinguish outstanding
feedback from changes already completed. Recording these notes does not initiate
figure edits or report regeneration.

Panel references below refer to the current figures unless explicitly marked
as proposed. Proposed S1 will explain periodic cue following using current
S3A–B; current S1 becomes a substantially redesigned S2. The spring-model
supplement and any remaining current S3 material will need renumbering or
consolidation once the final layout is settled.

## Figure 1

- Reintroduce periodically shifting cue following in the main figure. The
  intended four-test protocol is peak competition, gradient reversal, finite
  opposing-cue reversal with stimulus-free recovery, and periodic cue following.
- Exclude Otsuji and Goryachev from both direct cue-reversal analyses on
  analytical/prior-art grounds; record the specific arguments and citations
  when implementing this revision. Periodic cue following provides a fairer
  test of their response to changing cue direction. Keep cue following
  conceptually distinct from lasting reversal after cue withdrawal.
- Current 1F should include only wave-pinning and LEGI. Current 1H should
  include only wave-pinning and follow the presentation style of current 1B
  and 1E.
- Keep the scaling-control overlays out of the main figure pending further
  thought about their biological meaning. Assay-specific abundance scaling is
  not readily interpretable as a biological scaling rule.
- Current H's normalized first moment measures polarity. At L = 15, only wave-pinning
  shows lasting reversal under the tested protocol; this is not a universal
  restriction across models or lengths.
- The current example traces use different opposing-cue durations: approximately
  14.4 s for wave-pinning and 5,120 s for Otsuji and Goryachev. This matters when
  interpreting successful replacement versus failed challenges in the existing
  figure; the revised H will retain only wave-pinning.

## Proposed Supplementary Figure S1: periodic cue explanation

- Move the periodic shifting-cue explanation in current S3A–B here to support
  its return to Figure 1. Coordinate the main-figure results and supplementary
  explanation rather than retaining the old framing of periodic tracking as
  a demoted assay.

## Proposed Supplementary Figure S2: cross-assay scaling (current S1)

- Substantially re-imagine current S1. Re-showing the prior-art Chiou
  saturation result is no longer the objective; cite it as background where
  needed.
- Test whether mean-saturation rescaling is a defensible scaling rule. The
  concern is that the existing inactive-mean matching is assay-specific:
  initial abundance is adjusted to match the reference inactive mean at an
  assay-specific time. A separate match for each assay does not establish a
  common biological scaling rule.
- Make a good-faith attempt to find one mean-preserving abundance-versus-size
  rule for wave-pinning that works across all four Figure 1 assays. Choose
  wave-pinning because it works in all of these assays. Define the mean being
  preserved, the reference state, and the matching criterion explicitly;
  use the same rule across assays without assay-specific retuning.
- The working hypothesis is that a satisfactory shared rule will fail, but
  this must remain an empirical question. Report a successful shared rule if
  found; if none is found, describe the search and its limits rather than
  claiming impossibility.
- Assess both whether the intended mean is preserved and whether eligible
  polarized states and assay responses survive across lengths. Retain
  assay-specific matches as diagnostics of why a shared rule succeeds or
  fails, not as evidence that one has already been found.
- Prefer interpretable active and inactive spatial profiles at representative
  lengths (approximately L = 15 and L = 30) to the current cluttered D–F and
  normalized diagnostics in G–I. The final panel layout should serve the
  cross-assay question rather than repeat the old saturation demonstration.
- Clearly distinguish fixed-concentration scaling, inactive-mean matching,
  and physical peak preservation wherever compared. Physical peak
  preservation retains reference peak profiles while filling added domain
  with background concentrations; matching inactive mean alone does not
  preserve peak content.

## Additional models: common four-test screen

- Trial a two-species polarity model, with Model 3 of
  [Holmes et al. (2012), Modelling Cell Polarization Driven by Synthetic
  Spatially Graded Rac Activation](https://doi.org/10.1371/journal.pcbi.1002366)
  as the suggested example. This is the mutually inhibitory Rac–Rho model.
  Here “two species” means two GTPases, each with active/inactive pools;
  specify whether the optional phosphoinositide layer is included.
- Screen this model, the simple spring model, and the deBelly model with the
  exact same four-test protocol used in the revised Figure 1: peak
  competition, gradient reversal, finite opposing-cue reversal with
  stimulus-free recovery, and periodic cue following.
- Finalize the common protocol before screening, including conditioning,
  cue schedules, size sampling, eligibility criteria, readouts, and success
  criteria. Document any required mapping of model units or cue inputs;
  avoid model-specific assay substitutions or success definitions. Record
  ineligible states and failures as outcomes of the shared screen.

## Figure 2

- No additional revision was identified in this discussion; most of the
  previously requested restructuring is already present.
- Tension-feedback ablation is already shown as Rho feedback 0×. It removes
  tension-dependent Rho activation while retaining the other mechanical
  processes and MCA-dependent Rac feedback.

## Current Supplementary Figure S2: spring model (final number pending)

- Peak competition and kymographs already exist.
- Address E's apparent orientation flip. The underlying example never changes
  polarity sign: both stimulated peaks decay toward a nearly uniform state,
  after which left polarity develops. The visualization makes the early bias
  difficult to read.
- Make the visualization and interpretation clearer. The example is not a
  clean demonstration of direct peak-to-peak competition because it nearly
  depolarizes before selecting a winner.
- Add reference-length “moving kymographs” for both polarity reversal and
  peak competition, potentially replacing current G–H.
- Incorporate the common four-test screen above into the revised plan.

## Figure 3

- Consider a persistence-versus-noise cost/loss surface showing relative loss
  across parameter combinations, with the group-mean baseline and WGD
  parameter values marked.
- Add representative WGD− and WGD+ track plots to visualize their differences.
- Use annotations to show how tracks produce speed, MSD, path length and
  directness observables.

## Current Supplementary Figure S3: periodic tracking (to redistribute)

- Return periodic cue-following results to Figure 1 and move the explanation
  in A–B to proposed S1. Decide the placement of remaining material with the
  final layout.
- Wherever current C is retained, change its red example boxes to match
  current A–B: L = 15 with periods 100 and 400 s.
