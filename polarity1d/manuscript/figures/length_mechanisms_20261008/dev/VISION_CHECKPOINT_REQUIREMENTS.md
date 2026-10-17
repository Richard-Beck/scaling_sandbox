# Three-figure checkpoint: requirements and evidence

Active work resumed on 17 October 2026 after the user agreed the revised scope.
This is the coordination record, not a result or the next review artifact.
The author’s [manuscript vision](../../../../../MANUSCRIPT_VISION.md) governs
the scientific role. The [preceding checkpoint](../../../../../current%20manuscript%20figs/updated_checkpoint_20261017.md)
remains preserved.

## Required final artifact

The final embedded review package will be written to the exact user-requested
path `current manscript figs/vision_aligned_checkpoint.md`. A pointer from the
existing `current manuscript figs` collection will make it discoverable.
The three working figure numbers do not assign final manuscript figure numbers
or replace the reserved new WGD motility-results component.

## Completion evidence

| Requirement | Evidence required before completion |
| --- | --- |
| Figure 1, all five established linear 1D models | Corrected native basal provenance; incumbent organization; full retained gain/recovery size comparisons; actual finite long-cue responses across sizes, with declared native input, duration and post-withdrawal readout. |
| Figure 2, supported WP/Goryachev/Holmes rings | Visual interval/ring correspondence; the same gain/recovery definitions; recovery on basals and declared challenged states across circumference/location; justified long-cue response/outcome maps; highlighted exact native histories and matching global activity vectors. |
| Figure 3, replicated noisy rings | Frozen noise statistics and size-scaling convention; numerical pilot; independent seeds at several circumferences for all three models; deterministic mean/no-noise controls as applicable; vector correlation/persistence, magnitude and local-history diagnostics; physical site-motion diagnostics where identifiable. |
| Scientific presentation | Gain/recovery and finite competition both feature in Figures 1 and 2; model differences, nonmonotonicity and null results retained; native time/amplitude conventions explicit; no assumed velocity law or causal width ablation requirement. |
| Review package | Three readable rendered figure candidates, embedded in Markdown, each with message, rationale for panel/example selection, captions and corresponding experimental measurements. Supporting data, methods, provenance and fresh-context critiques linked. |
| Verification | Appropriate numerical checks cover new integrations, source/input/output provenance and selection. Root inspects actual final renders; final independent review addresses claims and figure coherence. Links/hashes verify; multi-GiB groups recorded in DATA_MANIFEST.md. |

Full completion requires actual results for all three figures. A design-only
noise figure or placeholder ring-recovery box cannot satisfy the resumed goal.

## Evidence accepted so far

- **Linear evidence complete:** the [all-five result report](vision_linear_20261017/RESULTS.md)
  contains all 97 eligible gain/recovery measurements and all 127 finite-cue
  branches (97 primary, 30 duration anchors). Thirteen unavailable preparations
  remain explicit. The [fresh actual-result review](vision_linear_20261017/supporting/metadata/FRESH_RESULT_CRITIQUE.md)
  independently reconstructs all 97 basal/function sources and 127 finite
  histories and finds no essential checkpoint blocker. Root has read the report
  and inspected the assembled all-five figure.
- **Ring evidence complete:** the [ring result report](vision_ring_20261017/RESULTS.md)
  contains 50 basal and 652 challenged-state recovery measurements, full
  separation/circumference finite outcomes, 652 matching activity-vector traces
  and prespecified duration/readout checks. The [fresh actual-result review](vision_ring_20261017/supporting/metadata/FRESH_RESULT_CRITIQUE.md)
  verifies the measurement definitions and retained histories. Root has read
  the report and inspected both continued figure sheets. Positive conditional
  recovery remains distinct from loss of activity in the evolving original
  region.
- **Noise work remains active:** the [reviewed numerical pilot](vision_noise_pilot_20261017/RESULTS.md)
  justified the separately frozen 243-case main ensemble: 216 independent
  noisy realizations, 18 mean controls and 9 zero-input controls. Actual main
  results, seedwise uncertainty, final plots and fresh actual-result review
  are still required. The analysis-only [numerical amendment review](vision_noise_20261017/supporting/metadata/FRESH_NUMERICAL_AMENDMENT_REVIEW.md)
  validates evaluation of the unchanged correlation estimator; it is not a
  substitute for review of the completed ensemble.

The final Markdown entrypoint, complete Figure 3, whole-package critique and
final provenance/navigation/storage audit remain outstanding at this stage.

## Ownership and staged execution

- `all_five_linear_assay_design`: full linear recovery and long-cue evidence,
  fresh `vision_linear_20261017` outputs.
- `chemical_ring_functional_campaign`: matched ring recovery and finite-cue
  evidence, full-grid vector descriptors and components,
  fresh `vision_ring_20261017` outputs.
- `vision_noise_campaign`: bounded numerical/noise pilot and frozen design,
  then replicated ensemble after root reviews the pilot; fresh
  `vision_noise_20261017` outputs.
- `compact_integrated_chemical_figure`: render-only composition and final
  embedded package, fresh `vision_aligned_checkpoint_20261017` outputs.
- Root: scientific scope, integration, user communication, final checks and
  navigation. Native cores and preceding campaigns remain immutable.

Compute-heavy new work uses SLURM and the existing `scaling_sandbox` environment.
Only required compact measurements and selected full histories should be
retained; storage policy applies to all newly generated output groups.
