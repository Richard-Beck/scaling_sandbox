# Selected polarity manuscript evidence and retirement inventory

Decision date: 2026-10-03. `polarity1d/manuscript/` is the single working home
for the corrected experiment, its figures, and selected intermediate data.
The active manuscript uses corrected basal establishment and its ensuing
reversal protocol. **The corrected five-model rerun completed and passed independent audit: 73,975 trials, 375 successful readouts, 207 censored, and 78 preparation-ineligible readouts.** Old linear/width
results are excluded from the main figure and supplements; they remain only
as selected experiments in the separate personal native-field report.

## Keep in the active manuscript core

| Retained material | Canonical destination | Scope |
| --- | --- | --- |
| Corrected basal establishment compact records, 2,476 selected native checkpoints and selected arrival histories | `data/establishment_basal/` | Authoritative current preparation evidence, including the LEGI negative result. |
| Corrected establishment model-grouped heatmaps and combined PDF | `figures/supplement_establishment/` | Current supplement. |
| Corrected protocol illustration and completed corrected reversal length trends | `figures/` and `code/` | Main figure; only the fresh matched comparison supplies reversal results. |
| Model equations, native solvers, corrected cue/readout definitions and frozen scales | `code/core/`, `data/inputs/fixed_scales.json` | Corrected execution and reproduction. |
| Original localized-patch native-field report, all 42 full native fields and metadata | `data/reversal_fields/`, `figures/diagnostic_reversal_fields/`, `code/reversal_report/` | Personal research reference, excluded from main/supplement inference. |

`data/provenance/data_sources.json` records retained immutable files' original
sources, canonical destinations, SHA256, byte count, and role. Copy hashes
match originals. The diagnostic's derived record subset explicitly records its
original source hash and selected scope. Generated figures have independent
hashes after deliberate updates. Original scientific campaigns remain unchanged.

## Minimal personal diagnostic support

`data/reversal_fields/original_campaign/` contains original native solver/cue
source code, protocol, campaign manifest/provenance, and a compact subset of
**42 selected cases and 42 selected trials**. Its record file is approximately
14 KB compressed, rather than a duplicate of the complete width campaign.
`data/reversal_fields/reconstruction_inputs/` preserves 54 unique original
checkpoint/trace files needed by those 42 experiments (about 102 MB). Original
repository-relative identities remain explicit for reconstruction.

This older diagnostic includes known zero-total-input LEGI artifacts and Otsuji
starting-reference ambiguity. Its presence does not validate those experiments
for the corrected comparison. Full native histories support rebundling and
inspection without rerunning the old campaign. Its frozen source files are
scientific provenance, not active manuscript campaign entry points.

The newly copied `data/linear_v11/`, `data/width_reversal/`, baseline-record
and normalized-v10-record duplicates were removed after extracting this minimal
support. No additional historical archive of these new duplicates was created.
`code/retire_legacy_data.py` records the extraction and authorized duplicate
cleanup; `data/provenance/legacy_slimming_plan.json` states the exact scope.
The initial legacy consolidation script was removed from the active core.

The corrected establishment already has physical local copies of the six
unchanged models' results. Selected checkpoints/histories are independent of
old-campaign links. `data/establishment_basal/support/ref15_relative20/` retains
compact predecessor provenance explaining the 2,772 unchanged non-LEGI rows.
The full corrected scalar trial records and withdrawal correction audit are
retained. Complete unselected duration-scan histories remain outside this
selected bundle; no claim of a full 105,947-trial re-audit is made.

## Deprecated and marked for later deletion review

Paths below are relative to the project root. Sizes are approximate allocated
storage measured with `du -sh` on 2026-10-03; they are not promised recoverable
space, especially where copies/hard links or independent workstreams exist.

| Path | Size | Retirement reason and prerequisite |
| --- | --- | --- |
| `tmp/manuscript_figures/{current_01_assay_explanation,pilot_02_width_and_response_definitions}.png` and their old suite source figures | review copies | Historical visual references only. The current main uses the corrected protocol; old width readouts are excluded from the main and supplements. |
| `polarity1d/four_assay/pilots/polarity_establishment_20261003/` | 30 GB | Earlier two-width establishment campaign; superseded by standardized-reference corrected-basal evidence. Preserve unique history until archive decision. |
| `polarity1d/four_assay/pilots/polarity_establishment_ref15_relative20_20261003/` | 42 GB | Incorrect zero-total-input LEGI withdrawal; superseded scientific interpretation. Its compact non-LEGI reuse provenance is retained. |
| `polarity1d/four_assay/pilots/polarity_establishment_basal_copy_attempt_20261003/` | 6.7 GB | Failed/intermediate copy attempt. Review unique files before pruning. |
| `polarity1d/four_assay/pilots/polarity_establishment_submission_attempt_20261003/` | 3.2 GB | Scheduler/submission development residue. |
| `polarity1d/four_assay/pilots/polarity_establishment_integration_20261003/` | 56 MB | Superseded development integration smoke campaign. |
| `polarity1d/four_assay/pilots/polarity_establishment_ref15_relative20_integration_20261003/` | 61 MB | Superseded reference-width integration smoke campaign. |
| `polarity1d/four_assay/pilots/polarity_establishment_smoke_20261003/` | 8.1 MB | Development smoke results; tests/source survive in core. |
| `polarity1d/four_assay/pilots/width_full_lengths_20261002_prelaunch/` | 16 MB | Prelaunch snapshot; authoritative completed frozen sources survive. |
| `polarity1d/four_assay/pilots/width_mean_slope_20261002/` | 19 GB | Earlier width campaign, excluded from current manuscript. Only selected personal-report reconstruction traces remain in the core. |
| `polarity1d/four_assay/pilots/full_reorientation_width_20261002/` | 74 MB | Earlier pilot, excluded from current manuscript. Only selected personal-report reconstruction traces remain in the core. |
| `polarity1d/four_assay/pilots/linear_shared_20261001/` | 274 MB | Earlier shared-linear exploratory comparison. |
| `polarity1d/four_assay/pilots/physical_cue_20261002/` | 1.1 GB | Earlier physical-cue exploratory comparison. |
| `polarity1d/four_assay/pilots/stage1_endpoint_audit_20261002/` | 3.1 MB | Prior starting-reference audit. Retain diagnostic conclusions in assumptions notes; old output is historical. |
| `polarity1d/four_assay/results/reversal_field_report_linear_wrong_campaign_archive_20261003/` | included in old results | Reconstructed the wrong scientific campaign; explicitly superseded. The selected localized-patch report survives. |
| `polarity1d/four_assay/results/length_response_profiles*` and `length_response_profiles_fields/` | included in old results | Earlier linear comparison companion outside selected figure plan. |
| `polarity1d/four_assay/results/length_response_and_settling.*`, `length_response_gradients.*`, and `index.html` | included in old results | Prior settled-linear suite presentation; excluded from current manuscript; no complete duplicate is retained in the canonical bundle. |
| `polarity1d/four_assay/work/` | 2.5 GB | Full revision-11 workers/diagnostics. Only the personal-report selected reconstruction inputs remain in the core; other payload is deprecated historical evidence. |
| `polarity1d/four_assay/reference/` | 1.7 GB | Previous preparation chains/continuation evidence. Frozen corrected scales survive; old linear chains are excluded from current manuscript. |
| `polarity1d/four_assay/pilots/width_full_lengths_20261002/` | 22 GB | Original old-protocol width campaign is deprecated and excluded from active figures. Selected personal-report reconstruction inputs survive; the newly copied full canonical duplicate was removed. |
| `polarity1d/four_assay/pilots/polarity_establishment_ref15_relative20_basal_20261003/` | 1.7 GB | Original authoritative corrected campaign becomes frozen provenance archive after canonical verification; retain before any source cleanup. |

All old `four_assay` execution entry points and exploratory plotters are deprecated
for manuscript work under the corrected manuscript scope.
The originals remain unchanged; retirement labels describe workflow selection,
not removal of numerical or scientific evidence.

## Historical and unrelated material protected

The following remain outside this focused cleanup. Their presence is not an
endorsement of old figures as current manuscript evidence.

- `polarity1d/four_assay/baseline/` (4 MB): explicitly immutable accepted v8 reference.
- `polarity1d/archive/four_test_suite_legacy/` (10 GB): historical suite, including
  the explicitly approved `publication_figures/assay_explanation_v2.*`, its
  generating scripts, reconstructed fields, and provenance. Do not delete that
  keeper while retiring surrounding history.
- `polarity1d/archive/four_assay_reproduction_20261001/` (2.8 GB),
  `four_assay_normalized_v9_20261001/` (1.4 GB),
  `four_assay_conditioning_v10_20261001/` (120 MB), and
  `four_assay_v8_pre_normalization_20261001/` (7.1 MB): independent/historical
  baseline reproduction evidence; review separately before deletion.
- `polarity1d/src/`, `scripts/`, `campaigns/`, `results/data/`, and
  `results/figures/`: broader WP1 qualification/deBelly manuscript work, with
  links from `reports/TRIAL_REPORT_GENERATING_TRAILS.md`; not automatically
  obsolete because the seven-model suite was streamlined.
- `reports/figure_revision/`, `tmp/manuscript_figures/main_03_migration_OU.png`,
  `tmp/manuscript_figures/supp_04_Hosny_sensitivity.png`, and
  `literature_review/migration_data/`: migration/OU and literature evidence;
  unchanged by this polarity cleanup.
- `M2072`, `LFCT`, `migration1d`, and all unrelated repository changes:
  outside scope. Existing root `RETENTION.md` still governs their retention.

## Live producer and stale-link audit

The SLURM queue inspected on 2026-10-03 contained only interactive container and
RStudio sessions; no polarity campaign producer was running. This observation
is a snapshot, not authorization to modify inputs of a later live campaign.

No symbolic links were found in the old pilot/result trees in the bounded link
audit. Existing records contain old relative or absolute path strings and some
old HTML reports link outside their directory. The canonical copy manifest and
diagnostic resolver preserve these identities and resolve selected personal-report inputs;
unselected raw-trial pointers remain archival. Old `reports/CURRENT_STATUS.md`,
`polarity1d/README.md`, and temporary review-set provenance require explicit
handoff links to the canonical bundle before users treat them as current.

The establishment supplement report may link to the canonical data directory for
records/manifest/checkpoints; copied HTML must not imply that every duration trial
is locally bundled. Report generation times remain visible, with no timed browser
self-refresh.

Before any later physical deletion: verify immutable canonical hashes, regenerate
the selected figures/report, check selected operational path resolution, and
record where unique unselected historical evidence will reside. Git ignore rules
or regeneration recipes alone are not an independent backup.
