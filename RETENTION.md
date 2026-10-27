# Repository retention and commit preparation

## Working policy

**Figure 2 exception, authorized and executed 2026-10-22:** the user requested
physical cleanup of Figure 2 data outside the evidence needed by
`manuscript_outline/revision_round3_20261022/`. The exact
[assessment](reports/figure2_cleanup_round3_20261022/assessment.json),
[deletion ledger](reports/figure2_cleanup_round3_20261022/delete_plan.jsonl) and
[execution record](reports/figure2_cleanup_round3_20261022/execution.json)
bound this exception. It supersedes earlier blanket preservation of the
inventoried unused native histories only; it grants no standing permission for
further deletion. Removed 33,821 unused NPZ files and 4,754 empty terminal logs.
Current width/train histories, all calibration and numerical checks, exact
inputs/basals, source code, scalar results and delivered figures remain.
Unselected historical parameter sweeps no longer have complete native-field
coverage; their summaries/provenance remain, and their complete raw histories
require rerunning. The current main and 101-page supplement were rebuilt and
verified after cleanup; see the
[post-cleanup verification](reports/figure2_cleanup_round3_20261022/after_verification.json).
Figures 3/4 and unrelated workstreams were outside the deletion scope.

Keep the active polarity1d and LFCT workstreams. Preserve deprecated code and
historical material while their unique scientific value is reviewed.
This policy generally changes Git eligibility only: it does not authorize
deleting data, retiring implementations, or promoting provisional results.

**M2072 exception, authorized 2026-10-17:** the user requested almost complete
removal, retaining a high-level record and a couple of important models. This
supersedes the earlier M2072 preservation rules. Keep its retired landing notes
and `literature_models/M2072/retirement_20261017/`; remove the old `tmp/`,
`results/`, `archive/` and `src/` trees after verifying the compact retained
files. This is deletion, not relocation to another archive. The polarity1d
storage audit does not itself authorize further deletion.

| Material | Git policy |
| --- | --- |
| Source, model assumptions, licences, tests, configurations, launchers, and acquisition scripts | Eligible; stage by workstream after review. |
| polarity1d/results/data and results/figures | Explicitly eligible curated evidence, including compressed figure source data and the pixel-comparison checkpoint. Review new additions for relevance and size. |
| polarity1d/results/campaigns | Ignored raw/regenerable campaign outputs. Reproduction instructions and campaign inputs remain eligible outside results. |
| Retired M2072 record | Keep the retirement README, two model XMLs, compact parameters/summaries and provenance manifest. The user authorized removal of the former campaigns and development pipeline on 2026-10-17. |
| migration1d source/configs and writer bundle | Eligible historical material; generated results and figure cache are ignored. Figures and prose remain candidates for later curation. |
| Root data and downloaded biosensor payloads | Ignored; preserve acquisition scripts and add compact source manifests deliberately outside ignored payload directories when needed. |
| LFCT | Keep code, config, source URL/hash manifest, and QC decisions. Its own .gitignore excludes annotation payloads, results, and delivery reports. |
| Generated top-level HTML reports, LFCT viewer assets, logs, caches, and builds | Ignored; analysis and report source files remain eligible. |

Eligibility does not mean that a file has been selected or validated for a
commit. Avoid blanket force-adds to bypass these boundaries. Add explicit,
narrow exceptions when selecting further scientific evidence. Ignoring a path
also does not remove files already tracked in Git.

## Reproducibility and storage

A source checkout will not contain every output linked by workstream READMEs.
M2072 now retains selected examples and summaries rather than a complete
campaign archive; its old trajectories and pipeline are intentionally discarded.
Raw LFCT inputs, segmentation masks, historical outputs, and expensive
simulation results need separate durable storage. No external backup location
has yet been verified. Do not delete local outputs on the assumption that Git
or a regeneration recipe is a backup.

For later evidence selection, retain the originating configuration, source and
input identities, software environment, seeds where applicable, checksums,
summary data, scientific limitations, and a reproduction entry point. Record
an actual archive location before replacing unique local data with a pointer.

## Commit sequence

1. Establish layout documentation and retention rules, and reconcile stale
   staging without changing scientific files.
2. Review polarity1d source and curated evidence together; reconcile its
   qualification caveats, links, and reproduction entry points.
3. Review the compact retired M2072 record and LFCT independently, including
   their provenance and selected evidence. Check earlier LFCT tools for unique
   functionality before consolidation.
4. Review historical material and relocated literature workflows. Pair old-path
   deletions with their replacement files where appropriate; do not commit the
   deletions as a standalone cleanup. Decide website retirement separately.

Step 1 stages only .gitignore, README.md, and this policy. Existing scientific
changes, environment changes, HPC notes, relocations, and deletions remain
unstaged for their own review. No commit is created by this preparation.

Before clearing the old index, its binary patch, index copy, HEAD identifier,
and short status were saved under `.git/local-recovery/step1-*/`. This local
recovery record includes staged files no longer present in the working tree;
it is not a backup of the full workspace. To inspect an old staged file without
affecting the current index, use `GIT_INDEX_FILE` pointing at that saved index
with `git show :path/to/file`. The binary patch is the portable recovery record
for staged changes relative to the recorded HEAD; applying it is a separate,
deliberate recovery operation.


## Publication curation, 2026-10-17

The root ignore rules now protect unselected historical `polarity1d/archive`
outputs, four-assay work/results/reference payloads and compacted trial archives,
`reports/figure_revision/data`, exploratory earlier `wgd_langevin` branches, and
`tmp`. These boundaries preserve local material. They do not authorize deletion.
Historical Markdown remains eligible, together with the approved assay-explanation
keeper and selected Figure 3/Hosny assets. Current manuscript evidence and any
new selected native input must be added through a narrow, reviewed exception.
The current migration analysis is `wgd_langevin/refined_ou_20261003`.

The authorized publication task supersedes the earlier three-file-only staging
step: stage reviewed source, tests, configurations, licences, compact evidence,
and the manuscript outline as coordinated workstreams. Pair verified relocations
with old-path deletions. Omit unexplained deletions until their replacement or
intentional retirement is recorded; files recoverable from Git history can be
preserved at their old paths without inventing a relocation.


## Stage-2 preservation decision, 2026-10-17

Preserve the reviewed earlier WGD scientific source/documentation through the
exact historical-source selector, with raw results and native histories still
ignored. Retain the primary compact v8 summary/figure exports once, the historical
writer workbook referenced by its methods/catalogue, and the compact source XML;
duplicate frozen baseline output copies remain local. Source preservation does
not promote old fits, earlier centered-endpoint deBelly histories, or diagnostic
reports to current manuscript evidence.

The **81 preexisting absent tracked originals will remain absent in this local
workspace**. Their intact copies remain published in commit `6199cdc` and the
current retained Git tree, and their paths/sizes/hashes are recorded in
[`publication_held_deletions_pre_recovery_20261017.json`](reports/publication_held_deletions_pre_recovery_20261017.json).
This bounded decision preserves the earlier local workspace state without
assuming its intent: no original files are recreated, and no deletion of those
81 paths is staged or published. The 31 verified relocations already committed
in V1 remain in effect. A clean status is not required to preserve source safely;
independent source-checkout validation must distinguish these known local
absences from the committed tree. Website retirement remains outside this change.

<!-- figure4-round3-compact-storage-20261022 begin -->

## Figure 4 round-3 retention exception, authorized 22 October 2026

The user authorized replacement of OLD Figure 4 noise native HPC outputs
**only after the new compact campaign is COMPLETE**. Its current state at
2026-10-23T03:16:44.915711+00:00 is **COMPLETE**, with
45,696 verified compact seeds.
The scoped old Figure4 noise native profiles were deleted after the COMPLETE audit; the separate retirement execution ledger records actual reclaimed storage. Preserve all preceding native arrays while any
allocated seed or final precision assessment is missing. This exception does
not cover Figures 2/3, exact qualified basals, source/protocol/audit records,
or repository round-2 compact CSVs and figures.

Keep the new `polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022/` compact NPZ files, immutable
manifest/frozen sources, exact 20 basal inputs, per-seed integrity and numerical
receipts, CSV index, 2000-resample whole-seed precision reports, allocation
decisions and actual SLURM receipts. Seeds 0/1 retain float16 examples bounded
at 400 spatial bins and 2000 frames. No full native profile history is retained
for new seeds. Unresolvable rates remain explicitly undefined; finite
SE(log k)>0.07 triggers batches of 96 to at most384. Optional largest/weak
extensions require a recorded condition-specific allocation and remain
separate from the original compact evidence.

[Current allocated storage record](reports/figure4_round3_storage_20261022.json).
The immutable initial manifest says small QoS; effective xxlarge throughput
resources are recorded in `results/qos_burst_upgrade.json`. The measured shared
filesystem free space is distinct from unknown personal quota headroom.

<!-- figure4-round3-compact-storage-20261022 end -->

<!-- figure4-round3-old-native-retirement-20261022 begin -->

## Figure 4 old native-profile retirement completed, 23 October 2026

Measured **2026-10-23T03:16:44.863818+00:00**. Status **EXECUTED_AFTER_COMPACT_COMPLETE**.
The exact [reviewable plan](reports/figure4_round3_retirement_20261022/summary.json)
and [path/stat ledger](reports/figure4_round3_retirement_20261022/delete_plan.jsonl)
cover 72 old Figure 4 noise campaign roots:
**241,709,301,760 allocated bytes**
(225.109 GiB), including directory
and symbolic-link allocation and counting hard-linked files once.
The 12,596 selected native profile/coupled-state paths
occupy **166,686,736,384 unique allocated bytes**;
estimated reclaim is **166,686,736,384 bytes**
(155.239 GiB), allowing for preserved hard links.

Selection is limited to `primary.npz`, `dense_primary.npz`, `full_states.npz`,
`final_state.npz`, and full native-field CSV shards beneath
`results/native_field_CSV_archives/` in those named roots. Exact qualified basal
inputs, frozen sources/protocols and all repository round-2 CSVs/figures remain.
The **75,022,565,376 allocated bytes** of remaining
scalar and descriptor evidence are measured scientific records: legacy
`trace.npz`/scalar histories support historical M1/M1initial and q1=P diagnostics,
matched-mean forced flags and 12-seed comparisons; compressed contour/track JSON
preserves previously measured region identities and sensitivity audits. These
are not retained full native profile histories and are not assumed negligible.

Execution completed after the new compact campaign's COMPLETE audit, final
precision receipt and successful final rebuild. The
[execution receipt](reports/figure4_round3_retirement_20261022/execution.json)
records **241,713,340,416 allocated bytes before** and **75,026,046,976 after**,
with **166,687,293,440 bytes (155.240 GiB) reclaimed** across 12,596 paths.
The execution ledger is authoritative for actual removals; the planning
measurements above are historical. No further deletion is authorized by this inventory.

<!-- figure4-round3-old-native-retirement-20261022 end -->

<!-- ring-r1r2-20261027 begin -->

## Ring figures R1 and R2, 27 October 2026

The R1/R2 brief authorizes two new campaigns without changing native kinetics.
Preserve the authoritative records on the HPC, with allocated storage measured
in `DATA_MANIFEST.md`. Neither Git eligibility nor this inventory authorizes
removal of Richard's existing event campaign or the round-3 compact histories.

| Group | Retain on the HPC | Version in Git |
| --- | --- | --- |
| `polarity1d/manuscript/data/runs/ring_r1_fast_cytosol_20261027/` | Every fast-cytosol compact trajectory and 16 native solver checks, exact input states, matched-noise hashes, tracked-region identities, numerical receipts, and decimated profiles for seeds 0/1 only. Full profiles exist transiently in memory for extraction and are not written for other seeds. | Immutable manifest, frozen executable sources, exact five small input states, launch/collection receipts and compact index; no per-seed trajectory archives or logs. |
| `polarity1d/manuscript/data/runs/poisson_ring_modes_20261027/` | Richard's existing per-perimeter Fourier records and all continuation/checkpoint information, unchanged. No new exemplar runs or output streams are added by the ring-figure task. | Required source/design provenance and frozen inputs; selected compact plotted summaries live in the R2 deliverable. Full Fourier records stay on the HPC. |
| `polarity1d/manuscript/data/runs/ring_r2_per_cell_20261027/` | Event logs, native-cadence polarity vector, active modes 0/1 at validated cadence, Spring circumference, and bounded continuation checkpoints. Middle-size records remain shared references to the existing arm. Keep cadence-validation comparisons and extension decisions. | Frozen code/manifest, source identities, cadence-validation tables and selection, storage forecasts, launch and completion receipts; no copied existing-arm output or bulk trajectory streams. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/ring_r1r2_20261027/` | Same small deliverables as Git. | R1/R2 PNG/PDF/SVG, every plotted value/curve/bound with uncertainty and lag, legends, scientific readmes, campaign provenance, request, storage and independent visual-review receipts. |

The R2 storage plan retains only active modes 0/1 and uses a cadence selected
against the unchanged full-record mode-1 estimator and its whole-seed bootstrap
intervals. A separate small native-cadence polarity stream supports the required
2/5/20-frame event windows. The initial forecast is 0.719 GiB, and the forecast
with every planned extension is 2.901 GiB, including event-log allowance; this is
7.15% of the existing arm's initial Fourier allocation. Forecasts are distinct
from allocated-byte measurements. Existing source data and all plotted bounds
remain preserved even when checks differ from the brief or a rate is unresolved.

<!-- ring-r1r2-20261027 end -->
