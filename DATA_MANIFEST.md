# Project data and output storage manifest

This is the high-level inventory of large **local** data and generated outputs
in `scaling_sandbox`. Paths are relative to the repository root. The inventory
groups related files under meaningful campaign or dataset parents; it is not a
file-by-file catalogue or a statement that a dataset may be deleted.

**Current post-cleanup workspace: 650.90 GiB (698.90 decimal GB)**,
measured on **2026-10-22 during 19:59:03–20:00:57 UTC**:
**698,898,661,376 allocated bytes**. `polarity1d/`
accounts for **90.83%** of the workspace;
manuscript `data/runs/` accounts for **549.77 GiB
(84.46%)**. The user-authorized Figure 2 cleanup reclaimed
**409.436 GiB** from its 39 campaign roots;
the whole-workspace net reduction from the fresh pre-cleanup scan was
**409.425 GiB**, including small audit/retention additions.
The [cleanup summary](reports/figure2_cleanup_round3_20261022/summary.json),
[post-cleanup directory measurements](reports/figure2_cleanup_round3_20261022/after.du.tsv)
and [redraw/verification receipt](reports/figure2_cleanup_round3_20261022/after_verification.json)
record the current footprint and retained evidence. Current Figure 2 artwork,
all plotted values, the complete 101-page supplement and other main figures
passed post-cleanup checks. This is a non-atomic filesystem traversal; small
audit/documentation records were completed after it. No external backup is implied.

**Historical pre-cleanup audit: 1,060.15 GiB (1.138 decimal TB)**, measured on
2026-10-22 during **17:25:03–17:26:59 UTC**: 1,138,332,401,664 allocated bytes.
The [original audit](reports/storage_audit_20261022.json) and
[directory measurements](reports/storage_audit_20261022.du.tsv) remain as the
pre-cleanup record. Its native-coverage claims are historical for retired
Figure 2 arrays. The fresh immediate pre-cleanup traversal at
19:47:01–19:55:14 UTC measured **1,060.325 GiB**.

**Historical post-cleanup workspace: 178.79 GiB (191.97 GB)**, measured on
2026-10-17 during **20:59:40–21:00:35 UTC**: 191,973,122,048 allocated bytes.
Native ring-noise workers and the statistics/report job were terminal, and the
lossless-export compression writer had completed before this traversal. Small
final-assembler/organizer records could still change during the measurement
interval, which is not an atomic filesystem snapshot. The native/noise-export
group and updated figure parent are measured explicitly below. The earlier
complete scan was approximately 381 GiB; the 17:29 post-cleanup measurement was
168.57 GiB before subsequent manuscript production and Git preparation. Both
authorized bulk cleanups are complete and their retained inputs/records have
passed post-cleanup checks.

## Measurement and maintenance

- Include every file or related output group of approximately **2 GiB or more**.
  Include smaller groups where they explain an important retained dataset or a
  retirement decision. A parent entry covers its related descendants; detail
  rows are components of their parent, not additional storage.
- Record **allocated storage in GiB** (1 GiB = 1,073,741,824 bytes), the path,
  purpose, importance/status and measurement date. Logical file size and quota
  usage can differ from allocated storage.
- Measure with `du -x -B1` or an equivalent allocated-block count. Do not follow
  directory symlinks or count external storage as local. `du` counts a hard-linked
  inode once within a traversal, so separate subtree measurements can overlap;
  never interpret the sum of independently measured rows as reclaimed space.
- Tasks that generate, substantially enlarge, move or remove large outputs must
  update this file before handoff, as required by [AGENTS.md](AGENTS.md). Prefer
  extending an existing campaign/group entry to listing thousands of files.
- Source control and a regeneration command are not backups. No external backup
  of all local outputs has been verified. Follow [RETENTION.md](RETENTION.md)
  and workstream-specific scientific/provenance rules before any removal.

## Workspace overview

This table comes from one complete `du -x -B1 --max-depth=5 .` traversal during
**19:59:03–20:00:57 UTC on 2026-10-22**, counting each hard-linked inode once
across the workspace and excluding external filesystems and symlink targets.
Detailed older inventory rows retain their original dates and coverage status;
parent/child rows overlap and must not be added.

| Root-relative area | Allocated GiB | Main purpose |
| --- | ---: | --- |
| `polarity1d/` | **591.22** | Manuscript native simulation evidence dominates: 549.77 GiB in `manuscript/data/runs/`, plus retained historical archives and qualification/reproduction data. |
| `ht1080biosensors/` | **29.96** | Biosensor source archive, extracted native movies and derived analyses. |
| `literature_models/` | **13.97** | Publication-specific implementations/reproduction data, dominated by the Marée reconstruction; M2072 remains a compact retired record. |
| `reports/` | **5.17** | Figure-revision experiment data, report sources/assets, LFCT viewers and compact audit receipts. |
| `data/` | **4.12** | LFCT source images, annotations and early segmentation variants. |
| `.git/` | **2.57** | Versioned scientific source, portable evidence and repository history; native campaigns remain outside Git. |
| `migration1d/` | **1.62** | Historical moving-interval framework and generated campaign results. |
| `lfct/` | **1.02** | Current segmentation labels, annotations, assignment/QC tables and delivery copies. |
| Remaining root areas and root files | **1.25** | Smaller source, literature/OU analyses, manuscript outlines, working figures and temporary assets. |
| **Whole workspace** | **650.90** | **698.90 decimal GB**; local allocation, not personal quota usage. |

### Completed Figure 2 cleanup, 22 October 2026

The user authorized physical retirement of Figure 2 outputs outside the evidence
needed for `manuscript_outline/revision_round3_20261022/`. The cleanup removed
**33,821 unused NPZ files and 4,754 empty terminal scheduler logs**. The sealed
[inode/file ledger](reports/figure2_cleanup_round3_20261022/delete_plan.jsonl)
accounts for **409.435 GiB of removed file allocation**, with no hard links
outside the deletion set. The measured campaign-parent reduction was
**409.436 GiB**; directory allocation changes and new retirement notices explain
the small difference. Nothing outside the inventoried Figure 2 worker/work/log
paths was deleted. All frozen sources, configurations, inputs/basals and
nonempty scalar/JSON records remain.

The joint 39-root measurement during **19:59:03–19:59:11 UTC** reduced the
Figure 2 campaign footprint from **624.477 GiB to 215.041 GiB**. The current
whole-workspace traversal separately measures `polarity1d/manuscript/` at
**567.438 GiB**, including **549.774 GiB**
in `data/runs/` and **11.844 GiB** in `figures/`.
These are nested components, not additional storage.

| Disjoint Figure 2 scope | Before GiB | After GiB | Purpose and scientific importance/status |
| --- | ---: | ---: | --- |
| 20 October Figure 2 signature campaigns | **510.639** | **174.174** | 26 `signature_*20261020*/` roots excluding Figure 3/4: preserve all 450 width workers, selected 180-condition temporal evidence, all calibration/focal-support/numerical campaigns, exact inputs and scalar results. Unselected train/noise/pilot full fields were retired. |
| 19 October superseded Figure 2 compensation/width | **73.262** | **0.293** | 12 exact `revision_five_model_*`, `revision_compensation_*` and `revision_width_*` roots listed in the summary: obsolete native sweeps retired; basals, frozen sources, configurations and scalar results remain. |
| 22 October current Holmes evidence | **40.576** | **40.574** | `signature_holmes_retention_round2_20261022/`: all 32 native train/noise/control workers preserved; only empty terminal logs removed. |

Retain native phase and combined-trace files for all **450 plotted width and
180 temporal conditions**, the **32** current Holmes workers including noise
controls, all calibration and numerical qualification evidence, and the **97**
canonical panel-A basal checkpoints. **8,757 native files** retain their original
inode, size and modification time; **350 source/delivery hashes** verify exactly.
The post-cleanup redraw reproduces the current PNG and all five plotted-data
CSVs byte for byte; the main PDF and all **101** supplement pages retain identical
decoded contents. The existing independent manuscript verifier passes all 450
width and 180 temporal coordinate checks and confirms unchanged other main
figures. No simulations or retuning were performed.

The remaining 215.041 GiB includes authoritative current native evidence,
calibration/numerical controls and historical scalar/provenance records; it is
not a minimal artwork-only package. Unselected old full-field sweeps require
rerunning for arbitrary native re-analysis. Complete historical supplements and
scalar results do not replace the retired fields. Each affected campaign has
an `OUTPUT_RETIREMENT_20261022.json` notice. See the
[assessment](reports/figure2_cleanup_round3_20261022/assessment.json),
[execution record](reports/figure2_cleanup_round3_20261022/execution.json),
[scope measurement](reports/figure2_cleanup_round3_20261022/after_scope.measurement.json)
and updated [retention policy](polarity1d/manuscript/RETENTION.md).

### Pre-cleanup manuscript campaign breakdown, 22 October 2026 (17:25 UTC)

At the pre-cleanup 17:25 UTC measurement, `polarity1d/manuscript/` allocated **976.719 GiB**;
its `data/` component allocates **964.941 GiB**, including **959.169 GiB**
in `data/runs/`. Its separate `figures/` component allocates **11.730 GiB**.
All values in this subsection are historical pre-cleanup measurements. The
post-cleanup overview and retirement record above supersede their current
footprints and native-coverage claims; the original measurements are retained
for growth and retirement accounting.

The following disjoint campaign families partition the **224 direct run
roots**. Their sum plus 98,304 bytes allocated by the run parent itself equals
the measured `data/runs/` total. The audit JSON gives every exact campaign path,
allocated size and family; **50 campaigns individually exceed 2 GiB**.
These rows are native/analysis roots only; figure/report exports are separate.

| Paths under `polarity1d/manuscript/data/runs/` | Allocated GiB | Purpose and scientific importance/status |
| --- | ---: | --- |
| `signature_*20261020*/`, excluding `signature_figure3_*` and `signature_figure4_*` (26 roots) | **510.639** | Figure 2 signature assays: calibration, width, cue trains, bounded noise, support assays and numerical companions. Full native coupled fields, phase/restart evidence, per-frame readouts and analysis generations are retained, including superseded pilots and failures. This preserves scientific re-audit capability but includes repeated storage representations. |
| `signature_figure4_*20261020*/` (62 roots) | **214.620** | Figure 4 noisy-ring stages, long horizons, region analyses, recovery attempts and native-field atlas data. Preserve seed identity, original failures and finite-window qualifications; these outputs do not establish a universal persistence law. |
| Other retained native campaigns (82 exact roots in audit JSON) | **90.520** | Earlier corrected comparison and mechanism evidence, mixed historical campaigns, preparations and independent controls. Includes the 28.617-GiB corrected deBelly no-protrusion experiment. Earlier centered-endpoint deBelly histories remain deprecated evidence; a retained history is not automatically expendable. |
| `revision_*20261019*/` (36 roots) | **87.422** | Figures 2–4 revision native fields/checkpoints, exact basals, full-length compensation, stationary-cue assays and noise replications. This scope excludes separate report exports and the compensation-review mapping root included in the historical combined revision row. |
| `signature_holmes_retention_round2_20261022/` and `signature_figure3_holmes_extensions_round2_20261022/` | **40.940** | Completed round-2 Holmes trains, bounded-noise seeds/controls and exact-endpoint cue-free continuations. Native fields, checkpoints and provenance are authoritative evidence for the revision. |
| `signature_figure3_*20261020*/` (16 roots) | **15.028** | Fixed-size distance assays, regional handovers, pulse/width trials, mirrors and numerical controls, with original rejected analyses retained as provenance. |

Largest individual campaign components, already included in the families above:

| Path under `polarity1d/manuscript/data/runs/` | Allocated GiB | Purpose/status |
| --- | ---: | --- |
| `signature_train_core_20261020/` | **194.257** | Main cue-train native histories and repeated phase/combined-trace representations; see file breakdown below. |
| `signature_noise_bounded_full_20261020/` | **114.864** | Full bounded-noise trajectory campaign and analysis evidence. |
| `signature_figure4_E4_horizon_20261020/` | **105.231** | Long-horizon noisy-ring native histories and controls. |
| `signature_holmes_retention_round2_20261022/` | **40.576** | Today's completed Holmes train/noise follow-up. |
| `signature_train_holmes_ready_20261020/` | **38.402** | Holmes cue-train extension with qualified preparations. |
| `revision_five_model_compensation_strength_v2_20261019/` | **36.478** | Additional-strength five-model compensation arm. |
| `signature_train_holmes_extended_20261020/` | **34.520** | Extended Holmes cue-train histories. |
| `revision_five_model_compensation_baseline_v2_20261019/` | **29.348** | Baseline compensation protocols across all qualified model/length states. |
| `debelly_no_protrusion_20261010/` | **28.617** | Unique corrected endpoint-protrusion ablation histories. |
| `signature_figure4_native_atlas_final_20261020/` | **27.264** | Derived native-field atlas data/export evidence; distinct from new integrations. |
| `signature_width_core_20261020/` | **24.215** | Main equal-commanded-area width assays. |
| `signature_calibration_core_20261020/` | **21.198** | Native cue calibration evaluations and retained audit evidence. |

A separate file-level read of the 194.257-GiB train campaign on 22 October
found **8,478 `phase_*.npz` files allocating 85.736 GiB**, **972 `trace.npz`
files allocating 87.031 GiB**, and **972 per-trial `readouts.json` files
allocating 15.054 GiB**. Other JSON/analysis files and small sources/logs account
for the remainder. In the sampled `workers/0000/trial/`, NumPy equality checks
verified that the combined trace repeats the phase fields after omitting
duplicate phase-boundary frames; both representations also store identical
`active` and `field_active` arrays. This identifies a storage-efficiency issue,
not a verified campaign-wide reclaimable-space estimate or deletion approval.
Any consolidation must preserve restart information and exact dependent source
identities, which were not exhaustively audited here.

## Retirements authorized on 2026-10-17

| Path | Former allocated storage | Current record and scientific consequence |
| --- | ---: | --- |
| `literature_models/M2072/` | 86.66 GiB | **Retirement complete: 0.91 MiB retained**, 19 files and 472,183 logical bytes. Compact historical notes, original V1–V6/selected leader parameters, two selected XML models and provenance survive in [the retired branch](literature_models/M2072/README.md). Bulk campaigns and development pipeline were deleted; this is not an external archive. Measured allocated-footprint reduction: 93,049,708,544 bytes. See [the completed action record](literature_models/M2072/retirement_20261017/cleanup_status.json). |
| `polarity1d/four_assay/` | 127.02 GiB | **Cleanup and post-cleanup checks complete: 0.88 GiB retained** (944,349,184 allocated bytes), a measured footprint reduction of 135,440,498,688 bytes (126.14 GiB). Compact metadata in 24 verified archives, scalar records, source/protocol provenance, accepted v8 material, final v11 profiles and selected illustration inputs survive. All 5,006 canonical selected-input hashes passed again; data-only regeneration reproduced the current main and width-supplement PNGs byte for byte. Full native coverage of old trials is intentionally retired. The active corrected manuscript workflow is separate. See [the execution/retention record](polarity1d/four_assay/OUTPUT_RETIREMENT_20261017.md), [human cleanup report](reports/storage_cleanup_four_assay_20261017.md) and [machine-readable ledger](reports/storage_cleanup_four_assay_20261017.json). |

## Inventory

**Historical coverage notice:** component rows measured before the 22 October
Figure 2 cleanup retain their original dates. Where affected by the exact
retirement ledger, their old full-native-archive claims are historical, not a
statement that every referenced file is still local. The current retirement
record and per-campaign notices above determine retained native coverage.

Initial component measurements below were taken on **2026-10-17,
approximately 16:37–16:44 UTC**. Updated manuscript parents and the ring-noise native/export
group were remeasured during **20:59:40–21:00:35 UTC**, with dates stated in their rows.
Unchanged component entries retain their earlier measurements. The two retired
trees are reported separately above; the retained datasets below are outside
that deletion scope. Sizes are rounded to two decimals.

### Polarity modelling and manuscript evidence

Read [current report status](reports/CURRENT_STATUS.md), the
[manuscript workflow](polarity1d/manuscript/README.md) and its
[retention inventory](polarity1d/manuscript/RETENTION.md) before interpreting
these outputs. “Corrected” in an early campaign name can refer to basal
preparation, rather than to the later endpoint correction. **deBelly histories
using `center_boundary_motion=True` are deprecated scientific evidence** even
when retained in a mixed-model campaign; the other models are not implicated.

| Path or grouped component | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/` | **97.82** (105,035,956,224 allocated bytes; remeasured 2026-10-17, 20:59:40–21:00:35 UTC) | Working home of the active five-model comparison, native simulation evidence, report/figure sources and development experiments. This is a working manuscript resource, not a completed manuscript. The following manuscript rows are components of this total. |
| `polarity1d/manuscript/data/runs/` | **90.03** (96,669,343,744 allocated bytes; remeasured 2026-10-17, 20:59:40–21:00:35 UTC) | Native histories, frozen execution inputs, checkpoints, scalar records and audits for the retained comparison and mechanism experiments. Individual large campaigns are identified below; smaller campaigns are covered by this parent. |
| `polarity1d/manuscript/data/runs/vision_noise_main_20261017/`, `vision_noise_pilot_20261017/`, `vision_noise_numerical_audit_20261017/`, and related `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_20261017/` exports | **3.31 combined** (3,556,212,736 allocated bytes; joint traversal 2026-10-17, 21:00:34 UTC) | Finite-time chemical-ring polarity-vector and local-pattern noise ensemble, exact retained WP/Goryachev/Holmes basal inputs at circumferences 20/40/60, frozen input/source hashes, native histories, scalar/statistical readouts and report exports. Main run: **2,802,688,000 bytes (2.61 GiB)**; 243 workers comprising 216 fluctuating-input runs (two mean fractions, twelve seeds per model/size), 18 matched-mean controls and nine basal controls. The 30-worker pilot allocates **81,870,848 bytes**; temporal-interpolation/ODE controls allocate **30,064,640 bytes**; related figure/report exports allocate **641,589,248 bytes (0.60 GiB)**. All native workers completed at 20:45:20 UTC; statistics/report SLURM `21126788` completed at 20:57:36 UTC; deterministic gzip-export writer completed before measurement. Seven large raw CSVs (993,719,748 logical bytes) remain alongside lossless gzip copies (221,166,834 logical bytes), recorded in [`publication_compression_20261017.json`](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_20261017/supporting/metadata/publication_compression_20261017.json). These are activity-polarity vectors under declared interpolated OU-knot/lognormal activation forcing, not measured cell velocities or migration trajectories; numerical completion does not establish the final biological interpretation. Preserve the full native run roots, frozen protocol/inputs, actual-time/eligibility records and statistical/figure provenance. Full run source: [`vision_noise_20261017.py`](polarity1d/manuscript/code/mechanisms/vision_noise_20261017.py); orchestration: [`vision_noise_main_workflow_20261017.py`](polarity1d/manuscript/code/mechanisms/vision_noise_main_workflow_20261017.py); main analysis/report: [`vision_noise_main_report_20261017.py`](polarity1d/manuscript/code/mechanisms/vision_noise_main_report_20261017.py). Native components are included in `data/runs/`, exports in `figures/`; this grouped row is not additional storage to add to either parent. |
| `polarity1d/manuscript/data/runs/chemical_ring_functional_v2_20261017/` | **3.67** (3,940,261,888 allocated bytes; measured 2026-10-17, 17:42 UTC) | User-requested systematic native periodic WP/Goryachev/Holmes size/location characterization: 867 planned cells at 17 sizes and 17 separations, exact shared basals, initial cue-gain readouts, matched inactive-withdrawal/return probes at predefined epochs, numerical references and sparse native primary histories. New development evidence linking conditional functional assays to pattern fates; not global stability or freely migrating-cell persistence. Preserve frozen execution sources, full matrix exclusions and native probe returns. A retained-control verification adds 50 matched no-cue controls and 1304 checked finite-time increments, without new scientific conditions. Sources and reviewed main-figure candidates are unchanged. |
| `polarity1d/manuscript/data/runs/chemical_ring_functional_20261017/` | 0.28 (297,435,136 allocated bytes; measured 2026-10-17, 17:33 UTC) | First 36-condition native ring functional pilot, all 12 model/size workers completed. Exact C20/30/40/50 preparations and matching response/probe results are referenced by the expanded v2 campaign; this is a retained input dependency, not a redundant disposable draft. |
| `polarity1d/manuscript/data/runs/debelly_no_protrusion_20261010/` | 28.62 | Unique histories from 430 corrected-model integrations disabling endpoint protrusion: preparation scans, exact native-branch continuations and small-seed diagnostics. Audited development evidence and linked figure; not included in the authorized historical cleanup. |
| `polarity1d/manuscript/data/runs/three_test_mechanisms_20261006/` | 6.68 | Original five-model mechanism replays, neighboring cue doses and diagnostic controls behind the earlier three-test reader. Mixed historical evidence: its centered-endpoint deBelly results are deprecated, while other-model evidence may remain relevant. |
| `polarity1d/manuscript/data/runs/debelly_corrected_three_test_20261008/` | 6.45 | Corrected independent-endpoint deBelly preparations and all three opposing-cue tests, including native trial histories and independent audits. Authoritative replacement for earlier deBelly comparison outputs. |
| `polarity1d/manuscript/data/runs/corrected_reversal_20261003/` | 5.04 | First clean comparison from corrected basal establishment; duration-scan trial data and controls. Retains unchanged non-deBelly evidence; deBelly interpretation must use the later independent-endpoint replacement. |
| `polarity1d/manuscript/data/runs/fixed_slope_reversal_20261006/` | 4.22 | Finite constant-physical-slope opposing-gradient scan, cue-duration brackets, controls and audits. Model-specific native evidence; earlier centered deBelly component is superseded. |
| `polarity1d/manuscript/data/runs/continuous_fixed_slope_20261006/` | 3.10 | Continuous-gradient trajectories and first-orientation-crossing readouts for 97 eligible starts. Different outcome from persistence after cue withdrawal; earlier centered deBelly component is superseded. |
| `polarity1d/manuscript/data/runs/inactive_sink_loss_20261008/` | 2.82 | Goryachev/Wave-pinning inactive-substrate sink experiments, native histories and audits. Retained chemical distance/competition mechanism evidence. |
| `polarity1d/manuscript/data/runs/debelly_mca_threshold_profiles_20261010/` | 2.70 | Corrected deBelly MCA profile/withdrawal development experiments supporting the mechanism report. Conditional trajectories do not establish a universal MCA shutdown threshold. |
| `polarity1d/manuscript/data/runs/length_debelly_review_20261009/` | 2.11 | Replayed corrected independent-endpoint histories and measurement exports for the first length-report review. Native evidence behind later revised figures remains retained. |
| `polarity1d/manuscript/data/runs/wave_pinning_slow_diffusion_20261008/` | 1.77 | Superseded variant reducing both active and inactive diffusivity. Retained historical native histories; the diffusion-contrast replacement is a different experiment. Below the main threshold, recorded because its status matters. |
| `polarity1d/manuscript/data/runs/wave_pinning_diffusion_contrast_20261008/` | 1.63 | Requested replacement with active diffusivity 0.01 and inactive diffusivity 100, verified variant preparations and matched challenges. Separate from the superseded diffusion variant. |
| `polarity1d/manuscript/data/archive/goryachev_exploration_20261008/` | 4.39 | Archived depletion/transport diagnostic arrays, frozen sources and earlier exploratory measurements. Excluded from the active Goryachev presentation; archive status alone does not make unique histories redundant. |
| `polarity1d/manuscript/data/establishment_basal/`, `reversal_fields/` and other non-run data | included in manuscript total | Selected authoritative preparations, 4,952 independently copied establishment NPZ inputs, 54 historical diagnostic reconstruction inputs and compact source manifests. These local copies decouple current inputs from retired four-assay arrays. |
| `polarity1d/manuscript/figures/` | **3.511** (3,770,040,320 allocated bytes; remeasured 2026-10-19, 22:09 UTC) | Report PDF/HTML, vector/PNG exports, revision archives, development candidates, final consolidation and source/artifact manifests. Includes the revision report component below; do not add that component again. Native roots are separate. Historical manuscript/workspace totals retain their original dates. |
| `polarity1d/manuscript/data/runs/*ring*20261017/`, `chemical_ring_distance_pilot_20261017/`, `holmes_distance_robustness_20261017/`, `holmes_size_bridge_20261017/`, `debelly_planar_loop_preparation_20261017/`, `debelly_two_arc_distance_20261017/` | approximately 1.88 combined | Circular preparation, distance, finer-mesh and finite-window robustness experiments for the current main-text framework discussion. Smaller retained run group, already included in the run total. Spring and deBelly constructions introduce declared geometric/closure assumptions; these are not universal stability results. |
| `polarity1d/archive/` | **14.11** | Independent historical suite/normalization/reproduction evidence outside the authorized `four_assay/` tree. The next two rows are components; no further removal is implied. |
| `polarity1d/archive/four_test_suite_legacy/` | 9.90 | Historical suite, sources, reconstructed fields and provenance. Includes the approved `publication_figures/assay_explanation_v2` keeper and its inputs. |
| `polarity1d/archive/four_assay_reproduction_20261001/` | 2.73 | Independent prior accepted-baseline reproduction evidence. The smaller v9 normalization archive (1.35 GiB) and v10 continuation archive are covered by the archive parent. |
| `polarity1d/results/` | **8.96** | Broader WP1 qualification, curated data/figures and native reproduction campaigns; not automatically obsolete when the seven-model manuscript workflow was retired. Campaigns occupy 8.85 GiB of this total. |
| `polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_base/` | 3.70 | Native Julia/standalone Figure 2 baseline reproduction campaign, including conditioning audits, task histories and summaries. Component of `results/`, retained as separate native-model/reproduction evidence. |
| `polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000/` | 1.71 | Longer native Figure 2 panel-A reproduction histories. Component of `results/`; related native evidence below the main inventory threshold. |
| `polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response/` | 1.58 | Native Figure 2 linear-response reproduction histories. Component of `results/`; other smaller campaigns are covered by the parent. |

The current chemical-reference figures reuse retained histories rather than
creating a second full campaign. Exact all-five initial-state metrics and
Goryachev/Wave-pinning substrate-removal/recovery evidence are smaller components
covered by `data/runs/`; their scientific importance is not proportional to size.

### Three-figure checkpoint additions, 2026-10-17

These completed checkpoint components are included in the current workspace
and manuscript parents above. Their smaller component measurements are retained
for scientific provenance; they are not extra quantities to add to those parent
totals. The replicated-noise native and statistical/report exports are grouped
in the updated polarity inventory above, after the final producer and
compression writer completed.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/vision_linear_20261017/` | 1.36 (1,456,201,728 bytes) | Completed all-five linear evidence: exact frozen basal/source inputs, 97 original common active-loss recovery archives and 127 finite-cue native histories, duration/readout sensitivities and independent audits. Authoritative inputs to the new checkpoint; no main archive was cancelled or substituted. Measured 2026-10-17. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_linear_20261017/` | 0.006 (6,193,152 bytes) | Embedded linear report, compact measurement exports, PNG/PDF/SVG figures, designs and fresh reviews. Derived presentation with linked native archives. Measured 2026-10-17. |
| `polarity1d/manuscript/data/runs/vision_ring_recovery_20261017/` | 0.76 (820,133,888 bytes) | Complete common assay on 50 ring basals and 652 exact challenged states, all amplitudes/native fields and independent reconstruction. Required new interval/ring comparison evidence. Measured 2026-10-17. |
| `polarity1d/manuscript/data/runs/vision_ring_duration_20261017/` | 0.13 (143,040,512 bytes) | Prespecified half/double-duration anchor histories and whole-grid finite-readout sensitivity. Distinguishes fixed finite outcomes from universal saturation or attractor claims. Measured 2026-10-17. |
| `polarity1d/manuscript/data/runs/vision_ring_moments_20261017/` | 0.025 (26,927,104 bytes) | Read-only moment extraction from all 652 retained finite histories: global vectors, magnitude, harmonics and availability flags. The underlying original chemical-ring histories remain separate retained dependencies. Measured 2026-10-17. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_ring_20261017/` | 0.006 (6,029,312 bytes) | Embedded ring report, compact sources, correspondence/outcome/vector figures and independent actual-result critique. Measured 2026-10-17. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_pilot_20261017/` | 0.013 (14,147,584 bytes; measured 2026-10-17 21:07 UTC) | Reviewed immutable single-seed pilot report, figures and diagnostics. Historical numerical/regime evidence; the complete replicated ensemble supplies uncertainty claims. This smaller report is additional to the native-pilot/report group listed above. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_aligned_checkpoint_20261017/` | 0.021 (22,880,256 bytes; measured 2026-10-17 21:08 UTC) | Final three-figure assembly on four review sheets, exact used-data/native-history exports, PDF/SVG/PNG files, identity manifest and accepted fresh whole-package critique. The [requested Markdown entrypoint](<current manscript figs/vision_aligned_checkpoint.md>) remains separate from these generated artifacts. |

The completed noise branch's later [final storage handoff](polarity1d/manuscript/data/runs/vision_noise_main_20261017/results/final_storage_footprint.json)
was measured at **21:07 UTC**: 3,570,515,968 allocated bytes (3.325 GiB) across
its five disjoint native-main, main-report, numerical-audit, native-pilot and
pilot-report groups. This extends the earlier joint noise measurement above
by including the small pilot report and final documentation/provenance writes;
it is not additional storage to add to that earlier group. The record gives
each exact path, allocated size, purpose and complete/retained status. All
native outputs, authoritative statistical CSVs and reviewed PNGs remain
unchanged; failed and corrected analysis provenance is retained separately.

### Spring ring development addition, 2026-10-19

These small additions are recorded for scientific continuity. The manuscript
parent measurements above predate them; this is a measurement of the three
named subtrees, not a fresh whole-workspace inventory.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/spring_ring_validation_20261019/`, `spring_ring_api_smoke_20261019/`, and related `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/spring_ring_validation_20261019/` | **0.033 combined** (35,389,440 allocated bytes; measured 2026-10-19, 14:25 UTC) | Declared Spring periodic/global-extension development analogue at resting circumference 2. Native validation: **33,210,368 bytes**, nine completed SLURM workers, 23 complete phase archives, frozen inputs/sources, fresh 200/400-cell preparations, independent solver replay, rotation, finite-pulse and conserved-perturbation checks, independent audits, source-backed geometry/dilution accounting and browser verification. API smoke: **532,480 bytes**, zero-active initialization and exact staged restart on the corrected public adapter. Scientific reader: **1,646,592 bytes**, interactive retained raw profiles, vector equations, compact evidence, provenance and runnable example. The user agreed to carry over the original stretching, dilution and diffusion terms; mean-force circumference dynamics remain an explicit circular adaptation. Bounded numerical qualification is complete, without a size law, migration-persistence or native 2D-mechanics claim. Retain raw evidence and frozen sources. This row is a component of the manuscript area, not storage to add to an updated parent total. |

### deBelly periodic-domain attempt retirement, 2026-10-19

The user requested consolidation and removal of artifacts from this attempt.
This is a measurement of the named retired outputs and retained package only;
the manuscript parent measurements above have not been refreshed. The
[retirement inventory](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_ring_review_20261019/cleanup_plan.json),
[execution receipt](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_ring_review_20261019/cleanup_receipt.json) and
[storage record](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_ring_review_20261019/supporting/metadata/storage_footprint.json)
record the exact paths and measurements. Inventory alone did not authorize deletion;
the user's explicit cleanup instruction did.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_ring_review_20261019/` | **0.122 retained** (131,530,752 allocated bytes; measured 2026-10-19, 20:12:00 UTC) | Single [executive report and seven-field reviewer](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/debelly_ring_review_20261019/index.html) for the incomplete periodic deBelly attempt. Retains all 18 event clocks from the original and two fixed-grid-face C50 challengers at native 0.5-s resolution in event windows, 375 full-precision whole-response context frames, four compact evidence groups, minimal frozen reference model/basal state, provenance and audits. Viewer regeneration and relocation pass. Physical cap/material/reference closure, causal travel/reversal mechanism and numerical convergence remain unresolved. Recommendation: omit deBelly from current ring assays and keep corrected linear deBelly. No migration-persistence or accepted ring size law resulted. |
| Superseded 19 October deBelly run/report/source/cache groups (exact list in inventory) | **0 remaining**; 7,185,711,104 workspace allocated bytes retired (6.692 GiB), plus 146,305,024 temporary bytes | Removed 29 dated run roots, the old development report topic, 123 standalone scripts, 49 caches, the dated status note and 101 task temporary entries after independent retention/dependency checks. Full native histories outside retained event windows, context samples and selected summaries were retired; historical audit receipts do not imply those full inputs still exist. Net workspace allocation reduction after retaining the review: **6.570 GiB**. Established corrected linear results, earlier independent diagnostics, `tmp/x50_mechanisms_20261003/debelly`, current figures and independent Spring work are protected. This component must not be added to a refreshed manuscript parent total. |

### Figures 2–4 revision additions, 2026-10-19

Final allocated storage was measured during **22:08:59–22:09:09 UTC on 19 October** over **38 disjoint roots**: the revision figure group, every direct `data/runs/revision_*20261019*` root (including suffixed roots), and the separate compensation-review mapping root. All scientific producers are terminal; archived failed attempts and unsubmitted strict payloads remain explicit. Only small final documentation/audit records can change after this measurement. The [per-path storage record](polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/supporting/metadata/storage_current.json) records allocated bytes, purposes, producer states and scheduler observations. No deletion was performed or authorized by this inventory.

This is a final measurement of the revision group, not a fresh whole-workspace inventory. Historical manuscript/workspace parent totals retain their original dates. The separately measured report parent includes the revision report component and must not be added to it.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/revision_*20261019*/`, `review_compensation_mapping_20261019/` and `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/` | **88.226 combined** (94,732,214,272 allocated bytes; 2026-10-19, 22:08:59–22:09:09 UTC) | Complete native fields/checkpoints, frozen inputs, independent replay and derived measurements for the requested Figures 2–4 revision. Figure 2 contains 3,492 core protocols at all 97 qualified model–length states across five models, plus two targeted corrected deBelly withdrawals. Figure 3 contains all 645 stationary-cue protocols; Figure 4 contains 936 core cases plus 26 separate transition cases. Full finite-tail failures, Spring sampled-duration limits and archived execution failures remain explicit. Main figures, comprehensive supplements and compact redraw inputs are retained alongside local native evidence. Native archives are local, while the compact figure/reconstruction package is intended for Git delivery. |

Large components of that combined group (included in its total, not extra storage):

| Component path | Allocated GiB | Purpose/status at measurement |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/revision_arc_narrow_grid_20261019/` | **3.162** (3,395,076,096 allocated bytes; 2026-10-19, 22:09 UTC) | Complete native stationary-cue experiment: 645 logical protocols with exact reuses; full profiles/restarts and finite-tail qualifications underlie Figure 3 and its complete supplement. |
| `polarity1d/manuscript/data/runs/revision_five_model_compensation_baseline_v2_20261019/` | **29.348** (31,512,027,136 allocated bytes; 2026-10-19, 22:09 UTC) | Complete five-model/all-97-length, three-width/five-withdrawal baseline: 1,552 logical outcomes including 97 controls; native arrays, exact basals and restart provenance retained. |
| `polarity1d/manuscript/data/runs/revision_five_model_compensation_strength_v2_20261019/` | **36.478** (39,167,950,848 allocated bytes; 2026-10-19, 22:09 UTC) | Complete four-additional-strength arm at middle width: 1,940 positive protocols with exact high-strength reuses; full native fields and scalar/finite-tail status retained. |
| `polarity1d/manuscript/data/runs/revision_noise_strength_broad_20261019/` | **2.904** (3,118,391,296 allocated bytes; 2026-10-19, 22:09 UTC) | Complete 576-case nine-strength pilot coverage and numerical companions; original failed Spring sampled-duration comparison remains explicit. |
| `polarity1d/manuscript/data/runs/revision_noise_strength_replication_scoped_20261019/` | **3.619** (3,885,432,832 allocated bytes; 2026-10-19, 22:09 UTC) | Complete 360-case native replication extension giving 936 core cases; independently qualified occupancy/vector/memory readouts, full seed histories and controls retained. |

### Figures 2–4 signature, distance and noise revision additions, 2026-10-20

Allocated storage was measured from **2026-10-21T03:41:16.734383+00:00** to **2026-10-21T03:45:44.389892+00:00** over **112 disjoint roots**.
All required producers have finished; small documentation/receipt changes can follow this measurement.
The [per-path record](docs/tasks/briefs/20261020_revision_storage_current.json) records allocated bytes, native versus export components and exact paths. Each inode is counted once in this joint traversal. Original failed or superseded work is retained. Inventory does not authorize deletion. Historical parent totals retain their dates.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/signature_*20261020*/` excluding Figure3/4 roots, and `dev/signature_assays_20261020/` exports | **510.669 combined** (548,326,498,304 allocated bytes; native/analysis 548,294,533,120, exports 31,965,184) | Five-model signature assays: exact qualified basals, frozen sources and full coupled fields; 36 calibration-worker audits covering 528 native evaluations; 450 equal-commanded-area width runs; 1,458 defined cue trains; 486 bounded-noise trajectories; supporting assays and numerical companions. The final main and 82-page supplement retain 162 calibration-unavailable Holmes cells, Wave-pinning L19 noise losses and deBelly timing/field sensitivities. Superseded pilots and failures remain preserved; native arrays stay local. |
| `polarity1d/manuscript/data/runs/signature_figure3_*20261020*/` and `dev/revision_program_20261020/figure3_distance/` exports | **16.118 combined** (17,306,558,464 allocated bytes; native/analysis 16,135,503,872, exports 1,171,054,592) | Four-model fixed-size distance assays: 12 qualified basals, 176 A1/mirror histories, corrected full-region handover roots, 545 A2 pulse trials, 381 B trials including mirrors, and mesh/time-step controls. Final geometry/direct-heading displays preserve all bounds and numerical exceptions. The 1,150-history atlas has pixel-identical lossless WebP exports; original PNGs and native arrays remain local. Rejected fixed-footprint and discontinuous-boundary analyses remain diagnostic provenance. |
| `polarity1d/manuscript/data/runs/signature_figure4_*20261020*/` and `dev/signature_figure4_20261020/` exports | **220.436 combined** (236,691,611,648 allocated bytes; native/analysis 230,440,755,200, exports 6,250,856,448) | Four-model noisy-ring revision: corrected E0 inventory and 110 condition analyses; delivered 34-page interim; strict optical-region measurements and independent native-frame checks. Expanded stages contain 864 E1, 468 E2 and 528 E3 native runs. The literal E4 rule selects 52 noisy conditions, giving 624 noisy realizations plus 64 matched mean/zero controls. Original attempts, failed strict prefix comparisons, exact-cap restart pilots and actual-failure-only replacements are retained; selected replacements do not add independent seeds. S11 full native-field CSV.gz archives preserve exact plotted values on HPC; portable atlas figures and compact inventories retain their paths and hashes. The completed final main, 1,291-page statistical core/index and all 13 native-atlas volumes (578 pages) passed independent source/display review. H1 and H3 remain inconclusive and H2 fails its declared overall criterion. Whole-seed uncertainty, finite-window/coverage qualifications and the original Spring sampling failure remain explicit. Superseded export attempts are retained locally and the exact portable publication allowlist is preserved. Superseded compact-basal, pairwise-only and antipodal-route analyses remain diagnostic provenance. Original 936+26 native sources belong to the separately inventoried 19 October group. |
| `manuscript_outline/revision_20261020/` and the `outline_v4` JSON, PDF, build manifest and validation record | **0.008 combined** (8,806,400 allocated bytes; native/analysis 0, exports 8,806,400) | Static integrated report exports and exact source/assembly identities; no native integrations are launched by report assembly. |
| `.git/` after approved revision-payload staging; measured 2026-10-21T03:50:52.293191+00:00 | **3.076 GiB** (3,303,309,312 allocated bytes) | Versioned scientific source and portable figure/evidence exports plus preceding repository history. These duplicate the exported evidence rather than adding independent experiments. Native archives stay outside Git. The [measurement record](docs/tasks/briefs/20261020_revision_git_storage.json) identifies this staging snapshot; small documentation/index/commit objects can follow. No deletion authorized. |

These are component measurements of the revision and repository metadata. Do not add them to historical manuscript/workspace parent totals as if those totals were current.

### Observable-region native pilot and diagnostic archives, 22 October 2026

Measured allocated storage on **2026-10-22, 17:51–17:52 UTC**, with
`du -x -B1 -s`, after native pilot execution and report/artifact finalization.
These named related groups are recorded explicitly; historical workspace and
manuscript parent totals above retain their original dates. Archive component
measurements below are included in the output parent, not additional storage.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/figure4_region_model_20261022/` | **0.988220 GiB** (1,061,093,376 allocated bytes) | Author-corrected native Wave-pinning pilot: 25 conditions, 300 immutable native histories, six detector variants, 150 process fits with 2,000 whole-seed resamples, and 26,800 reduced trajectories across 134 supported ensembles. Includes saved heading/magnitude/occupancy arrays, G1/×2 cadence checks, V0–V5, parameter and validation CSVs, F1–F4 PNG/PDF/SVG plates, 168-page supplement, static timestamped reader, source identities and independent audits. **Requested pilot review stop:** raw G1 flags 3/25, V0 passes 25/25, V1 passes 15/21 estimable conditions with four unsupported; magnitude prediction and observable attribute feasibility remain failed population-closure diagnostics. No other model was run. Synthetic recovery is retired as a gate and was not extended. Both diagnostic archives remain inside this parent: `archive/recovery_diagnostic_20261022/` measures **10,731,520 bytes** (97 original files verified unchanged), and `archive/observable_width_diagnostic_v1_20261022/` measures **450,207,744 bytes** (225 relocated files verified unchanged; initial support-boundary flaw). Valid native fits were extracted and reused unchanged for the corrected replay. Preserve negative findings, uncertainty qualifications, archives and exact native/source identities; no native histories were changed or deleted. |
| `polarity1d/manuscript/code/mechanisms/region_model_20261022/` | **0.005829 GiB** (6,258,688 allocated bytes) | Observable-region detector/tracker, condition fitting, corrected support geometry and simulation, validation, report and independent-audit sources; immutable executed gate/fit/report snapshots and diagnostic historical code. Current state retains no hidden region identities; touching arcs are unioned. Sources distinguish numerical consistency from remaining scientific closure failures. Preserve executed sources and archive provenance; prepared report snapshots that were superseded before execution remain identified separately. |

### Literature-model implementations and reproduction data

| Path or grouped component | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `literature_models/`, excluding retired `M2072/` | **13.97** (root measured 2026-10-17, 20:59:40–21:00:35 UTC; compact M2072 record is below 1 MiB) | Native/provenance implementations and publication-specific reproduction bundles. Most storage is Marée reconstruction output below; source bundles themselves are small. |
| `literature_models/Maree_2006/` | 13.65 | Standalone 2D Cellular Potts/reaction-diffusion reconstruction, source/specification and development output. Numerical and historical reconstruction limitations are documented in its [README](literature_models/Maree_2006/README.md); presence of outputs is not a claim of exact paper recovery. |
| `literature_models/Maree_2006/results/validation_01/` | 13.01 | Native development/validation campaign histories and spatial output. This is the dominant component of the Marée folder; smaller lambda/bias sweeps are covered by the parent. Not part of the M2072 retirement. |
| Other literature bundles, including `M7682/` and `debelly_2026_mechanochemical/` | covered by literature parent | Smaller original model/source/provenance bundles and selected demonstrations. The deBelly authors' implementation is a reference for current adaptations. M7682 contains a separate area-scaling campaign and does not belong to the retired M2072 branch. |

### HT1080 biosensor source data and analysis

| Path or grouped component | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `ht1080biosensors/` | **29.96** | Downloaded HT1080 Rac/Rho biosensor source data and derived segmentation, polarity, boundary-activity and trajectory-timescale analyses. The next two rows dominate this total. |
| `ht1080biosensors/GTPase4.rar` | 19.61 | Downloaded source archive; acquisition URL, expected size and MD5 are recorded in `download_GTPase4_aria2.slurm`. Compressed source payload, not a simulation cache. A recorded public URL does not by itself verify future availability. |
| `ht1080biosensors/sample_nd2/` | 9.33 | Selected extracted native ND2 movies across Rac and Rho recordings. Quantitative source channels used by segmentation/FRET/polarity extraction; overlap in content with the compressed archive is deliberate source/extracted representation, not an additional acquisition. |
| Remaining biosensor analysis folders | included in parent | Derived segmentation boundaries, extracted vectors and FRET data, QC/movie assets, fitted trajectory-timescale summaries and conditional-model simulations. Exact quantitative source channels cannot be replaced by contrast-normalized viewer PNGs. |

### LFCT data, reports and historical frameworks

| Path or grouped component | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `reports/` | **5.16** | Historical/working report sources, figures, native model experiment data and LFCT viewers. Includes the 4.48-GiB figure-revision group below; generated report assets are derived outputs. |
| `reports/figure_revision/` | 4.48 | Separate figure-only revision package: sources, assay manifest, control experiments and validation for the earlier report. Its `data/` is 4.45 GiB of full fields/case outcomes, including physical-peak controls, finite countercues, Spring and mechanochemical experiments; not covered by four-assay deletion. Follow its [README](reports/figure_revision/README.md). |
| `data/` | **4.12** | LFCT source images, TrackMate annotations and earlier segmentation/QC inputs. The MCF10A group below dominates; smaller MDA-MB test and 20× annotation groups are covered by this parent. |
| `data/lfct_mcf10a_fld102/` | 3.28 | MCF10A FLD_102 acquisition/source-manifest bundle, extracted phase/RFP images and earlier CPSAM segmentation variants. Image subfolder 1.97 GiB; downloaded source subfolder 1.24 GiB. Used by LFCT workflows; quantitative images and masks have distinct roles. |
| `migration1d/` | 1.62 | Deprecated moving-interval 1D framework retained for historical reproducibility. Almost all storage is generated single-/three-species campaign output under `results/` (1.59 GiB). New polarity work belongs in `polarity1d/`. |
| `lfct/` | 1.02 | Active compact resegmentation: compressed original annotations, durable cell/nuclear label stacks, assignment/QC tables, provenance and viewers. `results/` is 0.22 GiB; `reports/` is 0.66 GiB including delivery copies. Masks are quantitative evidence; viewer contours/PNGs are derived presentation. |
| `wgd_langevin/`, `current manuscript figs/` and `literature_review/` | below threshold | Scientifically important smaller empirical migration/OU analysis, working figure collection and source literature measurements. Start with [the figure README](<current manuscript figs/README.md>) and [refined OU source](wgd_langevin/refined_ou_20261003/README.md) for the WGD/enlargement motivation and interpretation limits. |

Other root directories are individually below the inventory threshold and are
primarily source, configuration, tools, Git metadata or small working assets.
This inventory has not independently revalidated every historical dataset.

## Round-2 follow-up retained outputs (22 October 2026)

Measured 2026-10-22T10:58:36.332652+00:00. These smaller groups are recorded for continuity; they do not update the older whole-workspace total above. Native inputs were not modified. No simulations were executed in this follow-up.

| Root-relative group | Allocated storage | Purpose and scientific status |
| --- | ---: | --- |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/` | 0.3937 GiB (422,731,776 allocated bytes) | Completed Figure 3/4 recompositions, post hoc/frozen H1 side-by-side, complete retained supplements, plotted CSVs, and Figure 2 in-progress review. |
| `polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022/` | 0.0206 GiB (22,069,248 allocated bytes) | Frozen source/checkpoint campaign for 18 Holmes deterministic trains plus L39 bounded-noise controls/replicates (32 workers); unsubmitted because scheduler socket access failed. |
| `manuscript_outline/revision_round2_20261022/` | 0.0284 GiB (30,474,240 allocated bytes) | In-progress review outline, validation, historical interim and recovery bundle; final outline_v5 integration deferred until Holmes measurements exist. |

## Round-2 completion retained outputs (22 October 2026)

Measured 2026-10-22T13:51:12.281884+00:00 using allocated blocks (`du -sB1`), after the scientific producers stopped. This updates the round-2 groups below; the older whole-workspace total remains an explicitly dated historical measurement. Original inputs were preserved and no scientific data were deleted.

| Root-relative group | Allocated storage | Purpose and scientific status |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022/` | 40.5760 GiB (43,568,185,344 allocated bytes) | Completed and independently audited native evidence: 18 Holmes L31–L39 four-pulse trains, 12 L39 bounded-input seeds and two controls. Exact checkpoints, frozen sources, all native phase fields and joined trajectories are scientifically authoritative; retain these raw outputs. |
| `polarity1d/manuscript/data/runs/signature_figure3_holmes_extensions_round2_20261022/` | 0.3641 GiB (390,971,392 allocated bytes) | Nine approved exact-endpoint Holmes cue-free continuations through 80,000; six resolutions and three lower bounds, with native fields, conservation and resolution audits. Required evidence for the revised decomposition. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/` | 0.3969 GiB (426,123,264 allocated bytes) | Completed derived Figure 2–4 artwork, complete retained supplements, plotted-value CSVs and scientific/display acceptance. Includes the HPC-only 164 MiB whole-seed bootstrap archive; do not mistake derived presentation for the underlying native trajectories. |
| `manuscript_outline/revision_round2_20261022/` | 0.4145 GiB (445,022,208 allocated bytes) | Final manuscript reader, main plates and legends, exact-source acceptance, archived review interim and historical recovery bundles. Final integration is complete; previous in-progress artifacts are retained history. |

## Figure 3 round-3 retained outputs (22 October 2026)

Measured 2026-10-22T19:50:17.492860+00:00 using allocated blocks (`du -x -B1 -s`), after the four native producers finished. These smaller, scientifically important groups are recorded for continuity; this does not replace the dated whole-workspace measurement above. No original scientific inputs or prior editions were removed. [Targeted storage record](reports/figure3_round3_storage_20261022.json).

| Root-relative group | Allocated storage | Purpose and scientific status |
| --- | ---: | --- |
| `polarity1d/manuscript/data/runs/figure3_round3_relaxation_v2_20261022` | 0.040138 GiB (43,098,112 allocated bytes) | Authoritative full native cue-on/cue-off fields for the four requested, fixed-parameter relaxation protocols; all stationarity stops, with iii withdrawing at 20,000. Frozen sources, exact basals, restart continuity and accepted-step invariant records retained. |
| `polarity1d/manuscript/data/runs/figure3_round3_relaxation_20261022` | 0.000572 GiB (614,400 allocated bytes) | Preserved import-only startup failure and frozen inputs; no integrations started and no scientific trajectories produced. Four successful protocols are in v2. |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3` | 0.003380 GiB (3,629,056 allocated bytes) | Derived round-3 main, publication legend, all 152 placement records, 120 native saved-field crossings/topologies, undefined-heading bands and 15-page supplement preserving every round-2 main/supplement page plus 90/99% heading-time pages. Scientific verification passed. |
| `manuscript_outline/revision_round3_figure3_20261022` | 0.001236 GiB (1,327,104 allocated bytes) | Two-page figure and legend, current main plates 2–4, offline reader and exact-section/scientific verification. Every other figure section preserves outline_v6. |
| `manuscript_outline/outline_v7.pdf` | 0.005630 GiB (6,045,696 allocated bytes) | Integrated outline_v7 PDF with Figure 3 round 3; all prior outline editions retained. |

## Figure 3 single-C layout retained outputs (22 October 2026)

Measured 2026-10-22T20:07:58.029705+00:00 using allocated blocks after report generation and verification. These small derived groups are recorded for continuity; no native campaigns were created, enlarged or removed, and the previous round-3 outputs remain intact. The older whole-workspace measurement is unchanged. [Targeted storage record](reports/figure3_single_C_layout_storage_20261022.json).

| Root-relative group | Allocated storage | Purpose and scientific status |
| --- | ---: | --- |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3_single_C` | 0.003983 GiB (4,276,224 allocated bytes) | Single-C landscape Figure 3: eight C=30 profiles at left, four-model/three-size B and C at right; exact copies of the audited source CSVs. The 17-page supplement retains the preceding four-placement main and legend, including C=47.6, and every round-2 page. No new simulations; scientific coordinates unchanged; layout/reuse/preservation checks passed. |
| `manuscript_outline/revision_round3_figure3_single_C_20261022` | 0.001205 GiB (1,294,336 allocated bytes) | Two-page figure and legend, revised main plates 2–4, reader and source/preservation verification. Prior round-3 and outline_v7 editions retained. |
| `manuscript_outline/outline_v8.pdf` | 0.005524 GiB (5,931,008 allocated bytes) | Integrated outline_v8 with the revised Figure 3 placement; all other figure sections and artwork preserve outline_v7. |

<!-- figure4-round3-compact-storage-20261022 begin -->

## Figure 4 round-3 compact campaign, 22 October 2026

Measured **2026-10-23T03:16:44.915711+00:00**: **28,756,893,696 allocated bytes
(26.781944 GiB)** under
`polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022/`. Status **COMPLETE**;
45,696 compact seeds complete, 24,960 initial
and at most 99,840 after precision-driven 96-seed additions. This targeted campaign
measurement does not replace the older whole-workspace totals above.
[Targeted storage and quota record](reports/figure4_round3_storage_20261022.json).

The scientifically authoritative new records retain post-burn-in float32
activity moments and P = M1/M0, contrast and exact round-1 tracked-region
descriptors, plus bounded float16 illustrations for seeds 0/1 only. Native
profiles exist in memory during extraction and are never written to durable
storage. Initial planning estimated 26.586 decimal GB, maximum 103.823 GB,
with a 207.647 GB reservation (2× margin); the shared filesystem exposed
34,251,346,018,304 free bytes at its dated probe.
Personal quota headroom could not be read and is explicitly unknown.
The initial small-QoS arrays were upgraded to xxlarge for independent-seed
throughput; `results/qos_burst_upgrade.json` records the effective scheduler
setting and capacity evidence, while the immutable manifest preserves the
original submission resources. The scoped old Figure4 noise native profiles were deleted after the COMPLETE audit; the separate retirement execution ledger records actual reclaimed storage.

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

<!-- figure4-scoring-revision-20261023 begin -->

## Figure 4 event-scoring revision, 23 October 2026

Measured **2026-10-23T13:05:36.887264+00:00** using allocated storage (`du -s -B1`). This reanalysis reuses the complete compact campaign above; it creates no new simulations and changes no retained simulation evidence.

| Path / output group | Allocated storage | Purpose and scientific status |
| --- | ---: | --- |
| `polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4` | 0.077888 GiB (83,632,128 allocated bytes) | Revised normalized-drift main and supplement, 1,300 event-based decomposition records, 520 reused detector-free comparisons, and 260 retained whole-seed bootstrap archives. Complete: 154 event-based closure ratios; 0/33 within 15% where multiple fronts occupy most frames. Scoring fractions and missing-angle bounds expose the conditional angle estimate. Preceding round-3 exports remain unchanged. |
| `manuscript_outline/revision_figure4_scoring_20261023` | 0.001755 GiB (1,884,160 allocated bytes) | Figure 4 with legend, combined main plates, reader, preservation and numerical/source checks. Complete; all other manuscript figure sections and artwork preserve outline_v9. |
| `manuscript_outline/outline_v10.pdf` | 0.007477 GiB (8,028,160 allocated bytes) | Integrated outline_v10, including the corrected Figure 4 interpretation and methods. Complete and validated. |

[Measurement receipt](reports/figure4_scoring_storage_20261023.json). All new output groups are below 2 GiB. The existing compact campaign was also remeasured on 23 October at **28,756,893,696 allocated bytes (26.781944 GiB)**, unchanged from its preceding storage receipt. It remains scientifically authoritative simulation evidence; this revision changes only its derived analysis and presentation.

<!-- figure4-scoring-revision-20261023 end -->
