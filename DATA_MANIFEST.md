# Project data and output storage manifest

This is the high-level inventory of large **local** data and generated outputs
in `scaling_sandbox`. Paths are relative to the repository root. The inventory
groups related files under meaningful campaign or dataset parents; it is not a
file-by-file catalogue or a statement that a dataset may be deleted.

**Current post-cleanup workspace: 178.79 GiB (191.97 GB)**, measured on
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

This table comes from one complete `du -x -B1 --max-depth=6 .` traversal during
20:59:40–21:00:35 UTC on 2026-10-17, counting each hard-linked inode once across
the whole workspace. Native/statistical producers and export compression were
terminal; small final assembly records could still change. Detailed rows below
carry their own measurement dates; they are components of these root totals,
not extra storage to add to the overview.

| Root-relative area | Allocated GiB | Main purpose |
| --- | ---: | --- |
| `polarity1d/` | **121.80** | Active corrected manuscript evidence and development experiments, independent historical archives, native qualification/reproduction data, and the compact retired four-assay record. |
| `ht1080biosensors/` | **29.96** | Biosensor source archive, extracted native movies and derived analyses. |
| `literature_models/` | **13.97** | Publication-specific model implementations and reproduction data, dominated by the Marée reconstruction; M2072 now occupies only 0.91 MiB. |
| `reports/` | **5.16** | Figure-revision model experiment data, report sources/assets and earlier LFCT viewers. |
| `data/` | **4.12** | LFCT source images, annotations and early segmentation variants. |
| `migration1d/` | **1.62** | Historical moving-interval framework and generated campaign results. |
| `lfct/` | **1.02** | Current segmentation labels, annotations, assignment/QC tables and report/delivery copies. |
| Remaining root areas and root files | **1.13** | Smaller source, literature/OU analyses, working figures, Git metadata and temporary assets. |
| **Whole workspace** | **178.79** | **191.97 decimal GB**; this is local allocated storage, not a personal quota report or an external-backup inventory. |

## Retirements authorized on 2026-10-17

| Path | Former allocated storage | Current record and scientific consequence |
| --- | ---: | --- |
| `literature_models/M2072/` | 86.66 GiB | **Retirement complete: 0.91 MiB retained**, 19 files and 472,183 logical bytes. Compact historical notes, original V1–V6/selected leader parameters, two selected XML models and provenance survive in [the retired branch](literature_models/M2072/README.md). Bulk campaigns and development pipeline were deleted; this is not an external archive. Measured allocated-footprint reduction: 93,049,708,544 bytes. See [the completed action record](literature_models/M2072/retirement_20261017/cleanup_status.json). |
| `polarity1d/four_assay/` | 127.02 GiB | **Cleanup and post-cleanup checks complete: 0.88 GiB retained** (944,349,184 allocated bytes), a measured footprint reduction of 135,440,498,688 bytes (126.14 GiB). Compact metadata in 24 verified archives, scalar records, source/protocol provenance, accepted v8 material, final v11 profiles and selected illustration inputs survive. All 5,006 canonical selected-input hashes passed again; data-only regeneration reproduced the current main and width-supplement PNGs byte for byte. Full native coverage of old trials is intentionally retired. The active corrected manuscript workflow is separate. See [the execution/retention record](polarity1d/four_assay/OUTPUT_RETIREMENT_20261017.md), [human cleanup report](reports/storage_cleanup_four_assay_20261017.md) and [machine-readable ledger](reports/storage_cleanup_four_assay_20261017.json). |

## Inventory

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
