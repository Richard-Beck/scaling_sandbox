# LFCT: resegmentation, visual QC, and isolated trajectories

Planning scope: the existing **25 fields at 20×**, across four cell lines, comprising **8,310 frames / 16,620 channel images**. Expansion to other magnifications would be a separate addition. This document proposes the workstream; no existing results or input images have been moved or deleted.

## Recommended design

Process one field at a time through download → segmentation → centroid assignment → isolation analysis → QC selection → compact output validation → temporary-file cleanup. The durable scientific outputs are lossless label masks and small audit tables. Raw image ZIPs and full image sequences are temporary inputs. Reports contain a small selected set of visualizations, referencing PNG files rather than embedding full movies.

Use the existing `cpose` environment and cached CPSAM weights. GPU work runs through SLURM, with the established CUDA node exclusions. Start with one field per GPU job, four CPU cores, 16 GiB RAM, and QoS `small`. Add concurrency only when useful; if a scratch-space budget requires a worker/concurrency limit, record that specific justification. Avoid separate scheduler jobs or permanent output files for every frame.

## Durable layout

Keep the workstream together under `scaling_sandbox/lfct/`, consistent with the repository's other self-contained workstreams:

```text
lfct/
  README.md                      # entry point, run instructions, current outputs
  PLAN.md                        # this plan
  config.json                    # segmentation, matching, isolation, QC settings
  scripts/
    prepare.py                   # source manifest, downloads, extraction
    segment.py                   # CPSAM cells/nuclei and compact mask writing
    analyze.py                   # assignments, branch runs, neighbor clearances
    report.py                    # small overview and targeted-QC viewer
  workflow/
    run_field.slurm              # field-level pipeline; resumes completed stages
  tests/                         # assignment, isolation, and browser checks
  data/
    sources.csv                  # cell line, field, URLs, frame mapping, hashes
    annotations/<cell_line>/*.xml.gz
  results/
    fields/<cell_line>/<field>/
      cells.tif                  # losslessly compressed, multipage label stack
      nuclei.tif                 # same frame order; independent nuclear labels
      detections.csv.gz          # IDs, positions, assignments, distances, flags
      masks.csv.gz               # every mask, including unmatched objects
      runs.csv.gz                # candidate uninterrupted isolated runs
      provenance.json            # settings, source/model/code hashes, stage status
    summary.csv
    qc_manifest.csv              # selected frames/clips, reasons, coordinates
  reports/
    index.html                   # small report; metadata/selected contours inline
    assets/
      overview/                  # 50 selected full-field phase/RFP PNG pairs
      tracks/                    # selected cropped phase/RFP clip PNG pairs
```

Use a cluster-approved scratch location outside the durable tree for ZIPs, extracted full sequences, temporary per-frame masks, and disposable logs. Do not duplicate model weights. Keep code, configuration, manifests, and documentation in version control; ignore generated masks, PNGs, downloaded annotations, and scratch files. Keep one canonical current report and its assets, without accumulating copies of prior full-stack viewers.

Retain source XML losslessly as `.xml.gz`, preserving its exact original information. Preserve Spot IDs, Track IDs, track names, positions, and links. Frame-local segmentation labels remain separate identifiers mapped through the tables.

## 1. Establish a compact pilot before scaling

Refactor the existing FLD_102 implementation into this layout rather than rewriting the working segmentation/viewer. First repack its existing masks into two compressed TIFF stacks and reproduce its assignment results from those stacks. This requires no segmentation rerun.

Then inspect early, middle, and late frames from one field of each cell line. Confirm source/channel/time alignment and that CPSAM settings transfer adequately. Start with the current phase settings and local RFP normalization, but allow a documented setting per cell line if the pilot supports it. Fix settings before full processing; avoid ad hoc frame-by-frame tuning. Do not assume weak fluorescence or original polygons are reliable nuclear/whole-cell ground truth.

## 2. Rerun segmentation and preserve compact scientific outputs

- Download both channels for each field and verify frame indices, dimensions, and checksums against the source manifest and XML.
- Run CPSAM independently on phase for cell masks and RFP for nuclear masks, retaining the original tracking graph.
- Reuse conservative one-to-one centroid assignment, preserving unmatched centroids and masks with explicit statuses. Keep a reproducible distinction between containing, nearby, nuclear-overlap fallback, and ambiguous assignments.
- Save **all cell and nuclear label masks**, not just masks assigned to tracks. Unassigned cells matter for contact detection.
- Store masks as two losslessly compressed uint16 multipage TIFF stacks per field. Keep integer labels exact. Avoid a parallel permanent representation as all-frame boundary JSON/SVG/PNG: derive boundaries from these masks when needed, and store only the small selected-viewer geometry.
- Consolidate per-frame JSON and metrics into field-level compressed tables and provenance. Keep concise failure records rather than large redundant logs.
- Validate mask counts, frame alignment, original identifiers, label mappings, and file readability before removing temporary raw images and intermediates. Write a completed-field marker only after validation. A failed field retains its recoverable working data.

Field processing should be resumable. A parameter change creates an explicitly identified replacement run, not silently mixed old/new masks.

## 3. Define and find isolated trajectory runs

The analysis unit is an **uninterrupted, nonbranching trajectory**, not a whole TrackMate connected-component ID. Break chains at splits, merges, missing temporal links, and frame gaps. Give each derived chain/run its own identifier while retaining original IDs and names.

A starting definition is **at least 10 consecutive valid image frames** with every cell boundary at least **10 pixels** from every other detected cell boundary. Make both values configurable. Also retain each frame's actual minimum clearance so results can be filtered at 0, 5, 10, or 20 pixels without rerunning segmentation. Until calibration is verified, report pixels and frames rather than inventing micrometers or elapsed minutes.

For each frame:

1. Measure cell-mask boundary clearance against **all other cell masks**, including unassigned masks and cells sharing the same TrackMate ID. Bounding boxes may accelerate candidate search but must not determine final contact distances.
2. Require a uniquely supported cell assignment. Shared/ambiguous masks or an unmatched centroid break a valid run. Record non-containing matches separately; initially exclude them from the strict primary candidate set.
3. Do not interpret missing segmentations as empty space. Nearby unmatched tracked detections or nuclear masks without a corresponding cell mask flag the local neighborhood as uncertain. Use a documented conservative guard distance, informed by observed cell sizes, and expose it in configuration.
4. Exclude cells too close to the image boundary to verify the requested clearance. A missing nuclear mask alone does not disqualify an otherwise well-supported phase-cell isolation result, but remains a QC flag.
5. Calculate contiguous qualifying runs; report their length, minimum/median clearance, assignment quality, and any exclusion reasons. Do not bridge invalid frames to reach the 10-frame minimum.

Call these **isolated at every observed frame**. The available imaging intervals cannot establish that cells never approached or contacted between observations. Automated segmentation errors can also hide neighbors, which is why targeted QC remains necessary.

## 4. Retain a small, deliberate visual-QC set

**Overview: 50 full-field frames total**, initially two per field. Use a fixed seed and balance selected temporal thirds across fields and cell lines. Sample independently of segmentation success so failures are not quietly filtered out. Keep the exact selection in `qc_manifest.csv`.

For those frames, retain phase and RFP as full-resolution, contrast-normalized **8-bit PNGs**, recording the normalization bounds. These are visual-QC images, not replacements for the original quantitative fluorescence measurements. No lossy JPEG is needed for this smaller selection.

**Targeted QC: initially up to 20 short clips**, approximately 10–14 consecutive frames each, distributed across cell lines. Include promising isolated runs and cases close to the contact threshold or with relevant assignment failures. Keep fixed spatial crops over each clip, with enough surrounding context to show nearby cells; store the crop offset in original image coordinates. Use full-field context when a tight crop would conceal a relevant neighbor.

Select and save provisional targeted crops while the field's raw images still exist, with a bounded per-field allowance. The final overview chooses the desired clips from this small pool. Unneeded provisional crops can be removed after report validation. Later QC requests can re-download only the relevant field, rather than requiring permanent retention of every raw image.

Use one report with an overview selector and a targeted-track selector. Preserve the useful current controls: phase/RFP/plain background, cell/nuclear/original outlines, tracked centroids, zoom/pan, frame slider, original-ID inspection, and issue flags. Include clearance versus time and assignment status for each targeted clip.

The HTML embeds only small metadata/selected contour data and loads PNGs from `assets/` on demand. It should open locally without a server; do not rely on fetching external JSON through `file://`. Share the report and assets as **one ZIP**, rather than embedding the PNGs as base64. Keep one canonical PNG copy; the delivery ZIP is disposable/rebuildable, not another permanent dataset. Do not export SVGs for all frames or repeat raw backgrounds inside other reports.

## Storage budget

Measurements from the current MCF10A field:

- Downloaded ZIPs: 1,129 MiB; extracted TIFFs: 1,629 MiB.
- Exact compressed cell/nuclear masks for all 385 frames: **8.8 MiB** combined, excluding redundant per-frame metadata.
- Source XML: 73.3 MiB uncompressed, approximately **8.5 MiB compressed**.
- A representative phase/RFP PNG pair: **1.53 MiB** as 8-bit display images, versus 2.8 MiB as original-valued 16-bit PNGs.
- All 25 current XMLs: 953 MiB uncompressed. Filesystem allocation differs from logical size because of server compression and per-file overhead; consolidate files and monitor both.

Initial persistent planning allowance for all 25 fields:

| Item | Allowance |
|---|---:|
| Compressed original annotations and source manifest | 0.1–0.25 GiB |
| All cell/nuclear mask stacks | 0.2–0.6 GiB |
| Assignments, masks, runs, provenance | 0.05–0.2 GiB |
| 50 overview PNG pairs | 0.08–0.15 GiB |
| Bounded targeted-QC crops | 0.1–0.3 GiB |
| Report, code, small logs, miscellaneous | under 0.1 GiB |

**Aim for roughly 0.5–1.6 GiB of logical durable data, with a 2 GiB planning target and about 3 GiB of headroom for filesystem allocation/denser fields.** These are projections, not measured all-dataset totals. Measure the first field from each cell line before committing the rest. If usage grows unexpectedly, reduce duplicate/QC storage before considering any loss of the all-frame masks.

Transient space can be budgeted separately, initially around **5–10 GiB per active field worker**. Free the downloaded ZIPs, extracted full sequences, temporary masks, unselected image previews, and rejected pilot outputs after durable-output validation. This scratch budget can govern field concurrency without changing the scientific outputs.

## Implementation order and cleanup

1. Create the compact FLD_102 pilot, small overview report, and durable layout; verify identical IDs and mask content after repacking.
2. Add isolation metrics and test their behavior on touching cells, separated cells, shared IDs, segmentation gaps, ambiguous neighbors, and field edges. Inspect candidate clips before scaling.
3. Check transfer on the other three cell lines; finalize per-line settings and storage estimates.
4. Run the remaining fields through the resumable field pipeline; assemble the 50-frame overview and bounded targeted-QC collection.
5. Verify report links, all-frame scientific outputs, hashes, ID preservation, and QC records. Then replace the current scattered exploratory LFCT outputs with the compact canonical workstream and remove superseded raw copies, full-stack embedded reports, and unused pilot variants.

The existing large viewer and source files remain available until a validated compact replacement is ready. No migration, deletion, or all-dataset compute is authorized merely by writing this plan.
