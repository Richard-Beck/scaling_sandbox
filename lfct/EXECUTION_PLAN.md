# Staged execution plan

Updated 2026-09-09. Core pipeline is implemented under `lfct/`; the smoke scope passed. The user subsequently accepted HEK293 limitations and checked U87 alignment, authorizing direct full processing. See [README](README.md) for storage layout and commands. Full field array **20236390** and dependent audit/report **20236391** are now submitted; the qualification/rehearsal stages below are superseded by this user decision. See `data/qc_decisions.json` and `results/full/submission.json`.

## Smoke scope and evidence

Four fields, **21 full-resolution frame pairs / 42 channel images**:

| Cell line / field | XML frames (zero based) | Cell assignments | Nuclear assignments |
|---|---|---:|---:|
| HEK293 / FLD_2 | 1, 115, 228 | 71/134 (53%) | 79/134 (59%) |
| MCF10A / FLD_102 | 1–12 | 46/48 (96%) | 22/48 (46%) |
| MDA-MB / FLD_1 | 1, 192, 383 | 60/71 (85%) | 64/71 (90%) |
| U87 / FLD_1 | 1, 144, 286 | 63/90 (70%) | 46/90 (51%) |

These are assignment coverage, **not segmentation accuracy**. Shared masks are retained and flagged rather than split arbitrarily. Weak RFP causes missing nuclear masks in early MCF10A. HEK293 has substantial underassignment and shared masks; qualify segmentation before applying isolation broadly. No 10-frame isolated run was found in the MCF10A test window (longest eligible stretch: 6). The other three sampled windows are nonconsecutive and cannot establish 10-frame isolation. The targeted clip is explicitly a negative control. A synthetic positive 10-frame run passes the tests.

Completed checks include original IDs/positions and detection counts, one-to-one label mappings, exact compressed-mask pixel round trips and TIFF frame labels, all-mask neighbor distances, gaps/branch separation, uncertain neighbors, missing nuclei, field edges, source mapping exceptions, and safe resume checks. Unit tests: 16. Browser checks exercise local PNG loading, frame slider/wheel/keyboard, zoom, layers and every field/clip selector. Machine-readable measurements live in `results/smoke/benchmark.csv`, field provenance and `report_size.json`.

Source issues encountered:

- Conda's urllib certificate chain failed on Box; secure system `curl` with normal certificate verification fixes downloads.
- U87 FLD_1 contains 289 source images against 288 XML frames. RFP suffix 1 is anomalous. Visual alignment at XML 0,1,144,287 supports skipping source suffix 1 in both channels. This is a **provisional explicit mapping assumption**, not dataset-author confirmation. Inspect adjacent phase/RFP frames during qualification; evidence is in `results/smoke/alignment/`. Any other unmatched index range stops the pipeline with a failure record.
- Segmentation-setting/source-mapping changes are refused for existing mask outputs: use a new run name. Analysis-only changes can reuse verified masks. Earlier failed-attempt logs remain for audit; the latest `completed.json` is the completion authority.

## Measured runtime and interpretation

CPSAM 4.0.6, PyTorch 2.7.1+cu126, NVIDIA A30 on `hpctpa3pc0064`, four CPU threads. The first successful inference pass took about **2.1–2.4 seconds per phase/RFP pair** on the first three fields. Initial end-to-end times: HEK293 29.6 s for 3 pairs; MCF10A 36.3 s for 12 pairs using existing source TIFFs; MDA-MB 40.6 s for 3 pairs including fresh source downloads. CPU analysis/QC took 1.6–2.7 s per field sample. Peak batch RSS in that allocation was approximately 1.9 GiB. U87 took 28.9 s for three freshly segmented pairs, including 8.1 s preparation and 4.4 s analysis/QC, with median inference 2.41 s/pair. Its process peak RSS was 2.17 GiB. Detailed per-field measurements are in the benchmark CSV.

The final rerun reused the first three fields' masks, so its wall times measure resume/reanalysis rather than fresh inference. Inference timings remain attached to the original segmentation job. The first allocation lasted 1:58 including the subsequent U87 failure; do not divide its duration by successful frames as a throughput estimate.

At 2.1–2.4 s/pair, 8,310 pairs imply **4.8–5.5 GPU-hours for inference alone**. Budget **6–10 allocated GPU-hours** for the full baseline, allowing startup, downloads, compression, denser scenes and analysis. This is a small-sample extrapolation, not a full-stack benchmark. Download volume is largely independent of retained frame count because each channel is supplied as a field ZIP. Queue delay and shared network contention are excluded from all estimates. Model variants add approximately proportional inference cost.

## Stages and gates

| Stage | Scope | Expected execution, excluding queue/QC time | Submission recommendation |
|---|---|---|---|
| 0: smoke — complete | 4 fields, 21 pairs | A few minutes; measured per-field values above | One serial allocation, A30 ×1, 4 CPUs, 16G, 30 min |
| 1: qualification | Same 4 fields, 30 pairs each: early/middle/late 10-frame windows | 2–5 min/field; 10–20 allocated GPU-min baseline | Four field tasks; A30 ×1/task, 4 CPUs, 16G, 20 min/task |
| 2: full-field rehearsal | Same 4 fields, all 1,288 pairs | 15–30 min/field; about 1–2 allocated GPU-hours | Four field tasks; A30 ×1/task, 4 CPUs, 16G, 60 min/task |
| 3: remaining fields | Other 21 fields, 7,022 pairs | About 5–8 allocated GPU-hours | One task per field, same resources and 60 min limit initially |
| 4: final report and audit | All 25 field outputs; 50 overview pairs, ≤20 targeted clips | Budget 5–15 CPU-min | 1–2 CPUs, 4G, 15 min; no GPU |

Stage 1 gate: visually check alignment and outlines in every cell line/time window, particularly U87 adjacent-frame alignment and HEK293 merged/missing cells. If needed compare a small number of phase segmentation settings in separate run names before selecting fixed settings. The current implementation uses independent phase/RFP segmentation; combined-channel experiments would require an explicit implementation change. Do not relax isolation thresholds to hide missing neighbors. Manual QC remains necessary; assignment coverage alone cannot approve segmentation.

Stage 2 gate: inspect candidate clips plus exclusion/near-threshold controls; confirm actual full-field runtime, peak RAM/scratch and compressed bytes per frame. Check all frame mappings and source identifiers, ensure scratch cleanup and successful resume, and revise estimates if dense fields exceed expectations. Expand the bounded QC selection logic if negative/near-threshold clips are needed in addition to the currently selected one clip per field.

Stage 3: submit only after those gates pass. Failures remain field-local. Re-run the same command after fixing recoverable download issues; changed scientific inputs require a named replacement run. Never merge masks from different configurations implicitly.

Stage 4: rebuild the external-asset HTML and ZIP; verify every asset and stack checksum and compare aggregate counts with the catalog. The current full profile chooses two overview frames per field, using a fixed seed across distinct temporal thirds (50 overall). It saves at most one provisional targeted clip per field and publishes at most 20; selection can be refined for cell-line balance before final delivery. Keep one canonical final run; remove pilot duplicates only after the final output is validated and cleanup is agreed. Existing exploratory data/reports are untouched.

## Concrete SLURM commands

Run from `/home/4473331/projects/scaling_sandbox`. `red` / QoS `small` was inspected and used successfully. All GPU scripts exclude the four known CUDA-broken nodes. No array concurrency throttle is applied; each task uses one GPU and local scratch. A future throttle should be tied to a measured storage/network/scheduler limit.

```bash
# Stage 1 (not yet submitted)
sbatch --array=0-3 --time=00:20:00 lfct/workflow/run_field.slurm --profile pilot
# After all four succeed:
/home/4473331/.conda/envs/cpose/bin/python lfct/scripts/report.py --profile pilot --zip

# Stage 2: submit one full job for each explicitly selected field.
for field in hek293/FLD_2 mcf10a/FLD_102 mdamb/FLD_1 u87/FLD_1; do
  sbatch lfct/workflow/run_field.slurm --profile full --field "$field"
done

# Stage 3: full catalog array. The four already-completed fields validate and exit.
sbatch --array=0-24 lfct/workflow/run_field.slurm --profile full
# Stage 4 after completion:
/home/4473331/.conda/envs/cpose/bin/python lfct/scripts/report.py --profile full --zip
```

The final two report commands are light CPU work; use a CPU batch allocation if measurements exceed the repository's login-node limits. Keeping download, inference and analysis in one field job simplifies recovery; the GPU remains allocated during download/analysis, which is included in the allocated-time budget. Do not interpret 6–10 GPU-hours as elapsed time: four available GPUs could ideally finish that work in roughly 1.5–2.5 h, but queue delays, task lengths and download contention determine reality.

## Storage

Retain every lossless cell/nuclear mask and audit table, plus compressed original XML. Keep only selected display PNGs. Measured smoke storage: HTML **0.91 MiB**, 48 PNGs **18.16 MiB**, and all 42 label-mask pages **0.26 MiB**. Including four compressed source XMLs, result tables, diagnostics and screenshot, the smoke occupies approximately **49 MiB logical / 67 MiB allocated**, excluding the optional approximately 18.4 MiB delivery ZIP. Both the original folder and independently extracted ZIP passed browser checks with all 48 images loading. All four completed fields also passed a repeat checksum-only invocation without downloads or inference.

Initial durable allowance remains **0.5–1.6 GiB**, with a 2 GiB target and 3 GiB headroom, pending full-field measurements. ZIP delivery copies are optional duplicates. Budget **5–10 GiB transient scratch per active field** (archives plus extracted images); processing concurrent fields multiplies this. Node `/tmp` was used successfully but does not provide recovery across job allocations. Mask stacks and completion records live in the durable repository workstream. Do not count existing large exploratory inputs/viewers as new pipeline storage or delete them before validated migration.
