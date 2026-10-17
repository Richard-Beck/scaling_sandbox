# Literature migration data assets

This collection supports replacing approximate publication readings with traceable measurements for the five studies in the current migration and OU analysis. Xuan 2018 retracking is the highest priority. Assets and extraction packages were generated on 3 October 2026. See [EXTRACTION_STATUS.md](EXTRACTION_STATUS.md) for products and remaining evidence limits. Downstream fitting has not been changed.

The [agreed OU fitting disposition](OU_FITTING_DISPOSITION.md) specifies upgraded Xuan MSD fits, refined estimates for existing measurements, and diagnostic checks using the Xuan movie trackplots, with the existing groups retained.

## Finding assets

Each study has its own directory: `xuan_2018`, `mallin_2023`, `jemaa_2017`, `jemaa_2023`, and `hosny_2026`.

| Directory | Contents |
| --- | --- |
| `source/` | Downloaded article snapshots, original supplements and archives. Mallin videos are publisher streaming renditions, explicitly identified in the manifest. |
| `figures/` | Downloaded paper figures for digitization. |
| `derived/` | Extracted text, embedded supplement images, sampled movie frames, and local viewing conversions. |
| `access/` | Link inventories, player snapshots and publication metadata. `failed_downloads/` contains HTML challenges, never usable media. |
| `extraction/` | Reproducible scripts, protocols, pixel and calibrated measurements, annotated overlays, QC and methods. |

[manifest.csv](manifest.csv) inventories source assets and extraction products with study, role, access status, source URL where recorded, derivation parent, byte count and SHA-256. Paths are relative to this directory. `source_url_precision` distinguishes exact URLs from article landing pages; a landing page is not an exact download locator. Modification times are filesystem metadata, not asserted acquisition times.

[studies.csv](studies.csv) summarizes access and limitations. [relocation.csv](relocation.csv) maps former temporary paths to current project paths. [audit/video_decode_validation.json](audit/video_decode_validation.json) records full decoding checks for the original Xuan video, both Mallin streaming videos and eight Jemaà 2017 videos. It does not establish compatibility with every viewer.

`audit/initial_records/` preserves the initial audit and correction notes unchanged. Those records contain former filenames and paths; the current manifest and this README govern the organized collection.

## Video sources and viewing

- **Xuan:** [original AVI](xuan_2018/source/xuan2018_movie.avi), 9,941,250 bytes; a fresh publisher download matched its SHA-256. It contains 74 frames at 1040 by 1392 pixels. [H.264 preview](xuan_2018/derived/xuan2018_movie_preview.mp4) is a lossy local viewing conversion. Retrack the AVI. The user's inspection identifies apparently one PACC/PGCC and multiple non-PACCs; verify identities during tracking.
- **Mallin:** [non-PACC movie](mallin_2023/source/mallin_movie1.mp4) and [PACC movie](mallin_2023/source/mallin_movie2.mp4) are H.264 streaming MP4s from the publisher's players. They are not the listed WMV originals, so their smaller sizes are expected. Original WMVs were not retrieved.
- **Jemaà 2017:** eight genuine original AVIs were extracted from the Europe PMC supplementary archive. Filenames retain archive member names: `s004` is Video 1 through `s011` as Video 8. The initial 1.8 KB `.avi` download was an HTML challenge and is now segregated under `access/failed_downloads/`.

Movie playback time is not experimental elapsed time. Xuan's supplement labels the movie 16 hours while the main migration methods describe 12 hours; that discrepancy remains open.

## Evidence and access limits

No public numerical trajectory dataset was identified for these five comparisons in this audit. Xuan and Mallin movies are candidates for recovering representative tracks, not established substitutes for complete experimental replicates. Jemaà 2017 movies do not cover the key diploid–tetraploid tracking comparison. Jemaà 2023's supplementary ZIP contains a two-page PDF, without movies or trajectory tables.

Hosny's preprint PDF, main figures and S1–S6 supplement are available here. The September 2026 journal version was identified, but its full text and supplements were not retrieved. Elsevier XML contains metadata only. A related Zenodo deposit lists oxygen-intensity images and code; no cell trajectories were identified there, and that dataset was not downloaded.

## Extraction and retention

Follow [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) before producing measurements. Keep source files immutable and record new derivatives and outputs with their source hashes. Existing literature summaries remain in the parent folder; current OU observations and fits have not been changed.

Downloaded and derived payloads are locally ignored by this directory's `.gitignore`; manifests, the extraction plan and compact extraction outputs remain eligible for later Git review. Local retention is not a verified external backup. Do not delete originals on the assumption that Git or a preview preserves them.
