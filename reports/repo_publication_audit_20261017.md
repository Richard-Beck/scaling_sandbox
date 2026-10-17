# Repository publication audit — 17 October 2026

This audit supports the authorized source/manuscript commit and push. It does
not authorize further data deletion. The root agent owns the Git index, final
commit, push and rendered manuscript audit. This audit changed ignore/retention
policy and repaired the three specifically reviewed relocated literature bundles;
it did not stage, commit, push, delete data, reset files, or alter the old website.

## Selected scope and index boundary

The original index contained only `.gitignore`, `README.md` and `RETENTION.md`.
The working tree contained 112 preexisting tracked deletions alongside source
changes and extensive untracked research material. Broad `git add` was unsafe:
Git initially exposed about 15.6 GB of logical untracked content, including
60,728 historical polarity archive files (10.11 GB) and 1,348 figure-revision
files (4.20 GB). Eligibility alone did not identify publication evidence.

[`publication_files_20261017.txt`](publication_files_20261017.txt) is an explicit
candidate selector for the reviewed root/workstream scope. It selects source,
configs, tests, licences, compact retirement/provenance records, the current
figure collection, reviewed literature extraction tables, refined OU evidence,
and the previously curated WP1 figure inputs. It excludes the separately owned
`polarity1d/manuscript` selector, old trial-report PDFs, unselected diagnostic
reports and older uncurated WGD branches. Its 1,429 entries at the final source
check included 1,397 existing files totalling 259,656,490 logical bytes, 31
paired deletion paths, and this report, which had not yet been written. Counts
will change when the final manuscript selector is combined and reviewed.

Only the 31 verified paired relocations are deletion candidates in that selector.
[`publication_relocations_20261017.json`](publication_relocations_20261017.json)
classifies every preexisting tracked deletion: 31 exact relocations, 14 modified
intended replacements and 67 paths without verified replacements. Exact identity
was checked against `HEAD` before the separate portability repairs described
below; SHA-256 identities are recorded. The remaining 81 old-path deletions are
excluded. On a published checkout, those originals therefore remain present
from `HEAD`, even though their preexisting working-tree deletions remain locally
unstaged. In particular, no website/Pages workflow retirement is selected.

The 14 modified replacements were examined in full diffs, retained in
[`publication_modified_sources_20261017.json`](publication_modified_sources_20261017.json).
They include a substantially reimplemented Julia axisymmetric example/project,
changed Julia launcher defaults and optogenetic hooks, expanded Taniguchi
allocation-free kernels and tests, and additions to the Chiou abundance adapter.
These are preexisting developments, not certified byte-identical moves. None of
their old-path deletions is selected, so unique original implementation details
are preserved in the published tree. The new sources are added separately.

Combine the root selector with the **final**
`manuscript_checkout_files_20261017.txt` and final outline/builder assets only
after the coordinated producer has completed its checkpoint. Stage using a
NUL-delimited pathspec file derived from those selectors, preserving filenames
with spaces. Inspect `git diff --cached --name-status`, all selected deletions,
`git diff --cached --check`, total staged bytes and largest individual objects
before committing. Do not replace this with a blanket force-add.

## Ignore and retention repairs

The root `.gitignore` now protects historical polarity archive payloads;
four-assay work/results/reference trees and compacted worker archives;
`reports/figure_revision/data`; unselected historical report figures; earlier
OU/mechanochemical experiment branches; scratch outputs; downloaded Morpheus and
compiled model binaries; and model-local native report outputs. Markdown history,
the approved assay-explanation keeper and selected Figure 3/Hosny assets remain
eligible. Obsolete M2072 campaign-promotion exceptions were removed after that
branch's separately authorized retirement. `RETENTION.md` explains the current
publication boundary and supersedes its earlier three-file staging limit.

`git check-ignore` confirmed protection for representative historical NPZ,
report native data, compacted trial archives, scratch files and both downloaded
binaries. It also confirmed the narrow current helper/keeper exceptions.
`git check-ignore --no-index --stdin` reported **zero ignored existing files**
in the explicit root candidate selector. No force-add is needed for this scope.

## Repaired reproduction entry points

The Giese and Taniguchi Python drivers now import their packages from bundle-local
`src/`. Shell submitters enter their own bundle and submit scripts at their new
locations. Their SLURM launchers accept submission from the repository root or
model bundle. Test discovery now resolves the relocated packages. Output paths
retain their model-local relative layout and are ignored. The README files give
current entry points and remove the false assertion that a `reorg/legacy` backup
exists. No new Conda environment or dependency was installed.

The Chiou Python drivers now resolve the existing shared `migration1d/src` or
Giese `src/bulk_surface` dependency. R scripts resolve their relocated bundle
under `Rscript` or `source()`/agent-runner invocation, and source its `src` and
`params-or-config` paths. The published MAT default and report source reference
use the retained `params-or-config/initial_state_1D.mat`. The missing original
49-line Chiou scenario was recovered from `HEAD` into that configuration folder;
its only adaptation changes the first source path. Original Git blob and SHA-256
provenance are in
[`params-or-config/RECOVERY_20261017.md`](../literature_models/chiou_2015_mcas/params-or-config/RECOVERY_20261017.md).
The root-level preexisting deletion remains unselected.

Focused verification passed:

- Existing Giese/bulk-surface and Taniguchi test suites: **41 passed in 17.93 s**.
- Bash syntax checks for all repaired shell/SLURM files: no failures.
- CLI imports/help for Giese native, Taniguchi and Chiou bulk-surface drivers: passed.
- Nine Chiou R files parsed; source imports, original MAT input and the original
  18-condition design passed in the existing `R-4.2.3` environment.
- `check_chiou_port.R` executed from `/tmp` using its absolute path and passed,
  including the short conservative numerical smoke check.
- `git diff --check` passed before final manuscript staging.

The current refined OU renderer imports
`wgd_langevin/persistence_only/campaign.py`, which itself requires
`wgd_langevin/hosny/data.py` and `wgd_langevin/data/observations.csv`.
These compact dependencies are selected through a narrow ignore exception;
older persistence-only campaign outputs remain excluded. Frozen extraction
inputs referenced by the refined OU preparation/validation scripts are eligible.
Full native literature reproduction campaigns were not rerun by this audit.
Other model bundles remain provenance resources with their documented limitations;
this report does not certify every historical launcher in the repository.

## Storage inventory and live production

A complete `du -x -B1 --max-depth=2 .` traversal completed around 20:30 UTC on
2026-10-17. It counted hard-linked inodes once across the workspace and did not
follow directory symlinks. The scan was a measurement interval while new
manuscript work was running, rather than an atomic filesystem snapshot.

| Area | Allocated bytes | Allocated GiB |
| --- | ---: | ---: |
| Whole workspace | 188,084,649,984 | 175.167 |
| `polarity1d/` | 127,301,320,704 | 118.559 |
| `polarity1d/manuscript/` | 101,551,104,000 | 94.577 |
| `polarity1d/archive/` | 15,152,185,344 | 14.112 |
| `polarity1d/results/` | 9,624,829,952 | 8.964 |
| `polarity1d/four_assay/` | 944,349,184 | 0.879 |
| `ht1080biosensors/` | 32,168,869,888 | 29.960 |
| `literature_models/` | 14,990,311,424 | 13.961 |
| `literature_models/Maree_2006/` | 14,652,014,592 | 13.646 |
| `reports/` | 5,538,316,288 | 5.158 |
| `reports/figure_revision/` | 4,810,629,120 | 4.480 |
| `data/` | 4,426,620,928 | 4.123 |

The root `DATA_MANIFEST.md` covers every approximately 2-GiB top-level output
group found in this traversal, with purpose and scientific status. It also
identifies the large native manuscript campaigns, archived source histories,
HT1080 archive/movies, Marée validation output and LFCT image group. M2072
retirement and four-assay cleanup records remain explicit and independently
referenced. The audit created only small source/report files and did not change
large-output retention.

**Required final refresh:** the manifest's earlier 168.57-GiB overview and
91.66-GiB manuscript measurement predate live production. After the coordinated
producer is terminal, remeasure the changed manuscript run groups and one
whole-workspace traversal, update the manifest's measurement time, current total,
root overview and affected component rows, and then stage the refreshed file.
Do not add separate subtree `du` values to estimate global storage or reclaimed
space. Root and component observations differ when hard links are counted in a
different traversal order.

## Publication checks and limitations

A redacted pattern scan of the existing selected text files below 2 MiB found
no private-key markers, GitHub/OpenAI/Google/AWS key-shaped strings, signed-cloud
URL signatures or credential filenames. The selected small files contained no
ELF executable headers. This is a bounded check, not a guarantee that all local
research data or every possible credential format has been inspected. The
ignored native/downloaded payloads were not scanned into chat output.

The initial explicit selector's largest individual source-evidence objects are
compressed curated NPZ files below 44 MB, not bulk native campaign histories.
The coordinated manuscript selector and newly generated outline must receive
the same final size, dependency and rendered-artifact review. Local raw data and
native histories remain outside Git; no external backup was verified here.

## Retained Pages workflow compatibility follow-up

The existing `HEAD` Pages workflow remains selected through the unchanged index
entry; its preexisting working-tree deletion is unstaged. It automatically runs
on a `main` push and executes the original simulation-path check, plot-data
exporter and report renderer. The exporter still has all three root report
sources and the committed WGD snapshot/analysis dependency in the staged tree.
The 31 paired relocations would, however, remove the original Holmes and
wave-pinning reports from the renderer's directory glob and its public output.
The newly selected LFCT Rmd would enter that glob despite requiring local
annotation tables/assets absent from a source checkout.

A small compatibility repair adds `scripts/public_report_sources.R`: it retains
the original nine-report catalogue scope and resolves missing old Holmes/WP
paths to their canonical literature-bundle sources. The existing safety checker
and renderer share this source list. Rendering uses the original `reports/`
knitting directory, preserving old asset and `../docs/data` export paths, and
an absolute `docs/reports` output directory. The workflow, site landing page,
scientific report contents and deployment behavior are preserved.

A separate 73-file site-dependency snapshot was materialized from the current
Git index under `/tmp/scaling-site-preflight-20261017-ahpf7_g6`, then overlaid
with these three proposed script changes. It contains the retained index versions
of preexisting working-tree deletions; no original worktree files were restored
and the index was not modified. In that snapshot, the explicit nine-source
selection, both canonical fallbacks, LFCT exclusion, simulation-path safety scan
and renderer parse all passed. `git diff --check` also passed. The existing
`R-4.2.3` environment lacks ReacTran, deSolve, Matrix and rootSolve, so a complete
local website render was not performed. The retained GitHub workflow explicitly
installs those dependencies. This preflight supports the relocation repair;
it does not claim that a full Pages build or deployment has run successfully.

Add the new source-list helper and both modified checker/renderer scripts to
the reviewed publication scope. No automatic workflow disabling or website
retirement is part of this change.


## Live noise output inventory refresh, 20:45 UTC

The requested interim inventory update is now in `DATA_MANIFEST.md`. A complete
workspace traversal during 20:44:23–20:45:26 UTC measured **191,129,616,384
allocated bytes (178.00 GiB)**; the updated polarity parent is 129,942,798,336
bytes, manuscript parent 104,192,581,632 bytes, and native-run parent
96,666,099,712 bytes. The manifest labels these as a measurement interval while
production was live and refreshes the affected overview/parent/figure rows.

The new noise main/pilot/numerical-control group is explicitly inventoried with
its purpose, scientific limits, full source paths and status. A joint traversal
at 20:45:34 UTC allocated **2,903,146,496 bytes (2.70 GiB)**: 2,791,211,008 for
the main run, 81,870,848 for the pilot and 30,064,640 for numerical controls.
The main manifest specifies 243 jobs: 216 fluctuating-input replicates, 18
matched-mean controls and nine basal controls. Direct SLURM observation confirmed
12 live tasks at 20:44:23; at 20:45:34 the queue was empty and authoritative
status recorded 243/243 complete, terminal at 20:45:20. Completion of worker
execution does not certify the final figures or biological interpretation.

The root agent was notified of this terminal transition. A final parent/group
measurement refresh remains appropriate after report generation before the
publication commit/push. This follow-up changed only inventory/audit documents;
it did not stage files or delete data.


## Final publication storage measurement, 21:00 UTC

The statistics/report job `21126788` completed at 20:57:36 UTC. The curator
confirmed deterministic gzip-copy session 48632 exited successfully before
measurement: seven copies contain 221,166,834 logical / 248,397,824 allocated
bytes, while all 993,719,748 logical raw-CSV bytes remain preserved. Its mapping
and source hashes are in the topic's `publication_compression_20261017.json`.

One complete `du -x -B1 --max-depth=6 .` traversal ran from
**2026-10-17T20:59:40.273693+00:00 to 2026-10-17T21:00:34.614765+00:00**,
counting hard-linked inodes once. It measured **191,973,122,048 allocated bytes
(178.79 GiB)** overall; polarity parent 130,786,172,928; manuscript parent
105,035,956,224; native-run parent 96,669,343,744; figure parent 2,133,884,928
bytes (1.99 GiB, approximately the manifest's 2-GiB inventory threshold).
A subsequent joint traversal at 21:00:34 UTC measured the native main/pilot/
numerical-control roots plus related noise report exports at **3,556,212,736
bytes (3.31 GiB)**, including 641,589,248 allocated bytes of figure/statistical
exports. `DATA_MANIFEST.md` now records these measurements, scientific roles,
source paths and completed producer/compression states. Native/statistical
production was terminal throughout; small final assembly/organizer records
could still change during the stated measurement interval. No data were
deleted and no index mutation was performed by this measurement.
