# Repository retention and commit preparation

## Working policy

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
