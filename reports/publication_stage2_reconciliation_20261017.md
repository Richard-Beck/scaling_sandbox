# Stage-2 repository reconciliation decision — 17 October 2026

V1 is published as `6199cdc`. This follow-up audited the remaining workspace
and implemented reviewed append-only ignore/source-landing rules. It did not
recover original files, modify the Git index, commit or push. New V2
renderer/assets and mechanism-review edits
belong to the root agent's ongoing scientific review and were left alone.

## Held tracked originals

There are **81 absent tracked paths**, containing **8,478,861 logical bytes** in
published `HEAD`. The pre-recovery identities are recorded in
[`publication_held_deletions_pre_recovery_20261017.json`](publication_held_deletions_pre_recovery_20261017.json),
including each current `HEAD` SHA-256 and the commit identity. They include the
public-site workflow/landing resources, old report snapshots, original VCell R
sources, earlier Julia/Taniguchi variants and legacy launchers. These 81 paths
were deliberately excluded from V1 deletion selection; all remain in the
published repository. They are distinct from the 31 verified relocations already
committed in V1.

**Accepted decision: keep all 81 originals absent in this local workspace,
while retaining their intact copies in published `6199cdc` and the current Git
tree.** No recovery is performed, and none of these 81 deletions is staged or
published. Recreating those files could reverse the earlier local workspace
arrangement without evidence of its intent. The prestate ledger preserves exact
path/size/hash evidence; a clean local status is not a publication requirement.
The 31 verified relocations already committed in V1 remain in effect.

This separates preservation in Git from local working-tree state. Independent
checkout/index-derived validation uses the retained original site resources,
report snapshots and original source variants. The active local workflow uses
the canonical bundles and current manuscript pipeline. Some old Julia scripts
refer to helpers moved by the 31 verified relocations; their preserved bytes
are historical evidence, and their portable entry points are under
`literature_models/debelly_2026_mechanochemical/`. Source retention does not
certify every older launcher as current or fully runnable.

The evidence supports preservation. Root `RETENTION.md` says to keep deprecated
code while its unique scientific value is reviewed, permits preserving
recoverable old paths, pairs verified relocations with deletion, and reserves
website retirement for a separate decision. The explicit authorized physical
retirements concern M2072 and large four-assay outputs. No direct user message or
user-attributed local record was found authorizing general removal of these 81
tracked originals. The old `docs/steering.txt` and `developer_note.txt` contain
scientific direction, but their authorship is not established; they cannot be
treated as direct user authorization for source deletion or website retirement.
The current user-provided AGENTS instructions specifically preserve the accepted
v8 baseline and distinguish historical provenance from current scientific results.

## Remaining scientific source and compact evidence

[`publication_historical_sources_20261017.txt`](publication_historical_sources_20261017.txt)
selects **146 historical source/documentation files, 907,783 logical bytes**:
86 Python files, 32 Markdown files and 28 SLURM launchers. Of these, 107 are
currently hidden by whole-branch ignore rules. Most are earlier WGD active
Brownian/OU/Poisson/feedback model variants, fitting code, scientific reports,
tests and launchers. They are inexpensive, unique scientific history; preserve
them as explicitly historical source while keeping their generated results
local. All selected Python files parse successfully. The source selection also
includes the old x50 diagnostic report builder, whose input paths depend on
local scratch and retired historical viewer arrays. Preserve that builder as
provenance, without promoting it as a current executable manuscript workflow.

[`publication_compact_history_20261017.txt`](publication_compact_history_20261017.txt)
selects **11 compact evidence files, 4,664,120 logical bytes**:

- The nine primary accepted-v8 baseline scalar/figure exports, selected once
  from `polarity1d/four_assay/baseline/`; repeated frozen copies stay local.
- `migration1d_writer_bundle/EXPERIMENT_SUMMARIES.xlsx`, which is an explicit
  dependency of the already tracked `METHODS.md` and `FULL_DATA_CATALOG.md`.
  Its metadata identifies openpyxl generation, so it is a descriptive historical
  summary rather than direct user-authored preference evidence.
- The retained Jilkine model XML in the corrected basal pilot's reference source
  folder, as a compact model-source/provenance artifact.

These are historical selections. Current manuscript source/figure authority
remains the corrected workflow, refined OU evidence and the root agent's V2
revisions. The new `scripts/render_manuscript_mechanisms.py`,
`render_manuscript_migration_v2.py`, `render_noise_interpretation_v2.py` and their
corresponding selected assets should enter through the root's final V2 review,
with their own data/source/hash closure. Their in-progress files were not copied
into this historical selector.

## Implemented ignore boundary

[`publication_stage2_ignore_proposal_20261017.txt`](publication_stage2_ignore_proposal_20261017.txt)
records the reviewed append-only rules now applied to root `.gitignore`.
It keeps a default-deny boundary for the older WGD branches while exposing only their exact reviewed source files. It
also protects nested historical Hosny `results/` trees, duplicate frozen
baseline records, old trial-report delivery PDFs, diagnostic preview/report
copies and known generated pilot outputs. It leaves the primary compact baseline
and workbook eligible. No native NPZ, movie, mask, downloaded archive or bulk
simulation tree is added by this source selector.

A temporary independent Git repository first tested these rules: all
**157 source/compact paths remained eligible**, and six representative raw or
duplicate payload paths were ignored. The same eligibility/exclusion checks
were then repeated against the actual root rules after application. No root
index mutation was performed. The WGD landing README now identifies historical
source preservation, current scientific authority and local-only output limits;
`RETENTION.md` records the bounded 81-path local-state decision. The existing
`persistence_only/campaign.py` helper used by current refined OU rendering
retains its narrow exception. No blanket force-adds are needed.

The actual-root implementation check at **21:15:36 UTC** verified all 157
selected paths eligible, all six raw exclusions enforced, all 81 original local
paths still absent and all 81 published HEAD hashes unchanged. The machine
record is [`publication_stage2_preservation_validation_20261017.json`](publication_stage2_preservation_validation_20261017.json).

## Lightweight validation

The first LFCT collection attempt in `scaling_sandbox` failed because that
modelling environment lacks `cv2` and `tifffile`. The LFCT README specifies the
existing `cpose` environment; rerunning there passed **16 LFCT core tests** and
**8 reassignment tests**. These cover unique/ambiguous mask assignment, occupied
masks, nuclear overlap fallback, out-of-bounds positions, exact mask clearance,
uncertain neighbors, source indexing, branch/gap handling, eligibility duration
and compressed-mask round trips. No GPU work, segmentation, downloads or package
installation ran. The historical root WGD model's **five invariant/moment tests**
also passed in the existing modelling environment. These focused checks do not
qualify all historical fitting campaigns.

## Bounded completion sequence

1. Keep the 81-path prestate ledger and the originals' published Git copies;
   preserve their local absences and all 31 committed relocations. Do not stage
   the 81 preexisting deletions.
2. The reviewed ignore rules and historical-source landing note are implemented.
   The root agent can stage the exact 146 historical source paths, 11 compact
   evidence paths and reviewed source/retention/audit edits. No old campaign
   reruns are needed to preserve their code.
3. Add the root agent's final V2 scientific sources/assets separately, inspect
   every remaining tracked change and untracked path, and account for any
   remaining local-only files by their defined output/scratch boundary.
4. Recheck source selectors, largest staged objects and archive exclusions;
   run the relevant focused V2 checks and the public-site source selection/safety
   check. Verify the held-path Git identities and preserve the actual CI result.
   Update any measured large-output groups if V2 materially enlarges them, then commit/push
   the complete reviewed scope.

This plan closes necessary source/evidence preservation without interpreting an
ignored output as deleted data or treating a tidy status as scientific validation.
