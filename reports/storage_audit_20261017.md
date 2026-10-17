# Polarity storage audit — 17 October 2026

This was the initial read-only audit of `polarity1d`. The user subsequently
authorized physical removal of large deprecated four-assay outputs; execution
status and measured results are in the
[cleanup ledger](storage_cleanup_four_assay_20261017.json) and
[retirement note](../polarity1d/four_assay/OUTPUT_RETIREMENT_20261017.md).
The estimates and narrower six-root proposal below describe the pre-cleanup
snapshot. Manuscript empty-log and superseded-experiment suggestions were not
included in the later four-assay-only authorization.

The complete repository scan measured approximately 238 GiB in `polarity1d`:
128 GiB in `four_assay`, 88 GiB in `manuscript`, 15 GiB in `archive`, and 9 GiB
in `results`. Rounded directory totals should not be added: some imported
archives use hard links.

## Recommended historical payload retirement

The precise proposed predicate is: **regular files with names ending in `.npz`,
recursively under the `work/` directories of these six roots**. Do not follow
symlinked directories. This proposal does not remove entire campaign roots.

All paths below are relative to the repository root:

1. `polarity1d/four_assay/pilots/width_full_lengths_20261002/`
2. `polarity1d/four_assay/pilots/width_mean_slope_20261002/`
3. `polarity1d/four_assay/pilots/polarity_establishment_20261003/`
4. `polarity1d/four_assay/pilots/polarity_establishment_ref15_relative20_20261003/`
5. `polarity1d/four_assay/pilots/polarity_establishment_submission_attempt_20261003/`
6. `polarity1d/four_assay/pilots/polarity_establishment_basal_copy_attempt_20261003/`

Preserve all non-NPZ records: JSON trial/case records, compressed scalar records,
CSV summaries, results and figures, manifests, provenance, frozen source trees,
launch records and logs. Before implementation, add explicit retirement notices
to distinguish preserved scalar evidence from removed native-array coverage.

| Disjoint contribution to joint retirement | Allocated GiB |
| --- | ---: |
| Two old width campaigns, exclusive payload | 38.768 |
| Two old establishment campaigns, exclusive payload | 64.870 |
| Two failed intermediate campaigns, exclusive payload | 2.098 |
| Payload shared between establishment and failed-copy trees, counted once | 6.027 |
| **Total actually reclaimable by retiring all six together** | **111.763** |

The inode audit counted **418,926 unique NPZ inodes**, allocating
**120,004,370,432 bytes**. For every inode, the number of file names encountered
inside the proposed six work trees equalled its link count: **zero hard links
remain outside this candidate set**. This is the joint allocated-space estimate,
not a sum of independent campaign `du` sizes. In particular, the failed copy
attempt's approximately 6.6 GiB footprint overlaps the old establishment tree.

These are retired scientific trial arrays, not uniformly redundant copies.
Deleting them would prevent full native-array re-audit or arbitrary retrospective
analysis of those old scans without rerunning. The old width workflows and
preparation scans are superseded; the standardized-reference predecessor also
contains the incorrect zero-total-input LEGI withdrawal protocol. Compact records
should retain those scientific limitations.

## Selected-input independence and historical symlink cost

The existing canonical retention inventory is
[`polarity1d/manuscript/RETENTION.md`](../polarity1d/manuscript/RETENTION.md).
The audit checked the following independent manuscript copies:

- All **4,952 selected establishment NPZ files** listed in
  `manuscript/data/provenance/data_sources.json` exist as local regular files,
  with **no symlinks**; their recorded logical size totals 462,557,720 bytes.
  This check verified existence and file type, not every establishment hash.
- All **54 selected diagnostic reconstruction NPZ files** listed in
  `manuscript/data/provenance/legacy_slimming_plan.json` passed their recorded
  **SHA-256 checks**, covering 102,195,863 bytes.
- The diagnostic builder explicitly reads the canonical
  `manuscript/data/reversal_fields/reconstruction_inputs/` copies.

The historical
`four_assay/pilots/polarity_establishment_ref15_relative20_basal_20261003/work/cases/`
tree has symlinks for six unchanged models into the old
`polarity_establishment_ref15_relative20_20261003` work tree. Removing the proposed
old arrays would therefore retire native-array coverage behind those historical
links. It would not remove the independent canonical selected manuscript copies.
This dependency must be described in retirement notices; do not imply that all
105,947 old duration trials remain locally re-auditable afterward. Canonical
hash verification and selected operational-path checks should precede any
subsequent physical deletion, as required by the retention inventory.

## Empty scheduler logs

Under `polarity1d/manuscript/data/runs/`, pre-17-October campaign log directories
contain **53,310 zero-byte regular `.out`/`.err` files**, all with one hard link.
They allocate **1,249.453 MiB = 1.2202 GiB** because of filesystem allocation
granularity. After confirming those completed campaigns have no running
producers, removing only these empty scheduler logs would retain all native
arrays, records and nonempty logs. The largest concentrations are in
`corrected_reversal_20261003`, `fixed_slope_reversal_20261006`,
`three_test_reversal_20261006`, `wave_pinning_diffusion_contrast_20261008`, and
`wave_pinning_slow_diffusion_20261008`.

The historical NPZ proposal plus these empty logs would free approximately
**112.98 GiB (about 113 GiB)**.

## Optional smaller evidence retirement

`polarity1d/manuscript/data/runs/wave_pinning_slow_diffusion_20261008/work/chunks/`
and `work/preparations/` hold approximately **1.66 GiB** of unique, single-link
arrays for the explicitly superseded experiment with active diffusivity 0.01 and
inactive diffusivity 1. The requested replacement is
`wave_pinning_diffusion_contrast_20261008/`, with inactive diffusivity 100.
Retirement could preserve the old summary/PDF, frozen inputs and provenance,
while sacrificing its independent native-array re-audit. Its empty logs overlap
the preceding log estimate and must not be counted again.

## Protected material

Keep the accepted v8 baseline; the approved `assay_explanation_v2` keeper with its
sources, reconstructed fields and provenance; canonical model sources; current
corrected comparison records and selected inputs; report PDFs and source native
arrays; and all 17-October circular preparation/distance/robustness runs and
figure inputs. Keep
`polarity1d/manuscript/data/runs/debelly_no_protrusion_20261010/` (approximately
29 GiB): it contains unique audited histories supporting the 430-integration
development experiment and linked figure. An archive label or reproduction
recipe alone does not make scientific payload expendable.

This report records a bounded storage audit. It does not authorize deletion or
resume the paused scientific figure work.
