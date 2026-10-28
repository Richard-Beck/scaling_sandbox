# Ring figures R1 and R2

[Figure R1](figureR1/README.md) · [Figure R2](figureR2/README.md) · [Original brief](request/brief_ring_r1r2_20261027.md)

Revision 1 retains every R1 A–H and R2 A–G panel at its original path. Sizes use the fixed orange/green/blue palette, kymographs use linear viridis, and model identity comes from facets. R1 C–D use the complete saved f = 0.01 post-burn records. R2 G reports event-triggered excess over random windows in the same runs, across w = 2–64 original frames.

The current R2 extraction remains a campaign snapshot. Pending estimates occupy a separate band below the α axes and are distinct from finite-record bounds. Both scientific readmes report data-generation time, producer state and Unexpected results. Existing simulation campaigns and frozen inputs are preserved; no new simulation campaign or retuning is introduced.

Portable rebuild from repository root, using the existing `scaling_sandbox` environment:

```bash
conda run -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figureR1_20261027.py
conda run -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render
```

The existing dependency finalizer follows both event campaigns through their extension checks. Once both report COMPLETE, it rebuilds R2 again, runs the final-data numerical/display checks, refreshes links and allocated-storage inventory, then commits and pushes only the authorized revision paths. It uses a private Git index to preserve unrelated staged and working changes; the current branch and validated source hashes are checked before publication. No additional simulation is submitted. A failed check or push is recorded for attention rather than retried as a new campaign.

[Scientific/display audit](supporting/metadata/revision1_audit.json) · [Delivery setup](supporting/metadata/revision1_delivery.json) · [Completion state](supporting/metadata/completion_audit.json) · [Final delivery record](supporting/metadata/final.txt) · [Release artifact audit](supporting/metadata/release_audit.json) · [Allocated storage](supporting/metadata/storage_measurement.json)
