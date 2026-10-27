# Ring figures R1 and R2

[Figure R1](figureR1/README.md) · [Figure R2](figureR2/README.md) ·
[Original brief](request/brief_ring_r1r2_20261027.md)

R1 contains the completed matched-noise fast-cytosol campaign and native solver
check. R2 is being rebuilt as Richard's per-perimeter campaign and the new
per-cell arm finish their existing extension checks. Execution-pending
conditions are identified explicitly and remain separate from statistical rate
bounds. The figure scientific readmes give the data-generation time, producer
state and every unexpected result. Native Figures 1–6 and the outline are
preserved.

Rebuild the figures from their versioned compact plotted records, using the
existing `scaling_sandbox` environment, from the repository root:

```bash
conda run -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figureR1_20261027.py
conda run -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render
```

The rebuild entry points also provide extraction commands for the authoritative
HPC records. Bulk extraction and evolution run under SLURM; no new environment
or model retuning is required. Full records stay on the HPC under the paths in
[the data manifest](../../../../../../DATA_MANIFEST.md) and
[retention policy](../../../../../../RETENTION.md).

[Allocated-storage measurement](supporting/metadata/storage_measurement.json) ·
[Release artifact and scope audit](supporting/metadata/release_audit.json)
