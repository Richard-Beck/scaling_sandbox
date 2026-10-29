# Ring figures R1 and R2

[Figure R1](figureR1/README.md) · [Figure R2](figureR2/README.md) · [Original R1/R2 brief](request/brief_ring_r1r2_20261027.md)

R1 retains revision 1. R2 is rebuilt in place as six paired-size cases: three null processes and three low-rate/high-amplitude event-model cases. Each has a kymograph pair and a matching heading-correlation plot. Its former exponent panels, event-triggered analysis, event-log regeneration and drift guide have been deleted. The [R2 readme](figureR2/README.md) begins with Unexpected findings and links the current legend, source records, validation and visual review log.

Portable rebuilds use the existing `scaling_sandbox` environment:

```bash
conda run -n scaling_sandbox python -B polarity1d/manuscript/code/mechanisms/rebuild_figureR1_20261027.py
conda run -n scaling_sandbox python -B polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render
```

Both original event campaigns are complete. The previous release finalizer is historical and is not the current display rebuild entry point. Scientific simulation sources and records remain unchanged. The minimum missing-profile replay is documented inside the R2 package.
