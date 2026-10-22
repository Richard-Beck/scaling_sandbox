# Figure 3: approved Holmes cue-free continuations

The nine cases explicitly approved in the 22 October 07:03 steering were continued from their exact archived 20,000 endpoints to a total cue-free horizon of 80,000. The cue remains zero. Native kinetics, mesh, BDF tolerances and the complete-region resolution definition are unchanged. Original campaign outputs are read only and retain their original hashes.

SLURM array **21279705**, small QoS, nine unthrottled tasks, 1 CPU and 4 GB per task, completed successfully (18–23 seconds each). All trajectories reach the full 80,000 horizon. [Frozen manifest](manifest.json), [submission receipt](submission.json), [collection](results/collection.json), [per-case CSV](results/extension_evidence.csv), [independent field audit](results/independent_validation.json).

Six cases resolve: all five at circumference 30 (sampled brackets from 29,420–29,440 through 32,320–32,340), plus one at circumference 37.8 (34,380–34,400). Three circumference-37.8 cases remain censored at 80,000. Their finite-field stationarity gates are respectively false, false and true; unresolved spatial fate and finite-field stationarity are separate findings.

At matched withdrawal gap 7.6088469922, the middle-size/small-size time ratio is at least 2.373819. The large-size ratio is only bounded below by 0.618429, which does not prove a positive size effect. Both extra-separation ratios remain indeterminate because their numerator and denominator are censored. Other Holmes trajectories retain the original 20,000 horizon; differing bounds must not be read as a decrease with gap. No horizon is substituted as a point estimate.

Full continuation fields and complete-region readouts are stored under `workers/`. Endpoint and original-status copies are immutable inputs. Original full tail paths and hashes remain in the manifest so the original 0–20,000 history can be combined with each new 20,000–80,000 continuation. Scientific importance: these runs resolve the smallest-size Holmes denominator and convert the middle-size matched-gap decomposition from indeterminate to a finite one-sided lower bound, while retaining the unresolved comparisons honestly.

Recollect or verify without new integration:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/figure3_holmes_extensions_20261022.py collect
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/verify_figure3_round2_20261022.py
```

No resubmission is needed. The frozen command is retained for provenance, and the producer rejects duplicate submission receipts.
