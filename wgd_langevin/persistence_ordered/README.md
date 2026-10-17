# Ordered persistence-time rerun

All four variants and all nine experimental pairs from `../persistence_only` are
refitted here. No polarity calibrations are used. Both groups share every parameter
except persistence times, and the ordered alternatives enforce, in hours:

```
0 < tp_reference < tp_WGD+ < 12
0 < tm_reference < tm_WGD+ < 12    (variants with a memory state)
```

There are 36 ordered alternatives and 36 separate equal-time null benchmarks. Each
experiment is fitted independently. The other parameter bounds, data, fitted and
excluded measurements, and covariance closures are retained from the original
campaign. A log-positive reference time plus a logistic fraction of the remaining
interval up to 12 enforces the inequalities by construction. The numerical lower
floor is 1e-6 hours; the reference upper cap is 12−1e-8 hours; gap logits span
−12 to 12. This permits arbitrarily small effects for practical purposes without
claiming that strict ordering excludes an almost-equal solution.

Eight starts, at most 120 least-squares evaluations per search start; two candidates
polished with up to 80 evaluations each. Time seeds start between 0.025 and 8 h;
these are starting values, not search bounds. Other starts perturb generic values.
Each ordered fit includes its corresponding new null with nearly equal times.
Initial, best search, polish and validation checkpoints are written atomically to
individual JSON files. An initial/search/polish record with `complete=false` is
provisional. A separate checkpoint publisher maintains the HTML report, figures and CSVs every
30 seconds while the fitting jobs run. The fit workers themselves do not publish reports.

Linear Xuan/Jemaà objectives use analytic speed and MSD. Other linear simulations
use exact integrated-velocity Gaussian covariance, with a relative numerical jitter
instead of a fixed absolute jitter at very small times. Nonlinear simulations use
the original stochastic Heun equations, with a common integration grid across both
groups. The maximum step is halved until it is below 0.2 times each relaxation time
and a self-activation stability scale. This refinement avoids the former fixed
step being unstable near zero. The independent-seed, maximum-step halving and burn-in
doubling checks are retained. If time-constant stability refinement already dictates
a finer step, lowering the maximum step may leave the actual step unchanged.

To keep the campaign bounded, candidates requiring more than 2,000,000 integration
steps per independent cell are flagged as numerically unresolved and excluded from
optimizer candidates via a large residual. Nonfinite simulations are treated the
same way. These are computational failures, never evidence of biological or
mathematical infeasibility. Their counts and any unresolved validation settings are
retained. An extremely short-time optimum may therefore require a dedicated solver
study before interpretation. JSON `complete=true` means the scheduled search and
checks finished, not that a fit is adequate or numerically resolved everywhere.

Search/polish/validation ensembles: linear 512/2048/8192 independent cells where
sampling is needed; nonlinear 96/512/2048. All objectives and acceptance checks
remain descriptive, based on the existing reading intervals or Hosny IQRs.

Run in the existing `scaling_sandbox` Conda environment. SLURM red/small, one CPU,
3 GB, four-hour per-task limit; no array concurrency throttle. Two arrays both use
indices 0–35. The second uses ORDERED_OFFSET=36 and `aftercorr` to begin each
ordered task after its matching null, without waiting for unrelated null fits.
The manifest and execution record are under `results/`.

```
python -m unittest discover -s wgd_langevin/persistence_ordered -p 'test_*.py' -v
python wgd_langevin/persistence_ordered/campaign.py --manifest
sbatch --array=0-35 wgd_langevin/persistence_ordered/run.slurm
sbatch --array=0-35 --export=ALL,ORDERED_OFFSET=36 --dependency=aftercorr:NULL_JOB_ID wgd_langevin/persistence_ordered/run.slurm
```

Source/data hashes are retained in every checkpoint. Completed results are reused
only when hashes agree. Incomplete tasks currently restart their seeded search;
the last atomic checkpoint remains available if a task hits its wall limit.

## Automatically updated report

Open [the live report](results/report.html). It follows the previous campaign
layout: experiment selector, model comparisons, persistence times, shared parameters,
observed-versus-predicted figures, excluded checks and downloadable CSVs. Provisional
initial/search/polish checkpoints are visibly distinguished from completed fits.
Validation is never inferred from absent or numerically unresolved checks.

The publisher reads atomic checkpoints without importing or modifying the fitting
code. It polls files and SLURM every 30 seconds, atomically replaces artifacts, and
marks the final snapshot when all fits complete or no fit jobs remain active.
Generation time, last checkpoint, scheduler check and publisher state are visible.
The browser has no timed refresh; opening or manually refreshing loads the latest
server-generated snapshot. Figures regenerate only when their inputs change.

Start a publisher with `sbatch wgd_langevin/persistence_ordered/report.slurm`. It
uses red/small, one CPU, 2 GB, a nine-hour limit covering the two four-hour fit
stages, and exits when production stops. For one snapshot, run
`python wgd_langevin/persistence_ordered/publish.py`. Live status and all checkpoint
stages are exported in `results/report_status.json` and `results/progress.csv`.
