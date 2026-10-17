# Refined migration OU refit — 3 October 2026

This standalone update implements
[the agreed fitting disposition](../../literature_review/migration_data/OU_FITTING_DISPOSITION.md)
for the five-study Gaussian OU hierarchy used by the selected manuscript migration
figure. The original experiments, model, population contrasts, five prior/pooling
settings, and Hosny sensitivity comparisons are retained. The earlier campaigns
and manuscript figures remain available at their original locations.

Updated manuscript exports:
[main migration figure](figures/main_03_migration_OU.png) and
[Hosny sensitivity supplement](figures/supp_04_Hosny_sensitivity.png), each also
available as PDF and SVG. Read [fit results](results/RESULTS.md),
[population MSD residuals](figures/xuan_population_MSD_diagnostics.png), and
[representative movie findings](diagnostics/FINDINGS.md) for adequacy and limits.
The [comparison with the earlier selected fit](results/comparison_to_previous.csv)
reports parameter changes only; altered targets make objective scores incomparable.

All 40 starts converge. The primary pooled persistence/noise folds are approximately
15.80/0.131; omitting Hosny and refitting gives 19.71/0.107. These conditional
summaries do not establish a mechanism. The full Xuan curves show systematic shape
residuals, and the smaller-nucleus movie tracks show prolonged slow episodes and
longer correlations than the conditional OU prediction. The data-only Jacobian
has rank 18/22, with the four remaining directions constrained by regularization.

## Where to look

- `data/`: frozen reviewed fitting inputs, reserved measurements, full digitized
  Xuan population MSD curves, and source/uncertainty documentation.
- `scripts/`: fitting, numerical validation, figure rendering, and movie diagnostics.
- `results/`: individual bounded multistart fits, selected scenario results,
  frozen-parameter simulation checks, and compact numerical exports.
- `figures/`: updated `main_03_migration_OU` and
  `supp_04_Hosny_sensitivity` manuscript exports, captions, and selected diagnostics.
- `diagnostics/`: representative Xuan movie comparisons and their numerical tables.
- `ARTIFACTS.csv`: final retained file sizes and SHA-256 hashes.
- `scripts/*.slurm` and `scripts/hosny_sensitivity/run.slurm`: SLURM entry
  points; execution logs and identifiers are retained with the results.

Read the input notes before interpreting numerical precision. Raster/calibration
allowances, SEM, and quartiles describe different quantities. The loss remains
descriptive regularized least squares, with explicitly declared scales; it is
neither a likelihood nor a source of statistical confidence intervals.

## Scope and interpretation

The hierarchy comprises Xuan 2018, Mallin 2023, Jemaà 2017 Cl1/Cl8,
Jemaà 2023 D1/T1, and the Hosny hypoxic Non-PACC/PACC baseline. Treatments,
additional clones, and representative movie nuclei are not population fit targets.
Five comparisons retain the primary 10-minute persistence prior, 5- and
20-minute priors, and stronger/weaker pooling. Relaxation times and velocity-noise
strength may differ between groups; signed drift is fitted only for Hosny.

Both published Xuan MSD curves now contribute across early, intermediate, and
late lags alongside sampled mean speed. Each curve is one balanced loss block;
digitized columns are not counted as independent biological observations. The
legacy early-MSD check overlaps the upgraded curve and is therefore a correlated
check, not independent validation.

The movie supplies eight visually reviewed full-window nuclear-centroid tracks:
one apparent giant and seven smaller nuclei. Their appearance does not establish
ploidy. No recovered acquisition timestamps resolve the movie's 12-hour/16-hour
timing discrepancy. Physical-time diagnostics are conditional on named timing
scenarios, and centroid localization allowances are sensitivity assumptions.

The retained Hosny approximation fits medians analytically, checks quartiles by
simulation, and profiles drift and leave-Hosny-out sensitivity. Its observations
come from the recovered preprint; equivalence to the journal version is unresolved.

## Environment and execution

Use the existing `/home/4473331/.conda/envs/scaling_sandbox/bin/python`.
No environment or dependencies were installed. Compute-heavy fits, validation,
and sensitivity calculations use SLURM, with one CPU per task and `small` QoS.
Set `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for local inspections/rendering.
Execution records and artifact provenance accompany the outputs. Intermediate
caches and temporary files are removed after final validation; source snapshots,
scientific tables, and execution logs are retained.

To rerun from the repository root, first read [the frozen input notes](data/README.md).
Regenerate and validate inputs before launching any fits, then keep the inputs
fixed until all dependent producers stop:

```sh
PY=/home/4473331/.conda/envs/scaling_sandbox/bin/python
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
$PY -B wgd_langevin/refined_ou_20261003/data/prepare_inputs.py
$PY -B wgd_langevin/refined_ou_20261003/data/validate_inputs.py
sbatch --array=0-39 wgd_langevin/refined_ou_20261003/scripts/run.slurm
sbatch --dependency=afterok:FIT_JOB_ID wgd_langevin/refined_ou_20261003/scripts/validate.slurm
sbatch --dependency=afterok:FIT_JOB_ID wgd_langevin/refined_ou_20261003/scripts/hosny_sensitivity/run.slurm
```

After validation and sensitivity have finished, render the main figure and movie
checks with the scripts documented in their captions and
[movie methods](diagnostics/README.md). The final
[integrity audit](scripts/audit_final.py) independently recomputes all 40 objectives,
checks frozen input/source hashes, verifies balanced curve weights and all six
manuscript exports, and rejects leftover temporary files.

```sh
$PY -B wgd_langevin/refined_ou_20261003/scripts/diagnostics.py
$PY -B wgd_langevin/refined_ou_20261003/scripts/compare_previous.py
$PY -B wgd_langevin/refined_ou_20261003/scripts/main_surface.py
$PY -B wgd_langevin/refined_ou_20261003/scripts/build_manuscript_figures.py
$PY -B wgd_langevin/refined_ou_20261003/scripts/curve_diagnostics.py
$PY -B wgd_langevin/refined_ou_20261003/scripts/movie_diagnostics.py --fit wgd_langevin/refined_ou_20261003/results/selected_fits.json
$PY -B wgd_langevin/refined_ou_20261003/scripts/audit_final.py
$PY -B wgd_langevin/refined_ou_20261003/scripts/inventory.py
```

The existing output files are a completed run. Preserve them before a deliberate
future refit; the launch commands above write into this directory.
