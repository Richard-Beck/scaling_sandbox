# Taniguchi 2013 excitable-wave model

`src/taniguchi/` contains the fixed-circle and deformable phase-field
implementations plus metrics and rendering support. Campaign-generation,
validation, plotting, and SLURM files are retained in `scripts/`; focused tests
are under `tests/`.

The 2026-10-17 publication pass updates imports to `src/` and SLURM paths
to this bundle. From the repository root, run
`conda run -n scaling_sandbox python -m pytest -q literature_models/taniguchi_2013_excitable_waves/tests`.
For a native campaign, change into this folder, create `reports/output/taniguchi/logs`,
and submit `sbatch scripts/run_taniguchi_fixed.slurm` or
`sbatch scripts/run_taniguchi_phase.slurm`. Native outputs stay local under
`reports/output/` and are excluded from Git. No `reorg/legacy/` backup exists
in the current workspace; this source bundle does not promise recovery of all
historical stochastic outputs.
