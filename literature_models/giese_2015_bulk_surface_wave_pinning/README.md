# Giese 2015 bulk-surface wave pinning

`src/giese_wp/` contains the model-specific solver and metrics.
`src/bulk_surface/` contains the shared finite-element bulk/membrane engine.
Reproduction, validation, plotting, and scheduler entry points are preserved
under `scripts/`; tests are under `tests/`.

The 2026-10-17 publication pass updates imports to `src/` and submission
paths to this model bundle. From the repository root, run
`conda run -n scaling_sandbox python -m pytest -q literature_models/giese_2015_bulk_surface_wave_pinning/tests`.
For native campaigns, change into this folder and run `scripts/submit_giese_native.sh`;
its launcher creates the expected log folder. Rebuilt native outputs remain local
under `reports/output/`; they are excluded from Git. Portability checks and focused
tests do not substitute for rerunning the full native reproduction campaign.
