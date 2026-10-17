# Revision-8 post-processing renderer

User-requested mode classification, OLS lines, combined preparation/response
plot and slope heatmap. Read revision8/results/README.md and
visualization_provenance.json for definitions, units and limitations.

Run with the existing scaling_sandbox environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python reports/four_test_suite/revision8_visualization/render.py

The frozen simulation campaign is preserved; raw data and original exports
remain available. Renderer inputs are original CSVs in before_mode_revision/,
plus raw competition records for peak ratios. Outputs replace current figure
exports and derived endpoint statuses, and add regression/flag CSVs.
