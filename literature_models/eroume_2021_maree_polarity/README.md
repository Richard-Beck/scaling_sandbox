# Eroumé 2021 / Marée polarity model

`src/maree_migration/` preserves the standalone nine-field implementation,
metrics, stimulus, rendering, and its validation notes. The related scripts
and tests are retained alongside it. The 1D moving-interval adapter lives in
`../../migration1d/src/migration1d/eroume.py` and shares the same reaction
equations.

The fixed-shape implementation is literature-derived; WGD geometry and
endpoint-migration experiments are repository extensions.
