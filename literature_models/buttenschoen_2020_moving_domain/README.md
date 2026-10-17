# Buttenschön 2020 moving-domain model

The Python implementation and renderer are in `src/buttenschoen/`, with the
spatial and WGD drivers in `scripts/` and a focused test in `tests/`.
`src/buttenschoen/README.md` records the implementation details.

Legacy drivers may still assume the former root-level `python/` layout; they
are preserved for a later path-normalization pass.
