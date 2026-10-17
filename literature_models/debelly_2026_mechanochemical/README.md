# De Belly 2026 mechanochemical polarity model

This bundle preserves the Rac-Rho-MCA mechanochemical polarity circuit and
runnable 1D Julia implementation associated with the 2026 paper under
`src/debelly/` and the explanatory Python 1D implementation under
`src/debelly_1d/`. Diagnostic, plotting, and SLURM files are grouped in
`scripts/`, with Python tests in `tests/`.

The membrane-cortex mechanics derive from the authors' 2023 work; that
mechanical lineage does not make the complete polarity model a standalone
2023 model.

For WP1 simulations, common protocols/phenotypes, and the performance-selected
conservative 1D port, use `../../polarity1d/`. This directory remains the
native/provenance bundle and selected cross-check implementation.

The Julia trees retain their own Project/Manifest files. Use those existing
environments rather than introducing another Conda environment.
