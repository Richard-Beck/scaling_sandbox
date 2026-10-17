# Chiou 2015 MCAS peak competition

This bundle preserves the R port, published initial state, validation and
reproduction scripts, and scientific report. The finite-strip migration
adapter is maintained in `../../migration1d/src/migration1d/chiou.py`.

The source is a transcription of the authors' supplementary MATLAB model.
The repository's length/abundance scaling experiments are extensions, not
claims about the source paper.


The 2026-10-17 publication pass repairs relocated R/Python entry points. R
scripts locate this bundle when invoked with `Rscript` or `source()`, then use
`src/chiou_mcas.R`, `params-or-config/initial_state_1D.mat`, and the recovered
published-size scenario documented in
[`params-or-config/RECOVERY_20261017.md`](params-or-config/RECOVERY_20261017.md).
Run the compact checks with an existing R environment:
`conda run -n R-4.2.3 Rscript literature_models/chiou_2015_mcas/scripts/check_chiou_port.R`
from the repository root. Prefer `agentRrunner.R` when it is available.
Native reproduction and exported report outputs stay local in this bundle's
`reports/output/` and `reports/assets/`; these directories are ignored by Git.
The full reproduction launcher uses parallel workers and should run through
SLURM with explicitly selected resource limits.
