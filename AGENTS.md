# Agent Notes

- Account for every large file or related output group (approximately 2 GiB or
  more) in the root [DATA_MANIFEST.md](DATA_MANIFEST.md). Record its path,
  measured allocated storage, purpose and scientific importance/status; related
  files may be grouped under a meaningful parent folder. Tasks that create,
  substantially enlarge, move or remove large outputs must update this manifest
  and its measurement date before handing off. Inventory does not imply
  permission to delete data.
- For the biological motivation of size/length effects, start with
  `current manuscript figs/README.md` and the migration figure's source in
  `wgd_langevin/refined_ou_20261003/`. Length represents WGD-associated enlargement;
  the empirical motivation includes greater migration persistence. Distinguish
  observed migration summaries, conditional OU persistence estimates and 1D
  polarity-reversal readouts. The figure collection is a working manuscript
  resource, not evidence that a completed manuscript exists.
- For report figures and the seven-model four-assay comparison, read `reports/CURRENT_STATUS.md` and `polarity1d/four_assay/README.md` before continuing or interpreting results. The active four-assay pipeline and accepted revised-v8 baseline are in `polarity1d/four_assay/`; old report campaigns are historical archives. Completed output generation does not mean the suite or its length-response interpretation is final.
- This repository is used on an HPC system. Check existing Conda environments before creating or installing a new one.
- Use `agentRrunner.R` for R work when available; it provides the expected execution environment for LLM-run R scripts.
- Use SLURM for compute-heavy work (roughly >4 cores, >10 minutes, GPU work, or >8 GB RAM). Select a QoS suited to the job shape rather than defaulting to `normal`.
- Do not add a SLURM array `%` concurrency throttle unless resource, filesystem, scheduler, or user constraints justify it.
- Avoid `hpctpa3pc0066`, `hpctpa3pc0067`, `hpctpa3pc0068`, and `hpctpa3pc0069` for CUDA work. On 2026-08-29, repeated jobs allocated `gpu:a30:1` on each of these nodes but failed with `CUDA_ERROR_NO_DEVICE`; exclude them until the cluster configuration is repaired.
- For M2072 work, read `literature_models/M2072/AGENTS.md` before choosing campaign inputs. It records the user's intended original V1–V6 lineage and the explicit development-only exception for existing refined-module runs.
- In reports and documents, use properly rendered mathematical formatting consistent with scientific publications. Put substantial equations on separate centred display lines, with equation numbers right aligned; use true mathematical subscripts, superscripts, fractions and derivatives. LaTeX/Markdown math is expected where the document renderer supports it. Avoid raw LaTeX/Markdown equations only in CLI chat responses, where they do not render.
- For HTML reports, do not use timed browser self-refresh. Keep generated data current while its producer runs so opening or manually refreshing the page loads the latest available report; show the data-generation time clearly, especially if the producer has stopped.
