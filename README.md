> **deBelly warning (2026-10-08):** the old centered two-end extension is flawed and deprecated. In-loop mean subtraction imposes opposite-end motion and changes internal tension; it is not a reference-frame change. Do not use simulations with `center_boundary_motion=True` for scientific inference or new manuscript figures. Use the independently commanded, reflection-symmetric endpoint adapter described in [the manuscript workflow](polarity1d/manuscript/README.md). Old data are retained only as explicitly deprecated provenance.

# Scaling sandbox

This research workspace studies how cell size affects polarity, mechanochemical
communication, and migration. The active work combines one-dimensional polarity
models, literature-based migration analysis, and microscopy analysis. Scientific
results should distinguish published model assumptions from new adaptations
and size-dependent experiments.

The [current manuscript vision](MANUSCRIPT_VISION.md) records the intended
components, bundles mechanistic explanations with their corresponding measurement
proposals, and reserves a prominent results section for the planned SUM-159 WGD
motility assays. The order of the middle components remains open.

The [manuscript outline](manuscript_outline/README.md) assembles the proposed main
figures, detailed legends, narrative conclusions and pending experimental work.
Its PDF versions preserve the initial integration and subsequent critical review.

The biological motivation is **whole-genome doubling (WGD), cell enlargement,
and increased migration persistence**. The working figure collection is
[current manuscript figs](<current manuscript figs/README.md>); it is a collection
for developing a manuscript, not a completed manuscript. Start with its
[migration figure](<current manuscript figs/main_03_migration_OU.png>) and the
[source analysis](wgd_langevin/refined_ou_20261003/README.md). The compiled
observations include reduced sampled speed but greater late-lag MSD in Xuan,
and reduced path length but greater directness in Mallin. The OU fits infer
longer velocity memory in the enlarged/polyploid groups, conditional on their
model, priors and pooling; these fitted times are distinct from measured
directness or polarity-reversal times.

The 1D length experiments ask whether enlargement, with unchanged local kinetic
parameters, can supply a physical explanation for this persistence. Length here
represents the size change associated with WGD, rather than transient elongation
and shortening during a turn. The
[current length-mechanism report](polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf)
provides candidate polarity mechanisms; their connection to spontaneous turning
and whole-cell migration remains to be established. The literature groups span
stable tetraploids, giant polyploids and treatment-induced PACCs, so they are not
interchangeable causal tests of WGD or size alone.

The M2072 Cellular Potts branch explored a 2D bridge from mechanochemical
polarity to migration and WGD size experiments. The user authorized its
retirement on 2026-10-17. A [compact historical record and two model examples](literature_models/M2072/README.md)
are retained; bulk campaigns and the development pipeline have been removed.
Its exploratory fits did not establish an independently qualified model ensemble
or a general explanation of WGD-associated persistence.

## Report figures: current status (2026-10-17)

Read [the report work status](reports/CURRENT_STATUS.md) before interpreting or
continuing this work. The current comparison is the
[five-model, three-finite-challenge manuscript workflow](polarity1d/manuscript/THREE_TESTS.md),
with corrected independent-endpoint deBelly results and the user-requested
Spearman presentation. The seven-model four-assay tree is deprecated provenance;
its retained summaries and approved historical figures are not the active suite.
The refined migration figure and its sensitivity analyses are a separate
evidence component. Newly assembled mechanism figures remain review checkpoints.

## Current layout

| Area | Role |
| --- | --- |
| [current manuscript figs](<current manuscript figs/README.md>) | Working manuscript figure collection; start with the migration figure for the empirical persistence motivation and follow its source links. |
| [wgd_langevin/refined_ou_20261003/](wgd_langevin/refined_ou_20261003/README.md) | Source of the current migration figure: reviewed literature measurements, conditional OU persistence/noise fits, sensitivity analyses and provenance. |
| [polarity1d/](polarity1d/README.md) | Active WP1 framework: canonical Jilkine models and a De Belly mechanochemical port, common assays, numerical qualification, and curated figures. |
| [literature_models/M2072/](literature_models/M2072/README.md) | Retired 2D Cellular Potts branch; compact historical notes, parameters and two model examples survive the authorized 2026-10-17 storage cleanup. |
| [lfct/](lfct/README.md) | Active microscopy segmentation, track assignment, isolation analysis, and QC workflows. |
| [literature_models/](literature_models/README.md) | Publication-specific provenance and native/reproduction implementations; relocated workflows still require path checks. |
| [migration1d/](migration1d/README.md) | Deprecated moving-interval framework retained for historical reproducibility. |
| migration1d_writer_bundle/ | Historical methods, manuscript material, and figures associated with migration1d. |
| reports/, R/, scenarios/ | Retained transport and intracellular-organization analyses, plus earlier LFCT analyses. |
| scripts/, tests/, workflow/ | Shared and earlier analysis tools, checks, and SLURM launchers; LFCT consolidation is pending. |
| [plan.md](plan.md) | Scientific programme and proposed experiments; a roadmap, not a completion record. |

## Running and reproducing work

Follow the relevant workstream README. For the curated polarity figures, start
with [retained-result reproduction](polarity1d/REPRODUCE_RETAINED_RESULTS.md).
Check existing Conda environments before installing anything: modelling work
uses `scaling_sandbox` and LFCT segmentation uses `cpose`. The root
[environment.yml](environment.yml) describes the modelling environment, not
all dependencies of every workstream. Follow [AGENTS.md](AGENTS.md) for HPC
execution, R runners, and SLURM resource guidance.

## Source and scientific evidence

[DATA_MANIFEST.md](DATA_MANIFEST.md) inventories large local datasets and
simulation outputs, their storage footprints and their scientific purpose.
Update it when generating, moving or removing large outputs, as required by
[AGENTS.md](AGENTS.md).

[RETENTION.md](RETENTION.md) defines what belongs in Git and what requires
separate storage. Source, assumptions, configuration, tests, provenance, and
selected compact evidence are the durable core. Raw images, downloads,
provisional runs, and bulk simulation output remain local and ignored.
Ignoring a file does not back it up or delete it.

The workspace is partway through a reorganization. The old public website is
not maintained by the current working tree. Some historical report assets and
linked local campaign outputs will not be present in a source-only checkout.
There is no `reorg/legacy/` backup in this workspace. Git history can recover
previously committed files; it cannot recover never-committed work. No broader
backup or external result archive has been verified.
