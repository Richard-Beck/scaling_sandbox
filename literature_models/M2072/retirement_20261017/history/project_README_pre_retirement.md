# M2072 deBelly CPM campaigns

## Anchor and motility workflow (authorized 2026-09-22)

[Implementation and protocols](src/anchor_motility/README.md): numerical pilot →
32 mechanical anchors → native perimeter-response identification → circular-PDE
biochemical initialization → native speed/MSD and relaxed deBelly refinement →
fresh-seed assessment. New campaigns do not optimize morphology matching.
The [numerical pilot](tmp/anchor_numerical_pilot_20260922/report/README.md)
records whether coarser updates can replace the reference configuration.
This supersedes the 2026-09-19 workflow below for new campaigns; historical
results and their actual objectives remain unchanged.

Deployed [anchor-motility campaign](tmp/anchor_motility_20260922/report/index.html)
uses the reference numerics after all 45 numerical-pilot runs completed; none
of the coarser configurations passed every comparison. Its frozen sources,
stage dependencies, parameters and current status are linked from the report.

Scratch cleanup authorized 2026-09-19: the [workflow provenance bundle](tmp/workflow_provenance_20260919/README.md)
retains feature-selection, pair/scoring evidence and compact campaign summaries.
Its manifest maps historical paths below to retained copies; obsolete trajectory
and gallery links may no longer resolve after cleanup. The active joint search
and small original-path compatibility inputs are protected and remain in place.

## Adopted parameterization workflow (2026-09-19)

1. Optimize the biochemical component separately against deBelly assays.
2. Develop/use the tension emulator and coupled diversity screen to select
   several mechanically distinct coupled parameter sets with predicted tension
   responses reasonably compatible with the biochemical target responses.
3. Use those sets directly as seeds for native-CPM joint refinement of
   biochemical and mechanical/coupling parameters against morphology (pair
   loss), baseline movement (speed/MSD), and migration transfer.
4. Validate selected leaders independently with fresh seeds and more replication.

The emulator is for seed selection, not a surrogate assumed valid after
biochemistry/coupling changes. Mechanical diversity does not guarantee diverse
or adequate cell morphologies. The intermediate large mechanical sweep,
eight-hour morphology/movement emulator, 48-candidate pool, and local emulator
screen are optional development branches, not prerequisites.

See [AGENTS.md](AGENTS.md) for authoritative workflow and lineage instructions,
and the [current joint-refinement runner](src/overnight_refinement/README.md).
The [original-six overnight report](tmp/joint_original6_overnight_20260918/report/index.html)
supports adopting this approach; its winners remain search results, not held-out
validation. "Original/unrefined" seeds already have upstream-fitted biochemistry;
they have not undergone the older coupled joint reoptimization.
The existing eight-hour/48-candidate branch instead used jointly refined V1–V6;
its historical provenance is preserved, not made the default for fresh runs.

## Project role and current direction

M2072 is the active 2D bridge between the project's polarity mechanisms and
experimentally observable cell migration. Its purpose is to produce a compact,
reproducible parameterization workflow and a small ensemble of CPM cells with
credible morphology, shape dynamics, migration, and mechanochemical responses.
The resulting ensemble will be used for a held-out WGD scaling experiment and
for tests of whether weakening long-range Rac-Rho coordination has a larger
effect in enlarged cells.

The existing campaigns are exploratory evidence toward that goal, not a
requirement to optimize every available morphology, biochemical, mechanical,
or emulator objective. In particular, the cells in
[`tmp/migration_volume_half_vs_one_recalibrated_20260910/`](tmp/migration_volume_half_vs_one_recalibrated_20260910/README.md)
do not provide an adequate target-like morphological baseline for final WGD
inference. Some candidates developed through
[`tmp/initial_population48_1h_20260917/`](tmp/initial_population48_1h_20260917/README.md)
are more promising, but remain provisional. The workflow above is now selected;
qualification evidence for the final ensemble remains to be established.

Remaining high-level work is to:

* consolidate the adopted workflow and establish ensemble qualification evidence;
* retain a small ensemble of plausible parameterizations selected without
  reference to their WGD outcomes;
* freeze and audit the WGD scaling rule, then apply it across that ensemble;
* compare the resulting migration phenotypes with literature-derived WGD
  observations at the level supported by those data;
* incorporate additional experimental constraints where they resolve a known
  adequacy problem, including dynamic perturbation evidence summarized in the
  [de Seze review](../../literature_review/deSeze_2025/experiments_and_outcomes.md);
* test practical mechanochemical ablations and the intended experimental design.

The experimental sources have different roles and biological contexts. De
Belly assays constrain mechanochemical Rac-Rho behavior, HT1080 observations
constrain free-migration morphology and dynamics, de Seze constrains local
Rho-pathway stimulation in RPE1 cells, and Jemaà/Xuan provide sparse polyploid
migration comparisons. They should not be merged into an assumed universal
cell target without evidence that the same parameterization should reproduce
all contexts. If no parameterization can jointly provide plausible target-like
cells and the required mechanochemical behavior, that is a model-adequacy
result and a reason to reconsider model structure or scope.

## Current promoted campaigns

Open the [V1–V6 joint functional reoptimization report](results/joint_functional_V1_V6_20260909/index.html).
This is the promoted joint fit of CPM mechanics, tension transfer, and all
identifiable biochemical coordinates for each of the six original regions. All
60 searches and six independent held-out assessments completed. The frozen
results, execution method, inputs, and provenance are retained together; see
the campaign [README](results/joint_functional_V1_V6_20260909/README.md).
Its compact retention policy keeps initial/final arrays and contours, archives
all search metadata, and discards intermediate trajectories and redundant raw
native outputs. The compaction audit records file hashes and retention decisions.

Open the [minimal CPM campaign report](results/debelly_minimal_campaign_20260909/index.html): emulator refitting, parameter-regime identification, and fully coupled validation of six diverse vector/module pairs. It preserves 18 resting calibrations, all 108 coupled trajectories, three-field snapshots and movies, reproducible inputs, and verification/checksums. See its [README](results/debelly_minimal_campaign_20260909/README.md) for limitations and future refinements. Temporary originals remain intact.

## Working-output policy

Use [`tmp/`](tmp/) as the default output location for all new M2072 work.
Results in `tmp/` are provisional and must remain there until the user explicitly
approves promotion elsewhere in `M2072`.

## Unrefined V1–V6 ensemble benchmark

The [provisional assessment report](tmp/ensemble_test_unrefined_V1_V6_20260914/report/index.html)
compares the original, unrefined vector/module pairs against sigma=8 morphology,
direct HT1080 motility, deBelly biochemical targets, and simple OU projections
at the literature migration protocols. All 180 simulations and the output/browser
audits completed. Perimeter-ratio terms from the stopped morphology optimizer
are excluded. See the [results and limitations](tmp/ensemble_test_unrefined_V1_V6_20260914/README.md)
and [reusable test runner](src/ensemble_test/README.md).

The [fixed Rac-noise optimization dashboard](tmp/ensemble_optimization_rac010_gauge26_20260914/report/index.html)
tracks the replacement 12-hour campaign launched on 2026-09-14 at 22:23 Eastern,
with a 1,000-core ceiling. It optimizes corrected morphology, motility and
biochemical losses around equivalent normalized V1–V6 models. The search has
26 coordinates: alpha0, beta0 and xi are fixed to one to remove three exact scale
redundancies. Nominal 10% Rac noise is fixed, and screening and full validation
run independently. The automatic deadline is 2026-09-15 at 10:23 Eastern.
See the [campaign details](tmp/ensemble_optimization_rac010_gauge26_20260914/README.md).

## Historical deBelly-only joint functional refinement

This older module is not the full morphology/movement/migration-transfer joint
screen in the adopted workflow; use the [overnight refinement module](src/overnight_refinement/README.md)
for that implementation. Historical promoted outputs below remain unchanged.

The [joint refinement module](src/joint_refinement/README.md) varies CPM mechanics,
tension transfer and biochemical parameters together against the original
biochemical scoring suite. It supports overlapping local mechanical regions.
The initial V1 pilot is provisional in
[`tmp/joint_functional_v1_20260909/`](tmp/joint_functional_v1_20260909/);
its completed report is `REPORT.md` in that directory.

The promoted V1–V6 campaign uses the
[budgeted scale controller](src/joint_refinement/scale/README.md). Its completed
HTML report and independent assessment are linked above; new refinements still
belong in `tmp/` until explicitly promoted.

The previous historical result for this model is the compact, self-contained
archive at [`archive/debelly_reduced_pipeline_20260909`](archive/debelly_reduced_pipeline_20260909/README.md).

Open the final report at
[`archive/debelly_reduced_pipeline_20260909/final_campaign/report/index.html`](archive/debelly_reduced_pipeline_20260909/final_campaign/report/index.html).

The archive preserves the selected biochemical candidates, compact emulator
fit, perimeter-to-tension calibration, all 36 final closed-loop CPM trajectories,
frozen source/model/binary inputs, report generator, documentation, and SHA-256
manifest. Superseded model variants, benchmark suites, broad screens, raw logger
output, and intermediate reports were removed on 2026-09-09.

## Prescribed-drive emulator pilot (2026-09-15)

[768 one-hour CPM training trajectories](tmp/drive_emulator_pilot_20260915/README.md)
cover 256 mechanical/drive configurations with three replicas each, no biochemistry.
Four prescribed-drive families, broad historical morphology mechanics, and
compact lossless cell masks every minute; emulator fitting is deferred.

The [emulator comparison harness](tmp/drive_emulator_models_20260915/README.md)
compares a joint model, grouped specialists, and separate output regressors
against a training-mean baseline on held-out configurations. It is queued to
train after the pilot CPM array and feature extraction complete successfully.

## Fixed ten-condition mechanical suite (2026-09-15)

[10,000 twenty-minute CPM trajectories](tmp/mechanical_suite_1000_20260915/README.md)
apply the same five protocols at two mean-drive levels to 1,000 mechanical settings,
with a 2,000-core cap. Compact cell geometry is saved every 10 seconds.
The [mechanics-only response-panel emulators](tmp/mechanical_suite_models_20260915/README.md)
are queued to train after successful simulation and extraction, with whole-setting
800/100/100 splits and joint/protocol/quantity specialist comparisons.
