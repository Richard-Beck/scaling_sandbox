# M2072 campaign lineage

## New workflow authorized 2026-09-22

The user authorized implementation and deployment of
`src/anchor_motility/`, superseding the upstream emulator and morphology-fit
workflow below for NEW campaigns. First perform and digest a numerical pilot;
adopt coarser settings only if the observed responses and motility support it.

Use 32 fixed scrambled Sobol mechanical anchors: coupling gain 10–250,
volume strength 0.5–2, surface strength 1–6 (log coordinates), and surface
target 1–1.20 (linear). Identify a causal perimeter-response model using native
CPM perturbations, then optimize the circular PDE without CPM calls, retaining
exactly one best finite biochemical seed per each of all 32 anchors. The perimeter-to-
tension coupling law remains explicit. No adaptive anchor expansion.

Jointly refine biochemical and mechanical parameters against literature-based
speed/MSD regimes and the existing relaxed deBelly transfer criterion. Use the
frozen evaluator in `tmp/joint_original6_overnight_20260918/snapshot/score.py`
as the authority; remove morphology matching, retain its response directions,
normalization, weights and fragmentation term. Keep relative Rac noise 0.10.
Long validation trajectories are capped at six hours, with fixed one-hour
burn-in and five-hour movement measurement; short forced branches share the
saved state. The joint search now has two queues: single-cell 60-minute
screens (80 minutes for hourly speed targets), and three-cell fresh-seed
six-hour validations. Split allocated CPU-hour budgets 50:50; initially allow
one validation per roughly 20 screens, then use measured costs. Rank a retained
backlog with fair anchor service; only full validations update up to three
separated parents per stream. Propose asynchronously across 32 balanced streams
using small/medium/broad steps and biochemical, mechanical/coupling, or joint
parameter blocks; occasionally return to the original anchor seed.
The authorized resource controls are 2,000 concurrent cores and 20,000 total
core-hours (2026-09-22), adjustable independently. No 640-screen cap.
Use independent asynchronous dispatch under the shared core ceiling, with
unused capacity lent between queues and no batch-completion barrier.
Use xxlarge QoS for these worker arrays (small is capped at 475 cores/user).
Carry forward existing work and the running validation; retry the four false
checkpoint failures with identical parameters/seeds after correcting centroid
comparison modulo the 600-site periodic domain. Keep field and mask checks strict.
Retention authorized 2026-09-22: retain full outputs only for the current
three validated parents per stream. Displaced historical parents retain only
seeds, complete parameter dictionaries and scores. Other completed candidates
retain only records needed by the active queues; no dense trajectories, masks,
checkpoints, XML or logs after consumption. Never remove active/unconsumed
outputs. Retain frozen scientific code and shared inputs. Background cleanup
uses atomically published dispatcher state and records reclaimed storage.
Langevin-derived targets are conditional compatibility curves, not independent
data or confidence intervals. Physical distance calibration is explicit.

Preserve the original six as biochemical multistart provenance, not as the new
mechanical anchors. Preserve all historical campaigns. Missing/failed execution must be reported. Perimeter-model fit quality and
biochemical compatibility checks are diagnostics, not anchor exclusion gates.
There is no serial native-pilot release gate. Fresh-seed validation remains required. See the new module README for
finite deployment budgets, exact protocols and numerical-pilot evidence.

## Scratch retention decision — 2026-09-19

The user authorized removal of obsolete scratch simulations/models after
verification of `tmp/workflow_provenance_20260919/`. Consult its README,
retained-file manifest and cleanup status for feature-selection, pair-metric,
movement/transfer scoring and comparison provenance. Historical full-run paths
below can be retired; this does not erase or change their documented lineage.
Keep `tmp/joint_original6_overnight_20260918/` untouched by this cleanup. Small
original-path compatibility islands (experimental masks, reference/normalization
inputs and setup snapshots) remain; do not remove them merely because their
parent scratch campaigns were retired. New workflow code should consult the
retention map rather than assume all old development directories still exist.

## Current intended workflow — agreed 2026-09-19

This supersedes the 2026-09-18 plan that required a mechanical sweep/emulator
and a 48-candidate intermediate population. The adopted workflow is:

1. Separately optimize the biochemical component against deBelly assays.
2. Develop/use a tension-response emulator (the minimal CPM perimeter emulator
   plus perimeter-to-tension coupling in the existing implementation). Use it
   in a coupled tension/diversity screen to identify several distinct coupled
   parameter sets with diverse mechanics and predicted tension responses
   reasonably compatible with the biochemical target responses.
3. Use those coupled parameter sets directly as seeds for full native-CPM
   joint refinement. Adjust biochemical and mechanical/coupling parameters
   together against morphology (pooled pair loss), baseline movement
   (speed/MSD), and migration transfer. Migration transfer replaces the older
   downstream deBelly-only joint objective, not the upstream biochemical fit.
4. Independently validate selected leaders using fresh seeds and greater
   replication before treating them as a qualified ensemble.

The emulator selects seeds; it is not assumed valid after joint refinement
changes biochemistry or coupling. Mechanical parameter diversity does not
guarantee morphological diversity or adequacy; native simulations establish
those properties. The large intermediate mechanical sweep, eight-hour
morphology/movement emulator, 48-candidate pool, and emulator-guided local
screen are optional development branches, not required workflow stages.

The user accepted the joint-refinement approach after the original-six
`tmp/joint_original6_overnight_20260918/` screen. Its current equal-weight,
frozen-scaled loss was considered tolerable: this workflow decision does not
authorize changing weights, adding transfer thresholds, restarting jobs, or
promoting search winners as independently validated results. Six seeds,
960 proposals/generation, and three baseline cells are campaign choices, not
universal requirements of the workflow.

Do not silently insert the older deBelly-only joint reoptimization between seed
selection and the new full joint screen. Original and jointly refined V1–V6
have the same labels but different parameters; labels do not establish lineage.
"Original/unrefined" here means before coupled joint refinement, not before
the upstream biochemical optimization.

## Authoritative original inputs

The promoted original campaign is
`results/debelly_minimal_campaign_20260909/`:

- `emulator/`: minimal CPM perimeter emulator and training provenance.
- `tension_screen_20260909/`: original screen and diversity selection.
- `coupled_diversity_20260909/panel.json`: original six selected pairs,
  including mechanical parameters, biochemical dictionaries and fitted kP.
- `coupled_diversity_20260909/resting_reference.json`: their associated native
  resting-perimeter calibration; consult the frozen builder/design for how
  panel-level kP/tau and biochemical defaults are combined.

Before submitting a new scratch campaign, record source paths and hashes and
verify the actual complete input parameter dictionaries against the intended
original panel, including coupling, tension parameters and calibration. Do not
copy an existing sweep's `candidates.json` or a builder's default source without
checking it. Freeze the resolved inputs in the new campaign manifest.

## Existing development exception

The existing `tmp/biomechanical_8h_V1_V6_20260915/` and
`tmp/biomechanical_emulator_8h_pair64_20260917/` actually inherit jointly refined
modules from `results/joint_functional_V1_V6_20260909/`, as do the resulting
`tmp/initial_population48_1h_20260917/` and
`tmp/migration_transfer48_20260918/` candidates.

The user accepts this existing branch for development. Preserve its true
provenance and results; do not relabel it as an original-panel campaign. This
exception does not make refined modules the default for future scratch runs.
An emulator conditional on refined biochemistry/coupling is not a validated
predictor for the original modules.

The 2026-09-18 concern about the 24,990-run eight-hour sweep is resolved at the
workflow level by making that branch optional. Preserve its historical inputs
and outcomes; do not reinterpret it as a required precursor to joint refinement.

## Live report and size diagnostics — 2026-09-22

The user reserved 2,000 core-hours for WGD size diagnostics WITHIN the existing
20,000 core-hour total and shared 2,000-core ceiling. Subtract the existing
48-hour controller reservation; screening and validation each receive 8,976
core-hours. Diagnostics use at most 198 worker cores (three per size job),
with unused capacity available to refinement. Pending jobs count against the
ceiling and wall-time reservations prevent budget oversubscription.

For each current retained parent, reuse the existing 1x validation and run
0.5x and 2x volume with the same three validation seeds and six-hour protocol
(one-hour burn-in). Follow the established volume^(1/3) geometry convention:
area target scales as s^2, volume strength as 1/s^2, initial radius and resting
perimeter as s; other mechanical coefficients, biochemical parameters,
noise amplitudes and numerical settings stay fixed. No size-specific fit or
resting calibration; size results never affect refinement scores. These are
size-only diagnostics, not additional gene-dosage perturbations.

Keep compact movement scalars and a fixed first-seed 50–60 minute mask excerpt
for diagnostics of current parents. Remove these and derived report assets
when their parent is displaced, after active jobs finish. Historical progress
retains scalar losses and compute coordinates, not raw trajectories.

The self-contained live report defaults to the lowest-loss current parent and
provides anchor/parent selectors, best validated loss against cumulative charged
completed-worker core-hours, tracks, speed, MSD, size ratios, GIF excerpts and
existing relaxed deBelly scores. Refresh report assets every minute and browser
every two minutes. Different target-family losses remain campaign objectives,
not comparable likelihoods. Pending sizes have explicit status, not predictions.
