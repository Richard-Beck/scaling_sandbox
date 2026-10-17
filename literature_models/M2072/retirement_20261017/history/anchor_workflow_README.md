# Retention policy authorized 2026-09-22

Full outputs persist only for current validated parents (up to three per stream,
96 total). Displaced parents retain a `result.json` with exactly `seeds`,
`parameters`, and `scores`. Other completed candidates lose intermediate outputs;
a backlogged screen temporarily retains request, calibrated candidate and loss
for promotion. Once no queue needs it, that candidate directory is removed.
Frozen code, shared scientific inputs, controller bookkeeping and cleanup ledger
remain. No morphology arrays or raw diagnostic examples are retained separately.

`retention.py` runs in the dispatcher's background thread. Only atomically
published, consumed results whose workers have exited are eligible. Active and
current-parent paths are protected. Parent history is recorded during selection,
including parents displaced in the same dispatcher tick. Replacement compact
records are written atomically before deleting raw files. Repeated/interrupted
cleanup is safe. Retired task manifests and worker stdout are removed once no
worker needs them. Explicitly allowlisted campaign output roots prevent touching
other campaigns. `retention_summary.json` reports reclaimed allocated bytes;
`retention_ledger.json` records candidate dispositions. Running jobs can still
use substantial temporary disk space; the ~12 GB parent collection is not a
peak campaign disk cap.

The active dispatcher uses a versioned controller revision; the scientific
worker snapshot and running simulations are unchanged.

---

# Active asynchronous deployment (2026-09-22 revision)

This section supersedes the batch dispatcher below. `async_queue.py` imports
completed screens and biochemical seeds from `anchor_motility_two_queue_20260922`,
accounts for its still-running validation, and dispatches independent queues.

User-authorized controls: **2,000 concurrent cores**, **20,000 total core-hours**,
50:50 per-queue allowances. Edit `control.json` atomically to change either limit.
The controller polls SLURM and individual result files every 15 seconds. Pending
and running workers both consume concurrency commitments. One core is included
for the dispatcher; worker QoS is `xxlarge` (3,000-core user limit). Unused queue
capacity is lent to the other queue. Lowering limits drains existing work without
cancelling simulations. `pause` stops new submissions; `stop` drains and exits.
No array throttle or generation barrier; the scheduler still controls allocation.

Dispatch reserves the full requested wall-time cost (screens: 1 core × 2 hours;
validation: 3 cores × 6 hours). A reservation is replaced by measured worker
core-hours when its result arrives. This bounds worker expenditure conservatively
and can reduce concurrency near the budget limit. Reserve 48 core-hours within
the total for up to four 12-hour controller allocations; unused controller
reservation is not spent. The remaining worker allowance is split equally. Screening earns validation credit using completed plus in-flight
estimated cost, initially 20:1 and then measured costs. Both cumulative allowances
remain equal; audit promotions consume the validation allowance. Unused budget is
not a reason to relax scientific screening or fabricate validation candidates.

Each of 32 streams has fixed normalization bounds and a fixed target family.
Proposals are allocated to the least-proposed stream. Parents are the best three
separated finite, eligible validated candidates (minimum normalized RMS parameter
distance 0.05); best loss always retained. Before validation, use the anchor seed;
afterwards choose it with probability 0.15, otherwise choose a validated parent
uniformly. Mutation blocks: biochemical 30%, mechanical/coupling 30%, joint 40%.
Normalized Gaussian scales: 0.025 / 0.06 / 0.15 with probabilities 50% / 35% / 15%.
Reflect at bounds. Log coordinates for positive parameters; linear surface target.
Each request records its exact parent, block, scale, axes and serial number.
Scientific objective, forcing and noise are unchanged from the two-queue design.

The four false checkpoint failures are resubmitted once with identical parameters
and seeds. The fix compares centroid positions modulo the 600-site periodic box;
membrane fields, other scalars, and checkpoint masks still receive strict checks.
Prior failed outputs and their compute cost remain recorded. There are no automatic
retries of other failed proposals or favorable-seed replacement.

`async_summary.json` reports current utilization and spend; `async_state.json`
records active jobs, backlog, parent sets and proposal counts. Submissions have
intent/receipt files for crash diagnosis. Missing worker outputs become explicit
failures after a propagation grace period and are charged their full reservation.
A submission with an ambiguous receipt requires reconciliation before restart.
The dispatcher checkpoints continuously and hands off to a successor before its
wall limit. Frozen manifests exclude only the intentionally mutable resource controls.

---

# Active two-queue deployment (2026-09-22 revision)

This section supersedes conflicting legacy deployment details below.
All 32 characterized mechanical anchors are retained. Each contributes its
best finite biochemical result across the six original V1–V6 multistarts.
Perimeter adequacy, biochemical sign checks and half-timestep checks are
recorded diagnostics; none excludes an anchor. No serial native pilot gate.

`two_queue.py` runs a finite campaign of up to 640 screens (20 per anchor),
with a 2000 allocated CPU-hour forecast budget split equally across screening
and validation. After each batch, validation slots use accrued screening spend
and measured validation cost (20 times screening cost until measured).
Runtime variation can cause one-batch budget overshoot. Both arrays can run
concurrently; a single dependent controller collects them before the next batch.
There is no array concurrency throttle. Each queue uses suitable small/medium QoS.

Screens: one cell, 10-minute burn-in plus 50-minute movement; hourly-cadence
speed-only targets use 70-minute movement, total 80 minutes. Rac-front and
Rho-front pulses plus matched control use existing amplitudes and 90-second
forcing, with 300-second recovery. Relative Rac noise stays at 0.10.
The single-cell front-only legacy-score adapter is preliminary triage, not a
replacement for the full relaxed deBelly criterion: its noise denominator uses
within-cell control variation only; between-cell variance is unavailable.
Wrong-sign z < -3 is
flagged, not rejected. Reject execution/numerical failure, disappearance,
persistent fragmentation (>5% average off-largest-component mass), or movement
more than a factor four outside every coherent target (speed and RMS distance).
Weak/ambiguous deBelly responses stay in the ranked backlog. No morphology fit.

Joint-MSD streams use each retained Langevin curve at 1, 3 and 10 minutes,
plus speed at the study's cadence. Jemaà streams use speed only at the original
cadence, with no invented MSD constraints. Short-lag predictions retain their
model-conditioned interpretation. Stream target families never change.

Validation: three fresh seeds, six-hour trajectory each, one-hour burn-in,
full inherited front/flank relaxed deBelly score and original speed/MSD targets.
Only these results update incumbents. Proposals perturb biochemical and
mechanical parameters around a validated incumbent, or the anchor seed until
one exists. Validation serves least-served anchor streams first, then cheap
rank; every tenth slot randomly audits an eligible unpromoted candidate within
that fairness policy. Audits draw from the same budget. Unpromoted candidates
remain in the backlog, not classified as failures. Seeds are fixed within each
queue and disjoint across queues; this is validation for optimization, not an
untouched final statistical holdout.

Status: `queue_summary.json`, full backlog/provenance: `queue_state.json`,
per-worker raw results: `screen/b*/` and `validation/b*/`. Existing response
fits and completed biochemical multistarts are reused with provenance; skipped
multistarts are run. Historical campaigns remain intact.

---

# Anchor-conditioned biochemical initialization and motility refinement

New workflow authorized on 2026-09-22. Historical campaigns remain unchanged.
The circular PDE uses a small causal perimeter-response model fitted to native
CPM perturbations; it never invokes CPM inside biochemical optimization. Joint
refinement subsequently evaluates the complete native model.

## Numerical pilot

`numerical_pilot.py` freezes three far-separated ORIGINAL V1–V6 dictionaries,
the original resting references, the latest joint campaign's evaluator and
all execution sources. Three common search seeds evaluate five configurations:

| Name | Membrane samples | Fast systems/projections | CPM sweep duration |
|---|---:|---:|---:|
| reference | 256 | 0.05 s | 0.01 s |
| membrane128 | 128 | 0.05 s | 0.01 s |
| chemistry100ms | 256 | 0.1 s | 0.01 s |
| mechanics20ms | 256 | 0.05 s | 0.02 s |
| combined | 128 | 0.1 s | 0.02 s |

The lattice remains 600×600 and volume target 4800. Tension retains its 0.25 s
update. Membrane diffusion coefficients scale by (N/256)^2 to preserve the
angular diffusion operator. Nodal Rac noise scales by sqrt(N/256), preserving
low-mode spatial white-noise power; native SDE timestep treatment is retained.
These adjustments themselves are part of the numerical comparison.

Each pilot cell runs two hours, with a fixed 600-s burn-in. Five matched branches
(control, Rac/Rho front/flank) run a 90-s pulse and 300-s recovery. All boundary
states and masks are checked at restart; MCA mass is checked in every branch.
The baseline continuation is unforced. Long baseline output is scalar every
10 s and connectivity masks every 300 s; assay output remains every 5 s.

Frozen comparison: speed and RMS displacement at 1, 3, 10 and 30 minutes,
transfer score, runtime, and basic cell validity. To adopt a coarser alternative,
all three settings must have movement changes within 25% and within two panel
standard errors (5% numerical floor), transfer change <=10 points and within
two panel standard errors (2-point floor), and >=1.1× median runtime speedup.
These are engineering conventions, not experimental confidence statements.
A passing alternative requires a separate six-hour confirmation before adoption.
Otherwise retain the reference configuration. Three seeds do not establish
formal equivalence. The decision JSON is required by campaign setup.

## Fixed anchor design and response identification

32 scrambled Sobol points, seed 20260922. Coupling gain 10–250, volume strength
0.5–2, surface strength 1–6 in log coordinates; surface target 1–1.20 linearly.
No adaptive additions. Mechanics-only CPM has the same selected membrane and
CPM resolution as the coupled simulations.

At each anchor, amplitudes 0.25 and 0.75 of bounded RacHill are applied both
uniformly and with the inherited localized kappa=20 profile, from 600–1200 s,
followed by release to 1800 s. A third amplitude, 0.5, is held out and pulsed
from 600–780 s under both profiles. Each protocol has three seeds. The uniform
input probes the resting-drive dependence needed when biochemistry changes.

Fit an absolute resting perimeter and a saturating equilibrium response to
the spatial mean and standard deviation of RacHill, with separate rise/fall
time constants. The first anchor pilot showed that mean-only loading failed
to distinguish uniform from localized forcing. The contrast term was added
before campaign freeze; identification logs perimeter every second to resolve
the fast response. This is a deliberately
small phenomenological closure, not a general emulator or a mechanistic law.
All training and held-out response errors must have normalized RMSE <=0.35,
using a denominator no smaller than 1% of resting perimeter, with positive gain
and <=5% mean fragmentation. Inadequate anchors are recorded and excluded,
not repaired by an expanding model family or silently accepted.

## Circular PDE initialization

Six original biochemical dictionaries supply independent starts per adequate
anchor. Mechanical coordinates stay fixed. The archived biochemical bounds
are used, holding KT, imposed input amplitudes, noise, and numerical settings
fixed; kP and tau_T are adjustable. Each start has at most 64 proposals and
90 minutes of compute. This finite first campaign is not a convergence claim.

The causal perimeter state is updated from the current spatial RacHill mean
and contrast,
then coupled through rectified relative perimeter strain to tension. Its
resting perimeter is calculated from the candidate's homogeneous RacHill.
Drive excursions above the largest training mean (0.75) are explicitly flagged
in biochemical results; native calibration/refinement is not replaced by that
approximation. Circular MCA inventory is conserved. Additive Rac noise and rectified Hill
inputs follow the native convention; occasional negative noisy Rac is not
silently clipped or rejected merely because it is negative.

The upstream score is explicitly an ADAPTER of the frozen signaling terms:
same four distal response signs, same pulse/preinput windows, same equal-cell
control-variance calculation, same saturation at z=3. No movement term is
invented for a fixed circle. Three seeds per evaluation. A candidate is eligible
as a seed only when all four mean signed effects are positive and the score and
signs survive a half-timestep check (score difference <=10). No exact fold-change
targets are restored. This seed screen is a modeling choice; native assessment
remains authoritative. Retain up to two starts separated by >=0.05 RMS distance
in normalized biochemical/coupling coordinates.

## Native objective and movement targets

Each candidate receives three uniform-drive native resting calibrations before
its coupled evaluation. This prevents the original panel's fixed Prest values
from being mistaken for calibration of newly fitted biochemical parameters.
Calibrations and complete resolved dictionaries are retained.

Three six-hour cells: one fixed hour of burn-in and five hours of unforced
movement. Dense records only cover the final ten minutes of burn-in and the
short forced branches. Every forced branch starts at the SAME one-hour state;
absolute checkpoint times are retained, while the scoring adapter shifts the
assay clock to the frozen 0–990 s convention. Stimulation times do not depend
on favorable polarity or stationarity. Ambiguous fronts retain fixed-axis labels.
All branches retain relative Rac noise 0.10, inherited input profiles/amplitudes,
and matched reseeding. The movement continuation contains no imposed forcing.

Native relaxed transfer uses the exact frozen score functions:
75% signaling + 25% movement, with the original 25% fragmentation penalty.
The shape penalty is zero; no feature-reference or morphology matching is run.
The original 33 nonoverlapping 90-s movement-noise windows are translated to
the new post-burn-in interval. Signal-noise windows and boundary-length-weighted
distal observations are unchanged. Fragmentation also controls basic eligibility.

Movement uses 10, 30, 60 and 120 minute lags, equal cell weights and periodic
centroid unwrapping. Physical units use an explicit um_per_site (initial proposal
0.2, inherited lengthscale) and 3600 seconds/hour. This is not fitted to target
speeds. Current target library contains 32 coherent curves for each Xuan and
Mallin group, reconstructed from joint compatible grid points. Mallin uses the
existing documented finite-track sensitivity, including its weakly constrained
parental regime. Jemaà conditions retain speed-only intervals at the reported
cadence; they do not acquire invented persistence/MSD constraints. Hosny steering
is not imposed on unforced trajectories. Target uncertainty is model-conditioned
and assumes stationarity, not an empirical confidence interval.

Every search stream retains one target family. This prevents switching to a
speed-only family to evade its MSD objective. For a joint family, score the
closest coherent curve across ALL lags: equal weight log2 speed discrepancy and
mean squared log2 RMS-displacement discrepancy, followed by a square root.
Native scalarization weights movement and (1-transfer/100) equally. Raw scores,
matched target, per-cell summaries and all branch components remain available.

## Deployment, bounds and retention

Use existing `scaling_sandbox`; no environment installation. Compute uses SLURM.
Prepare with `campaign.py freeze --out NEW --pilot DIGESTED_PILOT`, then
`campaign.py deploy --out NEW`. Setup refuses to overwrite. Sources, original
dictionaries, targets and numerical decision are hashed in the frozen package.
The first anchor is a native deployment pilot before releasing the other 31.

Dependent stage callbacks check completeness before submitting the next stage:
response identification → PDE fitting → initial native evaluation → three
proposal generations → independent assessment. At most 64 seed streams and
256 native search evaluations; native workers request three CPUs and 8 GiB,
medium QoS, maximum 12 hours. Response/PDE tasks use one CPU, 3 GiB and small
QoS. Arrays have no concurrency throttle. An additional conservative 2000
CPU-hour native budget, based on three allocated CPUs times elapsed evaluation
time, prevents starting another generation without a validation reserve.
This is a forecast admission limit, not a scheduler-enforced aggregate CPU cap;
the evaluation count and individual wall limits are the hard bounds.

Validation selects at most one leader from each of the eight target families
and evaluates six fresh cells, in two three-cell waves. Scores are independent
assessments, not automatic qualification of an ensemble. All native results,
failures, masks, compressed trajectories and checkpoint evidence are retained
for this finite first campaign; temporary native logger/TIFF files are removed
by the node-local runner after successful extraction. No historical files are
cleanup targets. Current reports refresh at stage boundaries, not every second.

`steering.json` stop/pause prevents the next stage; submitted jobs drain. A
stopped boundary records its pending stage and requires an explicit resume
command using that stage. Missing task files stop the chain with a diagnosis.
Individual failed native proposals remain visible and ineligible, never retried
under a favorable seed. A wholly unsuccessful initial population stops search.

Run scientific unit checks with:

    /home/4473331/.conda/envs/scaling_sandbox/bin/python -m unittest discover -s literature_models/M2072/src/anchor_motility -p 'test_*.py' -v

Native and circular-PDE execution checks are run separately on SLURM before
deployment. See `tmp/anchor_numerical_pilot_20260922/` and the implementation-check
directories for observed evidence, including failed preliminary checks and fixes.

## Deployed campaign

`tmp/anchor_motility_20260922/` was released after all 45 numerical pilot tasks
completed. None of the four alternatives passed every predeclared comparison,
so the original 256-sample, 0.05-s fast-system, 0.01-s CPM configuration is
retained. Membrane-only and chemistry-only coarsening accelerated by 1.77x and
1.82x; chemistry-only preserved all transfer comparisons but failed the movement
variability comparison at one setting. CPM-only coarsening gave only 1.03x speedup
and changed movement; the combined configuration gave 3.51x but changed movement.

The verified first anchor is reused byte-for-byte after source, numerical and
anchor-dictionary checks. The remaining 31 were submitted as array 20510906;
callback 20510907 admits the biochemical stage after checking all response fits.
The production initial population likewise releases one complete six-hour
evaluation first, then admits the rest only after its technical success.
Deployment release job: 20510904; initial controller: 20510905. All later job
IDs and stage outcomes are recorded in the campaign's state.json.

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
existing relaxed deBelly scores. Refresh report assets every minute while the dispatcher runs; the browser
loads the latest generated report only when opened or manually refreshed. Different target-family losses remain campaign objectives,
not comparable likelihoods. Pending sizes have explicit status, not predictions.

## Stopped-campaign report reconstruction — 2026-09-24

`python live_report.py --out OUT --completed-snapshot` reads completed worker
results left after a dispatcher failure into an in-memory copy of the state.
It does not change the dispatcher state, submit jobs or restart the campaign.
The resulting static HTML summarizes all finished screens, validations and
size diagnostics. Its generation timestamp and last worker-result time are
shown separately. The report remains static until explicitly regenerated.

Progress is faceted by reference and polyploid targets for Jemaà 2017,
Jemaà 2023, Xuan 2018 and Mallin 2023. Its x coordinate is unique short-screen
proposals submitted per anchor when each validation was launched; the four
identical bug retries are excluded. This count is reconstructed exactly from
submission receipts and the original two proposals per anchor. The curve is retrospectively ordered by the proposal count at validation
submission, since validations finish out of submission order. Circles mark
new best validations and short ticks mark non-improving validations.

The WGD summary uses the current lowest-loss validated parent in each anchor.
The MSD exponent is a log-log least-squares slope of the equal-cell mean MSD
at lags 10–120 minutes. Points show 1x/0.5x and 2x/1x exponent ratios for reference and polyploid
leaders, grouped by study, with joint loss encoded by color. A missing size run for that exact leader produces
no point. Earlier parents are not substituted and there is no imputation.

## Manual completion of missing WGD size runs — 2026-09-24

After the dispatcher stopped, use the finished validation results in a read-only
snapshot to select each anchor's lowest-loss validated parent. The five anchors
without a completed validation have no eligible WGD parent. Schedule only
missing 0.5x and 2x runs for current leaders, all in one unthrottled SLURM array
with three cores per task and the established six-hour walltime. Existing
completed size runs are reused. This array does not restart refinement or
change the 20,000-core-hour campaign ceiling. Its reserved core-hours count
against the previously authorized 2,000-core-hour WGD diagnostic allowance.

The manual registry stays outside dispatcher state. A dependent one-core report
job runs after all array tasks finish, incorporates available results into an
in-memory state copy, republishes the static report and records any missing or
failed size tasks. No browser auto-refresh is used. The scientific scaling
protocol, numerical settings and three validation seeds are unchanged.
