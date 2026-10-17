# Mechanism and measurement component: manuscript intake

Evidence first checked on 17 October 2026 at approximately 20:29 UTC, with
production updates below. This intake is
for the manuscript-outline publication task. It does not supersede the live
checkpoint producer or the author's figure-review workflow. The first outline
release should integrate the delivered checkpoint substantially as presented;
the subsequent release should incorporate a separate coherence review.

## Governing direction and manuscript role

The controlling scientific direction is [MANUSCRIPT_VISION.md](../MANUSCRIPT_VISION.md).
Enlargement changes polarity organization and the response to an opposing
stimulus; those changes motivate, but do not demonstrate, greater migration
persistence in WGD-associated enlarged cells. Mechanistic observations and
corresponding proposed cellular measurements belong together. The new SUM-159
chemically induced WGD+/− motility results require a separate prominent major
results section and figure placeholder. The stochastic model figure must not
fill that experimental gap.

The [current manuscript instructions](../polarity1d/manuscript/AGENTS.md),
[revised three-figure plan](../current%20manuscript%20figs/figure_revision_plan_20261017.md)
and [completion requirements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/VISION_CHECKPOINT_REQUIREMENTS.md)
record explicit author clarification: all five established linear models;
supported Wave-pinning/Goryachev/Holmes rings; gain/recovery **and** substantial
finite competition in both; and a separate replicated noisy-ring persistence
figure. Ring results cannot substitute for interval results. Corrected deBelly
is eligible in the linear comparison; its ring development remains a sidecar.
Working checkpoint figure numbers do not assign final manuscript numbers.

The author's scope is supported by the records above, but these are repository
records of clarification rather than original chat transcripts. The outline
should not describe an agent review or completion audit as user approval.

## Delivery and authoritative production evidence

The expected exact deliverable is
`current manscript figs/vision_aligned_checkpoint.md` (the spelling is
intentional); the canonical `current manuscript figs` collection should link to
it. Neither directory contained that final file when first inspected. A
[composition plan](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_aligned_checkpoint_20261017/supporting/metadata/COMPOSITION_PLAN.md)
and two actual assembled ring review sheets already existed. Do not replace
the final package with that plan or an older checkpoint while current work is
demonstrably live.

| Component | Authoritative evidence at intake | Publication implication |
| --- | --- | --- |
| Linear finite competition | `data/runs/vision_linear_20261017/finite/results/status.json`: 127/127 jobs, finished, verification passed; 97 primary branches plus 30 duration anchors. | Actual completed new evidence, available for integration. |
| Full linear recovery | SLURM array 21126305: original slot 78, actual job 21126371, was RUNNING on hpctpa3pc0029 from 16:18:03 EDT (20:18:03 UTC), with a 30-minute limit. Other new slots and 15 pilot imports existed. | A status file saying initialized is stale, not proof the producer stopped. |
| Recovery early-window amendment | `recovery/early_window_amendment/results/status.json`: job 21126520 completed; passed, 20:27:05 UTC. Two Spring sizes L₀0.8/0.9 retain the same already fixed measurement clock and omit unused later dynamics. | Check the final full-series collector records this amendment; do not silently call the unresolved original trajectory a successful complete integration. |
| Ring recovery | `data/runs/vision_ring_recovery_20261017/results/status.json`: finished, 702 states (50 basal, 652 challenged), 2106 amplitude rows, all numerical gates; independent read-only native-data audit. | Completed data, not a placeholder. |
| Ring finite/vector analysis | `vision_ring_moments_20261017/status.json`: complete; 652 histories, 2608 exact readouts; separate duration campaign complete. | Actual full-grid measurements and finite-duration checks exist. |
| Noisy-ring pilot | 30/30 pilot cases complete with raw numerical checks and actual-result review. | One seed per condition establishes feasibility/regime, not an ensemble trend. |
| Main noisy-ring ensemble | `vision_noise_main_20261017/status.json`: SLURM handle 21126522, submitted 20:27:54 UTC; `squeue` independently confirmed its array tasks RUNNING from 16:27:54/55 and 16:28:10 EDT. | This is a verified live wait; continue other outline preparation while it runs. |

The main ensemble is prospectively frozen at 12 independent seeds × three
models × C20/C40/C60 × two strengths (216 noisy cases), plus 27 deterministic
mean/no-input controls, horizon 32000 and maximum lag 8000 in native time.
No runtime ETA is established by the status record. The pilot's 0.79–119.59 s
worker range at horizon 4000 suggests minutes for many full workers but is not
a guarantee of full production or analysis time. Poll the same job handle
before deciding whether production has stopped. Scheduler timestamps here are
EDT; machine-readable campaign timestamps use UTC.

**20:32 UTC update:** the full linear recovery collector now reports finished,
all 97 cases complete and all fixed clocks qualified (generated 20:31:09 UTC).
SLURM independently reports original task 21126305_78 COMPLETED, exit 0:0,
elapsed 11:43. The main noise handle 21126522 still has 168 RUNNING and
75 COMPLETED tasks, with no reported failures; its live status agrees.
The final checkpoint remains absent. The initialized linear record described
in the intake table is historical and no longer the current state.

**20:48 UTC update:** noisy producer 21126522 finished all 243 jobs with exit
0:0, and the separate native/archive/noise verifier 21126765 finished with
all 243 cases passing. Frozen-source verification covers 22 executable hashes,
11 scientific plan copies and nine native inputs. Worst pool relative error
is 8.52 × 10⁻¹²; every native concentration remains positive. The group occupies
2,791,211,008 allocated bytes (2.60 GiB). Statistics handle 21126773 is now
RUNNING; final statistics, noisy figure and the final checkpoint remain
pending. No native integration is restarted. The later numerical-correlation
amendment preserves the frozen estimator and original source, correcting
tiny-tail cancellation with positive suffix energy sums and direct products.

## Proposed integration legend: all-five linear component

**Working message.** Enlargement increases incumbent excess-active material in
all five retained models, through model-dependent concentration and spatial
organization. The measured initial response and small-disturbance recovery
describe specified local functions; substantial finite stimulation reveals
whether challenger activation persists after withdrawal. Neither the common
organization change nor one local coefficient determines a universal
polarity-replacement law.

**Legend for final assembly, to be matched to its actual panel lettering.**
Incumbent organization uses the exact qualifying nominal relative-20% basal
checkpoints in the established linear models: Wave-pinning, Goryachev,
corrected independent-endpoint deBelly, Holmes 3 and Spring. Excess-active
amount integrates activity above the basal background in the original
half-contrast incumbent core; peak concentration and core width remain separate
observables. Spring uses its measured physical domain and fixed material
coordinates. Native concentration and time units are not cross-model
calibrations. Initial local-input gain measures the derivative of native
primary-activity response to the declared opposing input, rather than a
spontaneous event rate.

For the common recovery assay, a smooth compact taper of the basal incumbent
excess-active profile specifies a small removal direction. The primary
amplitude is 0.5%, with 0.1% and 1% sensitivity checks. Removed active material
is transferred pointwise to the corresponding explicit inactive pool in the
four models that possess one; deBelly changes Rac alone. Every native field
then evolves without a challenger cue, against an independently evolved sham.
In the fixed original core, the repaired fraction is the reduction in the
sham-minus-probe active deficit divided by its initial value; division by the
fixed model-specific measurement window gives the plotted finite-window
coefficient. Fixed clocks are 0.1778279410 for Wave-pinning/Holmes,
0.00562341325 for Goryachev, 5 s for deBelly and 0.01778279410 for Spring.
This is a response along one specified disturbance, not an eigenvalue or
generic stability measure. The final collector must retain qualified,
unresolved and failed statuses rather than choose clocks by size.

The standardized opposing cue uses endpoint raised-cosine support 3 (Spring
0.2), native peaks 0.1/6.4/3.6/1/4 in the model order above, duration 1000
native units (Spring 50), and 13000 cue-free units after withdrawal. Support,
peak and duration are fixed across sizes within each model. The top finite-
response row shows the state at withdrawal; the lower row shows the later
cue-free state. Signed excess in fixed incumbent/challenger cores is
normalized by the initial incumbent amount; negative excess denotes activity
below the original background, not negative concentration. Raw integrals,
native fields, spatial patterns and independent orientation/strength/field-
tail gates remain linked. Halved/doubled durations at prespecified size anchors
are separate sensitivities, not a threshold search or a cross-model dose match.

**Narrative points supported now.** Longest retained Wave-pinning L39 and
Holmes L50 retain two above-basal-half-contrast regions after withdrawal;
their incumbent/challenger core-excess ratios are approximately 0.581/0.451
and 0.766/0.248. Goryachev L39 activates a challenger during stimulation and
returns to the incumbent afterward. All retained deBelly and Spring primary
branches pass the original orientation gate, but that is distinct from full
strength restoration or field stationarity. At deBelly L18 the challenger
core excess falls from 7.99 times the initial incumbent at withdrawal to about
1.013 after 13000 cue-free units. Do not equate a strong cue-on signal with a
retained state, or a fixed-region amount with uniquely tracked site identity.

**Corresponding measurements.** Image incumbent physical extent, peak and
integrated excess reporter signal; interpret molecular amount only with
calibration. Under characterized local input, separate early challenger gain,
incumbent local loss/rebound, cue-on dynamics and post-withdrawal site survival.
Record spatial support, input amplitude, perturbation region and observation
window. A local signal can leave a fixed region by redistribution, so also
track the site itself and the surrounding boundary.

Completed source: [linear finite results](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_linear_20261017/supporting/metadata/FINITE_RESULTS.md).
The full recovery legend and final numeric conclusions require the pending
collector's actual report, not the preceding 15-state pilot.

**20:34 UTC update:** the completed [all-five report](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_linear_20261017/RESULTS.md)
now supplies every original full-span archive and the assembled main linear
figure. Near-doubling common-recovery folds are 1.02562 (Wave-pinning), 3.43942
(Goryachev), 1.10175 (deBelly), 1.00871 (Holmes) and 1.608 (Spring). The
original slow Spring trajectory completed; the early-window amendment remains
additional validation and is **unused** by the main collection. The final main
assembly displays one common final-stage regional scale and exact native
histories, resolving the earlier concern about using the full 0–8 cue-on range
for the small post-cue effects. The earlier pending-collector statement above
records the intake chronology and is now superseded.

## Proposed integration legend: supported-ring component

**Working message.** Supported ring versions reproduce the established local
organization and response closely at matched C=2L, then allow challenger
placement to vary inside one complete basal state. Initial input gain varies
little with separation; later activity patterns differ markedly. Conditional
recovery remains positive even while the actual early incumbent core loses
activity, separating a restoring contribution from net maintenance.

**Legend, two continued sheets.** A, independently prepared native interval
L15 and ring C30 profiles, with the interval reflected for display. All 13
matched C=2L pairs retain near-doubled ring width/core amount, approximately
unchanged peaks and matched gain, and common-assay recovery agreement within
1.02%. States and native meshes are independently prepared; this is neither
exact identity nor a continuum-convergence result. Matching an interval
endpoint support 3 requires reflected ring full support 6; finite ring
challenges below use full support 3 and are a different input.

B, initial cue-associated input gain versus circumference and shorter-route
fraction d/C, divided by each model's own C30 antipodal value, on one shared
0–3 scale. C, the same active-deactivation assay applied exactly 10 native
units after cue onset, with cue continuing identically in probe and sham,
divided by the corresponding basal repaired fraction at that circumference.
The original basal perturbation direction and original core remain frozen;
all fields evolve. Fixed clocks are identical to those used for the chemical
interval models. All 867 planned cells remain represented: 652 completed
challenges, 198 prespecified cue/core overlap exclusions and 17 unavailable
Holmes C80 preparations. Crosses and dark gray cells distinguish these statuses.

D, local-pattern geography after width-3 raised-cosine stimulation for 1000
units and 13000 cue-free units. Native peaks are 0.1, 6.4 and 1 for
Wave-pinning, Goryachev and Holmes. Geographic classes distinguish a sole
old-associated maximum with unresolved identity, a cue-associated peak,
distinct old/cue-associated peaks, one spanning region, and shifted/unresolved
patterns. They are finite observations rather than attractor classifications.
Exact native Holmes histories highlighted by the map labels show C30 d12
cue-associated survival, C40 d12 merging/bridging and C40 d16 distinct sites.
Only the two C40 cases change cue placement in the identical complete basal
state. The C30-to-C40 comparison also changes material, architecture, both
routes and cue angle. Kymographs show the first 4000 units; endpoint maps use
14000 units from onset. The original peak remains inside the periodic display
cut; gray padding lies outside the C30 domain.

The activity moment is the first circular harmonic divided by its own initial
complex moment; initially it is (1,0). Display its two components and magnitude
together. Moment-space trajectories are activity-state trajectories, not cell
positions. Direction is unavailable below relative magnitude 0.1, with
0.05/0.2 sensitivities retained. The vector is never altered by the gate.
Two opposed strong activity sites can cancel while total activity remains
substantial. Native local profiles, total activity and both first/second
harmonic ratios distinguish this from homogeneous loss of activity. A second
harmonic is not itself a unique site count.

**Narrative points.** Within-size gain ranges are below 0.09% across the
qualified placements, while finite site outcomes depend strongly on placement.
All 652 challenged conditional repaired fractions are positive; every actual
sham core loses active material over its measurement window. Relative to the
same-size basal response, challenged recovery spans 0.469–1.000 in
Wave-pinning, 0.798–1.000 in Goryachev and 0.921–0.999 in Holmes. Basal
recovery itself changes with size, so normalization must not hide those
curves. Nine bounded half/double-duration anchors retain their local endpoint
class, while vector magnitude/angle can change. Four of 652 local classes
still change between 8000 and 13000 cue-free units. The chosen final time is
a declared window, not a universal settling certificate.

**Corresponding measurements.** Follow initial and tracked activity sites;
record extinction, replacement, merging and two-site lifetimes against size,
physical center separation and edge gaps. In closed-boundary geometries record
both routes. Estimate a whole-boundary reporter moment and magnitude alongside
local activity, and image local protrusion separately to test its relation to
activity. Moving a cue changes both routes and local basal context; these
results do not isolate transport distance from organization.

Source: [complete ring report](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_ring_20261017/RESULTS.md).
Both assembled ring sheets were visually inspected at intake; they convey
their stated comparisons and expose unavailable conditions.

## Prospective integration legend: replicated noisy-ring component

**Status.** Keep this section prospective until actual ensemble output and
its fresh result critique arrive. The pilot is not a substitute. Final
legend must include actual coverage, uncertainty, outcome patterns and the
reproducible example-selection rule.

**Methods text ready for integration.** Native periodic Wave-pinning,
Goryachev and Holmes basals at C20/C40/C60 receive nonnegative lognormal
activation input through their original cue channel. Local Gaussian field
variance, physical spatial correlation width 1.5, temporal correlation time
50 and native mean fractions 0.001/0.01 of the established cue peak remain
fixed within each model. Stationary OU knots spaced by 5 are interpolated
with analytic marginal-variance normalization; the resulting law is a
declared interpolated process, not an exact continuous OU path. No realization
is demeaned, clipped or normalized around its own mean. Larger circumferences
therefore contain more fluctuating boundary and have different whole-ring
input modes, despite matched local statistics. Deterministic uniform-mean and
zero-input controls remain visible.

Twelve independent realizations per model/size/strength are followed for
32000 native units. Report initial-axis memory separately from uncentered
vector lag correlation. The latter pools within-seed lag products and
normalizes by both endpoint energies; seed resampling provides uncertainty,
not independent time origins. The signed restricted correlation integral
through common lag 8000 is an observation-window descriptor, not an
extrapolated fitted persistence time. Full-, first-half and last-half windows
expose continuing adaptation. Show vector magnitude, total activity,
harmonics and native site histories; direction-only results require adequate
magnitude and their eligible-pair fraction. Physical site displacement must
accompany angular summaries when identifiable. At fixed physical displacement,
angle changes less on a larger ring. Finite-cadence endpoint eligibility does
not prove uninterrupted site identity.

The pilot already shows a specific caution: high-mean Goryachev C20 loses
polarity even in its deterministic mean-input control. That forced-state
regime cannot be interpreted as a pure fluctuation-induced directional-memory
comparison around the original polarized branch. It must remain labelled in
the ensemble and the manuscript. Wave-pinning's pilot at high input also
shows stronger angular memory at C60 despite larger eligible physical site
displacement than at C20. This is one realization, not an ensemble conclusion;
it demonstrates why angular memory and physical drift must remain separate.

**Corresponding measurements.** Whole-boundary reporter correlation,
resultant magnitude, weak-moment intervals and local site lifetimes can test
activity-memory predictions. In irregular cells declare boundary registration,
physical arclength weighting and the vector embedding; equal ring angle must
not become equal physical distance without justification. Observationally
matched detected excursions estimate subsequent establishment, not equal-input
gain. Compare spontaneous event rates per cell-time and accessible-boundary-
length–time. A measured relationship between reporter-polarity memory and
migration displacement/directness is an experimental test, not a simulated
velocity law.

Sources: [frozen noise design](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_20261017/DESIGN.md),
[actual pilot results](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_pilot_20261017/RESULTS.md),
[experimental-readout note](../current%20manuscript%20figs/vision_checkpoint_experimental_readouts_20261017.md).

## Concrete proposed measurements, with inference limits

These are candidate follow-up measurements, not experiments already agreed
with the collaborating laboratory. Keep them alongside their motivating
mechanistic panel in the outline rather than presenting an unattached assay
wish list.

| Mechanistic observation | Candidate measurement and comparison | What it could establish, and what remains open |
| --- | --- | --- |
| More incumbent excess-active material through different width/peak changes | Boundary-resolved activity reporter: physical half-contrast extent, peak, background and excess integrated signal versus measured cell size. Retain the whole boundary and segmented shape. | Distinguishes broadening from increased concentration-like signal. Integrated fluorescence requires calibration before being called molecular amount; reporter expression, thickness and geometry can confound it. |
| Lower equal-input opposing gain in selected models | Characterize a local activation input's amplitude, physical support and duration, then measure the early incremental challenger response relative to matched sham. Compare sizes at matched local input. | Tests gain under a characterized input. Counts of spontaneous reporter excursions alone mix upstream event frequency, input strength and downstream establishment. |
| A smaller sham-relative local deficit after a small incumbent perturbation | A calibrated local reduction of active signaling, with a matched sham and whole-boundary imaging, followed through a fixed early window and longer recovery. An appropriate biochemical intervention must change activity; reporter photobleaching by itself does not implement the model's active-deactivation kick. | Tests recovery along one declared disturbance. Separate actual local loss/rebound from the shrinking paired deficit, and fixed-region loss from tracked-site redistribution. It does not measure generic stability. |
| Cue-on challenger activation can fail to persist after withdrawal | Image both sites during a fixed substantial local input and after removal; record local maxima/regions, integrated signals, time to loss and any later return. Retain censored event lifetimes. | Separates activation, sustained competition and retained replacement. One successful selected pulse does not define a universal threshold or establish that withdrawal has no effect. |
| Same-state ring competition changes with placement despite nearly equal initial gain | Vary stimulation location at fixed measured size and initial reporter organization; record both physical boundary routes, edge gaps and local input overlap. Repeat across sizes without treating their initial states as identical. | Tests placement as a modifier. Organization and local substrate/context remain possible mediators, so route length is not isolated automatically. |
| Strong opposing sites cancel the first moment | Calculate reporter-polarity vector magnitude/direction together with native boundary signal, total reporter activity, site counts and local protrusion. Withhold direction in weak-resultant intervals. | Distinguishes global asymmetry from local activity and identifies cancellation. The vector is not a unique reconstruction of the profile or a cell-motion law. |
| Noisy polarity has size-dependent angular and physical descriptors | Independent-cell or independent-recording activity-vector correlations, magnitude and local-site trajectories, reporting acquisition window and physical site displacement separately from angle. Compare migration tracks in parallel where feasible. | Tests activity-memory and its empirical relation to displacement/directness. The simulated input law is not established as cellular noise, and an activity lag integral is not the conditional OU velocity-memory estimate. |

For the planned WGD pilot, the agreed starting pair and collaboration are
established, while density, frame interval, field count and reporter feasibility
remain open. The outline should reserve substantial motility results and
describe these mechanism-linked reporter measurements as proposed follow-up
unless the final checkpoint or direct author evidence establishes a firmer
experimental commitment.

## Critical coherence review to apply after the first publication

1. **Lead with what the evidence explains.** The accepted five-model,
   three-challenge comparison establishes model/stimulus-dependent polarity
   replacement, including nonmonotone responses. The mechanism figures should
   explain that pattern or identify testable modifiers, not imply that all
   five models share one greater-stability law. Add explicit links between the
   original duration curves and the corresponding native long-cue histories.
2. **Preserve the deterministic/stochastic bridge.** Small-kick conditional
   recovery, finite challenger fate and noisy activity memory are related
   observables, but no established causal decomposition links their effects.
   Positive conditional recovery while the actual core declines is a useful
   measured distinction, not contradictory evidence to be hidden.
3. **Make model-specific results visible.** Goryachev emphasizes increasing
   peak/concentration and strong response changes; broadened fronts in
   Wave-pinning/Holmes can retain multiple sites; corrected deBelly and Spring
   can orient under the fixed long cue without demonstrating every field has
   settled. Show null/near-null gain and small recovery effects honestly.
4. **Protect geometry claims.** Ring comparison supports a placement effect
   within a declared geometry. It changes both routes and local context.
   Cross-size angular memory can partly reflect the circumference-to-angle
   mapping. A larger angular correlation alone cannot establish slower
   physical movement of an activity site.
5. **Keep the empirical inference conditional.** Observed displacement and
   directness support a persistence hypothesis; fitted OU memory is model-
   conditional. Neither is the same as a polarity-replacement cue duration
   or activity-vector lag integral. The outline should show these as linked
   hypotheses and measurements, not interchangeable persistence units.
6. **Control main-text density.** The ring component is already two continued
   review sheets. The first outline can include those delivered sheets; the
   second outline can retain readable main panels and move complete-grid
   diagnostics/methods to supporting pages while preserving every proposed
   main figure and detailed legend.
7. **Keep the WGD experiment substantial.** Give the planned SUM-159 WGD+/−
   contribution a conspicuous figure and results placeholder, with measured
   size, migration trajectories, displacement/directness and polarity–
   protrusion measurements labelled proposed rather than reported. Do not
   fabricate pilot results or settled acquisition settings.

The final post-delivery review must inspect actual assembled figures and
ensemble tables before concluding that these concerns are resolved. Existing
review files and green numerical audits qualify their stated scope; they do
not by themselves establish manuscript coherence.
