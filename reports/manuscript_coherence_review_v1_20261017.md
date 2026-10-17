# Independent manuscript coherence review of outline v1

17 October 2026. Reviewed artifact: [outline v1](../manuscript_outline/outline_v1.pdf),
31 pages, SHA-256 `59b70848912e674b546339fac7d3b11e37d4ce4bf5eb966508c451716facb2e1`,
published in commit `6199cdc33fea355133e4ede149131455dc5ccb91`.
Reviewer: independent manuscript-coherence agent, after the first publication.
This is an assessment of the supplied artifact, not author approval or an
endorsement inferred from previous agent reviews.

## Verdict and scientific argument

The outline preserves the intended manuscript components and is a defensible
first integrated checkpoint. Its central hypothesis is coherent: enlargement
changes polarity organization and replacement, and these changes may contribute
to migration persistence. It does not yet present a complete causal explanation
of migration. The reserved SUM-159 results section is prominent, has no invented
observations, and remains separate from mechanism-linked experimental proposals.

The second outline requires substantive repair. The noisy-ring main image hides
the physical-site comparison that most strongly qualifies its persistence
message. The mechanistic section presents useful new response measurements but
omits several existing explanations and evolved interventions. Merely renaming
those panels as descriptive phenotypes would leave the manuscript vision's
mechanistic component underdeveloped. Restore the existing explanations beside
specific measurements, distinguish their evidence types, and keep the common
active-loss assay as one contribution. The migration-B revision requested by the
author also remains outstanding in v1.

Existing completed evidence suffices for these repairs. No fresh native campaign
is necessary for a coherent, scientifically conditional second outline. A fresh
campaign becomes necessary only if the manuscript chooses stronger claims such
as intervention-specific threshold shifts, causal shares of reversal duration,
pure transport-distance effects or a calibrated polarity-to-velocity law.

## Evidence actually inspected

The review read [MANUSCRIPT_VISION](../MANUSCRIPT_VISION.md), root and
[manuscript agent instructions](../polarity1d/manuscript/AGENTS.md), the
[foundation audit](manuscript_foundation_audit_20261017.md), the delivered
[checkpoint](../current%20manscript%20figs/vision_aligned_checkpoint.md), the
outline's JSON and extracted full PDF text. Actual rendered PDF pages 3, 6, 9,
10, 11, 15, 18, 20 and 26 were inspected, covering the canonical comparison,
all mechanism/noise sheets, mathematical definitions, literature figure and
full duration supplement.

Scientific inspection included the linear/ring/noise actual-results reports;
the complete noise persistence table; full-window lag-8000 rows read directly
from the compressed physical-displacement table; the seed/direction table
definitions; the refined-OU actual results; the current mechanism PDF; corrected
deBelly endpoint states and decline-clock CSVs; and the actual resting-ablation
CSV/validation. The original three-test `cell_evidence.json` and producer were
inspected for unchanged Wave-pinning, Holmes and Spring interventions. Original
centered-deBelly controls were excluded.

The direct Oct-17 user messages at lines 124, 215, 307 and 337 of the session
identified in the foundation audit were independently read with role `user`.
They establish the theoretical-basis/persistence hypothesis, continued importance
of existing figures, integrated mechanisms/measurements, and the separate major
experimental gap. The direct Oct-13 messages at lines 148 and 214 were also
read: reserve lines for the model, show scraped observations as small points,
show curvature/exponent behavior, and avoid implying recovered experimental lag
counts. These are stronger evidence than agent-authored completion statements.

## Required repair 1: show physical-site motion and eligibility with angular memory

Figure 4, PDF page 15, shows initial-axis memory, magnitude, uncentered vector
correlation and selected local histories. It does not show physical displacement
or the eligible-pair fraction. The explanatory text on page 18 gives a generic
geometry caveat but omits the actual contrasting Holmes result.

Direct reading of
`vision_noise_20261017/supporting/data/vision_noise_site_displacement.csv.gz`
gives these high-input, full-window, lag-8000 values. RMS is the square root of
the stored mean squared shortest physical displacement; intervals use the same
monotone transform of the stored whole-seed bootstrap bounds.

| Model | C20 RMS [95% interval] | C40 RMS [95% interval] | C60 RMS [95% interval] | Eligible pair fractions C20/C40/C60 |
| --- | ---: | ---: | ---: | --- |
| Wave-pinning | 3.249 [2.714, 3.735] | 2.826 [2.399, 3.292] | 2.254 [2.032, 2.476] | 0.9923 / 0.9649 / 0.9461 |
| Holmes 3 | 4.568 [3.882, 5.174] | 6.334 [5.220, 7.279] | 5.539 [4.248, 6.695] | 0.9208 / 0.6640 / 0.3735 |
| Goryachev | One eligible pair; no population estimate | Unavailable | Unavailable | 0.00001736 / 0 / 0 |

The Goryachev C20 CSV contains a finite numerical RMS of 6.484, a collapsed
interval and bootstrap-defined fraction 0.621. That number must not be plotted
as a precise population curve: it is one pair out of 57,612, with zero eligible
pairs at C40/C60. Unavailability is not evidence of absent local activity.

Holmes' stronger angular memory with size does not correspond to a monotone
decrease in physical displacement. Its eligibility also falls sharply, so
eligible subsets differ across sizes. Equal physical displacement subtends a
smaller angle on a larger ring. Wave-pinning supports both increased angular
memory and lower conditional physical displacement in this ensemble; Holmes
supports increased angular memory with a different physical-site outcome.

Add a matched diagnostic sheet to main Figure 4 using the existing curves,
native physical units, shared scales, uncertainty and pair coverage. Retain
single-region/peak/contrast/magnitude eligibility definitions and state that
endpoint eligibility does not establish uninterrupted site identity. Do not
erase the Holmes exception by averaging it with Wave-pinning or relabel missing
Goryachev values as stable sites.

## Required repair 2: separate directional memory from forced-state differences

The sentence leading Figure 4's narrative, “larger rings retain more vector lag
memory in all three models,” is numerically true for the chosen restricted
descriptor, but scientifically overpromotes the Goryachev contrast. Its C20
high-mean deterministic control ends at relative resultant approximately
`3.24e-17`, whereas C40/C60 retain approximately 0.743/0.836. It is a changed
forced-state regime, not a clean comparison of noise-induced rotation around
the same polarized state. The current legend correctly discloses this; the
primary message should follow that distinction as well.

Lead with Wave-pinning and Holmes retaining more angular activity memory under
the declared forcing. Present Goryachev as the counterexample demonstrating that
positive input mean and state change can dominate a memory descriptor. Keep its
data, controls and local histories visible.

The estimator is uncentered and finite-window, with no mean subtraction. A
nonzero average axis contributes to its plateau. Sixteen of eighteen noisy
conditions do not reach `1/e` by lag 8000. Consequently restricted signed areas
are measured horizon-dependent descriptors, not resolved exponential persistence
times. Actual high-input first/full/last-half areas retain the size ordering
but are not identical; for example Holmes C20 is approximately 3976/4304/4481
and C60 7286/7103/6960. Show or summarize these window diagnostics rather than
assuming stationarity from a completed run. No stationary-correlation or
unlimited-time claim is needed to retain the useful result.

## Required repair 3: retain genuine mechanistic explanations beside measurements

The overview promises to explain the three-assay duration trends through gain,
recovery and competition. The new common assay removes a small incumbent
activity deficit and measures its shrinkage at a fixed model clock. It is not
an independently evolved challenger experiment, and its values do not allocate
reversal time. Wave-pinning/Holmes gain and near-doubling recovery are almost
unchanged despite strong patch-duration size effects. deBelly's Rac-only kick
initially leaves MCA/mechanics untouched. The new prolonged-cue sheet uses
fixed support, whereas the canonical comparison includes relative support and
constant slope with important opposite/nonmonotone responses.

There is nevertheless substantial existing explanatory evidence. The attached
[all-five mechanism crosswalk](../manuscript_outline/mechanism_crosswalk_v2.json)
provides compact table cells, detailed cue/region/clock/gate definitions, current
citations and the smallest additional source files needed for a portable sheet.

| Mechanistic account to retain | Actual supporting evidence | Corresponding measurement and scope |
| --- | --- | --- |
| Wave-pinning: wider conserved active domain and post-cue competition | L15 width-3 pulses of 20/60 collapse/retain the challenger; fourfold inactive diffusivity rescues L39 at the same native 5120-unit failed pulse, whereas remeshing alone does not | Record incumbent extent, challenger growth and both domains after withdrawal. This is an evolved fixed-dose rescue, not an intervention-specific threshold curve or constant-speed proof. |
| Goryachev: less challenger substrate and stronger response to specified incumbent damage | Native activation multiplies inactive concentration; current L19→39 gain falls to 0.654 and active-loss repair rises 3.439-fold. Separate inactive-bolus removal/return gives immediate reaction/deficit 0.9731→15.2854 at L9→39 | Measure substrate, early activation and localized sham-relative loss/rebound separately. The same conserved pool couples these effects; no additive causal percentages. |
| Holmes: coupled Rac/Rho spatial competition and support dependence | Fourfold inactive mobility does not rescue fixed-support L28/L50 at the cap; relative support passes orientation at all eligible sizes. Same-C40 ring placement produces merger versus distinct sites | Track both activity fields and site fates while changing support/placement. Relative support changes dose/coverage; transport is not a universal rescue. |
| Corrected deBelly: inherited old-end MCA depletion, then longest-cell refill delay | Corrected relative-20% endpoint CSVs: L23→31 basal old MCA 0.05348→0.04823 and 180-s levels 0.06451→0.06044; 10% old-Rac decline 123.49→142.97 s. L39 starts less depleted than L37 but refills much more slowly, with decline at 637.39 s | Measure old/new-end attachment, activity, flow and protrusion together. Local coupled budgets provide an explanation, not a universal threshold. A separate corrected tension-to-Rho resting ablation reduces displacement spread 50.4%, establishing a starting-state contribution rather than a reversal-time contribution. |
| Spring: extension-dependent deactivation and increasing ramp input | Setting only the extension-dependent deactivation coefficient to zero during cue defeats 4/4 native successful gradient pulses; matching ramp peak defeats 3/4. The longer two patch doses still succeed under that term's removal | Measure extension/activity through characterized inputs. These evolved, geometry-dependent fixed-dose controls retain motion/dilution and restore native recovery; they do not remove all mechanics. |

A readable companion sheet under Figure 2 can supply this crosswalk without
replacing the delivered three checkpoint candidates. Its explanatory text should
state where an evolved intervention exists and where evidence instead consists
of exact reaction laws, coupled source budgets or temporal association. Keep the
mechanism/readout pairing beside the result, as requested by the author.

For corrected deBelly, use only independent-endpoint evidence. The selected
relative-20% pulses at L23/L31/L37/L39 are
232.2449/285.2806/355.9949/809.0816 seconds; native-field flags are
false/true/true/false. The later fixed-support-1000 campaign's field qualification
must not be transferred to these selected pulses. The 0.068 MCA midpoint is a
smooth feedback parameter, not a withdrawal threshold. The full response curve
decreases over L18→21 and L31→33 even while the descriptive old-Rac-decline clock
increases. These limits prevent a misleading one-clock explanation.

Assay lineage also matters for Goryachev. The primary 0.5% whole-inactive removal
held for 0.1 units gives the 15.708-fold immediate reaction/deficit increase, but
its 90%-recovery clock is 18.87 at L9 and 2.64 at L39. A slower late-tail result
in the canonical older report concerns a different activation-suppression assay.
Do not combine the 15.708-fold number and a generic slower-tail sentence as
though they came from the same experiment.

## Required repair 4: implement the author's migration-B encoding request

V1's integrated migration image on page 20 retains connected observed curves
and a different line style for fitted OU curves. This precedes the direct
Oct-13 request. Make an outline-owned editorial reassembly, preserving the
selected fit and all scientific data. Show sparse, explicitly chosen display
samples of the digitized population curves as points; reserve continuous lines
for OU predictions. Add matched-window local log-slopes or another truthful
curvature display. Their purpose is to reveal the discrepancy, not hide it.

The chosen display samples are not recovered author-measured lag locations,
independent cells or a sample size. The observed descriptive crossover near
2.32 hours must remain distinct from the fitted crossover near 1.11 hours.
Fitted curves must remain unchanged at every original fitting lag; slopes must
use identical windows for digitized and model curves. The foundation audit
already supplies an actionable bounded renderer plan.

The rest of the literature interpretation is appropriately conditional and
should be retained: Xuan's lower sampled speed with greater late MSD, Mallin's
shorter path with greater directness, heterogeneous polyploid/PACC comparisons,
and the OU hierarchy's prior dependence/model discrepancy. Five conditional
positive persistence estimates do not constitute five independent direct
measurements of the same WGD effect. Speed-only Jemaà comparisons do not identify
their persistence times without hierarchy/priors. The new experimental gap must
not be filled with those fit estimates.

## Readability, definitions and source lineage

The equations on inspected pages 9 and 18 render as mathematical expressions,
with fractions, subscripts, centered displays and right-aligned equation
numbers. The vector definition retains magnitude; the uncentered correlation
states its realization weighting and horizon. Preserve that mathematical
precision. Define the measured physical-site statistic separately when adding
the diagnostic sheet, and retain endpoint gating/coverage.

The main sheets are dense but reviewable. Figure 2's lower kymograph annotations
and Figure 4's condition keys are small in portrait packaging. The two-sheet
ring presentation is justified; compressing it into one page would lose useful
information. Prefer additional readable companion sheets to tiny new insets.
Inspect any regenerated migration B at its final PDF size. Internal source
working numbers differ from proposed manuscript numbers; v1 labels this
honestly, but a new outline-owned render can remove that navigational friction.

Keep the canonical five-model Spearman comparison and full duration supplement.
They retain corrected deBelly, successful-range/censoring distinctions, negative
Spring gradient association, nearly null Holmes gradient orientation and
nonmonotone curves. Do not filter these in the name of a simpler story.
Older seven-model/four-assay campaigns and centered-deBelly controls are not
substitutes for the retained foundation.

For portability, a stored PNG and an in-place rebuild are different promises.
The revised renderer should read committed compact evidence and write only
outline-owned assets; independently verify it in a clean checkout without bulk
native histories. Preserve source hashes and numerical identities. The crosswalk
JSON names eight small additional source files; it does not require committing
native simulation archives. Existing manifests/audits support source lineage but
do not prove a stronger biological or causal claim.

## Prioritized action and fresh-campaign decision

1. Add physical displacement/eligibility beside noisy angular memory, retaining
   the Holmes exception and Goryachev unavailable statuses. Change the lead noise
   message to distinguish WP/Holmes angular memory from Goryachev forced-state
   loss. Existing tables suffice.
2. Add the all-five mechanism/evidence/measurement companion, using exact
   corrected sources and fixed-dose versus budget/association distinctions.
   Preserve mechanistic explanations rather than solving the mismatch through
   weaker labels alone. Existing tables and interventions suffice.
3. Revise migration B using truthfully labelled raster samples and curvature
   diagnostics; verify invariant selected-fit predictions and unchanged other
   scientific inputs. No new OU fit is indicated.
4. Build and inspect the second outline; verify clean-checkout regeneration and
   publication lineage. Preserve the first checkpoint and keep the experimental
   results placeholder prominent.

New simulations are **not required for these bounded claims**. If stronger
mechanistic claims become essential, specify the needed experiment first:
intervention-specific duration searches for threshold shifts; a corrected
deBelly counterfactual MCA-remodelling perturbation from exact corrected basals
to test a causal reversal-time contribution; or an explicitly matched transport
geometry/material design to isolate separation. None can be inferred from a
fixed-dose intervention or frozen-state source allocation. Likewise a direct
biological link between activity memory and migration requires the proposed
measurements, which remain pending. The outline can be complete as an outline
while those empirical outcomes remain explicit placeholders.

This review and its companion JSON create no large output groups and modify no
canonical figures, native results or Git index.
