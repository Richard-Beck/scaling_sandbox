# Length-mechanism figure drafts

**Current discussion:** the [two main-text candidates](main_text_figure_candidates_20261017.md)
present the revised substrate-capture figure and the new separation figure,
with fresh critiques. The drafts below retain the preceding exploration;
in particular, Draft 2's broad restoring-response interpretation is superseded
by the returned-state chemical analysis in current candidate 1.

These are discussion drafts for a single main-text figure. The biological
motivation is the evidence for greater migration persistence in WGD-associated
enlarged cells, summarized in the [current migration figure](main_03_migration_OU.png).
The current candidate structure has two parts: (1) enlargement changes the scale
and organization of the incumbent, potentially affecting its restoring response
and competitor activation; (2) greater separation between incumbent and challenger
may affect replacement independently of incumbent scale. These are candidate
mechanisms, with different levels of support below. The final message and panel
selection remain open.

The [supplied mechanism report](../polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf)
is the scientific starting point. The first three result drafts reuse retained
simulations. The later protocol schematic and polarization-feasibility pilot
develop the proposed distance experiment. None directly measures migration
persistence.

## Draft 1: challenger activation and incumbent persistence

**Proposed message:** Increasing length can change the relationship between
local responsiveness and the persistence of established polarity. Goryachev
shows progressively weaker early challenger buildup, whereas corrected deBelly
retains rapid, substantial challenged-end Rac activation alongside increasingly
persistent incumbent activity.

![Draft time-course comparison: challenged-end activity and incumbent retention across lengths in Goryachev and corrected deBelly.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/activation_incumbent_contrast_20261017/activation_incumbent_contrast.png)

[Vector PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/activation_incumbent_contrast_20261017/supporting/vectors/activation_incumbent_contrast.pdf)
· [Technical caption](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/activation_incumbent_contrast_20261017/supporting/metadata/caption.txt)
· [Measurements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/activation_incumbent_contrast_20261017/supporting/data/source_series.csv)

**Why this design serves the message.** The two model columns let the reader
follow challenged-end activity and incumbent retention separately. Length colors
connect the two signals within each model. Showing several lengths is essential:
the intended result concerns how enlargement changes the response sequence,
rather than simply whether two models have different intrinsic kinetics.
The displayed interval is the first 180 native time units within each model
(seconds for deBelly), entirely during cue application. Different selected pulse
durations therefore cannot create the displayed contrast. Raw challenged-end
concentrations expose response amplitude;
basal-normalized incumbent activity exposes retention of the starting signal.

**What supports it.** Goryachev uses the full spatial histories underlying report
Figure 7, re-expressed as endpoint trajectories. Corrected deBelly uses the native
Rac histories underlying Figures 8 and 10. The report independently shows that
the deBelly challenged-end Rac-dependent protrusion switch reaches at least 99%
activation within 3–5 s across all twelve eligible lengths. Its incumbent Rac
takes approximately 123, 143 and 637 s to decline by 10% at 23, 31 and 39 µm.
Those are descriptive response clocks, not reversal thresholds.

**What this draft does not establish.** Goryachev has an immediate small chemical
response; the directly supported claim is weaker early buildup, with no defined
front-establishment threshold.
Rac activation in deBelly is not proof of a new front that will survive cue
withdrawal. Native time units, signal magnitudes and cue amplitudes are distinct
between models, so absolute cross-model speeds or doses are not being calibrated.
The examples represent response patterns, not an exclusive division into
chemical and mechanical model classes; Wave-pinning also shows prompt,
substantial challenged-end activation alongside an active incumbent. These
trajectories do not establish a causal
partition of reversal time or spontaneous directional persistence.

## Draft 2: two chemical effects that can help retain polarity

**Proposed message:** Enlargement can reorganize a polarity system so that an
equal local challenge produces a weaker initial local activation rate while the incumbent
produces a stronger immediate restoring reaction. The effect is pronounced in
the selected Goryachev regime; Wave-pinning supplies a comparison with different
length dependence.

![Draft chemistry figure: Goryachev and Wave-pinning basal profiles, initial local cue activation, and immediate incumbent restoring gain across length.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemistry_main_text_candidate_20261017/chemistry_main_text_candidate.png)

[Vector PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemistry_main_text_candidate_20261017/supporting/vectors/chemistry_main_text_candidate.pdf)
· [Technical caption](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemistry_main_text_candidate_20261017/supporting/metadata/caption.txt)
· [Measurements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/chemistry_main_text_candidate_20261017/supporting/data/length_response_measurements.csv)

**Why this design serves the message.** The basal profiles give a physical
explanation for the lower-row measurements: Goryachev builds a taller incumbent
and lowers inactive concentration, whereas Wave-pinning primarily widens its active
region. Panel C then measures initial activation per unit local cue, and panel D
shows the response to an actual perturbation. The comparison tests whether these
effects accompany enlargement in both selected models. It is stronger support
for the proposed chemical mechanism than a reversal-duration curve alone.

**What supports it.** Report Figures 3–4 contain the underlying profiles,
reaction identity and perturbation results. Initial added activation equals local
cue strength multiplied by inactive concentration. The restoring measurement
follows removal of 0.5% of the whole inactive profile for 0.1 native time units,
then exact pointwise return, with a matched sham. It is excess basal reaction
within the fixed original active core divided by the active material lost there.
At length 39 relative to length 9, Goryachev has approximately 0.333 times the
initial local activation and 15.7 times the immediate restoring gain;
Wave-pinning has approximately 1.00 and 1.28 times, respectively.

**What this draft does not establish.** These coupled local diagnostics do not
measure the frequency of spontaneous front initiation or the probability of
takeover after a front forms. Per-unit-local-cue activation is distinct from
activation integrated over a length-dependent cue support. Immediate restoring
gain is distinct from late recovery time and global stability. The length range
is illustrative model evidence, not a calibrated size change caused by one WGD.

## Draft 3: how the incumbent itself grows

**Proposed message:** Enlargement changes the scale of an existing active site,
even before any additional site is challenged. Across all five selected models,
the incumbent contains more active material above background at every adjacent
eligible sampled length. The added material appears mainly as greater width in
Wave-pinning, Holmes, Spring and deBelly, and mainly as greater concentration in
Goryachev.

![Draft incumbent architecture comparison: basal Wave-pinning and Goryachev profiles, then width and peak concentration relative to their own length-9 states.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/incumbent_architecture_20261017/incumbent_architecture.png)

[Vector PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/incumbent_architecture_20261017/supporting/vectors/incumbent_architecture.pdf)
· [Technical caption](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/incumbent_architecture_20261017/supporting/metadata/caption.txt)
· [Measurements](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/incumbent_architecture_20261017/supporting/data/incumbent_architecture_measurements.csv)

The spatial profiles establish what width and concentration mean. Separate
normalized summaries then show which property changes with length, without
equating the models' native concentrations. All eligible sampled lengths are
included, and the constant/proportional guides are descriptive rather than
fitted laws. Goryachev uses an explicitly labelled incumbent-end view to resolve
its narrow site. [Report Figure 2](../polarity1d/manuscript/figures/length_mechanisms_20261008/revision_20261014/visual_palette/assets/initial_polarity_length.png)
provides the full five-model evidence.

| Model and illustrative lengths | Incumbent width | Peak active concentration |
| --- | --- | --- |
| Wave-pinning, 9 to 39 | 2.82 to 11.92 | 1.101 to 1.125 |
| Goryachev, 9 to 39 | 0.712 to 0.799 | 27.8 to 117.4 |
| deBelly, 23 to 39 | 7.83 to 12.22 | 0.889 to 0.902 |

Values are model-native. Width is the boundary-connected region above half its
own basal contrast; active material is integrated over that region. The uniform
five-model statement uses material above the field minimum. Including background
gives a small deBelly decrease from 37 to 39, despite its continuing increase in
width and above-background material. Absolute width should also be distinguished
from the fraction of the cell occupied by the site: a broader site need not
occupy a larger relative fraction.

**What the size comparison holds fixed.** Wave-pinning, Goryachev and Holmes
start with fixed mean conserved-pool concentrations, so their total available
material grows in proportion to domain length. Spring likewise starts with a
fixed mean active-plus-inactive concentration, with conserved amount proportional
to its reference length; subsequent extension changes its mean concentration.
Corrected deBelly maintains a fixed mean MCA density by its native uniform
correction, but its Rac and Rho obey source–decay reactions rather than conserved
pool constraints. Greater incumbent amount therefore needs to be interpreted
alongside these assumptions. Subtracting the inactive-end background does not
remove the effect of having a larger total pool.

The more informative contrast is spatial organization. Between lengths 9 and 39,
Wave-pinning's active width changes from 31.3% to 30.6% of cell length, whereas
Goryachev's changes from 7.91% to 2.05%. The latter maintains a narrow absolute
site while its peak concentration rises approximately in proportion to length.
These are approximate descriptions of the displayed range, not exact scaling
laws: eligible shorter lengths are included in the draft comparison, and
Goryachev width has a shallow minimum before increasing again.
The [basal measurements](../polarity1d/manuscript/data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv)
provide the raw widths, physical lengths and concentrations.

**Why this could be a main-figure message.** Spatial profiles and width,
concentration and amount measurements directly establish how extra cell size
changes an incumbent. They provide a common starting point across models before
introducing model-specific competition mechanisms. This would tell a reader
whether an enlarged cell maintains a broader active region or concentrates more
activity into a localized site, instead of attributing everything to the creation
of additional sites. Unchanged local kinetic parameters can accompany a changed
polarized starting state.

**What remains an interpretation.** An active signalling region can be treated
as a candidate protrusive region, but its 1D width and active integral are not
measured protrusion geometry or force. In deBelly, the native protrusion law adds
a specific functional link: between 23 and 39 µm, incumbent Rac changes only
slightly, from 0.889 to 0.902, across the native switch midpoint of 0.899.
The protrusion switch consequently rises from 11.6% to 65.2% activation, and its
outward cortical endpoint command rises from 0.0151 to 0.0939 µm/s (6.20-fold),
also affected by reduced tension inhibition.
The [native endpoint measurements](../polarity1d/manuscript/figures/length_mechanisms_20261008/revision_20261014/visual_palette/debelly/Figure_11_steady_state_mechanism_endpoint_metrics.csv)
provide these quantities. That command is distinct from membrane motion.
Larger active amount alone does not establish greater stability
or explain a reversal duration; the perturbation and transport evidence would
still be needed for those claims.

## Separating incumbent scale from challenger distance

The user's proposed second message concerns distance at fixed cell size and
incumbent state. The existing report supplies motivation and partial tests:

| Candidate claim | Current evidence | What remains to establish |
| --- | --- | --- |
| Enlargement reorganizes the incumbent | Report Figure 2; Draft 3; deBelly Figure 11 | Which measure best connects the architecture to a functional consequence |
| Incumbent maintenance becomes stronger | Immediate restoring reaction per unit damage in Figure 4 and Draft 2, particularly Goryachev | General stability and a causal link to the changed site architecture |
| Competitor activation is suppressed | Lower initial per-unit-local-cue activation in Goryachev, Figure 3 and Draft 2 | Generality: Wave-pinning has little change, while deBelly retains prompt substantial activation |
| Separation slows remote disruption | Figure 5 relocates the same inactive sink in Wave-pinning and Goryachev | A competing active site rather than irreversible removal, with boundary proximity controlled |
| Separation slows deBelly replacement | Spatial remodelling benchmark and MCA dynamics in Figures 1 and 9–12 | A controlled separation intervention; current lengths also change incumbent state and mechanics |

![Proposed distance protocol: move the same challenger cue around an intact chemical ring, or use the proposed deBelly two-arc construction around two unmodelled protrusive regions.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/distance_ring_protocol.png)

[Protocol PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/distance_ring_protocol.pdf)
· [Technical caption](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/caption.txt)

**Why this schematic serves the proposed message.** It makes the controlled
comparison visible: the incumbent and circumference stay fixed while the
challenger moves. It also shows both communication paths and distinguishes the
chemical ring from the deBelly construction that still needs arc coupling.
It illustrates the protocol, without drawing unmeasured response curves.

The proposed new experiment places the models on a circular connected domain
and varies challenger location at fixed circumference. Wave-pinning, Goryachev
and Holmes can retain their local reactions and use periodic diffusion. Each
requires a fresh verified cue-free ring checkpoint, reused exactly across cue
locations. Both arc lengths matter, and nearby challenges must be checked for
overlap or incumbent drift rather than distinct competitor formation.
The [three-model development adapter](../polarity1d/manuscript/code/mechanisms/ring_distance_prototype_20261017.py)
passes periodic flux, conserved-pool, rotation-covariance and spatial-readout
checks. The first preparation pilot now supplies a settled single-site checkpoint
for each of these three models; its new simulation evidence is shown below.

**What the distance panel would need to show.** After preparation, compare
challenger-site activity, incumbent-site retention and the complete spatial
profiles across separations. A distance dependence of early challenger activity
would implicate initiation; similar challenger activity with different incumbent
weakening would place the effect later in competition. Following a single moving
peak is different from forming a distinct competitor while the incumbent
persists. Finite-cue survival would be a separate readout after the maintained-cue
diagnostic. A weak distance effect would also be informative: it would limit the
distance explanation within the tested geometry and range, rather than establish
that the original length response is a transport delay.

For deBelly the user's proposed construction retains the unmodelled protrusive
regions and simulates the two intervening arcs. Their combined length stays
fixed as the challenger moves. Inspection of the authors' native 2D source
(commit `cd51f0f91e039929abac2b6ce466ff2edbb11b67`) confirms that it represents
the protrusive region by a cut boundary and a Rac/tension-dependent cortical
command. It also provides curved vector mechanics and projected membrane
tension. The released example is a truncated axisymmetric sphere with one
updated leading edge and an inactive opposite pole; it does not already
implement two movable protrusive gaps.
[Native driver and boundary conditions](https://github.com/VirtualEmbryo/mechanochemical_polarization/blob/cd51f0f91e039929abac2b6ce466ff2edbb11b67/mechanochemical_2D_axisymmetric/src/MechanochemicalAxisymmetricVector.jl)
· [Native surface weak forms](https://github.com/VirtualEmbryo/mechanochemical_polarization/blob/cd51f0f91e039929abac2b6ce466ff2edbb11b67/mechanochemical_2D_axisymmetric/src/WeakForms.jl).

Those equations and cap conditions are the reference for the proposed arc
construction. A planar circle needs its own geometric measure and curvature
terms, rather than retaining the radial area weight and azimuthal terms of
axisymmetry. The essential open choices are how the two arc ends belonging to
one protrusive region share activity and mechanical boundary values, how chemical
fields communicate across that unmodelled region, and how one MCA budget and
uniform correction apply across both arcs. The 2D leading-edge switch also uses
the opposite-pole Rac as its threshold reference; copying this rule separately
at two active sites would change the intended two-site activation behaviour.
The symmetric two-site reference must therefore be stated explicitly.

Preserving the incumbent before challenge is a separate control. An inactive
cut carries mechanical and transport conditions; placing it at different
locations can already alter the basal state. A transparent pre-cue interface
would preserve one shared state, but opening that interface subsequently changes
the model conditions as well as applying the Rac cue. The first implementation
retains independent local flank commands and a single MCA correction over both
arcs, with no diffusion across either cut. A matched cut-only control will
measure the effect of opening a new boundary. The asymmetric two-arc
implementation now passes its state-transfer and native-operator checks.
Its [first controlled result](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_two_arc_distance_result.md)
shows that unequal-cut geometry reorganizes the arcs even without an added
Rac cue. It remains construction evidence rather than the main distance
mechanism. Spring now has an explicit mean-force circular
analogue, summarized below.

**User direction on gap timing:** either cue-triggered opening or a common
Rac-triggered rule is acceptable. A Rac-triggered rule must also govern incumbent
preparation: begin with an intact ring, apply the polarizing cue, and let that
same rule form the incumbent gap. The first deBelly prototype will use the
cue-triggered version consistently for preparation and subsequent challenge.
The Rac-triggered alternative remains open. The user identifies polarization of
all five periodic models, with parameters similar to the interval case, as the
first milestone. That milestone now passes for all five declared constructions;
the bounded challenger-location pilots are the next step.

## First milestone: can the circular models retain one polarized site?

![Five circular constructions retain a polarized site after cue withdrawal: Wave-pinning, Goryachev, Holmes 3, the Spring mean-force analogue and the deBelly point-cut loop.](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/circular_polarization_overview_20261017/circular_polarization_overview.png)

[Five-model overview PDF](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/circular_polarization_overview_20261017/supporting/vectors/circular_polarization_overview.pdf)
· [Result and construction details](circular_polarization_result_20261017.md)
· [Fresh critique](reviews/circular_polarization_r1.md)
· [Frozen run inputs](../polarity1d/manuscript/data/runs/ring_preparation_pilot_20261017/manifest.json)

**What this establishes.** With unchanged native local reaction parameters,
Wave-pinning, Goryachev and Holmes each form one settled active site on a
circumference-40 periodic domain. The first 10-unit preparation pulse suffices
for all three. Each then passes two consecutive 1000-unit steady windows and a
10000-unit cue-free holdout, with conserved pools and no detected site movement
during holdout. Profiles are rotated for display only; the saved checkpoints
retain their original coordinates. This is new preparation evidence, not a
result from the supplied PDF or a challenger-distance measurement.

The preparation uses a raised-cosine cue of total support 3, with peaks 0.1,
0.1 and 1 for Wave-pinning, Goryachev and Holmes. These choices retain the native
chemical parameters but are not identical to the current interval relative-20%
protocol, including its Goryachev peak 6.4. All three also qualify on an 800-cell
grid with the same parameters and preparation. In the
[two-grid comparison](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/ring_preparation_mesh_20261017/ring_preparation_mesh_comparison.png),
the largest peak and half-contrast-width changes are Goryachev's −0.275% and
+0.311%; the corresponding Wave-pinning and Holmes changes are below 0.001%.
These are two-grid sensitivity measurements, not a proof of continuum convergence.

**Spring also polarizes under a declared circular analogue.** Its native reaction,
diffusion, local force, extension-feedback and dilution coefficients are retained.
The new mechanical rule replaces endpoint-force sampling with the mean native
force around the circle, preserving the complete uniform-field native derivative.
At resting circumference 2 with 200 cells, the first 10-unit cue (peak 4, support
0.2) produces one settled site. Full fields and circumference pass the same basal
windows and 10000-unit holdout. Final circumference is 2.04418. This is evidence
for the specified circular global-extension analogue, with no claim of equivalence
to the authors' full 2D elasticity.
The [Spring derivation and embedded result](../polarity1d/manuscript/SPRING_RING_CLOSURE_20261017.md)
give the exact new rule, source mapping and numerical checks.

deBelly's first test uses a cue-created point cut retained after withdrawal and
tangential mechanics on a fixed planar reference loop. Reflection symmetry of a
single front permits a half-loop computation with one shared full-loop MCA
budget. Its [reviewed preparation protocol](../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/distance_ring_protocol_20261017/debelly_single_incumbent_prototype.md)
states the assumptions and current manuscript parameters. At circumference 62
with 402 full-loop-equivalent cells, the first 60-s cue (peak 3.6, support 12.4)
produces one retained site. The original 20000-s qualification remained
unresolved because an early tension window narrowly missed its criterion.
A separate exact-state continuation passes the 18000–20000-s prewindows and
the prospective 20000–30000-s holdout without changing parameters or tolerances.
This qualifies the preparation within the reflection-symmetric subspace.
The subsequent asymmetric two-arc pilot is linked above; complete controller
transfer, matched cut-only controls and one combined MCA correction were
retained.

The independent review accepts this first milestone, while preserving the
distinction between three native periodic reaction–diffusion systems and two
declared mechanical adaptations. It also identifies complete deBelly boundary
controller transfer and pre-cut state identity as requirements for the next
experiment. Five successful preparations are useful development evidence;
they do not by themselves supply the main-text mechanism.

## What belongs in the eventual main figure

The current structure could connect incumbent architecture (Draft 3) to one or
both measured response changes (Draft 2), while retaining the distinction between
challenger activation and incumbent retention (Draft 1). The ring experiment
would test the separate distance message rather than infer it from varying
cell length. A deBelly mechanism would draw on selected MCA evidence from
report Figures 9–12 and the inspected native 2D construction.

The final figure should select a small number of these components. The proposed
causal chain from incumbent scale to resistance remains a hypothesis to assess.
The circular-domain pilot now supplies a concrete alternative: Holmes shows
similar early challenger activation followed by merger at separation 12 or
two retained sites at separation 16, at fixed circumference 40. See
[candidate 2](main_candidate_2_20261017.md) and its
[fresh critique](reviews/chemical_distance_r1.md). The
[Spring analogue result](../polarity1d/manuscript/SPRING_RING_DISTANCE_PILOT_20261017.md)
distinguishes near-site bridging from separated competitor growth followed by
old-site extinction. These are new trajectories, rather than findings already
contained in the original PDF.
Stronger claims about WGD migration persistence would additionally require a
link between these polarity responses and spontaneous turning or steering in
moving cells.
