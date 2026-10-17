Manuscript outline — size effects across progressively richer 1D polarity and migration models

Working title

Cell enlargement selectively disrupts polarity resolution across progressively richer models of cell migration

Alternative: Cell size preserves polarity competence but destabilizes winner-take-all polarity selection

Central biological message

Cell enlargement does not generically abolish the ability of a cell-polarity system to form a polarized state.

The reproducible size-sensitive phenotype is instead the resolution of competing polarity states:

enlarged systems more readily retain competing fronts;

enlarged systems are less able to erase an established pole and cleanly repolarize after an opposing cue.

This qualitative size-response fingerprint is already visible in a minimal mass-conserved wave-pinning model and persists in a richer Cdc42–Rac–Rho network.

Adding cell mechanics and tension-dependent inhibition changes whether chemical polarity is translated into sustained movement and can strongly restrict the range of conditions that remain motile.

The biological implication is that large cells may remain capable of front formation while becoming less decisive in front selection, repolarization, and directional migration.

Results

1. A common phenotype panel enables direct comparison across model complexity

Introduce the study as a progression through three levels of biological richness:

a minimal one-GTPase wave-pinning system;

a three-species Cdc42–Rac–Rho polarity network;

the same polarity modules coupled to a deformable moving interval with tension feedback.

Define the same core perturbation phenotypes wherever possible:

single-cue polarization — can one persistent front form?

opposing-front competition — does one front win or do two fronts coexist?

polarity reversal — can an established front be replaced after cue reversal?

migration output — once mechanics are present, does polarity generate sustained directional motion?

Treat calibration maps as model qualification rather than as the central biological result.

Establish that both chemical model classes contain broad responsive regimes in which localized cues can produce persistent polarity.

Biological conclusion: the comparison is made between models that are demonstrably capable of cell-like polarization, so subsequent size effects report changes in polarity behavior rather than simple loss of model function.

Figure/Supplement

Main-text schematic: progression from 1-species → 3-species → polarity + mechanics, with the same phenotype icons repeated underneath.

Supplementary regime maps:

single-species polarization across mu, Lambda, and reservoir depth h;

three-species polarization across kappa_R, Lambda, and inhibitory Hill exponent n.

2. In the minimal wave-pinning model, enlargement preserves polarity formation but impairs polarity replacement

Apply the concentration-preserving isometric size transformation:

volume multipliers 0.5×, 1×, and 2×;

linear dimension scaled by V^(1/3);

local concentrations and reaction/diffusion kinetics preserved;

membrane/cytoplasmic reservoir depth scaled consistently in the single-species model.

Single-cue response

Matched 1×→2× endpoint state is unchanged in 847/867 conditions (97.7%).

Overall polarized endpoint prevalence changes only modestly, from 9.34% at 1× to 10.27% at 2×.

Interpret this as strong preservation of basic polarity competence under enlargement.

Opposing-front competition

Persistent two-domain outcomes change only weakly in the minimal system: 4.24% at 1× versus 4.39% at 2×.

Thus simple geometric enlargement is not sufficient to cause a large multipolarity phenotype in this minimal chemistry.

Polarity reversal

Among runs that are left-polarized before the opposing cue, clean reversal falls from 71.9% at 1× to 48.6% at 2×.

Multiple persistent polarized domains rise from 24.6% to 51.4% of eligible cases.

The dominant size phenotype is therefore a failure of winner-take-all reorientation rather than a failure to generate polarity.

Biological conclusion: a minimal conserved polarity circuit can tolerate increased cell size during initial front formation, but enlargement makes an established polarity state harder to replace cleanly.

Main Figure 1

Panel A: model and WGD-null scaling schematic.

Panel B: paired endpoint-state retention at 1× versus 2×.

Panel C: two-domain persistence across 0.5×, 1×, 2×.

Panel D: clean reversal versus multi-domain outcome across 0.5×, 1×, 2×.

Representative trajectories: one size-insensitive single-cue case and one reversal case that switches from clean reversal at 1× to competing poles at 2×.

3. A Cdc42–Rac–Rho network reproduces the same size-response fingerprint and strengthens the multipolarity phenotype

Repeat the same concentration-preserving size transformation in the six-field Cdc42–Rac–Rho model.

Single-cue response

Matched 1×→2× endpoint state is unchanged in 487/507 conditions (96.1%).

Polarized endpoint prevalence increases only modestly, from 41.6% at 1× to 44.4% at 2×.

Basic polarity competence again remains largely intact.

Opposing-front competition

Persistent two-pole outcomes increase monotonically with size:

5.0% at 0.5×

8.6% at 1×

13.6% at 2×

This is a substantially stronger size-dependent coexistence phenotype than in the one-species model.

Polarity reversal

Clean reversal among eligible runs decreases:

62.8% at 0.5×

52.7% at 1×

40.9% at 2×

Multiple-domain outcomes rise from 36.1% at 0.5× to 50.2% at 2×.

Network crosstalk therefore does not remove the minimal-model size phenotype; it makes persistent front competition more prominent.

Biological conclusion: increased biochemical realism preserves the same qualitative consequence of enlargement—reduced decisiveness of polarity selection—while revealing a stronger tendency for large cells to support simultaneous competing fronts.

Main Figure 2

Panel A: three-species network and Rac-based common polarity readout.

Panel B: matched endpoint-state retention.

Panel C: two-pole persistence versus volume.

Panel D: reversal outcome versus volume.

Panel E: representative Cdc42/Rac/Rho profiles showing single-front resolution at smaller size and persistent competing-front organization after enlargement.

4. Size acts on polarity selection more consistently than on polarity competence

Directly compare the 1×→2× effect size for the same three observables in the two chemical models.

Single-cue polarization

one-species: +0.9 percentage points

three-species: +2.8 percentage points

conclusion: no evidence for a broad loss of the ability to polarize with enlargement.

Persistent opposing poles

one-species: +0.15 percentage points

three-species: +5.0 percentage points

conclusion: network crosstalk amplifies the capacity of larger domains to sustain competing fronts.

Clean reversal

one-species: −23.3 percentage points

three-species: −11.8 percentage points

conclusion: impaired repolarization is the most consistent directional size effect across model classes.

Use the paired 1×→2× reversal transitions as an additional robustness measure:

of conditions that reverse at 1×, only 68.3% remain reversed at 2× in the one-species model;

63.0% remain reversed in the three-species model.

State the cross-model synthesis explicitly:

size has little effect on whether a polarity circuit can establish a front;

size has a much larger effect on whether competing spatial states are resolved into a single winner;

the richer network strengthens the multipolarity component but does not change the overall direction of the phenotype.

Biological conclusion: cell enlargement preferentially compromises polarity decision-making rather than polarity generation.

Main Figure 3

Cross-model size-response fingerprint:

x-axis: single-cue polarization, persistent opposing poles, clean reversal;

y-axis: percentage-point change from 1× to 2×;

one-species and three-species shown side-by-side.

This should be the manuscript's synthesis figure and a candidate graphical-abstract core.

5. Mechanical coupling converts polarity defects into changes in migration and exposes strong tension control

Extend both chemical modules to a moving one-dimensional cell with:

Kelvin–Voigt endpoint mechanics;

active endpoint forces generated from polarity signals;

soft and stiff mechanical cases;

reference lengths from 5 to 50 µm;

global or spatial tension-dependent enhancement of GTPase/Rac deactivation.

Retain the same single, dual, biased-dual, and reversal perturbation logic so the moving-cell analysis remains connected to the fixed-domain phenotype panel.

No tension feedback

single-cue runs in the one-species moving model remain terminally single-polarized across all 20 length × mechanics combinations.

mean total travel across the length/mechanics sweep is approximately 11.2 µm over 500 s.

in the three-species moving model, mean total travel is approximately 9.37 µm and terminal polarity is more restricted across the length/mechanics range.

Global tension feedback

increasing feedback strength strongly reduces travel:

one-species: 11.2 → 4.82 → 0.78 µm for alpha_t = 0, 5, 20;

three-species: 9.37 → 2.91 → 1.95 µm.

terminal single-cue polarity is correspondingly lost across much of the length/mechanics design.

Spatial tension feedback

the same nominal feedback strengths produce even stronger suppression:

one-species: 11.2 → 0.45 → 0.16 µm;

three-species: 9.37 → 1.06 → 0.35 µm.

local propagation of a tension-linked inhibitory signal therefore strongly constrains the persistence of front activity and sustained migration.

Analyze the Tier-2 per-length outputs as the decisive final comparison:

plot each common phenotype against reference length for 1-species, 3-species, and mechanochemical models;

test whether the fixed-domain enlargement signature—preserved single-cue competence, increased competing fronts, reduced reversal—survives after cell shape can change;

overlay displacement or directional-travel output to show when a chemical polarity defect becomes a migration defect.

Biological conclusion: mechanics does not merely report the polarity state; mechanochemical feedback can actively suppress it, providing a route by which size-dependent front competition can be converted into reduced or less decisive cell migration.

Main Figure 4

Panel A: moving-cell model with endpoint force, elastic tension, and global/spatial feedback.

Panel B: total travel versus tension-feedback strength for both chemical models and both feedback modes.

Panel C: fraction of single-cue runs retaining a terminal polarized state across the tested length/mechanics range.

Panel D: per-length migration output from the Tier-2 merged files, stratified by model richness.

Panel E: representative trajectory showing polarity, cell length, center position, and tension before and after feedback engagement.

6. A conserved mechanistic interpretation emerges across model richness

Summarize the hierarchy:

minimal polarity circuit: enlargement mainly impairs replacement of an existing front;

Rho-family network: the same defect persists and competing fronts become more size-sensitive;

mechanochemical model: tension feedback can suppress front persistence and translate polarity indecision into weak or curtailed migration.

Frame the central mechanistic idea as competition for spatial dominance:

larger cells provide more spatial separation between candidate polarity sites;

increased separation makes mutual exclusion/coarsening less decisive;

established fronts therefore coexist more readily with newly stimulated fronts or resist complete replacement.

Connect to cell biology:

large/polyploid cells need not be intrinsically incapable of polarizing;

their vulnerability may instead appear during turning, repolarization, navigation through changing cues, and suppression of supernumerary fronts;

these are precisely the behaviors required for persistent migration through heterogeneous tissue environments.

Final biological conclusion: the most robust predicted consequence of cell enlargement is a loss of polarity decisiveness—not a loss of polarity itself.

Methods

1. Modeling strategy and progressive model hierarchy

Define the three levels of model complexity and the rationale for comparing them with a common phenotype panel.

State that all models represent the cell along a one-dimensional longitudinal axis.

Separate:

fixed-domain chemical experiments;

concentration-preserving isometric volume scaling;

moving-domain mechanochemical experiments.

2. Single-species mass-conserved polarity model

Active membrane GTPase u and rapidly diffusing inactive pool v.

Reaction–diffusion equations and reaction term.

Conserved amount M = integral[u + h v] dx.

Dimensionless controls rho, Lambda, kappa, Gamma, mu, h, and Hill exponent n.

Homogeneous resting-state initialization.

Localized Gaussian activation cues.

3. Three-species Cdc42–Rac–Rho model

Six fields: active/inactive Cdc42, Rac, and Rho.

Rho inhibits Cdc42; Cdc42 promotes Rac; Rac inhibits Rho.

Conserved total for each GTPase.

Parameterization of totals, activation controls, diffusion, coupling, and Hill exponent.

Rac stimulation used for matched perturbation experiments.

4. Common perturbation protocols

Single left cue.

Simultaneous equal and weakly biased opposing cues.

Two seeded initial domains for fixed-domain competition.

Sequential left/right cues for reversal.

Relative-position cue coordinates so stimulus geometry scales with the cell.

5. Concentration-preserving cell-size scaling

Volume multipliers 0.5, 1, 2.

Linear scale s = V^(1/3).

Scale physical length and Lambda by s.

Preserve concentrations, diffusivities, reaction kinetics, cue amplitude, cue timing, and cue fractional geometry.

Single-species reservoir depth h -> s h, with abundance recalculated to preserve resting concentrations.

Three-species concentrations/totals represented in the implemented one-dimensional convention and preserved under the stated transformation.

Pair each scaled simulation to the same baseline molecular condition.

6. Moving-domain mechanics

Affine moving interval with mechanical degrees of freedom at the two endpoints.

Kelvin–Voigt internal mechanics.

Soft and stiff mechanical cases, with stiffness and viscosity normalized by reference length.

Active endpoint forces:

one-species from active GTPase;

three-species from outward Rac and inward Rho contributions.

Reference-length sweep: 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40, 50 µm.

7. Tension-dependent polarity feedback

Global feedback from positive elastic spring tension.

Spatial feedback via a diffusing/relaxing local tension-associated inhibitory field.

Feedback strengths alpha_t = 0, 5, 20.

Added deactivation transfers material to the inactive pool and preserves total GTPase abundance.

8. Numerical solution

Conservative finite-volume discretization on the moving interval.

Amount-based integration so changing control-volume size is handled conservatively.

No-flux boundaries.

BDF integration with stated relative/absolute tolerances.

Node-density rules for fixed and moving experiments.

Conservation-error thresholds and trajectory rejection rules.

9. Common phenotype definitions

Primary active field:

u in the one-species model;

Rac for matched endpoint analyses in the three-species model.

Normalized polarization amplitude.

Connected high-domain counting.

Terminal labels: homogeneous, single-polarized, dual-polarized, multi-peak, transient, oscillatory.

Reversal eligibility based on a prechallenge single left-biased domain.

Reversal outcome classes and reversal-time definition.

Keep three-species composite Cdc42/Rac/Rho coherence metrics as secondary mechanistic annotations.

10. Migration and mechanics observables

Cell center and center velocity.

Net displacement.

Absolute, signed, and polarity-projected speed.

Cumulative leftward/rightward travel and direction classification.

Final/minimum/maximum length and strain.

Elastic tension and feedback-exposure summaries.

Prespecified time windows for cue, establishment, reversal, reorientation, and late behavior.

11. Analysis strategy

Use matched condition-level comparisons for 0.5×/1×/2× size experiments.

Emphasize effect direction and paired state transitions rather than raw prevalence alone.

For cross-model synthesis, report the 1×→2× change for the same phenotype definitions.

For moving-length experiments, preserve reference length as the x-axis rather than treating the ten lengths only as replicates; stratify or facet by mechanics and tension mode.

Use representative trajectories only after the population-level pattern is established.

Suggested Supplementary Results/Figures

Calibration/regime maps for both chemistries.

Single-species local bistability map.

Full 0.5×/1×/2× endpoint transition matrices.

Full competition maps across parameter grids and initial asymmetry.

Reversal-time distributions among successful reversals.

Three-species spatial-stability and composite Cdc42/Rac/Rho annotations.

Soft-versus-stiff mechanochemical comparisons.

Global-versus-spatial tension-feedback diagnostics.

Conservation and numerical-quality summaries.