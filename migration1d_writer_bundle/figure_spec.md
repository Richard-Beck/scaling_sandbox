# Figure specification

## Overall figure strategy

The main figures should follow the progression:

**minimal 1-species polarity → 3-species Rho-family polarity → polarity coupled to mechanics/tension → cross-model synthesis**

Across these increasingly rich models, retain a common set of biological observables wherever possible:

* ability to establish a single front after a localized cue;
* ability to resolve simultaneous competing fronts;
* ability to replace an established front after cue reversal;
* persistence of multiple polarity domains;
* after mechanics is introduced: cell displacement, directional motion, and the relationship between polarity state and movement.

Each main figure should contain several complementary kinds of evidence:

1. **Model/experimental schematic** — what is being perturbed and what is measured.
2. **Representative spatial or spatiotemporal behavior** — profiles and preferably kymographs, not just endpoint cartoons.
3. **Population-level behavior across parameter space** — prevalence or effect-size summaries.
4. **Matched-condition analysis** — show what happens to the *same model condition* when size is changed.
5. **A visually explicit biological takeaway** — preferably apparent from the data rather than a large explanatory cartoon.

The reader should be able to move from left to right through each figure and understand:

**model → perturbation → individual behavior → population behavior → biological conclusion.**

Use the same phenotype terminology and iconography throughout. For example:

* one front = single polarized;
* two spatially separated fronts = competing/dual polarity;
* reversal arrow = successful reorientation;
* failed reversal = retained old front or coexistence of old and new fronts.

Avoid making simple bar plots the dominant visualization. Whenever possible show the underlying matched conditions, distributions, state transitions, or parameter-space structure.

---

# Figure 1. Cell enlargement preserves polarization but impairs polarity replacement in a minimal wave-pinning system

### Central result

The minimal polarity model establishes that enlargement does **not** generally destroy the ability to polarize. Its major effect is instead on the resolution of competing polarity states, most clearly during reversal.

### Panel A — Minimal polarity model and size-scaling experiment

Compact mechanistic schematic showing:

* active membrane species `u`;
* inactive rapidly diffusing reservoir `v`;
* activation/inactivation exchange;
* mass conservation;
* localized external activation cue.

Immediately adjacent, show the size-scaling transformation:

**0.5× volume → 1× → 2× volume**

with corresponding longitudinal sizes proportional to `V^(1/3)`.

Explicitly indicate that concentrations and kinetic parameters remain fixed while geometry and the reservoir term are scaled according to the defined WGD/null transformation.

Purpose: make clear that the experiment asks what geometry alone does to an otherwise equivalent polarity circuit.

### Panel B — Single-cue polarization is largely size invariant

Use **matched representative kymographs** from one or several parameter conditions.

Columns:

* 0.5×;
* 1×;
* 2×.

Rows can represent two representative parameter regimes.

Plot active `u(x,t)`.

Show:

* cue position and cue interval;
* cell boundaries;
* eventual polarized front.

The important visual point should be that the same front forms at all three sizes despite the increasing physical distance across the domain.

Beside the kymographs, show a compact matched-state summary for the entire parameter sweep.

Preferred visualization:

* 1× endpoint state on the y-axis;
* 2× endpoint state on the x-axis;
* transition matrix with counts.

Annotate prominently:

**847/867 matched conditions retain the same endpoint state (97.7%).**

This is more informative than simply plotting polarized fractions.

### Panel C — Size has little effect on simultaneous two-front competition in the minimal chemistry

Show endpoint profiles or kymographs for simultaneous left/right stimulation at 0.5×, 1×, and 2×.

Then show the population-level result as either:

* paired condition traces connecting 0.5× → 1× → 2× for domain number; or
* three aligned parameter-space phenotype maps.

Preferred parameter maps use the same axes and ordering at every size, with each model condition colored by:

* homogeneous;
* one dominant domain;
* two persistent domains;
* > 2 domains/other.

The reader should be able to see that the maps are nearly unchanged with size.

Include the aggregate two-domain incidence as a small annotation:

* 1×: ~4.24%;
* 2×: ~4.39%.

This panel establishes that enlargement alone does not inevitably generate multipolarity in the minimal model.

### Panel D — Enlargement strongly interferes with reversal

This should be the visual centerpiece of Figure 1.

Use one **matched parameter condition** displayed at different sizes.

For each size:

1. initial left cue;
2. established left front;
3. opposing right cue;
4. subsequent dynamics.

Show active-species kymographs rather than only final profiles.

Ideal example:

* 1×: old front disappears and the right front wins;
* 2×: new right front develops without complete elimination of the old left front.

Cue timings should be indicated directly above the kymograph.

### Panel E — Matched reversal outcomes across the whole parameter ensemble

Show outcome transitions with size.

Preferred visualization:

**alluvial/transition diagram or matched categorical strip**

for eligible parameter conditions across:

**0.5× → 1× → 2×**

Outcome states:

* clean reversal;
* retained original polarity;
* dual/multiple persistent domains;
* depolarized/other.

The dominant visual feature should be migration of conditions away from clean reversal as size increases.

Annotate:

* clean reversal at 1×: **71.9%**;
* clean reversal at 2×: **48.6%**;
* multiple polarized domains: **24.6% → 51.4%**.

### Panel F — What actually changes during failed reversal

Use trajectory-derived mechanistic summaries rather than another endpoint percentage.

Recommended measures:

* old-front amplitude versus time;
* new-front amplitude versus time;
* number of detected domains versus time.

Compare successful 1× and failed/dual 2× versions of the **same parameter condition**.

This should show that enlargement principally affects the competitive elimination/coarsening stage rather than the initial generation of the new front.

### Figure 1 conclusion

**A minimal conserved polarity system remains able to generate a front after enlargement, but increased size makes established polarity harder to replace.**

---

# Figure 2. A Cdc42–Rac–Rho network reproduces the reversal defect and reveals stronger size-dependent multipolarity

### Central result

Increasing biochemical complexity does not eliminate the size effect identified in the minimal system. Instead, the same loss of polarity decisiveness persists and a stronger tendency toward simultaneous competing fronts emerges.

The architecture of Figure 2 should deliberately echo Figure 1 so the reader can compare them almost panel-for-panel.

### Panel A — Three-species polarity network

Show the core interactions:

* Cdc42;
* Rac;
* Rho;
* active/inactive pools.

Include the implemented regulatory relationships:

* Rho inhibits Cdc42;
* Cdc42 promotes Rac;
* Rac inhibits Rho.

Identify Rac as the principal common polarity readout used in the matched phenotype analysis.

Include the same 0.5× → 1× → 2× geometric scaling icon from Figure 1.

### Panel B — Single-cue polarity remains robust to enlargement

Show matched Rac kymographs at 0.5×, 1×, and 2×.

Under each Rac kymograph, show late-time line profiles of:

* Cdc42;
* Rac;
* Rho.

This adds biological information unavailable from the one-species model and shows that the front corresponds to an organized Rho-family signaling state rather than a Rac-only artifact.

Beside these examples, show the 1× versus 2× endpoint transition matrix.

Annotate:

**487/507 matched conditions retain the same endpoint class (96.1%).**

### Panel C — Enlargement increasingly permits competing fronts

This should be much more visually prominent than the equivalent Figure 1 panel because this is where the richer network behaves differently.

Show aligned endpoint phenotype maps for:

* 0.5×;
* 1×;
* 2×.

Use identical parameter axes.

Persistent two-front states should be visually distinct and easy to follow as their region expands.

Alongside the maps, give a compact quantitative trend:

**5.0% → 8.6% → 13.6% persistent two-domain outcomes.**

Do not make this a standalone three-bar chart. A small line or dot summary next to the parameter maps is sufficient.

### Panel D — Spatial organization of a large-cell two-front state

Take one parameter condition that resolves to a single pole when small but retains two domains when enlarged.

Display:

* Rac kymograph;
* corresponding Cdc42 and Rho kymographs or selected time-point profiles;
* late-time active-species profiles.

The biological purpose is to demonstrate that the second pole reflects a coherent signaling structure, not threshold noise in the Rac classifier.

### Panel E — Reversal again deteriorates with enlargement

Repeat the matched reversal presentation from Figure 1.

Use one strong illustrative example:

* pre-existing left front;
* right-sided challenge;
* successful smaller-cell repolarization;
* incomplete or multipolar larger-cell response.

Show Rac as the main kymograph and Cdc42/Rho late-time profiles as supporting traces.

### Panel F — Population reversal outcomes

Show categorical matched outcomes across 0.5× → 1× → 2×.

Annotate clean reversal:

* **62.8% at 0.5×**;
* **52.7% at 1×**;
* **40.9% at 2×**.

Also show the increase in multiple-domain outcomes:

**36.1% → 50.2% from 0.5× to 2×.**

### Figure 2 conclusion

**The Rho-family network reproduces the size-dependent loss of repolarization seen in the minimal system and additionally makes persistent competing fronts increasingly likely as cells enlarge.**

---

# Figure 3. Increasing model complexity preserves a common size-response fingerprint

### Central result

This is the main cross-model synthesis figure.

It should answer directly:

> Which observable consequences of enlargement are model-specific, and which survive increasing biological complexity?

Do not organize this as six unrelated summary plots. Organize it as a **phenotype-by-model matrix**.

### Panel A — Model-complexity ladder

A horizontal visual progression:

**1-species polarity → Cdc42/Rac/Rho polarity → polarity + mechanics**

Underneath, show the common phenotype challenges:

**single cue | competing cues | cue reversal | migration**

Use this as a map for the remainder of the paper.

### Panel B — Raw phenotype response to cell size

Three rows:

1. single-cue polarization;
2. persistent competing fronts;
3. clean reversal.

Two initial model columns:

* 1-species;
* 3-species.

For each mini-panel, plot outcome prevalence at:

**0.5×, 1×, 2×**

using points/lines rather than bars.

Where possible, add faint matched-condition trajectories or the denominator `n`.

The reader should immediately see:

* near-flat single-cue response;
* stronger increase in dual polarity for the 3-species model;
* decreasing reversal in both models.

### Panel C — Size-effect fingerprint

Convert the preceding observations into a common effect metric.

For each phenotype plot the **1× → 2× paired effect**, with one point for each model.

Example:

| phenotype                | 1-species | 3-species |
| ------------------------ | --------: | --------: |
| single-cue polarization  |  ~+0.9 pp |  ~+2.8 pp |
| persistent dual polarity | ~+0.15 pp |  ~+5.0 pp |
| clean reversal           | ~−23.3 pp | ~−11.8 pp |

Use a zero-centered forest/dot plot.

This is where the reader should see the conserved sign structure:

* essentially neutral polarity competence;
* neutral-to-positive multipolarity;
* strongly negative reversal.

### Panel D — What enlargement does to the same conditions

Show the paired categorical transition structure for reversal in both models side-by-side.

Rather than only asking how the population prevalence changes, emphasize:

**conditions that successfully reverse at 1× frequently cease to do so at 2×.**

Annotate retention of successful reversal:

* 1-species: **68.3%** of 1× reversing conditions remain reversing at 2×;
* 3-species: **63.0%** remain reversing.

This is particularly persuasive because it establishes a genuine condition-level size effect.

### Panel E — Mechanistic synthesis

Small schematic grounded in the preceding panels:

**larger spatial domain → greater separation/persistence of competing polarity regions → slower or incomplete winner-take-all resolution**

Then distinguish:

* **front generation:** largely retained;
* **front selection:** size sensitive;
* **front replacement:** consistently impaired.

Keep this schematic simple. It should summarize the data rather than introduce a new hypothetical pathway.

### Figure 3 conclusion

**Across distinct polarity architectures, enlargement affects the resolution of competing spatial states much more strongly than the capacity to generate polarity itself.**

This figure supplies the paper's central phrase:

### **Cell enlargement reduces polarity decisiveness rather than polarity competence.**

---

# Figure 4. Coupling polarity to cell mechanics determines whether size-dependent polarity defects become migration defects

### Central result

The chemical models establish the common size-response fingerprint. The mechanochemical model asks whether that fingerprint persists once polarity must deform and move a cell.

This analysis should use the **per-reference-length Tier-2 data**. Do not collapse the ten reference lengths into a single distribution.

### Panel A — Mechanochemical model

Show:

* moving cell boundaries;
* polarity field along the cell;
* front-generated active force;
* rear/inward force where appropriate;
* Kelvin–Voigt mechanical element;
* cell center displacement;
* reference length `L0`.

For the 3-species implementation, indicate the Rac/Rho contribution to force.

The panel should make clear that cell size is now an explicit mechanical variable rather than only a fixed reaction–diffusion domain.

### Panel B — Representative small and large migrating cells

For the same biochemical/mechanical parameter condition, show two or three reference lengths.

Each condition should contain vertically aligned:

1. polarity kymograph;
2. left and right cell-boundary trajectories;
3. cell-center trajectory;
4. cell length versus time.

This allows the reader to see the causal sequence:

**polarity dynamics → deformation → displacement.**

Use equivalent examples for the 1-species and 3-species models if the layout permits.

### Panel C — Common polarity phenotypes versus reference length

Use reference length explicitly on the x-axis.

Plot, as separate aligned rows:

* single-cue polarization success;
* persistent dual-domain frequency;
* reversal success.

Columns:

* moving 1-species model;
* moving 3-species model.

Soft and stiff mechanics should be shown as separate traces or facets rather than averaged together.

This is the crucial test of whether the fixed-domain size-response fingerprint persists once the geometry can move.

### Panel D — Migration output versus cell size

For the same reference-length axis, show biologically interpretable migration outputs:

* net displacement;
* cumulative directional travel;
* mean polarity-projected speed;
* directional persistence/reversal if supported by the trajectories.

Prefer distributions/condition-level points at each reference length with a population trend overlay.

Do not use total travel alone, because high nondirectional movement and effective directional migration are biologically different.

### Panel E — Link polarity phenotype directly to migration phenotype

Stratify the migration outputs according to the polarity state reached by each simulation.

For example compare:

* clean single-pole runs;
* persistent dual-pole runs;
* failed reversal runs.

Plot net directional displacement or polarity-projected travel for each category.

This panel is important because it demonstrates whether the signaling phenotype identified in Figures 1–3 actually explains the movement phenotype.

### Panel F — Cross-model size sensitivity

For each common observable, calculate a normalized size effect across increasing model richness:

**fixed 1-species → fixed 3-species → moving 1-species → moving 3-species**

This can be a compact effect-size matrix.

The emphasis should be on sign and consistency rather than ranking arbitrary numerical magnitudes between models with different implementations.

### Figure 4 conclusion

**The mechanochemical model tests whether the conserved loss of polarity decisiveness remains biologically consequential once the cell must convert signaling asymmetry into directed movement.**

---

# Figure 5. Tension feedback creates a mechanochemical constraint on polarity and migration

I would separate tension feedback from Figure 4 unless the final dataset shows that it contributes little. It represents a genuinely new mechanism, not merely another mechanical parameter.

### Central result

Mechanical tension can feed back onto the polarity circuit and restrict sustained front activity. The question is whether this feedback modifies or amplifies the size dependence identified in the preceding figures.

### Panel A — Feedback mechanism

Extend the Figure 4 model schematic.

Show:

**cell extension → positive tension → increased deactivation/inhibition of front activity**

Include both implemented alternatives:

* global tension feedback;
* spatially propagated tension-associated feedback.

Make clear that both act on polarity chemistry but differ in spatial organization.

### Panel B — Representative mechanochemical trajectories

Use a matched condition at:

* `alpha_t = 0`;
* intermediate feedback;
* strong feedback.

For each display:

* polarity kymograph;
* cell length/strain;
* tension;
* center position.

A reader should see directly that increasing feedback shortens or extinguishes persistent front activity and therefore reduces displacement.

### Panel C — Length × feedback phenotype maps

This should be one of the strongest panels.

For each model and feedback mode, construct two-dimensional maps with:

* x-axis = reference cell length;
* y-axis = feedback strength.

Separate maps for:

* terminal polarized state;
* reversal success;
* persistent competing poles;
* directional displacement.

This exposes whether large cells and strong feedback interact rather than treating feedback strength and cell size as independent one-dimensional sweeps.

### Panel D — Global versus spatial feedback

Direct matched comparison.

Plot size-dependent migration curves under:

* no feedback;
* global feedback;
* spatial feedback.

Prefer net directional travel or polarity-projected travel.

The previously observed aggregate suppression can appear as a small annotation, not the main panel:

**1-species mean total travel**

* no feedback ~11.2 µm;
* global `αt=20` ~0.78 µm;
* spatial `αt=20` ~0.16 µm.

**3-species**

* no feedback ~9.37 µm;
* global `αt=20` ~1.95 µm;
* spatial `αt=20` ~0.35 µm.

The length-resolved behavior matters more than these global averages.

### Panel E — Does tension amplify the size effect?

For each model compute a size-effect metric separately under:

* no feedback;
* global feedback;
* spatial feedback.

Do this for:

* polarization;
* reversal;
* dual polarity;
* directional migration.

Display as an effect-size matrix or forest plot.

This panel answers the biologically important question:

**Does mechanics merely suppress motility overall, or does it specifically make enlarged cells more vulnerable?**

### Panel F — Representative mechanistic chain

Use one matched example to connect:

**large size → altered polarity competition → increased mechanical/tension state → inhibition of front activity → reduced directional displacement**

Only include arrows supported by the simulated trajectory.

### Figure 5 conclusion

**Mechanical feedback provides a mechanism through which changes in cell geometry can be converted into stronger defects in front persistence and migration.**

---

# Figure 6. Integrated model-richness summary

This figure may be unnecessary if Figures 3 and 5 already provide a clean synthesis. If used, it should be extremely compact and conceptual.

### Panel A — Common phenotype matrix

Rows:

* establish one front;
* suppress competing fronts;
* reverse polarity;
* sustain directional migration.

Columns:

* minimal polarity;
* Rho-family polarity;
* polarity + passive mechanics;
* polarity + tension feedback.

Each cell should show the observed direction of the enlargement effect together with a tiny quantitative indicator.

The important pattern should visually emerge:

**front formation ≈ preserved**

while

**front competition/reversal → progressively vulnerable**

and

**migration becomes vulnerable once mechanochemical coupling is present.**

### Panel B — Biological interpretation

Large-cell cartoon contrasting two situations:

**Stable external cue**

* large cell can still form a front.

**Changing/competing cues**

* spatially separated polarity programs coexist longer;
* repolarization is less decisive;
* directional migration becomes less robust.

### Final manuscript-level conclusion

**Increasing cell size does not primarily compromise the biochemical capacity to polarize. Across progressively richer models, its more reproducible consequence is impaired winner-take-all selection between competing polarity states. Network crosstalk strengthens the emergence of competing fronts, while mechanical feedback determines whether this signaling defect propagates into altered cell migration.**

---

# Supplementary figure plan

### Figure S1 — Single-species calibration and numerical validation

* local bistability structure;
* polarization regime maps;
* dependence on `h`, `mu`, `Lambda`, and other principal controls;
* conservation error;
* discretization/numerical checks.

### Figure S2 — Full single-species size-response maps

* complete single-cue classifications at 0.5×/1×/2×;
* complete simultaneous competition classifications;
* full reversal classification maps;
* reversal-time distributions.

### Figure S3 — Three-species calibration

* complete network parameter/regime maps;
* Cdc42/Rac/Rho profile examples;
* composite network coherence metrics;
* conservation checks.

### Figure S4 — Full three-species size-response maps

* endpoint transitions;
* dual-polarity parameter maps;
* reversal maps;
* reversal-time distributions;
* secondary network phenotypes.

### Figure S5 — Moving-domain mechanical validation

* soft versus stiff mechanics;
* passive relaxation;
* length/strain ranges;
* force normalization with reference length;
* representative no-feedback trajectories.

### Figure S6 — Full reference-length phenotype atlas

For every reference length and mechanics condition:

* polarization;
* domain number;
* reversal;
* displacement;
* speed;
* length;
* strain.

This is the exhaustive counterpart to the selected curves in Figure 4.

### Figure S7 — Tension-feedback diagnostics

* global feedback exposure;
* spatial feedback-field profiles;
* tension versus deactivation;
* global versus spatial comparisons;
* `alpha_t` sensitivity.

### Figure S8 — Numerical and classification robustness

* resolution sensitivity;
* integration/conservation diagnostics;
* phenotype-threshold sensitivity;
* examples around classification boundaries.

---

# Visual consistency requirements

The same visual language should be retained throughout the paper.

* **Cell size:** always left-to-right small → baseline → large.
* **Time:** always left-to-right.
* **Cell position:** same front/back orientation throughout.
* **Polarity-state labels:** identical terminology across models.
* **Single/dual/reversal icons:** reused verbatim.
* **Model richness:** progressively add components to the same base schematic rather than drawing unrelated diagrams.
* **Matched parameter conditions:** explicitly labeled as matched so readers understand that size, not parameter choice, changed.
* **Kymographs:** use common normalization where scientifically defensible; otherwise clearly indicate per-panel scaling.
* **Population summaries:** show actual condition-level structure whenever possible.
* **Representative examples:** chosen because they illustrate the population result, with their location in parameter space marked on the corresponding population map.
* **Main figures:** emphasize biology and model behavior.
* **Supplements:** carry exhaustive parameter sweeps, calibration, classification details, and numerical diagnostics.

The dominant visual rhythm should therefore be:

**schematic → representative dynamics → parameter/population structure → matched size transition → biological synthesis.**
