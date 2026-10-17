# Plan: Build an Experiment-Ready Manuscript Linking Cell Size, Polarity, Mechanics, and WGD Migration

## Overall objective

Develop a manuscript-shaped study that uses established polarity and mechanochemical models to generate experimentally discriminating predictions about migration in enlarged/WGD cells, then carries those mechanisms into a full 2D migration model that can be used to design and simulate the intended experiment.

The central question is not whether cell size affects polarity. That is already established theoretically.

The study should instead determine:

1. **Which polarity behaviors are robustly altered by increasing cell size across established model classes?**
2. **Whether long-range mechanochemical Rac–Rho coupling changes or compensates for those size effects.**
3. **How these polarity and mechanical effects translate into whole-cell migration in two dimensions.**
4. **Whether enlarged/WGD cells are predicted to become particularly dependent on long-range mechanochemical coordination.**
5. **Which experimentally measurable migration phenotypes most cleanly distinguish these mechanisms.**

The immediate endpoint is an **experiment-ready nascent manuscript**: theoretical and simulation figures substantially complete, the experimental hypotheses prespecified, and the proposed wet-lab experiment represented by a clearly defined missing-result panel supported by synthetic-data validation.

---

# Scientific premise

## Established size dependence of cell-polarity systems

Jilkine and Edelstein-Keshet, *A Comparison of Mathematical Models for Polarization of Single Eukaryotic Cells in Response to Guided Cues* (2011), compared several canonical 1D polarity models and explicitly examined changes in cell size.

Their enlargement experiment already captures the relevant physical null:

* the cell domain becomes longer;
* local reaction kinetics are unchanged;
* physical diffusion coefficients are unchanged;
* increasing domain size is mathematically equivalent, after coordinate rescaling, to reducing effective diffusion relative to reaction;
* total molecular material increases with domain size when concentrations are held fixed.

They showed that larger domains can resolve competing polarity sites more slowly and, in some models, support persistent competing fronts.

That result should therefore be treated as **established theory and a validation target**, not as a novel finding.

The opportunity is to extend this literature systematically toward:

* reorientation and polarity replacement;
* mechanochemical compensation;
* whole-cell migration;
* WGD-associated enlargement.

---

## Long-range mechanochemical polarity control

De Belly et al. 2023 and 2026 provide a complementary mechanism for cell-wide communication.

The central biological picture is approximately:

**local Rac-driven protrusion
→ mechanically transmitted cell-wide signal
→ distal Rho activation**

with a reciprocal arm:

**local Rho/contractility
→ mechanical communication
→ distal Rac activation**.

This provides a plausible alternative to relying entirely on diffusion-mediated communication across an enlarging cytoplasmic domain.

The important unresolved scaling question is:

> **Does long-range mechanochemical communication buffer the deterioration in polarity coordination predicted by local reaction–diffusion systems as cell size increases?**

A further experimentally important possibility is:

> **If long-range mechanical communication compensates for increasing scale, enlarged/WGD cells may become more dependent on that pathway.**

This creates a potentially testable size-specific vulnerability.

---

# Empirical motivation from polyploid migration

Two studies provide useful but incomplete experimental constraints.

## Jemaà et al. 2023 — stable tetraploid cells

Jemaà and colleagues compare stable diploid and tetraploid cancer-cell clones and report increased migration in tetraploid cells, including increased single-cell migration speed.

No reusable raw single-cell trajectory dataset appears to be available; the supplementary archive contains supplementary figures rather than raw tracking coordinates.

The study nevertheless establishes an important phenotype:

> stable tetraploidization can increase migration speed.

Its major limitation for mechanistic interpretation is that migration is only weakly decomposed into dynamical components. The paper does not provide a detailed account of directional persistence, turning behavior, polarity stability, state transitions, or front–back dynamics.

---

## Xuan et al. 2018 — giant polyploid MDA-MB-231 cells

Xuan and colleagues examine much larger polyploid giant cancer cells.

Their migration phenotype is different:

* lower instantaneous migration speed;
* increased directional persistence;
* long-timescale displacement that remains substantial because increased persistence partially compensates for lower speed.

They also report:

* major actin reorganization;
* altered RhoA/ROCK-associated cytoskeletal behavior;
* reduced intracellular particle MSD consistent with increased local mechanical constraint.

Raw trajectory coordinates do not appear to be deposited, but figures and methods contain enough quantitative information to make the migration results useful as model targets.

The appropriate synthesis is therefore not:

> polyploidy increases migration.

Instead:

> **polyploidization can alter distinct components of migration differently, particularly propulsion and directional persistence.**

This motivates a mechanistic framework capable of separating:

* speed / propulsion;
* polarity stability;
* turning;
* persistence;
* reorientation;
* migration efficiency.

---

# Work Package 1 — Unified 1D polarity and mechanochemical modeling framework

## Goal

Build a common modeling framework containing established polarity models and a modern mechanochemical polarity model, then use identical phenotype definitions and perturbation protocols to determine which consequences of cell size are robust across mechanisms.

The framework should be explicitly literature-based.

It should not consist of a bespoke progression from a simple model toward increasingly elaborate custom models.

---

## 1A. Reproduce the canonical Jilkine model repertoire

Implement the models compared in Jilkine and Edelstein-Keshet 2011 as faithfully as practical, including their published parameterizations and stimulation protocols.

The core repertoire should initially include:

* wave-pinning;
* Goryachev-type polarity;
* Otsuji-type polarity;
* LEGI-type sensing/polarization.

For each model, reproduce the major published behaviors before performing extensions.

Important validation targets include:

* single localized stimulation;
* competing opposing stimuli;
* gradient or stimulus reversal where applicable;
* published cell-size effects.

### Deliverable

A reproduction/qualification figure demonstrating that each implementation behaves consistently with the original literature.

This establishes that subsequent predictions are extensions of canonical models rather than artifacts of new parameterizations.

---

## 1B. Add the De Belly mechanochemical model to the same framework

Implement the published De Belly 2026 1D model, preferably using or directly adapting the released implementation rather than constructing a loose approximation.

Represent its major dynamical components explicitly:

* Rac;
* Rho;
* local antagonism;
* Rac-driven mechanical signaling to distal Rho;
* Rho-driven mechanical signaling to distal Rac;
* the relevant mechanical/tension variable or variables.

The De Belly model should be subjected to the **same phenotype-scoring framework** as the Jilkine models wherever biologically appropriate.

The purpose is not to treat mechanochemical signaling as an afterthought or separate extension.

It is one of the principal candidate mechanisms within the common model comparison.

---

## 1C. Establish a common phenotype panel

Define model-independent outcome measures that can be extracted from all relevant 1D systems.

Core phenotypes:

### Initial polarization

* does a localized cue establish a persistent dominant front?
* time to front establishment;
* polarization amplitude.

### Competing-front resolution

* does one front win?
* number of persistent fronts;
* time to resolution;
* coexistence duration.

### Polarity replacement / reorientation

After establishment of one pole and presentation of a competing/opposing perturbation:

* successful replacement probability;
* reversal latency;
* old-front decay time;
* new-front establishment time;
* coexistence duration.

### Dynamic tracking

Where feasible, apply repeated or periodically changing perturbations and estimate:

* response lag;
* phase lag;
* frequency at which the polarity system ceases to track environmental changes reliably.

This can be interpreted as a **polarity-response bandwidth**.

---

## 1D. Extend the published size analysis rather than redefining it

For the Jilkine models, use the published cell-size formulation as the starting point.

There is no need to introduce a separate “WGD enlargement protocol” where the published model already explicitly defines how domain size is varied.

The extension should instead ask:

> How do the previously characterized models behave across cell size for dynamical phenotypes that were not comprehensively examined in the original comparison?

Particular emphasis should be placed on:

* reversal;
* polarity replacement;
* transient coexistence;
* resolution times;
* repeated reorientation.

For the De Belly model, vary its explicit cell/domain size parameter while otherwise retaining the published model formulation.

The goal is to compare **size-response fingerprints**, not to force every model into an artificial identical mathematical transformation.

---

## 1E. Distinguish slower dynamics from changed qualitative behavior

This is a critical analysis requirement.

Increasing domain size at fixed physical diffusion coefficients naturally lengthens diffusion-mediated communication times.

Therefore every apparent “failure of polarity resolution” must be separated into:

1. **slower resolution with the same eventual state**, versus
2. **a genuine change in stable or effectively persistent polarity organization**.

For each relevant model:

* run simulations long enough to determine long-time outcomes;
* quantify resolution/reversal times explicitly;
* test whether differences collapse under appropriate diffusive time rescaling;
* separately report final-state changes and kinetic delays.

This distinction should be central to interpretation.

---

## 1F. Mechanochemical ablation and dependency analysis

Within the De Belly model, systematically weaken or remove the long-range pathways.

Compare:

* intact mechanochemical coupling;
* reduced Rac→distal-Rho coupling;
* reduced Rho→distal-Rac coupling;
* loss of both long-range pathways.

Perform these perturbations across cell size.

The key quantity is therefore not simply:

[
\text{effect of mechanics}
]

but:

[
\text{cell size} \times \text{mechanochemical coupling}.
]

Questions:

* Does intact coupling flatten the relationship between size and polarity-resolution time?
* Does weakened coupling reveal a progressively larger defect as cells enlarge?
* Is there a size beyond which long-range coupling becomes insufficient?
* Are some polarity phenotypes rescued while others remain size-sensitive?

### Central candidate result

A particularly interesting outcome would be:

> local biochemical polarity becomes increasingly difficult to coordinate with size, whereas long-range mechanochemical coupling buffers this effect; removing that coupling disproportionately destabilizes polarity in larger cells.

This would generate a direct experimental prediction.

### WP1 deliverables

1. Reproduced canonical Jilkine models.
2. Reproduced De Belly 1D mechanochemical model.
3. Common phenotype-analysis pipeline.
4. Cross-model size-response figure.
5. De Belly intact-versus-ablated size-dependency figure.
6. Explicit list of experimentally discriminating predictions to carry into 2D.

---

# Work Package 2 — Full 2D mechanistic migration model

## Goal

Build a full 2D cell-migration model in which intracellular polarity dynamics generate deforming-cell movement.

This is a required component of the project.

The purpose is to connect the 1D polarity mechanisms to the actual observables produced by migration experiments.

The current implementation direction is the M2072 Cellular Potts framework.
Its role is not to identify a uniquely optimal or universally realistic cell.
It should provide a small ensemble of defensible model cells for testing whether
WGD-like enlargement changes migration and dependence on long-range
mechanochemical coordination. Broad calibration, emulator, and mechanics
campaigns are useful insofar as they establish that ensemble or diagnose a
specific threat to the size comparison; they are not an independent endpoint.

Several experimental sources constrain different parts of this model and should
not be collapsed into one assumed universal target. De Belly data constrain the
mechanochemical Rac-Rho response in their assay context. HT1080 data constrain
free-migration morphology and dynamics. The de Seze experiments provide spatial
and kinetic tests of local Rho-pathway stimulation in RPE1 cells. Jemaà and Xuan
provide sparse polyploid migration observations. A parameter may legitimately
vary between cell types or assay contexts when that choice is declared. Any
claim that one fixed parameterization represents all of these systems requires
separate evidence.

---

## 2A. Use an established 2D migration framework

The model should be literature-derived rather than invented from scratch.

The current working framework is the M2072 Cellular Potts implementation with
explicit Rac-Rho and mechanochemical coupling. The alternatives below record
the selection logic and remain relevant if M2072 fails the model-adequacy gates;
framework comparison is not itself an ongoing project objective.

Candidate foundations include:

### Camley-type phase-field migration model

Advantages:

* directly couples reaction–diffusion polarity to deforming-cell movement;
* explicitly supports crawling and turning;
* relatively close conceptually to wave-pinning-style polarity;
* natural framework for cell size and shape.

### Plazen-type Rac/Rho + Cellular Potts model

Advantages:

* explicit Rac/Rho dynamics;
* whole-cell movement;
* comparison against experimental migration quantities such as speed, MSD and directionality;
* potentially closer to the Rac/Rho mechanochemical biology of interest.

### Other literature-derived alternatives

Additional models may be considered if they provide:

* explicit 2D cell shape;
* mechanistic intracellular polarity;
* experimentally meaningful movement;
* computational tractability for parameter sweeps;
* a defensible way to vary cell size.

### Framework criterion

Choose the model that best supports the central question:

> **How does increased cell size interact with polarity and mechanochemical coordination to alter whole-cell migration?**

---

## 2B. Reproduce the chosen 2D model before modification

Before adding WGD/size questions:

* reproduce published morphology;
* reproduce migration behavior;
* reproduce speed/persistence ranges where reported;
* verify numerical stability;
* reproduce key published perturbations.

The model should have a clear validation baseline before being used as a hypothesis generator.

### Defensible parameterizations and stopping rule

A parameterization is suitable for the WGD experiment when it passes separate
model-adequacy checks rather than only obtaining a favorable aggregate loss.
At minimum, it should:

* produce recognizable cells within the broad morphological support of the
  target population, including plausible shape variation over time;
* produce credible free migration and responses to relevant polarity or
  mechanochemical perturbations;
* preserve the intended signs and spatial organization of Rac, Rho,
  protrusion, contraction, and long-range coupling;
* remain numerically stable, connected, and adequately resolved across the
  sizes and perturbations under study;
* be evaluated on held-out observations or simulations where practical.

Static shape agreement alone is insufficient because cell shape is part of the
CPM mechanics and directly affects migration and tension. Conversely, exact
contour matching is not required. Retain a small ensemble spanning plausible
baseline behaviors rather than selecting a single apparent optimum. Selection
and acceptance must not use the direction or magnitude of the later WGD result.

Parameterization work is sufficient when this ensemble passes the declared
baseline gates and further fitting is unlikely to change interpretation of the
primary size-by-coupling comparison. If no parameterizations jointly produce
plausible morphology, migration, and mechanochemical responses, reconsider the
model structure or intended scope instead of assuming a larger search will fix
the incompatibility.

---

## 2C. Introduce cell-size variation within the 2D framework

Vary reference cell size using the model's native geometry and conservation rules.

Preserve the biological quantities that the original model interprets as concentrations or local kinetic parameters unless a model-specific reason requires otherwise.

Freeze and audit this transformation before examining its biological outcome.
State which geometric, mechanical, biochemical, noise, and numerical quantities
change and which remain fixed. Size-specific recalibration is permissible only
when it has a stated biological interpretation; otherwise it risks removing or
creating the size effect being tested.

The purpose is to identify how enlargement changes:

* polarity;
* number and stability of protrusions;
* cell shape;
* speed;
* persistence;
* turning;
* displacement.

Importantly, 2D enlargement should allow phenomena impossible in 1D, including:

* lateral front formation;
* multiple spatially separated protrusions;
* asymmetric turning;
* branching;
* shape-dependent mechanical communication.

---

## 2D. Incorporate long-range mechanochemical feedback

The 2D model should include a mechanochemical feedback mechanism motivated as directly as possible by the De Belly Rac–Rho framework.

The aim is not necessarily to reproduce every molecular component of the De Belly finite-element model.

The required functionality is:

* Rac-associated protrusive activity creates a cell-scale mechanical signal;
* that signal contributes to distal Rho activation/inhibition of competing fronts;
* reciprocal Rho/contractile activity can influence Rac elsewhere in the cell;
* coupling strengths are independently manipulable.

Then repeat the size analysis under:

* intact coupling;
* partial inhibition;
* strong inhibition.

The key prediction remains a size × coupling interaction.

---

## 2E. Quantify whole-cell migration phenotypes

For every simulated trajectory, extract the same statistics intended for the experiment.

Minimum outputs:

* speed;
* path length;
* displacement;
* MSD;
* directional persistence;
* velocity autocorrelation;
* turning-angle distribution;
* orientation autocorrelation;
* protrusion number;
* protrusion lifetime;
* front-switching frequency;
* cell elongation;
* shape variability.

This provides the direct bridge from intracellular mechanism to experimental analysis.

---

## 2F. Test whether established polyploid phenotypes occupy different mechanistic regimes

Use Jemaà and Xuan as held-out quantitative or qualitative comparison targets.
Their WGD outcomes should not be used to select the baseline parameterizations
whose predictions are being assessed.

### Jemaà-like phenotype

Target:

* increased speed in stable tetraploid cells.

Ask which model parameter changes can produce higher speed without necessarily altering persistence.

### Xuan-like phenotype

Target:

* reduced speed;
* increased directional persistence;
* altered long-timescale MSD;
* increased mechanical/cytoskeletal constraint.

Ask whether enlargement combined with stronger mechanical coordination or altered contractility can generate this regime.

The objective is not necessarily to uniquely infer mechanisms from these papers.

Rather, ask:

> Can a common established migration framework explain how polyploid cells can occupy qualitatively different speed–persistence regimes?

If so, this provides a useful reconciliation of apparently divergent published observations.

---

## 2G. Identify experimentally discriminating predictions

The 2D model should ultimately reduce the candidate mechanisms to a small number of predictions.

For example:

### Hypothesis A — enlargement mainly affects propulsion

Prediction:

* speed changes strongly;
* persistence changes weakly;
* perturbing long-range coordination produces little size-specific effect.

### Hypothesis B — enlargement creates a polarity-coordination burden that is mechanically compensated

Prediction:

* basal migration may remain normal or increased;
* persistence/protrusion organization depends increasingly on mechanochemical feedback;
* inhibiting that feedback disproportionately affects larger cells.

### Hypothesis C — enlargement itself stabilizes mechanically persistent migration

Prediction:

* enlarged cells may move more slowly but maintain direction longer;
* inhibition of actomyosin/mechanochemical feedback selectively reduces persistence.

### WP2 deliverables

1. A streamlined, reproducible CPM parameterization and qualification workflow.
2. A small ensemble of defensible model cells selected without their WGD outcomes.
3. An explicit, audited WGD scaling rule.
4. Size-dependent migration and mechanochemical-ablation simulations across the ensemble.
5. Held-out comparison with literature-derived WGD and perturbation observations.
6. Speed/persistence/MSD/protrusion phenotype maps and a small set of experimentally distinguishable predictions.

---

# Work Package 3 — Explicit experimental design and synthetic validation

## Goal

Design the actual experiment before data collection and demonstrate that the intended inference pipeline can identify the model-predicted effects.

The primary experiment should be experimentally straightforward.

A rapid cue-switching system is optional rather than mandatory.

---

## 3A. Core factorial experiment

Use matched lower- and higher-ploidy/WGD populations under:

1. control treatment;
2. perturbation of mechanochemical polarity coupling.

Conceptual design:

|                        | Control | Mechanochemical perturbation |
| ---------------------- | ------- | ---------------------------- |
| Lower ploidy / smaller | tracks  | tracks                       |
| Higher ploidy / larger | tracks  | tracks                       |

Use multiple matched cell-line pairs where feasible.

---

## 3B. Mechanochemical perturbation

Candidate intervention points derived from De Belly include:

### mTORC2-dependent arm

Possible acute pharmacological pilot:

* KU-0063794.

More specific follow-up:

* RICTOR perturbation;
* mSIN1 perturbation.

### Contractile arm

Possible perturbations:

* Y-27632 / ROCK inhibition;
* myosin inhibition.

The pharmacological experiments are attractive as rapid pilots but should not be interpreted as uniquely blocking long-range communication because these pathways have broader effects on migration.

The crucial result is therefore not:

> drug reduces migration.

It is:

> **drug sensitivity depends on cell size.**

---

## 3C. Measure size directly

Ploidy should not be the sole mechanistic variable.

For individual tracked cells, obtain or infer physical-size measurements where feasible:

* projected cell area;
* nuclear area;
* calibrated volume where available.

The principal statistical analysis should retain measured size as a continuous predictor.

This enables direct testing of:

[
\text{migration phenotype}
\sim
\text{size}
\times
\text{perturbation}.
]

Ploidy remains an important grouping variable and biological context, but size is the variable most directly connected to the modeling hypothesis.

---

## 3D. Migration endpoints

The experiment should extract multiple independent dimensions of movement.

At minimum:

* instantaneous speed;
* average speed;
* path length;
* net displacement;
* MSD;
* directional persistence;
* velocity autocorrelation;
* turning-angle distributions;
* track duration.

If morphology allows:

* protrusion number;
* protrusion lifetime;
* cell elongation;
* orientation;
* shape fluctuations.

This is necessary to distinguish Jemaà-like higher-speed phenotypes from Xuan-like lower-speed/higher-persistence phenotypes.

---

## 3E. Synthetic end-to-end experiment

Use the full 2D model to generate synthetic tracks under the competing hypotheses.

Then mimic the intended imaging experiment:

* finite frame rate;
* finite observation duration;
* field-of-view loss;
* measurement noise;
* variable cell size;
* biological heterogeneity;
* multiple wells;
* multiple experimental replicates.

Run exactly the intended tracking and analysis pipeline on these synthetic data.

The objective is to establish:

1. whether the migration summaries recover the known simulated mechanisms;
2. which statistics best distinguish competing models;
3. how many tracks are required;
4. how many wells/biological replicates are required;
5. whether the size × perturbation interaction can be estimated reliably.

This creates a direct bridge:

[
\text{mechanistic model}
\rightarrow
\text{simulated experiment}
\rightarrow
\text{analysis pipeline}
\rightarrow
\text{statistical inference}.
]

---

## 3F. Statistical framework

Primary analyses should explicitly test interactions.

Candidate structure:

[
Y_{ijk}
=======

\beta_0
+
\beta_1,\text{size}
+
\beta_2,\text{perturbation}
+
\beta_3,(\text{size}\times\text{perturbation})
+
\text{line/experiment/well effects}.
]

Where appropriate, ploidy can be included independently:

[
Y
\sim
\text{size}
\times
\text{perturbation}
+
\text{ploidy}
+
\text{line}
+
\text{random effects}.
]

For trajectory-level quantities, hierarchical models should account for:

* cells nested within wells;
* wells nested within experiments;
* cell-line background.

The synthetic-data work should establish power and estimator performance before experiments begin.

---

# Work Package 4 — Build the experiment-ready manuscript

## Goal

Produce a manuscript that already makes a coherent scientific argument before the new wet-lab data arrive.

The experimental section should appear as a deliberately unresolved but precisely specified prediction rather than a vague future direction.

---

## Figure 1 — Canonical polarity models and size-dependent dynamics

Contents:

* reproduction of the Jilkine model repertoire;
* recovery of published competing-front size effects;
* explicit acknowledgement that these effects are established;
* extension to reversal, coexistence duration, resolution time and dynamic reorientation.

Main question:

> **Which polarity dynamics consistently become more difficult as cell size increases?**

---

## Figure 2 — Mechanochemical compensation in the unified 1D framework

Contents:

* De Belly model alongside the canonical biochemical models;
* cell-size sweeps;
* intact versus weakened long-range Rac–Rho coupling;
* resolution/reversal times;
* interaction between size and coupling.

Main question:

> **Does long-range mechanochemical communication buffer the polarity-scaling problem?**

This is likely to contain the principal new 1D modeling result.

---

## Figure 3 — Whole-cell consequences in a full 2D migration model

Contents:

* validated 2D migration model;
* small versus large cells;
* intact versus perturbed mechanochemical coupling;
* representative cell shapes and tracks;
* speed, persistence, MSD, turning and protrusion dynamics.

Overlay or otherwise relate the simulated phenotype space to the existing Jemaà and Xuan observations.

Main question:

> **How do size-dependent polarity and mechanical effects manifest as experimentally measurable migration phenotypes?**

---

## Figure 4 — Experimental prediction and validation

Before wet-lab data:

* experimental schematic;
* low/high ploidy × control/perturbation design;
* synthetic 2D trajectories generated from the mechanistic model;
* simulated image/track sampling if appropriate;
* recovered statistical effects;
* power analysis;
* explicit predicted outcome ranges.

After experimentation:

replace the synthetic-result panels with the real WGD migration data while preserving the same analysis structure.

Main question:

> **Are enlarged/WGD cells disproportionately dependent on mechanochemical polarity coordination?**

---

# Experiment-ready decision criteria

The project is ready for wet-lab execution when:

1. The canonical 1D literature models have been reproduced.
2. The De Belly mechanochemical model is integrated into the same comparison framework.
3. A robust size-dependent polarity phenotype has been identified.
4. The effect of long-range mechanochemical coupling on that size dependence is known.
5. A small ensemble of 2D parameterizations produces credible cell morphology,
   shape dynamics, migration, and mechanochemical responses without selection on
   WGD outcomes.
6. The audited WGD transformation generates experimentally measurable and
   parameterization-aware differences between the candidate mechanisms.
7. At least one experimentally practical mechanochemical perturbation exists.
8. Synthetic experiments demonstrate that the proposed assay can identify the predicted size × perturbation interaction.
9. Required sample size and replication have been estimated.
10. Figures 1–3 of the manuscript are substantially complete and Figure 4 contains a concrete experimental-design placeholder.

---

# Immediate implementation sequence

## Phase 1 — Common 1D framework

> **Frozen-reproduction decision (2026-08-29).** Exact recovery of the Jilkine
> figure trajectories is no longer a Phase 1 objective. The paper does not
> document a reproducible symmetry-breaking mechanism and omits other
> numerically consequential details. The current equation-level,
> parameter-documented qualification is therefore frozen and will be labelled
> as a qualitative reconstruction rather than tuned further. The optimized De
> Belly port is likewise accepted as the Phase 1 screening implementation; a
> separate native-Gridap trajectory cross-check is not required. Phase 1 effort
> now starts with phenotype hardening and production execution.

1. Freeze the documented Jilkine equation-level implementations and qualification assets.
2. Retain the four canonical models as literature-derived comparison classes without further figure tuning.
3. Treat the existing size qualification as qualitative historical context rather than an exact reproduction target.
4. Use the accepted optimized De Belly 1D port as the mechanochemical screening model.
5. Implement a common phenotype-scoring system.
6. Extend all relevant models to reversal, resolution-time and dynamic-reorientation assays.
7. Separate kinetic slowing from genuine changes in final polarity state.
8. Perform mechanochemical ablations across size.
9. Identify the strongest experimental prediction.

## Phase 2 — Full 2D migration model

> **Frozen Tao baseline decision (2026-08-29).** A prespecified 16-candidate
> calibration screen followed by an exact, disjoint-seed GPU pilot did not
> recover the V2 deformation target. The initial V3 check also failed under its
> provisional protocol, but a later audit of `main.pdf` showed that its
> 60-minute window and 10-second speed sampling were not source-supported; V3
> is therefore indeterminate rather than a quantitative failure. The full
> 140-run ensemble was stopped for futility because V2 failed independently.
> The Tao implementation is
> therefore frozen as a qualitative literature-derived baseline: it supports
> directionality and qualitative tension contrasts, but not claims about
> absolute speed or published V2/V3 distributions. Subsequent Phase 2 work must
> hold the selected missing protocol inputs fixed and use relative contrasts.
> See `migration2d/results/qualification/v2_v3/baseline_freeze.json`.

10. Use M2072 as the working 2D polarity–migration model while retaining model
    adequacy, rather than sunk effort, as the criterion for continuing with it.
11. Distill and preserve the M2072 parameterization and qualification workflow,
    including the distinct roles and limitations of its experimental targets.
12. Select a small ensemble with credible baseline morphology, dynamics,
    migration, and mechanochemical responses, without consulting WGD outcomes.
13. Freeze and audit the controlled WGD size transformation.
14. Test size, mechanochemical coupling, parameterization, and stochastic
    variation together using prespecified migration and morphology outputs.
15. Compare the held-out size results with Jemaà, Xuan, and other relevant
    literature-derived observations at the level those data support.
16. Determine whether failures reflect parameter uncertainty, incompatible
    experimental contexts, or inadequate model structure.
17. Identify the migration statistics and feasible perturbations that best discriminate mechanisms.

## Phase 3 — Simulated experiment

18. Formalize the actual control/perturbation × ploidy/size assay.
19. Simulate realistic experimental datasets using the 2D model.
20. Introduce observation noise, track censoring and biological heterogeneity.
21. Run the complete intended statistical analysis.
22. Estimate power and required replication.
23. Freeze primary hypotheses and endpoints.

## Phase 4 — Manuscript assembly

24. Produce publication-quality Figures 1–3.
25. Produce the provisional synthetic-data Figure 4.
26. Draft the manuscript around the distinction between established theory, new modeling synthesis and experimental prediction.
27. Use the resulting manuscript package to justify immediate experimental execution.

---

# Target scientific story

The intended manuscript should ultimately support a claim along the following lines:

> **Increasing cellular scale creates a long-range polarity-coordination problem predicted by established reaction–diffusion models. Rapid mechanochemical Rac–Rho communication can compensate for this scaling constraint, but in doing so may create an increasing dependence on cell-wide mechanical coordination as cells enlarge. A full two-dimensional migration model predicts how this interaction alters speed, persistence, turning and protrusion dynamics, providing a mechanistic framework for the distinct migration phenotypes reported in stable tetraploid and giant-polyploid cancer cells. These predictions motivate a direct test of whether WGD-associated enlargement increases sensitivity to perturbation of long-range mechanochemical polarity control.**

The novelty should therefore come from the integrated progression:

**established size-dependent polarity physics
→ mechanochemical compensation
→ full 2D migration consequences
→ testable WGD-specific dependency**

rather than from claiming discovery of the generic effect of cell size on polarity.
