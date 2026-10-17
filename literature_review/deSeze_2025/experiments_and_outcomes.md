# de Seze et al. — spatial, morphological, and kinetic experiment reference for CPM calibration

**Paper:** Jean de Seze, Maud Bongaerts, Benoit Boulevard, and Mathieu Coppey. *Optogenetic control of a GEF of RhoA uncovers a signaling switch from retraction to protrusion.* eLife 12:RP93180. Version of Record published 27 May 2025. DOI: 10.7554/eLife.93180; version-specific DOI: 10.7554/eLife.93180.4.

**Purpose:** Extract stimulus–response observations that can constrain a two-dimensional Cellular Potts Model (CPM), especially the sign and geometry of local boundary motion, coordinated movement elsewhere in the cell, and reporter kinetics and localization. Here, “translation” means cell-body displacement, not protein synthesis.

## Scope, provenance, and access actually achieved

This report follows the supplied de Belly experiment reference in separating experimental conditions, cell morphology, spatial signals, kinetics, and proposed screening criteria. It does not transfer the de Belly cells, light geometry, bleb statistics, reporter amplitudes, or mechanisms into this study.

**R** = reported in article text or caption. **P** = a number or label printed on a figure. **F** = my visual description or approximate reading of an image/plot; not digitized data. **D** = my calculation from an identified source value. **M** = a proposed model observable or screening rule, not an experimental result. Values prefixed by **≈** are approximate. Numerical ranges around F estimates describe reading precision, not confidence intervals, standard errors, or observed population ranges unless explicitly stated.

The final eLife full text, Methods, Appendix 1, and movie captions were accessible. I directly inspected **all six main figures in the preserved eLife reviewed-preprint v2 PDF**, including the actual microscopy montages and graphs. I additionally inspected the **final Figure 1—figure supplement 1** and **final Appendix 1—figure 1** as standalone images. Main-figure visual estimates are therefore explicitly labeled **v2**; final text supplies the principal interpretation and protocol. This is not a claim that every pixel or dataset in the final main figures was independently verified as unchanged. [S1–S4]

**The seven movies were not successfully retrieved for frame-by-frame inspection.** Nor were the numerical source-data ZIPs or raw microscopy files downloaded and analyzed. Other supplementary/appendix figures are used only to the extent their text, captions, or author discussion was accessible; they are not assigned invented image measurements. The earlier blanket access assessment overstated verified access to these assets. The article states that numerical data and plotting code accompany figures and that raw images are deposited under BioImage Archive accession **S-BIAD1842**; that availability statement is not equivalent to successful retrieval here. [S1, Data availability]

### Source key

**S1 — Final article, full text, Methods, Appendix, and captions:** https://elifesciences.org/articles/93180

**S2 — Visually inspected reviewed-preprint v2, DOI 10.7554/eLife.93180.2, 11 December 2024:** titled *Control of cell retraction and protrusion with a single protein*. PDF mirror: https://www.researchgate.net/publication/392029340_Control_of_cell_retraction_and_protrusion_with_a_single_protein/fulltext/683115cfd1054b0207f1cf90/Control-of-cell-retraction-and-protrusion-with-a-single-protein.pdf . Main figures occur on PDF pages 4, 6, 8, 10, 12, and 13; Figure 6 caption is on page 14. References to PDF pages below are one-based.

**S3 — Visually inspected final Figure 1—figure supplement 1:** https://iiif.elifesciences.org/lax%3A93180/elife-93180-fig1-figsupp1-v1.tif/full/1500%2C/0/default.jpg

**S4 — Visually inspected final Appendix 1—figure 1:** https://iiif.elifesciences.org/lax%3A93180/elife-93180-app1-fig1-v1.tif/full/1500%2C/0/default.jpg

**S5 — Final-paper full-text mirror, including mathematical values omitted by the eLife HTML parser:** https://www.researchgate.net/publication/392128140_Optogenetic_control_of_a_GEF_of_RhoA_uncovers_a_signaling_switch_from_retraction_to_protrusion . Its final PDF exceeded the browser's file-size limit; its parsed text was accessible.

**S6 — Author responses embedded in S2, especially PDF pages 27–28 and 33–34.** These are identified separately from the main paper and figures.

---

# 1. Experimental system and the observation model

## 1.1 Cells, baseline state, and the exact comparison

**R.** The principal system is isolated, adherent **hTERT RPE1** cells, transiently expressing optogenetic constructs. HeLa cells provide a secondary demonstration of the two morphological responses. There is **no induced polyploidy comparison, matched ploidy-clone experiment, or controlled cell-size manipulation** in this paper. The central comparison is low versus high basal expression of the PRG DH–PH construct, and additional experiments manipulate its domains, downstream signaling, or recruitment waveform. [S1, Results; Materials and methods]

**R.** Experiments begin at least **30 hours after DNA transfection**. RPE1 cells are cultured in DMEM/F-12 with 10% serum and imaged at 37°C and 5% CO₂. Cells are scattered at low density on coverslips rather than confined as round suspended cells. [S1, Cell culture; Transfection; Imaging; Cell finder]

**F, v2.** Initial outlines are heterogeneous: compact rounded/ovoid cells, spread irregular cells, and substantially elongated cells all occur. Representative compact cells in Figure 1G/I have visible dimensions of roughly **30–50 µm**, judged against 10-µm scale bars; other examples are longer. This is a soft image-scale description, not a cell-size distribution. Some cells have pre-existing narrow extensions or ruffled edges before activation. A universally circular 10–15-µm initial condition would not reproduce these images. [S2, Fig. 1G/I and Fig. 3A/D/G]

**R.** High-optoPRG/protrusion-associated cells have approximately **1.5 times the initial projected area** of low-optoPRG/retraction-associated cells and greater spontaneous membrane activity before stimulation. This is an expression-associated baseline difference, not proof of a direct size effect or a 1.5-fold diameter/volume difference. [S1, “Cell phenotype … depends on the cytosolic concentration”; Fig. 2—figure supplement 1]

## 1.2 What light controls

**R.** Membrane-associated **mVenus–iLID–CAAX** recruits a cytosolic **PRG/ARHGEF11 DH–PH–SspB** construct upon blue illumination. Most experiments track the actuator with RFPt; an iRFP version is used when needed for spectral compatibility. The two parts are expressed from a P2A-containing construct, approximately at a one-to-one expression ratio. Recruitment changes local GEF availability; it is **not a direct clamp of RhoA activity**. [S1, first Results section; Plasmids]

**R/P.** Stimulation uses a DMD-generated rectangular/square region intersecting one pole of the cell. The illumination wavelength is **440 ± 10 nm**. The iLID–SspB dissociation timescale is approximately **20 s**. Many displayed experiments use repeated pulses every **30 s**; the frequency-switch experiment compares **30-s and 15-s periods**. Observation windows are generally **30–60 min**, while early reporters are compared over seconds to minutes. [S1, Fig. 1A/C; Fig. 5J; Fig. 6F; Imaging]

**F, v2.** The illuminated rectangle typically spans a substantial portion of one end of the cell. It is not a 1–2-µm circular spot. Its overlap with the cell changes as the boundary retracts or extends. Exact physical ROI dimensions and pulse energy should be recovered from each experiment rather than assigned universally from the montage. [S2, Figs. 1, 3–6]

**R.** The paper reports approximately **3–10-fold accumulated recruitment** in its initial Results and an experimentally accessible **1.1–3-fold single-pulse input increase** in the model section. The author response clarifies the former as accumulation after repeated activation, around five minutes with 30-s spacing, and the latter as the effect of an individual pulse. These must not be combined into a single dose scale or interpreted as 3–10-fold RhoA activation. [S1, first Results and model sections; S6, p. 28]

## 1.3 What the images and plotted quantities actually measure

**R.** The imaging combines DIC and fluorescence with a **TIRF configuration**, using 60×, NA 1.45 objectives. Fluorescence therefore samples the substrate-proximal region rather than a complete three-dimensional membrane surface. [S1, Imaging]

**M.** For a CPM, distinguish the projected cell mask, the basal membrane-associated reporter field, and any abstract perimeter-only signaling field. A broad bright basal patch must not automatically be encoded as a narrow cortical crescent or an entire circumferential ring. Changes in cell–substrate geometry may also affect TIRF intensity; intensity is not a direct force measurement.

**R.** The Methods define “surface displacement” as the **area of intersection between the activated region and the segmented cell**. Thus, with cell mask C(t) and illumination mask M:

`A_ROI(t) = area[C(t) ∩ M]`

The normalized plots divide this by the initial mean area in that region. This is an area/occupancy measure, **not centroid displacement, cell speed, three-dimensional membrane area, or a measurement of whole-cell growth**. [S1, Data processing; Fig. 1E/F]

**R.** Biosensors are background-subtracted and processed as the mean fluorescence in the activated part of the cell divided by the mean in its non-activated part, then normalized to the pre-stimulus value:

`B(t) = {mean(I_background-corrected, C(t) ∩ M) / mean(I_background-corrected, C(t) \ M)} / baseline_ratio`

Consequently, **B = 1.2 means a 20% increase in this normalized local-to-rest ratio**, not necessarily 20% more active GTPase, 20% more total protein, or a membrane/cytosol ratio of 1.2. A changing denominator can contribute. [S1, Data processing]

**M.** Apply these same measurement operations to simulated images/fields before fitting amplitudes. Preserve the experimental ROI in laboratory coordinates within each stimulation segment rather than continually recentering it on the moving tip. For long retractions, shrinking ROI occupancy also requires a minimum-valid-area criterion; that criterion is a modeling/data-QC choice, not a reported threshold.

---

# 2. Local optoPRG at low basal expression: stimulated-side retraction

## What is done

**R.** Recruit the PRG DH–PH construct repeatedly to one pole of a low-expressing RPE1 cell. In the concentration-stratified analysis, approximately **95% of cells below about 40 arbitrary fluorescence units** retract in the activation region. This threshold is specific to the authors' fluorescence measurement and cannot be transferred directly into another microscope or molecular-concentration unit. Phenotype is assessed using normalized ROI area at **five minutes**. [S1, Fig. 2A and corresponding Results]

## Morphology and spatial organization

**F, v2.** The stimulated edge withdraws from the fixed rectangle; the local boundary becomes less extended and can flatten into a compact trailing edge. Dense actuator signal remains near the withdrawing cell edge rather than implying a new protrusive front. The cell is not necessarily globally round: the opposite side can remain irregular, lobulated, or extended while the stimulated end shrinks. [S2, Fig. 1G; Fig. 3D/G]

**R.** Blebs **or** protrusions can appear at the non-activated pole. The text does not establish that all such extensions are pressure-driven blebs. **F.** The images support heterogeneous opposite-side outlines, not one invariant smooth hemispherical bleb. They do not provide a reliable bleb number, neck width, lifetime, angular distribution, or reporter coating pattern for those distal structures. [S1, first Results section; S2, Fig. 1G]

## Movement and quantitative area targets

**F, final supplement.** The activated region loses projected cell area while the non-activated region gains it. Approximate plotted central-trace values are:

| Time after onset | Activated-region area change | Non-activated-region area change |
|---|---:|---:|
| 5 min | ≈−15 to −20 µm² | ≈+5 to +15 µm² |
| 10 min | ≈−35 to −45 µm² | ≈+25 to +35 µm² |
| 15 min | ≈−55 to −65 µm² | ≈+50 to +60 µm² |

These are visual readings of the plotted central curves, **not per-cell acceptance intervals**. [S3, panel C]

**R/F.** The spatially opposed changes support coordinated movement away from the stimulated end, rather than simple uniform loss of projected cell area. Their approximate balance is consistent with relatively little net area change over this interval; it does not establish exact area or volume conservation. [S1, first Results section; S3, panel C; S6, p. 33]

**F, v2.** Across Figure 1E, normalized ROI area declines with substantial cell-to-cell variation. Some trajectories approach zero after roughly 10–20 minutes; others retract only partly. Zero means the cell has largely left the selected region, not that the cell has disappeared. [S2, Fig. 1E]

## Minimal CPM criterion

**M.** Require a sustained negative change in area at the stimulated pole, usually accompanied by area gain elsewhere, and local contractile-reporter recruitment. Do not count whole-cell collapse or a purely local inward indentation with no coordinated opposite-side response as equivalent to the full phenotype. An occasional opposite-side bleb is permitted but is not a mandatory, numerically calibrated outcome for every cell.

---

# 3. Local optoPRG at high basal expression: stimulated-side protrusion

## What is done

**R.** The same DH–PH recruitment procedure is applied to cells with higher basal optoPRG. Approximately **85% above the roughly 40-a.u. threshold** extend into the illuminated region in the five-minute analysis. [S1, Fig. 2A]

## Morphology

**R.** The authors describe **ruffles and filopodia** at the activation site. **F, v2.** The stimulated sector expands through an irregular, frilled/digitated outline, sometimes with several small extensions along a broader protrusive base. It is not adequately represented as a single smooth pressure bleb. The overall cell can become more polarized, but pre-existing shape varies and there is no universal circular-to-teardrop transition. [S1, first Results section; S2, Fig. 1I; Fig. 3A; Fig. 5J]

**F, v2.** The red actuator signal is enriched in the illuminated, substrate-proximal protrusive region. In Figure 1 this is **the optogenetic construct**, not Rac, Cdc42, or RhoA reporter fluorescence. Reporter localization must instead be inferred from the explicitly labeled experiments in Figures 3 and 5. [S2, Fig. 1 caption and panels G/I]

**R/F.** The opposite, non-activated region commonly retracts. The useful target is an expanding activated end coordinated with withdrawal elsewhere, not unrestricted isotropic spreading. [S1, first Results section; S3, panel B]

## Quantitative area targets

**F, final supplement.** Approximate central-trace changes are:

| Time after onset | Activated-region area change | Non-activated-region area change |
|---|---:|---:|
| 5 min | ≈+60 to +70 µm² | ≈−55 to −70 µm² |
| 10 min | ≈+85 to +95 µm² | ≈−75 to −90 µm² |
| 15 min | ≈+100 to +110 µm² | ≈−80 to −100 µm² |

The dominant pattern is reciprocal redistribution of projected occupancy. There can be residual net area gain; do not impose exact equality of the two curves. [S3, panel B]

**F, v2.** The five-minute scatter in Figure 2A includes many protruding cells with normalized ROI area around **1.1–1.7**, while the population time courses extend higher at later times. This is an approximate visible range among many points, not a reported population interval. The normalized area can exceed two at later times because the denominator is only the initially occupied illuminated region. [S2, Fig. 2A; Fig. 1E]

## Minimal CPM criterion

**M.** Generate a finite-width, irregular protrusive sector at the stimulus, with positive local area change and frequently negative area change in the rest of the cell. The accompanying pathway should show early Cdc42 recruitment and transiently suppressed or delayed accessible RhoA, followed by delayed Rac1 recruitment. It must not require immediate strong Rac1 activation or a smooth bleb at the stimulated end.

---

# 4. Mixed responses, baseline ruffling, and repolarization

## 4.1 Mixed is not the same as inactive

**R.** Mixed cells can form local ruffles/filopodia while also appearing to contract and developing distal blebs, yet show little clear net movement. **F, final supplement.** The example retains a broadly similar position over 0–25 minutes, but its outline remains active and irregular. The local actuator signal stays within the lower illuminated region. [S1, first Results section; S3, panel A]

**M.** Score biochemical recruitment, local contour activity, signed area change, and centroid translation separately. A cell with moving ruffles but near-zero net displacement is not equivalent to a static cell with no signaling response. Figure 1's combined “mixed/no movement” category cannot be automatically split into separate numerical classes without the underlying annotations.

## 4.2 Baseline membrane activity differs before light

**R.** Future protruding/high-expression cells show more spontaneous membrane displacement than retracting/low-expression cells before stimulation. **F/P, v2.** Figure 2E's plotted central markers are approximately **2–3 versus 4–5**, with an axis labeled **µm²/min**. The authors call this membrane mean speed, but those units are not linear velocity. The figure supports greater basal membrane activity, not a claim that the cells migrate at 2–5 µm/min. [S1, Results accompanying Fig. 2E; S2, Fig. 2E]

**M.** Use this as a baseline contour-activity comparison. Do not substitute these values for centroid speed or assume that the exact Fig. 2E algorithm is fully specified by its label.

## 4.3 Stimulating the other pole does not usually reverse the local response class

**R.** After recovery, approximately **one hour after the first activation**, the opposite pole is illuminated. Almost all clearly protruding cells protrude again; clearly retracting cells retract again. Some change between clear and mixed phenotypes. **P, v2.** The associated Sankey diagram is labeled **N = 88**. [S1, Fig. 1G–I and corresponding Results; S2, Fig. 1H]

**F, v2.** The selected rectangle switches from one side of the cell to the other. A protrusion-following cell consequently changes which end advances; a retraction-following cell withdraws whichever end is activated. This distinguishes the **sign of the local response** from the **direction of motion in laboratory coordinates**. [S2, Fig. 1G/I]

**M.** Do not hard-code the original cell rear as permanently contractile or original front as permanently protrusive. Under the standard waveform, the same expression-conditioned response should be inducible at either end. This experiment is not a seconds-scale cue-reversal measurement: it includes a long recovery interval, and no precise turning angle or repolarization latency is extracted here.

---

# 5. Local actin, myosin, and RhoA signals: morphology-linked kinetics

## 5.1 What was measured

**R.** Separate reporter experiments use **Lifeact-iRFP**, **MRLC-iRFP** for myosin regulatory light chain, and a relocation biosensor for accessible active RhoA, labeled **RBD-2xTdTomato** in Figure 3. These are not simultaneous three-color trajectories in the same cell. The RhoA experiment uses an iRFP actuator because of spectral overlap and yields a cohort biased toward retracting cells because high expression is harder to obtain. [S1, “Two distinct signaling pathways …”; Fig. 3]

**R.** Differences in actin and myosin are already detectable in the first post-activation frame, approximately **30 s** after onset, whereas the morphological response becomes recognizable after approximately **1–2 min**. This is a sampling-limited statement, not proof of a precise 30-s biochemical latency. Pairwise comparisons are made at **60 s for Lifeact and RBD**, and **90 s for MRLC**. [S1, Fig. 3 and corresponding Results]

## 5.2 Retraction: a local actomyosin-rich contractile region

**R.** Actin polymerization rises promptly in the activation region, myosin recruitment follows, and the edge retracts over the following minutes. RBD rises strongly and early. [S1, Fig. 3]

**F, v2.** Lifeact is visibly organized into elongated/bundled basal structures, with enhanced signal near the stimulated end. MRLC develops a dense local patch/band adjacent to the withdrawing boundary. RBD becomes enriched over a substantial stimulated-side basal region, often with a brighter band near the edge. None of these examples establishes a uniform whole-cell ring or a narrow, isolated one-pixel rear spot. The reporter-bright region and contracting edge remain spatially associated as the outline moves. [S2, Fig. 3A/D/G]

## 5.3 Protrusion: early local relaxation, then actin-rich edge growth

**R.** Both Lifeact and MRLC initially decrease locally. Lifeact subsequently rises as the protrusion develops, whereas MRLC remains below baseline. The authors interpret the initial decrease as likely stress-fiber disassembly; that interpretation is not a direct measurement of filament turnover. RBD shows little or no significant increase during the first minute, followed by a delayed rise. [S1, Fig. 3]

**F, v2.** The protruding Lifeact example transitions from a broader basal network to enhanced signal near the new, irregular advancing margin. The signal is not confined to a spherical bleb dome. The MRLC-rich local structures fade rather than becoming the dominant bright region at the growing end. The RBD signal can later become visible within the protrusive region; therefore, a protruding cell should not be forced to have permanently absent RhoA activity. [S2, Fig. 3A/D/G]

## 5.4 Approximate kinetic targets from the actual graphs

All numbers below are **F, v2**, read from the plotted mean curves, with baseline normalized to 1. They are **local-to-rest reporter ratios**. The bounds are deliberately broad visual read-off ranges, not experimental errors. [S2, Fig. 3B/E/H]

| Reporter and phenotype | ≈1 min | ≈2 min | ≈5 min | Main temporal signature |
|---|---:|---:|---:|---|
| Lifeact, retraction | 1.10–1.15 | 1.15–1.20 | 1.13–1.20 | Early rise, then modest plateau/fluctuation |
| Lifeact, protrusion | 0.98–1.02 | 1.04–1.10 | 1.18–1.25 | Small early dip, delayed rise |
| MRLC, retraction | 1.05–1.15 | 1.20–1.35 | 1.70–1.90 | Sustained accumulation |
| MRLC, protrusion | 0.85–0.95 | 0.75–0.85 | 0.75–0.85 | Sustained local depletion |
| RBD, retraction | 1.25–1.40 | 1.75–1.95 | 2.70–2.90 | Strong early and continuing increase |
| RBD, protrusion | 0.98–1.05 | 1.15–1.25 | 1.50–1.65 | Delayed, smaller increase |

**M.** The joint discriminator is more informative than a single endpoint: early actin increase plus growing MRLC and RBD in retraction, versus early local relaxation, continued low MRLC, and delayed actin/RBD accumulation in protrusion. These curves should not be stitched into a supposedly simultaneous single-cell time series.

---

# 6. Expression-state manipulation: background PRG changes the response to the same local recruitment

**R.** The amount of recruited optoPRG overlaps across response classes: some low-expression retracting cells attain substantial absolute recruitment. Thus, peak recruited amount alone does not cleanly determine the response. [S1, Fig. 2C]

**R.** Adding a **non-recruitable cytosolic PRG DH–PH construct** while keeping optoPRG expression low raises the protruding fraction from approximately **0% to 45%** in the authors' comparison. Greater non-recruitable PRG expression is associated with the protruding response. [S1, Results accompanying Fig. 2D; S2, Fig. 2D]

**M.** Separate a cell's pre-stimulus PRG-conditioned state from the amplitude of the acute membrane input. A model that treats optoPRG as a context-free, positive RhoA source will miss this experiment. An effective baseline-state variable is acceptable for screening; the data do not require the modeler to assert that baseline differences arise exclusively from one biochemical pool rather than possible longer-term cellular adaptation.

---

# 7. Fast pulse responses and isolated PH recruitment

## 7.1 optoPRG can produce positive, biphasic, or negative RBD responses

**R/F, v2.** Figure 4A shows three representative pulse-response classes. At low expression, RBD follows recruitment with a delayed positive response. At intermediate/high expression, recruitment can cause an initial RBD dip followed by recovery and a slower increase. In a high-expression example, RBD repeatedly drops below baseline, approximately into the **0.85–0.95** range, rather than accumulating strongly above 1. [S1, Results accompanying Fig. 4A; S2, Fig. 4A]

**M.** Preserve the distinction between **total GTP-bound RhoA** and **RhoA available to bind effectors/reporters**. The paper's mechanism is sequestration of active RhoA by PRG's PH domain. A decrease in RBD does not by itself show conversion of all active RhoA to GDP, global loss of RhoA, or movement of RhoA to the opposite pole. [Mechanistic basis: S1, Fig. 4D; Fig. 6A and Appendix 1]

## 7.2 Isolated optoPH lowers the local RBD signal

**R.** Recruiting the PH domain alone can produce an immediate local decrease in the RhoA reporter. The mean effect is significant but the response is not conspicuous in every cell, consistent with dependence on the pre-existing RhoA state. [S1, Fig. 4B–D; Fig. 4—figure supplement 1 discussion]

**F, v2.** In the strong representative example, normalized optoPH recruitment rises to roughly **1.7–1.9** around 15–25 s, while RBD falls to approximately **0.65–0.75** around 25–40 s and recovers partway, toward **0.85–0.95** around 70–80 s. These are example-trace values, not a population mean. [S2, Fig. 4B]

**F, final appendix.** A different displayed fitting example has a smaller input peak, about **1.35–1.4**, and a shallower RBD minimum near **0.87–0.90** around 25–45 s. This directly illustrates why the strong Figure 4 example should not be imposed as the response of every cell. [S4]

**F, v2.** The optoPH example is an elongated/irregular adherent cell, and the rectangular stimulus covers one end. The corresponding signal decrease is local to that region in the image and kymograph. No reliable quantitative protrusion, opposite-side reporter increase, or cell-translation trajectory accompanies this reporter-only demonstration. [S2, Fig. 4C]

## 7.3 PH-mediated inhibition alone does not create a robust protrusion

**R.** In cells overexpressing background PRG, recruitment of isolated PH does not reproducibly generate either clear protrusion or clear retraction. **F, v2.** Five-minute normalized ROI areas cluster near 1, with less dispersion than optoPRG; individual deviations are still present. The result is not “every boundary is motionless.” [S1, Fig. 4H/I; S2, Fig. 4H/I]

**M.** A model should allow a local RBD decrease without automatically creating a front. Inhibiting the contractile branch is not sufficient to replace active engagement of the protrusive branch.

---

# 8. PH-binding mutant: recruitment without the usual protrusion switch

**R.** The PH-binding mutations **F1044A and I1046E** prevent PRG PH binding to RhoA-GTP. In the ordinary mutant-optoPRG experiment, the authors observe retraction rather than the usual high-expression protrusion. [S1, Fig. 4E/F]

**F, v2.** The representative high-expression mutant cell retains a compact cell body and withdraws the stimulated end over the montage; the intensely recruited construct does not produce the expanding ruffled front of wild-type high-expression optoPRG. Bright fluorescent structures at the edge are actuator signal and should not themselves be scored as an advancing lamellipodium. [S2, Fig. 4F]

**R.** To separate preconditioning from acute PH action, the authors also overexpress non-recruitable wild-type PRG DH–PH and then recruit mutant optoPRG. Most cells still retract, although a few at very high background expression show small ruffles. [S1, Fig. 4G]

**P/F, v2.** Figure 4G compares **N = 24 wild-type and N = 15 mutant** cells in the high-background-PRG comparison. The bars are approximately **80–85% protruding with wild type**, versus **10–15% protruding and about 65–70% retracting with the mutant**. These are visual bar estimates, not manually reconstructed exact counts. [S2, Fig. 4G]

**M.** The informative rule is loss of the strong high-expression protrusion preference when PH binding is disabled, including after wild-type background PRG preconditioning. Do not impose an absolute universal ban on any ruffle in the separate high-background mutant condition, because residual ruffling is reported.

---

# 9. Cdc42 and Rac1: an early local Cdc42 response, followed by delayed Rac1

## 9.1 Reporter identity and localization

**R.** Figure 5 uses **PBD-iRFP**, sensitive to both Cdc42 and Rac1, plus more selective **mCherry-WaspGBD** for Cdc42 and **mCherry-3xp67Phox** for Rac1. PBD is not an exclusive Rac reporter. Measurements are made in separate experiments, not simultaneous Cdc42/Rac/Rho imaging of the same protrusion. [S1, Fig. 5; Plasmids]

**F, v2.** In protruding cells, PBD and Cdc42 reporter enrichment occupies a finite illuminated/protrusive basal region, with brighter areas near the advancing margin. The signal is not demonstrably restricted to the extreme apex of one filopodium. Rac1 enrichment develops later at the stimulated protrusive pole; one illustrated cell is narrow and elongated, so its signal geometry should not be generalized into a broad front cap in every cell. [S2, Fig. 5A/D/G]

**F/M.** The montages support stimulus-aligned local enrichment but do **not** establish an exact angular domain width, a reporter decay length beyond the light boundary, complete coating of every filopodium, or enrichment around an entire distal bleb. They also do not supply a paired quantitative measurement of the opposite pole's reporter activity. Those would require spatial profiles or raw image analysis.

## 9.2 Temporal ordering

**R.** PBD increases early in both response classes. It continues rising in protruding cells but shows little further increase after about one minute in retracting cells. Cdc42 is activated immediately in both phenotypes, with much stronger continued activation in protruding cells. Rac1 activation is delayed and occurs in the protruding phenotype. [S1, Results accompanying Fig. 5]

**P/R.** Figure 5 compares PBD at **60 s**, and the specific Cdc42 and Rac1 reporters at **30 s**. Cdc42 is already increased at that early point, whereas Rac1 is not significantly increased. Thus, “Rac activates after Cdc42” is supported, but the exact biochemical delay is not resolved as a universal number of seconds. [S1, Fig. 5C/F/I]

## 9.3 Semi-quantitative mean-curve targets

All entries are **F, v2**, visual ranges around plotted mean traces, normalized to 1 at baseline. [S2, Fig. 5B/E/H]

| Reporter and phenotype | ≈1 min | ≈2 min | ≈5 min | Interpretation |
|---|---:|---:|---:|---|
| PBD, retraction | 1.06–1.10 | 1.04–1.09 | 1.02–1.08 | Small early rise without sustained growth |
| PBD, protrusion | 1.06–1.12 | 1.15–1.20 | 1.30–1.36 | Progressive enrichment |
| Cdc42, retraction | 1.00–1.04 | 1.00–1.04 | 1.05–1.09 | Early small response; modest late drift |
| Cdc42, protrusion | 1.06–1.11 | 1.12–1.17 | 1.20–1.26 | Early and sustained enrichment |
| Rac1, retraction | 0.98–1.02 | 0.96–1.00 | 0.98–1.02 | No sustained positive response |
| Rac1, protrusion | 1.00–1.04 | 1.06–1.10 | 1.16–1.22 | Delayed, progressive enrichment |

The modest late Cdc42 increase in the retracting mean trace is worth retaining: the correct target is **much weaker continued Cdc42 recruitment than in protrusion**, not the artificial requirement that Cdc42 equal exactly 1 at every later time.

**M.** Fold changes from different reporters are not a shared concentration or force scale. For example, a late RBD ratio of 1.5 and a Cdc42 ratio of 1.2 do not establish that the contractile branch dominates. Each reporter has its own baseline, response properties, and observation operator; compare like-with-like trajectories before constructing any cross-pathway balance.

## Minimal CPM criterion

**M.** A retraction-compatible parameter set must tolerate an early Cdc42/PBD transient; “any Cdc42 increase implies protrusion” would fail. A protrusion-compatible set should sustain local Cdc42, recruit Rac1 later, and couple that pattern to an irregular advancing edge. Do not require the Rac1 signal to precede the first visible morphological change.

---

# 10. PAK inhibition with IPA3: loss of protrusion and occasional phenotype reversal

## Protocol

**R.** Cells are first stimulated for **30 min**, then allowed **30 min after the end** of that activation before medium replacement with **5 µM IPA3**. They are stimulated again after **5 min of drug exposure**. Figure 5J depicts 30-s optogenetic pulse spacing. [S1, Drug assay; Fig. 5J]

## Morphology and movement

**R.** Many cells become round; none in the tested set protrude again. Retracting cells often retain their ability to retract, and some previously protruding cells now retract. **F, v2.** The representative cell shows a ruffled, expanding illuminated margin before drug and retreat of the illuminated boundary after drug. This is a change in the sign of local movement, not merely a dimmer version of the same advancing front. [S1, Fig. 5J/K; S2, Fig. 5J]

**P/R/F, v2.** The transition diagram is labeled **N = 28**, with no orange protrusion category after treatment. Approximately **55% retract and 45% are mixed/no movement** in the post-treatment bar. The robust constraint is the absence of renewed protrusion in this sample, not a universal exact 55:45 probability. [S2, Fig. 5K; S1, corresponding Results]

## Interpretation for a CPM

**M.** Suppress the PAK-dependent protrusive output while allowing contractile signaling and some retraction to remain. This is not a clean Cdc42 deletion experiment: the perturbation targets a downstream kinase, and the paper does not provide post-IPA3 Cdc42 and Rac1 reporter traces that could justify setting those activities to zero. Some rounding/non-response should be allowed rather than forcing every previously protruding cell to become a strong migrator in the opposite direction.

---

# 11. Temporal-waveform switching in selected individual cells

## What is controlled

**R.** At intermediate expression, the authors compare **stronger, less frequent pulses every 30 s** with **weaker, more frequent pulses every 15 s**, keeping total integrated illumination the same. The former can drive retraction, the latter protrusion. This changes both frequency and per-pulse strength; it is not a comparison of frequency alone at identical pulse amplitude. [S1, Fig. 6F and model Results]

**R.** The search targets cells with measured fluorescence **14–35 a.u.**. Only **3 of 90 screened cells** show repeatable switching between clear phenotypes under the tested protocol; most show mixed responses regardless of waveform. The approximately 40-a.u. split in the earlier cohort should not be forced into an identical universal numerical boundary. **D.** 3/90 is approximately **3.3%** of this selected screen. [S1, model Results]

## Shape and local area traces

**F, v2.** Figure 6F shows an elongated cell undergoing successive stimulation at labeled top/bottom positions. Several retraction episodes reduce normalized ROI area to approximately **0.5–0.8** over about 5–10 minutes, with one stronger episode approaching **0.3**. The clearer protrusion episode reaches approximately **1.5–1.6**; another is much smaller, near **1.1**. These are selected episodes in an illustrative switching cell, not population means. [S2, Fig. 6F]

**F, v2.** In the retraction montage, the stimulated end withdraws from the rectangular region while the remaining body retains an elongated geometry. In the protrusion montage, a rough/ruffled extension grows into the selected region. The plot resets the local time/region for successive episodes; it is not one continuous monotonic centroid trajectory. [S2, Fig. 6F]

## Strong screening discriminator

**M.** Permit a narrow parameter/state region in which recruitment waveform changes the sign of the morphological response. Do not demand that every cell switch, and do not treat the three switching cells as the entire population. Reproducing only this cherry-picked cell while making the remaining population uniformly switchable would miss an important observation.

**R, author response.** The authors also state that changing pulse intensity changed response strength without switching phenotype in their experiments. The main-paper intensity sweep in Figure 6—figure supplement 1B is explicitly a **model result**. These are different evidence levels: the author response supports a qualitative statement, but there is no independently inspected dose-response dataset here from which to assign a universal experimental no-switch probability. [S1, model Results; S6, p. 28]

---

# 12. Secondary systems and GEF specificity

## 12.1 optoLARG and optoGEF-H1 are not interchangeable with optoPRG

**R/P/F, v2.** Figure 1D compares optoPRG (**N = 74**), optoLARG (**N = 14**), and optoGEF-H1 (**N = 35**). optoLARG has no protrusion category in this sample; approximately 80% are retracting and the rest mixed/nonmoving. optoGEF-H1 has many round or unresponsive cells; only about 40% show distinct morphological responses, mostly retraction, with a smaller protruding subgroup. [S1, first Results section; S2, Fig. 1D]

**R.** The authors note greater difficulty obtaining healthy/highly expressing cells with optoLARG and optoGEF-H1. The proportions therefore include expression/selection effects and should not be treated as a pure GEF-specific kinetic experiment with identical background states. [S1, Discussion]

**M.** Preserve construct identity. A model need not make every nominal “RhoGEF stimulus” behave identically. In particular, PRG's ability to favor protrusion at high expression does not imply the same robust behavior for LARG.

## 12.2 HeLa: the same response classes, different representative geometry

**R.** HeLa cells also show optoPRG-induced protrusion and retraction. **P/F, final supplement.** The displayed bar is **N = 20**, with approximately **60% protruding, 20% retracting, and 20% mixed/no movement**. These are broad bar readings for this sample, not expression-conditioned probabilities. [S1, first Results section; S3, panels D/E]

**F, final supplement.** The retracting example begins as a spread oval cell with a broad stimulated margin; the lower boundary withdraws progressively over 10–30 minutes. The protruding example is already irregular and branched, with long thin processes; the activated side extends further as a thin/irregular protrusive region. It does not resemble a spherical suspended cell producing one large round bleb. [S3, panel D]

**M.** Treat HeLa as a qualitative cross-cell-line check of stimulus-aligned movement, not a requirement to reuse RPE1's exact outline, area changes, or reporter amplitudes. No separate HeLa reporter-kinetics dataset is extracted here.

---

# 13. Intrinsic reporter kinetics and what the authors' model does—and does not—supply

## 13.1 RBD is a kinetic filter, not an instantaneous activity readout

**R.** Appendix 1 explicitly models the reporter's own binding/unbinding dynamics. In normalized form the observation equation is:

`db/dt = k_off,reporter × (r − b)`

where r is the accessible active-RhoA signal and b the normalized biosensor response. The authors estimate the reporter rate using PH recruitment, then estimate the RhoA deactivation parameter from low-optoPRG pulse responses. [S1, Appendix 1; S5, Appendix 1]

**R/P/F.** The reported central reporter value is **k_off ≈ 0.08 s⁻¹**. The final appendix plot shows individual fitted values approximately **0.04–0.18 s⁻¹**, from **three cells with repeated measurements**, and a central marker near 0.08. The printed uncertainty in the paper is **±0.4 s⁻¹**, which is inconsistent with the much narrower plotted spread. **This report does not silently change that to ±0.04.** Treat the uncertainty as a source-quality flag, not a reliable parameter prior. [S4; S5, Appendix 1; S2, p. 14]

**R.** The RhoA deactivation parameter is reported as **k₂ = 0.014 ± 0.003 s⁻¹**, fitted from three low-expression cells. **D.** The corresponding nominal first-order timescales are **1/0.08 = 12.5 s** for reporter relaxation and **1/0.014 ≈ 71 s** for the RhoA deactivation parameter. These are inverse-rate calculations, not observed time-to-peak values. [S1, Appendix 1 captions; S5, Appendix 1]

**M.** Use these as optional kinetic priors with their provenance and sample-size limits, not as exact mechanical CPM parameters. Do not assign the RBD off-rate to Lifeact, MRLC, Cdc42, or Rac1 reporters. The optogenetic actuator also has its own kinetics, distinct from the reporter.

## 13.2 Useful mechanistic structure, but not an independent morphology dataset

**R.** The model tracks free active RhoA, its PRG-bound complex, and Cdc42. PRG can activate both GTPases, while PH binding removes RhoA-GTP from the freely signaling pool without necessarily destroying GTP-bound RhoA. The model is local and non-spatial at the activation region; it does not calculate a cell outline, membrane curvature, substrate traction, or centroid velocity. [S1, Fig. 6; Appendix 1]

**R.** The model's phenotype indicator compares integrated RhoA and Cdc42 responses during the **first two minutes**, with an arbitrarily chosen relative weighting and a mixed-response band. Equal Cdc42 and RhoA deactivation rates are assumed for simplification. [S1, Fig. 6D/E; Appendix 1]

**M.** These are modeling choices, not newly measured experimental constraints. A CPM should be screened against the measured phenotype, reporter traces, and geometry—not awarded a pass merely because it duplicates the authors' chosen indicator or threshold. The experimental data permit mechanisms beyond the minimal model, which does not fit every observed RhoA trajectory. [S1, Fig. 6—figure supplement 1 discussion]

---

# 14. Compact CPM benchmark set

The following **M** criteria summarize what the experiments can discriminate. Their biological basis is the experiment-specific evidence above; they are proposed screening rules rather than author-defined pass/fail cutoffs.

| Condition | Shape/movement target | Local reporter target | Important failure to avoid |
|---|---|---|---|
| Low basal optoPRG | Stimulated edge retracts; occupancy gains elsewhere | Early RBD/Lifeact rise; accumulating MRLC | Whole-cell collapse counted as coordinated retraction |
| High basal optoPRG | Irregular/ruffled stimulated protrusion; withdrawal elsewhere | Early Cdc42; depressed MRLC; delayed Rac1 and RBD | Smooth bleb substituted for the ruffled front |
| Mixed response | Active ruffling but little net translation | Recruitment can occur despite weak net motion | Mixed class replaced by an inert cell |
| Opposite-pole re-stimulation | Same local response sign at the new pole | Response follows the new stimulus | Permanently fixed front/back identities |
| Extra non-recruitable PRG | More protrusion despite low recruitable actuator | Baseline-state dependence | Classification solely by peak recruited amount |
| Isolated PH | No reproducible directed protrusion | Local RBD suppression allowed | RBD inhibition automatically creates a front |
| PH-binding mutant | Strong shift toward retraction | Loss of acute sequestration function | High baseline expression always forces protrusion |
| PAK inhibitor | Protrusion lost; some retraction retained | Post-drug GTPase amplitudes not measured here | Every reporter forced to zero |
| 15-s versus 30-s pulses | Clear switching possible in a minority of intermediate cells | Input waveform must be temporally resolved | Universal switching across the population |
| HeLa validation | Same response classes with different starting shapes | No independent kinetic targets here | Copying RPE1 or HL-60 geometry by assumption |

## 14.1 Suggested observables

**M. Local area and translation.** Retain A_ROI(t), area outside the ROI, whole-cell projected area, and signed centroid displacement along the stimulus axis as separate outputs. Only the first two have the directly extracted numerical time courses supplied here. For new raw-data analysis, also calculate cell-body displacement, path length, reversals, and persistence; do not populate those fields with the article's area-based values.

**M. Shape.** Record projected area, aspect ratio, perimeter, circularity, stimulated-edge normal motion, protrusion width, branch/filopodium count, and coarse edge roughness. Existing images justify qualitative irregular-versus-smooth distinctions but not universal numerical cutoffs for these metrics. Roughness measures will depend on optical/lattice resolution.

**M. Spatial signals.** Compare mean reporter values within the stimulated basal region, local edge-associated signal, and spatial overlap with the moving boundary. Keep signal on the opposite side as a separate variable rather than inferring it from a local-to-rest ratio. Domain width, beyond-ROI spread, and perimeter coverage remain unquantified in this extraction.

**M. Timing.** Keep actuator recruitment, reporter response, and visible morphology distinct. Prioritize the ordering of early Cdc42, actomyosin divergence, shape change, and delayed Rac1; avoid sub-frame latency claims. Distinguish early transient RBD suppression from late RBD accumulation in protruding cells.

## 14.2 How to use the numerical ledger

The accompanying JSON records each numerical target's metric, time, source, version, and evidence type. **Visual lower/upper values are read-off bounds, not biological acceptance bands.** Fit broad experimental features first; use the approximate amplitudes as soft targets until numerical source files are retrieved. Do not treat successive points on the same mean trace as independent experimental replicates.

Keep at least one mechanistically discriminating perturbation—such as isolated PH, the binding mutant, or IPA3—out of parameter fitting for validation. This is a suggested modeling workflow, not a procedure performed in this paper. When comparing to both this paper and the de Belly reference, retain cell type, adhesion state, optical geometry, stimulation mask, reporter identity, and baseline-expression state as experiment-specific metadata.

---

# 15. Source-quality flags and quantities deliberately not invented

**Response percentages.** The main prose gives approximately 35% retracting, 40% protruding, 10% nonresponding, and 25% mixed for optoPRG; these sum to 110%. Figure 1D also combines mixed/nonmoving responses. Do not use those prose values as an exhaustive probability vector or silently renormalize them. The concentration-stratified 95%/85% comparisons and directly labeled cohort sizes are separate observations. [S1, first Results; S2, Fig. 1D]

**Movie construct labels.** Several final movie captions refer to optoLARG even when their linked Results concern optoPRG, including the protrusion and reporter experiments. These captions should not be used alone to relabel the experimental construct. The report uses the experiment-specific main text and figure captions and makes no claim to movie-frame inspection. [S1, Videos 3, 5–7 and linked Results]

**Reporter-rate uncertainty.** The printed “0.08 ± 0.4 s⁻¹” does not match the plotted spread; the central value and the problematic uncertainty are kept separate. [S4; S5]

**Area-based quantities.** Figure 2A includes an area-normalization label alongside area units, and Figure 2E calls an area-rate measure “speed.” This extraction preserves normalized ROI area as dimensionless and distinguishes the printed µm²/min measure from linear translation speed. [S2, Fig. 2A/E; S1, Data processing]

**Model transcription.** The displayed second reaction in the final Appendix uses a GDP label, whereas the surrounding text and differential equation describe PH binding to active RhoA-GTP. The mechanistic interpretation here follows that explicitly stated active-RhoA binding process; a coding implementation should not blindly copy the isolated reaction label. [S1, Appendix 1, reaction list and following ODEs]

**Not available as verified numerical targets in this extraction:** centroid velocity or persistence distributions; exact turn angles; complete spatial reporter profiles; cap width or angular coverage; filopodial dimensions/counts; bleb counts, sizes, lifetimes, or membrane-coating patterns; traction/contractile force; membrane tension; adhesion strength; absolute GTPase concentrations; and frame-resolved movie-derived delays. The paper may enable additional extraction through its raw images, but that reanalysis was not performed here.

## Overall experimental signature

The strongest combined benchmark is not simply “Rho causes contraction.” It is **the same locally recruited PRG DH–PH module producing opposite boundary motions depending on baseline state and recruitment dynamics**: low-state cells retract with an early, strong accessible-RhoA/actomyosin response; high-state cells protrude through an irregular edge with early Cdc42, reduced local myosin, and delayed Rac1/RBD increases. Sequestration alone can reduce the RhoA reporter without generating a front, and the PH-binding mutant and PAK inhibition disrupt the protruding outcome in different ways. Those joint morphology–kinetics–perturbation distinctions are the most useful constraints for a CPM. [S1, Figs. 2–6; S2, Figs. 3–6; S3]