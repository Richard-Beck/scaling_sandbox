Below I use the title rather than the misattributed author name where the two disagree. Thus **study 1 is Nanda et al. 2023, not Mahlandt et al.**, and **study 4 is Inaba, Miao & Nakata 2021, not Katsura et al.** All four title-matched papers were retained because I could access the complete main text and main figures. For supplementary material, I distinguish material I could actually inspect from material that was merely linked.

**Provenance tags throughout:** **R** = value/detail explicitly reported in main or supplementary text/caption; **P** = value/text printed directly on a displayed figure; **F** = my approximate read from a plotted curve, bar, scatter, or image. `~` always denotes an **F** estimate, never an exact value.

# 1. Nanda et al., 2023 — Rac→Rho crosstalk through ARHGEF11/12

**Correct citation.** Suchet Nanda, Abram Calderon, Arya Sachan, Thanh-Thuy Duong, Johannes Koch, Xiaoyi Xin, Djamschid Solouk-Stahlberg, Yao-Wen Wu, Perihan Nalbant & Leif Dehmelt. “Rho GTPase activity crosstalk mediated by Arhgef11 and Arhgef12 coordinates cell protrusion-retraction cycles.” *Nature Communications* **14**, 8356 (2023). DOI 10.1038/s41467-023-43875-y. The paper specified by title is therefore **Nanda et al.**, not Mahlandt et al. ([Nature][1])

## Experimental system and exact comparison

The study uses **Neuro-2a, A431, HeLa, NIH3T3 and U2OS** cells. There is **no ploidy manipulation, polyploid state, or clone-matched ploidy comparison**: those fields are N/A. The causal experiments instead manipulate Rho-family signaling acutely by chemically induced membrane recruitment of constitutively active Rac1/Cdc42/RhoA, by PA-Rac1 photoactivation, by ARHGEF11/12 overexpression, and by siRNA depletion. These are perturbation/control comparisons in the same cell-line background, but **not isogenic clone comparisons** in the genetic sense. ([PubMed Central (PMC)][2])

Culture/imaging conditions matter considerably. Neuro-2a were grown in MEM + 10% FBS; A431 in DMEM + 10% FBS; HeLa, NIH3T3 and U2OS in DMEM + 10% FBS. Neuro-2a and HeLa were imaged on uncoated glass, NIH3T3 on poly-L-lysine, U2OS on 10 µg/ml collagen I, and A431 on **10 µg/ml fibronectin**. For peripheral GEF/sensor enrichment, cells were replated on fibronectin-coated glass prepared at 10 µg/ml for 16 h at 4 °C. ([PubMed Central (PMC)][2])

### Figure 1 — establishing rapid reversible Rho-family perturbation

**1A.** Schematic, not quantitative: constitutively active Rac1-Q61L, Cdc42-Q61L or RhoA-Q63L lacking their CAAX motif is recruited from cytoplasm to membrane by SLF′-TMP-induced FKBP′/eDHFR heterodimerization; trimethoprim reverses recruitment. **1B** shows Neuro-2a TIRF images with mCherry-actin. Images are shown **30 s before, 24 min during and 24 min after** activation; scale 10 µm, 0.26 µm/pixel. At least three independent repetitions and ≥40 cells/condition were used. ([PubMed Central (PMC)][2])

The printed response frequencies are approximately **86% Rac1-induced protrusion, 74% Cdc42-induced protrusion, and 84% RhoA-induced retraction [P]**; phenotypic reversal after removal of activation is approximately **94%, 92%, and 75%, respectively [P]**. The exact per-condition cell counts printed in the panels are approximately n=46, 77 and 63 for initial activation and n=40, 61 and 51 for reversibility [P]. The important distributional point is that the response is common but not universal; especially RhoA reversal has a visibly lower success fraction than Rac/Cdc42 reversal.

**Qualification:** these are binary phenotype scores rather than migration speeds or displacement measurements.

### Figure 2 — pairwise Rho-family crosstalk

**2A** defines perturbation + simultaneous TIRF sensor measurement. Sensors are translocation reporters: p67^phox^ GBD for the Rac subgroup, WASP GBD for the Cdc42 subgroup and Rhotekin GBD for the Rho subgroup; critically, these reporters do **not** distinguish individual close paralogues such as RhoA versus RhoB/C. ([PubMed Central (PMC)][2])

**2B–C** show raw Neuro-2a images and a representative Rac1→Rho experiment. The authors explicitly correct activity signals with a cytosolic volume reporter because membrane/cell-volume changes themselves alter TIRF intensity. This correction is consequential: in the example, the Rho reporter rises while the volume-control signal falls, so an uncorrected trace underestimates the response. ([PubMed Central (PMC)][2])

**2D–E** show time courses and pre-versus-5-min changes. Selected sample sizes are **n=48 RhoA→Cdc42, n=40 RhoA→Rac, n=51 Rac1→Rho, n=28 Rac1→Cdc42 [P/Figure]**. Visually, the 5-min responses are roughly **+10% for RhoA→Cdc42, −5% for RhoA→Rac, +7% for Rac1→Rho and +7% for Rac1→Cdc42 [F]**. All these estimates have substantial cell-to-cell spread; the interesting point is sign consistency at the mean, not separation into non-overlapping populations. Statistical comparisons are two-sided Student t-tests, with the displayed crosstalk effects reaching ****p<0.0001 where marked. ([PubMed Central (PMC)][2])

**2F** is the resulting influence diagram. It should not be read as a kinetic model; it encodes statistically detectable pairwise responses at the chosen 5-min readout.

### Figure 3 — Rac→Rho is reproduced across cell lines

This is the first especially useful figure for model benchmarking.

**3A** switches to mCerulean-PA-Rac1, illuminated under TIRF at 445 nm. Neutral-density attenuation was deliberately strong (1,000–10,000× plus ~30% built-in filter) to avoid saturating PA-Rac1; activation light remains on during the activation intervals except while fluorescence images are acquired. ([PubMed Central (PMC)][2])

**3B/F.** Rho response is measured in NIH3T3, HeLa, N2a and U2OS as the difference between **five measurements immediately before and five after activation**, with each measurement representing a 25-s frame **[R]**. Approximate figure means are NIH3T3 **~+9%**, HeLa **~+24%**, N2a **~+18%**, U2OS **~+7% [F]**; matched controls cluster close to zero. Sample sizes are approximately n=43, 55, 56 and 49 [P]. The major result is therefore not a fixed gain but **a positive response with a roughly threefold range in magnitude across cell backgrounds**. ([PubMed Central (PMC)][2])

**3C–E. A431 is strongly heterogeneous.** Among **49 cells from three independent experiments [R]**, **67% showed a positive reversible Rho response**: 29% had sustained activation and 38% a transient pulse; **25% showed no response and 8% a negative response [R]**. This is unusually useful information because the mean masks genuinely different dynamic regimes. ([PubMed Central (PMC)][2])

**3D** is the cell-by-cell heat map; visually the responses span roughly −25% to +30% [F]. **3E** averages them, which consequently smooths the transient responders and nonresponders together. **3G** compares peripheral versus central attachment regions; there is no strong spatial bias in this global PA-Rac1 response.

### Figure 4 — spontaneous Rac/Rho dynamics during protrusion–retraction cycles

A431 cells are now observed **without imposed Rac activation**. This figure is highly relevant to phenomenological migration constraints.

**4A–C** use improved tandem Rac or Rho sensors and a modified ADAPT edge-tracking pipeline. Rac dataset: **n=23 cells**, Rho dataset **n=21**, each from ≥3 independent experiments [R]. Cell-edge velocity is mapped around the circumference through time. Red is local protrusion and blue local retraction. ([PubMed Central (PMC)][2])

The analysis explicitly classifies **protrusion as edge velocity > +0.075 µm/min and retraction as < −0.075 µm/min [R]**. Sensor enrichment is the local reporter intensity divided/normalized against the cell-wide time-averaged signal, rather than a biochemical concentration.

**4D/F.** A representative local Rac burst accompanies positive edge velocity; visually, individual Rac enrichment peaks can reach several-fold above the normalized baseline while edge speed reaches approximately **0.4–0.5 µm/min [F]**. These are local excursions, not mean migration speed.

**4G/H.** Cross-correlations place Rac essentially at the protrusive phase and Rho at the retractile phase, with little detectable lag of the maxima/minima. Rough visual enrichment maxima are about **+40% Rac in protrusions versus near-baseline/small enrichment in retractions [F]**, while Rho is approximately **+35–40% during retraction versus ~15–20% during protrusion [F]**. The distributions overlap; these are preferential enrichments, not exclusive states. ([PubMed Central (PMC)][2])

### Figure 5 — ARHGEF11/12 temporally bridge protrusion and retraction

**5A** is the proposed Rac→Lbc-GEF→Rho mechanism.

**5B.** A focused A431 screen compares seven Lbc-family GEFs. Printed n values are **ARHGEF1 12, ARHGEF2 13, ARHGEF11 14, ARHGEF12 12, AKAP13 13, ARHGEF18 10, ARHGEF28 16 [P]**, three independent experiments. ARHGEF11 and 12 are the conspicuous members enriched around both protrusion and retraction, whereas several others remain comparatively flat. ([ResearchGate][3])

**5C–D.** Higher-sensitivity measurements use **n=37 ARHGEF11 and n=22 ARHGEF12 [P]**. Both show broad enrichment over the protrusion–retraction sequence. Approximate maximal enrichment is of order **~35–45% above the control-normalized baseline [F]**.

**5E** is particularly useful temporally. Relative to protrusion onset, active Rac peaks first, ARHGEF11 at roughly **160 s**, ARHGEF12 around **210 s**, and Rho considerably later, beyond **~470 s [R/figure-defined maxima]**. Thus this is not simply synchronous Rac/GEF/Rho activation.

### Figure 6 — causal recruitment of ARHGEF11/12 by Rac

**6B–F.** A431 cells coexpress PA-Rac1 and full-length/mutant ARHGEF11/12. Three independent experiments are used. Recruitment is quantified in an **early 25-s window during photoactivation** and a **late window 1 min afterwards [R]**. Wild-type ARHGEF11 and ARHGEF12 both rapidly recruit after Rac activation. ([PubMed Central (PMC)][2])

The PH-domain mutants do not simply eliminate recruitment: visually, ARHGEF11 WT and PH mutant both show small early gains (~3–4% [F]), whereas ARHGEF12 WT recruitment is stronger (~6% versus ~2% mutant [F]). Thus the mechanism is not reducible to “PH domain on/off.”

**6E–F** additionally test ARHGEF11's F-actin-binding region. The ΔFAB mutant particularly compromises the later response [P/F], supporting more than one route to membrane/cortical association.

**6G–H.** siRNA against either GEF reduces Rac-induced Rho activation. Approximate 25-s Rho response means are **control ~14%, siARHGEF11 ~8%, siARHGEF12 ~9% [F]**. Knockdown efficiencies reached roughly **83±15% for ARHGEF11 and 75±11% for ARHGEF12 [R]**. ANOVA with Holm–Sidak post-tests is used here. ([PubMed Central (PMC)][2])

### Figure 7 — morphodynamic and whole-cell migration consequences

This is the paper's most directly migration-comparable figure.

**7A–B, overexpression.** A431 cells imaged at **1.5 frames/min**. There are three independent experiments and **>26 cells per condition**. P–R cycle duration is defined from edge-velocity dynamics, from onset of a protrusive episode to onset of the following retractile episode. Approximate central values are **control ~20 min, ARHGEF11 ~15–16 min, ARHGEF12 ~17 min [F]**. The distributions strongly overlap, but both overexpressors shift toward shorter cycles. ([PubMed Central (PMC)][2])

**7C–D, knockdown.** Three independent experiments and **>105 cells/condition**. Approximate means are **control ~14–15 min, siARHGEF11 ~21–22 min, siARHGEF12 ~17–18 min [F]**. The ARHGEF11-depleted population has a conspicuous long tail, with some cycles extending toward **~80 min [F]**; ARHGEF12's tail is shorter, about ~40 min. The paper explicitly warns that the control means in 7B and 7D should not be compared directly because the protocols differed (single versus dual transfection). ([PubMed Central (PMC)][2])

**7E defines the migration metrics exactly.**

* **Distance** = cumulative trajectory/path length.
* **Displacement** = Euclidean straight-line distance from initial to final location.
* **Directionality** = displacement / distance.

This is **not** persistence time, MSD exponent, instantaneous speed or mean velocity. ([PubMed Central (PMC)][2])

**7F–H.** A431 nuclei are tracked at **1 frame/min for 4 h**, three independent experiments, **>491 cells per condition**. The migration assay was unusually well specified: 96 h after siRNA, cells were seeded on **10 µg/ml fibronectin at ~143 cells/mm²**; 6 h later they were stained with SPY-650-DNA for 1.5 h, washed and rested for 1 h before imaging. Cells leaving the FOV were excluded. ([PubMed Central (PMC)][2])

Approximate distributions/means:

| Metric               |         Control |      siARHGEF11 |      siARHGEF12 |
| -------------------- | --------------: | --------------: | --------------: |
| 4-h **distance**     | ~175–180 µm [F] | ~145–150 µm [F] | ~145–150 µm [F] |
| 4-h **displacement** |   ~13–15 µm [F] |   ~13–15 µm [F] |   ~13–15 µm [F] |
| **directionality**   | ~0.075–0.08 [F] | ~0.10–0.105 [F] | ~0.09–0.095 [F] |

Distance distributions overlap substantially and extend to roughly 350–400 µm [F]. Displacement distributions overlap almost completely. Directionality has a large right tail in every condition, so the mean shift is modest relative to single-cell variability. The key decoupling is therefore **less accumulated movement without appreciably less endpoint displacement**, producing a higher displacement/path-length ratio. ANOVA + Holm–Sidak is used. ([PubMed Central (PMC)][2])

### Figure 8

Schematic model only: more ARHGEF11/12 → greater Rac/Rho crosstalk → shorter protrusion–retraction cycles; depletion does the reverse. No new quantitative data.

### Supplementary material

The PMC record exposes a **2.2-MB Supplementary Information PDF, seven supplementary movies, reporting summary and source-data spreadsheet**. I could inspect the supplement inventory and main-text references to individual supplementary results, but the supplementary PDF/movie frames themselves were not rendered through the current interface. I therefore **do not report F-type graphical values from supplementary figures**. Main-text-described supplement findings include dose-dependence of PA-Rac1→Rho response, sensor/volume controls, individual-cell versions of the enrichment analyses and independent siRNA migration results. ([PubMed Central (PMC)][2])

## Internal/provenance issues

The strongest qualifications are that the “Rho,” “Rac,” and “Cdc42” sensors are subgroup rather than isoform specific; TIRF signals require explicit morphology/volume correction; the Figure 7 control groups are not interchangeable; and the migration assay's “directionality” is simply endpoint displacement/path length and should **not** be translated into persistence time. The strong A431 heterogeneity in Fig. 3 also means a single mean response hides sustained, transient, null and even inverse Rac→Rho responses.

---

# 2. Wu et al., 2009 — photoactivatable Rac and directed motility

**Citation.** Yi I. Wu, Daniel Frey, Oana I. Lungu, Angelika Jaehrig, Ilme Schlichting, Brian Kuhlman & Klaus M. Hahn. “A genetically encoded photoactivatable Rac controls the motility of living cells.” *Nature* **461**, 104–108 (2009). DOI 10.1038/nature08241. ([Nature][4])

## Experimental system and exact comparison

The main motility experiments use **mouse embryonic fibroblasts (MEFs)** stably expressing PA-Rac1; HeLa cells are used extensively for morphology/RhoA measurements, with HEK293 and MET cells appearing in supplementary experiments. There is **no ploidy/polyploid manipulation and no ploidy clone structure**. The strongest causal comparison is acute light-on versus dark/local-light location within the same cell, backed by light-insensitive C450A/C450M and dominant-negative T17N variants. ([PubMed Central (PMC)][5])

Live-imaging cells were plated on coverslips coated with **5 µg/ml fibronectin**, in phenol-red-free Ham's F-12K containing **2% FBS [R]**. For the principal MEF localization experiments, cells were additionally described as **serum-starved** to suppress spontaneous activity. Imaging used an Olympus FV1000 confocal, typically a 60×/1.42-NA oil objective. Photoactivation was at 458 or 473 nm; spot diameter was typically 10 or 20 µm depending on experiment. ([Gale][6])

### Figure 1 — engineering and spatial localization of PA-Rac1

**1A–B.** Structural/biochemical design. LOV2 is fused to Rac1, with the final PA-Rac1 carrying Q61L/E91H/N92H and linkage around LOV residue 546 to Rac1 residue 4. C450A is a light-insensitive control and I539E mimics the lit state. These panels establish caging; they do not measure migration. ([PubMed Central (PMC)][5])

**1C.** Whole-cell 458-nm activation in HeLa rapidly generates lamellipodial protrusions and ruffling. A separate supplementary kymograph analysis found maximum protrusions approximately **4× longer** for PA-Rac1 than LOV-alone or C450A controls **[R, reported from Supplementary Fig. S4 in main text]**. ([Gale][6])

**1D.** Serum-starved PA-Rac1 MEFs receive a **20-µm circular illumination every 60 s** at one edge. The outlined cell at baseline versus 10 min shows strongly localized extension near the spot and little movement elsewhere; the kymograph shows repeated ruffles followed by protrusion. This is local membrane displacement, not centroid speed. ([PubMed Central (PMC)][5])

Useful physical characterization reported alongside/supplementary to Fig. 1: PA-Rac1 lateral diffusion **D=0.55 µm²/s**, FRAP **t½=12.1 s**, PA-GFP tracking **t½=14.6 s**, and dark-state recovery **t½=43 s at room temperature [R]**. A simulation predicted at most **7.5% spillover activation** into an adjacent unilluminated 10-µm spot. These numbers set the approximate spatial/temporal resolution of the actuator. ([PubMed Central (PMC)][5])

### Figure 2 — Rac is sufficient to impose cell polarity and motion

**2A.** After one local pulse, **all 50/50 PA-Rac1 MEFs** protruded near the illuminated site and retracted at the opposite side **[R]**. Conversely, **all 25/25 PA-Rac1-T17N cells** retracted locally and protruded elsewhere. The panel is a protrusion/retraction map, not a trajectory plot. This is unusually clean categorical causality. ([PubMed Central (PMC)][5])

**2B.** Repeated edge activation at **2-min intervals** drives directional movement. Mean movement is **0.8 µm per pulse, n=6 [R]**. Interpreted literally, that is pulse-locked centroid/forward movement; it should not be reported as a continuous speed without also stating the 2-min stimulation interval. The implied pulse-averaged rate would be 0.4 µm/min, but that conversion is ours rather than an author-defined metric and therefore should not replace the reported 0.8 µm/pulse. ([PubMed Central (PMC)][5])

**2C.** Local PA-Rac1 activation is repeated with **1 µM ML-7, 1 µM blebbistatin, or 10 µM Y-27632**. Each condition has **n>25**, values displayed as means ±95% CI. ML-7 and blebbistatin leave local protrusion comparatively intact, whereas Y-27632 reduces it visibly. Exact y-axis values are not recoverable at adequate resolution from the displayed figure, so I do not turn the bars into false precision. ([PubMed Central (PMC)][5])

**2D. Metric definition is critical.** The plotted “polarity index” is the **cosine of the angle θ** between (i) the vector from illuminated spot to initial cell centroid and (ii) the vector from initial to final centroid. Hence +1 means movement in the optically specified direction, 0 orthogonal movement, −1 the reverse. It is **not** displacement/path-length directness, migration persistence, or velocity. n>25, means ±95% CI. ML-7 and blebbistatin markedly reduce this directional alignment despite preserving protrusion, whereas the ROCK inhibitor has a smaller/different effect. ([PubMed Central (PMC)][5])

A caption error is worth retaining: the Figure 2 caption says “throughout figure 3 irradiation at 458 nm,” plainly referring to the displayed Figure 2 experiment. That is a typographical internal inconsistency rather than a biological result. ([PubMed Central (PMC)][5])

### Figure 3 — local Rac inhibits RhoA

**3A.** HeLa cells express a RhoA FRET biosensor plus PA-Rac1 or C450M. A **single 473-nm pulse in a 10-µm circle** is applied. n=3 cells per construct. In active PA-Rac1 cells, the irradiated region remains at reduced RhoA activity after the initial photobleaching and an adjacent unbleached region gradually drops to a similar level. In C450M controls, the bleached region recovers and the neighboring region remains stable. Thus the RhoA suppression propagates beyond the optical ROI. ([ResearchGate][7])

**3B** contrasts naturally occurring versus PA-Rac1-induced pseudopods. The induced protrusion uses a **20-µm 473-nm circle**.

**3C.** Quantification is **percent increase in FRET/CFP ratio within 1 µm of the cell edge relative to the mean of a flat region to the left of the line scan**, with **18 line scans from six cells per bar**, means ±95% CI. This is a spatial edge-enrichment metric, not total RhoA activity. ([ResearchGate][7])

### Figure 4 — structure, not a migration assay

**4A–B** show dark-state PA-Rac1 crystal structure and the LOV–Rac interface. **4C** uses the key Rac Trp interaction to improve LOV inhibition of Cdc42: introducing **F56W** enables caging; a lit-state mutation restores effector binding. This establishes generalizability of the actuator design but adds no quantitative migration phenotype. ([ResearchGate][8])

### Supplementary figures and movies

This study has unusually extensive supplementation: **Supplementary Figs S1–S14, Tables S1–S5 and Movies S1–S16**. The Nature supplement index and supplementary methods/text were accessible; I did **not** visually quantify the graphs in the supplementary PDF itself, so no supplementary `F` estimates are reported. ([Nature][4])

Several supplement results are nonetheless explicitly reported in accessible text and are therefore legitimate **R** values: S4 gives the ~4× protrusion-length difference; S8 gives **D=0.55 µm²/s, FRAP t½=12.1 s, PA-GFP t½=14.6 s, dark recovery 43 s, ≤7.5% adjacent-spot activation**; the supplementary methods specify the 5 µg/ml fibronectin/2% FBS imaging condition and define protrusion length as the post-irradiation maximal edge position minus the pre-irradiation mean edge position along a kymograph. ([Gale][6])

The movie inventory is informative but I did not perform frame-by-frame quantitative analysis: S9 shows PA-Rac1-induced polarization, S10–11 polarity reversal with T17N, S12 directional MEF movement, S13 directional HEK293 movement, S14 propagating RhoA downregulation, and S15 ruffle translocation. ([Nature][4])

## Internal/provenance issues

This is a strong causal migration paper, but many commonly quoted “migration” results are actually **local protrusion or angular alignment**, not trajectory statistics. The 0.8 µm/pulse datum is based on only **six MEFs**. The MEF motility behavior also does not generalize trivially to HeLa: HeLa generated local protrusions but generally failed to generate the coordinated opposite-side retraction needed for whole-cell motion. ([PubMed Central (PMC)][5])

---

# 3. de Seze et al., 2025 — one RhoA GEF can drive protrusion or retraction

**Citation.** Jean de Seze, Maud Bongaerts, Benoit Boulevard & Mathieu Coppey. “Optogenetic control of a GEF of RhoA uncovers a signaling switch from retraction to protrusion.” *eLife* **12**, RP93180 (Version of Record, 27 May 2025). DOI 10.7554/eLife.93180.4; umbrella article DOI 10.7554/eLife.93180. ([eLife][9])

## Experimental system and exact comparison

The main system is **hTERT-RPE1**, with HeLa used as an independent cell background. Again there is **no ploidy or polyploid comparison** and no matched clone structure. Cells are transiently transfected with iLID/SspB constructs recruiting the DH–PH region of **LARG/ARHGEF12, GEF-H1/ARHGEF2, or PDZ-RhoGEF/ARHGEF11 (PRG)**. The striking comparison is therefore **low versus high transient expression within one cell line**, followed by identical local light recruitment. This is mechanistically causal with respect to membrane recruitment but **not a clean randomized causal manipulation of baseline GEF abundance** because expression level covaries with the pre-activation cell state. ([PubMed Central (PMC)][10])

Cells were deliberately plated at **low density** so individual transfected cells could be found across a 25-mm coverslip. Experiments were performed ≥30 h after DNA transfection, at 37 °C/5% CO₂, using 60×/1.45-NA optics, DIC + fluorescence/TIRF and patterned **440±10-nm** illumination from a DMD. No extracellular chemotactic gradient is used. The paper does not specify a matrix coating for these assays. ([PubMed Central (PMC)][10])

### Figure 1 — the same local recruitment yields opposite cell-polarity responses

**1A–C** define iLID/SspB membrane recruitment and the three RhoA GEFs. iLID/SspB dissociates in darkness with an approximately **20-s off rate/time scale [R]**. Cells are transiently transfected ≥30 h, then observed for **30–60 min** while pulse duration, intensity and frequency are varied. ([PubMed Central (PMC)][10])

**1D.** LARG is essentially retractile. GEF-H1 often produces round/nonresponsive cells, with only **~40% showing a distinct response [R]**; those responders are mostly retractile but include some protrusions. PRG is the robust bidirectional case. The text gives **~35% retracting, ~40% protruding, ~10% no clear response and ~25% mixed [R]**. These numbers are internally problematic: they sum to **~110%**, not ~100%. That is either rough rounding or a textual/categorization error and should not be silently normalized away. ([PubMed Central (PMC)][10])

**1E–F.** The primary quantitative morphology readout is area change **inside the activated region**, normalized to the mean initial area. Protruders increase area while retractors decrease it; unilluminated controls show far smaller motion. The time course extends about 15 min in the displayed plots. This is **local projected cell area**, not path length or cell speed. ([eLife][11])

**1G–I.** Cells are activated on one pole, allowed to recover for **1 h**, then stimulated on the opposite pole. Many reproduce the same phenotype at the newly chosen side, supporting light-defined polarity rather than pre-existing front/back identity.

**Figure 1—supplement 1** was accessible through eLife. It shows mixed cells that first ruffle/protrude and later contract, opposite-side area changes demonstrating whole-cell polarization, and reproduces both protruding and retracting responses in HeLa. ([eLife][11])

One significant caption inconsistency: the eLife video legends repeatedly refer to examples as being activated with **“optoLARG”** even where the surrounding Results discussion concerns the optoPRG bidirectional phenotype. This makes the movie labels unsafe for tool-identity inference without cross-checking the main text. ([PubMed Central (PMC)][10])

### Figure 2 — baseline optoPRG abundance predicts phenotype

**2A.** The y-axis is normalized local area after **5 min**; x-axis is initial cytosolic optoPRG fluorescence in arbitrary units. Below roughly **40 a.u., ~95% of cells retract; above it, ~85% protrude [R]**. This is one of the clearest threshold-like relationships across the four papers. However, the two clouds are not defined by an absolute protein concentration; 40 a.u. is microscope/expression-system specific. ([PubMed Central (PMC)][10])

**2B** gives representative low-expression retractor, intermediate/mixed cell and high-expression protruder.

**2C** plots absolute recruitment. Retractors can achieve high absolute recruited fluorescence, arguing against “more recruited GEF at the illuminated membrane” alone explaining the switch.

**2D.** Adding a separately expressed, **non-recruitable PRG DH–PH domain** changes phenotype. Among low optoPRG cells, the protruding fraction rises from **0% to 45% [R]**, demonstrating that pre-existing global PRG abundance/state matters independently of optically recruited dose. ([PubMed Central (PMC)][10])

**2E.** “Membrane mean absolute displacement” before stimulation is higher in eventual protruders, Mann–Whitney **p<0.0001 [R]**. It is a local membrane-motility measure, not centroid displacement.

**Figure 2—supplement 1.** High-expressing/eventual-protruding cells are about **1.5× larger [R]** before stimulation and have higher basal RhoA-GTP staining. This is a major confound if one wants to interpret expression level as the only changing variable. ([PubMed Central (PMC)][10])

### Figure 3 — biochemical divergence precedes visible morphology

**3A–C, Lifeact.** In retractors, actin signal rises immediately after optoPRG recruitment. In protruders, it initially falls, then later rises as protrusion develops. Pairwise t=0 versus **60 s** comparisons are shown.

**3D–F, MRLC/myosin.** Retractors recruit/increase myosin after activation; protruders show an early decline that remains below baseline longer. Pairwise endpoint is **90 s**.

**3G–I, RhoA reporter.** Retractors show strong, immediate RhoA-GTP reporter recruitment. Protruders show a much smaller, delayed signal and by 60 s no significant net increase. The plotted group means have clearly overlapping scatter but opposite early trends. Tests are described as Wilcoxon for within-cell before/after and independent t-test between phenotypes. ([eLife][11])

An important methods inconsistency exists here: the Methods call the paired-timepoint test a **“Wilcoxon rank sum”** test. Rank-sum is conventionally an unpaired test; for paired observations one would normally expect Wilcoxon signed-rank. The figure legend simply says “Wilcoxon test.” The exact implementation therefore remains ambiguous from the prose. ([PubMed Central (PMC)][10])

### Figure 4 — PRG's PH domain can suppress available RhoA-GTP

**4A.** Three exemplar RhoA reporter traces show qualitatively different dynamics with increasing optoPRG expression: low expression produces straightforward RhoA activation; intermediate/high expression produces an early dip followed by recovery/activation. Thus the population split is already visible at signaling level.

**4B–C.** Recruitment of the isolated PRG PH domain decreases the RhoA reporter. The corresponding supplement averages **26 cells [R]** and shows PH recruitment rising while RBD signal falls. ([PubMed][12])

**4E–G.** A PH mutant designed to disrupt RhoA-GTP binding largely converts high-expression cells back toward the retractile phenotype. A few very high-expressing cells still form ruffles, so the mutation is not an absolute “protrusion off” switch.

**4H–I.** Recruiting PH alone into cells globally overexpressing PRG does **not** generate a clear protrusion or retraction. Thus PH-dependent RhoA inhibition is **necessary but not sufficient** for the protrusive response. Variance in final normalized area differs at p<0.05 by Levene's test. ([PubMed Central (PMC)][10])

### Figure 5 — Cdc42, rather than Rac1, is the immediate protrusive branch

**5A–C.** A PAK-binding-domain biosensor that reads Rac/Cdc42 rises immediately in both phenotypes. After about **1 min**, signal continues increasing in protruders but plateaus in retractors.

**5D–F.** The more Cdc42-selective WASP-GBD reporter rises immediately in both groups but remains elevated preferentially in protruders.

**5G–I.** The Rac reporter (3×p67^phox^) rises later and predominantly in protruders. Thus the visually recoverable ordering is **Cdc42 first, Rac later**, rather than simultaneous Rac/Cdc42 activation. Pairwise endpoints are 60 s for PBD and **30 s** for the more specific Cdc42/Rac reporters. ([PubMed Central (PMC)][10])

**5J–K.** Cells are first phenotyped, allowed a 30-min interval, then exposed to **5 µM IPA-3**; a second stimulation begins **5 min after drug-containing medium replacement [R]**. The Methods say **20–30 cells** are selected for the drug experiment. No previously protruding cells can protrude again under IPA-3; some instead switch to retraction, whereas many original retractors retain retraction. ([PubMed Central (PMC)][10])

There is another textual typo: the Figure 5 caption calls this the “IP3 experiment” although the reagent is **IPA-3**, a PAK inhibitor.

### Figure 6 — minimal dynamical model and within-cell switching

**6A–D.** The minimal model has PRG activating RhoA, PRG PH binding/sequestering RhoA-GTP, and PRG activating Cdc42. Most parameters are fixed/fitted separately, with **G_eq_/K_b_** serving as the key cell-to-cell parameter associated with basal GEF level. The phenotype classifier is an artificial balance, γ, between **integrated RhoA and Cdc42 activity over the first 2 min**. Positive γ is assigned retraction; negative γ protrusion. The scale factor between RhoA and Cdc42 is chosen to reproduce the observed transition and the model includes an **arbitrarily defined gray/mixed zone [R]**. ([PubMed Central (PMC)][10])

**6E.** At intermediate G_eq_/K_b_ (~0.5), model behavior changes with **pulse frequency**, whereas changing pulse intensity alone cannot cross the phenotype boundary in this linear construction.

**6F.** Experimentally, cells in an intermediate fluorescence range of **14–35 a.u.** were screened. Only **3 of 90 cells [R]** could be driven repeatedly between retraction and protrusion by alternating **low-frequency/high-power pulses every 30 s** with **high-frequency/lower-power pulses every 15 s**, while keeping total integrated light comparable. That ~3.3% success rate is important: the model prediction is experimentally realizable, but only in a very small selected subset. ([PubMed Central (PMC)][10])

**6G–H** are conceptual summaries rather than new measurements.

**Figure 6—supplement 1** explicitly shows RhoA dynamics the minimal model does not capture and demonstrates in the model that pulse intensity alone cannot switch phenotypes. This is a useful negative result rather than just a cosmetic supplement. ([PubMed Central (PMC)][10])

### Metric definitions that should not be conflated

The image-analysis definition of **surface displacement** is area of overlap/intersection between the activation region and segmented cell. “Persistence” elsewhere in their image-processing code is defined as **distance from the initial point divided by the sum of absolute displacements**. That is a geometric directness-like ratio, **not persistence time** and not a directional autocorrelation decay. Biosensor intensities are locally background-subtracted, divided by signal in the nonactivated cell, then normalized to the prestimulation value. ([PubMed Central (PMC)][10])

### Supplementary figures/movies

Unlike the other papers, a meaningful fraction of the supplementary figure material was integrated and inspectable through eLife HTML: Fig. 1-supplement 1, Fig. 2-supplement 1, Fig. 4-supplement 1, Fig. 6-supplement 1, and three appendix fitting figures. The appendix fits use **N=3 cells** for several parameter-estimation steps. eLife also exposes **seven videos**, whose captions I inspected; I did not perform numerical frame extraction from the actual movies. ([PubMed Central (PMC)][10])

## Internal/provenance issues

The largest are: the response percentages sum to ~110%; several movie legends call the construct optoLARG where the Results imply optoPRG; “IP3” appears where IPA-3 is meant; the statistical prose ambiguously says Wilcoxon rank-sum for paired measurements; expression is in acquisition-specific arbitrary units; and baseline expression is entangled with **cell size, basal RhoA state and pre-activation membrane dynamics**. Consequently, Fig. 2 is extremely useful phenomenologically but should not be treated as a clean dose-response in which only PRG concentration changes.

---

# 4. Inaba, Miao & Nakata, 2021 — RhoA–calcium coupling

**Correct citation.** Hironori Inaba, Qianqian Miao & Takao Nakata. “Optogenetic control of small GTPases reveals RhoA mediates intracellular calcium signaling.” *Journal of Biological Chemistry* **296**, 100290 (2021). DOI 10.1016/j.jbc.2021.100290. This title is **not a Katsura et al. paper**. ([PubMed Central (PMC)][13])

## Experimental system and exact comparison

Cell lines are hTERT-RPE1, HeLa, MDCK and HEK293T. No ploidy manipulation, polyploid state or clone comparison is present. RPE1 were maintained in DMEM/F12 +10% FBS; the other lines in DMEM +10% FBS. Experimental comparisons are acute iLID-mediated recruitment of GEF domains for **RhoA, Rac1, Cdc42, Ras, Rap or Ral**, plus inactive/control switches and pharmacological/siRNA perturbations. These are causal optogenetic signaling tests, not migration assays in the usual trajectory sense. ([PubMed Central (PMC)][13])

For live imaging, medium was exchanged to Leibovitz L-15 or defined Ringer's solution, and cells were imaged **serum-starved at 37 °C without CO₂**. Imaging was confocal; activation was usually 458 nm every **10 s**. Pixel size was 0.31 µm in higher-resolution signaling panels and 1.24 µm in the population calcium assays. No seeding density or matrix coating is reported for these assays. ([PubMed Central (PMC)][13])

### Figure 1 — validating opto-RhoA

**1A–B.** iLID is CAAX anchored; SspB carries LARG-DH. Blue light recruits the GEF to membrane. **1C–D** locally illuminate HeLa cells from **300–600 s every 20 s**. Seven repeats are reported. LARG-DH recruitment rises rapidly to about **1.3× baseline [F]** and the Rhotekin RBD reporter follows more gradually to roughly **1.3–1.35× [F]**; both decline after illumination. The unilluminated region changes much less. ([PubMed Central (PMC)][13])

The key limitation is that this is local reporter recruitment/activity, not a calibrated RhoA-GTP concentration.

### Figure 2 — specificity and kinetics of six opto-GTPases

**2A–B** define the six switches and ddRFP biosensors.

**2C–H.** HeLa cells are imaged in L-15. For RhoA/Rac/Cdc42/Ral, activation is from **60–230 s**, every 10 s; Ras/Rap run to **350 s**. Each group has >5 cells, with individual cell numbers printed as points. Means ±SD are plotted. One-way ANOVA/Dunnett is used except the Ras panel, which uses Brown–Forsythe ANOVA + Dunnett T3. The reported omnibus p values are <0.0001 for several Rho-family panels but much weaker for Ras. ([PubMed][14])

An important supplementary characterization is that **opto-Rap membrane recruitment itself is slower**, taking about **50 s to maximum [R]**, whereas most other GEF constructs translocate with broadly similar kinetics. Thus differences in downstream time course are partly actuator-specific. ([ResearchGate][15])

### Figure 3 — RhoA uniquely induces acute calcium transients in the screen

This is the most quantitatively useful main figure.

**3A–B.** RPE1 and HeLa cells express R-GECO1 and one of the six opto-GTPases. Activation is **458 nm every 10 s for 2 min**. A calcium “responder” is defined as **F/F₀ > 2 [R]**. >150 cells/photoswitch from >3 independent experiments are scored. Only opto-RhoA produces a substantial responder fraction. Approximate bar heights are **RPE1 ~30–35% and HeLa ~23–25% [F]**; the other five GTPases are near zero to a few percent [F]. ([PubMed Central (PMC)][13])

The GEF-dead Y940A control is much cleaner numerically: **1/136 RPE1 cells (0.74%) and 0/133 HeLa cells [R]** respond. ([PubMed Central (PMC)][13])

**3C** shows the population images.

**3D.** Only responding cells enter the kinetic average, and each cell is rescaled so baseline =0% and its own maximum =100%. Therefore curves **cannot be used to compare calcium amplitude among cell types**. RPE1, HeLa and HEK293T begin increasing within ~10 s and peak around **20 s [R]**; MDCK begins at **20–30 s** and reaches maximum only around **2 min [R]**. HeLa has visibly broader temporal variability and a more sustained tail. ([PubMed Central (PMC)][13])

Supplementary Fig. S5 adds exact responder rates for the other two cell lines: **MDCK 30% (n=50)** and **HEK293T 41% (n=52) [R]**. HEK293T cells already show membrane blebbing before illumination, which the authors interpret as dark-state leakiness of opto-RhoA. ([PubMed Central (PMC)][13])

### Figure 4 — the coupling mechanism differs by cell type

Cells are tested with **50 µM Y-27632, 20 µM blebbistatin, 10 µM SKF96365, 100 µM 2-APB or 10 µM U73122**. Each condition has **>200 cells from three independent experiments**; responder definition remains F/F₀>2. One-way ANOVA with Dunnett comparisons is used. ([PubMed][14])

**4A–B, RPE1/HeLa.** ROCK or myosin-II inhibition has comparatively little effect, while 2-APB suppresses the response. U73122 appears ineffective in L-15 but abolishes the calcium response when tested in Ringer's solution; the authors attribute the medium dependence to chemical inactivation of U73122 by components of culture medium. Calcium-free Ringer substantially reduces responders, but **~10% remain [R]**, supporting an intracellular-store component. ([PubMed Central (PMC)][13])

**4C–D, MDCK/HEK293T.** Y-27632 and blebbistatin strongly suppress calcium responses, as does SKF96365, producing the opposite mechanistic profile. The bars are visibly near baseline under effective inhibitors, but I do not assign exact `F` percentages because the small plotted replicate count and SD make such readings unstable. ([PubMed Central (PMC)][13])

This figure is important for migration-model interpretation because “RhoA activation” does **not** imply a universal downstream mechanical/signaling program across cell backgrounds.

### Figure 5 — PLCε is necessary in RPE1/HeLa

**5A** validates three PLCε siRNAs by Western blot.

**5B.** More than 150 cells/condition, three experiments. In both RPE1 and HeLa, all three PLCε siRNAs collapse the opto-RhoA calcium-response fraction. Visually: RPE1 control **~30%**, siRNAs roughly **0–4% [F]**; HeLa control **~30–35%**, siRNAs roughly **0–4% [F]**. ANOVA p<0.0001. ([PubMed][14])

**5C–E** perform rescue with wild-type and mutant mouse PLCε in PLCε-depleted HeLa. Approximate responder fractions are **WT ~35%; lipase-dead F702A ~30–35%; H1433L ~0–3%; R2130L ~35–40%; ΔYins ~0–3% [F]**, depending on exact mutant annotation in the panel. The central conclusion from the authors' mutational logic is that the RhoA-responsive Y-insert/lipase-related functions are required, whereas some other PLCε regulatory domains are dispensable. The exact bar heights should not be promoted to reported values. ([PubMed][14])

### Figure 6 — plasma-membrane RhoA and phosphoinositides

**6A–C.** RhoA is retargeted either to the plasma membrane, TGN or Golgi. In RPE1 and HeLa, plasma-membrane opto-RhoA induces calcium in roughly **40–50% [F]**, whereas TGN/Golgi versions are approximately **0–10% [F]**. >150 cells/photoswitch, three experiments; ANOVA p=0.0005 and p=0.0003 for the two cell types. The site of activation therefore matters strongly. ([PubMed][14])

**6D–F.** An optogenetic 5-phosphatase is the positive control for PI(4,5)P₂ depletion. mCherry-PLCδ-PH rapidly redistributes with the 5-phosphatase, whereas opto-RhoA produces only a small change; the authors therefore describe effects on bulk phosphoinositide dynamics as marginal. Images are acquired every 10 s. ([PubMed][14])

### Figure 7 — the calcium pulse propagates to NFAT

Approximately **30% of opto-RhoA cells show NFAT nuclear translocation [R]**, comparable with the fraction showing calcium transients. PLCε knockdown eliminates the response, as do BAPTA-AM and cyclosporin A. Reported n values are approximately **61, 58, 72 cells** in the siRNA comparison and **69, 44, 50** in the drug comparison, across three experiments [P/caption]. The nuclear/cytosolic NFAT metric is calculated as the ratio at 5 min divided by the ratio at 0 min. ([PubMed Central (PMC)][13])

### Figure 8

Summary diagram only: RPE1/HeLa use RhoA→PLCε→IP₃R/ER Ca²⁺; MDCK/HEK293T instead use a ROCK–myosin pathway that the authors propose may couple contractile force to a mechanosensitive Ca²⁺ channel. ([PubMed Central (PMC)][13])

### Supplementary material

The PMC record exposes **Figures S1–S6 plus five videos**. The support PDF itself was not graphically rendered in this pass, so I do not claim figure-read values from it. However, integrated main-text descriptions provide several explicit supplementary results: S1 shows light-dose-dependent translocation with **2% 458-nm laser power sufficient for maximal recruitment [R]** and reversal within ~1 min; S2 shows broadly similar translocation kinetics except slow opto-Rap; S5 provides the exact MDCK/HEK responder fractions and prelight HEK blebbing; S6 supports the PLCε rescue experiments. ([PubMed Central (PMC)][13])

## Internal/provenance issues

The paper's responder fraction is thresholded at **F/F₀>2**, so modest but real calcium changes are deliberately counted as negative. Fig. 3D then includes **only responders** and renormalizes every cell to its own maximum, making it unsuitable for comparing response amplitude or population penetrance. HEK293T shows pre-illumination blebbing, indicating actuator leak. U73122 changes apparent effectiveness depending on buffer chemistry. Finally, the paper contains **no actual migration trajectory assay**, so it should constrain a migration model only through signaling/mechanical side effects of RhoA, not through speed/directness data.

---

# Important provenance / measurement qualifications

Across the four studies, the cleanest directly comparable quantities are **not all “migration” in the same sense**. Nanda gives 4-h trajectory **path length, Euclidean displacement and path directness**; Wu gives pulse-induced movement and an **angular polarity index**; de Seze mainly gives **local segmented area/membrane displacement and phenotype probability**; Inaba gives **signaling-response incidence and timing**, not migration. Converting these into one generic “motility” variable would throw away experimentally important definitions.

Several numbers above are intentionally tagged **F** because I read them from plotted means/distributions. They should be used as approximate model-screening targets or order-of-magnitude constraints, not as digitized ground truth. For fitting, the better path for Nanda and de Seze would be to use their released source data rather than these visual estimates: Nanda provides source-data files and analysis code, while de Seze provides per-figure source data/code and raw imaging through the BioImage Archive. ([PubMed Central (PMC)][2])

Supplement coverage is uneven. de Seze's supplementary figures were substantially inspectable through the integrated eLife presentation. For Nanda and Inaba I verified the supplementary file inventories and main-text descriptions but did not visually render the PDF supplements, and for Wu I inspected the supplementary index/methods and text-reported supplementary measurements but not the supplementary graphs themselves. I therefore did not manufacture F-type supplemental values.

---

# Authors’ interpretations and hypotheses

**Nanda et al.** interpret Rac→Rho activation as a mechanism coupling a protrusive Rac phase to a subsequent contractile/retractive Rho phase. ARHGEF11 and ARHGEF12 are proposed as Rac-responsive intermediates. In this view, stronger coupling shortens protrusion–retraction cycles and promotes frequent direction changes/exploratory migration, whereas depletion makes trajectories less extensive but straighter. They further speculate that delayed Rho self-inhibition/negative feedback may explain adaptation and permit a new protrusive cycle. ([PubMed Central (PMC)][2])

**Wu et al.** conclude that local Rac activity is sufficient to establish cell polarity and, in a permissive cell background, directional migration. Their myosin experiments lead them to distinguish local protrusion generation from whole-cell directional coordination: Rac can still induce protrusion when myosin function is impaired, but directionality deteriorates. They interpret spatial RhoA suppression as one component of Rac-mediated polarity. ([PubMed Central (PMC)][5])

**de Seze et al.** propose that PRG is capable of multiplexing two antagonist outputs. At low PRG abundance, optogenetic recruitment predominantly activates RhoA and retraction; at high abundance, PRG's PH domain sequesters active RhoA while PRG also activates Cdc42, shifting the outcome toward protrusion. Their minimal model interprets pulse-frequency sensitivity as a consequence of fast RhoA sequestration versus slower RhoA activation. They explicitly connect this abundance-dependent signaling behavior to situations of unusually high protein expression, including cancer, but that cancer relevance is a hypothesis rather than a demonstrated disease mechanism here. ([eLife][9])

**Inaba et al.** interpret acute RhoA→Ca²⁺ coupling as broadly present but mechanistically cell-type dependent. In RPE1/HeLa they propose plasma-membrane RhoA→PLCε→IP₃R/ER Ca²⁺ release, whereas in MDCK/HEK293T they propose RhoA→ROCK→myosin-II→mechanosensitive-channel coupling; they mention TRP-family channels as candidates because the pharmacology is inconsistent with Piezo being the sole route. ([PubMed Central (PMC)][13])

---

# Migration-relevant quantitative core

| Study            | System / perturbation                           | Exact metric                                  | Quantitative core                                                                                              |
| ---------------- | ----------------------------------------------- | --------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| **Nanda 2023**   | A431, spontaneous migration ± ARHGEF11/12 siRNA | **4-h path length**                           | control ~175–180 µm; si11/si12 ~145–150 µm **[F]**                                                             |
|                  |                                                 | **4-h Euclidean displacement**                | ~13–15 µm all groups **[F]**; comparatively little effect                                                      |
|                  |                                                 | **directionality = displacement/path length** | control ~0.075–0.08; si11 ~0.10–0.105; si12 ~0.09–0.095 **[F]**                                                |
|                  | A431 morphodynamics                             | **P–R cycle duration**                        | OE: control ~20, GEF11 ~15–16, GEF12 ~17 min **[F]**; KD: control ~14–15, si11 ~21–22, si12 ~17–18 min **[F]** |
|                  | Rac/Rho spontaneous cycles                      | **local edge velocity classification**        | protrusion >+0.075 µm/min; retraction <−0.075 µm/min **[R]**                                                   |
|                  | PA-Rac1 across A431                             | response heterogeneity                        | 67% positive: 29% sustained + 38% transient; 25% none; 8% negative **[R]**                                     |
| **Wu 2009**      | PA-Rac1 MEFs                                    | **pulse-induced movement**                    | **0.8 µm/pulse**, 2-min pulse interval, n=6 **[R]**                                                            |
|                  | Single local Rac pulse                          | polarized morphology                          | 50/50 PA-Rac1 cells protruded locally/retracted opposite; 25/25 T17N showed reverse polarity **[R]**           |
|                  | Myosin inhibition                               | **polarity index = cos θ**                    | strong decrease with ML-7/blebbistatin despite relatively preserved protrusion; n>25 **[R/P]**                 |
|                  | PA-Rac1 actuator                                | spatial/temporal resolution                   | D=0.55 µm²/s; FRAP t½=12.1 s; dark recovery t½=43 s; predicted adjacent-spot activation ≤7.5% **[R]**          |
| **de Seze 2025** | RPE1 optoPRG                                    | local protrusion versus retraction            | below ~40 a.u.: ~95% retract; above ~40 a.u.: ~85% protrude **[R]**                                            |
|                  | optoPRG + nonrecruitable PRG                    | protruding fraction in low-optoPRG cells      | **0% → 45% [R]**                                                                                               |
|                  | high- vs low-expression states                  | prestimulation projected area                 | high/protruding state ~**1.5× larger [R]**                                                                     |
|                  | pulse-frequency switch                          | within-cell phenotype control                 | **3/90 cells** successfully switched repeatedly; 30-s versus 15-s pulse intervals **[R]**                      |
| **Inaba 2021**   | RPE1/HeLa opto-RhoA                             | Ca²⁺ responder, **F/F₀>2**                    | roughly ~30–35% RPE1 and ~23–25% HeLa **[F]**; GEF-dead control 1/136 and 0/133 **[R]**                        |
|                  | MDCK/HEK293T                                    | Ca²⁺ responder                                | **30% of 50 MDCK; 41% of 52 HEK293T [R]**                                                                      |
|                  | responding cells only                           | Ca²⁺ latency                                  | RPE1/HeLa/HEK: starts ≤10 s, peak ~20 s; MDCK starts 20–30 s, peak ~2 min **[R]**                              |
|                  | downstream signaling                            | NFAT nuclear translocation                    | ~**30% [R]** after opto-RhoA; PLCε/Ca²⁺/calcineurin dependent                                                  |

For a migration-model screening suite, the **highest-value hard constraints from this set** are probably: Nanda's path-length/directness decoupling and P–R cycle times; Wu's 0.8 µm/pulse plus the near-deterministic local/opposite-side polarization; de Seze's abundance-dependent sign switch and its low within-cell switching success rate; and Inaba's warning that identical upstream RhoA activation can engage quite different downstream pathways and timescales across cell types.

[1]: https://www.nature.com/articles/s41467-023-43875-y "https://www.nature.com/articles/s41467-023-43875-y"
[2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10724141/ "https://pmc.ncbi.nlm.nih.gov/articles/PMC10724141/"
[3]: https://www.researchgate.net/figure/Identification-of-Arhgef11-and-Arhgef12-as-Rac-effectors-in-local-cell_fig5_376556871 "https://www.researchgate.net/figure/Identification-of-Arhgef11-and-Arhgef12-as-Rac-effectors-in-local-cell_fig5_376556871"
[4]: https://www.nature.com/articles/nature08241 "https://www.nature.com/articles/nature08241"
[5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2766670/ "https://pmc.ncbi.nlm.nih.gov/articles/PMC2766670/"
[6]: https://go.gale.com/ps/i.do?asid=f2139121&id=GALE%7CA207829295&it=r&p=AONE&u=googlescholar&v=2.1 "https://go.gale.com/ps/i.do?asid=f2139121&id=GALE%7CA207829295&it=r&p=AONE&u=googlescholar&v=2.1"
[7]: https://www.researchgate.net/figure/nhibition-of-RhoA-by-PA-Rac1a-HeLa-cells-expressing-RhoA-biosensor-and-either-PA-Rac1_fig1_26753666 "https://www.researchgate.net/figure/nhibition-of-RhoA-by-PA-Rac1a-HeLa-cells-expressing-RhoA-biosensor-and-either-PA-Rac1_fig1_26753666"
[8]: https://www.researchgate.net/figure/Crystallization-and-structural-modelling-of-PA-Rac1a-Dark-state-crystal-structure-of_fig2_26753666 "https://www.researchgate.net/figure/Crystallization-and-structural-modelling-of-PA-Rac1a-Dark-state-crystal-structure-of_fig2_26753666"
[9]: https://elifesciences.org/articles/93180 "https://elifesciences.org/articles/93180"
[10]: https://pmc.ncbi.nlm.nih.gov/articles/PMC12113259/ "https://pmc.ncbi.nlm.nih.gov/articles/PMC12113259/"
[11]: https://elifesciences.org/articles/93180/figures "https://elifesciences.org/articles/93180/figures"
[12]: https://pubmed.ncbi.nlm.nih.gov/40421714/ "https://pubmed.ncbi.nlm.nih.gov/40421714/"
[13]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7949103/ "https://pmc.ncbi.nlm.nih.gov/articles/PMC7949103/"
[14]: https://pubmed.ncbi.nlm.nih.gov/33453281/ "https://pubmed.ncbi.nlm.nih.gov/33453281/"
[15]: https://www.researchgate.net/publication/348449609_Optogenetic_control_of_small_GTPases_reveals_RhoA_mediates_intracellular_calcium_signaling "https://www.researchgate.net/publication/348449609_Optogenetic_control_of_small_GTPases_reveals_RhoA_mediates_intracellular_calcium_signaling"
