# Jemaà et al. 2017 — structured data summary

**Jemaà M et al.** *Heterogeneity in sarcoma cell lines reveals enhanced motility of tetraploid versus diploid cells.* **Oncotarget** 8:16669–16689 (2017; published online Dec 27, 2016). ([Oncotarget][1])

* **R** = explicitly reported in text/caption.
* **P** = numerical value printed directly on the figure.
* **F** = approximate value read visually from the figure.

I keep measurements separate from the authors' mechanistic interpretation.

---

## Experimental system and provenance of the diploid/tetraploid comparison

The study uses two human **undifferentiated pleomorphic sarcoma** cell lines, **MFH152** and **MFH137**, both established from resected metastatic tumors. The parental cultures are highly heterogeneous mixtures of diploid, tetraploid and aneuploid cells. IMR90 fibroblasts serve as a non-transformed reference. All were grown in DMEM + 10% FCS. ([Oncotarget][2])

This is importantly **not an experimentally induced WGD system**. The authors FACS-sorted cells from the pre-existing heterogeneous MFH137/MFH152 cultures on the basis of **FSC/SSC cell size/granularity**, deposited single cells into 96-well plates and clonally expanded them. For each parental line they sorted **1,000 “small” and 1,000 “big” single cells**. Most surviving clones were aneuploid; they ultimately obtained **>20 stable diploid clones per line**, **5 stable tetraploid MFH152 clones**, and **3 stable tetraploid MFH137 clones**. Selected clones retained stable DNA content for up to ~3 months. ([Oncotarget][2])

Therefore the correct comparison is:

> **clonally isolated diploid versus tetraploid subpopulations from the same heterogeneous tumor-derived cell line**

rather than a demonstrated diploid → WGD → tetraploid lineage. The authors explicitly acknowledge later that transcriptomic similarity is consistent with a shared origin but **does not prove that the tetraploid clones derive from the diploid clones**. ([Oncotarget][2])

This is consequently a weaker causal WGD design than Jemaà 2023, although the replication across multiple independent clones and across **two sarcoma lines** is useful.

---

# Figure 1 — baseline motility/invasion of the heterogeneous parental sarcoma lines

These experiments precede clone separation. They establish that the original MFH152/MFH137 populations are unusually motile relative to IMR90.

## 1A — scratch-wound migration

**Assay:** confluent monolayers scratched with a **1,000-μL pipette tip** and imaged **hourly**. “Grouped cell speed” is calculated from wound-closure rate. ([Oncotarget][2])

**F:**

| Cell line | Grouped migration speed |
| --------- | ----------------------: |
| IMR90     |            **~21 μm/h** |
| MFH152    |            **~46 μm/h** |
| MFH137    |            **~61 μm/h** |

MFH152 and MFH137 are both *** versus IMR90.

Thus the parental sarcoma populations migrate approximately **2.2× and 2.9× faster**, respectively, in this collective assay.

---

## 1B–C — Oris™ migration

Cells were seeded at **2×10⁴ cells/well** in 96-well plates around silicone stoppers and grown **overnight**. After stopper removal they migrated into the initially empty circular region for **24 h**. ([Oncotarget][2])

### 1C — filled area after 24 h

**F:**

| Cell line |        Filled area |
| --------- | -----------------: |
| IMR90     |      **~0.34 mm²** |
| MFH152    |      **~0.54 mm²** |
| MFH137    | **~0.66–0.67 mm²** |

Both sarcoma lines *** versus IMR90.

Panel 1B is useful qualitatively for density: the cells form a dense annulus immediately outside the stopper-defined migration region, i.e. this is a **collective/high-density migration assay**, not sparse-cell tracking. ([Oncotarget][2])

---

## 1D — single-cell migration on fibronectin tracks

Cells were plated at approximately **20 cells/mm²**, specifically to minimize cell-cell contacts, on fibronectin-coated linear CYTOO™ tracks. Imaging was **every 15 min for 24 h**. ([Oncotarget][2])

**F:**

| Cell line | Mean individual-cell speed |
| --------- | -------------------------: |
| IMR90     |               **~12 μm/h** |
| MFH152    |            **~28–29 μm/h** |
| MFH137    |            **~70–72 μm/h** |

Error bars are large, particularly for MFH137 (~±20 μm/h). Both sarcoma lines *** versus IMR90.

This is a genuinely high-frequency single-cell assay, but it compares the **heterogeneous parental cultures**, not diploid versus tetraploid clones.

---

## 1E–F — 3D collagen invasion

Cells were suspended at **10⁵ cells/mL** in **2.3 mg/mL bovine collagen**, 100 μL/well. After polymerization, serum was added above the gel to **5% final concentration**. After 24 h, cells were counted at the bottom of the well (**0 μm**) and at **50 μm depth**; invasion index = cells at 50 μm / cells at 0 μm. ([Oncotarget][2])

### 1F

**F:**

| Cell line | % invasive cells / plotted invasion measure |
| --------- | ------------------------------------------: |
| IMR90     |                                **~0.2–0.3** |
| MFH152    |                                      **~6** |
| MFH137    |                                     **~11** |

Both *** versus IMR90.

The y-axis is labelled **“% invasive cells”**, although Methods define the underlying metric as the **50 μm/0 μm cell-count ratio**. ([Oncotarget][2])

---

# Figure 2 — heterogeneity of the parental MFH137/MFH152 populations

This figure explains why the authors subsequently isolate clones.

## 2A — cell-cycle duration

Cells were followed for up to **72 h**.

Exact means are printed:

| Cell line | Mean cell-cycle duration |  n |
| --------- | -----------------------: | -: |
| IMR90     |               **17.2 h** | 18 |
| MFH152    |               **32.3 h** | 80 |
| MFH137    |               **31.9 h** | 22 |

**R:** MFH152 individual cycles range **21–50 h** and MFH137 **22–46 h**. ([Oncotarget][2])

---

## 2B — mitosis duration

Exact means printed:

| Cell line | Mean mitosis duration |   n |
| --------- | --------------------: | --: |
| IMR90     |          **26.6 min** |  77 |
| MFH152    |         **93.15 min** | 177 |
| MFH137    |          **70.5 min** | 102 |

**R ranges:** MFH152 **24–432 min**; MFH137 **16–368 min**. ([Oncotarget][2])

The distributions have conspicuous long upper tails, especially MFH152.

---

## 2C–D — DNA-content distributions

PI-flow cytometry shows:

* a conventional 2n/4n profile in IMR90;
* substantial populations with **≥4n DNA** in both sarcoma cultures;
* a high-DNA “shoulder” beyond the conventional 4n peak. ([Oncotarget][2])

Panel 2D quantifies these distributions, but the principal point is that the parental sarcoma populations cannot be treated as one ploidy class.

---

## 2E–F — chromosome counts

Metaphase spreads directly show near-diploid and near-tetraploid cells.

Sample sizes:

* IMR90: **n=17**
* MFH152: **n=140**
* MFH137: **n=141**

IMR90 is tightly clustered at approximately **46 chromosomes**.

MFH152 and MFH137 contain:

* a cluster around **~40–50 chromosomes**;
* cells around **~90–100 chromosomes**;
* numerous intermediate chromosome counts.

The authors classify these as diploid, tetraploid and aneuploid subpopulations. ([Oncotarget][2])

---

# Figure 3 — characterization of isolated MFH152 diploid/tetraploid clones

For MFH152 the principal exemplar clones are:

* **Cl1, Cl2** = diploid;
* **Cl7, Cl8** = tetraploid.

The broader motility experiment later includes Cl1–4 and Cl5–9.

## 3A–B — DNA content

Cl1 shows the expected diploid **2n G1 / 4n G2-M** profile.

Cl7 shows the corresponding tetraploid **4n G1 / 8n G2-M** profile.

The authors selected clones whose G2/M populations were **<20%**, in order to minimize ambiguity between a diploid G2 cell and a tetraploid G1 cell. ([Oncotarget][2])

---

## 3C — cell-cycle duration

Individual points are shown.

**R/P/F:**

* Cl1 diploid: mean approximately **22 h**, **n=36**
* Cl7 tetraploid: mean approximately **26 h**, **n=32**
* difference **p<0.05**

The distributions overlap strongly, but Cl7 is shifted upward by roughly **4 h (~18%)**. ([Oncotarget][2])

---

## 3D — mitosis duration

* Cl1: central value roughly **~80 min**, **n=61**
* Cl7: roughly **~65–70 min**, **n=71**
* **n.s.**

Both distributions have strong right tails; Cl1 includes values near ~200 min and Cl7 one around ~240 min. ([Oncotarget][2])

---

## 3E–F — phospho-H3 mitotic fraction

**F:**

| Clone | pH3+ fraction |
| ----- | ------------: |
| Cl1   |     **~2.7%** |
| Cl2   | **~2.4–2.5%** |
| Cl7   |     **~1.9%** |
| Cl8   |     **~1.8%** |

The tetraploid clones carry * markers relative to Cl1, but absolute differences are <1 percentage point. ([Oncotarget][2])

---

## 3G–H — spindle-assembly checkpoint

After **16 h** with nocodazole (200 nM) or paclitaxel (100 nM):

**F, nocodazole:**

* Cl1: ~62%
* Cl2: ~82%
* Cl7: ~49%
* Cl8: ~45%
* parental MFH152: ~22%

**F, paclitaxel:**

* Cl1: ~60%
* Cl2: ~59%
* Cl7: ~39%
* Cl8: ~46%
* MFH152: ~31%

Thus the diploid clones initially arrest more strongly.

For prolonged nocodazole, the text gives explicit **50% escape times**:

| Clone          | Approx. time for 50% SAC escape |
| -------------- | ------------------------------: |
| Cl1 diploid    |                        **63 h** |
| Cl2 diploid    |                        **46 h** |
| Cl7 tetraploid |                        **36 h** |
| Cl8 tetraploid |                        **30 h** |

([Oncotarget][2])

---

# Figure 4 — the central diploid-versus-tetraploid migration data

This is the key figure for our purposes.

## 4A–B — Oris collective migration, MFH152

Four diploid clones (**Cl1–4**) and five tetraploid clones (**Cl5–9**) were tested independently.

Assay:

* **2×10⁴ cells/well**
* overnight growth around silicone stopper;
* stopper removed;
* **24 h migration**
* DAPI segmentation;
* authors deliberately count **nuclei entering the migration zone**, rather than occupied area, to avoid confounding by the larger size of tetraploid cells.
* **n=3 experiments.** ([Oncotarget][2])

### 4B — cells inside migration zone

**F:**

| Ploidy     | Clone | Cells inside migration zone |
| ---------- | ----- | --------------------------: |
| Diploid    | Cl1   |                    **~25%** |
|            | Cl2   |                    **~25%** |
|            | Cl3   |                    **~27%** |
|            | Cl4   |                    **~15%** |
| Tetraploid | Cl5   |                    **~39%** |
|            | Cl6   |                    **~44%** |
|            | Cl7   |                    **~38%** |
|            | Cl8   |                    **~39%** |
|            | Cl9   |                    **~34%** |

The five tetraploid clones therefore cluster around approximately **34–44%**, whereas the four diploid clones span approximately **15–27%**.

The authors report ***p<0.001**, with the caption stating **each tetraploid clone compared with each diploid clone** by ANOVA. ([Oncotarget][2])

The mean across the visually read values is approximately:

* diploid ≈ **23%**
* tetraploid ≈ **39%**

or roughly **1.7× more cells** within the migration zone, reasonably consistent with the manuscript's description of “roughly twice as many.”

---

## MFH137 replication — Supplementary Fig. S2E–F

The authors repeat the same Oris experiment using independently isolated **diploid and tetraploid MFH137 clones** and report the same direction: tetraploid clones enter the migration zone substantially more frequently. ([Oncotarget][2])

I do **not** have a sufficiently reliable rendering of Supplementary Fig. S2 to figure-read its individual bars quantitatively, so I would record this as a **qualitative/author-reported independent replication**, not invent numerical values.

---

# Figure 4C–D — low-density individual-cell migration

This is probably the single most valuable experiment in the paper for us.

### Assay conditions

For the diploid/tetraploid clone comparison:

* cells seeded at **30,000 cells per 6-well plate**
* explicitly described as **very low density / non-confluent**
* time-lapse imaging every **30 min**
* nominal duration **24 h**
* 20× objective
* tracks reconstructed using ImageJ. ([Oncotarget][2])

There is a small reporting discrepancy: the Results say the individual cells were followed over a **20-h time range**, whereas the Figure 4 caption and Methods describe **24 h**. ([Oncotarget][2])

This is notably better temporal resolution than the **1-h cadence in Jemaà 2023**, although still coarser than Xuan's 10-min movies.

### 4C — example trajectories

The panel shows reconstructed tracks for:

* diploid Cl1
* tetraploid Cl8.

Tetraploid tracks are visibly longer and often spatially extensive.

However, these are illustrative tracks. No coordinate data, displacement, persistence or turning-angle statistics are presented.

### 4D — individual-cell speed

**Scatter plot**, y-axis 0–80 μm/h.

**F mean/central values:**

| Clone          | Approx. mean speed |
| -------------- | -----------------: |
| Cl1 diploid    |       **~16 μm/h** |
| Cl2 diploid    |    **~15–16 μm/h** |
| Cl7 tetraploid |    **~24–25 μm/h** |
| Cl8 tetraploid |       **~27 μm/h** |

Thus the average shift is roughly **1.5–1.8×**.

Statistical annotations:

* Cl7 versus diploid reference: **p<0.05**
* Cl8: ***p<0.001**
* caption gives **n=3 experiments**. ([Oncotarget][2])

The individual-cell distributions contain more information than the means.

**Cl1/Cl2:**

* dense bulk roughly **10–22 μm/h**
* most observations <25 μm/h
* upper points around **27–29 μm/h**.

**Cl7:**

* bulk roughly **15–35 μm/h**
* several observations near **35–38 μm/h**.

**Cl8:**

* bulk approximately **15–35 μm/h**
* several cells around **45–55 μm/h**
* extreme observations around **~60 μm/h**.

So the tetraploid phenotype includes both an upward shift of the central distribution and a **high-speed tail** largely absent from the diploid clones.

The paper does **not specify whether “individual-cell speed” is total track length / observation time or an average of framewise instantaneous speeds**.

---

# Figure 4E — Matrigel invasion

This is a genuine ECM invasion assay, unlike the somewhat ambiguous transwell terminology in Jemaà 2023.

### Assay conditions

* serum deprivation for **24 h** beforehand;
* CIM-Plate membrane coated with **300 μg/mL Matrigel**;
* **20,000 cells/well**
* serum-free medium in upper chamber;
* invasion monitored hourly for **24 h** with xCELLigence;
* **n=4**;
* plotted “Cell Index Invasion” = cell index of **Matrigel-coated / uncoated** wells at the corresponding time. ([Oncotarget][2])

Only:

* parental MFH152;
* diploid **Cl1**;
* tetraploid **Cl7**

are compared.

### Approximate trajectory

| Time | Cl1 diploid | MFH152 parent | Cl7 tetraploid |
| ---: | ----------: | ------------: | -------------: |
| 10 h |   **~0.48** |         ~0.59 |      **~0.78** |
| 12 h |       ~0.50 |         ~0.65 |      **~0.86** |
| 14 h |       ~0.58 |         ~0.72 |      **~0.93** |
| 16 h |       ~0.65 |         ~0.76 |      **~0.97** |
| 18 h |       ~0.68 |         ~0.84 |      **~1.00** |
| 20 h |       ~0.73 |         ~0.86 |      **~1.00** |

**F from Fig. 4E.**

The rank order is stable throughout:

**tetraploid Cl7 > parental heterogeneous MFH152 > diploid Cl1**. ([Oncotarget][2])

---

# Figure 5 — molecular differences between diploid and tetraploid clones

## 5A — parental sarcoma cytoskeleton / KIF11

Mitotic MFH137/MFH152 cells show prominent membrane blebbing, abnormal attachment to fibronectin micropatterns and strong KIF11 staining. These are representative fluorescence images (**n>5**) rather than quantified measurements. ([Oncotarget][2])

---

## 5B — RNA-seq similarity

RNA-seq was performed on:

* diploid MFH152 Cl1/Cl2
* tetraploid Cl7/Cl8
* IMR90.

Reported correlations:

| Comparison            |          r |
| --------------------- | ---------: |
| diploid Cl1 vs Cl2    | **0.9901** |
| tetraploid Cl7 vs Cl8 | **0.9939** |
| diploid vs tetraploid | **0.9701** |
| diploid vs IMR90      | **0.9114** |
| tetraploid vs IMR90   | **0.9046** |

Thus the diploid/tetraploid clones are transcriptionally substantially closer to one another than either is to fibroblasts. The authors use this to argue that their large phenotypic difference is occurring against a broadly similar sarcoma transcriptional background, while explicitly noting that it **doesn't establish ancestry**. ([Oncotarget][2])

---

## Table 1 — CINSARC transcripts

A striking point is that the increased protein abundance described below is **not simply mirrored by uniformly higher transcript abundance in tetraploids**.

Examples, normalized to IMR90 = 1:

| Gene  | D Cl1 | D Cl2 |   T Cl7 |   T Cl8 |
| ----- | ----: | ----: | ------: | ------: |
| AURKA |   1.8 |   2.1 | **1.6** | **1.5** |
| AURKB |   3.7 |   4.4 | **3.0** | **3.1** |
| BUB1B |   1.8 |   2.2 | **1.7** | **1.7** |
| KIF11 |   2.2 |   2.1 | **1.9** | **1.7** |
| KIF23 |   2.6 |   2.8 | **3.1** | **3.0** |
| PRC1  |   2.0 |   2.1 | **2.2** | **1.9** |
| TPX2  |   4.2 |   4.5 | **4.0** | **3.5** |
| TTK   |   2.7 |   3.0 | **3.3** | **3.2** |

([Oncotarget][2])

Many CINSARC transcripts are similarly expressed in diploids and tetraploids, and several are actually somewhat lower in the tetraploid clones.

---

## 5C — protein abundance

Western blots compare **three diploid clones (Cl1–3)** with **three tetraploid clones (Cl6–8)**.

The authors report visibly higher protein abundance in tetraploids for:

* **BUBR1**
* **KIF11**
* **TTK**
* **TPX2**
* **KIF23**
* **PRC1**
* **PLK1**
* **AURKA**
* **CCNB1**
* **CDK1**

while several non-CINSARC proteins such as **EB1, CFL1, JNK1/JNK2** are comparatively similar. CCNE1 is relatively higher in the diploids. ([Oncotarget][2])

There is **no densitometric quantification**, so these differences remain visually assessed western-blot differences.

---

## 5D — KIF11 localization

Interphase immunofluorescence shows more prominent cytoplasmic KIF11 in tetraploid Cl7/Cl8 than Cl1/Cl2. The paper particularly notes KIF11 staining in **lamellipodial/peripheral regions**.

Again, no fluorescence-intensity quantification is presented. ([Oncotarget][2])

## 5E — AURKA localization

Tetraploid mitoses show unusually prominent AURKA staining at the **midzone/midbody**, whereas the diploid examples do not. No numerical quantification. ([Oncotarget][2])

---

# Table 2 / Figure 6 — pharmacological perturbation of migration

The authors initially screened **22 mitotic/cytoskeletal inhibitors**, each at **0.1, 1 and 10 μM**, excluding treatments producing >30% cytotoxicity. Migration was measured using the same Oris nuclear-count assay. ([Oncotarget][2])

Only two compounds met their prespecified criterion of **≥40% migration inhibition with acceptable toxicity** at 10 μM:

| Drug, 10 μM   | Diploid migration inhibition | Tetraploid migration inhibition |
| ------------- | ---------------------------: | ------------------------------: |
| **Reversine** |                     **0.40** |                        **0.47** |
| **SP600125**  |                     **0.44** |                        **0.47** |

Other usable compounds generally produced ≤0.25 inhibition under the conditions tested. ([Oncotarget][2])

---

## Figure 6A — toxicity

Three diploid clones (Cl1–3) and three tetraploid clones (Cl6–8) are tested with **10 μM reversine or SP600125 for 24 h**.

PI positivity + low DiOC6 indicate dead/dying cells.

**F:** total affected fractions remain roughly **~10–24% across all groups**, consistent with the authors' <30% selection criterion. There is no obvious preferential catastrophic toxicity in the tetraploid clones. ([Oncotarget][2])

---

## Figure 6B — migration after broad kinase inhibition

Percentage of cells entering the migration zone:

| Clone |     Control | Reversine | SP600125 |
| ----- | ----------: | --------: | -------: |
| D Cl1 |    **~25%** |      ~19% |     ~15% |
| D Cl2 |    **~25%** |      ~14% |     ~14% |
| D Cl3 | **~26–27%** |      ~11% |     ~12% |
| T Cl6 |    **~44%** |      ~17% |  ~11–12% |
| T Cl7 |    **~38%** |      ~19% |  ~15–16% |
| T Cl8 |    **~39%** |      ~17% |     ~16% |

**F**, n=3. ([Oncotarget][2])

Thus the untreated tetraploid/diploid difference is large, but after either broad inhibitor the groups collapse into roughly the same **~11–19% migration range**.

The figure marks all treated comparisons significant within clone and reports **### p<0.001** for greater inhibition of tetraploid versus diploid clones. ([Oncotarget][2])

---

# Specific-target follow-ups in the supplement

This part is important because it substantially qualifies the mechanistic claim.

The authors targeted **AURKA, AURKB, TTK and KIF11 individually** using siRNA. They report:

* near-complete protein depletion with **5 nM siRNA**;
* <30% toxicity even at **100 nM**;
* **no significant migration inhibition** at either 5 or 100 nM for any individual target. ([Oncotarget][2])

They then used highly selective inhibitors of the same proteins:

* MLN8237 / alisertib — AURKA
* AZD1152-HQPA / barasertib — AURKB
* SB743921 — KIF11
* Mps-BAY2a — TTK.

At concentrations producing the expected cell-cycle phenotypes while maintaining <30% cytotoxicity, **none significantly inhibited migration in either diploid or tetraploid clones**. ([Oncotarget][2])

So experimentally:

> broad multikinase drugs suppress migration strongly; individual depletion/inhibition of the nominated targets does not reproduce that effect.

---

# Figure 7 — 3D spheroid invasion

This experiment does **not compare diploid and tetraploid cells**.

The authors attempted spheroids from both, but report that tetraploid MFH152 cells failed to form stable spheroids even after **48–72 h**, instead producing loose clumps. The experiment was therefore restricted to a **diploid MFH152 clone**. ([Oncotarget][2])

### Assay

* **1,000 cells** used to form spheroid;
* 2.4% methylcellulose;
* 24 h spheroid formation;
* embedded in collagen;
* only **0.1% FCS** to minimize proliferation;
* ±10 μM reversine or SP600125;
* imaged **hourly**;
* full Z-stack each time;
* measurements through 24 h. ([Oncotarget][2])

## 7A

Representative frames at:

* 0 h
* 6 h
* 12 h
* 18 h
* 24 h.

Untreated cells progressively radiate outward from the spheroid into collagen. Both drug-treated spheroids remain compact.

## 7B — invasion area

**F:**

| Treatment |    Invasion area |
| --------- | ---------------: |
| Control   | **~4.6–4.7 mm²** |
| Reversine |     **~0.1 mm²** |
| SP600125  |   **near 0 mm²** |

Both *** versus control. ([Oncotarget][2])

This is therefore strong evidence that these compounds inhibit 3D invasion in the diploid clone, but it provides **no direct tetraploid/diploid 3D comparison**.

---

# Supplementary movies — what they actually contain

The paper has **eight AVI supplementary movies** plus the supplementary PDF/tables. ([PubMed Central (PMC)][3])

The assignments from the main text are:

| Movie(s)  | Experiment                                                  |
| --------- | ----------------------------------------------------------- |
| **S1**    | parental MFH137/MFH152/IMR90 scratch migration              |
| **S2**    | parental cells migrating individually on fibronectin tracks |
| **S3–S5** | long-term parental cell-cycle/mitosis imaging               |
| **S6**    | untreated diploid MFH152 spheroid invasion                  |
| **S7–S8** | reversine/SP600125 spheroid invasion                        |

([Oncotarget][2])

An important limitation for our purposes is that **there does not appear to be a supplementary movie corresponding to the key diploid-versus-tetraploid single-cell comparison in Fig. 4C–D**. The paper provides example reconstructed trajectories, but not a published movie or numerical trajectory dataset for those clone tracks.

---

# Important provenance and measurement qualifications

### 1. Tetraploidy was selected, not experimentally generated

Unlike Jemaà 2023, these tetraploid cells were already present within long-cultured, genomically heterogeneous sarcoma populations. Sorting by cell size and subsequent clone selection produced stable near-diploid and near-tetraploid lines, but does not establish that WGD itself caused the motility phenotype. The authors explicitly propose a future experiment using experimentally generated tetraploids to answer that question. ([Oncotarget][2])

### 2. There is substantial independent-clone replication

The strongest part of the design is Fig. 4B:

* **4 independent diploid MFH152 clones**
* **5 independent tetraploid MFH152 clones**

all showing separated migration ranges.

The same qualitative result is additionally reported in **MFH137 clones**. ([Oncotarget][2])

So the phenotype is not driven solely by one unusually motile tetraploid clone.

### 3. Individual-cell imaging is reasonably high frequency

The key clone experiment is **30-min imaging for ~20–24 h**, at low density (**30,000 cells/6-well**) and explicitly non-confluent. ([Oncotarget][2])

That is enough to contain considerably richer migration information than the authors extract.

However, they publish only:

* reconstructed example trajectories;
* average cell speed.

They do **not** report:

* net displacement;
* path length separately;
* persistence;
* turning-angle distribution;
* velocity autocorrelation;
* pause/move state occupancy;
* directional changes;
* MSD;
* cell morphology;
* protrusion dynamics.

### 4. Tetraploid cells are visibly larger, but the authors deliberately avoid area-based migration quantification

For collective migration they count **nuclei entering the migration region** rather than area coverage specifically because tetraploid cells are larger. ([Oncotarget][2])

The paper unfortunately does not provide a quantitative cell-area or volume comparison between the isolated diploid/tetraploid clones.

### 5. The 3D invasion result in Fig. 4E is genuine ECM invasion

Unlike some simple transwell assays, the xCELLigence membrane is explicitly coated with **300 μg/mL Matrigel**, and the plotted invasion index is normalized to migration through an uncoated membrane. ([Oncotarget][2])

### 6. The drug experiments do not identify a single causal protein

This is important enough not to lose in a cross-paper database: broad inhibitors strongly suppress migration, but **specific siRNA and selective inhibition of AURKA, AURKB, TTK or KIF11 individually do not**. ([Oncotarget][2])

---

# Authors' interpretations and hypotheses — separate from measurements

The authors interpret the increased motility/invasion of the tetraploid clones as being associated with their higher steady-state levels of a suite of **mitotic proteins**, particularly proteins that can also interact with cytoskeletal machinery outside mitosis. ([Oncotarget][2])

They propose that inappropriate persistence/activity of mitotic regulators during **interphase** alters actin/microtubule organization and thereby increases migration. KIF11 localization at cell protrusions and AURKA persistence at the midbody are used as supporting observations. ([Oncotarget][2])

Because neither siRNA nor selective single-target inhibitors reproduce the migration suppression caused by reversine/SP600125, the authors ultimately hypothesize that **several deregulated mitotic proteins act together**, rather than one of AURKA/AURKB/TTK/KIF11 being individually responsible. ([Oncotarget][2])

Most importantly, they explicitly **do not claim to have established that tetraploidization itself causes increased migration**. They state that experimentally generated tetraploids from a non-CINSARC background would be required to address that question. ([Oncotarget][2])

---

## Migration-relevant quantitative core

| Measurement                    |         Diploid |      Tetraploid | Context                               |
| ------------------------------ | --------------: | --------------: | ------------------------------------- |
| Oris migrated cells, MFH152    |      **15–27%** |      **34–44%** | 4 D vs 5 T clones, 24 h               |
| Individual speed               | **~15–16 μm/h** | **~24–27 μm/h** | low density; 30-min cadence; ~20–24 h |
| Individual speed upper tail    |        ~30 μm/h | **~50–60 μm/h** | Fig. 4D                               |
| Matrigel invasion index, 10 h  |       **~0.48** |       **~0.78** | Cl1 vs Cl7                            |
| Matrigel invasion index, 20 h  |       **~0.73** |       **~1.00** | Cl1 vs Cl7                            |
| Cell-cycle length              |       **~22 h** |       **~26 h** | Cl1 vs Cl7                            |
| Oris migration after reversine |         ~11–19% |         ~17–19% | three clones each                     |
| Oris migration after SP600125  |         ~12–15% |         ~12–16% | three clones each                     |

([Oncotarget][2])


[1]: https://www.oncotarget.com/article/14291/text/?utm_source=chatgpt.com "Heterogeneity in sarcoma cell lines reveals enhanced motility of tetraploid versus diploid cells | Oncotarget"
[2]: https://www.oncotarget.com/article/14291/text/ "Heterogeneity in sarcoma cell lines reveals enhanced motility of tetraploid versus diploid cells | Oncotarget"
[3]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5369993/?utm_source=chatgpt.com "Heterogeneity in sarcoma cell lines reveals enhanced motility of tetraploid versus diploid cells - PMC"
