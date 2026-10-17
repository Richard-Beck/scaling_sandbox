# Jemaà et al. 2023 — structured data summary

**Jemaà M et al.** *Tetraploidization Increases the Motility and Invasiveness of Cancer Cells.* International Journal of Molecular Sciences 24, 13926 (2023). ([PubMed Central (PMC)][1])


* **R** = explicitly reported in text/caption.
* **P** = numerical value printed directly on the figure.
* **F** = approximate value read visually from the figure.

I keep experimental observations separate from the authors' mechanistic interpretation.

---

## Experimental system and how the diploid/tetraploid clones were generated

This study uses two cancer-cell systems:

**RKO human colon carcinoma.** This is the system used for essentially all of the new **in-vitro migration and invasion experiments** in Figures 1–4.

**MFH152 malignant fibrous histiocytoma / soft-tissue sarcoma.** These clones are characterized partly in the supplement, included in the proliferation experiment, and used for the **in-vivo tail-vein metastasis experiment**. ([PubMed Central (PMC)][1])

The RKO system is substantially cleaner with respect to WGD than the Xuan PGCC experiments. Parental RKO cells were exposed to **600 ng/mL cytochalasin D for 48 h**, blocking cytokinesis. After washout, cells were stained with **1 μg/mL Hoechst 33342** and sorted by DNA content. Because diploid G2/M cells and tetraploid G1 cells both have a 4n-equivalent DNA content, the authors instead sorted **diploid G1 cells (2n DNA)** versus **tetraploid G2/M cells (8n-equivalent DNA)** and established single-cell clones. The study ultimately uses **2 diploid RKO clones (D1, D2) and 2 tetraploid RKO clones (T1, T2)**. ([Semantic Scholar][2])

For MFH152, the parental population already contained heterogeneous ploidies. Cells were sorted according to **FSC/SSC size and granularity**, with single large and small cells cloned. The authors state that they sorted **1000 “small” + 1000 “big” clones**, that the majority of surviving clones were actually aneuploid, and that they isolated stable diploid/tetraploid clones from the survivors. They use **3 diploid and 4 tetraploid MFH152 clones**, preferentially at early passage. ([Semantic Scholar][2])

So, for provenance:

> **RKO:** experimentally induced cytokinesis failure → clonally derived stable near-diploid and near-tetraploid populations from the same parental line.

But it is still **not a single clone measured immediately before and after WGD**: D1/D2 and T1/T2 are independently isolated clones, so clone-to-clone variation is present alongside ploidy.

---

# Figure 1 — characterization of RKO near-diploid and near-tetraploid clones

## 1A — forward-scatter distributions

**Plot:** flow-cytometry histogram of FSC, used as a proxy for cell size.

The near-tetraploid distribution is visibly shifted rightward relative to the near-diploid distribution. No individual numerical values are printed on this histogram. 

## 1B — mean FSC

**Bar plot**, y-axis 0–500 arbitrary units.

**F:**

| Clone |          FSC |
| ----- | -----------: |
| D1    | **~275–280** |
| D2    | **~285–290** |
| T1    | **~330–335** |
| T2    | **~345–350** |

T1 and T2 are each marked ******* relative to D1. Caption: data are means ± SEM; **n=5**. 

Thus the tetraploid clones are only about **~20–25% larger by FSC**, rather than displaying the enormous size difference seen in Xuan's PGCC populations.

---

## 1C — cell/nuclear morphology

Representative confocal images:

* phalloidin = actin;
* DAPI = DNA;
* scale bar **10 μm**;
* one nucleus outlined white;
* one cell outlined yellow.

The shown tetraploid cell is visibly larger in both cell and nuclear cross-section. This panel itself is qualitative. 

---

## 1D — nuclear area

**Scatter plot:** nucleus size, y-axis 0–200 μm²; individual nuclei plus mean.

**F approximate mean values:**

| Clone | Mean nuclear area |
| ----- | ----------------: |
| D1    |    **~32–34 μm²** |
| D2    |    **~31–33 μm²** |
| T1    |    **~41–43 μm²** |
| T2    |    **~44–46 μm²** |

T1 and T2 are each *** vs D1. 

The distributions are worth recording because the means do not capture them well:

* D1/D2: most nuclei approximately **20–45 μm²**, occasional values around **60–70 μm²**.
* T1: most roughly **25–65 μm²**, with several values extending to approximately **80–90 μm²**.
* T2: most ~30–60 μm², with several large nuclei around **80–100 μm²** and an extreme point around **150–160 μm²**.

Thus the tetraploid distributions overlap substantially with diploids but have larger means and a pronounced upper tail. 

---

## 1E–F — chromosome counts

Metaphase spreads are shown in 1E; 1F quantifies chromosome number.

The numbers of metaphases are printed:

* D1 **n=20**
* D2 **n=21**
* T1 **n=20**
* T2 **n=20**

**F approximate means:**

| Clone | Chromosomes |
| ----- | ----------: |
| D1    |  **~42–44** |
| D2    |  **~42–44** |
| T1    |  **~93–95** |
| T2    |  **~90–92** |

The individual values are extremely tight around those means compared with the nuclear-area distributions. 

This directly validates that T1/T2 are genuinely **near-tetraploid** rather than merely large/polyploid cells of uncertain DNA content.

---

## 1G–H — DNA-content/cell-cycle distributions

1G shows representative PI histograms, with:

* diploid G1 around **2n**, G2/M around **4n**;
* tetraploid G1 around **4n**, G2/M around **8n**.

1H shows percentage of cells in G1, S, and G2/M.

**F:**

| Clone |          G1 |       S |    G2/M |
| ----- | ----------: | ------: | ------: |
| D1    |     ~51–52% |    ~21% | ~27–28% |
| D2    |     ~52–53% | ~23–24% | ~23–25% |
| T1    | **~63–64%** | ~26–27% |  ~9–10% |
| T2    | **~64–66%** | ~24–25% |  ~9–10% |

The G1 fractions of T1/T2 carry *** markers. Caption gives **n=5**. 

The authors subsequently use this G1 accumulation as one indication that the tetraploid clones proliferate more slowly. 

---

# Figure 2 — proliferation

This experiment is important because both wound healing and Oris assays run for **48 h**, giving proliferation enough time to influence area closure.

Cells were seeded at **2,000 cells/well in 96-well plates**, cultured for up to **5 days**, fixed at each time point, crystal-violet stained, and quantified by absorbance at **595 nm**. The plotted quantity is a normalized **“proliferation index.”** ([Semantic Scholar][2])

## 2A — RKO

**Line plot**, 0–120 h.

**F approximate values:**

|  Time |             D1 |             D2 |           T1 |           T2 |
| ----: | -------------: | -------------: | -----------: | -----------: |
|   0 h |              1 |              1 |            1 |            1 |
|  24 h |           ~2.0 |       ~1.8–1.9 |         ~1.5 |     ~1.5–1.6 |
|  48 h |           ~4.0 |       ~3.4–3.6 |         ~2.7 |     ~2.4–2.5 |
|  72 h |       ~5.7–5.9 |           ~5.0 |         ~3.7 |         ~3.9 |
|  96 h |       ~9.2–9.5 |       ~7.8–8.0 |         ~5.2 |         ~5.2 |
| 120 h | **~14.0–14.3** | **~13.7–14.0** | **~9.4–9.6** | **~8.8–9.1** |

At the endpoint both tetraploid clones are substantially below both diploid clones. T1/T2 are marked *** relative to D1. **n=5.** 

At 120 h, therefore, the tetraploid proliferation index is roughly **~65–68% of D1**.

---

## 2B — MFH152

The same assay is performed across 3 diploid and 4 tetraploid MFH152 clones.

At 120 h, visually:

* diploid clones: approximately **~5.2–5.7**
* tetraploid clones: approximately **~3.0–3.6**

Again, all tetraploid clones are below the diploid group, with *** indicated. 

The migration advantage described below therefore occurs despite a measurable **proliferation disadvantage** in the same clone system.

---

# Figure 3 — RKO migration

This is the most important figure for our specific question because it contains **three separate migration assays**, including true single-cell time-lapse tracking.

---

## Figure 3A–B — scratch/wound-healing assay

### Assay conditions

The cells are grown as **confluent monolayers**.

A scratch is generated using a sterile **200-μL pipette tip** and monitored for **up to 48 h**.

Because initial wound widths vary, each measurement is normalized as:

$$
\mathrm{Remaining\ free\ area}(t)/\mathrm{free\ area}(0)
$$

so time 0 is **100%** in every clone. ([Semantic Scholar][2])

No numerical seeding density is given for this assay; the relevant density information is simply that cells were explicitly **confluent** before scratching.

### 3A — representative fields

Images at 0 and 48 h, scale **1 mm**.

The initial scratch is visibly wide relative to individual cells and the monolayers appear densely confluent. Yellow dashed lines mark the free area.

By 48 h, both tetraploid wounds appear substantially narrower than either diploid wound. 

### 3B — remaining cell-free area

The figure visually overlays the 0-, 24-, and 48-h values in each bar:

* white boundary = 0 h = 100%;
* gray boundary = 24 h;
* black boundary = 48 h.

**F:**

| Clone | 24 h remaining area | 48 h remaining area |
| ----- | ------------------: | ------------------: |
| D1    |            **~55%** |         **~31–32%** |
| D2    |         **~69–70%** |         **~46–48%** |
| T1    |         **~35–36%** |           **~3–5%** |
| T2    |         **~29–31%** |         **~13–15%** |

T1 and T2 are marked *** compared with D1. Data are means ± SEM; **n=5**. 

This is a large phenotype. At 48 h:

* D1 still has roughly one-third of the original wound open;
* D2 nearly one-half;
* T1 is essentially closed;
* T2 retains only ~14%.

Because tetraploids proliferate **more slowly**, the increased wound closure cannot be explained by faster population growth in these clones.

---

# Figure 3C–D — Oris™ 2D migration assay

This provides a second collective-migration geometry without mechanically scratching cells.

### Assay conditions

Cells were plated at **4 × 10⁴ cells/well** in a **96-well Oris plate** containing a silicone stopper and grown overnight. Removing the stopper creates an initially cell-free circular region; cells then migrate into that region for up to 48 h. ([Semantic Scholar][2])

The Results and figure caption describe the cells as **confluent** when the stopper is removed. ([PubMed Central (PMC)][1])

### 3C — representative images

Images show the circular migration zone at 0 and 48 h; scale **1 mm**.

The T1 example shows visibly greater occupation of the originally empty circle than the D1 example. 

### 3D — remaining cell-free area

**F:**

| Clone |        24 h |        48 h |
| ----- | ----------: | ----------: |
| D1    | **~93–95%** | **~57–60%** |
| D2    | **~97–99%** | **~64–67%** |
| T1    | **~84–86%** | **~49–51%** |
| T2    | **~81–83%** | **~33–36%** |

T1 is marked ** and T2 *** relative to D1. 

The phenotype is less dramatic than in the scratch assay but has the same direction in both independent tetraploid clones.

T2 is particularly strong: by 48 h its remaining open area is approximately **half that of D2**.

---

# Figure 3E–F — individual-cell migration

This is the part of the paper that most directly addresses migration independent of collective monolayer dynamics.

### Assay conditions

Cells were seeded at **low concentration specifically to avoid cell-cell contacts**.

The paper does **not give a numerical plating density** for this experiment.

Cells were imaged by time-lapse microscopy:

* **every 1 h**
* for **24 h**
* using an inverted Zeiss Axio Observer Z1.

Trajectories were reconstructed in ImageJ. ([Semantic Scholar][2])

So temporal resolution is much coarser than Xuan 2018/2020: **24 one-hour steps per trajectory**, rather than 10-min imaging.

### 3E — example trajectories

Representative phase images are shown with colored tracks overlaid; scale bar **50 μm**.

The tetraploid fields visibly contain several much longer individual tracks than the diploid fields. 

The raw track geometry is shown only illustratively; the paper does **not quantitatively extract**:

* displacement;
* path length;
* straightness;
* turning angle;
* persistence time;
* MSD;
* arrest fraction.

Only **average cell speed** is quantified.

### 3F — individual-cell speed

**Scatter plot**, y-axis **individual cell speed (μm/hour)**, range 0–40.

Each dot appears to represent an individual tracked cell; a horizontal bar shows the central/mean value. Caption gives **n=5**, apparently referring to experimental replicates rather than the number of displayed tracks. The total number of tracked cells is not explicitly stated. 

**F approximate central values:**

| Clone | Mean individual speed |
| ----- | --------------------: |
| D1    |       **~2.5–3 μm/h** |
| D2    |           **~3 μm/h** |
| T1    |       **~7–7.5 μm/h** |
| T2    |       **~6.5–7 μm/h** |

So both tetraploid clone means are approximately **2.3–2.8× D1**.

T1 and T2 are both *** relative to D1. 

The distribution itself contains additional information:

**D1**

* majority roughly **1–4 μm/h**
* upper observations ~5–6 μm/h.

**D2**

* majority roughly **1–5 μm/h**
* several ~6–9 μm/h.

**T1**

* dense bulk approximately **4–10 μm/h**
* several tracks approximately **10–16 μm/h**
* one conspicuous outlier at about **34 μm/h**.

**T2**

* bulk approximately **3–10 μm/h**
* several tracks ~12–18 μm/h
* two high observations around **24–26 μm/h**.

Thus the difference is not only a small shift of the entire speed distribution: the tetraploid populations also contain a **long high-speed tail** that is almost absent from the diploid clones. 

The paper does not report whether these speed values are:

* total path length / 24 h,
* average of frame-to-frame speeds,
* or some other exact ImageJ tracking-derived definition.

The Methods simply state that the time lapses were analyzed in ImageJ. ([Semantic Scholar][2])

---

# Figure 4 — transmembrane migration / “invasion”

## Figure 4A–B — Boyden chamber assay

### Assay conditions

For RKO:

* **1,000,000 cells** added to the upper chamber;
* upper medium is **serum-free**;
* lower compartment contains growth medium + **10% FBS** as chemoattractant;
* assay duration **24 h**;
* cells reaching the lower side of the filter are fixed with 4% PFA, crystal-violet stained and counted. ([Semantic Scholar][2])

Figure 4A shows representative lower-filter fields, scale **100 μm**.

### 4B — migrated cell count

Y-axis is **cell count (10³)**.

**F:**

| Clone |    Migrated cells |
| ----- | ----------------: |
| D1    | **~1.0–1.1 ×10³** |
| D2    | **~0.9–1.0 ×10³** |
| T1    |     **~3.1 ×10³** |
| T2    | **~3.2–3.3 ×10³** |

Approximate SEMs:

* D1 ~0.45–0.5×10³
* D2 ~0.4×10³
* T1 ~0.6×10³
* T2 ~0.7×10³.

T1/T2 are both ***. Caption: **n=3**. 

Thus approximately **3× as many tetraploid RKO cells** reach the underside of the membrane after 24 h.

### Important assay-definition qualification

The paper repeatedly describes this as an **“invasion”** assay. However, the printed Methods do **not describe coating the Boyden membrane with Matrigel, collagen, or another extracellular matrix barrier**. The documented protocol therefore establishes serum-directed **transmembrane migration**; whether an unreported matrix coating was present cannot be determined from the method as published. ([Semantic Scholar][2])

---

# Figure 4C — real-time xCELLigence assay

This experiment compares only:

* **D1**
* **T1**

rather than all four RKO clones.

### Assay conditions

CIM-16 plates are used.

* **2 × 10⁴ cells/well**
* serum-free upper compartment;
* complete medium + **10% serum** below;
* impedance measured continuously;
* impedance is converted to a dimensionless **cell index** proportional to the number of cells reaching the electrode-containing underside;
* monitored for **24 h**; caption says mean impedance measured for **26 h**.
* **n=3**. ([Semantic Scholar][2])

### Quantitative trajectory

The T1 curve separates from D1 very early and continues to rise.

**F approximate values:**

|     Time |   D1 cell index |  T1 cell index |
| -------: | --------------: | -------------: |
|      0 h |              ~0 |             ~0 |
|     ~6 h | ~−0.04 to −0.05 | **~0.04–0.05** |
|    ~12 h |          ~−0.05 | **~0.10–0.12** |
|    ~18 h |          ~−0.05 | **~0.27–0.30** |
|    ~21 h |          ~−0.05 |      **~0.40** |
| ~24–25 h |          ~−0.05 | **~0.55–0.60** |

The figure marks the overall comparison ***. 

The T1 curve also develops visibly increasing SEM at later times, reaching roughly ±0.05–0.07 cell-index units by ~24 h.

Again, although the paper calls this a **3D invasion assay**, the printed Methods do not describe coating the CIM membrane with an ECM gel. The documented measurement is impedance generated by cells crossing a membrane toward serum. ([Semantic Scholar][2])

---

# Figure 5A–B — in-vivo lung colonization after tail-vein injection

This experiment switches from RKO to **MFH152 sarcoma clones**.

### Assay

* NSG mice, **8 weeks old**
* **n=12 total**
* **6 mice/group**
* **1 × 10⁶ cells/mouse**
* intravenous tail-vein injection
* diploid or tetraploid MFH152 clone
* animals sacrificed after **8 weeks**
* lungs fixed and histologically assessed. ([Semantic Scholar][2])

### 5A — histology

The diploid example shows approximately normal lung parenchyma.

The tetraploid example shows a region described histologically as:

* cellular proliferation within lung parenchyma;
* dense cellular infiltration;
* cellular atypia.

Several magnification levels are shown, with **200, 50, 10 and 20 μm** scale bars. 

### 5B — incidence

**P:**

* diploid-injected: **1/6 mice**
* tetraploid-injected: **4/6 mice**

Thus observed lesion incidence is approximately:

* **16.7% diploid**
* **66.7% tetraploid**.

No significance test or confidence interval is shown for this 1/6 versus 4/6 comparison. 

Also, unlike the RKO in-vitro work, this experiment appears to compare **one selected diploid MFH152 clone against one selected tetraploid MFH152 clone**, not all 3/4 MFH clones.

---

# Figure 5C–D — human sarcoma karyotype data

These panels use the **Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer**, not the experimental RKO/MFH cultures. ([Semantic Scholar][2])

## 5C — matched primary versus metastatic sarcoma karyotypes

**Violin plot**, approximate ploidy on a log-spaced axis labelled 1n, 2n, 4n, 8n.

The figure states:

* **n=61**
* **p=0.0035**.

The text, however, says they analyzed **63 matched pairs**. This is an internal sample-size discrepancy worth retaining. 

Visually:

**Primary tumors**

* strongly concentrated around **~2n**
* most values roughly **1.6–3n**
* upper values close to ~4n.

**Metastases**

* still contain many ~2n cases;
* much broader upper distribution;
* numerous cases around **3–4n**;
* some approaching approximately **5–6n**.

Importantly, the caption states that the statistical comparison is a **Fisher test of variances**, not a test of mean/median ploidy. The reported **p=0.0035 therefore directly supports a difference in distributional variance**, not by itself a higher average ploidy. 

The authors' text interprets the broader metastatic distribution as enrichment for **near-triploid and near-tetraploid** karyotypes. ([Semantic Scholar][2])

---

## 5D — primary → first metastasis → second metastasis

A smaller dataset, **n=12**, compares three sequential sites:

* P = primary
* M1 = first metastatic site
* M2 = second metastatic site.

Y-axis is labelled **Log2 number of chromosome per clone**, with reference lines corresponding approximately to:

* **2n ≈46 chromosomes**
* **4n ≈92 chromosomes**.

**F distributions:**

**Primary**

* roughly **35–80 chromosomes**
* strongly centered near ~45–50.

**M1**

* roughly **38–85**
* central distribution still approximately 45–55.

**M2**

* extremely broad;
* values from roughly **20–110 chromosomes**
* several values around or above the **4n ~92-chromosome line**.

Figure reports **p=0.0196**, again using a **Fisher variance test**. 

So the strongest directly demonstrated feature is the **increased karyotypic spread/variance**, particularly at the second metastatic site.

---

# Supplementary figures explicitly referenced by the main article

The publisher/PMC record lists a **387.2-kB supplementary ZIP**. ([PubMed Central (PMC)][1])

The main article identifies at least:

**Figure S1A:** FSC/SSC gating used to isolate small versus large MFH152 single-cell clones.

**Figure S1B:** Hoechst/DNA-content gating strategy used after cytochalasin-D treatment to isolate RKO diploid versus tetraploid clones.

**Figure S2A:** MFH152 cell-size/FSC comparison; tetraploid clones are larger.

**Figure S2B:** MFH152 nuclear-area comparison; tetraploid clones have larger nuclei.

**Figure S2C:** MFH152 metaphase chromosome counts; diploid clones are approximately diploid and tetraploid clones approximately doubled.

**Figure S2D:** MFH152 cell-cycle distributions, again showing greater G1 accumulation in tetraploid clones. ([Semantic Scholar][2])

I have **not numerically figure-read those supplementary panels here**, because the ZIP itself is not rendering into inspectable figures through the browser interface in this pass. I therefore wouldn't fabricate approximate values from S1/S2. The main-paper RKO figures above have been directly inspected.

---

# Important provenance / measurement qualifications

### This is much closer to a WGD experiment than Xuan 2018/2020

The RKO tetraploids were explicitly generated by **cytokinesis failure**, then chromosome-counted at approximately 90–95 chromosomes. That makes the ploidy manipulation substantially cleaner than classifying naturally occurring giant cells by size/DNA content. ([Semantic Scholar][2])

At the same time, D1/D2/T1/T2 are **independently selected clonal derivatives**. The design therefore compares replicated diploid versus tetraploid clones rather than observing the exact same clone immediately before and after WGD.

### The single-cell time series is real but temporally coarse

Individual cells are imaged **once per hour for 24 h**. ([Semantic Scholar][2])

That is sufficient for estimating gross motility/speed but poorly suited to resolving:

* protrusion/retraction cycles;
* short-lived pauses;
* rapid turning;
* polarity reversals;
* front competition;
* sub-hour migration states.

The paper does not calculate any of those quantities.

### Cell density differs radically among the three Figure 3 assays

* Scratch: **confluent monolayer**
* Oris: **4×10⁴ cells/well**, grown overnight, described as confluent
* Single-cell tracking: deliberately **low density / non-contacting**

The tetraploid migration advantage is observed in all three contexts. ([Semantic Scholar][2])

### Proliferation acts in the opposite direction

Both RKO tetraploid clones grow substantially more slowly than diploid clones in the five-day crystal-violet assay, yet they close wounds and Oris gaps more rapidly. 

### “Invasion” is not fully documented as ECM invasion

For both the Boyden and xCELLigence experiments, the Methods describe cells crossing porous membranes toward serum but do **not mention an ECM/Matrigel coating**. ([Semantic Scholar][2])

Therefore the directly documented quantity is transmembrane chemotactic migration. The paper's terminology is “invasion.”

### The migration analysis is overwhelmingly velocity/endpoint based

Despite showing individual trajectories, the paper does not report:

* persistence;
* net displacement;
* accumulated path length;
* turning angles;
* MSD;
* autocorrelation;
* morphology;
* aspect ratio;
* protrusions;
* Rac/Rho activity;
* migration-state occupancy.

The most detailed single-cell quantitative result is simply **speed**, approximately ~3 μm/h in diploids versus ~7 μm/h in tetraploids. 

---

# Authors' interpretations and hypotheses — separate from the measured data

The authors' primary interpretation is that **tetraploidization confers enhanced migratory, invasive and metastatic capacity despite a proliferation disadvantage**. They regard the consistency of the scratch, Oris, individual-cell, transwell, xCELLigence and mouse experiments as support for this overall conclusion. ([PubMed Central (PMC)][1])

They do **not identify a mechanism experimentally in this paper**.

In the Discussion, they instead cite several hypotheses from prior literature:

* supernumerary centrosomes might promote invasion;
* polyploid cells may exhibit enhanced EMT programs;
* polyploidy has been associated elsewhere with SNAI, TWIST and N-cadherin;
* expression changes in migration/invasion-related genes may contribute;
* metabolic reprogramming may contribute;
* reduced expression of membrane-associated genes could facilitate detachment;
* polyploidy-associated inflammatory signaling may maintain a pro-metastatic state. ([Semantic Scholar][2])

Those are **literature-based proposed mechanisms**, not mechanistic results from the Jemaà 2023 experiments.

---

## Migration-relevant quantitative core

If we eventually reduce this paper to the handful of quantities most useful for the WGD migration comparison, I would preserve the following as the raw record:

| Measurement                |         Diploid |        Tetraploid | Context                        |
| -------------------------- | --------------: | ----------------: | ------------------------------ |
| FSC size                   |        ~280–290 |          ~330–350 | RKO clones                     |
| Nuclear area               |      ~32–33 μm² |        ~42–45 μm² | RKO                            |
| Chromosome number          |             ~43 |            ~91–94 | RKO                            |
| Proliferation index, 120 h |             ~14 |            ~9–9.5 | RKO                            |
| Scratch free area, 48 h    |      **31–48%** |        **~4–14%** | confluent collective migration |
| Oris free area, 48 h       |     **~58–66%** |       **~35–50%** | dense collective migration     |
| Individual speed           | **~2.5–3 μm/h** | **~6.5–7.5 μm/h** | low-density, 1-h cadence, 24 h |
| Boyden migrated cells      |      **~1×10³** |      **~3.2×10³** | 24 h, serum gradient           |
| xCELLigence index at 24 h  |          ~−0.05 |    **~0.55–0.60** | D1 vs T1 only                  |
| Lung lesion incidence      |         **1/6** |           **4/6** | MFH152 tail-vein model         |


[1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10531202/?utm_source=chatgpt.com "Tetraploidization Increases the Motility and Invasiveness of Cancer Cells - PMC"
[2]: https://pdfs.semanticscholar.org/9f68/913f4aaa76af8b05e1243711f659f28d027e.pdf "Tetraploidization Increases the Motility and Invasiveness of Cancer Cells"
