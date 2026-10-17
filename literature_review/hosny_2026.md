### Tracking / observation protocol

PC3-GFP cells were imaged for **16 h at 30-min intervals**. Tiled GFP and oxygen-sensor images were acquired with a Nikon Eclipse TE2000-E, **4× objective**, Andor Zyla camera; GFP exposure was 300 ms and mCherry oxygen-sensor exposure 500 ms. Segmentation/tracking used **Cellpose-SAM in TrackMate/Fiji**, yielding cell position and area at each frame. 

Live-cell PACCs were defined retrospectively as cells with **area ≥1500 μm² at 16 h and ≥3-fold area increase from 0–16 h**. Cells failing either criterion were classified as non-PACCs. Normoxic PACCs were excluded. 

Fixed-cell calibration gave mean areas of **487 ± 3 μm²** for non-PACCs and **1490 ± 48 μm²** for PACCs (mean ± SEM). 

Cells analyzed for migration were selected from **regions of the steepest O₂ gradient**. The system generated a radial field with approximately **\(2\times10^{-3}\) bar O₂** at the core and **0.2 bar O₂** at the oxygen-rich periphery.  

### Fig. 3 — baseline motility

Definitions: **Accumulated Distance** = sum of traveled distance; **Euclidean Distance** = net displacement; **Directness** = Euclidean/Accumulated Distance, ranging from 0 to 1. 

\(n=60\) cells/group: 20 cells from each of 3 experiments. Data are median ± IQR. 

Approximate direct figure reads:

| Condition        | Accumulated distance (μm) | Euclidean distance (μm) |             Directness |
| ---------------- | ------------------------: | ----------------------: | ---------------------: |
| Hypoxic PACC     |       **~175 [140, 223]** |       **~82 [60, 124]** | **~0.51 [0.37, 0.74]** |
| Hypoxic non-PACC |        **~105 [85, 150]** |        **~55 [34, 95]** | **~0.31 [0.26, 0.49]** |
| Normoxic control |       **~210 [170, 285]** |        **~42 [20, 67]** | **~0.20 [0.12, 0.26]** |

Equivalent accumulated-distance/16-h values are approximately **10.9, 6.6, and 13.1 μm/h**, respectively. At the 30-min sampling interval, these correspond to mean accumulated increments of approximately **5.5, 3.3, and 6.6 μm/frame**.

### Fig. 4 — directional migration / aerotaxis

For each trajectory, the oxygen-gradient direction was approximated as the radial vector from the system center to the **midpoint between the trajectory's initial and final positions**. Tracks were rotated so this local gradient direction aligned with \(+\hat y\). 

The Rayleigh test was applied to each trajectory's **final displacement angle** after this rotation. Exact reported \(p\)-values:

| Condition        | Rayleigh \(p\) |
| ---------------- | -------------: |
| Normoxic control |      **0.673** |
| PACC             |     **0.0254** |
| non-PACC         |     **0.0390** |



Approximate Fig. 4C/D reads:

| Condition |            Parallel FMI |       Perpendicular FMI |
| --------- | ----------------------: | ----------------------: |
| PACC      |  **~0.25 [0.05, 0.50]** | **~0.03 [-0.33, 0.21]** |
| non-PACC  | **~0.10 [-0.17, 0.39]** | **~0.02 [-0.21, 0.31]** |
| Normoxia  | **~0.00 [-0.29, 0.17]** | **~0.00 [-0.18, 0.20]** |

\(n=60\) cells/group; median ± IQR. Fig. 4E/F instead reports mean ± SEM by radial bin. 

Approximate Fig. 4F radial profile:

| Radius from center (μm) | Mean parallel FMI ± SEM |
| ----------------------: | ----------------------: |
|                  ~2,000 |      **~0.145 ± 0.065** |
|                  ~4,000 |      **~0.075 ± 0.070** |
|                  ~4,800 |      **~0.100 ± 0.062** |
|                  ~5,600 |      **~0.108 ± 0.067** |
|                  ~6,400 |      **~0.350 ± 0.052** |
|                  ~7,500 |      **~0.312 ± 0.066** |
|                  ~8,650 |      **~0.085 ± 0.073** |
|                 ~11,400 |      **~0.105 ± 0.074** |

The authors report that time-resolved **per-step parallel FMI** increased during the hypoxia time course in PACCs and declined again at later times; confidence bands were broad because of the smaller PACC sample. No numerical S6 values are available in the attached PDF. 

One reporting detail: the printed FMI equation in Methods is parsed as a projection of net displacement onto the gradient axes, whereas the plotted FMI values are dimensionless and approximately bounded by ±1. The attached PDF does not explicitly resolve that normalization detail. 

### Fig. 6 — HIF-1α knockdown

Cells were transfected with non-targeting control siRNA or HIF-1α siRNA and tracked under hypoxia for 16 h. \(n=60\) cells/group. Distance/directness and FMI are median ± IQR. 

Approximate direct figure reads:

| Metric                    |            PACC siNTC |             PACC siHIF |         non-PACC siNTC |         non-PACC siHIF |
| ------------------------- | --------------------: | ---------------------: | ---------------------: | ---------------------: |
| Accumulated distance (μm) |    **~198 [154,258]** |     **~166 [121,207]** |      **~121 [86,175]** |       **~81 [52,141]** |
| Euclidean distance (μm)   |      **~87 [62,120]** |        **~54 [43,92]** |        **~60 [39,78]** |        **~30 [17,50]** |
| Directness                | **~0.49 [0.33,0.62]** |  **~0.32 [0.21,0.49]** |  **~0.36 [0.25,0.50]** |  **~0.26 [0.13,0.36]** |
| Parallel FMI              | **~0.21 [0.00,0.49]** | **~0.06 [-0.14,0.25]** | **~0.13 [-0.06,0.30]** | **~0.04 [-0.14,0.17]** |
| Perpendicular FMI         |                 ~0.03 |                 ~−0.03 |                  ~0.04 |                  ~0.02 |
| Rayleigh \(p\)            |           **0.00352** |              **0.114** |             **0.0152** |              **0.416** |

Accumulated-distance/16-h values are approximately **12.4 → 10.4 μm/h** for PACCs and **7.6 → 5.1 μm/h** for non-PACCs.

The authors report significant reductions in accumulated distance, Euclidean distance, directness, and parallel FMI after HIF-1α depletion; perpendicular FMI remained approximately zero. 

### Fig. 7 — RhoA knockdown

Same 16-h tracking design; \(n=60\) cells/group. Distance/directness and FMI are median ± IQR. 

Approximate direct figure reads:

| Metric                    |            PACC siNTC |            PACC siRhoA |         non-PACC siNTC |        non-PACC siRhoA |
| ------------------------- | --------------------: | ---------------------: | ---------------------: | ---------------------: |
| Accumulated distance (μm) |    **~182 [148,240]** |     **~164 [117,210]** |     **~132 [111,215]** |       **~94 [75,133]** |
| Euclidean distance (μm)   |      **~87 [52,115]** |        **~59 [32,93]** |        **~63 [39,85]** |        **~38 [16,57]** |
| Directness                | **~0.51 [0.27,0.65]** |  **~0.28 [0.18,0.45]** |  **~0.38 [0.22,0.52]** |  **~0.25 [0.15,0.42]** |
| Parallel FMI              | **~0.22 [0.00,0.36]** | **~0.04 [-0.23,0.19]** | **~0.13 [-0.08,0.27]** | **~0.03 [-0.24,0.15]** |
| Perpendicular FMI         |                    ~0 |                     ~0 |                     ~0 |                     ~0 |
| Rayleigh \(p\)            |            **0.0264** |              **0.286** |             **0.0351** |              **0.457** |

Accumulated-distance/16-h values are approximately **11.4 → 10.2 μm/h** for PACCs and **8.2 → 5.9 μm/h** for non-PACCs.

The authors report reduced accumulated and Euclidean distances, reduced directness, loss of significant directional bias by Rayleigh test, and reduced parallel FMI after RhoA knockdown; perpendicular FMI remained near zero. 

### Additional trajectory-associated observations

Across PACC and non-PACC populations, the authors report **no significant Spearman correlation** between cell area and accumulated distance, Euclidean distance, directness, parallel FMI, or perpendicular FMI. Exact coefficients are in supplementary figures not present in the attached PDF. 

For the supplementary subpopulation analysis, motility classes were defined from the **control-PACC Euclidean-distance distribution**: high motility = above the 75th percentile, low motility = below the 25th percentile, intermediate = between those thresholds. Cells were additionally divided by the sign of parallel FMI. HIF-1α and RhoA knockdown reduced the fraction of positive-FMI/high-motility PACCs and redistributed cells toward lower-motility and negative-FMI classes. 

The attached PDF does **not** report numerical MSD-versus-lag curves, velocity autocorrelation functions, turning-angle distributions, step-length distributions, or explicit persistence times.
