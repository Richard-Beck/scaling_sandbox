# Figure R2. Process-dependent size effects on activity-heading persistence

The activity heading is the direction of the first spatial Fourier moment of activity around the ring. Slower heading-correlation decay means greater persistence. This simulated activity readout is distinct from measured cell migration and conditional OU velocity-persistence estimates.

Per-perimeter supply means a fixed event density per unit circumference; per-cell supply means a fixed arrival intensity per cell. Goryachev is a biochemical polarity model; Spring couples polarity chemistry to a deforming domain. Model events supply localized, exponentially decaying activation pulses at uniform material positions.

Each row pairs an illustrative activity history (A, C, E, G, I, K) with a heading-correlation plot (B, D, F, H, J, L). The horizontal rule separates nulls from models. Circumference C is labelled C20/C60, or native resting size C1/C3 for Spring. Orange solid curves denote the smaller size; blue dashed curves the larger. Kymograph time increases upward. Physical arclength, time, lag and activity retain their native units: nulls use their defined length/time and unit peak height; each model uses its own uncalibrated native length, time and active-concentration units. Values across models are not calibrated against one another.

Kymograph widths follow circumference; Spring's edge follows its instantaneous circumference, with white space outside the domain. Within each row, both sizes share a time window and a viridis activity scale. The relocation scale is 0–1; all four Spring facets share one scale. White lines in G, I and K place the first-moment direction at its equivalent arclength, which need not follow a local activity maximum, where first-mode coherence is at least 0.1; lines break at the periodic seam or excluded frames. No profile is recentered. Finite Fourier reconstructions in G and I are **clipped at zero for display only**; source values and heading calculations remain signed. K uses full native profiles. Fixed first-replicate model examples start at the first post-burn sample and span 40 event decay times (displayed window lengths: 2,000 Goryachev; 200 Spring time units), uniformly sampled at 801 times. Ensemble correlations use all four original histories per size and the full common half-record lag range, independently of those short examples.

**A–B, uniform relocation with constant success:** enlargement reduces persistence. Relocations choose a uniform heading; per-perimeter supply gives

$$
\lambda(C)=\lambda_0\frac{C}{40},\qquad \lambda_0=\frac{1}{7000},\qquad
R(\ell)=\exp[-\lambda(C)\ell]. \tag{1}
$$

A reuses round-7 paths over 24,000 time units at interval 20, displaying a Gaussian peak of SD 1.5. B is the closed form.

**C–D, uniform relocation with success proportional to inverse circumference squared:** enlargement increases persistence. Per-perimeter supply with success proportional to $1/C^2$ gives $\lambda(C)=\lambda_0(40/C)$; D shows equation (1) with this rate. C reuses the existing fixed-seed R2 paths over the same duration and cadence as A. Retained relocation paths use the original finite-step Bernoulli arrivals, whereas B and D show the continuous-time closed forms.

**E–F, patch field with uniform births:** the many-patch prediction is size-independent. Gaussian patches (SD 1.5, unit peak) are born at uniform random positions in a Poisson process of density $b$ per unit length per unit time. Each remains fixed for an independent exponential lifetime of mean $1/r$; activity is the sum of live patches. At both sizes $r=1/7000$ and $b=r/6$. E is a new stationary sample path. F shows

$$
R(\ell)=\frac{\pi}{4}\rho\,{}_2F_1\!\left(\frac12,\frac12;2;\rho^2\right),
\qquad \rho=e^{-r\ell}. \tag{2}
$$

Orange solid and blue dashed curves coincide. This **many-patch approximation** is particularly approximate at C20 (mean 3.33 live patches, versus 10 at C60), rather than a fit to E. A finite path can lose all patches and have no instantaneous heading.

**G–H, Goryachev with per-perimeter events:** enlargement increases persistence. The existing low-rate/high-amplitude campaign has arrival density $10^{-4}$ per unit arclength per unit time, mean exponential event amplitude 19.2 (three times reference 6.4), Gaussian event SD $1.5/\sqrt{2}$, and exponential event decay time 50. G uses all retained active Fourier modes.

**I–J, Spring with per-perimeter events:** enlargement reduces persistence. The existing low arrival density is 0.01 per unit instantaneous circumference per unit time, mean exponential amplitude 12 (three times reference 4), Gaussian SD $0.15/\sqrt{2}$, and exponential event decay time 5. I uses all retained active Fourier modes.

**K–L, Spring with per-cell events:** enlargement increases persistence. Both outer sizes have constant intensity 0.0204418072951 per cell per unit time, fixed from the exact size-2 basal circumference. Other event parameters and local dynamics match I–J. The original arm saved only modes 0/1, so K uses one full-profile replay per size with unchanged seeds and event sequences; L uses the original four histories for inference.

H, J and L show the pooled mean $\langle\cos[\theta(t+\ell)-\theta(t)]\rangle$, with both endpoints satisfying $|P|\geq0.1$, $P=\overline{M_1}/M_0$ for the stored Fourier convention. Grey dotted lines mark correlation 0.5. Bands are pointwise 95% whole-seed bootstrap intervals (n = 4 histories per size; 2,000 resamples of full histories). These lag axes are logarithmic, spanning 64,000 Goryachev or 6,400 Spring units; exact zero-lag correlation 1 is omitted. B, D and F use linear lag axes with no uncertainty bands.

**Separate mode-1 analysis; these are not fits to H, J or L.** The unchanged campaign exponent α is the negative three-size log-rate/log-mean-circumference slope, using first mode-1 half-height crossings:

| Model / supply | α | 95% whole-seed interval |
| --- | ---: | --- |
| Goryachev / per-perimeter | ≥ 5.1654 | Lower endpoint 5.0579; upper unbounded |
| Spring / per-perimeter | −0.2389 | [−0.3038, −0.1788] |
| Spring / per-cell | 0.8695 | [0.8396, 0.8945] |

Goryachev C60 remains censored. The per-cell exponent shares the original size-2 trajectory with its circumference-dependent hazard. [Source records and methods](../../../README.md) document selection, replay checks and native conventions.
