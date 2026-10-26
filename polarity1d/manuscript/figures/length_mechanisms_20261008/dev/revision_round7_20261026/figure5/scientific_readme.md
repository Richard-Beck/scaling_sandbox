# Unexpected

At strong input in Goryachev C60, the median |P| is 0.047. The raw phase is therefore poorly localized even though it is mathematically defined; its broad displacement density cannot be identified uniquely as relocation of a coherent dominant domain. The distributions are shown without a cut or a pass criterion. The new panel does not by itself prove that Goryachev heading loss is exclusively jump-mediated.

## Raw heading-change method

Every retained compact seed for each displayed condition enters (96–384 seeds). Heading is arg(P). Changes are principal signed differences, so the readout is the shortest wrapped displacement over a lag, not cumulative path length. Physical circumference is averaged between endpoints. There is no |P|, contrast, region-count or track-identity mask. Only nonfinite numerical differences would be excluded, and their counts are exported. Variable physical-width bins are symmetric asinh spaced with scale C×10⁻⁴; centers and bounds are in D_raw_heading_distributions.csv.gz. Density normalization uses all finite pairs and the exact bin widths. Rescaling by √lag changes both coordinates and density. The linear signed x axis makes the core shape and broad tails visible without the previous logarithmic stretching. All 48 histograms retain every finite pair; no tail was removed to improve visual collapse.

## Coherence distribution per displayed condition

Quantiles are over all retained frames and seeds; these are descriptive frame distributions, not independent-sample confidence intervals.

|Model|f|C label|Seeds|Mean coherence|5%|25%|Median|75%|95%|Fraction coherence<0.1|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|Wave-pinning|0.03|20|96|0.486|0.457|0.484|0.493|0.499|0.506|0.001|
|Wave-pinning|0.03|60|384|0.462|0.207|0.477|0.510|0.516|0.522|0.011|
|Wave-pinning|0.3|20|384|0.341|0.104|0.247|0.367|0.444|0.509|0.046|
|Wave-pinning|0.3|60|96|0.164|0.043|0.100|0.155|0.217|0.315|0.248|
|Holmes|0.02|20|96|0.387|0.364|0.388|0.395|0.400|0.407|0.007|
|Holmes|0.02|60|192|0.328|0.117|0.265|0.392|0.405|0.411|0.037|
|Holmes|0.3|20|192|0.277|0.094|0.207|0.289|0.352|0.426|0.057|
|Holmes|0.3|60|96|0.161|0.043|0.102|0.156|0.214|0.296|0.243|
|Spring|0.005|1|96|0.128|0.120|0.126|0.128|0.130|0.136|0.003|
|Spring|0.005|3|96|0.159|0.154|0.159|0.161|0.163|0.165|0.013|
|Spring|0.3|1|384|0.180|0.068|0.136|0.177|0.220|0.299|0.117|
|Spring|0.3|3|384|0.128|0.034|0.080|0.122|0.170|0.240|0.367|
|Goryachev|0.005|20|192|0.612|0.188|0.500|0.678|0.765|0.838|0.015|
|Goryachev|0.005|60|96|0.926|0.905|0.919|0.927|0.935|0.945|0.000|
|Goryachev|0.3|20|96|0.166|0.044|0.102|0.157|0.221|0.321|0.241|
|Goryachev|0.3|60|384|0.050|0.013|0.030|0.047|0.066|0.097|0.958|

## Lag relative to counted relocation interval

The following diagnostics use the frozen round6 dominant-track event detector after removal of 5-frame flicker. They are separate from the ungated heading distributions. The mean waiting interval is reciprocal pooled counted-switch rate; it is a detector-dependent proxy, not a threshold-free biological relocation time. We call 16τ not well below this interval when 16τ/mean wait ≥ 0.1, purely as a caution about overlapping events. Unscored switches and low coherence limit this comparison, especially strong-input Goryachev C60. No counted switch gives no finite waiting estimate, not proof of no relocations. The companion CSV also reports raw consecutive-frame track changes to expose detector sensitivity.

|Model|f|C label|16τ|Counted mean wait|16τ / mean wait|Not well below|
|---|---:|---:|---:|---:|---:|---|
|Wave-pinning|0.03|20|800|12094.5|0.066|no / unresolved|
|Wave-pinning|0.03|60|800|10072.1|0.079|no / unresolved|
|Wave-pinning|0.3|20|800|223.8|3.575|yes|
|Wave-pinning|0.3|60|800|115.7|6.912|yes|
|Holmes|0.02|20|800|9974.0|0.080|no / unresolved|
|Holmes|0.02|60|800|4426.5|0.181|yes|
|Holmes|0.3|20|800|377.1|2.122|yes|
|Holmes|0.3|60|800|110.3|7.252|yes|
|Spring|0.005|1|80|30720.0|0.003|no / unresolved|
|Spring|0.005|3|80|5907.7|0.014|no / unresolved|
|Spring|0.3|1|80|35.9|2.230|yes|
|Spring|0.3|3|80|10.3|7.768|yes|
|Goryachev|0.005|20|800|2221.7|0.360|yes|
|Goryachev|0.005|60|800|unresolved (zero counted)|0.000|no / unresolved|
|Goryachev|0.3|20|800|272.6|2.935|yes|
|Goryachev|0.3|60|800|1090.8|0.733|yes|

## What the raw distributions show

Wave-pinning (weak f = 0.03): the median lag-scaled absolute change grows from about 0.031 to 0.064 at C20 and 0.030 to 0.054 at C60 from τ to 16τ. Its outer width changes less between 4τ and 16τ. There is no full three-lag core collapse.

Holmes (weak f = 0.02): both the scaled core and outer quantiles broaden. The median rises from 0.038 to 0.078 at C20 and 0.054 to 0.096 at C60; the 95th percentile also increases. This panel alone does not resolve a separate relocation tail from the continuously broadening distribution.

Spring (weak f = 0.005): the scaled median increases from 0.0080 to 0.0164 at C1 and 0.0089 to 0.0176 at C3. The outer quantile is nearly stable between 4τ and 16τ, while the core continues to broaden.

Goryachev (weak f = 0.005): at C20, the scaled median remains narrow (0.0177, 0.0198, 0.0137) while the raw probability beyond 1.5 physical units rises from 4% to 13% to 22%. At C60 there is no probability beyond 1.5, median coherence is 0.927 and no takeover is counted. This is consistent with suppression during the observed window, rather than proof that jumps exclusively carry heading loss or a measured shift of the relocation-onset input. No onset-versus-input curve is presented.

At strong f = 0.3, Wave-pinning, Holmes and Spring spread across most of the bounded wrapped-displacement range by 4τ–16τ. Their lag-scaled outer widths then narrow, which reflects bounded angular displacement rather than a Brownian collapse. Goryachev also has broad raw phase changes: C20 approaches the wrapped range by 16τ while C60 continues to spread across lags. The low coherence of strong-input Goryachev C60 prevents identifying those phase changes uniquely as relocations. All statements are descriptive; no pass criterion was imposed. Quantiles and raw fixed-distance tail masses for all 48 displayed histograms are in D_descriptive_lag_quantiles.csv.

## Drift, α and partial closure

A uses native deterministic response with shared translation, no gain recalibration to noisy data. B pools the eligible point set across f but does not average away input dependence. The exclusion inventory is retained in B_alpha_bridge.csv; Goryachev has no eligible bridge estimate. C uses the original complete condition inventory and identifies missing drift as a nonnegative, unbounded quantity. Its intervals are angular-loss bounds, not sampling uncertainty. Domain-translation closure fails in many strong-input conditions; a numerical interval containing one does not validate the event mechanism.

## New strong-input campaign

### Strong-input findings

The new campaign uses 64 seeds per condition. Numerical values are reported from the frozen protocol without tuning. The forcing is the unchanged repository NoiseInput on the native grid, retaining its spatial zero mode and normalized OU-knot interpolation. The geometric mean of per-seed half-correlation rates is the declared estimator.

The small Holmes fast-cytosol effect is now resolved: the mode-1 size exponent falls by about 0.08 at both strong inputs, with paired seed intervals excluding zero. Wave-pinning shows a larger decrease, about 0.54 at input fraction 0.1 and 0.44 at 0.3. The fast-cytosol Wave-pinning heading exponent at 0.3 is slightly negative. Thus disappearance of the positive size effect is supported, while an exactly zero exponent is not imposed.

Patch count grows approximately with circumference, but the local statistics are not size independent. At threshold 0.65 and input fraction 0.1, C20 to C60 density falls 20.5%, mean width rises 28.2%, and RMS matched-patch speed rises 23.4%. At input fraction 0.3 the changes are −9.6%, +11.1%, and +18.7%. The same directions occur at the other thresholds. The measured speeds, about 0.04–0.05 at fraction 0.1 and 0.085–0.10 at fraction 0.3, exceed the pilot's approximate 0.03 and 0.05. Our declared centroid-overlap estimator and unchanged repository forcing differ from the pilot's unspecified overlap details and standalone forcing, so these differences are reported rather than corrected.

### Exponents and pilot comparisons

| Model | f | Quantity | Native alpha (95% CI) | Fast alpha (95% CI) | Native minus fast (95% CI) | Pilot difference |
|---|---:|---|---|---|---|---:|
| Wave-pinning | 0.1 | mode1 | 0.902 [0.848,0.957] | 0.358 [0.317,0.400] | 0.544 [0.500,0.592] | 0.49 |
| Wave-pinning | 0.1 | heading | 0.717 [0.657,0.777] | 0.211 [0.166,0.257] | 0.505 [0.455,0.555] | 0.65 |
| Wave-pinning | 0.3 | mode1 | 0.458 [0.434,0.482] | 0.020 [-0.000,0.039] | 0.438 [0.425,0.452] | 0.30 |
| Wave-pinning | 0.3 | heading | 0.406 [0.376,0.435] | -0.027 [-0.052,-0.004] | 0.433 [0.417,0.450] | 0.41 |
| Holmes | 0.1 | mode1 | 0.819 [0.759,0.881] | 0.738 [0.681,0.799] | 0.081 [0.046,0.117] | 0.05 |
| Holmes | 0.1 | heading | 0.676 [0.607,0.747] | 0.592 [0.526,0.664] | 0.084 [0.042,0.128] | 0.10 |
| Holmes | 0.3 | mode1 | 0.238 [0.209,0.267] | 0.158 [0.128,0.186] | 0.081 [0.076,0.086] | 0.10 |
| Holmes | 0.3 | heading | 0.129 [0.093,0.164] | 0.055 [0.019,0.089] | 0.074 [0.066,0.081] | 0.18 |

A confidence interval containing zero leaves that native/fast exponent difference unresolved. The 64-seed design is fixed; no seeds were selected by outcome. The local exponent here is the negative least-squares slope across three log sizes and geometric mean rates, rather than a differential exponent at one circumference.

### Threshold sensitivity

All patch statistics use native Wave-pinning fields, with active concentration above 0.65, and components spanning at least 10 cells (0.5 length units). Checks repeat the identical analysis at 0.55 and 0.75; no births or deaths are counted. Density is mean count divided by circumference; width pools all patch widths within each seed; speed is the RMS signed physical centroid velocity of mutual-best-overlap matches. Sampling is every 5 time units. Ambiguous splitting/merging and nonoverlapping regions supply no velocity. Values below are means across 64 seeds.

| f | Threshold | Density C20,C40,C60 | Width C20,C40,C60 | RMS speed C20,C40,C60 |
|---:|---:|---|---|---|
| 0.1 | 0.55 | 0.05809, 0.04822, 0.0465 | 6.598, 7.945, 8.244 | 0.04134, 0.04936, 0.05002 |
| 0.1 | 0.65 | 0.05742, 0.0474, 0.04563 | 6.095, 7.493, 7.812 | 0.04008, 0.04876, 0.04947 |
| 0.1 | 0.75 | 0.05662, 0.04661, 0.04478 | 5.565, 6.994, 7.337 | 0.04003, 0.05069, 0.05124 |
| 0.3 | 0.55 | 0.07334, 0.06711, 0.06634 | 6.324, 6.874, 6.957 | 0.09183, 0.1047, 0.1072 |
| 0.3 | 0.65 | 0.07289, 0.06667, 0.06588 | 5.679, 6.224, 6.31 | 0.0849, 0.09854, 0.1008 |
| 0.3 | 0.75 | 0.07202, 0.06613, 0.06549 | 5.016, 5.522, 5.596 | 0.08054, 0.09446, 0.09637 |

Threshold changes alter absolute widths and detected density. These finite-threshold summaries describe the patch field; they do not establish a threshold-independent point process or identify patch creation mechanisms.

### Numerical and statistical protocol

Wave-pinning runs last 32,000 time units and discard the first 16,000; Holmes runs last 20,000 and discard 4,000. Both use spatial step 0.05 and time step 0.05; all recorded frames are 5 time units apart. Native inactive diffusivities are 10 and 50 respectively; the fast control is 1000. Reaction terms, total-pool concentrations, initial frozen fields and all other kinetics are unchanged. Exact frozen basal inputs, worker source hashes, job IDs, and per-seed compact results are retained under ../noise/. The solver applies explicit native reaction and implicit periodic finite-difference diffusion, with a uniform roundoff correction restoring each species' diffusion-preserved zero mode.

The maximum relative conserved-pool deviation was 2.87e-12; minimum concentration was 0.000701. Native/fast forcing hashes agreed in all 768 pairs.

Four paired time-step-halving controls were run for each model and diffusivity at C40 and both input fractions 0.1 and 0.3, with the identical stochastic input. These controls assess discretization sensitivity of stochastic observables, not long-trajectory pointwise agreement. See EF_step_controls.csv. The native spatial mesh is retained; no stochastic mesh-refinement campaign was performed. Seed intervals therefore exclude spatial-discretization and estimator systematic uncertainty.

| Quantity | Median absolute relative change | Maximum absolute relative change |
|---|---:|---:|
| mode1_rate | 0.064% | 1.542% |
| heading_rate | 0.051% | 4.416% |
| P_abs_mean | 0.012% | 0.810% |

Confidence intervals resample 64 independent seed blocks 10,000 times, jointly across sizes and diffusivities. Rates are extracted per seed, so intervals summarize independent-run variability without treating correlated frames as independent. Plotted mode rates are geometric means. Mode autocorrelations are uncentered because nonzero spatial modes have zero ensemble mean by rotational symmetry. First-frame decay gives a rate lower bound of log(2)/5; nondecay by the longest lag gives an upper bound. The full censor status is retained per seed.

All plotted quantities derive from full-precision fields; float16 kymographs are illustrations only. The post-burn-in Fourier coefficients and coherence are float32/complex64 compact records. Neither heading nor modes are conditioned on coherence. Heading persistence is an activity-vector readout; it is distinct from simulated polarity replacement, observed cell migration, and conditional OU velocity-memory inference.

### Mathematical definitions

$$a_n(t)=\frac{1}{C}\int_0^C u(x,t)e^{-2\pi i n x/C}\,\mathrm{d}x,\qquad q_n=\frac{2\pi n}{C}.\tag{1}$$

$$R_n(s)=\frac{\operatorname{Re}\left\langle a_n(t+s)a_n(t)^*\right\rangle_t}{\left\langle |a_n(t)|^2\right\rangle_t},\qquad \gamma_n=-\frac{\log R_n(s_{1/2})}{s_{1/2}}.\tag{2}$$

The selected lag is the sampled lag closest to correlation 0.5 immediately before or at its first crossing. Mode-1 heading uses the unit complex coefficient in the same autocorrelation. Its coherence is the magnitude of mode 1 divided by the mean mode. Low coherence is retained; it makes heading a less sharply defined activity readout.

$$\alpha=-\frac{\sum_{j=1}^3(\log C_j-\overline{\log C})\,\overline{\log\gamma_j}}{\sum_{j=1}^3(\log C_j-\overline{\log C})^2},\qquad C_j\in\{20,40,60\}.\tag{3}$$

The overbar on log rate averages independent seeds. The native-minus-fast difference uses the same resampled seed blocks.

### Reproduction

Run the producer prepare and submit commands with the existing scaling_sandbox Python, then the separate strongest-input control prepare command. The immutable manifests record their SLURM commands. After every worker has stopped, run round7_noise_analysis_20261026.py analyze, then its validate command. Plot helpers are imported by the Figure 5 renderer. Reproduction writes into a fresh output location; the delivered frozen run must not be overwritten.

Generated 2026-10-26T14:58:37.682854+00:00. Producer COMPLETE; no timed browser refresh.


## Reproduction and provenance

`round7_persistence_20261026.extract()` writes D from frozen compact archives; run it with SLURM for bulk file reads. `round7_persistence_render_20261026.py --figure5` assembles all six panels after the new strong-input campaign analysis completes. Exact compact source paths and hashes are in supporting/metadata/D_compact_sources.csv, indexed by the immutable complete-campaign CSV hash. No dropped round6 panel is reproduced as a supplement.

## Prospective cellular assays and predicted signatures

Apply a weak boundary pulse and track both lateral edges in matched WGD− and WGD+ cells. If input statistics are size independent, the squared WGD+/WGD− ratio of permanent edge displacement per dose predicts the physical drift-diffusivity ratio. Lower common boundary motion in WGD+ would favour Wave-pinning-like drift suppression; similar physical motion with longer heading memory would favour Holmes/Spring-like geometry. Goryachev has a small deterministic edge response but no validated quantitative drift bridge.

Measure raw heading changes at several lags and resolved secondary-region takeovers, reporting coherence alongside both. A narrow high-coherence core with fewer large displacements in WGD+ is consistent with relocation suppression; low-coherence phase changes cannot identify domain jumps. To distinguish Goryachev-like pool depletion from slow drift, combine this with the Figure2 opposite-site first-rise assay: a delayed challenger at the same local input in WGD+ is the predicted depletion signature.

Estimate inactive-pool mobility, for example by FRAP with a diffusion/exchange model, and combine effective diffusivity with perimeter to estimate C/C*. If WGD-associated enlargement moves cells toward the conduction crossover, WGD+ should have smaller edge gain and greater persistence than geometry alone predicts; when C is much smaller than C*, nearly unchanged physical drift and geometric persistence are expected. Mobility/perimeter alone does not establish Goryachev-like competitor suppression; the first-rise assay is also needed.

At strong input, compare full-perimeter Fourier-mode correlations against physical wavenumber, with coherence and activity-field class retained. Increased inactive mobility is predicted to raise the low-wavenumber rates preferentially and reduce the Wave-pinning-like size advantage. The Holmes perturbation reduces its strong-input mode-1 exponent by about 0.08; a substantial mode-1 size effect remains and is unexplained. These tests concern orientation of whole-ring activity imbalance and do not assume a cell-velocity law.
