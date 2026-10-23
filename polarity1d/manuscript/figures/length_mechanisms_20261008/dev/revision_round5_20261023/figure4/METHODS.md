# Native-amplitude heading accounting used in the main figure

The original normalized point-region reconstruction has a unit-modulus moment for a single region even when the native field has a broad active profile and $|P|<1$. Its unscaled motion fraction is therefore approximately $1/|P|$, with an offsetting residual. This is a normalization artefact, not evidence that a distinct biological process cancels motion. The main panel instead expresses the same three Shapley changes at the native endpoint mean moment amplitude. It does not change the native profiles, region detection, endpoints, lag, moment eligibility, births, takeover records or measured total heading loss.

Let

$$
r_t=|P_t|,\qquad \widetilde Q_t=r_tQ_t,\qquad \bar r=\frac{r_0+r_1}{2}. \tag{A.1}
$$

For each original symmetric Shapley change of $Q$, use

$$
\Delta P_j=\bar r\,\Delta Q_j,\qquad j\in\{\mathrm{motion,amount,membership}\}. \tag{A.2}
$$

The remainder closes the measured native moment change exactly:

$$
\Delta P_{\mathrm R}=\Delta P-\bar r\,\Delta Q
=\Delta(P-rQ)+\bar Q\,\Delta r,\qquad\bar Q=\frac{Q_0+Q_1}{2}. \tag{A.3}
$$

It contains native amplitude variation, within-region shape, unresolved activity and background; it is not a pure background term. Project the four components with equation 4.5 below. The same exact signed identity holds. In a single rigidly translating domain, $P_t=r_te^{\mathrm i\theta_t}$ and $Q_t=e^{\mathrm i\theta_t}$, so the radial remainder is orthogonal to the unit-heading chord: motion contributes exactly one and remainder zero, even when $r_t$ varies. A synthetic test varies both heading and amplitude at three lags; multiregion and full-campaign checks also close numerically. The corrected campaign contains the same 170 eligible conditions and 45,696 source histories. It was reanalysed through heading-only SLURM job 21333382 without regenerating events or noisy model runs.

`heading_amplitude_condition_metrics.csv` and `heading_amplitude_seed_metrics.csv` contain the main accounting. Original `birth_heading_condition_metrics.csv` and per-condition `seeds.csv` remain unchanged. The reviewed main and 40-page supplement are saved under `supporting/historical_unscaled_20261023/` and appended unchanged after the new supplementary pages. This convention is descriptive and not uniquely causal.

For the time-horizon check, use certified births with at least 80τ of remaining record, pool f>0.01, and hold this cohort fixed at every horizon from 1τ to 80τ. The supplementary curves show observed five-frame-dominance incidence for the same cohort, not lifetime probability. Extending a finite horizon cannot by itself identify lifetime success or remove certification selection.

## Original unscaled accounting retained as a diagnostic

# Original event and unscaled-heading methods

The unchanged completed campaign contains 260 designed conditions and 45,696 independent histories (96–384 per condition). No new noisy trajectories were generated. New analysis reads every compact archive, not the bounded illustration profiles. Source cohort, hashes, cadence, controls, drift and correlation estimates are inherited from round 4. Chemical saved cadence is 5; Spring cadence is 0.5. A five-frame hold therefore spans five samples (first-to-fifth elapsed time 20 or 2); the historical switch convention is retained.

## Perimeter-time exposure and persistent secondary onsets

An incident identity has its first saved appearance after the start of the analysis window. It is persistent if present in each of its first five consecutive frames. Its preceding frame must contain a resolved dominant incumbent for the secondary-onset cohort. Initial identities and onsets with fewer than five observed consecutive frames are excluded. Finite-history exposure uses trapezoidal actual circumference:

$$
\mathcal E=\sum_{j=0}^{N-2}\frac{C_j+C_{j+1}}{2}(t_{j+1}-t_j). \tag{4.1}
$$

The nominal size is a facet/condition label; the actual Spring circumference varies. Report rates in native inverse length–time units. A condition with no certifiable events has a measured zero lower bound; it does not establish that the true birth rate is zero.

## Why exact births and splits are only partially identifiable

The frozen exporter retained `region_id`, excess-weighted contour `region_heading`, complete native-basin `region_amount`, contour-cell `region_width`, amount availability, current contrast, moments and offsets. It discarded masks, endpoints, parent IDs, per-region peaks and topology events. Width plus an excess-weighted centroid does not identify a region interval; asymmetric positive within-contour weights can move the centroid while leaving its geometric support unchanged. The old aggregate region audit includes identity onsets and a different two-frame takeover rule; it is not the requested event source. Exact split counts and exact all-birth probabilities remain unavailable in all 260 conditions. Bounded float16, spatially binned and temporally thinned exemplars cannot recover discarded full-campaign native overlap masks.

Every qualified half-contrast component is contiguous on the periodic native material-cell grid. If its cell-quadrature width is less than half a ring, the positive-weight circular centroid lies in that component’s shortest support arc. All occupied material cells, including half-cell boundaries, lie inside the conservative angular enclosure centred on the centroid with radius equal to the **full width**, not the half-width. A strict 2×10⁻⁶-radian margin protects float32 storage. A new region is a certified no-parent birth only if its enclosure misses every preceding enclosure; any component at least a half-ring wide or with a missing centroid prevents certification. Spring uses angular material-cell support so changing physical circumference does not invalidate the overlap test. Certified persistent secondary births give a lower bound; all persistent secondary onsets give an upper bound, since the latter include unknown split descendants. Neither uncertified onsets nor missing parents are called splits. The enclosure was checked against 2,000 asymmetric positive-weight contiguous supports including seam crossings.

## Distance and finite-horizon takeover

The incumbent is the preceding raw largest-amount qualified track. The global high-minus-low contrast is an explicit estimator of the incumbent plateau excess height; the compact source does not contain its local plateau. Using preceding-frame incumbent amount and height, distance is

$$
d=\frac{C_{\rm birth}}{2\pi}\left|\operatorname{wrap}(\theta_{\rm birth}-\theta_{\rm inc})\right|-\frac{A_{\rm inc}}{2h},\qquad h=\max u-\min u. \tag{4.2}
$$

Half-width includes the full-basin amount, not the truncated half-contrast contour amount. This plateau approximation and any negative distance remain explicit. Distance bins are $d\leq0$, $0<d\leq C_{\rm birth}/6$ and $d>C_{\rm birth}/6$. They are scaled relative to actual birth circumference; this avoids mixing the different native length scales.

A takeover requires that the newborn identity become globally dominant for at least five consecutive saved frames. The event time is the fifth frame. A common follow-up of 20 noise correlation times is used: 1,000 chemical time units or 100 Spring units. Births with shorter remaining observation are excluded. Successes must occur by this horizon. The result is an observed common-horizon probability conditional on certification, not lifetime probability or the probability for all births. Pool designed inputs $f>0.01$ without balancing their event counts. Report each input separately in `takeover_by_condition_distance.csv`; short-followup exclusions, surviving/censored tracks and median available follow-up are recorded. Independent histories, including histories with zero events, are resampled together for 2,000 whole-history pointwise bootstrap draws. Bootstrap intervals do not resolve systematic certification selection.

## Exact endpoint moment convention

Use the existing condition-specific analysis lag, chosen by the inherited eligible correlation closest to 0.5. Every pair of exact saved endpoints at that lag is considered. The native normalized first moment and its point-region reconstruction are

$$
P_t=\frac{M_{1,t}}{M_{0,t}},\qquad
Q_t=\frac{\sum_{i\in I_t} A_{i,t}e^{\mathrm i\theta_{i,t}}}{\sum_{i\in I_t} A_{i,t}},\qquad R_t=P_t-Q_t. \tag{4.3}
$$

The complete-basin amounts and excess-weighted half-contrast centroids do not reconstruct each region’s first moment exactly. Therefore $R$ includes within-region shape and amplitude, centroids defined by a different support, untracked activity, and diffuse background. It is a **moment residual**, not solely background. In a broad single domain, $|Q|=1$ while $|P|<1$; substantial cancelling motion/residual contributions are consequently expected under this unit-normalized point-region convention.

For a lag pair, form three binary endpoint operations. Motion changes the positions of shared identities. Coexisting amount changes the amounts of shared identities. Membership removes identities present only at the first endpoint and inserts identities present only at the second, with their observed endpoint positions/amounts. Membership includes unclassified splits, mergers and threshold appearances/disappearances, not only the certified no-parent births of panel E(i). At each of the eight hybrids, recompute the normalization in $Q$. Average the contribution of each operation over the six possible orders (the three-factor Shapley convention):

$$
\Delta Q_j=\sum_{S\subseteq\{M,A,N\}\setminus\{j\}}
\frac{|S|!(2-|S|)!}{3!}\left[Q(S\cup\{j\})-Q(S)\right],\quad j\in\{M,A,N\}. \tag{4.4}
$$

Together with $\Delta R$, these changes add exactly to $\Delta P$. Shared-region normalized weight changes induced by birth/loss normalization are shared symmetrically across operations; this is a declared accounting convention rather than a unique process separation. Events entirely within the lag that leave the same endpoint identity set are not individually resolved by this endpoint analysis.

Let $H=P/|P|$. For each complex component $\Delta P_j$, define its signed contribution to the measured endpoint heading decorrelation by

$$
L_j=\frac{\operatorname{Re}\!\left[\Delta P_j\,\overline{(H_1-H_0)}\right]}{|P_0|+|P_1|},\qquad
\sum_j L_j=\frac{|H_1-H_0|^2}{2}=1-\cos(\theta_1-\theta_0). \tag{4.5}
$$

The identity follows from $\Delta P=\overline r\,\Delta H+\overline H\,\Delta r$ and $\operatorname{Re}(\overline H\,\overline{\Delta H})=0$ for unit endpoint headings. Here overbars on $r,H$ in that explanatory identity denote arithmetic endpoint averages; complex conjugation in equation 4.5 is explicit. The projection preserves signed cancellation and does not divide by an arbitrary sum of absolute component sizes. Pool by $F_j=\sum L_j/\sum(1-\cos\Delta\theta)$. Fractions can be negative or greater than one. They are descriptive and convention dependent, not independent causal shares. Whole-history resampling preserves every component’s covariance and uses the same numerator/denominator draw for each fraction.

Eligibility requires both native $|P|\geq0.1$, finite positive native total activity, finite available endpoint region amounts/centroids and positive totals in every hybrid. A hybrid with no domains cannot be normalized and excludes the pair; it is not assigned a zero region moment. The CSV records exclusions for low coherence, amount/centroid availability and empty hybrids, with potential/eligible pair counts. No eligible analysis lag remains unavailable. Scientific tests establish complex closure to 9.2×10⁻¹⁶ and heading-loss closure to 9.2×10⁻¹⁵ across the actual complete campaign; no inference of mechanistic sufficiency follows from this algebraic closure.

## Drift and historical closure

Physical drift is inherited without refitting from five-frame persistent dominant-region tenures. Normalize by the smallest nominal size at the same input. Fit the log–log input exponent at every size using all available $f\leq0.01$ points; do not change it to agree with the deterministic Figure 3 prediction. Panel D’s old drift-plus-replacement point closure is shown only at $f\leq0.01$ and multiple-region frame fraction $\leq0.05$:

$$
k_{\rm pred}=\left(\frac{2\pi}{\overline C}\right)^2D_s+\lambda_{\rm rep}\,\langle1-\cos\Delta\theta\rangle_{\rm scored}. \tag{4.6}
$$

The angular loss of unscored replacements is represented by the scored-event mean in that **conditional point estimate**. Missing-angle bounds in the historical diagnostic can change the closure conclusion. Among 33 multiple-front-majority conditions, 0 point estimates lie within 15%, but bounds permit the 15% band in 24 and span unity in 18. Goryachev has 37 finite historical point ratios 0.00575–0.344; its missing-angle bounds span unity in 11. The present endpoint allocation is defined independently of any putative failure of equation 4.6. Detector-free fits, full-input closure, effective replacement rate, multiple-region fraction and null-mixture exponents are supplementary. All preceding supplement pages are appended unchanged.
