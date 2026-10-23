"""Publication legend and explicit identification/attribution methods."""
from __future__ import annotations
import html,json,math
from datetime import datetime,timezone
from pathlib import Path
import figure4_round5_analysis_20261023 as A
import figure4_round3_render_20261022 as R

def documents(legend_only=False):
    out=A.OUT;data=out/'supporting/data';meta=out/'supporting/metadata'
    rows=A.read_csv(data/'birth_heading_condition_metrics.csv');fits=A.read_csv(data/'weak_input_drift_exponents.csv')
    receipt=json.loads((meta/'round5_analysis_receipt.json').read_text())
    exponents={m:(min(float(r['exponent']) for r in fits if r['model']==m),max(float(r['exponent']) for r in fits if r['model']==m)) for m in ('wave_pinning','holmes_model3','spring_ring')}
    exponent_text='; '.join(R.NAMES[m]+f' {v[0]:.2f}–{v[1]:.2f}' for m,v in exponents.items())
    panels=[
      ('A','Drifting peaks, same-rate null relocations and perimeter-rate null relocations give size exponents 2, 0 and −1. Null kymograph widths follow circumference; time increases upward. Autocorrelations C(ℓ) use solid C20 and dotted C60 curves. The null relocations are illustrative whole-peak moves, distinct from measured dominant-region replacements and birth-to-takeover events.'),
      ('B','Larger rings retain heading across stronger fluctuating input. Local size exponents α and noise-tolerance folds use model hues and light-to-dark circumference ramps. Open markers identify forced conditions; vertical intervals are pointwise 95% whole-history bootstrap intervals. Boundary triangles retain off-scale values. 125/260 local neighborhoods resolve, including 11 Goryachev neighborhoods. Largest-size tolerance folds are Wave-pinning 2.373, Holmes 2.155, Spring 3.088 and Goryachev 4.122. f₁/₂ is the first input crossing of C(ℓ)=0.5 at reference lag 2000 (Spring 200); Spring circumference is multiplied by 20 only in the fold panel.'),
      ('C','Each model kymograph is rolled to place its initial activity-moment heading at the ring centre, to the nearest retained spatial bin. The two sizes share the existing selected input and replicate; the time window and raw evidence are unchanged. Activity scales are shared within Wave-pinning, Holmes and Spring. Goryachev uses a separate scale at each size: 0–9.2 at C20 and 0–83 at C60. Black outside ticks indicate persistent dominant-region replacements. Ensemble C(ℓ) below retains defined small- and large-size correlations. A low-contrast initial heading may have weak biological meaning even though the numerical roll is defined.'),
      ('D','Physical within-tenure drift Dₛ is normalized by the smallest size at the same input. Unity gives the equal-physical-drift geometric expectation α=2. Weak-input log–log fits for f≤0.01 give exponents '+exponent_text+'. The second row shows the existing event-based predicted/observed decay only where f≤0.01 and multiple regions occupy at most 5% of frames. 45/47 eligible point comparisons lie within 15% (Wave-pinning 16/16; Holmes 15/16; Spring 14/14; Goryachev 0/1). These point predictions assign the scored-event mean angular loss to unscored replacements; this remains an assumption. Goryachev uses separate logarithmic ordinate ranges: 10⁻⁷–2 for normalized drift and 0.003–3.3 for decay ratio; other models use 0.2–8 and 0.3–3.3. The actual mean circumference, C̄, enters physical-to-angular conversion. Figure 3D independently predicts this measured drift from deterministic edge responses.'),
      ('E','i, Certified secondary-domain birth rates per actual perimeter-time. A new track must persist at least five saved frames, have an incumbent in the preceding frame, and have a conservative support enclosure disjoint from every preceding region. The curves are lower bounds, not estimates of all no-overlap births: original masks and parent IDs were discarded. All persistent secondary onsets, including unresolved splits, supply upper bounds in the supplement. Exact split counts are unavailable for all 260 conditions. ii, Probability of five-frame global dominance within a common follow-up of 20 noise blocks (1000 chemical units; 100 Spring), conditional on certified no-parent births. Births with shorter available follow-up are excluded; pointwise 95% whole-history bootstrap intervals are shown. Designed inputs f>0.01 are pooled without input reweighting. Distance d is the shortest centroid arc minus Aᵢₙ꜀/(2h), where incumbent full-basin excess amount Aᵢₙ꜀ and global contrast h from the preceding frame estimate a plateau half-width. This is an explicit plateau approximation. Curves distinguish d≤0, 0<d≤C/6 and d>C/6; circumference bins use nominal size while exposure uses actual circumference. The Figure 2B–C asymmetry motivates measuring initiation and replacement separately. The observed takeover decline applies to this finite-horizon, certified subset; unequal certification prevents a population-level birth/takeover claim. iii, Exact signed allocation of eligible endpoint heading loss over each condition’s existing analysis lag: M, shared-region motion; A, shared-region amount; B, track membership changes; R, moment residual. Membership includes unclassified splits and mergers. Three-factor symmetric Shapley changes of the normalized point-region moment plus the measured-moment residual sum exactly to the observed moment change. Their projections sum to 1−cos Δθ. Signed pooled fractions sum to one and may be negative or exceed one; cancellations are retained. The residual includes within-region shape and amplitude, unresolved activity and diffuse background. The symmetric-log color scale is linear for |fraction|≤0.1 and spans −10 to 10; black dots mark values outside that range. Gray cells lack an eligible lag or valid endpoint/hybrid moments. These are convention-dependent descriptive allocations, not causal fractions.')
    ]
    narrative=[
      'Enlargement extends heading memory and raises noise tolerance, while the size advantage weakens as input strengthens. The weak-input drift in Wave-pinning, Holmes and Spring is close to quadratic in input, connecting Figure 4 to the deterministic edge-response test in Figure 3.',
      'Direct event accounting distinguishes persistent secondary onsets from certifiable no-parent births. Certified takeover probability declines with size in each model at a common finite follow-up, but different certification fractions and pooled inputs limit population inference. Per-perimeter birth direction varies by model; it is not a universal consequence of enlargement.',
      'Signed heading-loss accounting closes numerically while showing substantial cancellation between point-region changes and the moment residual. This exposes the limitations of a point-domain reconstruction; it does not establish that one mechanism causes a specified fraction of native heading loss.'
    ]
    limitations=[
      'Exact no-overlap birth counts and split counts cannot be recovered for the full campaign because compact records discard native masks/parent IDs. Lower/upper birth bounds and conditional takeovers are reported instead. The global-contrast plateau estimate is not a recovered per-region plateau.',
      'Native clocks are not calibrated between models. A takeover is five-frame global dominance within 20 noise blocks, not a lifetime success probability. Full-history right censoring, short-followup exclusions and per-condition counts are retained in CSV. Inputs f>0.01 are pooled without reweighting; size cohorts can differ in their event-weighted input composition.',
      'Heading attribution uses only eligible lag endpoint pairs with measured |P|≥0.1, finite tracked amounts/centroids and all nonzero hybrid totals. Goryachev at f=0.03 has pair coverage 0.666 at C20 but only 0.0109 at C60. Endpoint membership does not reveal events occurring and reversing inside the lag. This algebraic allocation depends on the chosen convention and includes normalized point-region approximation error in its residual.',
      'Historical multiple-front-majority closure point predictions are conditional on scored angles representing unscored replacements: none of 33 point estimates is within 15%, but missing-angle bounds permit agreement within 15% in 24/33 and span unity in 18/33. This cannot generally distinguish omitted heading processes from unscored large-angle replacements. Goryachev’s 37 historical finite point ratios span 0.00575–0.344, with limited dominant-region exposure and angle scoring; missing-angle bounds span unity in 11/37. These diagnostics remain in the supplement and do not motivate the independent exact endpoint allocation.'
    ]
    sources=[]
    for label,path in [('Round5 complete campaign analysis',meta/'round5_analysis_receipt.json'),('Scientific identity/enclosure validation',meta/'scientific_validation.json'),('Birth and heading-loss conditions',data/'birth_heading_condition_metrics.csv'),('Takeovers with common follow-up',data/'takeover_by_size_distance.csv'),('Original campaign manifest',A.CAMPAIGN/'manifest.json')]:
        sources.append(dict(label=label,path=str(path.relative_to(A.REPO)),sha256=R.sha(path)))
    protocol=('Input fraction f sets the ensemble mean to f times the native reference '
      '(Wave-pinning 0.1; Holmes 1; Spring 4; Goryachev 6.4). Input is positive lognormal, '
      'mean × exp(σZ−σ²/2), with unit-variance Gaussian driving field Z and σ=1. '
      'The spatial parameter ℓ and OU correlation time τ are 1.5 and 50 in chemical models; '
      'Spring uses basal-physical ℓ=0.15 and τ=5, with physical ℓ scaling with its changing circumference. ')
    revised=[]
    for p,t in panels:
        t=t.replace('20 noise blocks','20 noise correlation times')
        if p=='B':
            t=t.replace('Local size exponents α',protocol+'Local size exponents α')
            t=t.replace('Open markers identify forced conditions;',
                'Open symbols in the exponent and drift panels identify conditions flagged as forced by '
                'the frozen matched-mean/lower-noise controls; open diamonds in the fold panel identify '
                'Goryachev irrespective of forcing;')
        if p=='E':
            t=t.replace('Distance d is the shortest centroid arc minus Aᵢₙ꜀/(2h), where incumbent full-basin excess amount Aᵢₙ꜀ and global contrast h from the preceding frame estimate a plateau half-width.',
                'Distance d is the shortest centroid arc minus the estimated incumbent half-width: '
                'the incumbent’s preceding-frame full-basin excess amount divided by twice the preceding-frame '
                'global contrast h.')
        revised.append((p,t))
    panels=revised
    limitations=[t.replace('20 noise blocks','20 noise correlation times') for t in limitations]
    equations=[
      dict(math=r'P_t=\frac{M_{1,t}}{M_{0,t}}=Q_t+R_t,\qquad Q_t=\frac{\sum_{i\in I_t}A_{i,t}e^{\mathrm{i}\theta_{i,t}}}{\sum_{i\in I_t}A_{i,t}}',number='4.1',
        description='A is complete region excess amount; θ is the tracked contour centroid angle. R is the measured moment residual, including region shape/amplitude, untracked activity and background.'),
      dict(math=r'\Delta P=\Delta P_{\mathrm{M}}+\Delta P_{\mathrm{A}}+\Delta P_{\mathrm{B}}+\Delta R',number='4.2',
        description='M, shared-track motion; A, coexisting amount; B, endpoint track membership including unresolved split/merge changes. The three changes in Q use the symmetric three-factor Shapley convention.'),
      dict(math=r'L_j=\frac{\mathrm{Re}\!\left[\Delta P_j\,\overline{(H_1-H_0)}\right]}{|P_0|+|P_1|},\quad H_t=\frac{P_t}{|P_t|},\quad\sum_j L_j=1-\cos(\theta_1-\theta_0)',number='4.3',
        description='Pooled signed fractions are sums of L_j divided by total measured heading loss. They retain cancellation, can exceed one, and depend on the stated accounting convention; they are not causal shares.'),
      dict(math=r'k(\ell)=-\frac{\ln C(\ell)}{\ell},\qquad k\propto\bar C^{-\alpha}',number='4.4',
        description='C(ℓ) is activity-heading autocorrelation; C̄ is actual mean circumference. Endpoint eligibility uses |P|≥0.1. Drift closure is shown only at f≤0.01 and multiple-region fraction≤0.05.')
    ]
    legend=dict(id='figure4',title='Enlargement extends heading memory as fluctuating input changes domain dynamics',
      role='Size-dependent heading persistence, edge drift and secondary-domain competition under fluctuating input.',
      status='complete_existing_campaign_with_partial_identification',
      assets=[dict(path=str((out/'figure4_main.png').relative_to(A.REPO)),caption='',page_size=[180/25.4,420/25.4],minimum_printed_width_inches=180/25.4,
        native_pdf_page=str((out/'figure4_main.pdf').relative_to(A.REPO)),page_title='')],
      legend=[letter+'. '+text for letter,text in panels],legend_panels=[dict(panel=p,text=t) for p,t in panels],
      narrative=narrative,equations=equations,measurements=['Measure lateral-edge fluctuations, secondary-protrusion births per actual perimeter-time, finite-horizon takeover probability and leading-edge displacement after boundary-straddling photoactivation. Preserve separate WGD-associated cell-size, observed migration and simulated polarity-heading readouts.'],
      limitations=limitations,sources=sources,main_pdf=str((out/'figure4_main.pdf').relative_to(A.REPO)),
      supplement_pdf=str((out/'figure4_supplement.pdf').relative_to(A.REPO)),reader=str((out/'index.html').relative_to(A.REPO)),campaign_completed=True)
    (out/'legend.json').write_text(json.dumps(legend,indent=2,ensure_ascii=False)+'\n')
    methods=r'''# Figure 4 methods: existing-campaign event and heading analysis

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

A takeover requires that the newborn identity become globally dominant for at least five consecutive saved frames. The event time is the fifth frame. A common follow-up of 20 independent noise-block durations is used: 1,000 chemical time units or 100 Spring units. Births with shorter remaining observation are excluded. Successes must occur by this horizon. The result is an observed common-horizon probability conditional on certification, not lifetime probability or the probability for all births. Pool designed inputs $f>0.01$ without balancing their event counts. Report each input separately in `takeover_by_condition_distance.csv`; short-followup exclusions, surviving/censored tracks and median available follow-up are recorded. Independent histories, including histories with zero events, are resampled together for 2,000 whole-history pointwise bootstrap draws. Bootstrap intervals do not resolve systematic certification selection.

## Exact endpoint moment convention

Use the existing condition-specific analysis lag, chosen by the inherited eligible correlation closest to 0.5. Every pair of exact saved endpoints at that lag is considered. The native normalized first moment and its point-region reconstruction are

$$
P_t=\frac{M_{1,t}}{M_{0,t}},\qquad
Q_t=\frac{\sum_{i\in I_t} A_{i,t}e^{\mathrm i\theta_{i,t}}}{\sum_{i\in I_t} A_{i,t}},\qquad B_t=P_t-Q_t. \tag{4.3}
$$

The complete-basin amounts and excess-weighted half-contrast centroids do not reconstruct each region’s first moment exactly. Therefore $B$ includes within-region shape and amplitude, centroids defined by a different support, untracked activity, and diffuse background. It is a **moment residual**, not solely background. In a broad single domain, $|Q|=1$ while $|P|<1$; substantial cancelling motion/residual contributions are consequently expected under this mandated normalization.

For a lag pair, form three binary endpoint operations. Motion changes the positions of shared identities. Coexisting amount changes the amounts of shared identities. Membership removes identities present only at the first endpoint and inserts identities present only at the second, with their observed endpoint positions/amounts. Membership includes unclassified splits, mergers and threshold appearances/disappearances, not only the certified no-parent births of panel E(i). At each of the eight hybrids, recompute the normalization in $Q$. Average the contribution of each operation over the six possible orders (the three-factor Shapley convention):

$$
\Delta Q_j=\sum_{S\subseteq\{M,A,N\}\setminus\{j\}}
\frac{|S|!(2-|S|)!}{3!}\left[Q(S\cup\{j\})-Q(S)\right],\quad j\in\{M,A,N\}. \tag{4.4}
$$

Together with $\Delta B$, these changes add exactly to $\Delta P$. Shared-region normalized weight changes induced by birth/loss normalization are shared symmetrically across operations; this is a declared accounting convention rather than a unique process separation. Events entirely within the lag that leave the same endpoint identity set are not individually resolved by this endpoint analysis.

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
'''
    methods=methods.replace('B_t=P_t-Q_t','R_t=P_t-Q_t').replace('Therefore $B$','Therefore $R$').replace('$\\Delta B$','$\\Delta R$')
    methods=methods.replace('20 independent noise-block durations','20 noise correlation times').replace('20 noise blocks','20 noise correlation times')
    if not legend_only:(out/'METHODS.md').write_text(methods)
    table=['| Model | Small birth lower–upper | Large birth lower–upper | Small→large certification |','|---|---:|---:|---:|']
    for m in R.MODELS:
        rr=sorted([r for r in rows if r['model']==m and float(r['f'])==.03],key=lambda r:float(r['C_label']))
        lo,hi=rr[0],rr[-1]
        bounds=lambda r:f"{float(r['certified_birth_rate_lower_bound']):.3g}–{float(r['persistent_onset_rate_upper_bound']):.3g}"
        table.append(f"| {R.NAMES[m]} | {bounds(lo)} | {bounds(hi)} | {float(lo['certification_fraction']):.3f}→{float(hi['certification_fraction']):.3f} |")
    results='''# Figure 4 results

Enlarged rings retain heading across stronger fluctuating input. Largest-size noise-tolerance folds are 2.373 (Wave-pinning), 2.155 (Holmes), 3.088 (Spring) and 4.122 (Goryachev). The size advantage weakens as input strengthens. Weak-input physical-drift exponents are '''+exponent_text+'''.

The new event analysis covers every one of the 45,696 independent histories, yielding 2,064,445 persistent secondary track onsets and 1,195,029 certified no-parent births. Exact split counts remain unavailable because parent overlaps were not retained. Certification varies strongly with model, size and input; the certified curves cannot establish a universal birth-rate increase with enlargement.

At input fraction 0.03, native per-perimeter-time birth identification bounds and certification are:

'''+ '\n'.join(table)+'''

At f=0.1, Spring’s largest-size birth lower bound (0.01723) exceeds the smallest-size upper bound (0.01439), identifying an increase for that condition. Goryachev instead has a small-size lower bound of 0.0002540 exceeding the large-size upper bound 0.000002642, identifying a decrease. No universal secondary-protrusion birth direction should be preregistered from these models. At f=0.03 the Wave-pinning and Holmes small/large bounds overlap.

For certified births with d>C/6 and a common 20-noise-block follow-up, the observed takeover probability declines from the smallest to the largest size: Wave-pinning 0.549→0.223, Holmes 0.503→0.247, Spring 0.755→0.419, Goryachev 0.388→0.0718. These probabilities pool f>0.01 and apply to the selected certified subset. The data support measuring births and takeovers separately, as motivated by Figure 2B–C; they do not identify the population-wide takeover probabilities.

The exact endpoint heading-loss allocation preserves large signed cancellations. For example, Spring at f=0.03 has motion fractions 6.83→3.88 and moment residual fractions −6.36→−5.10 from smallest to largest size. Point-region moments have larger amplitudes than the native field moment, so these values are not physical percentages or evidence of negative causal mechanisms. Membership and amount contributions increase in some conditions, but their interpretation remains tied to the stated endpoint convention. Goryachev f=0.03 pair coverage falls from 66.6% at C20 to 1.09% at C60, preventing a representative whole-trajectory interpretation there. The source tables retain every negative, large and unavailable value.

The restricted historical drift/replacement closure has 47 finite point comparisons and 45 within 15%. The Goryachev comparison is 0.283; one Holmes comparison is 0.843. The existing all-input event-angle scoring limitation remains: missing-angle bounds permit agreement within 15% in 24/33 multiple-front-majority conditions, so the unscored angles prevent a general claim that an additional heading-loss process is required. The new exact allocation does not depend on that claim.
'''
    results=results.replace('20-noise-block','20-noise-correlation-time')
    if not legend_only:(out/'RESULTS.md').write_text(results)
    availability=[dict(condition_id=r['condition_id'],model=r['model'],C_label=r['C_label'],f=r['f'],
        exact_no_parent_births='partially_identified_only',exact_splits='unavailable_missing_parent_masks',
        exact_local_plateau='unavailable_global_contrast_estimator_used',
        heading_allocation='available' if float(r['eligible_pairs'])>0 else 'unavailable_no_eligible_analysis_lag_or_endpoints',
        analysis_lag=r['analysis_lag'],pair_coverage=r['pair_coverage']) for r in rows]
    if not legend_only:A.write_csv(data/'requested_measure_availability.csv',availability)
    readme='''# Figure 4: enlargement and heading memory

[Main figure](figure4_main.pdf) · [Legend/source contract](legend.json) · [Complete supplement](figure4_supplement.pdf) · [Methods](METHODS.md) · [Results](RESULTS.md)

'''+ '\n\n'.join(narrative)+'''

The complete existing campaign was analyzed: 260 conditions, 45,696 independent histories. No new noisy model runs were generated. Birth rates are certified lower bounds; exact no-parent birth and split classification was not retained in the compact schema. The requested strict all-birth/takeover measures therefore remain partly unidentifiable. `requested_measure_availability.csv` lists all cases. Finite-horizon takeovers use a common 20-noise-block follow-up. Every historical supplementary page is appended unchanged after the eight new pages.

Rebuild with the existing `scaling_sandbox` Conda Python:

```bash
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_validation_20261023.py
sbatch polarity1d/manuscript/code/mechanisms/figure4_round5_analysis_20261023.slurm
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_analysis_20261023.py collect
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_render_20261023.py
python -B polarity1d/manuscript/code/mechanisms/figure4_round5_documents_20261023.py
```

Analysis used SLURM job 21332720 (260 single-core condition tasks, small QoS, no array throttle). Data are generated once and opening/manually refreshing the reader loads those saved results; there is no browser timer. The analysis does not mutate raw campaign data, basals, historical figures, reviews or author text.

Key CSVs: `birth_heading_condition_metrics.csv`, `birth_heading_seed_metrics.csv`, `takeover_by_condition_distance.csv`, `takeover_by_size_distance.csv`, `takeover_seed_counts.csv`, `weak_input_drift_exponents.csv`, `kymograph_centering.csv`, `plotted_D.csv`, and `requested_measure_availability.csv`. Per-condition retained event rows live under `supporting/condition_analysis/` and cover the entire campaign.
'''
    readme=readme.replace('20-noise-block','20-noise-correlation-time')
    if not legend_only:(out/'README.md').write_text(readme)
    generated=datetime.now(timezone.utc).isoformat()
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure 4 — enlargement and heading memory</title><style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1080px;margin:2rem auto;padding:0 1rem;color:#222}}img{{width:100%;height:auto}}a{{color:#17619c}}li{{margin:.5rem 0}}small{{color:#666}}</style><h1>Enlargement extends heading memory</h1><p>Generated {generated}. Analysis complete; reload manually to read newly generated files.</p><p><a href="figure4_main.pdf">Main PDF</a> · <a href="figure4_supplement.pdf">Full supplement</a> · <a href="METHODS.md">Methods</a> · <a href="RESULTS.md">Results</a> · <a href="supporting/data/requested_measure_availability.csv">Missing cases</a></p><img src="figure4_main.png" alt="Figure4 panels A–E, with size-dependent heading memory, drift, bounded births, finite-horizon takeover and signed heading-loss decomposition"><h2>Legend</h2>'''
    page+=''.join('<p><b>'+letter+'.</b> '+html.escape(text)+'</p>' for letter,text in panels)
    page+='<h2>Limitations</h2><ul>'+''.join('<li>'+html.escape(t)+'</li>' for t in limitations)+'</ul></html>'
    (out/'index.html').write_text(page)
    (meta/'round5_document_receipt.json').write_text(json.dumps(dict(generated_utc=generated,title=legend['title'],
        update_mode='legend_reader_only' if legend_only else 'all_documents',
        document_source_sha256=R.sha(Path(__file__)),noise_protocol_verified_from='original campaign manifest and per-worker noise_audit.json',
        noise_sigma=1,noise_chemical=dict(ell=1.5,tau=50),noise_spring=dict(basal_ell=.15,tau=5),
        requested_measure_cases=260,exact_split_unavailable=260,analysis_receipt=receipt,
        retained_historical_supplement=json.loads((meta/'supplement_preservation.json').read_text())),indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--legend-only',action='store_true')
    documents(parser.parse_args().legend_only)
