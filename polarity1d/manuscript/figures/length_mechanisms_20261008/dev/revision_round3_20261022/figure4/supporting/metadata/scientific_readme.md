# Figure 4 round 3: scientific record

## Unexpected findings and unresolved conditions

**Polarity normalization differs from the historical scalar gate.** Retained round-2 trace P was normalized by the initial first moment, whereas this brief explicitly defines P = M₁/M₀. The completed native compact records store the requested P = M₁/M₀; the matched-definition historical comparison reconstructs this quantity from retained round-2 moments. The changed magnitude gate can alter endpoint eligibility, coverage, selected lags, local exponents and the displayed points relative to the appendix sketch. Historical round-2 figures and S0 records are preserved unchanged.

**Switch counting and displacement availability require separate masks.** Earlier round-3 preview code removed an accepted dominant identity whenever its centroid was absent, then counted switches only across the remaining adjacent segments. This added an unintended switch-count criterion. In the checked new-campaign Wave-pinning C=20, f=0.3, seed-0 history, track 162 is replaced by track 166 for two frames (times 31,425 and 31,430), then track 167 persists for at least five frames from time 31,435. Every frame contains observed regions with finite headings; only the previously accepted track 162 disappears at time 31,430. The two-frame flicker is ignored, so the persistent 162→167 replacement counts. Its missing centroid interrupts displacement estimation without cancelling the observed switch. Final tracked summaries use persistence-filtered dominance for switch counts, valid exposure and median tenure, with a separate centroid-valid mask for displacement. Genuine missing-dominance intervals and seed boundaries remain breaks. Historical provisional preview records retain their earlier method and must not be substituted for the rebuilt final tracked summaries.

**Completed compact campaign.** Model measurements use the completed from-scratch compact records identified in the campaign manifest.

135 of 260 local size neighborhoods have no reported exponent. Every failed neighborhood remains a row in `local_alpha.csv` and `omitted_estimates.csv`; the omission reasons identify missing members, inadequate endpoint coverage, or no shared lag with all neighborhood correlations strictly between 0.03 and 0.995. A nearly empty Goryachev facet is retained as an estimation result. No absent size or input is interpolated into a neighborhood.

Available local exponents span 0.05786–3.435. The plotting range is −1.5 to 4: outward triangles preserve off-scale estimates and interval tails. Wild exponents are retained, not clipped in the data or smoothed into the main figure.

The completed native model ensembles contain 96, 192, 288, 384 independent histories per condition. A 2,000-resample bootstrap does not increase the number of independent histories. Sparse original round-2 ensembles appear only in the explicitly labelled historical comparisons; their actual independent-history counts remain in the comparison CSV.

Tracked diffusivity is finite in 255 exported conditions; a finite predicted/observed ratio is available in 42. Zero counted switches give an empirical switch contribution of zero; a log axis omits that zero rather than displaying a positive substitute. `jump_no_event_upper95 = 3 / valid_time` is a Poisson zero-count reference bound under a constant-hazard assumption, not an assumption-free confidence bound for these models. Whole-seed bootstrap intervals are separate and can be [0, 0] when no event was observed. Insufficient tenures do not imply a fitted diffusivity, and undefined measures are never imputed. Every unavailable decomposition quantity has a status/reason in `drift_decomposition.csv` or `omitted_estimates.csv`.

Null simulations can be too stable at weak input or decorrelate before the shortest declared lag at strong input. Such cells have no eligible local exponent; analytic formulas are used only for validation, never as replacements in the main panels. The simulated null analytic-check diagnostics and precision shortfalls are retained even if they are unfavorable.

The frozen new campaign contains 260 model–size–input conditions. The current source contains 260 observed conditions, leaving 0 planned cells without source observations. The final local-neighborhood census uses the complete fixed five-size grid at every designed model input. [Every planned condition](campaign_manifest.json).

The corrected switch diagnostic is available for 260/260 primary tracked conditions (five-frame hold, MSD range one quarter of median dominance tenure). Across their actual retained ensembles, 104381 counted replacements in 115 conditions lack the accepted centroid immediately before or after the event; the earlier centroid-adjacency gate would have suppressed these observed replacements. This diagnostic counts events, not independent histories. [All condition-level counts and actual ensemble sizes](../data/drift_decomposition.csv), column `n_switches_with_missing_adjacent_centroid`. Every affected primary condition is listed below as (circumference label, input, counted replacements with a missing adjacent centroid).

Wave-pinning: (20.0, 0.05, 1), (20.0, 0.1, 8), (20.0, 0.3, 440), (20.0, 1.0, 1295), (30.0, 0.05, 2), (30.0, 0.1, 4), (30.0, 0.3, 456), (30.0, 1.0, 1253), (40.0, 0.05, 5), (40.0, 0.1, 25), (40.0, 0.3, 394), (40.0, 1.0, 1865), (50.0, 0.05, 5), (50.0, 0.1, 33), (50.0, 0.3, 515), (50.0, 1.0, 2562), (60.0, 0.01, 1), (60.0, 0.05, 6), (60.0, 0.1, 51), (60.0, 0.3, 680), (60.0, 1.0, 12168).

Holmes: (20.0, 0.02, 1), (20.0, 0.05, 8), (20.0, 0.1, 29), (20.0, 0.3, 512), (20.0, 1.0, 3471), (30.0, 0.03, 3), (30.0, 0.05, 26), (30.0, 0.1, 128), (30.0, 0.3, 900), (30.0, 1.0, 14731), (40.0, 0.03, 9), (40.0, 0.05, 50), (40.0, 0.1, 277), (40.0, 0.3, 1768), (40.0, 1.0, 5348), (50.0, 0.02, 9), (50.0, 0.03, 23), (50.0, 0.05, 76), (50.0, 0.1, 399), (50.0, 0.3, 5115), (50.0, 1.0, 6659), (60.0, 0.02, 8), (60.0, 0.03, 17), (60.0, 0.05, 100), (60.0, 0.1, 528), (60.0, 0.3, 3568), (60.0, 1.0, 7553).

Spring: (1.0, 0.05, 16), (1.0, 0.1, 17), (1.0, 0.3, 671), (1.0, 1.0, 596), (1.5, 0.03, 1), (1.5, 0.05, 8), (1.5, 0.1, 250), (1.5, 0.3, 877), (1.5, 1.0, 1514), (2.0, 0.02, 2), (2.0, 0.03, 8), (2.0, 0.05, 62), (2.0, 0.1, 622), (2.0, 0.3, 896), (2.0, 1.0, 2459), (2.5, 0.01, 1), (2.5, 0.02, 7), (2.5, 0.03, 20), (2.5, 0.05, 65), (2.5, 0.1, 542), (2.5, 0.3, 5598), (2.5, 1.0, 3568), (3.0, 0.01, 2), (3.0, 0.02, 8), (3.0, 0.03, 26), (3.0, 0.05, 117), (3.0, 0.1, 424), (3.0, 0.3, 7556), (3.0, 1.0, 4477).

Goryachev: (20.0, 0.005, 2), (20.0, 0.007, 13), (20.0, 0.01, 16), (20.0, 0.015, 37), (20.0, 0.02, 32), (20.0, 0.025, 46), (20.0, 0.03, 53), (20.0, 0.05, 54), (20.0, 0.1, 46), (20.0, 0.3, 31), (20.0, 1.0, 18), (30.0, 0.01, 11), (30.0, 0.015, 32), (30.0, 0.02, 34), (30.0, 0.025, 43), (30.0, 0.03, 43), (30.0, 0.05, 46), (30.0, 0.1, 51), (30.0, 0.3, 33), (30.0, 1.0, 16), (40.0, 0.015, 19), (40.0, 0.02, 28), (40.0, 0.025, 18), (40.0, 0.03, 25), (40.0, 0.05, 38), (40.0, 0.1, 19), (40.0, 0.3, 6), (40.0, 1.0, 5), (50.0, 0.02, 13), (50.0, 0.025, 9), (50.0, 0.03, 11), (50.0, 0.05, 13), (50.0, 0.1, 11), (50.0, 0.3, 4), (60.0, 0.025, 6), (60.0, 0.03, 2), (60.0, 0.05, 1), (60.0, 0.1, 1).

A separate [initial-cohort seed-0 switch-contract audit](../data/seed0_switch_contract_audit.csv) checks 260 condition archives and records 737 restored replacements across 47 affected conditions. Its scope is replicate 0 from the initial 260-condition campaign, not the complete final ensembles or a bootstrap analysis. The CSV retains every audited condition, event times, dominance and centroid exposure/tenure diagnostics, archive paths and per-archive SHA-256 hashes. Its limited-scope event counts must not be substituted for the final primary-condition counts above.

Model rate-precision qualification uses four disjoint groups: 128 conditions meet SE(ln k) ≤ 0.07 with at least 97.5% valid bootstrap draws; 36 have adequately defined bootstrap rates but exceed that SE target (36 at the 384-history cap); 8 finite point estimates have fewer than 97.5% valid bootstrap draws; and 88 conditions have no eligible single-condition rate. A conditional SE computed from the finite subset of an invalid bootstrap is not a qualified finite-rate precision assessment. [Every condition, actual seed count, selected lag, SE and bootstrap-valid fraction](../data/analysis_lag_metrics.csv).

Conditions at the seed cap that still exceed the finite-rate precision target: Wave-pinning C=20.0, f=0.001, SE=0.07311, Wave-pinning C=20.0, f=0.005, SE=0.07668, Wave-pinning C=20.0, f=0.01, SE=0.07073, Wave-pinning C=30.0, f=0.003, SE=0.07368, Wave-pinning C=30.0, f=0.01, SE=0.07476, Wave-pinning C=40.0, f=0.003, SE=0.07756, Wave-pinning C=40.0, f=0.005, SE=0.0786, Wave-pinning C=40.0, f=0.02, SE=0.08247, Wave-pinning C=50.0, f=0.003, SE=0.07922, Wave-pinning C=50.0, f=0.005, SE=0.07822, Wave-pinning C=50.0, f=0.02, SE=0.07701, Wave-pinning C=60.0, f=0.01, SE=0.07217, Wave-pinning C=60.0, f=0.02, SE=0.07626, Wave-pinning C=60.0, f=0.03, SE=0.07242, Holmes C=20.0, f=0.003, SE=0.07485, Holmes C=30.0, f=0.001, SE=0.07127, Holmes C=30.0, f=0.003, SE=0.07422, Holmes C=30.0, f=0.005, SE=0.0752, Holmes C=40.0, f=0.001, SE=0.07757, Holmes C=40.0, f=0.003, SE=0.07202, Holmes C=40.0, f=0.005, SE=0.07628, Holmes C=50.0, f=0.005, SE=0.07117, Holmes C=50.0, f=0.01, SE=0.07774, Holmes C=60.0, f=0.003, SE=0.07138, Holmes C=60.0, f=0.01, SE=0.08005, Spring C=1.0, f=0.0001, SE=0.07162, Spring C=1.0, f=0.0003, SE=0.0708, Spring C=1.5, f=0.0001, SE=0.07441, Spring C=1.5, f=0.0003, SE=0.07287, Spring C=1.5, f=0.001, SE=0.07136, Spring C=2.0, f=0.001, SE=0.07128, Spring C=2.5, f=0.0003, SE=0.07028, Spring C=2.5, f=0.001, SE=0.07776, Goryachev C=20.0, f=0.004, SE=0.1705, Goryachev C=40.0, f=1.0, SE=0.1636, Goryachev C=60.0, f=0.3, SE=0.223.

Conditions with fewer than 97.5% valid bootstrap rate draws: Wave-pinning C=60.0, f=1.0, n=384, valid fraction=0.8925, Spring C=1.0, f=0.05, n=384, valid fraction=0.9655, Spring C=1.0, f=0.3, n=384, valid fraction=0.8505, Spring C=1.5, f=0.1, n=384, valid fraction=0.851, Spring C=2.0, f=0.1, n=384, valid fraction=0.932, Spring C=2.5, f=0.3, n=384, valid fraction=0.858, Spring C=3.0, f=0.3, n=384, valid fraction=0.9025, Goryachev C=50.0, f=1.0, n=384, valid fraction=0.6545.

**Sparse estimable jump contributions do not imply sparse switching.** The primary tracked ensembles contain 915,757 counted replacements, of which 321,423 have valid five-frame whole-P mean headings on both sides. Effective jump contributions comprise 108 measured zeros, 1 positive estimate and 151 undefined event-bearing conditions. A zero means no counted replacement in the observed ensemble; an undefined value means the complete jump contribution cannot be estimated under the declared scoring rule. Applying the shared polarity threshold to every frame of a five-frame jump window is an explicit conservative extension of the autocorrelation gate; the brief states this threshold directly for correlation endpoints, rather than separately specifying jump-window eligibility. Requiring every counted event to be scored avoids substituting the conditional angle distribution of scorable events for the unobserved angles. The resulting omissions describe conditional unscorability, not absence of replacements or evidence of agreement with the closure.

Among the 42 estimable prediction/observation comparisons, ratios span 0.9054896–1.148153; 42 lie within 15% of unity, and 41 have measured zero switch contributions. Goryachev contributes no estimable comparison. This agreement is conditional on estimability and does not validate jump closure generally or establish closure in conditions with undefined ingredients.

**Tracked and detector-free replacement estimates disagree at measured zeros.** The linear-axis cross-check contains 108 defined nonnegative pairs: 107 tracked zeros and 1 positive tracked estimate. In 107 of those zero pairs, the detector-free fit is positive. [Every matched condition, actual values and forcing flag](../data/supplement_scatter_lambda_eff.csv). 1 additional defined tracked condition has no defined matched detector-free estimate: Goryachev C=60.0, f=1.0: tracked λ_eff=0, detector status unavailable_too_few_valid_increments, eligible increments 0. [All detector-free conditions and increment counts](../data/detector_free_mle.csv). Undefined estimates remain omitted; the linear axes retain the measured zeros without altering either estimator.

**Shared grayscale suppresses small-ring Goryachev contrast.** Both selected seed-0 fields are intact and finite. Their shared linear activity scale ends at 83.3125; the small-ring maximum is 9.21875, only 11.065% of that scale, whereas the large-ring narrow peak reaches 83.3125. The apparently pale small-ring panel reflects the mandated common activity scale, not a missing history. Its native exemplar equality and source hashes are checked independently in the final display audit.

**The null B curves cannot cover the full analytic sketch under the declared lag rule.** At weak input, correlation can remain above 0.995 at the maximum lag of 16,000; at strong input, it can already be below 0.03 at the minimum lag of 250. The permitted local-neighborhood lag window therefore truncates the simulated curves, including possible asymptotic regimes. More independent histories do not repair this time-resolution limit. The shared-estimator and no-imputation requirements take precedence over drawing a complete analytic-looking curve.

[Every omitted estimate, with exact model/size/input and reason](../data/omitted_estimates.csv) · [Complete local-neighborhood census, including unselected lags](../data/local_alpha.csv) · [Every observed lag, correlation, coverage and interval](../data/condition_lag_metrics.csv). These inventories include each failed condition; no failed condition is removed from its designed grid.

Wave-pinning — unavailable_no_shared_lag_with_C_0.03_0.995_and_coverage_0.5; exact (circumference label, input) pairs: (20.0, 0.0001), (30.0, 0.0001), (40.0, 0.0001), (50.0, 0.0001), (60.0, 0.0001), (20.0, 0.0003), (30.0, 0.0003), (40.0, 0.0003), (50.0, 0.0003), (60.0, 0.0003), (20.0, 0.001), (30.0, 0.001), (40.0, 0.001), (50.0, 0.001), (60.0, 0.001), (20.0, 0.3), (30.0, 0.3), (40.0, 0.3), (20.0, 1.0), (30.0, 1.0), (40.0, 1.0), (50.0, 1.0), (60.0, 1.0).

Holmes — unavailable_no_shared_lag_with_C_0.03_0.995_and_coverage_0.5; exact (circumference label, input) pairs: (20.0, 0.0001), (30.0, 0.0001), (40.0, 0.0001), (50.0, 0.0001), (60.0, 0.0001), (20.0, 0.0003), (30.0, 0.0003), (40.0, 0.0003), (50.0, 0.0003), (60.0, 0.0003), (50.0, 0.001), (60.0, 0.001), (20.0, 1.0), (30.0, 1.0), (40.0, 1.0), (50.0, 1.0), (60.0, 1.0).

Spring — unavailable_no_shared_lag_with_C_0.03_0.995_and_coverage_0.5; exact (circumference label, input) pairs: (1.5, 0.0001), (2.0, 0.0001), (2.5, 0.0001), (3.0, 0.0001), (2.5, 0.03), (3.0, 0.03), (1.0, 0.05), (1.5, 0.05), (2.0, 0.05), (2.5, 0.05), (3.0, 0.05), (1.0, 0.1), (1.5, 0.1), (2.0, 0.1), (2.5, 0.1), (3.0, 0.1), (1.0, 0.3), (1.5, 0.3), (2.0, 0.3), (2.5, 0.3), (3.0, 0.3), (1.0, 1.0), (1.5, 1.0), (2.0, 1.0), (2.5, 1.0), (3.0, 1.0).

Goryachev — unavailable_no_shared_lag_with_C_0.03_0.995_and_coverage_0.5; exact (circumference label, input) pairs: (20.0, 0.0001), (30.0, 0.0001), (40.0, 0.0001), (50.0, 0.0001), (60.0, 0.0001), (20.0, 0.0003), (30.0, 0.0003), (40.0, 0.0003), (50.0, 0.0003), (60.0, 0.0003), (20.0, 0.001), (30.0, 0.001), (40.0, 0.001), (50.0, 0.001), (60.0, 0.001), (20.0, 0.003), (30.0, 0.003), (40.0, 0.003), (50.0, 0.003), (60.0, 0.003), (20.0, 0.004), (30.0, 0.004), (40.0, 0.004), (50.0, 0.004), (60.0, 0.004), (20.0, 0.005), (30.0, 0.005), (40.0, 0.005), (50.0, 0.005), (60.0, 0.005), (20.0, 0.007), (30.0, 0.007), (40.0, 0.007), (50.0, 0.007), (60.0, 0.007), (30.0, 0.01), (40.0, 0.01), (50.0, 0.01), (60.0, 0.01), (40.0, 0.015), (50.0, 0.015), (60.0, 0.015), (40.0, 0.02), (50.0, 0.02), (60.0, 0.02), (40.0, 0.025), (50.0, 0.025), (60.0, 0.025), (40.0, 0.03), (50.0, 0.03), (60.0, 0.03), (30.0, 0.05), (40.0, 0.05), (50.0, 0.05), (60.0, 0.05), (30.0, 0.1), (40.0, 0.1), (50.0, 0.1), (60.0, 0.1), (20.0, 0.3), (30.0, 0.3), (40.0, 0.3), (50.0, 0.3), (60.0, 0.3), (20.0, 1.0), (30.0, 1.0), (40.0, 1.0), (50.0, 1.0), (60.0, 1.0).

Null B maximum absolute relative rate deviation from its analytic check is 0.02704 (same_rate, C=50.0, f=0.01584893192461114, selected lag=16000.0); signed deviation 0.02704. Actual replicate counts span 512–33830; 0 reported precision-met cells exceed SE(ln k)=0.01. [Every null condition and check](../data/null_B_analytic_agreement.csv).

Null A maximum absolute simulated-minus-analytic autocorrelation deviation is 0.05476 across the displayed endpoint/lag records. [All autocorrelations and analytic checks](../data/null_A_acf.csv).

Null A has 8 negative finite-ensemble autocorrelation points; the minimum is -0.038088 for drift, C=20.0, lag=22000.0. They are retained on the displayed slightly negative ordinate and in CSV. At long lags, a 1,000-history ensemble can fluctuate around a small positive expected correlation; no history or seed was replaced to suppress this behavior.

Spring physical-coordinate sensitivity: 60/60 condition ratios are finite. Literal fixed-origin coordinate / intrinsic-arc diffusivity spans 0.99993232–1.000183, a maximum point-estimate change of 0.0183%. Both definitions use identical tenures and the fixed primary MSD range; all available primary intrinsic-arc estimates are independently reproduced. The alternative is a sensitivity check and does not replace or tune the primary ingredient. [Every ratio, interval and coordinate convention](../data/spring_coordinate_sensitivity.csv).

For B, the literal Bernoulli probability λΔt exceeds one at strong inputs. The simulator therefore uses the exact sampled Brownian-plus-Poisson-reset transition, with reset probability 1 − exp(−λΔt); every reset is uniform. A retains its requested Bernoulli-step experiment and validates against the corresponding discrete-step rate. This protocol qualification is recorded, rather than silently clipping a probability.

**Historical comparison only: tracking cadence.** Earlier provisional round-2 reanalysis used chemical cadences 5 or 40 and Spring cadences 0.5 or 4. Its five-frame holds therefore spanned 25 or 200 native chemical time units, and 2.5 or 20 Spring time units. The completed native compact campaign uses cadence 5 for chemical models and 0.5 for Spring, giving five-frame hold times 25 and 2.5. `source_cadence` and `flicker_hold_time` identify every current value in [the decomposition CSV](../data/drift_decomposition.csv). The earlier cadence qualification applies to the historical provisional evidence; no mixed cadence is pooled within a current condition.

Endpoint eligibility is abs(P)>=.1; tracked detector eligibility is distinct.

Forcing decisions are resolved for all 260 completed native conditions. Exact historical-control flags and the newly collected matched-mean/lower-spread controls are joined without changing the round-2 rule; the final analysis receipt reports no unknown forcing decisions.

The campaign failure/status receipt records scheduler availability and the exact ready submission command. No old Figure 4 native archive is deleted before the complete new compact campaign passes the stated retention gate. Round-2 S0/H1/H2/H3/A6 history remains unchanged; the current analysis has no hypothesis gates.

## Panel key numbers

A: three point-null families; the displayed circumferences are 20 and 60, with fixed physical width 1.5, drift diffusivity 0.002, jump rate 1/7000, total duration 24,000 and step 20. The perimeter-rate variant multiplies the jump rate by C/40. The requested autocorrelation ensemble is 1,000 independent histories per family and size; actual simulated counts and uncertainty are recorded in the null receipt.

B: 41 logarithmically spaced designed inputs from 10⁻⁴ to 1 and five fixed circumferences (20, 30, 40, 50, 60). The target simulation precision is SE(ln k) ≤ 0.01. Actual stopping counts and any unresolved conditions are recorded; they must not be inferred from the number of bootstrap draws.

C: 125/260 local estimates resolved: Wave-pinning 37, Holmes 43, Spring 34, Goryachev 11. Largest/smallest fixed-reference tolerance folds: Wave-pinning 2.373, Holmes 2.155, Spring 3.088, Goryachev 4.122. All folds, including censored and undefined values, remain in `noise_tolerance_folds.csv`.

D — Wave-pinning: f = 0.03; C_small(reference)<0.4 and C_large(reference)>0.6; nearest log tolerance midpoint among qualifying designed inputs. Requested/actual window: 10469.461975971391 / 10469.461975971391.

D — Holmes: f = 0.02; C_small(reference)<0.4 and C_large(reference)>0.6; nearest log tolerance midpoint among qualifying designed inputs. Requested/actual window: 8219.595058279263 / 8219.595058279263.

D — Spring: f = 0.005; C_small(reference)<0.4 and C_large(reference)>0.6; nearest log tolerance midpoint among qualifying designed inputs. Requested/actual window: 535.3699300272787 / 535.3699300272787.

D — Goryachev: f = 0.01; C_small(reference)<0.4 and C_large(reference)>0.6; nearest log tolerance midpoint among qualifying designed inputs. Requested/actual window: 1347.214579681161 / 1347.214579681161.

E: 260 condition records; 255 finite tracked-drift estimates; 42 finite predicted/observed comparisons. The panel retains every defined descriptive multiple-front fraction.

E — Wave-pinning: 16 finite ratios; minimum 0.9055 at C=50.0, f=0.005; maximum 1.077 at C=20.0, f=0.001.

E — Holmes: 13 finite ratios; minimum 0.9254 at C=30.0, f=0.005; maximum 1.148 at C=50.0, f=0.001.

E — Spring: 13 finite ratios; minimum 0.9446 at C=2.0, f=0.0003; maximum 1.089 at C=1.5, f=0.0003.

E — Goryachev: no finite decomposition comparison.

## Methods and estimands

Length represents WGD-associated enlargement. The empirical motivation includes greater migration persistence; observed migration summaries, conditional OU velocity-memory estimates, simulated activity-heading memory and one-dimensional polarity-reversal readouts are distinct evidence. The present analysis measures simulated activity heading rather than cell displacement.

The complex activity vector is the first spatial moment divided by total activity. Heading autocorrelation averages the cosine of the angular change over saved endpoint pairs with polarity magnitude at least 0.1 at both endpoints. Endpoint coverage is the eligible-pair count divided by all possible endpoint pairs within the declared observation window.

$$
C(\ell)=\left\langle\cos[\theta(t+\ell)-\theta(t)]\right\rangle_{|P(t)|,|P(t+\ell)|\geq0.1},\qquad k(\ell)=-\frac{\ln C(\ell)}{\ell}.
\tag{1}
$$

Target lags use square-root-of-two steps from 250 to 16,000 for chemical models and 25 to 1,600 for Spring, clipped to the available observation window. The compact lag table records both target and actual saved-endpoint lags; the analysis receipt states its endpoint-rounding rule. No temporal interpolation is hidden.

Each local exponent uses the complete fixed five-size grid: the focal size and its immediate neighbors, with two members at the ends. Actual mean circumference is used in the fit. A shared lag requires correlation strictly between 0.03 and 0.995 and coverage at least 0.5 in every member. Among eligible lags, the selected lag minimizes the worst distance of ln(−ln C) from ln(ln 2). Lag selection is held fixed in bootstrap draws.

$$
\alpha=-\frac{\operatorname{Cov}(\ln \bar C,\ln k)}{\operatorname{Var}(\ln \bar C)}.
\tag{2}
$$

The single-condition analysis lag is the eligible correlation closest to 0.5 in absolute correlation distance. Noise tolerance uses the first downward crossing of absolute correlation 0.5 at reference lag 2,000 for chemical models and 200 for Spring, interpolated linearly in logarithmic input. Forced conditions are included and encoded by open markers. Reference-lag sensitivities (1,000/4,000, or 100/400 for Spring) remain in the compact records. Missing crossings, left/right censoring and gaps are retained explicitly.

Confidence intervals are pointwise 95% whole-seed bootstrap intervals from 2,000 resamples. Every time point, lag and tracked event belonging to a resampled independent history is resampled together. Size ensembles are independent. Undefined logarithms or inadequate valid-draw fractions remain undefined or unbounded according to the estimator receipt.

A dominance tenure is a maximal interval with a constant dominant largest-amount track ID after requiring a replacement to persist for five frames; shorter flickers retain the previously accepted identity. Switches, valid exposure and the median tenure use these persistence-filtered identities before checking centroid availability. Exposure includes saved time intervals whose two endpoints have nonmissing filtered identities. A genuinely missing dominant ID resets the accepted identity; neither an event nor exposure crosses that gap or a seed boundary. A missing centroid of the accepted identity splits only the displacement segment, without deleting that identity from the switch/exposure record or shortening the dominance tenure used to choose the MSD range.

The primary physical coordinate integrates wrapped angular centroid increments multiplied by midpoint circumference, within each continuous centroid-valid segment of a dominance tenure. No missing centroid is interpolated or replaced by the short-lived instantaneous winner. A uniform dilation at fixed material angle therefore contributes no heading-changing drift. This is an intrinsic arclength increment convention, rather than the literal product of unwrapped angle and instantaneous circumference. Spring coordinate sensitivity is retained separately where available.

$$
\Delta s_i=\frac{C_i+C_{i+1}}{4\pi}\,\operatorname{wrap}(\theta_{i+1}-\theta_i),\qquad s_j=\sum_{i<j}\Delta s_i.
\tag{3}
$$

The primary diffusivity fits MSD = 2Dₛτ through the origin through one quarter of the median dominance tenure, pooling only displacement pairs within the continuous centroid-valid segments. The jump rate is counted persistent replacements divided by valid dominance exposure. Mean jump cosine uses the whole-activity P heading averaged over the five saved frames immediately before and after each counted replacement; the windows can include another identity but cannot cross an invalid-P interval or a seed boundary. Both windows require every polarity magnitude to be at least 0.1. Switch counting does not depend on these heading windows: if any counted switch lacks a valid jump-angle score, the effective jump rate remains undefined rather than discarding that switch or imputing an angle.

$$
\lambda_{\mathrm{eff}}=\lambda\left(1-\langle\cos\Delta\theta\rangle\right),\qquad k_{\mathrm{pred}}=\left(\frac{2\pi}{\bar C}\right)^2D_s+\lambda_{\mathrm{eff}}.
\tag{4}
$$

The detector-free mixture uses wrapped-Gaussian angular increments with variance 2Dθδ plus uniform relocations of weight 1 − exp(−λδ), at the shortest cadence and twice that cadence. Physical diffusivity is DθC²/(4π²). Declared tracked sensitivities use persistence thresholds of 2 and 10 frames, and MSD fit limits of one eighth and one half of median tenure. The supplementary PDF shows the cross-method comparisons, measured zeros, sensitivity checks and diagnostics of undefined ingredients.

## Rebuild, portable records and provenance

From the repository root, run:

```bash
conda run --no-capture-output -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figure4_round3_20261022.py
```

The default rebuild reads compact records only. It does not start a campaign, read native full-field archives or delete old outputs. Use the explicitly named preparation options to regenerate scalar/track summaries from retained round-2 evidence or to rerun null simulations. The manifest and hashes bind each figure to its actual source.

- `local_alpha.csv`, `end_to_end_alpha.csv` and `analysis_lag_metrics.csv`: estimates, selected lags, intervals, valid fractions and statuses.
- `condition_lag_metrics.csv`: every saved correlation/rate curve point, coverage and seed count.
- `noise_tolerance.csv` and `noise_tolerance_folds.csv`: every reference lag, crossing/fold, interval and censoring status.
- `drift_decomposition.csv`, within-tenure and sensitivity CSVs: directly measured ingredients and qualifications.
- `null_A_*` and `null_B_*`: newly simulated null records, distinct from model evidence.
- `plotted_*.csv`: exact compact curve/point sources consumed by the main renderer.
- `omitted_estimates.csv`: unavailable/off-domain estimates with their original reasons.
- `model_D_selection.csv` and the kymograph index: deterministic displayed-input, replicate and time-window selection.
- `analysis_receipt.json`, `null_analysis_receipt.json`, `render_receipt.json`, campaign manifest/status and this documentation: source hashes, estimator rules and provenance.

The main PDF is 180 mm wide and 500 mm high. It contains no overarching title or subtitle. Axes retain only left/bottom spines, time increases upward, model facets remain black, and each model has a local circumference key. No table appears in the current figure or its new supplementary diagnostic pages.
