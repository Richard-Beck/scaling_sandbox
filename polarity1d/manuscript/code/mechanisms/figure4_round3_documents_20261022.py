"""Scientific documentation bound to the portable Figure 4 round-3 records."""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path
from datetime import datetime, timezone

from figure4_round3_render_20261022 import (OUT, ROOT, REPO, MODELS, NAMES,
    read_csv, n, finite, primary, first, equivalent_size, write_json, sha)


def fmt(value):
    if value is None:return 'unavailable'
    try:value=float(value)
    except (TypeError,ValueError):return 'unavailable'
    if math.isnan(value):return 'unavailable'
    if math.isinf(value):return '−∞' if value<0 else '∞'
    return f"{value:.4g}"


def span(rows, key):
    values=[n(row,key) for row in rows if finite(n(row,key))]
    return f"{fmt(min(values))}–{fmt(max(values))}" if values else "unavailable"


def create(out=OUT):
    data,meta=out/'supporting/data',out/'supporting/metadata'
    path=meta/'analysis_receipt.json'
    receipt=json.loads(path.read_text()) if path.exists() else {}
    provisional=not receipt.get('completed_campaign',False)
    state='provisional_round2_pending_new_campaign' if provisional else 'complete_compact_campaign'
    alpha=primary(read_csv(data/'local_alpha.csv'))
    folds=primary(read_csv(data/'noise_tolerance_folds.csv'))
    tols=primary(read_csv(data/'noise_tolerance.csv'))
    lag=primary(read_csv(data/'condition_lag_metrics.csv'))
    single=primary(read_csv(data/'analysis_lag_metrics.csv'))
    end=primary(read_csv(data/'end_to_end_alpha.csv'))
    decomp=[r for r in primary(read_csv(data/'drift_decomposition.csv')) if n(r,'flicker_frames',default=5)==5 and math.isclose(n(r,'MSD_range_fraction',default=.25),.25)]
    omissions=read_csv(data/'omitted_estimates.csv')
    selection=read_csv(data/'model_D_selection.csv')
    null_agreement=read_csv(data/'null_B_analytic_agreement.csv')
    null_receipt_path=meta/'null_analysis_receipt.json'
    null_receipt=json.loads(null_receipt_path.read_text()) if null_receipt_path.exists() else {}
    seed_counts=sorted({int(n(r,'n_seeds')) for r in lag if finite(n(r,'n_seeds'))})
    planned_path=meta/'campaign_manifest.json'
    planned=json.loads(planned_path.read_text()).get('conditions',[]) if planned_path.exists() else []
    observed_keys={(r.get('model'),n(r,'C_label'),n(r,'f')) for r in single}
    pending_design=[r for r in planned if (r.get('model'),n(r,'size'),n(r,'fraction')) not in observed_keys]
    unresolved=[r for r in alpha if not finite(n(r,'alpha'))]
    tracked_count=sum(finite(n(r,'Ds')) for r in decomp)
    predicted_count=sum(finite(n(r,'kpred')) and n(r,'kobs')>0 for r in decomp)
    local_count=sum(finite(n(r,'alpha')) for r in alpha)
    n_by_model={m:sum(r.get('model')==m and finite(n(r,'alpha')) for r in alpha) for m in MODELS}
    main_lag=lambda m:200 if m=='spring_ring' else 2000
    main_tols=[r for r in tols if n(r,'reference_lag')==main_lag(r.get('model'))]
    main_folds=[r for r in folds if n(r,'reference_lag')==main_lag(r.get('model'))]
    fold_numbers=[]
    for model in MODELS:
        row=next((r for r in main_folds if r.get('model')==model and math.isclose(equivalent_size(r),60)),{})
        fold_numbers.append(f"{NAMES[model]} {fmt(n(row,'fhalf_fold','fold'))}")
    prefix='The model panels are a provisional reanalysis of retained round-2 histories; the newly requested 96–384-seed campaign has not completed. ' if provisional else ''
    d_selected=[r for r in selection if finite(n(r,'f'))]
    d_qualified=0
    for selected in d_selected:
        model=selected.get('model');ref=n(selected,'reference_lag',default=200 if model=='spring_ring' else 2000)
        refrows=[r for r in lag if r.get('model')==model and math.isclose(n(r,'f'),n(selected,'f')) and math.isclose(n(r,'lag'),ref)]
        small=next((r for r in refrows if math.isclose(equivalent_size(r),20)),{})
        large=next((r for r in refrows if math.isclose(equivalent_size(r),60)),{})
        d_qualified+=n(small,'C')<.4 and n(large,'C')>.6
    if not d_selected:
        d_observation='No model input has an available selection under the declared correlation or noise-tolerance-midpoint rules; exact reasons remain in the selection CSV.'
    elif d_qualified==len(d_selected):
        d_observation='At every selected input, small-ring reference-lag correlation is below 0.4 and large-ring correlation exceeds 0.6; the enlarged rings retain heading longer.'
    else:
        d_observation=f'{d_qualified}/{len(d_selected)} selected inputs satisfy the small-ring correlation below 0.4 and large-ring correlation above 0.6 criterion. The other selections use the declared nearest log-midpoint noise-tolerance fallback; exact reasons remain in the selection CSV.'
    ratios=[n(r,'kpred')/n(r,'kobs') for r in decomp if finite(n(r,'kpred')) and n(r,'kobs')>0]
    ratio_observation=(f'The prediction/observation ratio spans {min(ratios):.2g}–{max(ratios):.2g}, exposing conditions where this closure misses measured decay.' if ratios else 'The tracked decomposition remains unavailable where retained records do not support its ingredients.')
    paragraphs=[
        ('A', 'Simulated drifting and relocating peaks on small and large rings. Grayscale kymograph widths follow circumference and time increases upward; each pair shares its declared seed. Ensemble heading autocorrelations below use small-ring solid and large-ring dotted curves. Fixed physical drift, same-rate relocation and perimeter-rate relocation give size exponents 2, 0 and −1, respectively.'),
        ('B', 'Simulated mixtures illustrate how relocation erodes a geometric size advantage. Input fraction is f; the illustrative dependences are Dₛ = 8.8f² and λ = 135f⁴, with perimeter-rate relocation multiplying λ by C/40. The decay rate is k(ℓ) = −ln C(ℓ)/ℓ; α is minus the local slope of ln k against ln circumference. Noise tolerance f½ is the first downward crossing of correlation 0.5 at the fixed reference lag. Gray shades encode increasing circumference. Line styles distinguish mixtures in the tolerance-fold panel; its dashed guide is proportional to circumference.'),
        ('C', ('Provisional model estimates reuse round-2 histories; the new campaign is pending. ' if provisional else '')+'Local exponents use model hues and light-to-dark circumference ramps. Markers show designed inputs, open markers identify forced conditions, vertical intervals show pointwise 95% uncertainty, and outward triangles mark off-scale values or tails. The right panel shows tolerance folds with model-specific markers; Spring circumference is multiplied by 20. '+f'{local_count}/{len(alpha)} neighborhoods resolve; Goryachev contributes {n_by_model["goryachev"]}. Largest-size tolerance folds are '+', '.join(fold_numbers)+'.'),
        ('D', 'Matched first-replicate model histories use the selected designed inputs, with physical widths and shared grayscale activity scales within each model. Black outside ticks identify persistent dominant-region replacements. The ensemble curves below include end-to-end exponents and 95% intervals where estimable. '+d_observation),
        ('E', 'Columns distinguish models; rows show physical drift diffusivity, effective replacement rate, predicted/observed decay-rate ratio and multiple-front fraction. Circumference uses the same ramps. The prediction is (2π/C)²Dₛ + λ(1 − mean cos Δθ), with unity as the ratio reference. '+ratio_observation+' Undefined values and empirical zeros on logarithmic axes are omitted with reasons recorded. Multiple-front fractions remain descriptive. Time and length retain each model’s native units throughout.'),
    ]
    title='Null models for size-dependent heading persistence'
    legend='# '+title+'\n\n'+'\n\n'.join(f'**{letter}.** {text}' for letter,text in paragraphs)+'\n'
    (meta/'main_legend.md').write_text(legend)

    text=['# Figure 4 round 3: scientific record','',
          '## Unexpected findings and unresolved conditions','',
          '**Polarity normalization differs from the historical scalar gate.** Retained round-2 trace P was normalized by the initial first moment, whereas this brief explicitly defines P = M₁/M₀. The present analysis reconstructs P = M₁/M₀ from retained moments. The changed magnitude gate can alter endpoint eligibility, coverage, selected lags, local exponents and the displayed points relative to the appendix sketch. Historical round-2 figures and S0 records are preserved unchanged.','',
          '**Switch counting and displacement availability require separate masks.** Earlier round-3 preview code removed an accepted dominant identity whenever its centroid was absent, then counted switches only across the remaining adjacent segments. This added an unintended switch-count criterion. In the checked new-campaign Wave-pinning C=20, f=0.3, seed-0 history, track 162 is replaced by track 166 for two frames (times 31,425 and 31,430), then track 167 persists for at least five frames from time 31,435. Every frame contains observed regions with finite headings; only the previously accepted track 162 disappears at time 31,430. The two-frame flicker is ignored, so the persistent 162→167 replacement counts. Its missing centroid interrupts displacement estimation without cancelling the observed switch. Final tracked summaries use persistence-filtered dominance for switch counts, valid exposure and median tenure, with a separate centroid-valid mask for displacement. Genuine missing-dominance intervals and seed boundaries remain breaks. Historical provisional preview records retain their earlier method and must not be substituted for the rebuilt final tracked summaries.','',
          ('**Provisional model evidence.** The new from-scratch compact campaign has not completed. All displayed model measurements explicitly reuse retained round-2 observations, with the new round-3 estimators and required polarity normalization. Null simulations are newly generated and separately identified. This is a complete review figure with provisional model evidence; it does not claim new high-replication simulation results.' if provisional else '**Completed compact campaign.** Model measurements use the completed from-scratch compact records identified in the campaign manifest.'),'',
          f'{len(unresolved)} of {len(alpha)} local size neighborhoods have no reported exponent. Every failed neighborhood remains a row in `local_alpha.csv` and `omitted_estimates.csv`; the omission reasons identify missing members, inadequate endpoint coverage, or no shared lag with all neighborhood correlations strictly between 0.03 and 0.995. A nearly empty Goryachev facet is retained as an estimation result. No absent size or input is interpolated into a neighborhood.','',
          f'Available local exponents span {span(alpha,"alpha")}. The plotting range is −1.5 to 4: outward triangles preserve off-scale estimates and interval tails. Wild exponents are retained, not clipped in the data or smoothed into the main figure.','',
          f'The retained model ensembles contain {", ".join(map(str,seed_counts)) or "unavailable"} independent histories per condition. A 2,000-resample bootstrap does not increase the number of independent histories. Sparse original ensembles remain identified in the compact files.','',
          f'Tracked diffusivity is finite in {tracked_count} exported conditions; a finite predicted/observed ratio is available in {predicted_count}. Zero counted switches give an empirical switch contribution of zero; a log axis omits that zero rather than displaying a positive substitute. `jump_no_event_upper95 = 3 / valid_time` is a Poisson zero-count reference bound under a constant-hazard assumption, not an assumption-free confidence bound for these models. Whole-seed bootstrap intervals are separate and can be [0, 0] when no event was observed. Insufficient tenures do not imply a fitted diffusivity, and undefined measures are never imputed. Every unavailable decomposition quantity has a status/reason in `drift_decomposition.csv` or `omitted_estimates.csv`.','',
          'Null simulations can be too stable at weak input or decorrelate before the shortest declared lag at strong input. Such cells have no eligible local exponent; analytic formulas are used only for validation, never as replacements in the main panels. The simulated null analytic-check diagnostics and precision shortfalls are retained even if they are unfavorable.','']
    if planned:text.extend([f'The frozen new campaign contains {len(planned)} model–size–input conditions. The current source contains {len(observed_keys)} observed conditions, leaving {len(pending_design)} planned cells without source observations. The local-neighborhood census uses the represented input union and can include missing neighborhood members; it must not be confused with completed new-campaign conditions. [Every planned condition](campaign_manifest.json).',''])
    precision_unmet=[r for r in single if n(r,'SE_log_k')>.07]
    undefined_rate=[r for r in single if not finite(n(r,'k'))]
    low_boot=[r for r in single if finite(n(r,'k')) and n(r,'bootstrap_valid_fraction')<.975]
    capped=[r for r in precision_unmet if n(r,'n_seeds')>=384]
    text.extend([f'Model rate-precision qualification: {len(precision_unmet)} observed analysis-lag rates have SE(ln k) > 0.07; {len(capped)} are at the 384-history cap. {len(undefined_rate)} conditions have no eligible single-condition rate, and {len(low_boot)} finite point estimates have fewer than 97.5% valid bootstrap rate draws. '+('These are retained sparse historical ensembles; the new campaign precision target has not yet been assessed on completed new records. ' if provisional else '')+'[Every condition, actual seed count, selected lag, SE and bootstrap-valid fraction](../data/analysis_lag_metrics.csv).',''])
    if capped:text.extend(['Conditions at the seed cap that still exceed the finite-rate precision target: '+', '.join(f'{NAMES.get(r.get("model"),r.get("model"))} C={r.get("C_label")}, f={r.get("f")}, SE={fmt(n(r,"SE_log_k"))}' for r in capped)+'.',''])
    text.extend(['**The null B curves cannot cover the full analytic sketch under the declared lag rule.** At weak input, correlation can remain above 0.995 at the maximum lag of 16,000; at strong input, it can already be below 0.03 at the minimum lag of 250. The permitted local-neighborhood lag window therefore truncates the simulated curves, including possible asymptotic regimes. More independent histories do not repair this time-resolution limit. The shared-estimator and no-imputation requirements take precedence over drawing a complete analytic-looking curve.',''])
    text.extend(['[Every omitted estimate, with exact model/size/input and reason](../data/omitted_estimates.csv) · [Complete local-neighborhood census, including unselected lags](../data/local_alpha.csv) · [Every observed lag, correlation, coverage and interval](../data/condition_lag_metrics.csv). These inventories include each failed condition; no failed condition is removed from its designed grid.',''])
    for model in MODELS:
        failures=[r for r in unresolved if r.get('model')==model]
        reasons=sorted({str(first(r,'reason','status',default='No shared eligible lag')) for r in failures})
        for reason in reasons:
            selected=[r for r in failures if str(first(r,'reason','status',default='No shared eligible lag'))==reason]
            keys=', '.join(f"({r.get('C_label')}, {r.get('f')})" for r in selected)
            text.extend([f'{NAMES[model]} — {reason}; exact (circumference label, input) pairs: {keys}.',''])
    outliers=read_csv(data/'display_outliers.csv')
    if outliers:
        text.extend(['Off-scale finite point estimates are retained with outward edge triangles. [Exact outlier records and intervals](../data/display_outliers.csv):',''])
        for row in outliers:
            text.extend([f'- Panel {row.get("panel")}, {NAMES.get(row.get("model"),row.get("model"))}, C={row.get("C_label")}, f={row.get("f")}: '+(f'α={fmt(n(row,"alpha"))}, 95% interval [{fmt(n(row,"alpha_ci_lower"))}, {fmt(n(row,"alpha_ci_upper"))}]' if row.get('panel')=='C' else f'predicted/observed ratio={fmt(n(row,"value"))}')+'.'])
        text.append('')
    finite_checks=[r for r in null_agreement if finite(n(r,'relative_deviation'))]
    if finite_checks:
        worst=max(finite_checks,key=lambda r:abs(n(r,'relative_deviation')))
        counts=[n(r,'n_seeds') for r in null_agreement if finite(n(r,'n_seeds'))]
        unmet=[r for r in null_agreement if r.get('status')=='precision_met' and n(r,'SE_log_k')>.01]
        text.extend([f'Null B maximum absolute relative rate deviation from its analytic check is {abs(n(worst,"relative_deviation")):.4g} ({worst.get("model")}, C={worst.get("C_label")}, f={worst.get("f")}, selected lag={worst.get("selected_lag")}); signed deviation {n(worst,"relative_deviation"):.4g}. Actual replicate counts span {min(counts):g}–{max(counts):g}; {len(unmet)} reported precision-met cells exceed SE(ln k)=0.01. [Every null condition and check](../data/null_B_analytic_agreement.csv).',''])
    null_a=read_csv(data/'null_A_acf.csv')
    deviations=[abs(n(r,'C')-n(r,'analytic_C')) for r in null_a if finite(n(r,'C')) and finite(n(r,'analytic_C'))]
    if deviations:text.extend([f'Null A maximum absolute simulated-minus-analytic autocorrelation deviation is {max(deviations):.4g} across the displayed endpoint/lag records. [All autocorrelations and analytic checks](../data/null_A_acf.csv).',''])
    negatives=[r for r in null_a if n(r,'C')<0]
    if negatives:
        lowest=min(negatives,key=lambda r:n(r,'C'))
        text.extend([f'Null A has {len(negatives)} negative finite-ensemble autocorrelation points; the minimum is {n(lowest,"C"):.5g} for {lowest.get("model")}, C={lowest.get("C_label")}, lag={lowest.get("lag")}. They are retained on the displayed slightly negative ordinate and in CSV. At long lags, a 1,000-history ensemble can fluctuate around a small positive expected correlation; no history or seed was replaced to suppress this behavior.',''])
    spring_coordinate=read_csv(data/'spring_coordinate_sensitivity.csv')
    spring_ratios=[n(r,'literal_over_arc_Ds') for r in spring_coordinate if finite(n(r,'literal_over_arc_Ds'))]
    if spring_ratios:
        effect=max(abs(v-1) for v in spring_ratios)*100
        text.extend([f'Spring physical-coordinate sensitivity: {len(spring_ratios)}/{len(spring_coordinate)} condition ratios are finite. Literal fixed-origin coordinate / intrinsic-arc diffusivity spans {min(spring_ratios):.8g}–{max(spring_ratios):.8g}, a maximum point-estimate change of {effect:.4g}%. Both definitions use identical tenures and the fixed primary MSD range; all available primary intrinsic-arc estimates are independently reproduced. The alternative is a sensitivity check and does not replace or tune the primary ingredient. [Every ratio, interval and coordinate convention](../data/spring_coordinate_sensitivity.csv).',''])
    text.extend(['For B, the literal Bernoulli probability λΔt exceeds one at strong inputs. The simulator therefore uses the exact sampled Brownian-plus-Poisson-reset transition, with reset probability 1 − exp(−λΔt); every reset is uniform. A retains its requested Bernoulli-step experiment and validates against the corresponding discrete-step rate. This protocol qualification is recorded, rather than silently clipping a probability.','',
                 'Provisional tracking cadence varies with the retained observation cohort: chemical records include cadences 5 and 40, and Spring records include 0.5 and 4. A five-frame hold therefore spans 25 or 200 native chemical time units, and 2.5 or 20 Spring time units. The new compact campaign uses the scalar analysis cadence (5 and 0.5), giving hold times 25 and 2.5. `source_cadence` and `flicker_hold_time` identify the exact value for every tracked condition in [the decomposition CSV](../data/drift_decomposition.csv). This temporal-resolution difference qualifies provisional ingredients and switch ticks; no mixed cadence is pooled within a condition.',''])
    for limitation in receipt.get('limitations',[]):
        text.extend([str(limitation),''])
    text.extend(['The campaign failure/status receipt records scheduler availability and the exact ready submission command. No old Figure 4 native archive is deleted before the complete new compact campaign passes the stated retention gate. Round-2 S0/H1/H2/H3/A6 history remains unchanged; the current analysis has no hypothesis gates.','',
                 '## Panel key numbers','',
                 'A: three point-null families; the displayed circumferences are 20 and 60, with fixed physical width 1.5, drift diffusivity 0.002, jump rate 1/7000, total duration 24,000 and step 20. The perimeter-rate variant multiplies the jump rate by C/40. The requested autocorrelation ensemble is 1,000 independent histories per family and size; actual simulated counts and uncertainty are recorded in the null receipt.','',
                 'B: 41 logarithmically spaced designed inputs from 10⁻⁴ to 1 and five fixed circumferences (20, 30, 40, 50, 60). The target simulation precision is SE(ln k) ≤ 0.01. Actual stopping counts and any unresolved conditions are recorded; they must not be inferred from the number of bootstrap draws.','',
                 f'C: {local_count}/{len(alpha)} local estimates resolved: '+', '.join(f'{NAMES[m]} {n_by_model[m]}' for m in MODELS)+'. Largest/smallest fixed-reference tolerance folds: '+', '.join(fold_numbers)+'. All folds, including censored and undefined values, remain in `noise_tolerance_folds.csv`.',''])
    for row in selection:
        text.extend([f'D — {NAMES.get(row.get("model"),row.get("model"))}: f = {row.get("f")}; '+str(first(row,'selection_reason','reason',default='selection in compact record'))+'. Requested/actual window: '+str(first(row,'window_requested',default='recorded in index'))+' / '+str(first(row,'window_actual',default='recorded in index'))+'.',''])
    text.extend([f'E: {len(decomp)} condition records; {tracked_count} finite tracked-drift estimates; {predicted_count} finite predicted/observed comparisons. The panel retains every defined descriptive multiple-front fraction.','',
                 '## Methods and estimands','',
                 'Length represents WGD-associated enlargement. The empirical motivation includes greater migration persistence; observed migration summaries, conditional OU velocity-memory estimates, simulated activity-heading memory and one-dimensional polarity-reversal readouts are distinct evidence. The present analysis measures simulated activity heading rather than cell displacement.','',
                 'The complex activity vector is the first spatial moment divided by total activity. Heading autocorrelation averages the cosine of the angular change over saved endpoint pairs with polarity magnitude at least 0.1 at both endpoints. Endpoint coverage is the eligible-pair count divided by all possible endpoint pairs within the declared observation window.','',
                 '$$\nC(\\ell)=\\left\\langle\\cos[\\theta(t+\\ell)-\\theta(t)]\\right\\rangle_{|P(t)|,|P(t+\\ell)|\\geq0.1},\\qquad k(\\ell)=-\\frac{\\ln C(\\ell)}{\\ell}.\n\\tag{1}\n$$','',
                 'Target lags use square-root-of-two steps from 250 to 16,000 for chemical models and 25 to 1,600 for Spring, clipped to the available observation window. The compact lag table records both target and actual saved-endpoint lags; the analysis receipt states its endpoint-rounding rule. No temporal interpolation is hidden.','',
                 'Each local exponent uses the complete fixed five-size grid: the focal size and its immediate neighbors, with two members at the ends. Actual mean circumference is used in the fit. A shared lag requires correlation strictly between 0.03 and 0.995 and coverage at least 0.5 in every member. Among eligible lags, the selected lag minimizes the worst distance of ln(−ln C) from ln(ln 2). Lag selection is held fixed in bootstrap draws.','',
                 '$$\n\\alpha=-\\frac{\\operatorname{Cov}(\\ln \\bar C,\\ln k)}{\\operatorname{Var}(\\ln \\bar C)}.\n\\tag{2}\n$$','',
                 'The single-condition analysis lag is the eligible correlation closest to 0.5 in absolute correlation distance. Noise tolerance uses the first downward crossing of absolute correlation 0.5 at reference lag 2,000 for chemical models and 200 for Spring, interpolated linearly in logarithmic input. Forced conditions are included and encoded by open markers. Reference-lag sensitivities (1,000/4,000, or 100/400 for Spring) remain in the compact records. Missing crossings, left/right censoring and gaps are retained explicitly.','',
                 'Confidence intervals are pointwise 95% whole-seed bootstrap intervals from 2,000 resamples. Every time point, lag and tracked event belonging to a resampled independent history is resampled together. Size ensembles are independent. Undefined logarithms or inadequate valid-draw fractions remain undefined or unbounded according to the estimator receipt.','',
                 'A dominance tenure is a maximal interval with a constant dominant largest-amount track ID after requiring a replacement to persist for five frames; shorter flickers retain the previously accepted identity. Switches, valid exposure and the median tenure use these persistence-filtered identities before checking centroid availability. Exposure includes saved time intervals whose two endpoints have nonmissing filtered identities. A genuinely missing dominant ID resets the accepted identity; neither an event nor exposure crosses that gap or a seed boundary. A missing centroid of the accepted identity splits only the displacement segment, without deleting that identity from the switch/exposure record or shortening the dominance tenure used to choose the MSD range.','',
                 'The primary physical coordinate integrates wrapped angular centroid increments multiplied by midpoint circumference, within each continuous centroid-valid segment of a dominance tenure. No missing centroid is interpolated or replaced by the short-lived instantaneous winner. A uniform dilation at fixed material angle therefore contributes no heading-changing drift. This is an intrinsic arclength increment convention, rather than the literal product of unwrapped angle and instantaneous circumference. Spring coordinate sensitivity is retained separately where available.','',
                 '$$\n\\Delta s_i=\\frac{C_i+C_{i+1}}{4\\pi}\\,\\operatorname{wrap}(\\theta_{i+1}-\\theta_i),\\qquad s_j=\\sum_{i<j}\\Delta s_i.\n\\tag{3}\n$$','',
                 'The primary diffusivity fits MSD = 2Dₛτ through the origin through one quarter of the median dominance tenure, pooling only displacement pairs within the continuous centroid-valid segments. The jump rate is counted persistent replacements divided by valid dominance exposure. Mean jump cosine uses the whole-activity P heading averaged over the five saved frames immediately before and after each counted replacement; the windows can include another identity but cannot cross an invalid-P interval or a seed boundary. Both windows require every polarity magnitude to be at least 0.1. Switch counting does not depend on these heading windows: if any counted switch lacks a valid jump-angle score, the effective jump rate remains undefined rather than discarding that switch or imputing an angle.','',
                 '$$\n\\lambda_{\\mathrm{eff}}=\\lambda\\left(1-\\langle\\cos\\Delta\\theta\\rangle\\right),\\qquad k_{\\mathrm{pred}}=\\left(\\frac{2\\pi}{\\bar C}\\right)^2D_s+\\lambda_{\\mathrm{eff}}.\n\\tag{4}\n$$','',
                 'The detector-free mixture uses wrapped-Gaussian angular increments with variance 2Dθδ plus uniform relocations of weight 1 − exp(−λδ), at the shortest cadence and twice that cadence. Physical diffusivity is DθC²/(4π²). Declared tracked sensitivities use persistence thresholds of 2 and 10 frames, and MSD fit limits of one eighth and one half of median tenure. The supplementary PDF shows the available comparisons and explicitly identifies diagnostics that await compact records.','',
                 '## Rebuild, portable records and provenance','',
                 'From the repository root, run:','',
                 '```bash\nconda run --no-capture-output -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figure4_round3_20261022.py\n```','',
                 'The default rebuild reads compact records only. It does not start a campaign, read native full-field archives or delete old outputs. Use the explicitly named preparation options to regenerate scalar/track summaries from retained round-2 evidence or to rerun null simulations. The manifest and hashes bind each figure to its actual source.','',
                 '- `local_alpha.csv`, `end_to_end_alpha.csv` and `analysis_lag_metrics.csv`: estimates, selected lags, intervals, valid fractions and statuses.\n- `condition_lag_metrics.csv`: every saved correlation/rate curve point, coverage and seed count.\n- `noise_tolerance.csv` and `noise_tolerance_folds.csv`: every reference lag, crossing/fold, interval and censoring status.\n- `drift_decomposition.csv`, within-tenure and sensitivity CSVs: directly measured ingredients and qualifications.\n- `null_A_*` and `null_B_*`: newly simulated null records, distinct from model evidence.\n- `plotted_*.csv`: exact compact curve/point sources consumed by the main renderer.\n- `omitted_estimates.csv`: unavailable/off-domain estimates with their original reasons.\n- `model_D_selection.csv` and the kymograph index: deterministic displayed-input, replicate and time-window selection.\n- `analysis_receipt.json`, `null_analysis_receipt.json`, `render_receipt.json`, campaign manifest/status and this documentation: source hashes, estimator rules and provenance.','',
                 'The main PDF is 180 mm wide and 500 mm high. It contains no overarching title or subtitle. Axes retain only left/bottom spines, time increases upward, model facets remain black, and each model has a local circumference key. No table appears in the current figure or its new supplementary diagnostic pages.',''])
    e_details=[]
    for model in MODELS:
        selected=[r for r in decomp if r.get('model')==model and finite(n(r,'kpred')) and n(r,'kobs')>0]
        if not selected:
            e_details.extend([f'E — {NAMES[model]}: no finite decomposition comparison.','']);continue
        low=min(selected,key=lambda r:n(r,'kpred')/n(r,'kobs'));high=max(selected,key=lambda r:n(r,'kpred')/n(r,'kobs'))
        e_details.extend([f'E — {NAMES[model]}: {len(selected)} finite ratios; minimum {n(low,"kpred")/n(low,"kobs"):.4g} at C={low.get("C_label")}, f={low.get("f")}; maximum {n(high,"kpred")/n(high,"kobs"):.4g} at C={high.get("C_label")}, f={high.get("f")}.',''])
    marker=text.index('## Methods and estimands')
    text[marker:marker]=e_details
    (meta/'scientific_readme.md').write_text('\n'.join(text))
    render=json.loads((meta/'render_receipt.json').read_text()) if (meta/'render_receipt.json').exists() else {}
    campaign=list(meta.glob('*campaign*manifest*.json'))
    sources=[dict(label=path.stem.replace('_',' '),path=str(path.relative_to(REPO)),sha256=sha(path)) for path in [meta/'analysis_receipt.json',meta/'null_analysis_receipt.json',*campaign] if path.exists()]
    summary=dict(local_alpha_resolved=local_count,local_alpha_total=len(alpha),local_alpha_resolved_by_model=n_by_model,tracked_diffusivity_available=tracked_count,prediction_comparisons=predicted_count,seed_counts=seed_counts,bootstrap_resamples=2000,width_mm=180,height_mm=500)
    decomposition_by_model={}
    for model in MODELS:
        values=[n(r,'kpred')/n(r,'kobs') for r in decomp if r.get('model')==model and finite(n(r,'kpred')) and n(r,'kobs')>0]
        decomposition_by_model[model]=dict(n=len(values),minimum=min(values) if values else None,median=float(statistics.median(values)) if values else None,maximum=max(values) if values else None)
    summary['decomposition_by_model']=decomposition_by_model
    ratio_text=', '.join(f'{NAMES[m]} {fmt(v["minimum"])}–{fmt(v["maximum"])} (n={v["n"]})' for m,v in decomposition_by_model.items())
    legend_payload=dict(id='figure4',title=title,role='main',status=state,
        assets=[dict(path=str((out/'figure4_main.png').relative_to(REPO)),caption='',page_size=[180/25.4,500/25.4],minimum_printed_width_inches=180/25.4,native_pdf_page=str((out/'figure4_main.pdf').relative_to(REPO)),page_title='')],
        legend=[f'{letter}. {paragraph}' for letter,paragraph in paragraphs],
        narrative=[prefix+'Simple drift and relocation nulls expose distinct geometric size exponents; unchanged polarity models are measured with the same local estimators, and the directly measured decomposition retains all failures.',
                   'Largest/smallest fixed-reference noise-tolerance folds are '+', '.join(fold_numbers)+f'. Local estimation resolves {local_count}/{len(alpha)} fixed neighborhoods; the fixed magnitude gate uses M₁/M₀.',
                   'Available tracked predicted/observed decay-rate ratios span '+ratio_text+'. These condition-specific comparisons describe the closure; they do not establish multiple fronts as a causal explanation.'],
        measurements=['Measure activity-heading autocorrelation, centroid motion within persistent active regions, and persistent region replacements in WGD+ versus WGD− cells. Compare enlargement-associated changes in drift and relocation without equating activity heading to migration persistence.'],
        limitations=([prefix.strip()] if provisional else [])+['Native model clocks are not calibrated across models. Whole-seed intervals are pointwise. Missing eligible lags and tracked ingredients remain undefined; multiple-front fractions are descriptive. The required M₁/M₀ gate differs from historical round-2 normalization.'],sources=sources,
        main_pdf=str((out/'figure4_main.pdf').relative_to(REPO)),supplement_pdf=str((out/'figure4_supplement.pdf').relative_to(REPO)),reader=str((out/'index.html').relative_to(REPO)),code_commit=render.get('code_commit'),campaign_completed=not provisional)
    legend_payload['legend_panels']=[dict(panel=letter,text=text) for letter,text in paragraphs]
    write_json(out/'legend.json',legend_payload)
    write_json(meta/'panel_summary.json',summary)
    state_text='Provisional round-2 model evidence; newly simulated nulls. New compact campaign pending.' if provisional else 'Completed compact campaign.'
    html=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title><style>body{{font:16px/1.55 system-ui,sans-serif;max-width:1050px;margin:2rem auto;padding:0 1rem;color:#161616}}img{{width:100%;height:auto}}a{{color:#17619c}}.status{{padding:1rem;background:#f5f5f5;border-left:3px solid #555}}time{{font-size:.85rem}}h1{{font-size:1.6rem}}</style><h1>{title}</h1><p class="status">{state_text}</p><p><time>Model data generated: {receipt.get('generated_utc','unavailable')}. Null data generated: {null_receipt.get('generated_utc','unavailable')}. Plate rebuilt: {render.get('generated_utc','unavailable')}.</time></p><p><a href="figure4_main.pdf">Main PDF</a> · <a href="figure4_supplement.pdf">Complete new supplement</a> · <a href="supporting/metadata/main_legend.md">Legend</a> · <a href="supporting/metadata/scientific_readme.md">Scientific record and failures</a> · <a href="supporting/metadata/render_receipt.json">Hashed render receipt</a></p><img src="figure4_main.png" alt="Figure 4 null simulations, local size exponents, matched model histories and drift–jump decomposition"><p>{prefix}</p></html>'''
    (out/'index.html').write_text(html)
    return legend_payload


if __name__=='__main__':
    create()
