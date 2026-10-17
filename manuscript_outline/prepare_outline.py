"""Assemble the manuscript outline from reviewed local evidence; no simulation."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'manuscript_outline'


def section(identifier, title, role, assets=(), legend=(), narrative=(),
            measurements=(), limitations=(), sources=()):
    return dict(id=identifier, title=title, role=role, assets=list(assets),
                legend=list(legend), narrative=list(narrative),
                measurements=list(measurements), limitations=list(limitations),
                sources=[dict(label=label, path=path) for label, path in sources])


def foundation():
    return section(
        'figure1', 'Figure 1. Enlargement changes the cue duration required for lasting polarity replacement',
        'Comparative theoretical foundation: five established models, three finite challenges, and two distinct scalar replacement gates.',
        assets=[dict(path='polarity1d/manuscript/figures/main_polarity.png',
                     caption='Current canonical comparison, 10 October Spearman revision; corrected independent-endpoint deBelly.')],
        legend=[
            'A, Example wave-pinning preparation from a homogeneous state using a finite left patch. The cue ends at 0.773 native time units, whereas the subsequent basal relaxation establishes the polarized profile later (reported establishment arrival about 116). Coloured profiles and orientation samples distinguish input duration from autonomous development. The grey late interval is verification, not additional stimulation. Profiles are scaled for display.',
            'B, Three rightward finite challenges restart the identical qualifying nominal-amplitude, relative-20% basal checkpoint at each model and length: fixed physical support 3 (Spring 0.2 native length units), relative support 20%, and a linear ramp with constant physical slope. Patch peaks remain fixed within each model; a fixed-slope ramp has a larger peak in a longer domain. Goryachev patches use the independently selected, length-independent multiplier K=64 (peak 6.4), which makes this a response-permissive comparison rather than a calibrated equal biochemical dose across models.',
            'C, Spearman rank correlation between length and the shortest tested passing pulse duration, separately for orientation replacement and strength replacement, across Wave-pinning, Goryachev, corrected deBelly, Holmes 3 and Spring. Each cell reports the correlation and included successful-length count. All audited scalar successes enter, including those failing the separate native-field tail criterion. At least four distinct successful lengths and nonconstant durations are required; censored and preparation-ineligible cases receive no invented time. The companion duration curves retain all lengths and field-status markers.',
            'After each positive pulse, input returns to the native basal value for 10,000 native time units. Orientation requires all 201 final-500 samples to have scalar orientation above 0.05 and contrast at least 0.05, with orientation range below 0.005. Strength additionally requires every tail orientation to be within 1% of the absolute initial leftward scalar orientation; it does not compare individual pole heights or amounts. The separate native-field tail flag requires pointwise range and projected 500-unit drift within 0.1% of frozen field scales. The reported duration excludes preparation and recovery. It is the shortest passing tested duration, with sampled fail/pass boundaries refined to 1%; narrower unobserved intervals are not excluded.',
            'deBelly uses independently commanded corrected endpoints, rather than the deprecated in-loop centered extension. Spring length and time retain native units and are not silently converted to micrometres or seconds. Model-specific axes and stimulus definitions prevent interpreting cross-model amplitudes or clocks as biologically calibrated equivalents.'
        ],
        narrative=[
            'Fixed-width and relative-width patches show perfectly increasing rank order among observed successes in Wave-pinning, Goryachev, Holmes 3 and Spring; enlargement often increases the finite stimulation required for a retained response.',
            'The response depends on the cue geometry and model. Spring has a negative constant-slope association (Spearman correlation -0.874 for both gates), while Holmes orientation under the constant-slope cue is nearly unassociated with length (0.047).',
            'Corrected deBelly curves are nonmonotone despite moderate positive rank correlations. The full duration curves and eligibility ranges are essential to the result.',
            'These are deterministic polarity-replacement experiments. They motivate, but do not measure, spontaneous polarity memory or migration persistence.'
        ],
        measurements=[
            'Mechanism-linked cellular follow-up would match local input support, amplitude and duration, track the original and challenged activity regions through stimulation and withdrawal, and distinguish orientation replacement from recovered activity strength. This validation is proposed and has no experimental results in this outline.'
        ],
        limitations=[
            'Correlations are conditional on observed scalar success, rather than an unbiased summary of all possible sizes or a universal scaling law. Ineligible starts are distinct from responses unobserved by the finite horizon.',
            'A retained scalar response can coexist with a changing native field. Only 13 of the 72 corrected deBelly scalar successes pass the separate native-field tail flag. Resolution checks from the earlier centered-deBelly campaign cannot qualify the corrected replacement thresholds.'
        ],
        sources=[('Three-test protocol and Spearman rule', 'polarity1d/manuscript/THREE_TESTS.md'),
                 ('Correlation table', 'polarity1d/manuscript/data/derived/length_correlations.csv'),
                 ('Corrected endpoint method', 'polarity1d/manuscript/DEBELLY_BOUNDARIES.md'),
                 ('Current duration curves', 'polarity1d/manuscript/figures/supplement_width_response.png')])


def migration():
    return section(
        'figure5', 'Figure 5. Published migration summaries support a persistence hypothesis with conditional OU interpretation',
        'Existing empirical support, with observed migration summaries kept separate from fitted velocity-memory parameters.',
        assets=[dict(path='wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU.png',
                     caption='Latest integrated refined-OU figure. The original image retains its source label “Figure 3”; its proposed manuscript number here is 5.')],
        legend=[
            'A, Published/digitized enlarged or polyploid-group to reference-group ratios for nine migration summaries, alongside primary Ornstein–Uhlenbeck (OU) fit predictions. Xuan contributes sampled mean speed and population mean-squared displacement (MSD) at the last shared lag, approximately 10.95 h; Hosny summaries are group medians. Mallin central-bar statistic is not specified by the source and retains the historical mean approximation. These distinct source summaries are not interchangeable measures of persistence.',
            'B, Full digitized Xuan population MSD traces and fitted ensemble OU predictions on logarithmic axes. The digitized traces cross descriptively near 2.32 h; the fitted OU curves cross near 1.11 h, so their crossings are not identical. Digitized curve samples are correlated raster reconstructions of published curves, not independently measured biological observations or a recovered experimental lag grid. Both full curves are fitting targets. The early-MSD comparison overlaps these targets and is not independent validation.',
            'C, Representative Xuan primary-fit model-generated tracks, with sampled speed, total path length and directness (net displacement divided by path length). Both groups use identical spatial scales. The tracks are simulations, not published individual-cell trajectories; a slower sampled trajectory can nevertheless have greater endpoint directness. Dashed chords illustrate the contribution of finite-time displacements.',
            'D, With reference-group parameters fixed, vary enlarged-group persistence and velocity-noise ratios. Separate surfaces show change in the descriptive data-misfit score and in the total regularized objective; the selected fit and equal-group reference are marked. The heuristic score is not a calibrated likelihood, and the surfaces are not confidence regions.',
            'E, Primary paired persistence and velocity-noise estimates for Xuan 2018, Mallin 2023, Jemaà 2017, Jemaà 2023 and the Hosny hypoxic baseline. F, Enlarged/reference persistence folds under the primary 10-minute prior, 5- and 20-minute alternatives and stronger/weaker pooling, with the pooled geometric mean. Priors and pooling define conditional sensitivity scenarios, not uncertainty intervals. The Jemaà comparisons constrain speed alone, leaving persistence/noise tradeoffs dependent on hierarchy and priors.',
            'The fitted population contrasts span stable tetraploids, giant polyploids and treatment-induced polyaneuploid cancer cells (PACCs). The inherited WGD+/WGD− labels in the artwork are comparison shorthand and do not establish matched WGD-only interventions in every study. Hosny observations derive from the recovered preprint; equivalence to the journal version remains unresolved.'
        ],
        narrative=[
            'Xuan combines lower sampled speed (0.64-fold) with greater late-lag population MSD (2.08-fold); Mallin combines shorter path length (0.52-fold) with greater directness (3.12-fold). Persistence-related behaviour therefore deserves analysis separately from speed.',
            'Hosny baseline PACCs have greater median path length (1.77-fold) and directness (1.67-fold). Signed migration index is a separate directional summary.',
            'All five primary OU comparisons infer longer velocity memory. The pooled primary persistence fold is about 15.8, and remains above one across the declared prior/pooling scenarios; this is a conditional model result.',
            'The fit does not reproduce the complete Xuan curve shape within its raster envelope, and data-only sensitivity has rank 18 of 22. Optimization completion does not establish model adequacy or uniquely identified persistence times.'
        ],
        measurements=[
            'The new WGD assay should measure matched trajectories directly, with declared sampling and observation duration, and assess speed, displacement, directness and lag-dependent MSD before considering fitted persistence. An activity-polarity memory measurement would complement these quantities without equating activity and velocity.'
        ],
        limitations=[
            'Neither these heterogeneous group comparisons nor the OU fits establish that enlargement alone causes persistence. Cell identity, induction/treatment and acquisition differences require explicit treatment in the new experiments.',
            'The 13 October user feedback requested clearer observed-versus-model encoding and visible MSD curvature/exponent behaviour. The integrated figure shown here precedes that requested panel-B revision; this is recorded for the second editorial assessment.'
        ],
        sources=[('Refined OU analysis', 'wgd_langevin/refined_ou_20261003/README.md'),
                 ('Original panel legend', 'wgd_langevin/refined_ou_20261003/figures/main_03_migration_OU_caption.md'),
                 ('Observed folds', 'wgd_langevin/refined_ou_20261003/figures/main_03_observation_folds.csv'),
                 ('Fit results and adequacy', 'wgd_langevin/refined_ou_20261003/results/RESULTS.md'),
                 ('Hosny sensitivity supplement', 'wgd_langevin/refined_ou_20261003/figures/supp_04_Hosny_sensitivity.png')])


def planned_experiment():
    return section(
        'figure6', 'Figure 6. RESERVED: new WGD motility-assay results',
        'A substantial new experimental results section and main figure. No observations are available; the panel layout is a proposal, not an experimental result.',
        assets=[dict(path='manuscript_outline/assets/planned_WGD_results.png',
                     caption='Explicit placeholder. Every panel is pending acquisition, analysis and author review; no mock observations are displayed.')],
        legend=[
            'A, Proposed experimental context: the SUM-159 chemically induced WGD+/− pair is the agreed starting pilot with the Speros collaboration. Imaging should first be optimized for that pair before expansion to additional cell lines. Ploidy/group identity, physical size and morphology need to be documented so the intended comparison is clear.',
            'B, Reserved for representative experimental fields and trajectories under the selected acquisition conditions. The eventual panel should state sampling interval, usable duration, spatial scale, track eligibility and exclusions. Number of independent preparations and biological replicates remains to be specified.',
            'C, Reserved for empirical motility summaries: sampled speed, accumulated path length, endpoint displacement/directness and lag-dependent MSD. The final readouts, follow-up window and analysis plan are to be settled before interpreting the data; this placeholder does not imply that all listed measurements will be technically feasible.',
            'D, Reserved for the measured relationship between WGD-associated enlargement and migration behaviour, including size distributions and controls needed to interpret the group contrast. Any activity/protrusion measurements arising from the mechanism proposals require a separately agreed assay and should be distinguished from the basic motility pilot. Statistical models, uncertainty reporting and sample-size targets remain open.'
        ],
        narrative=[
            'This major results section will test whether the planned paired experimental system exhibits the migration behaviour motivating the theoretical work.',
            'Its conclusions must follow the acquired observations, including null or opposing results. Simulation output and literature fits cannot fill this gap.',
            'Agreement between increased polarity resistance or activity memory and migration persistence would support the proposed connection; causal attribution requires appropriate perturbations and controls.'
        ],
        measurements=[
            'Pending decisions include plating density, acquisition duration and cadence, field/condition throughput, tracking and exclusion rules, replicate structure, confirmation of WGD status, size quantification and the statistical analysis. The known throughput versus spatial/temporal resolution trade-off should be resolved by the pilot.',
            'Mechanism validation can subsequently examine incumbent extent, early challenger growth, local loss/rebound, site fates and global activity-polarity correlation under characterized perturbations. These activity experiments are proposals alongside the theoretical mechanisms and are separate from the reserved motility-results contribution.'
        ],
        limitations=['No pilot acquisition, sample size, effect estimate, significance claim or experimental outcome is supplied by the available evidence.'],
        sources=[('Author manuscript vision and supplied collaboration context', 'MANUSCRIPT_VISION.md')])


def base_document():
    return dict(
        title='Cell enlargement, polarity competition and migration persistence',
        subtitle='Vision-aligned manuscript outline with proposed main figures and detailed legends',
        date='17 October 2026', stage='Initial integrated checkpoint',
        overview=[
            'Working hypothesis: WGD-associated enlargement changes how established polarity is organized and replaced; these changes may provide a physical basis for greater migration persistence. The theoretical contribution tests size under unchanged local kinetic parameters. Length differs between model cells to represent WGD-associated enlargement, rather than elongation during a turn. Existing migration observations supply empirical motivation; new paired WGD motility results retain a substantial, explicitly reserved contribution.',
            'Proposed argument: first establish model- and stimulus-dependent size effects in finite polarity replacement, then explain them through local input response, recovery and site competition. Supported rings test circumference and separation, and stochastic ring experiments examine activity-polarity memory directly. Published migration summaries are evaluated separately, with OU interpretation clearly conditional. The new WGD assay closes the empirical gap when results arrive.',
            'This is an outline, not a completed manuscript or a claim that all figures are publication-final. The order of the mechanism and existing-migration components remains an editorial choice. Six proposed main figures are numbered here for navigation; source artwork may retain earlier working numbers. Mechanism-linked experimental validation is proposed, and Figure 6 remains wholly unmeasured.',
            'Suggested introduction: distinguish migration persistence from speed; describe why WGD-associated enlargement motivates a geometric polarity question; introduce established biochemical and mechanochemical models as complementary hypotheses; explain the need to connect mechanisms to measurable activity/protrusion behaviour without assuming a cell-velocity law.',
            'Suggested discussion: retain nonmonotone, null and opposite responses as informative model differences; separate deterministic cue resistance, noisy activity memory and velocity memory; identify concrete experiments; discuss finite horizons, stimulus geometry, conservation conventions and native times. The conclusions should remain conditional until the new motility and mechanism-validation experiments are obtained.'
        ],
        sections=[foundation(), migration(), planned_experiment()],
        supplementary=[
            dict(id='supplement1', title='S1. Full cue-duration response curves',
                 role='Essential curve context for the main comparison heatmap.',
                 assets=[dict(path='polarity1d/manuscript/figures/supplement_width_response.png', caption='Canonical corrected five-model duration curves; source working number Supplement 1.')],
                 legend=['Five model rows and three finite-cue columns retain the complete sampled length dependence. Coloured curves distinguish orientation from strength replacement. Open success markers retain scalar successes that fail the separate native-field tail flag; triangles indicate no passing pulse within the tested range; crosses indicate ineligible preparations. Missing values are not zero response times. Each duration excludes preparation and the subsequent 10,000-unit basal recovery. Spring uses native length and time. See Figure 1 for exact gates and stimulus conventions.'],
                 narrative=['Corrected deBelly is nonmonotone across all three geometries. The positive rank coefficients summarize association rather than its full curve shape.', 'Constant-slope Spring shows the opposite association to its patch responses. Holmes orientation under that cue varies with little monotone association. Censoring and native-field flags restrict interpretation even where rank correlation is strong.'],
                 sources=[dict(label='Current figure legend and gate definitions', path='polarity1d/manuscript/figures/CAPTIONS.md')]),
            dict(title='S2. Numerical qualification and preparation provenance', path='polarity1d/manuscript/THREE_TESTS.md', text='Exact common starting states, duration-boundary audits, corrected deBelly endpoint controls and selected mesh diagnostics; source histories remain distinguishable from compact committed summaries.'),
            dict(title='S3. Mechanism and stochastic sensitivity', path='polarity1d/manuscript/figures/length_mechanisms_20261008/dev/VISION_CHECKPOINT_REQUIREMENTS.md', text='Complete gain/recovery grids, clock qualification, finite-cue definitions, outcome tables, site/vector diagnostics and numerical/noise controls accompany the mechanism checkpoint.'),
            dict(id='supplement4', title='S4. Migration-model adequacy and Hosny sensitivity',
                 role='Supporting conditional sensitivity analysis; it is not an independent biological replication.',
                 assets=[dict(path='wgd_langevin/refined_ou_20261003/figures/supp_04_Hosny_sensitivity.png', caption='Existing Hosny sensitivity supplement from the refined OU analysis.')],
                 legend=['The retained supplement explores signed-drift and leave-Hosny-out sensitivity within the declared OU hierarchy. Its fitted changes and descriptive-score profiles are conditional on the model and supplied summaries, not confidence intervals. Data-only identifiability, Xuan population-curve residuals and prior/pooling assumptions remain part of the adequacy assessment. Representative movie diagnostics are separate: acquisition timing is uncertain and nuclear appearance does not establish ploidy.'],
                 sources=[dict(label='Original Hosny supplement caption', path='wgd_langevin/refined_ou_20261003/figures/supp_04_Hosny_sensitivity_caption.md'), dict(label='Sensitivity methods and results', path='wgd_langevin/refined_ou_20261003/scripts/hosny_sensitivity/METHODS.md')])
        ],
        references=[
            dict(label='Reviewed primary literature and extraction map', path='literature_review/migration_data/studies.csv'),
            dict(label='Scientific scope and author direction', path='MANUSCRIPT_VISION.md'),
            dict(label='Direct-user preference evidence and foundational figure audit', path='reports/manuscript_foundation_audit_20261017.md'),
            dict(label='Mechanism checkpoint intake and live-producer evidence', path='reports/manuscript_mechanism_intake_20261017.md'),
            dict(label='Source checkout, evidence selection and reproduction', path='reports/manuscript_checkout_audit_20261017.md')
        ])


def main():
    import argparse
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--mechanism-sections', type=Path)
    parser.add_argument('--output', type=Path, default=OUT / 'outline_v1.json')
    args = parser.parse_args()
    document = base_document()
    if args.mechanism_sections:
        mechanism = json.loads(args.mechanism_sections.read_text())
        document['sections'][1:1] = mechanism['sections'] if isinstance(mechanism, dict) else mechanism
    args.output.write_text(json.dumps(document, indent=2, ensure_ascii=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
