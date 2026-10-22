"""Package the replicated F4 signatures, full supplements and local HTML.

No integration or new statistical estimator. Readout definitions match Figure3.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import html
import importlib.util
import json
import os
import shutil

import numpy as np
import pandas as pd
from pypdf import PdfReader,PdfWriter

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript'
GROUP=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019'
FIG=GROUP/'figure4'
RUN=MS/'data/runs/revision_noise_strength_replication_scoped_20261019'
OUT=FIG/'signature_assembly_20261019'
HELPER=MS/'code/mechanisms/revision_figure4_html_explorer_20261019.py'
HELPER_SHA='6e88c60938a0fa2308ffe17a3f7728436dca2e65b93175de168cfe72c9ab2d3f'
LABELS=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring_ring='Spring ring')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(p):return str(Path(p).resolve().relative_to(ROOT))
def link(p,out):return Path(os.path.relpath(p,out)).as_posix()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def bundle(target,sections):
    document=PdfWriter()
    for title,path in sections:
        document.append(PdfReader(path),outline_item=title)
    pages=len(document.pages)
    with Path(target).open('wb') as handle:document.write(handle)
    document.close();return pages


def basal_context(out):
    cache=json.loads((RUN/'results/partial/case_cache.json').read_text())
    regions=np.load(out/'supporting/data/postburn_region_samples.npz')
    records=[];cells=[]
    for case_id,case in cache.items():
        meta=case['descriptors'][0]
        if meta['kind']!='basal':continue
        with np.load(case['profile']) as z:
            u=z['primary'][0];C=float(np.broadcast_to(z['circumference'],len(z['time']))[0]);q=z['material_positions'] if 'material_positions' in z else z['x']/C
        M0=float(C*np.mean(u));M1=C*np.mean(u*np.exp(2j*np.pi*q));M2=C*np.mean(u*np.exp(4j*np.pi*q))
        with np.load(case['trace']) as z:
            np.testing.assert_allclose(M0,z['M0'][0],rtol=0,atol=1e-10)
            np.testing.assert_allclose(abs(M1)/M0,z['q1'][0],rtol=0,atol=1e-12)
        records.append(dict(model=meta['model'],size=meta['size'],Cbasal=C,M0_initial=M0,
            M1_initial_real=float(M1.real),M1_initial_imag=float(M1.imag),M1_initial_magnitude=float(abs(M1)),
            q1_initial=float(abs(M1)/M0),q2_initial=float(abs(M2)/M0),contrast_initial=float(np.ptp(u)),
            basal_halfcontrast_largest_physical_width=float(regions[case_id+'__largest_physical'][0]),
            basal_halfcontrast_total_physical_width=float(regions[case_id+'__total_above_fraction'][0]*C),
            profile=case['profile'],profile_sha256=case['profile_sha256']))
        cells.extend(dict(model=meta['model'],size=meta['size'],material_fraction=float(x),primary=float(y),Cbasal=C) for x,y in zip(q,u))
    assert len(records)==12
    pd.DataFrame(records).to_csv(out/'supporting/data/basal_physical_context.csv',index=False)
    pd.DataFrame(cells).to_csv(out/'supporting/data/basal_primary_profiles.csv',index=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=OUT);a=p.parse_args();out=a.out
    generation=json.loads((out/'supporting/metadata/native_generation_snapshot.json').read_text())
    assert generation['complete'] and generation['scientific_cases_available']==936
    data=out/'supporting/data';metadata=out/'supporting/metadata';vectors=out/'supporting/vectors'
    shutil.copy2(data/'condition_window_summaries.csv',data/'condition_window_snapshot.csv')
    assert sha(HELPER)==HELPER_SHA
    spec=importlib.util.spec_from_file_location('unchanged_f4_explorer',HELPER);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    explorer,explorer_proof=helper.render(out)
    basal_context(out)
    main_pages=bundle(out/'figure4_main.pdf',[
        ('A–D · Native local organization',vectors/'figure4_local_dynamics.pdf'),
        ('E–H · Replicated signatures versus circumference',vectors/'figure4_signature_readouts.pdf')])
    sections=[('Temporal memory, usable coverage and full-vector correlation: all windows',vectors/'all_models_sizes_strengths_temporal_windows.pdf'),
        ('Organization, amount/coherence, controls and conditional widths: all windows',vectors/'all_models_sizes_strengths_organization_windows.pdf'),
        ('Every model/size/strength condition: native second-half medoid atlas',vectors/'all_conditions_second_half_native_atlas.pdf'),
        ('Late-quarter observed-window sensitivity',vectors/'late_quarter_sensitivity.pdf')]
    supplement_pages=bundle(out/'figure4_comprehensive_supplement.pdf',sections)
    d=pd.read_csv(data/'condition_window_bootstrap.csv');m=pd.read_csv(data/'memory_window_bootstrap.csv')
    display=[]
    for model,name in LABELS.items():
        size=3. if model=='spring_ring' else 60.
        for f in (.01,.03,.1):
            s=d[(d.model==model)&(d['size']==size)&(d.fraction==f)&(d.kind=='noise')&(d.window=='last_half')].iloc[0]
            t=m[(m.model==model)&(m['size']==size)&(m.fraction==f)&(m.kind=='noise')&(m.window=='last_half')&(m.lag==(800 if model=='spring_ring' else 8000))].iloc[0]
            value=f'{t.direction_cosine_reported:.3f}' if np.isfinite(t.direction_cosine_reported) else 'unavailable'
            display.append(f'<tr><th>{name}</th><td>{f:g}</td><td>{int(s.seeds)}</td><td>{100*s.qualified_multiregion_fraction:.1f}%</td><td>{100*s.low_contrast_fraction:.1f}%</td><td>{value}</td><td>{100*t.direction_coverage:.1f}%</td><td>{s.mean_R:.3f}</td><td>{s.mean_A:.3f}</td><td>{s.mean_Q:.3f}</td></tr>')
    caption=[
        'Figure 4. Circumference effects on activity-heading retention depend on the input regime. Four models retain their own native parameters, physical basals and forcing correlations. Whole positive input mean and fluctuation amplitude scale together. The main summaries use twelve independent seeds at each displayed f=0.001,0.003,0.01,0.03,0.1 condition and discard the first half: chemical times16000–32000 and Spring1600–3200. This is an observed-window choice, not a stationarity claim.',
        'A–D: largest-size native activity histories at f=0.1 in the fixed window T/2–5T/8, with linked profiles at9T/16. Basal, noisy and deterministic ensemble-mean control profiles share native concentration units. Each noisy trajectory is the declared seven-descriptor second-half medoid; no snapshot was selected for extra regions. Region counts use the half-current-contrast threshold and must be read with basal-relative contrast. Goryachev includes a native-scale enlargement. Spring uses material position for its changing circumference.',
        'E1–E4: pooled cosine of the activity-heading difference at native lag8000 for chemicals or800 for Spring, conditional on both endpoint relative first-moment magnitudes R≥0.1. F1–F4: eligible pairs divided by all available lag pairs. Sparse groups with fewer than100 eligible pairs or less than1% coverage remain unavailable. Whiskers are pointwise95% percentile intervals from2000 whole-seed resamples, paired across readouts, windows and lags; no individual time samples are treated as replicates.',
        'G1–G4: fraction of all scalar times with more than one periodic connected arc above the current minimum plus half the current contrast, additionally requiring contrast≥10% of the same case’s basal contrast. H1–H4: fraction below that contrast threshold. These descriptors distinguish organization and loss of contrast; they do not assign competition, site ancestry or exclusive phases. The physical vector P=M₁(t)/M₁(0), active amount M₀ and coherence q₁=|M₁|/M₀ are the same definitions used for Figure3. R=AQ exactly per frame; actual Spring C(t) enters physical quadrature. Activity heading is distinct from cell migration.',
        'The comprehensive supplement retains all nine strengths0.0001–1, all three sizes, every noisy seed and matched mean/no-input control in linked tables and native archives, full/half temporal memory and usable coverage, full-vector correlation, activity amount/coherence, contrast, region-width geometry and quarter-window sensitivity. Five strengths have twelve seeds; the four outer unreplicated levels retain two and are exploratory. Conditional widths state frame/seed weighting and qualification coverage. Original physical-site/angular-remapping and increment evidence remains linked supporting material. Spring’s original selected f=1 run-q90 observation gate remains failed; duration and event-rate inference is outside the accepted primary scope.'
    ]
    (out/'CAPTION.md').write_text('\n\n'.join(caption)+'\n')
    scope_link=link(RUN/'protocol/STABLE_ESTIMAND_SCOPE.json',out)
    failed_link=link(GROUP/'supporting/metadata/FIGURE4_BROAD_NUMERICAL_SPRING_RING_INDEPENDENT_CHECK.json',out)
    acceptance_link=link(RUN/'results/STABLE_ESTIMAND_ACCEPTANCE.json',out)
    baseline_link='supporting/data/basal_physical_context.csv'
    generated=datetime.now(timezone.utc).isoformat()
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure4 · Size, input regime and activity memory</title>
<style>body{{font:17px/1.6 system-ui;max-width:1180px;margin:2rem auto;padding:0 1.2rem;color:#20333e}}h1,h2{{line-height:1.2}}h2{{margin-top:2.5rem}}a{{color:#086c78}}img{{max-width:100%;height:auto;border:1px solid #dce5e9}}.stamp,.caption{{font-size:.9rem;color:#526573}}.notice{{padding:1rem;background:#eef6f7;border-left:4px solid #25818a}}.limit{{padding:1rem;background:#fff6e4;border-left:4px solid #bf8c2d}}nav{{display:flex;flex-wrap:wrap;gap:1rem}}table{{border-collapse:collapse;width:100%;font-size:.9rem}}th,td{{padding:.5rem;border-bottom:1px solid #dce3e7;text-align:right}}th:first-child{{text-align:left}}.scroll{{overflow-x:auto}}{helper.CSS}</style></head><body>
<h1>Figure4: circumference, input regime and activity-heading memory</h1>
<p class="stamp">Report generated {generated}. Native data generated {html.escape(generation['generated_at'])}; last observed {html.escape(generation['data_observed_at'])}: {html.escape(generation['producer'])}. All360 new cases and936 total scientific cases collected; no native failures. Collector stopped after terminal production. This static report has no timed refresh.</p>
<nav><a href="#main">Main evidence</a><a href="#results">Measured context</a><a href="#explorer">Full experiment explorer</a><a href="#readouts">Definitions</a><a href="#supplement">Complete supplements</a><a href="#limits">Scope and limits</a></nav>
<p class="notice">Two complementary signatures are retained: circumference-dependent heading memory across input regimes, and the local organization accompanying memory loss. Larger rings can retain more heading memory at intermediate input while also supporting more threshold regions; stronger input can remove that retention. A weak global vector does not identify activity depletion, uniformity or a unique local pattern.</p>
<section id="main"><h2>Short native histories and replicated size summaries</h2><p>All main curves use the second half after discarding the first half, without assuming a stationary distribution. The temporal and quarter-window supplements retain the observed sensitivity to that choice. The two main sheets keep twelve independent seeds per displayed condition and whole-seed uncertainty.</p>
<p><a href="figure4_main.pdf">Two-sheet main PDF</a> · <a href="CAPTION.md">Legend</a> · <a href="supporting/data/main_signature_points.csv">Exact points/intervals</a></p>
<figure><a href="figure4_local_dynamics.png"><img src="figure4_local_dynamics.png" alt="Four short native histories and exact activity profiles with matched mean controls"></a><figcaption class="caption">Actual unshifted native amplitudes, fixed time selection, declared medoids. Spatial patterns and weak resultants require separate interpretation.</figcaption></figure>
<figure><a href="figure4_signature_readouts.png"><img loading="lazy" src="figure4_signature_readouts.png" alt="Directional memory, eligible pair coverage and contrast-qualified organization versus circumference"></a><figcaption class="caption">Activity-heading memory and usable coverage are paired. Counts and low contrast remain separate. Native clocks and forcing correlations differ; these input factors are not calibrated cross-model doses.</figcaption></figure></section>
<section id="results"><h2>Largest-size measured context</h2><p>These rows use the same second-half windows and native lags as the main curves. Read direction jointly with coverage, and R jointly with active amount A and coherence Q. Their means are not an additive or multiplicative causal partition: mean R is the mean of paired A×Q.</p><div class="scroll"><table><thead><tr><th>Model</th><th>f</th><th>Seeds</th><th>Qualified multi</th><th>Low contrast</th><th>Direction cosine</th><th>Pair coverage</th><th>Mean R</th><th>Mean A</th><th>Mean Q</th></tr></thead><tbody>{''.join(display)}</tbody></table></div></section>
{explorer}
<section id="readouts"><h2>Same physical activity vector as Figure3</h2>{helper.MATH}
<p>The linked <a href="{baseline_link}">twelve-basal context table</a> retains physical M₀, |M₁|, raw q₁/q₂, contrast and half-contrast arc widths, so normalization does not hide initial size differences. M₀ is active-only amount, distinct from the native conserved total pool. Spring scalar quantities are saved at cadence0.5, chemical quantities at5; native profiles are saved at4 and40 respectively.</p>
<p>A threshold region is a circular connected set of cells above current min+half current contrast. Spatial widths join seam-crossing regions and use actual C(t). Largest-region means condition on amplitude-qualified frames; second-largest means condition on qualified multiregion frames. Condition means average available seed descriptors equally. Qualification fractions, exact frame counts and contributing seeds accompany these widths. No region birth or lifetime classifier is implied.</p></section>
<section id="supplement"><h2>Full underlying experiment and supporting evidence</h2><p><a href="figure4_comprehensive_supplement.pdf">Complete {supplement_pages}-page current-experiment supplement</a> contains every model, size and strength, all full/half temporal curves and coverage, organization/amount/coherence/contrast/mean controls, instantaneous widths, a108-condition native atlas and late-quarter sensitivity. Every individual seed remains in the tables; native arrays remain linked by exact archive hashes.</p>
<ul><li><a href="supporting/data/postburn_native_source_index.csv">All936 native case identities/archive hashes</a> and <a href="supporting/data/worker_inventory.csv">all360 new producer records</a></li>
<li><a href="supporting/data/case_window_descriptors.csv">All full/half case descriptors</a>, <a href="supporting/data/case_window_memory.csv">case temporal sufficient statistics</a> and <a href="supporting/data/condition_window_bootstrap.csv">condition uncertainty</a></li>
<li><a href="supporting/data/postburn_spatial_case_windows.csv">Spatial width/count/coverage case tables</a> and <a href="supporting/data/late_quarter_sensitivity.csv">paired late-quarter differences</a></li>
<li><a href="supporting/data/all_condition_native_atlas_selection.csv">Every medoid/source selection</a>, <a href="supporting/data/local_medoid_scores.csv">all seed selection scores</a> and <a href="supporting/data/basal_primary_profiles.csv">native basal profiles</a></li>
<li><a href="{link(FIG/'figure4_with_spring_geometry_patterns.png',out)}">Retained physical-site/common-geometry evidence at the original two strengths</a>; <a href="{link(FIG/'supporting/metadata/content.json',out)}">original complete supporting inventory</a></li>
<li><a href="{link(FIG/'broad_strength_campaign_20261019/review_complete_scientific/report_accepted_20261019T201733.html',out)}">Preserved accepted576-case interactive snapshot</a> and <a href="{link(FIG/'postburn_signature_candidate_20261019/README.md',out)}">preserved first signature draft</a></li></ul></section>
<section id="limits"><h2>What the evidence establishes and leaves open</h2>
<p class="limit">The original Spring C₀=3,f=1 run-q90 observation check remains <strong>failed</strong>:12 coarse/dense duration comparisons failed, while ODE and the selected occupancy/amount/coherence/memory measurements passed. No source, forcing or solver was retuned. A <a href="{scope_link}">separate ten-metric scope</a> was <a href="{acceptance_link}">explicitly accepted</a>; both raw/qualified run-q90 descriptors are excluded from that qualified scope. The <a href="{failed_link}">original passed:false receipt</a> remains unchanged. All run spans are cadence-dependent sampled descriptors, not continuous lifetimes or event rates. Instantaneous widths remain descriptive geometry, outside that numerical acceptance.</p>
<p>Controls must remain visible: Goryachev C20,f=0.01 already loses its first moment under the smooth matched ensemble-mean input; larger high-input Goryachev cases also often become low contrast under the mean control. This is a forced-state response, not zero-mean wandering about an unchanged incumbent. Other models can retain contrast while losing heading memory and supporting multiple threshold regions.</p>
<p>Second-half and late-quarter summaries remain finite observations, not equilibrium certification or discrete phase boundaries. Whole-seed intervals condition on exact common basals, native parameters, forcing family and fixed windows. Four selected upper-input numerical screens do not prove convergence of every possible condition. Outer strengths remain exploratory two-seed outcomes. Multiple regions alone do not establish competition or lineage.</p>
<p>The biological bridge is a measurable analogy: activity-heading persistence can be imaged alongside free migration, whereas a quantitative migration autocorrelation requires a measured activity-to-motion coupling. Circumference changes angular geometry as well as organization; the present heading summaries do not isolate physical site displacement or imply a velocity law.</p></section></body></html>'''
    (out/'index.html').write_text(page)
    narrative=['Circumference-dependent activity-heading retention changes with the declared input regime; global heading memory and its usable coverage must be read together.',
               'Memory loss can accompany contrast-rich multiple-region organization or a low-contrast forced state. Native profiles, active amount and coherence separate these observable outcomes without assigning competition or a migration law.']
    assets=[]
    for stem,label in [('figure4_local_dynamics','Short second-half native histories and linked actual profiles with matched mean controls.'),('figure4_signature_readouts','Replicated heading memory, usable coverage and spatial organization versus circumference.')]:
        assets.append(dict(kind='main',path=str(out/(stem+'.png')),pdf=str(vectors/(stem+'.pdf')),svg=str(vectors/(stem+'.svg')),caption=label))
    content=dict(title='Figure 4 · Input regime changes the relationship between circumference and activity-heading memory',role='stochastic activity signatures',
        claim=narrative[0],narrative=narrative,assets=assets,main_asset=assets[0]['path'],legend=caption,
        measurements=['Conditional heading memory with usable-pair coverage','Contrast-qualified threshold-region occupancy','Low-contrast occupancy','Native field organization, active amount and coherence','Instantaneous periodic region widths and observed-window sensitivity'],
        limitations=['Native clocks/input references are not cross-model calibrations','No stationary-distribution, exclusive-phase, site-lineage or competition classifier','Activity heading is distinct from cell migration','Spring sampled run-duration inference is excluded; width geometry remains descriptive'],
        sources=[dict(label='Frozen broad-strength campaign',path=rel(MS/'data/runs/revision_noise_strength_broad_20261019/manifest.json')),
                 dict(label='Frozen replicated native cases',path=rel(RUN/'manifest.json')),
                 dict(label='Explicit stable measurement scope',path=rel(RUN/'protocol/STABLE_ESTIMAND_SCOPE.json')),
                 dict(label='Original failed Spring observation receipt',path=rel(GROUP/'supporting/metadata/FIGURE4_BROAD_NUMERICAL_SPRING_RING_INDEPENDENT_CHECK.json')),
                 dict(label='All native sources and frame coverage',path=rel(data/'postburn_native_source_index.csv')),
                 dict(label='Joint whole-seed uncertainty protocol',path=rel(metadata/'BOOTSTRAP_PROTOCOL.json'))],
        supplementary=[dict(title='Full current underlying experiment',path=str(out/'figure4_comprehensive_supplement.pdf'),caption=f'All nine strengths, sizes, controls, temporal/window and organization suites; {supplement_pages} pages.'),
                       dict(title='Interactive full experiment',path=str(out/'index.html'),caption='All model/size/strength/window readouts with actual seed support and archived source links.'),
                       dict(title='Retained physical-site and angular-remapping evidence',path=str(FIG/'figure4_with_spring_geometry_patterns.png'),caption='Original two-strength geometry diagnostics; folded endpoint displacements are not travel.')])
    write(metadata/'content.json',content)
    readme=f'''# Figure4: replicated activity signatures

[Main two-sheet PDF](figure4_main.pdf) · [Scientific HTML report](index.html) · [Legend](CAPTION.md) · [Full {supplement_pages}-page supplement](figure4_comprehensive_supplement.pdf)

This package uses all936 scientific cases, with all360 new native replicas collected at {generation['generated_at']}. All producers are terminal. Main readouts use the second-half windows and12 independent seeds at each of five common strengths. Other four strengths retain2 seeds and remain exploratory. The accepted earlier pilot and signature drafts remain preserved.

The two main signatures are activity-heading memory across circumference and the contrasting local organizations accompanying memory loss. The same physical M₁/P and active-only M₀ definitions as Figure3 are used. Native model units/noise references differ; no cell-velocity law is assumed.

The comprehensive PDF covers all9strengths×3sizes×4models, full/half direction/full-vector/coverage curves, organization/contrast/amount/coherence with matched mean controls, conditional physical region widths, all108condition-native medoids and late-quarter sensitivity. Every seed and control is retained in the linked data tables and exact archive index. Region widths condition on qualified frames and then average available seed summaries; contributors and profile coverage are exported. They are descriptive geometry, not lifetimes.

The original failed Spring run-q90 receipt remains intact. The accepted primary scope explicitly excludes duration inference, and the README/report do not relabel the whole original gate as passed. Sampling windows do not certify stationarity. Confidence intervals resample whole independent seeds, jointly across readouts/windows/lags, conditional on the exact basals and forcing law.

Main selection, uncertainty and sources are detailed in [content.json](supporting/metadata/content.json), [bootstrap protocol](supporting/metadata/BOOTSTRAP_PROTOCOL.json), [every native archive](supporting/data/postburn_native_source_index.csv), [baseline physical context](supporting/data/basal_physical_context.csv), [late-quarter audit](supporting/metadata/LATE_QUARTER_AUDIT.json) and [medoid scores](supporting/data/local_medoid_scores.csv).
'''
    (out/'README.md').write_text(readme)
    receipt=dict(generated_at=generated,source=rel(__file__),source_sha256=sha(__file__),native_integrations=0,scientific_cases=936,
        main_pdf_pages=main_pages,supplement_pdf_pages=supplement_pages,html_sha256=sha(out/'index.html'),explorer=explorer_proof,
        independent_final_review='pending',original_spring_whole_gate_still_failed=True)
    write(metadata/'FINAL_REPORT_RECEIPT.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
