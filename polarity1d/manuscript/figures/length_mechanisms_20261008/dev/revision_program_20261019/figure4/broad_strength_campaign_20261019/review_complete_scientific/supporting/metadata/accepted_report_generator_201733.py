"""Build the first browsable Figure 4 progress report from reviewed exports."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import html
import importlib.util
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GROUP = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019'
OUT = GROUP / 'figure4/broad_strength_campaign_20261019/review_complete_scientific'
LABELS = {'wave_pinning': 'Wave-pinning', 'goryachev': 'Goryachev', 'holmes_model3': 'Holmes', 'spring_ring': 'Spring'}


def main():
    helper_path = ROOT / 'polarity1d/manuscript/code/mechanisms/revision_figure4_html_explorer_20261019.py'
    spec = importlib.util.spec_from_file_location('figure4_html_explorer', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    explorer, explorer_receipt = helper.render(OUT)
    limits = json.loads((OUT / 'supporting/metadata/SPRING_SAMPLING_LIMITS.json').read_text())
    failure_rows = ''.join(f"<tr><th>{r['window'].replace('_',' ')}</th><td>{r['reference']:g}</td><td>{r['comparison']:g}</td><td>{r['absolute_difference']:g}</td><td>{r['tolerance']:g}</td></tr>"
                           for r in limits['failures'] if r['check']=='observation_standard' and r['metric']=='multiregion_run_q90')
    source = OUT / 'supporting/data/largest_size_memory_context.csv'
    data = pd.read_csv(source)
    selected = data[data.fraction.isin([.01, .03, .1])]
    rows = []
    for model in LABELS:
        for r in selected[selected.model == model].sort_values('fraction').to_dict('records'):
            cosine = f"{r['direction_cosine_reported']:.3f}" if pd.notna(r['direction_cosine_reported']) else 'unavailable'
            rows.append(f"<tr><th>{LABELS[model]}</th><td>{r['fraction']:g}</td><td>{int(r['available_cases'])}</td>"
                        f"<td>{100*r['qualified_multiregion_fraction']:.1f}%</td><td>{100*r['low_contrast_fraction']:.1f}%</td>"
                        f"<td>{cosine}</td><td>{100*r['direction_coverage']:.1f}%</td></tr>")
    generated = datetime.now(timezone.utc).isoformat()
    page = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Figure 4 · broad noise experiment and directional memory</title>
<style>
body{font:17px/1.6 system-ui,sans-serif;color:#20333e;background:#fff;max-width:1180px;margin:2rem auto;padding:0 1.2rem}
h1,h2,h3{line-height:1.22}h1{max-width:950px}h2{margin-top:2.7rem}a{color:#086c78}img{max-width:100%;height:auto;border:1px solid #dce5e9}
.stamp,.caption{font-size:.9rem;color:#526573}.notice{background:#eef6f7;padding:1rem 1.2rem;border-left:4px solid #25818a}
.qualification{background:#fff7e8;padding:1rem 1.2rem;border-left:4px solid #be8a28}.links{display:flex;flex-wrap:wrap;gap:.65rem 1.2rem}
.cards{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1.2rem}.card{background:#f6f8f9;padding:1rem 1.2rem}
table{border-collapse:collapse;width:100%;font-size:.93rem}td,th{padding:.5rem .65rem;border-bottom:1px solid #d9e3e7;text-align:right}th:first-child{text-align:left}
.scroll{overflow-x:auto}figure{margin:1.5rem 0}details{margin:1.2rem 0}summary{cursor:pointer;font-weight:650}nav{position:sticky;top:0;background:#ffffffed;padding:.7rem 0;border-bottom:1px solid #dce5e9;z-index:1}
@media(max-width:760px){.cards{grid-template-columns:1fr}body{font-size:16px}}
__EXPLORER_CSS__
</style></head><body>
<h1>Figure 4: what changes as the whole noise signal becomes stronger?</h1>
<p class="stamp">Report generated __GENERATED__. Scientific plots/data snapshot: 19 October 2026, 19:36 UTC.
Native production last observed terminal at 20:02:32 UTC: 256/256 new outputs, 576 scientific cases, no native worker failures.
This page updates only when regenerated and manually refreshed. No timed refresh.</p>
<nav class="links"><a href="#findings">Findings</a><a href="#explorer">Interactive comparison</a><a href="#overview">All strengths and sizes</a><a href="#profiles">Native profiles</a><a href="#memory">Directional memory</a><a href="#readouts">Readout guide</a><a href="#qualification">Sampling qualification</a><a href="#next">Next steps</a></nav>
<p class="notice"><strong>The broader experiment is complete.</strong> Four models, three sizes each, nine common strength factors.
Wave-pinning, Holmes and Spring reach frequent multiple-region organization at stronger inputs.
Large Goryachev cells often instead become low-contrast relative to their initial incumbent.
Local organization and directional memory provide complementary information.</p>
<p class="qualification"><strong>The original Spring numerical gate did not pass in full.</strong> In the prespecified largest-size, strongest-input path, 12 comparisons of sampled multiple-region run durations fail when observations are saved ten times more densely. All66 canonical ODE comparisons and120 other observation comparisons pass in that screen. These failures do not invalidate the independently checked native data; they limit duration interpretation. The failed receipt is preserved. <a href="#qualification">Exact differences and scope</a>. Replication remains prepared and unsubmitted.</p>
<p class="links"><a href="broad_strength_review.pdf">Two-page PDF</a><a href="README.md">Detailed scientific account</a><a href="index.html">All sizes: raw and contrast-qualified tables</a><a href="supporting/data/condition_window_snapshot.csv">Exact measurements</a></p>

<section id="findings"><h2>What the new range adds</h2>
<div class="cards"><article class="card"><h3>Multiple regions also occur in WP, Holmes and Spring</h3>
<p>The earlier two strengths did not span comparable changes in all models. At larger tested inputs, these three models spend substantial time with multiple activity regions and retain strong spatial contrast. Size changes how early and how frequently this happens.</p></article>
<article class="card"><h3>A high region count can accompany a weak pattern</h3>
<p>In large Goryachev cells, stronger input often produces small spatial variations relative to the original tall incumbent. Counting regions relative to each frame's own small contrast can therefore suggest many peaks even when the field has little basal-relative contrast.</p></article></div>
<p>Small Goryachev cells also show substantial contrast-qualified multiple-region activity. The data support size-dependent differences within the sampled input family, rather than a model-wide binary classification.</p>
<p>Every model keeps its original positive-input noise law and correlations. The factors are 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3 and 1. Scaling the whole signal increases its mean and fluctuation amplitude together. Model-native forcing and clocks are not calibrated across models.</p>
</section>

__EXPLORER__

<section id="overview"><h2>The complete strength sweep</h2>
<p>Read down the model rows and across the three readouts. Colours identify sizes within each row. Left: time with multiple regions, including the basal-contrast qualification. Middle: global activity-vector magnitude. Right: active amount and spatial coherence, which together determine that magnitude.</p>
<figure><a href="broad_strength_exploratory.png"><img loading="lazy" src="broad_strength_exploratory.png" alt="All four models, three sizes and nine noise factors, with region occupancy, vector magnitude, activity amount and coherence"></a>
<figcaption class="caption">Click for the full-resolution image. New levels have two seeds; the original 0.001 and 0.01 levels have twelve. Whiskers show observed seed ranges. Connecting segments guide the eye across sampled factors. The image header retains its19:36 snapshot (three numerical companions then complete); all scientific curves were already complete and unchanged. The final fourth check and its duration limitation are reported above.</figcaption></figure>
</section>

<section id="profiles"><h2>What those changes look like locally</h2>
<p>These profiles use the largest original size of each model, common factors 0.01, 0.03 and 0.1, the first saved seed and a fixed midpoint time. They preserve native concentrations and remain unshifted. They illustrate actual states; the preceding time fractions measure prevalence across the full trajectories.</p>
<figure><a href="broad_strength_transition_profiles.png"><img loading="lazy" src="broad_strength_transition_profiles.png" alt="Native activity profiles for four models across three common noise strengths"></a>
<figcaption class="caption">Goryachev's inset enlarges its small high-input variations while retaining native units. Spring uses material position because its physical circumference changes.</figcaption></figure>
<details><summary>A complementary two-region example</summary><p>At the smallest chemical size and factor 1, both selected WP and Goryachev profiles contain two regions, but their global resultants differ markedly. The WP snapshot was selected conditionally to demonstrate this possibility; it is not a typical-state estimate.</p><a href="broad_strength_native_profiles.png"><img loading="lazy" src="broad_strength_native_profiles.png" alt="Two-region WP and Goryachev profiles with different global activity-vector magnitudes"></a></details>
</section>

<section id="memory"><h2>Organization changes and directional memory can separate</h2>
<p>In WP and Holmes, increasing multiple-region occupancy accompanies substantial long-lag memory loss in the largest cells. Spring already loses appreciable memory before multiple regions become frequent: at factor 0.003 its largest size has only 0.66% qualified multiple-region occupancy, magnitude near its basal value, and lag-800 direction cosine 0.543.</p>
<p>Large Goryachev cells require a different reading. As the resultant weakens, fewer direction pairs remain usable. A low conditional correlation with poor coverage should be read alongside that loss of directional definition.</p>
<p><strong>Largest size per model, full observation windows:</strong> direction cosines use lag 8000 in chemical native units and 800 in Spring native units. These clocks are not cross-model calibrations.</p>
<div class="scroll"><table><thead><tr><th>Model</th><th>Factor</th><th>Seeds</th><th>Qualified multiple regions</th><th>Low contrast</th><th>Direction cosine</th><th>Usable direction pairs</th></tr></thead><tbody>__ROWS__</tbody></table></div>
<p class="caption">Both direction endpoints require R ≥ 0.1. Pair coverage is relative to all available lag pairs. Means and memory estimates are finite-window observations; new-level seed counts remain small.</p>
</section>

<section id="readouts"><h2>A short guide to the measurements</h2>
__MATH__
<p><strong>Multiple regions:</strong> more than one connected arc above the current minimum plus half the current contrast. The qualified fraction additionally requires contrast at least 10% of that case's initial basal contrast. Both fractions use every saved time in their denominator.</p>
<p><strong>Activity vector:</strong> sum local activity with the ring's outward-normal directions. Its direction measures the global activity heading. Its magnitude R is relative to its own initial value and can exceed one; it is not cell speed.</p>
<p><strong>Amount and coherence:</strong> A is total active amount relative to basal; Q is first-harmonic coherence relative to basal. R = AQ at each frame. Their separate changes explain why high R need not imply a single unchanged peak, and low R need not mean little total activity.</p>
<p><strong>Memory:</strong> the mean cosine between usable global directions separated by a declared lag. The report retains eligible-pair coverage and full/first-half/last-half estimates, since organization can evolve during the finite experiment.</p>
</section>

<section id="qualification"><h2>Exactly what the Spring sampling check limits</h2>
<p>The selected Spring C₀=3, f=1 path was solved with the same finite forcing and seed using standard and tighter ODE tolerances. Canonical observations at cadence0.5 are exact subsets of dense observations at cadence0.05. The raw and contrast-qualified multiregion-run90th percentiles give the same values here:</p>
<div class="scroll"><table><thead><tr><th>Window</th><th>Cadence0.5</th><th>Cadence0.05</th><th>Difference</th><th>Frozen tolerance</th></tr></thead><tbody>__FAILURE_ROWS__</tbody></table></div>
<p>These three differences recur for raw and qualified spans in both standard and tight solves:12 failed comparisons. All66 canonical ODE checks pass; the remaining120 observation checks for occupancy, weak-resultant/low-contrast fractions, amount/coherence/magnitude and direction/coverage also pass. This is observation-cadence sensitivity, rather than evidence of a changed forcing realization or an ODE error in those tested estimands.</p>
<p>The result directly establishes a limit for that tested high-input path. It does not demonstrate failure at lower strengths. Every strength's spans remain descriptors at their actual sampling cadence, with edge censoring, rather than continuous site lifetimes or event rates. Quantitative cadence-invariant span comparisons at the failed scope are withheld. Any acceptance limited to stable estimands must be recorded explicitly; the original passed:false receipt will not be relabelled.</p>
<p>Mean controls are separate native outcomes. At Goryachev C20,f=0.01 the smooth ensemble-mean control already removes the first moment; that condition must be read as a forced-regime response, not clean spontaneous wandering about an unchanged incumbent.</p>
<p class="links"><a href="supporting/metadata/SPRING_SAMPLING_LIMITS.json">Exact scope and terminal status</a><a href="supporting/data/spring_sampling_failed_checks.csv">All12 failed rows</a><a href="../../../supporting/metadata/FIGURE4_BROAD_NUMERICAL_SPRING_RING_INDEPENDENT_CHECK.json">Preserved independent failed receipt</a><a href="../../../supporting/metadata/FIGURE4_BROAD_NEW_NATIVE_PROGRESS_INDEPENDENT_CHECK.json">Independent new-native data audit</a></p></section>

<section id="next"><h2>What remains open</h2>
<p>The current candidate story combines the existing dynamics panels with the broader organization sweep and native examples. The main-figure composition remains open; these pages are a review workspace.</p>
<p>A concrete replication proposal adds ten seeds at each of three common factors, 0.003, 0.03 and 0.1, across every model and size: 360 additional trajectories. This would estimate variability around early memory changes and the observed organization changes. It has not been submitted; the scope of Spring duration-sensitive readouts is being resolved first.</p>
<p>The two-seed pilot does not establish precise transition boundaries or exclude another organization outside the tested range. Multiple regions alone do not establish functional competition or identify site lineage.</p>
<p>Useful experimental readouts would measure local contrast and multiple-region occupancy alongside total active signal and the global activity heading. Sampling cadence should be checked before interpreting persistence of identified regions. These are proposals for measurable observables, not a classifier for molecular mechanisms or a cell-velocity law.</p>
<p class="links"><a href="supporting/data/memory_all_conditions_windows_lags.csv">All memory and coverage estimates</a><a href="supporting/data/organization_memory_seed_window_sensitivity.csv">Seed and window sensitivity</a><a href="supporting/metadata/REPLICATION_PROPOSAL.json">Replication proposal</a><a href="../review/README.md">Preserved earlier partial review</a><a href="report_initial_20261019.html">Preserved initial HTML snapshot</a></p>
</section></body></html>'''
    page = page.replace('__GENERATED__', html.escape(generated)).replace('__ROWS__', ''.join(rows))
    page = page.replace('__EXPLORER_CSS__',helper.CSS).replace('__EXPLORER__',explorer).replace('__MATH__',helper.MATH).replace('__FAILURE_ROWS__',failure_rows)
    path = OUT / 'report.html'
    temporary = path.with_suffix('.html.tmp')
    temporary.write_text(page)
    temporary.replace(path)
    receipt = {'generated_at': generated, 'html_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'source_csv_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'displayed_rows': len(rows),
               'native_integrations': 0, 'timed_browser_refresh': False,
               'scope': 'Expanded author-requested local HTML report; native/art sources retained, interactive controls use audited aggregate exports, Spring duration gate remains failed.',
               'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'explorer':explorer_receipt,
               'initial_html_sha256':hashlib.sha256((OUT/'report_initial_20261019.html').read_bytes()).hexdigest(),
               'spring_limits_sha256':hashlib.sha256((OUT/'supporting/metadata/SPRING_SAMPLING_LIMITS.json').read_bytes()).hexdigest()}
    (OUT / 'supporting/metadata/HTML_PROGRESS_DELIVERY.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
