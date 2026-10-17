"""Build the concise replacement reader from audited current figures."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import html
import json
import shutil
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'figures/length_mechanisms_20261008'
DATA = ROOT / 'data/runs/length_mechanisms_20261008'
TITLE = 'Mechanisms of cell-length-dependent polarity reversal'

CHAPTERS = [
    dict(id='responses', title='Length dependence of reversal cue duration', image='length_responses.png', shape=(11.7, 15.5),
         paragraphs=[
             'The minimum reversal cue duration depends on both cell length and stimulus geometry. Fixed-width patches produce pronounced increases in Wave-pinning and Holmes, whereas Goryachev slows under all three stimuli. Constant-slope gradients show weaker length dependence in several models, while their peak increases with cell length.',
             'Cue duration is the shortest tested pulse that produces reversed polarity after withdrawal and basal recovery. The orientation criterion requires persistent reversal of direction; the strength criterion additionally requires recovery of the original polarity magnitude within 1%. Redistribution may continue after withdrawal. The mechanistic analyses below relate these responses to initial pole material, chemical reaction kinetics and cortical-flow memory.'],
         caption='Figure 1. Minimum reversal cue duration in five models under fixed-width, relative-width and constant-slope stimuli. All tests share the same basal-polarized starting state at each length. Patch peaks are fixed; support is 3 (Spring 0.2) or 20% of length. Blue denotes the orientation criterion and orange the strength criterion. Open markers indicate failed native-field stability flags, triangles indicate tested duration limits, and crosses indicate ineligible preparations. Axes use model-native units.',
         sources=[('../supplement_width_response.pdf', 'Full vector duration figure'), ('../../data/runs/debelly_corrected_three_test_20261008/results/index.html', 'Current comparison and audit')]),
    dict(id='material', title='Larger active poles increase the redistribution requirement', image='initial_polarity_length.png', shape=(11.7, 12.0),
         paragraphs=[
             'Across all five models, longer cells contain more active material in the pre-existing pole. Reversal requires sufficient redistribution to establish a new pole that maintains polarity after countercue withdrawal. A length-independent reversal time therefore requires redistribution capacity to increase with the initial material requirement.',
             'The additional material is accommodated by different pole geometries. Goryachev develops a taller peak of nearly fixed width. Wave-pinning, deBelly, Holmes and Spring predominantly increase pole width. These differences motivate separate analyses of inactive availability and restoring kinetics in Goryachev, and persistent domain expansion in Wave-pinning.'],
         caption='Figure 2. Initial basal-polarized states at all 97 eligible model–length combinations. Integrated active material (A), half-contrast pole width (B) and pole concentration (C) distinguish expansion in width from expansion in height. Opposite-end inactive concentration is shown in D. Solid amount includes background; dashed amount subtracts the field minimum. Pole regions are connected to the original boundary and exceed half the active contrast. Amount is a one-dimensional concentration integral.',
         sources=[('../initial_polarity_length/initial_polarity_length_manuscript.pdf', 'Vector starting-state figure'), ('../../data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv', 'Corrected starting-state measurements')]),
    dict(id='initiation', title='Inactive depletion reduces countercue activation in Goryachev', image='comparative_revision/equal_cue_initiation.png', shape=(11.7, 15.5),
         paragraphs=[
             'The taller Goryachev pole is accompanied by lower inactive concentration throughout the cell. Its unsaturated autocatalytic activation remains strong at high active concentration despite inactive depletion: the zero-reaction contour falls toward lower inactive concentration as active concentration increases. Wave-pinning’s activation feedback saturates, while deactivation continues to increase with active concentration. This kinetic contrast helps explain the taller Goryachev peak and the nearly fixed Wave-pinning plateau, where additional material is accommodated through domain widening.',
             'In both models, the added cue contribution is cue strength multiplied by inactive concentration. At the opposite end, an equal local cue therefore produces only about one third as much initial activation in Goryachev at L = 39 as at L = 9; the Wave-pinning response is nearly unchanged. Increased cue support or the increasing peak of a constant-slope gradient can partly compensate for this initiation penalty.'],
         caption='Figure 3. Active (A–B) and inactive (C–D) profiles at L = 9, 21 and 39, and basal net-reaction maps (E–F). Red denotes net activation and blue net deactivation as functions of active concentration u and inactive concentration v; diffusion is excluded. Solid curves trace spatial profile states, circles mark the opposite end, triangles mark the initial pole, and dashed curves indicate zero net reaction. Insets show added cue activation s × v, normalized to L = 9. Scales are model-native.',
         sources=[('comparative_revision/equal_cue_initiation.pdf', 'Individual figure PDF'), ('../../data/runs/length_mechanisms_20261008/comparative_revision/provenance.json', 'Reaction maps and measurements')]),
    dict(id='restoration', title='Goryachev restoring responses increase with cell length', image='comparative_revision/incumbent_restoration.png', shape=(11.7, 15.5),
         paragraphs=[
             'Resistance of the pre-existing pole was assessed by removing 0.5% of the entire inactive profile for 0.1 native time units and returning the removed material pointwise. Reactions and diffusion were unchanged. Active loss during the perturbation and the restoring reaction after return distinguish susceptibility to depletion from the subsequent recovery response.',
             'Between L = 9 and 39, Goryachev’s fractional active loss decreases by approximately half, while its immediate restoring reaction per unit damage increases 15.7-fold. Wave-pinning’s fractional loss remains nearly constant and its corresponding restoring response increases only 1.3-fold. The taller Goryachev pole is thus less susceptible to this perturbation and develops a stronger immediate restoring response relative to the material lost.'],
         caption='Figure 4. Active-core loss during depletion and recovery (A–B), restoring reaction after inactive return (C–D), and fractional damage and immediate restoring response across 18 lengths (E–F). The shaded interval denotes 0.5% inactive removal for 0.1 native time units; exact pointwise return occurs at time zero. Measurements use the fixed original half-height core and subtract a matched sham. Reaction is normalized by the active deficit at return and excludes diffusion. Scales are model-native.',
         sources=[('comparative_revision/incumbent_restoration.pdf', 'Individual figure PDF'), ('../Goryachev_wave_pinning_trial/inactive_pool_perturbation/index.html', 'Full perturbation comparison'), ('../../data/runs/inactive_pool_figure_exports_20261008/README.md', 'Perturbation methods and audit')]),
    dict(id='sink', title='Pole separation and local resistance limit material loss', image='sink_distance_and_resistance.png', shape=(11.7, 10.2),
         paragraphs=[
             'Inactive material was continuously removed by a width-1 sink located either at the opposite cell end or at a fixed gap of 3 beyond the initial pole edge. Opposite-end separation increases with cell length. The fixed-gap sink accelerates pole erosion, with the largest effect in long cells, consistent with a limitation imposed by communication distance through the inactive pool.',
             'A length effect remains when separation is held fixed. Over L = 9–39, fractional active-material loss slows 4.9-fold in Goryachev and 3.0-fold in Wave-pinning. In conjunction with the depletion-and-return experiment, this supports a contribution from Goryachev’s length-dependent pole resistance. Wave-pinning can lose active material through edge retreat while largely preserving plateau height.'],
         caption='Figure 5. Peak-height (A–B) and active-amount (C–D) loss rates in the saturated strong-sink regime, k = 1000. Both sinks have width 1 and remain fixed throughout each trajectory. Rates average 10–30% loss relative to the initial peak or core amount. The fixed-gap active-amount loss rate decreases 3.0-fold in Wave-pinning and 4.9-fold in Goryachev over L = 9–39. Separate axes retain model-native scales. No reversal cue is applied.',
         sources=[('sink_distance_and_resistance.pdf', 'Vector plateau figure'), ('../inactive_sink_loss/index.html', 'Full strength sweep, morphology and trajectories'), ('../../data/runs/length_mechanisms_20261008/sink/plateau_measurements.csv', 'Selected plateau measurements')]),
    dict(id='interface', title='Wave-pinning retains active material through domain expansion', image='support/wave_pinning_interface.png', shape=(11.7, 10.5),
         paragraphs=[
             'Persistent Wave-pinning reversal involves expansion of the high-active domain. Cue-induced increases in new-pole height relax after withdrawal, whereas expansion through interface advance retains active material at the new end. The original domain narrows as its interface retreats.',
             'Longer cells require a wider new domain to retain a comparable share of active material. A fixed-width cue does not increase its support with this requirement; a relative-width cue improves coverage, but substantial interface displacement is still required. Unlike Goryachev’s taller peak, the wider Wave-pinning pole can lose material through edge retreat without uniformly reducing plateau concentration.'],
         caption='Figure 6. Spatial profiles during the relative-20% cue (A) and after withdrawal (B), pole widths during stimulation (C), and post-withdrawal polarity for adjacent effective and ineffective doses (D), at L = 39. After the 615-unit pulse, the new peak relaxes from approximately 1.56 to 1.12 while the high-active domain widens and recovers initial polarity strength. The 202- and 204-unit pulses straddle the orientation criterion. Width uses a fixed active-concentration threshold of 0.6.',
         sources=[('support/wave_pinning_interface.pdf', 'Vector interface figure'), ('../../data/runs/length_mechanisms_20261008/comparative_support/wave_pinning_interface_measurements.csv', 'Native trajectory measurements')]),
    dict(id='growth', title='Interface advance increases less than the material requirement', image='support/wave_pinning_length_growth.png', shape=(11.7, 8.8),
         paragraphs=[
             'Under equal-peak relative-20% cues, interface advance increases with length less strongly than the initial pole width. Between L = 21 and 39, the original high-active domain is 1.83 times wider, whereas average new-interface advance over cue times 25–60 is only 1.33 times faster.',
             'Normalized to the initial pole width, domain expansion is approximately 28% slower in the longer cell. Greater interface displacement is therefore required without proportionately faster expansion. The post-withdrawal trajectories establish why transient increases in concentration cannot substitute for this retained domain growth.'],
         caption='Figure 7. New-pole profiles at common cue times (A) and domain width normalized to the initial pole width (B), at L = 21 and 39. Cue peak is identical and support is 20% of length. Profiles are aligned to the new cell end; both widths use the active-concentration threshold of 0.6. The 25–60 interval measures average interface advance. Normalized growth is approximately 28% slower at L = 39.',
         sources=[('support/wave_pinning_length_growth.pdf', 'Vector length-growth figure'), ('../../data/runs/length_mechanisms_20261008/comparative_support/wave_pinning_length_growth.json', 'Interface measurements')]),
    dict(id='mca', title='MCA bias predicts withdrawal outcomes in the polarized regime', image='debelly_memory_landmark.png', shape=(11.7, 10.2),
         paragraphs=[
             'In deBelly, membrane–cortex attachments (MCA) retain mechanical polarity memory. Rac and Rho respond rapidly to the gradient, whereas cortical flow progressively redistributes the attachment bias. Attachments favor the rear, opposite the Rac front. The first zero crossing of bound-MCA polarity therefore marks removal of the original attachment bias.',
             'Withdrawal near this landmark produces persistent reversal at L = 23 and 39. At L = 18, the initial state is polarized, but the new bias remains weak after zero-crossing withdrawal. Small-perturbation analysis identifies a short-cell amplification regime in which strong polarity is more difficult to regenerate. MCA zero crossing and polarity speed remain informative measures of memory redistribution, while amplification limits the short-cell withdrawal response.'],
         caption='Figure 8. Rac and bound-MCA responses to withdrawal near MCA zero crossing at three lengths, and length-dependent spontaneous amplification. All withdrawal trajectories begin from qualifying basal-polarized states. Colored traces indicate the withdrawal offset relative to continuous-cue MCA zero. No-cue small-perturbation growth changes sign near L = 19.53. A cue-established polarized branch can exist below this onset.',
         sources=[('debelly_memory_landmark.pdf', 'Vector landmark figure'), ('../deBelly_mechanism/debelly_zero_crossing.pdf', 'Full corrected withdrawal and polarizability analysis')]),
    dict(id='flow', title='Longer MCA displacement and opposing flow delay memory reversal', image='debelly_flow_length.png', shape=(11.7, 10.2),
         paragraphs=[
             'The time to MCA zero crossing increases from 289 to 686 s over L = 23–39 µm and is approximately linear across these nine lengths (R-squared = 0.985). The required nominal MCA centroid displacement increases from 1.36 to 2.49 µm, while mean net centroid speed decreases from approximately 0.0047 to 0.0036 µm/s. Increased redistribution distance is therefore not compensated by faster average progress.',
             'A fixed-profile mechanical control shows that a longer cortex permits a more gradual velocity return and increases isolated Rho-driven flow, until attachment friction screens the coupling. In the native trajectories near MCA zero, however, prescribed-endpoint advection increasingly opposes Rho-driven redistribution. The resulting reduction in late net progress, combined with the longer required displacement, accounts for the observed timing trend.'],
         caption='Figure 9. Continuous-cue MCA-zero times (A), required centroid displacement (B), mean net centroid speed (C), fixed-Rho-profile domain-extension control (D–E), and instantaneous MCA polarity-rate budget at zero (F). Gradient slope is 0.24. Mean speed is displacement divided by zero-crossing time; the mechanical control and rate budget provide separate mechanistic evidence. Endpoint commands are independent. MCA-zero times are distinct from the minimum finite-cue durations in Figure 1.',
         sources=[('debelly_flow_length.pdf', 'Vector length and flow figure'), ('../deBelly_mechanism/debelly_advection.pdf', 'Full corrected flow decomposition'), ('../../data/runs/length_mechanisms_20261008/debelly/figure_audit.json', 'Timing fit and source audit')]),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def esc(value):
    return html.escape(value, quote=True)


def wrapped(fig, text, y, width=130, size=10.5, color='#263d4b', left=.065):
    lines=textwrap.wrap(text, width)
    fig.text(left, y, '\n'.join(lines), va='top', fontsize=size, color=color, linespacing=1.4)
    return y - len(lines)*size*1.4/72/fig.get_figheight()


def verify_sources():
    comparison=ROOT/'data/runs/debelly_corrected_three_test_20261008/results/status.json'
    status=json.loads(comparison.read_text())
    assert status['complete'] and status['producer']=='finished' and status['execution_errors']==0
    initial=ROOT/'data/runs/initial_polarity_length_corrected_20261008'
    provenance=json.loads((initial/'report_provenance.json').read_text())
    assert provenance['all_checkpoint_hashes_pass'] and provenance['qualified_preparations']==97
    assert provenance['source_manifest_sha256']==sha(ROOT/'data/runs/debelly_corrected_three_test_20261008/manifest.json')
    old=list(csv.DictReader((ROOT/'data/runs/initial_polarity_length_20261008/steady_state_metrics.csv').open()))
    new=list(csv.DictReader((initial/'steady_state_metrics.csv').open()))
    assert len(old)==len(new)==97
    for a,b in zip(old,new):
        assert (a['model'],a['domain_length'])==(b['model'],b['domain_length'])
        if a['model']!='debelly': assert a==b
    for checkpoint in provenance['checkpoints']:
        assert sha(ROOT/checkpoint['path'])==checkpoint['sha256']
    comparative=DATA/'comparative_support/provenance.json'
    evidence=json.loads(comparative.read_text())
    for path,digest in evidence['sources'].items():
        assert sha(ROOT/path)==digest
    for path,digest in evidence['figures'].items():
        assert sha(FIG/'support'/path)==digest
    revision=json.loads((DATA/'comparative_revision/provenance.json').read_text())
    assert revision['producer']=='finished' and revision['integration']=='none'
    for path,digest in revision['sources'].items():
        assert sha(ROOT/path)==digest
    for path,digest in revision['figures'].items():
        assert sha(FIG/'comparative_revision'/path)==digest
    assert json.loads((DATA/'comparative_revision/independent_observable_audit.json').read_text())['passed']
    sink=DATA/'sink/audit.json'
    assert json.loads(sink.read_text())['passed']
    debelly=DATA/'debelly/figure_audit.json'
    d=json.loads(debelly.read_text())
    assert d.get('passed', d.get('all_passed', False)), d
    for path,digest in d['source_sha256'].items():
        assert sha(ROOT/path)==digest
    for path,digest in d['figure_sha256'].items():
        assert sha(FIG/path)==digest
    return dict(comparison=status['generated_at'], initial=provenance['generated_at'],
                comparative=revision['generated_at'], sink=json.loads(sink.read_text())['generated_at'],
                debelly=d['generated_at'])


def pdf_report(stamp):
    plt.rcParams.update({'font.family':'DejaVu Sans', 'pdf.fonttype':42})
    path=FIG/'length_mechanisms.pdf'
    with PdfPages(path) as pdf:
        pdf.infodict().update(Title=TITLE, Author='Polarity manuscript working report',
                             Subject='Material redistribution, chemical interfaces, and cortical flow')
        for number,c in enumerate(CHAPTERS, 1):
            w,h=c['shape'];fig=plt.figure(figsize=(w,h))
            fig.text(.065, .975, f'POLARITY REVERSAL · FIGURE {number}', color='#52636b', fontsize=9, weight='bold', va='top')
            title_lines=textwrap.wrap(c['title'], 76)
            fig.text(.065, .947, '\n'.join(title_lines), fontsize=17, weight='bold', va='top', color='#153b4b')
            y=.947-len(title_lines)*22/72/h-.035
            for paragraph in c['paragraphs']:
                y=wrapped(fig, paragraph, y, width=132, size=10.4)-.018
            caption_lines=textwrap.wrap(c['caption'], 145)
            caption_top=.065+len(caption_lines)*9.1*1.4/72/h
            top=y-.010;bottom=caption_top+.035
            available_h=(top-bottom)*h; available_w=.89*w
            with Image.open(FIG/c['image']) as im:ratio=im.width/im.height
            display_w=min(available_w,available_h*ratio)
            display_h=display_w/ratio
            assert display_h>2.7, (c['id'],display_h)
            ax=fig.add_axes([(.5-display_w/w/2), bottom+(available_h-display_h)/2/h, display_w/w, display_h/h])
            ax.imshow(plt.imread(FIG/c['image']));ax.axis('off')
            wrapped(fig, c['caption'], caption_top, width=145, size=9.1, color='#51616b')
            fig.text(.065,.020, 'Generated '+stamp[:19]+' UTC · all source producers finished · sources in the HTML report', fontsize=8, color='#667983')
            fig.text(.94,.020, str(number), ha='right', fontsize=9, color='#29717b')
            pdf.savefig(fig, dpi=180);plt.close(fig)
    return path


def html_report(stamp, source_times):
    nav=''.join(f'<a href="#{c["id"]}">{i:02d} {esc(c["title"])}</a>' for i,c in enumerate(CHAPTERS,1))
    chapters=[]
    for i,c in enumerate(CHAPTERS,1):
        paragraphs=''.join(f'<p>{esc(p)}</p>' for p in c['paragraphs'])
        sources=' · '.join(f'<a href="{esc(url)}">{esc(label)}</a>' for url,label in c['sources'])
        chapters.append(f'<section id="{c["id"]}" aria-labelledby="title-{c["id"]}"><div class="number">{i:02d}</div>'
                        f'<h2 id="title-{c["id"]}">{esc(c["title"])}</h2><div class="prose">{paragraphs}</div>'
                        f'<figure><a href="{esc(c["image"])}" target="_blank" rel="noopener" aria-label="Open full-size figure: {esc(c["title"])}">'
                        f'<img src="{esc(c["image"])}" alt="{esc(c["caption"])}" loading="{"eager" if i==1 else "lazy"}"></a>'
                        f'<figcaption>{esc(c["caption"])}</figcaption></figure><p class="sources">{sources}</p></section>')
    times=''.join(f'<li>{esc(k)}: {esc(v)}</li>' for k,v in source_times.items())
    document=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{TITLE}</title><style>
:root{{--ink:#193847;--muted:#5a707a;--accent:#236e79;--rule:#dbe6e8;--paper:#fff;--wash:#f3f7f7}}
*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;background:var(--paper);color:var(--ink);font:17px/1.65 system-ui,sans-serif}}
main{{max-width:1190px;margin:auto;padding:48px 32px 64px}}header{{padding-bottom:28px;border-bottom:2px solid var(--ink)}}
.eyebrow,.number{{font-size:12px;letter-spacing:.07em;color:var(--muted);font-weight:650;text-transform:uppercase}}h1{{font-size:clamp(28px,4vw,38px);line-height:1.2;margin:12px 0 20px;max-width:1000px}}
h2{{font-size:clamp(23px,3vw,31px);line-height:1.25;margin:10px 0 22px;max-width:1000px}}p{{margin:0 0 16px}}a{{color:var(--accent);text-underline-offset:3px}}a:hover{{color:#12394d}}
.intro{{max-width:850px;font-size:19px}}.metadata{{font-size:13px;color:var(--muted)}}.download{{display:inline-block;margin:10px 12px 16px 0;background:var(--ink);color:white;padding:9px 16px;border-radius:5px;text-decoration:none;font-weight:650}}.download:hover{{color:white;background:var(--accent)}}
nav{{display:grid;grid-template-columns:1fr 1fr;gap:7px 24px;background:var(--wash);padding:20px 24px;margin:24px 0;font-size:14px}}nav a{{text-decoration:none}}
section{{padding:36px 0 32px;border-bottom:1px solid var(--rule);scroll-margin-top:20px}}.prose{{max-width:925px}}figure{{margin:26px 0 14px}}img{{display:block;width:100%;height:auto;background:white}}figcaption{{font-size:14px;line-height:1.55;color:var(--muted);max-width:1040px;margin-top:15px}}
.sources{{font-size:13px;color:var(--muted);line-height:1.7}}details{{padding:20px 24px;background:var(--wash);margin-top:30px}}summary{{cursor:pointer;font-weight:650}}details p{{margin-top:14px;font-size:14px}}details ul{{font-size:13px;padding-left:20px}}footer{{margin-top:24px;font-size:13px;color:var(--muted)}}
@media(max-width:700px){{main{{padding:28px 18px}}nav{{grid-template-columns:1fr;padding:16px}}body{{font-size:16px}}section{{padding-top:28px}}}}
@media(prefers-reduced-motion:reduce){{html{{scroll-behavior:auto}}}}
@media print{{nav,.download,details{{display:none}}main{{padding:0}}section{{break-before:page}}a{{color:inherit}}}}
</style></head><body><main>
<header><div class="eyebrow">Polarity manuscript · 8 October 2026</div><h1>{TITLE}</h1>
<p class="intro">Polarity reversal requires redistribution of active material and reorganization of polarity memory. The analyses below relate cell length dependence to pole geometry, reaction kinetics, interface motion and cortical flow.</p>
<a class="download" href="length_mechanisms.pdf">Report PDF ({len(CHAPTERS)} pages)</a> <a href="../index.html">Figure collection</a>
<p class="metadata">Generated {esc(stamp)}. Simulation and figure production completed. Full-resolution figures are linked from each image.</p></header>
<nav aria-label="Report chapters">{nav}</nav>{''.join(chapters)}
<details><summary>Measurements, current data and reproduction</summary>
<p>The full duration comparison retains five models and both outcome gates. The focused mechanism chapters explain Goryachev, Wave-pinning and deBelly. All current deBelly figures use independent local endpoint commands; the deprecated centered extension contributes no current evidence here. Endpoint constraints and quasistatic mechanics remain the corrected model assumptions.</p>
<p>The new compact figures are derived from audited native arrays, with no new simulation campaign. Sink geometry, inactive removal and fixed-profile flow controls are explicitly different perturbations. Equal-cue activation is a local initiation measure; immediate restoring gain is a reaction response per small damage. Neither assigns additive fractions of reversal time.</p>
<p><a href="../../THREE_TESTS.md">Three-test protocol</a> · <a href="../../ASSUMPTIONS.md">Model assumptions</a> · <a href="../../data/runs/length_mechanisms_20261008/report_provenance.json">Report hashes and verification</a> · <a href="../../data/runs/length_mechanisms_20261008/debelly/figure_audit.json">Corrected MCA and flow audit</a></p>
<p>Source-data generation times (UTC):</p><ul>{times}</ul>
<p>Rebuild with the existing scaling_sandbox Python environment: <code>code/mechanisms/length_mechanisms_report_20261008.py</code>. The original diagnostic set remains available through the <a href="../three_test_mechanisms/historical_index_20261006.html">historical reader</a>; its older deBelly results are deprecated.</p>
</details><footer>Five-model finite-cue comparison and mechanistic analyses of polarity reversal.</footer>
</main></body></html>'''
    (FIG/'index.html').write_text(document)


def compatibility():
    old=ROOT/'figures/three_test_mechanisms'
    for current,historical in [('index.html','historical_index_20261006.html'),
                               ('three_test_mechanisms.pdf','historical_three_test_mechanisms_20261006.pdf')]:
        if not (old/historical).exists():
            shutil.copy2(old/current, old/historical)
    history=old/'historical_index_20261006.html'
    content=history.read_text()
    content=content.replace('href="three_test_mechanisms.pdf"', 'href="historical_three_test_mechanisms_20261006.pdf"')
    banner='<div style="padding:18px;background:#fff0d8;border:1px solid #d7ba88">Historical 6 October reader. Its centered deBelly results are deprecated. <a href="../length_mechanisms_20261008/index.html">Read the current replacement report</a>.</div>'
    if banner not in content:
        content=content.replace('<body>', '<body>'+banner, 1) if '<body>' in content else content.replace('<h1>', banner+'<h1>', 1)
        history.write_text(content)
    target='../length_mechanisms_20261008/index.html'
    (old/'index.html').write_text(f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{TITLE}</title><style>body{{max-width:850px;margin:55px auto;padding:0 24px;font:18px/1.65 system-ui;color:#193847}}a{{color:#236e79}}iframe{{width:100%;height:85vh;border:1px solid #dbe6e8}}</style><h1>{TITLE}</h1><p><a href="{target}">Open the new report</a> · <a href="three_test_mechanisms.pdf">{len(CHAPTERS)}-page replacement PDF</a></p><p>The replacement follows material redistribution, Goryachev versus Wave-pinning, and corrected deBelly cortical flow. <a href="historical_index_20261006.html">Historical reader and diagnostics</a>.</p><iframe src="{target}" title="Current polarity length mechanism report"></iframe></html>')
    shutil.copy2(FIG/'length_mechanisms.pdf', old/'three_test_mechanisms.pdf')


def main():
    FIG.mkdir(parents=True, exist_ok=True);DATA.mkdir(parents=True, exist_ok=True)
    stamp=datetime.now(timezone.utc).isoformat()
    source_times=verify_sources()
    for src,name in [(ROOT/'figures/supplement_width_response.png','length_responses.png'),
                     (ROOT/'figures/initial_polarity_length/initial_polarity_length.png','initial_polarity_length.png')]:
        shutil.copy2(src, FIG/name)
    for c in CHAPTERS:assert (FIG/c['image']).exists(), c['image']
    pdf=pdf_report(stamp);html_report(stamp,source_times);compatibility()
    images={c['image']:sha(FIG/c['image']) for c in CHAPTERS}
    sources={str(p.relative_to(ROOT)):sha(p) for p in [
        ROOT/'data/runs/debelly_corrected_three_test_20261008/results/status.json',
        ROOT/'data/runs/initial_polarity_length_corrected_20261008/report_provenance.json',
        DATA/'comparative_support/provenance.json', DATA/'comparative_revision/provenance.json',
        DATA/'comparative_revision/independent_observable_audit.json',
        DATA/'sink/audit.json', DATA/'debelly/figure_audit.json']}
    record=dict(generated_at=stamp,producer='finished',pages=len(CHAPTERS),renderer_sha256=sha(Path(__file__)),
        source_data_times=source_times,source_hashes=sources,figure_hashes=images,
        html_sha256=sha(FIG/'index.html'),pdf_sha256=sha(pdf),
        revision='Scientific exposition and chemical-profile/restoration figure revision',
        checks=dict(corrected_initial_checkpoint_hashes=True,all_85_other_model_rows_unchanged=True,
                    audited_comparative_source_hashes=True,corrected_debelly_only=True),
        chapters=[dict(number=i,id=c['id'],title=c['title'],image=c['image']) for i,c in enumerate(CHAPTERS,1)])
    (DATA/'report_provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(html=str(FIG/'index.html'),pdf=str(pdf),pages=len(CHAPTERS),generated_at=stamp),indent=2))


if __name__=='__main__':
    main()
