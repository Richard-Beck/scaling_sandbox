"""Build the concise replacement reader from audited current figures."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import html
import json
import shutil
import io
import re
import sys
from pypdf import PdfReader, PdfWriter, Transformation
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'figures/length_mechanisms_20261008'
DATA = ROOT / 'data/runs/length_mechanisms_20261008'
TITLE = 'Mechanisms of cell-length-dependent polarity reversal'
sys.path.insert(0,str(ROOT/'code/mechanisms'))
from length_dev_output_layout import dev_output_path


def figure_artifact(image, filename):
    """Find a figure's companion export in its declared development layout."""
    image=Path(image)
    return dev_output_path(image.parent,filename) if image.parent==FIG/'dev' else image.parent/filename


def manifest_artifact(path, record, filename):
    """Resolve exports relative to the manifest's explicitly declared artifact root."""
    candidate=Path(filename)
    if candidate.is_absolute():return candidate
    local=path.parent/record.get('artifact_root','.')/candidate
    return local if local.exists() else FIG/candidate

CHAPTERS = [{'id': 'responses',
  'title': 'Cell length changes the work needed to reverse polarity',
  'image': 'review_20261009_quality/length_responses.png',
  'shape': (11.7, 20.5),
  'paragraphs': ['Minimum reversal cue duration depends on length and cue geometry (Fig. 1A–O). Shaded bands '
                 'mark decreases between adjacent successful samples. Fixed-width patches strongly slow '
                 'Wave-pinning and Holmes (Fig. 1A,J); Goryachev slows under all three cues (Fig. 1D–F). The '
                 'decreasing gradient arms in Wave-pinning, Holmes and Spring therefore need explanation '
                 '(Fig. 1C,L,O).',
                 'A constant physical slope increases cue peak with length. Holding peak at its '
                 'reference-length value prevents reversal at 2/4 sampled native doses in Wave-pinning, 2/3 '
                 'in Holmes and 3/4 in Spring. Removing Spring extension feedback prevents 4/4. These '
                 'controls support contributions from input growth and mechanics; they do not locate new '
                 'duration thresholds. For deBelly, we follow the length-dependent resting mechanics, '
                 'test its tension feedback, and track how a successful cue rewrites attachment memory '
                 '(Fig. 1G–I; Figs. 8–12).',
                 'The mean earth mover distance between each initial active distribution and its reflected '
                 'target grows with length in all five models (Fig. 1P–T). Larger cells thus have a longer '
                 'spatial reorganization target as well as more active material. Amount times distance '
                 'measures the burden of physically transporting that distribution; conversion and turnover '
                 'can redistribute activity through other routes.',
                 'Each pulse is followed by 10000 basal recovery units. Orientation reversal requires a '
                 'normalized active first moment above 0.05, contrast at least 0.05 and orientation range '
                 'below 0.005 throughout the final 500 units. The strength gate additionally restores the '
                 'original scalar polarity magnitude within 1%. Open markers flag a separate native-field '
                 'stability failure. Ineligible starts and unreached thresholds remain visible.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/length_responses.pdf', 'Vector figure'),
              ('review_20261009_quality/length_responses.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../../data/runs/length_mechanisms_20261008/review_20261009_quality/negative_intervals.csv',
               'Every decreasing interval'),
              ('../../data/runs/three_test_mechanisms_20261006/results/cell_evidence.json',
               'Matched-peak and feedback controls')]},
 {'id': 'material',
  'title': 'Two starting-state differences distinguish Goryachev',
  'image': 'review_20261009_quality/initial_polarity_length.png',
  'shape': (11.7, 12),
  'paragraphs': ['Longer cells contain more active pole material in every model (Fig. 2A). Goryachev '
                 'accommodates it mainly through increased height at approximately fixed width; the other '
                 'four predominantly widen their poles (Fig. 2B,C). Goryachev also depletes inactive '
                 'concentration as length increases, whereas the explicit inactive pools of Wave-pinning, '
                 'Holmes and Spring remain approximately constant (Fig. 2D). deBelly has no explicit '
                 'inactive Rac field.',
                 'We examine how these two differences—active geometry and inactive availability—affect '
                 'new-pole growth, incumbent resistance and lasting reversal. Wave-pinning is the simplest '
                 'comparable two-field mass-conserved model among the other four. Both models face '
                 'increasing material and reflected-target distance (Fig. 1P,Q), but height growth and '
                 'domain widening impose different kinetic demands.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/initial_polarity_length.pdf', 'Vector figure'),
              ('review_20261009_quality/initial_polarity_length.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../../data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv',
               'Unchanged starting-state measurements')]},
 {'id': 'initiation',
  'title': 'Inactive depletion and activation kinetics hinder new Goryachev peaks',
  'image': 'review_20261009_quality/chemical_hurdles.png',
  'shape': (11.7, 20.5),
  'paragraphs': ['Goryachev’s taller incumbent accompanies inactive depletion throughout the cell (Fig. 3A). '
                 'Opposite-end inactive concentration falls from 0.263 to 0.088 between L = 9 and 39. Added '
                 'activation equals cue strength times inactive concentration, so an equal local cue '
                 'initially activates only one third as much material. Wave-pinning’s inactive availability '
                 'changes little (Fig. 3B).',
                 'The rate maps retain absolute native reaction amplitudes; the overlaid curves trace paired '
                 'active/inactive concentrations across each spatial pole (Fig. 3C,D). Goryachev’s '
                 'unsaturated feedback supports higher active concentrations despite inactive depletion. '
                 'Peak basal net activation along its incumbent rises from 0.212 to 0.639 over L = 9–39; '
                 'Wave-pinning’s remains near 0.0265 and moves outward with the interface (Fig. 3E,F). '
                 'Diffusion balances these local reaction rates in the stationary profile.',
                 'For a competing peak, freeze inactive concentration at its initial opposite-end value and '
                 'leave the incumbent unchanged. Without cue, Goryachev’s unstable local activation '
                 'threshold then rises from 18.1 to 40.6 to 58.4 at L = 9, 21 and 39 (Fig. 3G). Its '
                 'negative-rate region extends to higher active concentration, and growth at a given '
                 'concentration weakens. Wave-pinning’s curves and threshold near 0.602 nearly coincide '
                 '(Fig. 3H). Inactive depletion therefore hinders both cue initiation and subsequent '
                 'cue-free growth. This approximation isolates local chemistry; spatial diffusion and pool '
                 'consumption still govern actual nucleation and survival.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/chemical_hurdles.pdf', 'Vector figure'),
              ('review_20261009_quality/chemical_hurdles.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../../data/runs/length_mechanisms_20261008/review_20261009_quality/chemical_measurements.csv',
               'Reaction rates and local thresholds')]},
 {'id': 'restoration',
  'title': 'The taller incumbent also has a stronger immediate restoring response',
  'image': 'review_20261009_quality/incumbent_restoration.png',
  'shape': (11.7, 14),
  'paragraphs': ['The strong incumbent reaction in Fig. 3E motivates a damage-and-return experiment. Remove '
                 '0.5% of the entire inactive profile for 0.1 native units, then return it pointwise; native '
                 'reactions and diffusion continue throughout. Active-core loss measures susceptibility '
                 '(Fig. 4A,B); the immediate reaction after return measures rebound per amount lost (Fig. '
                 '4C,D).',
                 'From L = 9 to 39, Goryachev’s fractional damage decreases by approximately half, while '
                 'Wave-pinning’s stays nearly constant (Fig. 4E). Immediate restoring reaction per unit '
                 'damage grows 15.7-fold in Goryachev and 1.3-fold in Wave-pinning (Fig. 4F). Together with '
                 'Fig. 3G,H, this distinguishes stronger incumbent resistance from weaker competitor growth. '
                 'Immediate rebound can coexist with slower late recovery; the coupled effects do not assign '
                 'additive fractions of reversal duration.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/incumbent_restoration.pdf', 'Vector figure'),
              ('review_20261009_quality/incumbent_restoration.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../Goryachev_wave_pinning_trial/inactive_pool_perturbation/index.html',
               'Full removal-and-return experiment')]},
 {'id': 'sink',
  'title': 'Separate fixed-gap erosion from the extra separation penalty',
  'image': 'review_20261009_quality/sink_distance_and_resistance.png',
  'shape': (11.7, 20.5),
  'paragraphs': ['At t = 0, switch on a width-1 inactive sink for the entire 20000-unit integration, without '
                 'reversal cue. Its equation is dv/dt = Dv × d²v/dx² − f(u,v) − k × χsink(x) × v; removed '
                 'material never returns. Compare the opposite-end support [L−1,L] with a support beginning '
                 '3 units beyond the initial half-height pole edge. Both remain fixed during erosion.',
                 'Use k = 1000 inverse native time units to approach the measured strength plateau: '
                 'increasing k from 100 to 1000 changes selected loss rates by less than 7%. This makes sink '
                 'chemistry less limiting and exposes transport and pole response. It is the largest tested '
                 'strength, rather than a universal saturation constant. Peak-height and active-amount loss '
                 'rates share model-comparison axes (Fig. 5A,B).',
                 'At fixed gap, the 10–30% active-loss interval grows from 15.8 to 77.6 in Goryachev and '
                 '0.646 to 1.92 in Wave-pinning over L = 9–39 (4.9-fold and 3.0-fold; Fig. 5C). At L = 39, '
                 'the opposite-end sink extends these intervals to 836 and 8.74: extra separation adds 759 '
                 'and 6.82 native time units, respectively (Fig. 5D). Active/inactive snapshots show how '
                 'Goryachev loses height while Wave-pinning retreats its interface (Fig. 5E–H), consistent '
                 'with Fig. 3E,F. The paired difference quantifies the extra-separation penalty alongside '
                 'the residual fixed-gap effect.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/sink_distance_and_resistance.pdf', 'Vector figure'),
              ('review_20261009_quality/sink_distance_and_resistance.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../inactive_sink_loss/index.html', 'Full sink-strength sweep and audited histories'),
              ('../../data/runs/length_mechanisms_20261008/review_20261009_quality/sink_distance_contrast.csv',
               'Fixed-gap interval and extra interval'),
              ('../../data/runs/length_mechanisms_20261008/review_20261009_quality/sink_saturation.csv',
               'Plateau check'),
              ('../../data/runs/length_mechanisms_20261008/review_20261009_quality/sink_snapshots.csv',
               'Every profile snapshot time')]},
 {'id': 'interface',
  'title': 'A fixed-width seed must expand to retain Wave-pinning reversal',
  'image': 'review_20261009_quality/wave_pinning_withdrawal.png',
  'shape': (11.7, 12.8),
  'paragraphs': ['At L = 15, a fixed physical width-3 cue seeds a new Wave-pinning peak narrower than the '
                 'eventual polarized domain. A 20-unit pulse raises its height while leaving much of the '
                 'incumbent intact; a 60-unit pulse allows further interface advance and incumbent retreat '
                 '(Fig. 6A,B). The nearly fixed plateau and interface-localized reactions in Fig. 3B,F '
                 'explain why height alone cannot replace domain growth.',
                 'After the short pulse, the new peak collapses and the original pole returns (Fig. 6C). '
                 'After the longer pulse, height relaxes while the new domain widens beyond cue support and '
                 'reaches a fully relaxed reversed polarized state (Fig. 6D). Both final fields match the '
                 'original or reflected basal state within 0.1% of fixed scales, with negligible '
                 'final-window drift. These examples bracket the 35.5-unit strength threshold.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/wave_pinning_withdrawal.pdf', 'Vector figure'),
              ('review_20261009_quality/wave_pinning_withdrawal.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../../data/runs/length_revision_20261009/independent_audit.json',
               'Fresh checkpoint and relaxed-endpoint audit')]},
 {'id': 'growth',
  'title': 'Goryachev delays redistribution, then switches mass abruptly',
  'image': 'review_20261009_quality/half_mass_dynamics.png',
  'shape': (11.7, 13),
  'paragraphs': ['The active half-mass coordinate x50 bisects the spatial active integral. Under a '
                 'continuously maintained relative-20% cue, Goryachev shows a length-dependent delay '
                 'followed by abrupt redistribution; Wave-pinning moves more gradually (Fig. 7A,B). Halfway '
                 'crossings toward the reflected basal position occur at 70.5, 410.5 and 2222 in Goryachev '
                 'at L = 9, 21 and 39, versus 20.5, 51 and 135 in Wave-pinning.',
                 'Initial half-mass velocity falls from 0.0356 to 0.0120 across these Goryachev lengths, but '
                 'rises from 0.0878 to 0.314 in Wave-pinning (Fig. 7C,D). Goryachev’s delay ends in a narrow '
                 'velocity pulse; Wave-pinning has broader motion, later slowing and secondary transients. '
                 'The difference connects the competitor kinetics in Fig. 3G,H and incumbent response in '
                 'Fig. 4 to persistent domain expansion in Fig. 6.',
                 'For two separated narrow poles, x50 rapidly traverses the low-concentration gap when the '
                 'new pole acquires half the active mass. This geometric median shift can precede '
                 'commitment: at L = 39, Goryachev crosses halfway at 2222, but requires a 2680-unit '
                 'orientation pulse (Fig. 7A). Velocity is a 0.5-unit window mean, rather than an interface '
                 'speed; all plotted cues remain on.'],
  'caption': '__LOAD_CAPTION__',
  'sources': [('review_20261009_quality/half_mass_dynamics.pdf', 'Vector figure'),
              ('review_20261009_quality/half_mass_dynamics.manifest.json',
               'Figure manifest, legend, hashes and reproduction'),
              ('../../data/runs/length_mechanisms_20261008/review_20261009_quality/half_mass_measurements.csv',
               'Half-mass timing and velocity measurements'),
              ('../../data/runs/length_revision_20261009/manifest.json', 'Exact replay inputs')]},
 {'id': 'mca', 'title': 'MCA polarity tracks the slow reversal timescale'},
 {'id': 'flow', 'title': 'MCA travel distance and speed explain the continuous-cue clock'},
 {'id': 'screening', 'title': 'A spring/dashpot picture makes mechanical screening visible'},
 {'id': 'patch_turnovers', 'title': 'Patch-duration dips reflect competing redistribution trends'}]


# Each specialist renderer supplies the scientific text associated with its new panels.
for override_path in [DATA/'quality_20261009_debelly/chapter_overrides.json',
                      DATA/'quality_20261009_mechanics/chapter_overrides.json']:
    overrides=json.loads(override_path.read_text())
    for chapter in CHAPTERS:
        if chapter['id'] in overrides:
            chapter.update(overrides[chapter['id']])

# Preserve the manuscript figure numbers while adding readable accompaniment pages.
for number,chapter in enumerate(CHAPTERS,1):chapter['figure_label']=str(number)
field_overrides=json.loads((DATA/'quality_20261009_field_accompaniment/chapter_overrides.json').read_text())
expanded=[]
for chapter in CHAPTERS:
    if chapter['id'] in field_overrides:chapter.update(field_overrides[chapter['id']])
    if chapter['id']=='flow':chapter['paragraphs']=chapter['paragraphs'][:2]
    if chapter['id']=='mca':
        chapter['paragraphs'][-1]+=' The complete L = 29 relative-patch field and endpoint histories accompany these traces in Fig. 8S1–8S2.'
    expanded.append(chapter)
    if chapter['id']=='mca':expanded.extend(field_overrides[k] for k in ['mca_fields','mca_globals'])
    if chapter['id']=='flow':expanded.append(field_overrides['flow_budget'])
CHAPTERS=expanded

# The narrative revision replaces the deBelly section as a complete sequence.
# Earlier override files above remain the provenance of the archived reader.
narrative={}
for section in ['resting_ablation','mechanics','memory','fields']:
    narrative.update(json.loads((DATA/'narrative_20261009'/section/'chapter_overrides.json').read_text()))
CHAPTERS=CHAPTERS[:7]
for key,label in [('resting','8'),('screening','9'),('ablation','10'),
                  ('mca','11'),('mca_recovery','11S1'),('mca_fields','11S2'),('mca_globals','11S3'),
                  ('flow','12'),('flow_budget','12S1')]:
    if key not in narrative:
        assert key=='flow_budget',key
        continue
    chapter=dict(narrative[key],id=key,figure_label=label)
    if key=='mca_recovery':chapter['paragraphs']=chapter['paragraphs'][:1]
    if key=='flow':
        chapter['paragraphs']=[
            'Under the uninterrupted constant-slope cue, increasing length from 23 to 39 µm raises the initial MCA centroid distance to the midpoint from 1.36 to 2.49 µm. Its mean distribution speed falls from 4.71 to 3.63 × 10⁻³ µm/s, so the spatial-zero clock grows from 289 to 686 s. The longer attachment timescale thus combines a greater redistribution distance with a slower mean rate. The endpoint-MCA clock grows from 258 to 491 s, exposing a distinct local exchange timescale.',
            'Across the same lengths, the shortest strength-passing gradient pulse increases from 150 to 293 s. It can end before the MCA crossing because redistribution continues after withdrawal. The attachment clocks therefore describe a slower memory process accompanying reversal; the finite pulse and subsequent basal recovery together rebuild the reversed state. The paired source budgets show how endpoint forcing and Rho-driven cortical flow contribute during stimulation.']
    CHAPTERS.append(chapter)

# Load the standalone scientific legends for every revised figure.
for chapter in CHAPTERS:
    caption_path=figure_artifact(FIG/chapter['image'],Path(chapter['image']).stem+'_caption.txt')
    if caption_path.exists() and not chapter.get('caption_override'):chapter['caption']=caption_path.read_text().strip()
    assert chapter['caption']!='__LOAD_CAPTION__',chapter['id']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def esc(value):
    return html.escape(value, quote=True)


def verify_sources():
    """Verify the active data and attached figures, independent of historical renderers."""
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
    for previous,current in zip(old,new):
        assert (previous['model'],previous['domain_length'])==(current['model'],current['domain_length'])
        if previous['model']!='debelly': assert previous==current
    for checkpoint in provenance['checkpoints']:
        assert sha(ROOT/checkpoint['path'])==checkpoint['sha256']
    for run in ['length_revision_20261009','length_debelly_review_20261009']:
        assert json.loads((ROOT/'data/runs'/run/'independent_audit.json').read_text())['passed']
    for audit in [DATA/'review_20261009_quality/quality_audit.json',
                  DATA/'quality_20261009_debelly/independent_audit.json',
                  DATA/'quality_20261009_debelly/execution_source_audit.json',
                  DATA/'quality_20261009_mechanics/audit.json',
                  DATA/'quality_20261009_field_accompaniment/endpoint_story_audit.json']:
        result=json.loads(audit.read_text())
        assert result.get('passed',result.get('frozen_sources_verified',False)),audit
    helper=ROOT/'code/mechanisms/length_dev_output_layout.py'
    verified={str(helper.relative_to(ROOT)):sha(helper)};active_manifests=[]
    for chapter in CHAPTERS:
        stem=(FIG/chapter['image']).with_suffix('')
        path=figure_artifact(stem,stem.name+'.manifest.json');record=json.loads(path.read_text())
        assert record.get('producer','finished')=='finished',path
        assert record.get('passed',True),path
        renderer_name=chapter.get('renderer') or {'review_20261009_quality':'length_revision_20261009_figures.py',
            'quality_20261009_debelly':'length_quality_debelly_traces_20261009.py',
            'quality_20261009_mechanics':'length_quality_mechanics_20261009.py'}[path.parent.name]
        renderer=ROOT/'code/mechanisms'/renderer_name
        if 'renderer_sha256' in record:assert sha(renderer)==record['renderer_sha256'],renderer
        verified[str(renderer.relative_to(ROOT))]=sha(renderer)
        exports=record.get('sha256',record.get('figure_sha256',{}))
        assert exports,path
        for filename,digest in exports.items():
            candidate=manifest_artifact(path,record,filename)
            assert candidate.exists() and sha(candidate)==digest,(path,candidate)
            verified[str(candidate.relative_to(ROOT))]=digest
        for filename,digest in record.get('source_sha256',record.get('sources_sha256',record.get('sources',{}))).items():
            candidate=Path(filename)
            if not candidate.is_absolute():candidate=ROOT/candidate
            assert candidate.exists() and sha(candidate)==digest,(path,candidate)
            verified[str(candidate.relative_to(ROOT))]=digest
        for filename,digest in record.get('data_sha256',{}).items():
            candidate=manifest_artifact(path,record,filename) if 'artifact_root' in record else DATA/path.parent.name/filename
            if not candidate.exists():candidate=path.parent/filename
            if not candidate.exists():candidate=ROOT/filename
            if not candidate.exists():candidate=DATA/path.parent.relative_to(FIG)/filename
            assert candidate.exists() and sha(candidate)==digest,(path,candidate)
            verified[str(candidate.relative_to(ROOT))]=digest
        for destination,_ in chapter['sources']:
            if destination.startswith(('http:','https:')):continue
            candidate=FIG/destination
            assert candidate.is_file(),('missing chapter source',destination)
            verified[str(candidate.relative_to(ROOT))]=sha(candidate)
        assert figure_artifact(stem,stem.name+'.pdf').exists(),stem
        active_manifests.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
    archived=FIG/'archive_20261008'
    original=json.loads((archived/'report_provenance.json').read_text())
    assert sha(archived/'length_mechanisms_report_20261008.py')==original['renderer_sha256']
    assert sha(archived/'length_mechanisms.pdf')==original['pdf_sha256']
    return dict(comparison=status['generated_at'],initial=provenance['generated_at'],
                active_manifests=active_manifests,verified_hashes=verified)


def pdf_report(stamp):
    """Assemble text pages and original vector plots without rasterizing the plots."""
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
    path=FIG/'length_mechanisms.pdf';writer=PdfWriter();layouts=[]
    for number,chapter in enumerate(CHAPTERS,1):
        image=FIG/chapter['image'];vector=figure_artifact(image,image.stem+'.pdf')
        source=PdfReader(vector)
        assert len(source.pages)==1
        figure_page=source.pages[0]
        fw,fh=float(figure_page.mediabox.width),float(figure_page.mediabox.height)
        page_w=842.4;margin=34.;text_w=page_w-2*margin;plot_w=text_w
        title_lines=textwrap.wrap(chapter['title'],72)
        body_lines=[textwrap.wrap(paragraph,132) for paragraph in chapter['paragraphs']]
        legend_lines=textwrap.wrap(chapter['caption'],145)
        measured_widths=[]
        for lines,size,weight in [(title_lines,15,'bold'),([line for paragraph in body_lines for line in paragraph],10.4,'normal'),(legend_lines,9,'normal')]:
            font=FontProperties(family='DejaVu Sans',weight=weight)
            for line in lines:
                width=TextPath((0,0),line,size=size,prop=font).get_extents().width
                assert width<=text_w,(chapter['id'],'text exceeds page width',width,text_w)
                measured_widths.append(float(width))
        top_h=34+15+19*len(title_lines)+18+sum(14.6*len(lines)+10 for lines in body_lines)
        plot_h=plot_w*fh/fw
        legend_h=12.6*len(legend_lines)
        page_h=top_h+20+plot_h+20+legend_h+48
        fig=plt.figure(figsize=(page_w/72,page_h/72))
        def put(text,y,size,color='#263d4b',weight='normal',spacing=1.4):
            fig.text(margin/page_w,1-y/page_h,text,va='top',fontsize=size,
                     color=color,weight=weight,linespacing=spacing)
        put(f'POLARITY REVERSAL · FIGURE {chapter["figure_label"]}',25,8.5,'#52636b','bold')
        put('\n'.join(title_lines),49,15,'#153b4b','bold',1.25)
        y=49+19*len(title_lines)+17
        for lines in body_lines:
            put('\n'.join(lines),y,10.4);y+=14.6*len(lines)+10
        plot_top=y+10;plot_bottom=plot_top+plot_h
        put('\n'.join(legend_lines),plot_bottom+20,9.0,'#51616b')
        fig.text(margin/page_w,17/page_h,'Generated '+stamp[:19]+' UTC · producers finished · methods and manifests in HTML',
                 fontsize=7.5,color='#667983',url='index.html#'+chapter['id'])
        fig.text((page_w-margin)/page_w,17/page_h,str(number),ha='right',fontsize=8,color='#29717b')
        buf=io.BytesIO();fig.savefig(buf,format='pdf');plt.close(fig);buf.seek(0)
        page=PdfReader(buf).pages[0]
        scale=plot_w/fw;x=margin-float(figure_page.mediabox.left)*scale
        ypdf=page_h-plot_bottom-float(figure_page.mediabox.bottom)*scale
        page.merge_transformed_page(figure_page,Transformation().scale(scale).translate(x,ypdf))
        writer.add_page(page)
        assert plot_bottom+20+legend_h<page_h-25
        layouts.append(dict(page=number,figure=chapter['figure_label'],page_width_points=page_w,page_height_points=page_h,
            figure_bounds=[margin,page_h-plot_bottom,margin+plot_w,page_h-plot_top],
            figure_vector_source=str(vector.relative_to(ROOT)),
            body_lines=sum(map(len,body_lines)),legend_lines=len(legend_lines),
            maximum_text_width_points=max(measured_widths),text_area_width_points=text_w))
    writer.add_metadata({'/Title':TITLE,'/Author':'Polarity manuscript working report',
        '/Subject':'Pole geometry, activation kinetics, redistribution and cortical mechanics'})
    with path.open('wb') as stream:writer.write(stream)
    (DATA/'quality_report_layout.json').write_text(json.dumps(dict(generated_at=stamp,
        vector_figures=True,pages=layouts),indent=2)+'\n')
    return path


def html_report(stamp, source_times):
    nav=''.join(f'<a href="#{c["id"]}">{esc(c["figure_label"])} {esc(c["title"])}</a>' for c in CHAPTERS)
    chapters=[]
    for i,c in enumerate(CHAPTERS,1):
        paragraphs=''.join(f'<p>{esc(p)}</p>' for p in c['paragraphs'])
        links=c['sources']+[(Path(c['image']).stem+'.manifest.json','Attached report figure manifest')]
        sources=' · '.join(f'<a href="{esc(url)}">{esc(label)}</a>' for url,label in links)
        chapters.append(f'<section id="{c["id"]}" aria-labelledby="title-{c["id"]}"><div class="number">Figure {esc(c["figure_label"])}</div>'
                        f'<h2 id="title-{c["id"]}">{esc(c["title"])}</h2><div class="prose">{paragraphs}</div>'
                        f'<figure><a href="{esc(c["image"])}" target="_blank" rel="noopener" aria-label="Open full-size figure: {esc(c["title"])}">'
                        f'<img src="{esc(c["image"])}" alt="{esc(c["caption"])}" loading="{"eager" if i==1 else "lazy"}"></a>'
                        f'<figcaption>{esc(c["caption"])}</figcaption></figure><p class="sources">{sources}</p></section>')
    times=''.join(f'<li>{esc(k)}: {esc(v)}</li>' for k,v in source_times.items() if k in ['comparison','initial'])
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
<header><div class="eyebrow">Polarity manuscript · figure and mechanics revision, 9 October 2026</div><h1>{TITLE}</h1>
<p class="intro">Cell length changes both the polarized state a cue must overturn and the process of building its replacement. Pole geometry and reaction kinetics explain chemical resistance; in deBelly, tension feedback amplifies resting mechanical differences, and MCA attachment redistribution supplies a slower memory of polarity.</p>
<a class="download" href="length_mechanisms.pdf">Report PDF ({len(CHAPTERS)} vector pages)</a> <a href="../index.html">Figure collection</a>
<p class="metadata">Report and figure revision generated {esc(stamp)}. Core comparison data generated {esc(source_times['comparison'])}; all producers finished. Figure-specific native-data sources and generation records accompany each image.</p></header>
<nav aria-label="Report chapters">{nav}</nav>{''.join(chapters)}
<details><summary>Measurements, current data and reproduction</summary>
<p>The full duration comparison retains five models and both outcome gates. The focused mechanism chapters explain Goryachev, Wave-pinning and deBelly. All current deBelly figures use independent local endpoint commands; the deprecated centered extension contributes no current evidence here. Endpoint constraints and quasistatic mechanics remain the corrected model assumptions.</p>
<p>Figures 1–7 retain their audited chemical evidence. Figures 8–10 establish the resting deBelly mechanical differences and the contribution of tension feedback. Figures 11–12 follow successful strength-passing pulses and uninterrupted-cue clocks. The accompanying pages retain all spatial fields and explain endpoint relaxation. Attached manifests identify checkpoints, histories, field definitions and reproduction commands. Equal-command mechanical controls hold chemistry fixed; the tension-to-Rho ablation independently evolves the coupled model. Instantaneous MCA budgets decompose the source at a saved state. These experiments answer complementary questions.</p>
<p>The demonstrated mechanical feedback contribution and the measured attachment redistribution provide a positive account of several length effects. Their precise contribution to each minimum pulse duration, including the patch-duration dips, remains an open question. Scalar strength replacement and the separate native-field stability flag retain their protocol definitions.</p>
<p><a href="../../THREE_TESTS.md">Three-test protocol</a> · <a href="../../ASSUMPTIONS.md">Model assumptions</a> · <a href="../../data/runs/length_mechanisms_20261008/report_provenance.json">Report hashes and verification</a></p>
<p>Core comparison and starting-state generation times (UTC; figure-specific times are in the attached manifests):</p><ul>{times}</ul>
<p>Rebuild with the existing scaling_sandbox Python environment using the commands in each figure manifest, then run <code>code/mechanisms/length_mechanisms_report_20261008.py</code>. <a href="report_text.txt">Plain-text report</a> · <a href="archive_20261009_before_narrative_revision/index.html">Archived preceding report</a>. The original diagnostic set remains available through the <a href="../three_test_mechanisms/historical_index_20261006.html">historical reader</a>; its older deBelly results are deprecated.</p>
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


def audit_report(pdf):
    """Check page content, attached hashes and local HTML destinations."""
    reader=PdfReader(pdf);assert len(reader.pages)==len(CHAPTERS)
    for number,(page,chapter) in enumerate(zip(reader.pages,CHAPTERS),1):
        content=page.extract_text()
        assert f'FIGURE {chapter["figure_label"]}' in content,(number,'missing title')
        assert len(content)>500,(number,'missing searchable text')
        assert '/Contents' in page and '/Font' in page['/Resources']
    document=(FIG/'index.html').read_text()
    assert 'http-equiv="refresh"' not in document.lower()
    destinations=re.findall(r'(?:href|src)="([^"#]+)"',document)
    for destination in destinations:
        if destination.startswith(('http:','https:','mailto:')):continue
        target=FIG/destination.split('#',1)[0]
        assert target.exists(),('broken report link',destination)
    manifests=set();figure_files=set()
    for destination in destinations:
        if not destination.endswith('.manifest.json'):continue
        manifest_path=FIG/destination;manifest=json.loads(manifest_path.read_text());manifests.add(str(manifest_path))
        for filename,digest in manifest.get('sha256',manifest.get('figure_sha256',{})).items():
            candidate=manifest_artifact(manifest_path,manifest,filename)
            assert candidate.exists() and sha(candidate)==digest,(manifest_path,candidate)
            figure_files.add(str(candidate.resolve()))
    return dict(passed=True,pages=len(reader.pages),vector_figures=True,
                searchable_text=True,local_links_checked=len(destinations),
                linked_manifests_checked=len(manifests),figure_artifacts_checked=len(figure_files),no_timed_refresh=True)


def main():
    FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now(timezone.utc).isoformat();source_times=verify_sources()
    for number,chapter in enumerate(CHAPTERS,1):
        image=FIG/chapter['image'];stem=image.with_suffix('')
        renderer_manifest=figure_artifact(image,stem.name+'.manifest.json')
        original=json.loads(renderer_manifest.read_text())
        artifacts=[p for p in [image,figure_artifact(image,stem.name+'.pdf'),figure_artifact(image,stem.name+'.svg')] if p.exists()]
        attached=dict(number=chapter['figure_label'],page=number,description=chapter['title'],legend=chapter['caption'],
            paragraphs=chapter['paragraphs'],generated_at=stamp,producer='finished',
            figure_sha256={str(p.relative_to(FIG)):sha(p) for p in artifacts},
            source_links=chapter['sources'],renderer_manifest=str(renderer_manifest.relative_to(FIG)),
            renderer_manifest_sha256=sha(renderer_manifest),
            reproduction=original.get('reproduction',original.get('reproduce',{})),
            report_reproduction='/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_mechanisms_report_20261008.py',
            corrected_deBelly_revision=number>=8)
        (FIG/(stem.name+'.manifest.json')).write_text(json.dumps(attached,indent=2)+'\n')
    pdf=pdf_report(stamp);html_report(stamp,source_times)
    text=[TITLE,'Generated '+stamp+' · producers finished','']
    for chapter in CHAPTERS:
        text.extend(['Figure '+chapter['figure_label']+' — '+chapter['title'],'',
                     '\n\n'.join(chapter['paragraphs']),'',chapter['caption'],''])
    (FIG/'report_text.txt').write_text('\n'.join(text))
    audit=audit_report(pdf);compatibility()
    (DATA/'quality_report_audit.json').write_text(json.dumps(dict(generated_at=stamp,**audit),indent=2)+'\n')
    record=dict(generated_at=stamp,producer='finished',pages=len(CHAPTERS),
        renderer_sha256=sha(Path(__file__)),source_data_times={k:source_times[k] for k in ['comparison','initial']},
        source_hashes=source_times['verified_hashes'],active_manifests=source_times['active_manifests'],
        figure_hashes={c['image']:sha(FIG/c['image']) for c in CHAPTERS},
        html_sha256=sha(FIG/'index.html'),pdf_sha256=sha(pdf),text_sha256=sha(FIG/'report_text.txt'),
        revision='deBelly resting mechanics, tension feedback ablation, strength-selected reversal histories and attachment redistribution; prior reader archived; obsolete turnover chapter removed',
        checks=dict(corrected_initial_checkpoint_hashes=True,all_85_other_model_rows_unchanged=True,
                    active_figure_and_source_hashes=True,corrected_debelly_only=True,**audit),
        chapters=[dict(page=i,number=c['figure_label'],id=c['id'],title=c['title'],image=c['image']) for i,c in enumerate(CHAPTERS,1)])
    (DATA/'report_provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(html=str(FIG/'index.html'),pdf=str(pdf),pages=len(CHAPTERS),generated_at=stamp,audit=audit),indent=2))


if __name__=='__main__':
    main()
