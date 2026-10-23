#!/home/4473331/.conda/envs/scaling_sandbox/bin/python
"""Apply the Claude follow-up presentation edits using terminal cached evidence.

The executed analysis and its receipts remain immutable. This entry point changes
prose, supplement labels/page order and a presentation audit, without estimating
new quantities, writing numerical exports or running a native model.
"""
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json,html,sys
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parents[1]
ROOT=next(p for p in OUT.parents if (p/'polarity1d/manuscript/code/mechanisms').is_dir())
META=OUT/'supporting/metadata'
DATA=OUT/'supporting/data'
ARCHIVE=OUT/'supporting/reviewed_terminal_028d9f5'
REVIEW_HASHES={
 'reviews/claude/20261023_manuscript_round6.md':'6dd620c7cd406625a473c248f80da028c6de15f2143c916fadf7e8325c17d041',
 'reviews/claude/20261023_manuscript_round6_2.md':'69e8c13963c31d912c1645c3fb85efa18c35d09db7aae96a43a9c05db7b4f5e2'}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text())
def dump(path,value):
    temporary=Path(str(path)+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    temporary.replace(path)
def verify_preserved():
    expected=load(ARCHIVE/'supporting/metadata/data_export_hashes.json')
    assert len(expected)==62
    current={name:sha(DATA/name) for name in expected}
    assert current==expected,'A numerical export differs from the reviewed terminal state'
    assert sha(META/'data_export_hashes.json')==sha(ARCHIVE/'supporting/metadata/data_export_hashes.json')
    assert sha(META/'stable_artifact_handshake.json')==sha(ARCHIVE/'supporting/metadata/stable_artifact_handshake.json')
    assert sha(META/'final_analysis_source_freeze.json')==sha(ARCHIVE/'supporting/metadata/final_analysis_source_freeze.json')
    for rel,digest in REVIEW_HASHES.items():assert sha(ROOT/rel)==digest,rel
    for rel,digest in load(ARCHIVE/'FILE_HASHES.json').items():assert sha(ARCHIVE/rel)==digest,rel
    freeze=load(META/'final_analysis_source_freeze.json')
    for record in freeze['files']:
        assert sha(META/'final_analysis_source_snapshot'/record['snapshot'])==record['sha256'],record['snapshot']
        if record['snapshot']!='rebuild_figure4.py':assert sha(ROOT/record['source'])==record['sha256'],record['source']
    receipt=load(ARCHIVE/'supporting/metadata/stable_artifact_handshake.json')
    for rel in ['figure4_main.pdf','figure4_main.png']:
        assert sha(OUT/rel)==receipt['artifact_sha256'][rel]
    assert sha(OUT/'supporting/vectors/figure4_A_J.pdf')==receipt['artifact_sha256']['figure4_main.pdf']
    return current

def write_documents():
    verify_preserved()
    section=load(ARCHIVE/'legend.json')
    section['role']='Size-dependent heading persistence and limits of region-based diagnostic accounting.'
    section['legend']=[s for s in section['legend'] if not s.startswith('The strong-input regime')]
    for index,text in enumerate(section['legend']):
        if text.startswith('I,'):
            section['legend'][index]=('I, Signed remainder share versus input, using point regions with native-amplitude scaling. '
              'R includes amplitude, within-region shape, unresolved activity and background contributions; its sign and magnitude depend on the accounting convention. '
              'This convention retains measured native coherence, avoids the unit-modulus single-region normalization artefact and gives a common comparison across inputs and sizes. '
              'The decomposition is descriptive. Exact own-mask and normalization sensitivities appear in S4.5; complete-mask reconstruction is an algebraic identity, not causal validation.')
        elif text.startswith('J,'):
            section['legend'][index]=('J, Probability that a certified persistent secondary birth becomes globally dominant within 20 input correlation times, at each model’s labelled input. '
              'Filled circles: 0 < d ≤ C/6; open squares: d > C/6, where d is estimated edge-to-edge birth–incumbent separation. '
              'Probabilities concern the certified positive-distance subset; a [0,1] interval is uninformative. '
              'Lower finite-horizon incidence can reflect slower or rarer takeover; lifetime takeover probability is unidentified. '
              'Input selection, certification denominators and interval procedures are in Methods.')
    section['narrative']=[s for s in section['narrative'] if not s.startswith('Fourier rate-collapse')]
    section['measurements'][1]=('Track the physical lateral displacement of both front boundaries. '
       'Size-suppressed physical drift would favour Wave-pinning-like or Goryachev-like susceptibility; this assay does not separate them. '
       'The Fig. 2 first-rise assay separates these alternatives, because initiation slows with size only in Goryachev. '
       'Nearly size-independent physical drift with slower angular wandering would favour Holmes/Spring-like geometry but does not separate those models. '
       'The Goryachev quantitative drift bridge remains unreliable.')
    section['measurements'][4]=('Measure perimeter Fourier spectra and mode correlations against physical wavenumber in both ploidy groups. '
       'A common rate curve would support a detector-free description; disagreement would reject that description without requiring a domain detector. '
       'This assay tests that description rather than distinguishing mechanisms of the size effect (S4.6–S4.8).')
    section['limitations']=[s for s in section['limitations'] if s!='A universal strong-input mechanism is not identified.']
    section['legend'].append('Scope/limitations: The strong-input regime remains unexplained.')
    dump(OUT/'legend.json',section)
    (OUT/'legend.md').write_text(section['title']+'\n\n'+'\n\n'.join(section['legend'])+'\n')

    methods=(ARCHIVE/'methods_provenance.md').read_text()
    methods=methods.replace('T is inherited eligible dominant-tenure exposure','T is eligible dominant-tenure exposure')
    methods=methods.replace('Main I keeps the original round 5 native-amplitude convention; these changes of convention are a sensitivity, not silent redefinition.',
      'Main I uses the original point-region/native-amplitude convention because it retains measured native coherence, avoids the unit-modulus single-region normalization artefact and compares all inputs and sizes consistently. The sign and magnitude of R depend on the accounting convention. Alternative mask and normalization conventions are sensitivities, not causal decompositions.')
    methods=methods.replace('Intervals resample whole independent histories including zero-event histories. These probabilities concern the certified subset over finite follow-up.',
      'The selected inputs are f = 0.05 for Wave-pinning, Holmes and Spring, and f = 0.025 for Goryachev; minimum total certified births across the five sizes are 590, 231, 580 and 1064, respectively. Total certified births and usable positive-distance denominators are exported separately in J_takeover_selection.csv and J_takeover_points.csv. Intervals resample whole independent histories, including zero-event histories, 2000 times. If fewer than 97.5% of resamples have a defined probability, the interval is replaced by the uninformative probability bounds [0,1]. These probabilities concern the certified subset over finite follow-up.')
    methods=methods.replace('K–L enter the main figure only when the declared completed-data screens support the hypothesis.',
      'The completed rate and slope diagnostics are supplemental S4.7 and S4.8. Their fixed four-model promotion criterion is not met; the primary and stencil-sensitivity results remain in the scientific readme and data tables. Historical K_/L_ CSV filenames are retained to preserve the executed numerical exports.')
    methods+='''\n## Presentation after the completed analysis\n\nThe Claude follow-up review at commit 028d9f5 led to presentation changes only: a shorter J caption, explicit convention dependence and rationale for I, corrected assay discrimination, and supplemental labels S4.7–S4.11. The main A–J artwork is byte-identical to that reviewed commit. The supplement reads previously exported points and guides without pooling, fitting, bootstrapping or writing numerical files. All 62 numerical data exports and both Claude reviews remain byte-identical; supporting/metadata/postreview_presentation_audit.json records the per-file proof and prior/current artifact hashes. The original stable_artifact_handshake.json, final_analysis_source_freeze.json, final_analysis_source_snapshot and native receipts are preserved. The current rebuild entry point and two presentation scripts were written after the executed analysis; their hashes document presentation sources, not a claim that they produced the native results. The reviewed terminal artifacts are archived in supporting/reviewed_terminal_028d9f5.\n'''
    (OUT/'methods_provenance.md').write_text(methods)

    readme=(ARCHIVE/'scientific_readme.md').read_text()
    readme=readme.replace('supplement S4.7b','supplement S4.10')
    readme=readme.replace('## Current targeted campaign','## Completed targeted campaign')
    readme=readme.replace('The independent precision controller updates mode evidence and rebuilds artifacts after the authorized additions. The source round 5 deliverables are unchanged.',
      'All authorized additions and precision analyses are complete; the native controller and live producer have stopped. The source round 5 deliverables are unchanged.')
    readme=readme.replace('These descriptive screens are not significance tests, and the table remains provisional while allocated cohorts are incomplete.',
      'These completed-cohort descriptive screens are not significance tests. The fixed primary rule and its stencil sensitivities are shown below; the corresponding plots are S4.7–S4.8.')
    start=readme.index('This one entry point reads source data, collects completed analyses, exports every plotted point and rebuilds')
    stop=readme.index('\n\nBiological scope:',start)
    readme=readme[:start]+('This entry point performs a presentation rebuild from the completed cached evidence. It verifies all 62 numerical export hashes, preserves the main PDF/PNG, renders the supplement from exported points and guides, updates prose and writes a separate presentation audit. It performs no native integrations, estimates, pooling or bootstrap draws. The executed analysis source and original completion receipt remain archived in supporting/metadata/final_analysis_source_snapshot and stable_artifact_handshake.json.\n\n'
       'Supplement contents: S4.1–S4.5 show certified births, heading curves, region shares, finite-horizon takeovers and mask sensitivities; S4.6 shows all-input spectra; S4.7–S4.8 show mode rates and the slope comparison; S4.9–S4.10 show phase diagnostics; S4.11 shows normalized complex correlations. Historical K_/L_ data filenames remain unchanged. The main caption and conclusions concern A–J; the strong-input limitation is stated in Scope/limitations.')+readme[stop:]
    (OUT/'scientific_readme.md').write_text(readme)
    timestamp=datetime.now(timezone.utc).isoformat()
    original_html=(ARCHIVE/'index.html').read_text()
    begin=original_html.index('<!-- LIVE BEGIN -->');end=original_html.index('<!-- LIVE END -->')+len('<!-- LIVE END -->')
    live=original_html[begin:end]
    (OUT/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Figure 4: enlargement and polarity heading</title></head><body>'
      '<h1>Enlargement and polarity heading</h1><p>Numerical data generation: 2026-10-23T22:58:48.104839+00:00. Presentation generation: '+html.escape(timestamp)+'.</p>'
      '<p><a href="figure4_main.pdf">Main A–J PDF</a> · <a href="figure4_supplement.pdf">Supplement PDF</a> · <a href="scientific_readme.md">Scientific readme</a> · <a href="methods_provenance.md">Methods</a> · <a href="legend.md">Legend</a></p>'
      '<p>Detector-free mode diagnostics: <a href="figure4_supplement.pdf#page=28">S4.7, mode rates</a> and <a href="figure4_supplement.pdf#page=29">S4.8, local slope</a>. Numerical results and precision limitations are in the scientific readme.</p>'
      '<p>Scope/limitations: The strong-input regime remains unexplained.</p><img src="figure4_main.png" style="max-width:100%" alt="Figure 4 A through J">'+live+'</body></html>\n')
    with (ARCHIVE/'supporting/metadata/panel_remap.csv').open(newline='') as handle:rows=list(csv.DictReader(handle))
    for row in rows:
        if row['quantity']=='Mode rate versus wavenumber':row['current_round6']='supplement S4.7'
        if row['quantity']=='Size exponent versus local mode slope':row['current_round6']='supplement S4.8'
    with (META/'panel_remap.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['previous_round6','current_round6','quantity']);writer.writeheader();writer.writerows(rows)
    return timestamp

def write_audit(started):
    numeric=verify_preserved()
    from pypdf import PdfReader
    supplement=PdfReader(OUT/'figure4_supplement.pdf')
    text=[page.extract_text() or '' for page in supplement.pages]
    assert len(text)==35
    expected_labels={1:'S4.1',2:'S4.1b',3:'S4.2',4:'S4.3',5:'S4.4',6:'S4.4',7:'S4.5',12:'S4.6',27:'S4.6',28:'S4.7',29:'S4.8',30:'S4.9',31:'S4.10',32:'S4.11',35:'S4.11'}
    for page,label in expected_labels.items():assert label in text[page-1],(page,label)
    receipt=load(META/'stable_artifact_handshake.json')
    current={name:sha(OUT/name) for name in receipt['artifact_sha256']}
    sources=['rebuild_figure4.py','scripts/presentation_postreview.py','scripts/cached_supplement.py']
    archive_hashes=load(ARCHIVE/'FILE_HASHES.json')
    audit={
      'stable':True,'no_more_writers':True,'presentation_only':True,
      'presentation_started_utc':started,'stable_utc':datetime.now(timezone.utc).isoformat(),
      'reviewed_commit':'028d9f5','prior_artifact_sha256':receipt['artifact_sha256'],'artifact_sha256':current,
      'numeric_data_exports_unchanged':True,'numeric_data_export_count':len(numeric),
      'per_file_sha256':numeric,'numeric_data_export_sha256':numeric,
      'data_export_hash_manifest_sha256':sha(META/'data_export_hashes.json'),
      'original_stable_handshake_sha256':sha(META/'stable_artifact_handshake.json'),
      'original_stable_handshake_preserved':True,'final_analysis_source_freeze_preserved':True,
      'final_analysis_source_freeze_sha256':sha(META/'final_analysis_source_freeze.json'),
      'final_analysis_snapshot_files_verified':10,
      'unchanged_executed_analysis_sources':9,
      'source_differences':{
         'rebuild_figure4.py':{'executed_analysis_sha256':sha(META/'final_analysis_source_snapshot/rebuild_figure4.py'),'current_sha256':sha(OUT/'rebuild_figure4.py'),'change':'Presentation-only cached entry point written after executed analysis; original source remains immutable in the final analysis snapshot and reviewed archive.'},
         'scripts/presentation_postreview.py':{'current_sha256':sha(OUT/'scripts/presentation_postreview.py'),'change':'New post-analysis prose, preservation checks and artifact audit; no numerical estimation or data export writes.'},
         'scripts/cached_supplement.py':{'current_sha256':sha(OUT/'scripts/cached_supplement.py'),'change':'New post-analysis renderer reading existing points and guides. Replaces pooling/guide calculation with cached CSV reads and changes supplemental labels/page order.'}},
      'presentation_source_sha256':{name:sha(OUT/name) for name in sources},
      'review_sha256':REVIEW_HASHES,
      'reviewed_archive':str(ARCHIVE.relative_to(ROOT)),
      'reviewed_archive_file_sha256':archive_hashes,
      'reviewed_archive_manifest_sha256':sha(ARCHIVE/'FILE_HASHES.json'),
      'main_artwork_byte_identical':all(current[name]==receipt['artifact_sha256'][name] for name in ['figure4_main.pdf','figure4_main.png']),
      'main_pdf_pages':len(PdfReader(OUT/'figure4_main.pdf').pages),
      'supplement_pdf_pages':len(text),'supplement_label_page_checks':{str(k):v for k,v in expected_labels.items()},
      'supplement_only_mode_panels':['S4.7','S4.8'],
      'changes':['Shortened main J caption; selection and interval procedures moved to Methods.',
                 'Main I explicitly states convention dependence and why point/native-amplitude is displayed.',
                 'Removed supplemental mode-letter references and duplicate Fourier-screen conclusions from main legend/narrative; retained strong-input limitation.',
                 'Reordered spectra before S4.7–S4.8 and numbered later supplemental diagnostics S4.9–S4.11.',
                 'Lateral-edge assay now states Wave-pinning/Goryachev ambiguity and first-rise discrimination; Fourier assay tests the detector-free description.'],
      'no_new_simulations':True,'no_new_numerical_analysis':True,'no_numeric_data_file_writes':True,
      'native_receipts_and_initial_snapshots':'Preserved without modification; no native campaign or producer was restarted.',
      'historical_identifiers':'K_/L_ CSV filenames and identifiers in executed analysis receipts remain unchanged for provenance; current displayed mode panels are S4.7–S4.8.'}
    dump(META/'postreview_presentation_audit.json',audit)
    dump(META/'artifact_hashes.json',{name:current[name] for name in ['figure4_main.pdf','figure4_main.png','figure4_supplement.pdf','legend.json']})
    return audit

def main():
    started=write_documents()
    from cached_supplement import render_cached_supplement
    render_cached_supplement()
    audit=write_audit(started)
    print(json.dumps({'stable':True,'numeric_exports_unchanged':audit['numeric_data_export_count'],'main_artwork_byte_identical':audit['main_artwork_byte_identical'],'supplement_pages':audit['supplement_pdf_pages'],'audit':str(META/'postreview_presentation_audit.json')},indent=2))
if __name__=='__main__':main()
