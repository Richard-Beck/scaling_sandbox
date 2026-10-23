"""Verify the single-C landscape layout and exact reuse of audited evidence."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
DEV=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022'
FIG=DEV/'figure3_single_C'
PREVIOUS=DEV/'figure3'
DEST=ROOT/'manuscript_outline/revision_round3_figure3_single_C_20261022'


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    layout=read(FIG/'supporting/metadata/layout.json')
    assert layout['profile_sizes']==[30.]
    assert layout['profile_labels']==['i','ii']
    assert layout['native_integrations']==0
    a,b,c=(layout[key] for key in ('profile_axes','time_axes','speed_axes'))
    assert [len(a),len(b),len(c)]==[8,4,4]
    assert max(x+w for x,y,w,h in a)<min(x for x,y,w,h in b+c)
    assert min(y for x,y,w,h in b)>max(y+h for x,y,w,h in c)
    assert abs(min(y for x,y,w,h in a)-min(y for x,y,w,h in c))<1e-10
    assert abs(max(y+h for x,y,w,h in a)-max(y+h for x,y,w,h in b))<1e-10
    width,height=layout['figure_size_inches']
    assert all(w*width>1.6 and h*height>2 for x,y,w,h in b+c)
    assert all(w*width>1.8 and h*height>1 for x,y,w,h in a)
    reuse=read(FIG/'supporting/metadata/layout_source_preservation.json')
    for relative,digest in reuse['source_data_sha256'].items():
        assert sha(PREVIOUS/relative)==digest
        assert sha(FIG/relative)==digest
    assert read(FIG/'supporting/metadata/render_coordinates.json')==read(PREVIOUS/'supporting/metadata/render_coordinates.json')
    upstream=read(PREVIOUS/'supporting/metadata/verification.json')
    assert upstream['passed'] and upstream['observed_crossings']==120
    assert upstream['all_four_marker_topologies_correct']
    assert sha(ROOT/reuse['previous_main'])==reuse['previous_main_sha256']
    assert read(FIG/'supporting/metadata/build_receipt.json')['layout_revision_native_integrations']==0
    legend=read(FIG/'legend.json')
    assert 'at C = 30' in legend['legend'][0]
    assert '47.6' not in legend['legend'][0]
    assert 'i–ii markers' in legend['legend'][1]
    assert 'iii:' not in legend['legend'][3] and 'iv:' not in legend['legend'][3]
    main_pdf=PdfReader(FIG/'figure3_round3.pdf')
    assert len(main_pdf.pages)==1
    text=main_pdf.pages[0].extract_text()
    assert '(i)' in text and '(ii)' in text
    assert '(iii)' not in text and '(iv)' not in text
    assert 'Cue placement relative' not in text
    # Every retained page is identical after concatenation, including the
    # preceding round-3 main/legend with both enlarged-ring examples.
    preserved=read(FIG/'supporting/metadata/supplement_preservation.json')
    supplement=PdfReader(FIG/'figure3_round3_supplement.pdf')
    offset=0
    for source in preserved['sources']:
        path=ROOT/source['path'];assert sha(path)==source['sha256']
        for page in PdfReader(path).pages:
            assert page.get_contents().get_data()==supplement.pages[offset].get_contents().get_data()
            offset+=1
    assert offset==len(supplement.pages)==17
    publication=read(DEST/'figure3_with_legend.build_manifest.json')
    assert sha(ROOT/publication['pdf']['path'])==publication['pdf']['sha256']
    pages=PdfReader(ROOT/publication['pdf']['path']).pages
    assert len(pages)==2
    assert pages[0].get_contents().get_data()==main_pdf.pages[0].get_contents().get_data()
    prior=read(ROOT/'manuscript_outline/outline_v7.json')
    current=read(ROOT/'manuscript_outline/outline_v8.json')
    integration=read(DEST/'integration.json')
    assert sha(ROOT/integration['baseline'])==integration['baseline_sha256']
    for ident,record in integration['protected_sections'].items():
        assert next(x for x in prior['sections'] if x['id']==ident)==next(x for x in current['sections'] if x['id']==ident)
        for relative,digest in record['asset_sha256'].items():assert sha(ROOT/relative)==digest
    assert read(ROOT/'manuscript_outline/outline_v8.validation.json')['passed']
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,
        profile_sizes=[30.],displayed_profiles=8,highlighted_placements=['i','ii'],
        A_left_of_B_and_C=True,B_above_C=True,all_three_curve_sizes_preserved=True,
        all_source_data_identical=True,native_integrations=0,
        upstream_scientific_audit=upstream['generated_at'],
        complete_supplement_pages=offset,previous_main_and_legend_preserved=True,
        all_other_sections_preserved=list(integration['protected_sections']))
    for path in (DEST/'verification.json',FIG/'supporting/metadata/verification.json'):
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
