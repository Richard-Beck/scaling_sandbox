"""Compact portable S11 reader; full native CSV/NPZ remain explicit HPC sources."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import html
import json
import shutil

NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring_ring='Spring')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path);args=p.parse_args();run=args.run.resolve()
    source=run/'results/export';target=args.output.resolve() if args.output else run/'results/portable_export'
    if target.exists():raise FileExistsError('Portable atlas packaging uses a fresh immutable directory')
    atlas=read(source/'complete.json');csv=read(source/'page_data/complete.json');pdf=read(source/'pdf/lossless_pdf_audit.json')
    independent=read(run/'results/independent_native_atlas_audit.json');assert atlas['passed'] and csv['passed'] and pdf['passed'] and independent['passed']
    assert csv['source_atlas_complete_sha256']==sha(source/'complete.json')==independent['atlas_complete_sha256']
    for name in ['webp','pdf','page_data']:shutil.copytree(source/name,target/name)
    metadata=target/'metadata';metadata.mkdir()
    files={
        source/'native_atlas_index.json':'native_source_index.json',source/'candidate_features_and_selection.csv':'all_candidate_features_and_selection.csv',
        source/'selection_rules.json':'selection_rules.json',source/'frozen_atlas.py':'frozen_atlas.py',source/'complete.json':'original_atlas_complete.json',
        run/'manifest.json':'source_manifest.json',run/'results/independent_native_atlas_audit.json':'independent_native_atlas_audit.json',
        run/'results/independent_numeric_diagnostic.json':'independent_numeric_diagnostic.json',run/'results/failed_numeric_audit_receipt.json':'failed_numeric_audit_receipt.json'}
    optional={'independent_numeric_diagnostic.json','failed_numeric_audit_receipt.json'}
    for path,name in files.items():
        if path.exists():shutil.copy2(path,metadata/name)
        elif name not in optional:raise FileNotFoundError(path)
    for name in ['physical_identity_binding.json','canonical_case_inventory.json','logical_coverage_freeze.json','identity_binding_freeze.json',
                 'resource_pilot_gate.json','production_resource_receipt.json','pipeline_freeze.json']:
        path=run/'inputs'/name
        if path.exists():shutil.copy2(path,metadata/name)
    if (run/'results/caption_layout_crosswalk.json').exists():
        shutil.copy2(run/'results/caption_layout_crosswalk.json',metadata/'caption_layout_crosswalk.json')
    if (run/'results/text_layout_audit.json').exists():
        shutil.copy2(run/'results/text_layout_audit.json',metadata/'text_layout_audit.json')
    index=read(source/'native_atlas_index.json');pages=read(source/'page_data/page_data_index.json')['pages'];entries=[]
    for group,page in zip(index['groups'],pages):
        assert group['condition']==page['condition']
        c=group['condition'];title=f"{c['experiment']} · {NAMES.get(c['model'],c['model'])} · size {c['C_label']:g} · f={c['f']:g} · {c['kind']} · σ={c['sigma']:g} · horizon={c['native_horizon']:g} · cadence={c['profile_cadence']:g}"
        files=[dict(label='PDF volume, page '+str(page['pdf_page']),url='pdf/'+page['pdf_volume']+'#page='+str(page['pdf_page'])),
            dict(label='Per-page selection and native-source identities CSV',url='page_data/'+Path(page['selection_CSV']['path']).name)]
        for name,item in page['companions'].items():files.append(dict(label=name.replace('_',' ')+' · full plotted values CSV.gz',url='page_data/'+Path(item['path']).name))
        archives=[dict(path=item['path'],sha256=item['sha256'],bytes=item['bytes'],rows=item['rows']) for item in page['native_field_CSV_shards']]
        entries.append(dict(title=title,image='webp/'+Path(group['artifact']['webp']).name,files=files,archives=archives,cells=group['artifact']['displayed'][0]['native_cells']))
    css='body{font:17px/1.5 system-ui;max-width:1400px;margin:32px auto;padding:0 22px;color:#14232b}img{width:100%;border:1px solid #ddd}input,select{font:16px system-ui;margin:8px;padding:8px;max-width:100%}code{font:12px monospace;overflow-wrap:anywhere}.note{background:#eef4f4;padding:15px;border-left:4px solid #46817a}li{margin:6px 0}'
    document='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure4 S11 E0 native atlas</title><style>'+css+'</style><h1>S11 · E0 native activity atlas</h1><p>Generated '+html.escape(csv['generated_utc'])+'. The E0 atlas is complete; this page has no timed refresh. '+str(atlas['groups'])+' condition pages show '+str(atlas['selected_displays'])+' seed histories selected from '+str(atlas['independent_native_cases'])+' independent native cases. Every model, available size, input strength, control and forcing arm remains explicit.</p><p>Each condition shows the predeclared summary-feature medoid and a deterministic distinct random case; singleton controls show one case. The seven-feature distance uses the condition IQR, a nonzero-range fallback and stable case-ID ties. The first half is shaded and excluded from selection statistics. No main-figure example is selected here.</p><p>All saved native active fields use physical arc. Spring uses its actual circumference at every saved time; only mesh boundary geometry between saved times is interpolated. Four snapshots use native samples nearest 0, 0.5, 0.75 and 1 horizon, with actual times printed. |P| uses basal M₁ normalization. Source concentration values and all retained frames are preserved.</p><div class="note"><strong>Native source archives:</strong> Every full-field per-page CSV.gz exists in the HPC workspace and passes 17-significant-digit float64 roundtrip checks, including time and actual circumference. These large CSVs, original NPZs and original PNGs are retained as local/HPC source archives outside this compact package. Their exact paths, hashes, row counts and grids are recorded in the tracked <a href="page_data/native_field_archive_inventory.csv">native-field archive inventory</a> and per-page selection CSVs. Paths shown below are provenance text, not links to absent GitHub files.</div><p><a href="pdf/index.html">Five lossless PDF volumes</a> · <a href="pdf/page_selection_inventory.csv">PDF page inventory</a> · <a href="metadata/all_candidate_features_and_selection.csv">All 962 candidate features and selection scores</a> · <a href="metadata/selection_rules.json">Predeclared selection rules</a> · <a href="metadata/independent_native_atlas_audit.json">Independent selection/geometry/pixel audit</a> · <a href="pdf/lossless_pdf_audit.json">Lossless PDF proof</a> · <a href="page_data/page_data_index.json">Complete numeric archive/companion index</a></p><input id="search" placeholder="Model, size, strength, control, stage or horizon"><select id="page"></select><h2 id="title"></h2><img id="image" alt="Native physical-arc kymographs, saved profiles and scalar traces"><div id="files"></div><details><summary>Full-field CSV source archive paths (HPC-only)</summary><p>Wide CSV columns: case index, native frame, native time, actual circumference, then every native active cell. The q grid is in the linked companion CSV.</p><div id="archives"></div></details><script>const pages='+json.dumps(entries,separators=(',',':')).replace('</',r'<\/')+';const byId=id=>document.getElementById(id);let selected=[];function filter(){const query=byId("search").value.toLowerCase();selected=pages.filter(p=>p.title.toLowerCase().includes(query));byId("page").replaceChildren(...selected.map((p,i)=>{const o=document.createElement("option");o.value=i;o.textContent=p.title;return o}));show()}function show(){const p=selected[byId("page").value];if(!p)return;byId("title").textContent=p.title;byId("image").src=p.image;const list=document.createElement("ul");for(const f of p.files){const li=document.createElement("li"),a=document.createElement("a");a.href=f.url;a.textContent=f.label;li.append(a);list.append(li)}byId("files").replaceChildren(list);const archives=document.createElement("ul");for(const item of p.archives){const li=document.createElement("li"),code=document.createElement("code");code.textContent=item.path+" · SHA256 "+item.sha256+" · "+item.rows+" native rows × "+p.cells+" cells";li.append(code);archives.append(li)}byId("archives").replaceChildren(archives)}byId("search").addEventListener("input",filter);byId("page").addEventListener("change",show);filter();</script></html>'
    context=atlas['stage']
    document=document.replace('S11 E0 native atlas','S11 '+html.escape(context)+' native atlas').replace('S11 · E0 native activity atlas','S11 · '+html.escape(context)+' native activity atlas').replace('The E0 atlas is complete','The '+html.escape(context)+' atlas is complete')
    document=document.replace('Five lossless PDF volumes',str(len(pdf['volumes']))+' lossless PDF volumes').replace('All 962 candidate features','All '+str(atlas['independent_native_cases'])+' candidate features')
    document=document.replace('seed histories selected from '+str(atlas['independent_native_cases'])+' independent native cases',
        'native histories selected from '+str(atlas['independent_native_cases'])+' physical runs ('+str(atlas.get('stochastic_realizations','unavailable'))+' noisy realizations; '+str(atlas.get('deterministic_controls','unavailable'))+' deterministic controls)')
    document=document.replace('if(!p)return;','if(!p){byId("title").textContent="No matching condition";byId("image").removeAttribute("src");byId("files").replaceChildren();byId("archives").replaceChildren();return;}')
    (target/'index.html').write_text(document)
    shutil.copy2(__file__,metadata/'frozen_package.py')
    content=dict(purpose='Native-field supplement with predeclared medoid and distinct random seed selection; no main-example choice.',analysis_arm=context,
        stochastic_realizations=atlas.get('stochastic_realizations'),deterministic_controls=atlas.get('deterministic_controls'),
        groups=atlas['groups'],displayed_histories=atlas['selected_displays'],independent_candidate_cases=atlas['independent_native_cases'],
        candidate_count_scope='Physical runs; deterministic controls contribute zero stochastic independent seeds. Cross-sigma pairing and E4 retained-prefix dependence remain explicit in the final coverage master.',
        PDF_volumes=pdf['volumes'],source_fields='Exact native active concentrations, active Rac for Holmes3; no source resampling or quantization.',
        geometry='Physical arc from fixed material origin; actual Spring circumference at every saved time. Only boundary mesh geometry is interpolated.',
        snapshots='Nearest native saved samples to 0, 0.5, 0.75 and 1 of the horizon; actual times printed and indexed.',
        selection='Last-half seven-feature L1 medoid with IQR/range scaling, stable ties, plus deterministic distinct SHA256 random seed. Singleton controls have one case.',
        numeric_data='All field/time/actualC CSV.gz files exist and roundtrip bit-exactly in HPC source archives; compact package carries full snapshots/scalar curves/q grids and archive/selection inventories.',
        source_atlas_complete_sha256=sha(source/'complete.json'),source_native_CSV_complete_sha256=sha(source/'page_data/complete.json'),
        independent_native_audit_sha256=sha(run/'results/independent_native_atlas_audit.json'),main_example_selection=False,native_integrations=0)
    write(metadata/'content.json',content)
    assets={str(path.relative_to(target)):sha(path) for path in sorted(target.rglob('*')) if path.is_file()}
    assert not any('native_field_CSV_archives/' in path or path.endswith('.png') or path.endswith('.npz') for path in assets)
    assert all((target/file['url'].split('#')[0]).is_file() for entry in entries for file in entry['files'])
    assert all((target/entry['image']).is_file() for entry in entries)
    receipt=dict(passed=True,created_utc=datetime.now(timezone.utc).isoformat(),analysis_arm=context,groups=len(entries),
        all_portable_links_exist=True,no_HPC_archive_href=True,fullfield_CSVs_distribution='HPC-only; tracked inventory plus source selection hashes',
        original_PNG_and_NPZ_retained=True,maximum_portable_file_bytes=max(path.stat().st_size for path in target.rglob('*') if path.is_file()),
        portable_logical_bytes=sum(path.stat().st_size for path in target.rglob('*') if path.is_file()),source_packager_sha256=sha(__file__),assets_sha256=assets)
    assert receipt['maximum_portable_file_bytes']<100*1024**2
    write(target/'package_audit.json',receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k!='assets_sha256'}))


if __name__=='__main__':main()
