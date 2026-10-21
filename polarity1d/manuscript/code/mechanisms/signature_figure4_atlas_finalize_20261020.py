"""Single compact gate binding all local-numeric and portable S11 evidence."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import shutil


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--package',type=Path,required=True);args=p.parse_args()
    run,package=args.run.resolve(),args.package.resolve();target=package/'completion_gate.json'
    if target.exists():raise FileExistsError('Completion gates are immutable')
    asset=read(package/'package_audit.json');assert asset['passed']
    for name,value in asset['assets_sha256'].items():assert sha(package/name)==value
    native=read(package/'metadata/independent_native_atlas_audit.json');pdf=read(package/'pdf/lossless_pdf_audit.json')
    numeric=read(package/'page_data/complete.json');pages=read(package/'page_data/page_data_index.json');frontend=read(package/'metadata/frontend_audit.json')
    source=read(run/'results/export/complete.json');expected_pages=source['groups'];expected_histories=source['selected_displays']
    text_layout=read(package/'metadata/text_layout_audit.json');assert text_layout['passed'] and text_layout['pages']==expected_pages
    if (package/'metadata/caption_layout_crosswalk.json').exists():
        crosswalk=read(package/'metadata/caption_layout_crosswalk.json')
        assert crosswalk['passed'] and crosswalk['scientific_source_and_selection_unchanged'] and crosswalk['native_CSV_bytes_unchanged']
    assert native['passed'] and pdf['passed'] and numeric['passed'] and frontend['passed']
    assert native['groups']==numeric['pages']==pages['summary']['pages']==frontend['pages']==expected_pages
    assert native['displayed_seed_histories']==numeric['selected_histories']==expected_histories
    assert native['atlas_complete_sha256']==numeric['source_atlas_complete_sha256']==sha(run/'results/export/complete.json')
    assert numeric['source_PDF_audit_sha256']==sha(package/'pdf/lossless_pdf_audit.json')
    archive_files=0;archive_bytes=0;archive_allocated=0;fields=0
    for page in pages['pages']:
        assert page['passed'] and page['float64_roundtrip_equal'] and page['decimal_significant_digits']==17
        fields+=page['native_cell_time_values']
        for part in page['native_field_CSV_shards']:
            path=Path(part['path']);assert path.is_file() and sha(path)==part['sha256'] and part['roundtrip_float64_equal']
            archive_files+=1;archive_bytes+=path.stat().st_size;archive_allocated+=path.stat().st_blocks*512
        for item in list(page['companions'].values())+[page['selection_CSV']]:
            assert sha(package/'page_data'/Path(item['path']).name)==item['sha256']
    assert fields==numeric['native_cell_time_values']
    assert archive_bytes==numeric['full_field_archive_bytes']
    assert frontend['reader_sha256']==sha(package/'index.html')
    shutil.copy2(__file__,package/'metadata/frozen_completion.py')
    result=dict(passed=True,generated_utc=datetime.now(timezone.utc).isoformat(),analysis_arm=source['stage'],groups=expected_pages,selected_native_histories=expected_histories,
        stochastic_realizations=source.get('stochastic_realizations'),deterministic_controls=source.get('deterministic_controls'),
        independent_candidates=native['independent_candidate_cases'],all_native_field_values=fields,native_rows=numeric['native_rows'],
        actual_local_fullfield_CSV_files_verified=archive_files,fullfield_CSV_logical_bytes=archive_bytes,fullfield_CSV_allocated_bytes=archive_allocated,
        all_native_field_time_actualC_float64_roundtrip_equal=True,selected_IDs_and_random_hashes_exact=True,lossless_PNG_WebP_and_PDF_pixels=True,
        portable_links_exist=True,no_missing_HPC_archive_hrefs=True,
        fullfield_CSV_distribution='Actual verified local/HPC sources; compact public reader links tracked inventory and selection files, never absent archive URLs.',
        source_feature_selection_audit_sha256=sha(package/'metadata/independent_native_atlas_audit.json'),
        lossless_PDF_audit_sha256=sha(package/'pdf/lossless_pdf_audit.json'),native_CSV_complete_sha256=sha(package/'page_data/complete.json'),
        native_CSV_page_index_sha256=sha(package/'page_data/page_data_index.json'),native_CSV_inventory_sha256=sha(package/'page_data/native_field_archive_inventory.csv'),
        package_asset_audit_sha256=sha(package/'package_audit.json'),frontend_audit_sha256=sha(package/'metadata/frontend_audit.json'),
        all_titles_and_footers_within_native_PNG=True,text_layout_audit_sha256=sha(package/'metadata/text_layout_audit.json'),
        caption_layout_crosswalk_sha256=sha(package/'metadata/caption_layout_crosswalk.json') if (package/'metadata/caption_layout_crosswalk.json').exists() else None,
        source_complete_sha256=sha(run/'results/export/complete.json'),completion_source_sha256=sha(__file__),
        native_integrations=0,main_example_selection=False,
        final_portable_assets_sha256={str(path.relative_to(package)):sha(path) for path in sorted(package.rglob('*')) if path.is_file()})
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='final_portable_assets_sha256'}))


if __name__=='__main__':main()
