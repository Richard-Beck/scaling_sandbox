"""Bind unchanged exact numeric archives to a caption-only atlas revision."""
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
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--revised',type=Path,required=True);args=p.parse_args()
    original,revised=args.original.resolve(),args.revised.resolve();old=read(original/'manifest.json');new=read(revised/'manifest.json')
    assert old['jobs']==new['jobs'] and old['rules_sha256']==new['rules_sha256']
    assert old['inventory_sha256']==new['inventory_sha256'] and old['verification_sha256']==new['verification_sha256']
    assert read(original/'results/export/page_data/complete.json')['float64_roundtrip_equal_every_field_time_actualC']
    out=revised/'results/export/page_data';out.mkdir(exist_ok=False);workers=revised/'results/csv_workers';workers.mkdir(exist_ok=False)
    shutil.copy2(original/'inputs/frozen_csv_export.py',revised/'inputs/frozen_csv_export.py')
    crosswalk=[]
    for job in old['jobs']:
        slot=job['slot'];source=original/'workers'/f'{slot:04d}.json';destination=revised/'workers'/f'{slot:04d}.json'
        a,b=read(source),read(destination)
        for key in ['cases','features','feature_names','selection','condition']:
            assert a[key]==b[key],f'Caption revision changed scientific {key}: slot{slot}'
        assert a['artifact']['displayed']==b['artifact']['displayed']
        if len(a['artifact']['displayed'])==2:
            assert a['artifact']['png_sha256']==b['artifact']['png_sha256'] and a['artifact']['webp_sha256']==b['artifact']['webp_sha256']
        old_csv=original/'results/csv_workers'/f'{slot:04d}.json';record=read(old_csv)
        assert record['source_atlas_PNG_sha256']==a['artifact']['png_sha256'] and record['selected_cases']==[x['case_id'] for x in b['artifact']['displayed']]
        for item in list(record['companions'].values())+[record['selection_CSV']]:
            path=Path(item['path']);assert sha(path)==item['sha256'];target=out/path.name;shutil.copy2(path,target)
            item['original_numeric_export_path']=str(path);item['path']=str(target)
        # Full-field source files remain the exact original verified archive.
        # The new receipt explicitly preserves its lineage and only rebinds
        # the figure-image identity; no numeric source is rewritten.
        relation=dict(slot=slot,original_numeric_receipt_sha256=sha(old_csv),
            original_PNG_sha256=a['artifact']['png_sha256'],revised_PNG_sha256=b['artifact']['png_sha256'],
            source_native_records_identical=True,features_and_selected_IDs_exact=True,
            snapshots_and_actual_geometry_identical=True,two_column_pixels_identical=len(a['artifact']['displayed'])==2,
            original_run=str(original),original_source_manifest_sha256=sha(original/'manifest.json'))
        record['source_atlas_PNG_sha256']=b['artifact']['png_sha256'];record['caption_only_rebind']=relation
        write(workers/f'{slot:04d}.json',record);crosswalk.append(relation)
    for name in ['independent_numeric_diagnostic.json','failed_numeric_audit_receipt.json']:
        shutil.copy2(original/'results'/name,revised/'results'/name)
    receipt=dict(passed=True,created_utc=datetime.now(timezone.utc).isoformat(),original_run=str(original),revised_run=str(revised),
        groups=len(crosswalk),scientific_source_and_selection_unchanged=True,native_CSV_bytes_unchanged=True,
        changes='Wrap singleton-control footer text and correct singular case wording. Every native input, feature, selected ID, snapshot and physical geometry is identical; all two-column PNG/WebP bytes are identical.',
        original_native_CSV_complete_sha256=sha(original/'results/export/page_data/complete.json'),
        original_independent_native_audit_sha256=sha(original/'results/independent_native_atlas_audit.json'),
        original_source_manifest_sha256=sha(original/'manifest.json'),revised_source_manifest_sha256=sha(revised/'manifest.json'),
        native_integrations=0,new_scientific_selections=0,source_sha256=sha(__file__),crosswalk=crosswalk)
    write(revised/'results/caption_layout_crosswalk.json',receipt)
    shutil.copy2(__file__,revised/'inputs/frozen_caption_rebind.py')
    print(json.dumps({k:v for k,v in receipt.items() if k!='crosswalk'}))


if __name__=='__main__':main()
