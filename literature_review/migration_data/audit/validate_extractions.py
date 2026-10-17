"""Validate extraction outputs and append/update their provenance inventory."""
from pathlib import Path
import csv, hashlib, json, datetime, math

BASE=Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with (BASE/'manifest.csv').open(newline='') as f:
    reader=csv.DictReader(f); fields=reader.fieldnames; inventory={r['path']:r for r in reader}
results={}; new_sources={}
hosny=BASE/'hosny_2026/extraction/output_inventory.csv'
if hosny.exists():
    with hosny.open() as f:
        for r in csv.DictReader(f):new_sources[r.get('path','')]=r
for study in ['xuan_2018','mallin_2023','jemaa_2017','jemaa_2023','hosny_2026']:
    directory=BASE/study/'extraction'
    for required in ['protocol.json','METHODS.md','digitized_points.csv','measurements.csv']:
        assert (directory/required).is_file(),(study,required)
    protocol=json.loads((directory/'protocol.json').read_text())
    counts={}
    for p in directory.rglob('*.csv'):
        with p.open(newline='') as f:
            rr=list(csv.DictReader(f));counts[str(p.relative_to(directory))]=len(rr)
    results[study]={'csv_rows':counts,'protocol_valid_json':True}
    source_parents=[r['path'] for r in inventory.values() if r['study']==study and r['role'] not in ['extraction_output','extraction_code'] and '/extraction/' not in r['path']]
    for p in (BASE/study).rglob('*'):
        if not p.is_file() or 'work' in p.relative_to(BASE/study).parts or '__pycache__' in p.parts:continue
        rel=str(p.relative_to(BASE))
        if '/extraction/' not in rel and rel in inventory:continue
        old=inventory.get(rel,{})
        details=new_sources.get(rel,{})
        role='extraction_code' if p.suffix=='.py' else 'extraction_output'
        if '/extraction/' not in rel:role=details.get('role','new_source_asset')
        inventory[rel]={key:old.get(key,'') for key in fields}
        parent=details.get('parent_asset') or (';'.join(source_parents) if p.name=='protocol.json' else f'{study}/extraction/protocol.json') if '/extraction/' in rel else details.get('parent_asset','')
        inventory[rel].update(path=rel,study=study,role=role,status=old.get('status') or 'available',source_url=old.get('source_url') or details.get('source_url',''),source_url_precision=old.get('source_url_precision') or ('exact' if details.get('source_url') else 'not_applicable'),parent_asset=old.get('parent_asset') or parent,bytes=p.stat().st_size,sha256=digest(p),asset_modified_utc=datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat(),notes=old.get('notes') or 'Extraction package; protocol.json and METHODS.md give precise evidence, source hashes and limitations.')

with (BASE/'xuan_2018/extraction/tracks.csv').open() as f:tracks=list(csv.DictReader(f))
keys=[(r['cell_id'],r['frame_index']) for r in tracks];assert len(keys)==len(set(keys))
with (BASE/'xuan_2018/extraction/equal_window_comparison.csv').open() as f:comp=list(csv.DictReader(f))
assert len(comp)==8 and all(r['first_frame']=='0' and r['last_frame']=='73' and r['n_frames']=='74' for r in comp)
assert all(float(r['path_length'])>=float(r['net_displacement']) and 0<=float(r['directness'])<=1 for r in comp)
with (BASE/'xuan_2018/extraction/msd.csv').open() as f:
    assert all(int(r['n_pairs'])==74-int(r['lag_frames']) for r in csv.DictReader(f))
results['xuan_2018'].update(unique_track_frame_keys=True,equal_window_tracks=8,msd_pair_counts_valid=True)
# Validate old assets too; checksum changes require an explicit new asset/version.
for rel,r in inventory.items():
    p=BASE/rel
    assert p.exists(),rel
    assert digest(p)==r['sha256'],rel
with (BASE/'manifest.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(inventory[k] for k in sorted(inventory))
results['generated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
results['manifest_assets']=len(inventory)
(BASE/'audit/extraction_validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
