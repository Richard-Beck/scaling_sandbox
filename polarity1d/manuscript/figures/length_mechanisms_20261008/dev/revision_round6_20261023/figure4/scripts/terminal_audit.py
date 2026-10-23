"""Audit the allocated compact evidence without repeating native integrations."""
from common import *
from datetime import datetime,timezone
from collections import Counter

def audit():
    manifest=read(RUN/'manifest.json');ledger=read(RUN/'allocation_targets.json')
    rows=[];failures=[];runtime=Counter();seeds=set()
    for ci,c in enumerate(manifest['conditions']):
        target=ledger['targets'][str(ci)]
        for rep in range(target):
            folder=RUN/'workers'/f'c{ci:03d}'/f'seed{rep:03d}';status_path=folder/'status.json';archive=folder/'modes.npz'
            status=read(status_path) if status_path.exists() else {}
            expected_seed=202610230600000+c['condition_index']*1000+rep
            valid=bool(status.get('complete')) and archive.exists() and archive.stat().st_size==status.get('bytes') and status.get('seed')==expected_seed
            if not valid:failures.append(str(folder))
            if expected_seed in seeds:failures.append('duplicate seed '+str(expected_seed))
            seeds.add(expected_seed)
            runtime[status.get('runtime_manifest_sha256','initial/runtime launch hash not recorded')]+=1
            checks=status.get('validation',{})
            rows.append(dict(condition_id=c['condition_id'],replicate=rep,seed=expected_seed,valid=valid,archive_bytes=archive.stat().st_size if archive.exists() else 0,
                recorded_minimum_concentration=checks.get('minimum_concentration'),recorded_maximum_pool_relative_error=checks.get('maximum_pool_relative_error'),
                runtime_manifest_sha256=status.get('runtime_manifest_sha256','not recorded at initial launch')))
    assert len(list(RUN.glob('workers/c*/seed*/modes.npz')))==len(rows),'Unexpected extra mode archives outside allocation ledger'
    source=read(SOURCE/'manifest.json');source_checks=[]
    for kind in ('source_sha256','input_sha256'):
        for rel,expected in source[kind].items():
            ok=digest(SOURCE/rel)==expected;source_checks.append(ok)
            if not ok:failures.append(str(SOURCE/rel))
    frozen=OUT/'scripts/frozen_native_v2';frozen_manifest=read(frozen/'runtime_manifest.json')
    for name,expected in frozen_manifest['files'].items():
        if digest(frozen/name)!=expected:failures.append(str(frozen/name))
    snapshot=OUT/'supporting/review_snapshot_initial';snapshot_manifest=read(snapshot/'SNAPSHOT_HASHES.json')
    for name,expected in snapshot_manifest['files'].items():
        if digest(snapshot/name)!=expected:failures.append(str(snapshot/name))
    assert digest(SOURCE/'manifest.json')==manifest['source_sha256']
    assert digest(CODE/'figure4_round3_campaign_20261022.py')==manifest['source_code_sha256']
    assert not failures,failures
    csvwrite(META/'terminal_seed_receipts.csv',rows)
    result=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,allocated_conditions=len(ledger['targets']),allocated_histories=sum(ledger['targets'].values()),
        complete_valid_receipts=len(rows),unique_seeds=len(seeds),native_profiles_saved=False,source_hash_checks=len(source_checks),source_hashes_all_match=all(source_checks),
        immutable_snapshot_files_verified=len(snapshot_manifest['files']),frozen_late_runtime_files_verified=len(frozen_manifest['files']),runtime_provenance_counts=dict(runtime),
        archive_validation='Final precision rereads all four arrays in every compressed archive, including ZIP CRC verification. This audit independently checks completion, bytes, seed identity and retained numerical receipts; SHA256 values remain in per-seed status records.',
        limitation='Initial full-wrapper launch hash was not recorded; recovered source history is labelled as recovery evidence. Native kernel and SOURCE manifest were contemporaneously hash-recorded.')
    write(META/'terminal_evidence_audit.json',result)
    return result
if __name__=='__main__':print(json.dumps(audit(),indent=2))
