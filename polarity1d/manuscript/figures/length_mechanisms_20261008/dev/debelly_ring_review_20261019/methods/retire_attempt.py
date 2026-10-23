"""Execute this attempt's reviewed, exact-path retirement inventory once.

Default is a read-only preflight. Execution requires the current plan's SHA256
and refuses to run again after writing a receipt. No glob expands deletion scope.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, collections, hashlib, importlib.util, json, os, shutil

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = next(p for p in PACKAGE.parents if (p / 'polarity1d/manuscript/AGENTS.md').exists())
PLAN = PACKAGE / 'cleanup_plan.json'
RECEIPT = PACKAGE / 'cleanup_receipt.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def resolve(label):
    p = Path(label)
    assert '..' not in p.parts, label
    return p if p.is_absolute() else ROOT / p

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-sha256', required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    assert sha(PLAN) == args.plan_sha256, 'Inventory changed; review the new plan first.'
    assert not RECEIPT.exists(), 'A retirement receipt already exists. Refusing another execution.'
    d = json.loads(PLAN.read_text())
    spec = importlib.util.spec_from_file_location('retirement_inventory', PACKAGE / 'cleanup_inventory.py')
    inventory = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inventory)
    candidates = [(resolve(r['path']), r) for r in d['candidate_delete']]
    protected = [resolve(r['path']) for r in d['protected']]
    for p, r in candidates:
        assert p != PACKAGE and PACKAGE not in p.parents and p not in PACKAGE.parents, p
        assert not any(p == q or p in q.parents or q in p.parents for q in protected), ('protected overlap', p)
        assert not any(p != q and (p in q.parents or q in p.parents) for q, _ in candidates), ('candidate overlap', p)
        assert r['tracking']['tracked_files'] == 0, p
        assert p.exists() or p.is_symlink(), p
        if r['category'] == 'system_temp':
            assert p.parent == Path('/tmp') and 'debelly' in p.name.lower() and p.lstat().st_uid == os.getuid(), p
        elif r['category'] == 'repository_temp':
            raise AssertionError('Unscoped repository temporary paths are not authorized by this inventory.')
        else:
            assert ROOT in p.parents and '20261019' in p.name, p
        if 'sha256' in r:
            assert sha(p) == r['sha256'], ('candidate file changed', p)
        for label, digest in r.get('provenance_hashes', {}).items():
            assert sha(resolve(label)) == digest, ('candidate provenance changed', label)
        actual = inventory.inventory(p)
        assert actual['allocated_bytes'] == r['allocated_bytes'], ('candidate size changed', p)
    hashes = {}
    for r in d['protected']:
        if 'sha256' in r:
            hashes[r['path']] = sha(resolve(r['path']))
        for label in r.get('protected_hash_samples', {}):
            hashes[label] = sha(resolve(label))
    retained_before = {str(p.relative_to(PACKAGE)): sha(p) for p in PACKAGE.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    summary = dict(plan_sha256=sha(PLAN), candidates=len(candidates),
                   categories=dict(collections.Counter(r['category'] for _, r in candidates)),
                   candidate_allocated_bytes=sum(r['allocated_bytes'] for _, r in candidates),
                   protected_hash_samples=len(hashes), retained_files=len(retained_before), execute=args.execute)
    if not args.execute:
        print(json.dumps(summary, indent=2))
        return
    rows = []
    for p, r in candidates:
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        else:
            p.unlink()
        assert not p.exists() and not p.is_symlink(), p
        rows.append(dict(path=r['path'], category=r['category'], allocated_bytes=r['allocated_bytes'], absent=True))
    assert all(sha(resolve(label)) == digest for label, digest in hashes.items()), 'Protected sample changed during removal.'
    assert all(sha(PACKAGE / label) == digest for label, digest in retained_before.items()), 'Retained file changed during removal.'
    summary.update(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        mode='EXECUTED_USER_AUTHORIZED_EXACT_PATH_RETIREMENT', removed=rows,
        protected_samples_unchanged=True, retained_files_unchanged=True,
        old_workspace_allocated_bytes=sum(r['allocated_bytes'] for _, r in candidates if r['category'] != 'system_temp'),
        old_system_temporary_allocated_bytes=sum(r['allocated_bytes'] for _, r in candidates if r['category'] == 'system_temp'),
        full_native_histories_retired=True,
        retained_coverage='All18 event clocks, native0.5s windows,375 whole-response context frames,selected summaries/figures,minimal historical model and basal restart.',
        measurement_note='Allocated storage is measured from filesystem blocks; it is not a filesystem free-space estimate. Net workspace reduction subtracts the retained package separately.')
    RECEIPT.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'removed'}, indent=2))

if __name__ == '__main__':
    main()
