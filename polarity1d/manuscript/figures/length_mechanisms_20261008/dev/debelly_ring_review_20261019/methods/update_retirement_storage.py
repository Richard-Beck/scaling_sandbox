"""Refresh only this completed retirement's storage accounting and report date."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, importlib.util, json, re

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = next(p for p in PACKAGE.parents if (p / 'DATA_MANIFEST.md').exists())
spec = importlib.util.spec_from_file_location('inventory', PACKAGE / 'cleanup_inventory.py')
inventory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory)
GIB = 1024 ** 3

def main():
    receipt_path = PACKAGE / 'cleanup_receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['passed'] and all(not inventory.ROOT.joinpath(r['path']).exists() if not Path(r['path']).is_absolute() else not Path(r['path']).exists() for r in receipt['removed'])
    timestamp = datetime.now(timezone.utc)
    iso = timestamp.isoformat()
    target = PACKAGE / 'supporting/metadata/storage_footprint.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    index = PACKAGE / 'index.html'
    report = index.read_text()
    report = re.sub(r'<p class="meta">Consolidated review.*?</p>',
        f'<p class="meta">Consolidated review · Generated {timestamp:%Y-%m-%d %H:%M:%S} UTC · Scientific development stopped · Static report; no automatic refresh</p>', report, count=1)
    index.write_text(report)
    for _ in range(6):
        measured = inventory.inventory(PACKAGE)
        retained = measured['allocated_bytes']
        net = receipt['old_workspace_allocated_bytes'] - retained
        receipt.update(retained_package_allocated_bytes=retained, retained_package_measured_utc=iso,
            net_workspace_allocated_reduction_bytes=net,
            retained_coverage='All 18 event clocks, native 0.5 s windows, 375 whole-response context frames, selected summaries/figures, minimal historical model and basal restart.')
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        storage = dict(measured_utc=iso, scope='This dated attempt only; no whole-workspace or parent refresh.',
            retained=dict(path=str(PACKAGE.relative_to(ROOT)), **measured),
            removed_workspace_allocated_bytes=receipt['old_workspace_allocated_bytes'],
            removed_system_temporary_allocated_bytes=receipt['old_system_temporary_allocated_bytes'],
            net_workspace_allocated_reduction_bytes=net,
            retirement_inventory='../../cleanup_plan.json', retirement_receipt='../../cleanup_receipt.json',
            purpose='Closed-out unsuccessful periodic deBelly extension; viewer and selected evidence preserved. Corrected linear, earlier evidence and independent Spring/figure work protected.',
            source_sha256={str(Path(__file__).relative_to(PACKAGE)): hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
        target.write_text(json.dumps(storage, indent=2) + '\n')
        paragraph = (f'<p id="storage-status">Measured {timestamp:%Y-%m-%d %H:%M:%S} UTC: '
            f'<strong>{retained / 1024**2:.1f} MiB retained</strong>. '
            f'{receipt["old_workspace_allocated_bytes"] / GIB:.3f} GiB of superseded workspace artifacts and '
            f'{receipt["old_system_temporary_allocated_bytes"] / 1024**2:.1f} MiB of task temporary files were retired. '
            f'Net workspace allocated-storage reduction after retaining this package: <strong>{net / GIB:.3f} GiB</strong>. '
            'These are measured filesystem allocations, not a filesystem free-space estimate.</p>')
        report = re.sub(r'<p id="storage-status">.*?</p>', paragraph, index.read_text(), count=1)
        index.write_text(report)
        if inventory.inventory(PACKAGE)['allocated_bytes'] == retained:
            break
    else:
        raise RuntimeError('Storage accounting did not stabilize.')
    manifest = ROOT / 'DATA_MANIFEST.md'
    text = manifest.read_text()
    start = text.index('### deBelly ring development addition, 2026-10-19') if '### deBelly ring development addition, 2026-10-19' in text else text.index('### deBelly periodic-domain attempt retirement, 2026-10-19')
    end = text.index('\n### ', start + 4)
    relative = str(PACKAGE.relative_to(ROOT))
    section = f'''### deBelly periodic-domain attempt retirement, 2026-10-19

The user requested consolidation and removal of artifacts from this attempt.
This is a measurement of the named retired outputs and retained package only;
the manuscript parent measurements above have not been refreshed. The
[retirement inventory]({relative}/cleanup_plan.json),
[execution receipt]({relative}/cleanup_receipt.json) and
[storage record]({relative}/supporting/metadata/storage_footprint.json)
record the exact paths and measurements. Inventory alone did not authorize deletion;
the user's explicit cleanup instruction did.

| Path or related group | Allocated GiB | Purpose and importance/status |
| --- | ---: | --- |
| `{relative}/` | **{retained / GIB:.3f} retained** ({retained:,} allocated bytes; measured {timestamp:%Y-%m-%d, %H:%M:%S} UTC) | Single [executive report and seven-field reviewer]({relative}/index.html) for the incomplete periodic deBelly attempt. Retains all 18 event clocks from the original and two fixed-grid-face C50 challengers at native 0.5-s resolution in event windows, 375 full-precision whole-response context frames, four compact evidence groups, minimal frozen reference model/basal state, provenance and audits. Viewer regeneration and relocation pass. Physical cap/material/reference closure, causal travel/reversal mechanism and numerical convergence remain unresolved. Recommendation: omit deBelly from current ring assays and keep corrected linear deBelly. No migration-persistence or accepted ring size law resulted. |
| Superseded 19 October deBelly run/report/source/cache groups (exact list in inventory) | **0 remaining**; {receipt['old_workspace_allocated_bytes']:,} workspace allocated bytes retired ({receipt['old_workspace_allocated_bytes'] / GIB:.3f} GiB), plus {receipt['old_system_temporary_allocated_bytes']:,} temporary bytes | Removed 29 dated run roots, the old development report topic, 123 standalone scripts, 49 caches, the dated status note and 101 task temporary entries after independent retention/dependency checks. Full native histories outside retained event windows, context samples and selected summaries were retired; historical audit receipts do not imply those full inputs still exist. Net workspace allocation reduction after retaining the review: **{net / GIB:.3f} GiB**. Established corrected linear results, earlier independent diagnostics, `tmp/x50_mechanisms_20261003/debelly`, current figures and independent Spring work are protected. This component must not be added to a refreshed manuscript parent total. |
'''
    manifest.write_text(text[:start] + section + text[end:])
    print(json.dumps(dict(measured_utc=iso, retained_allocated_bytes=retained, net_workspace_reduction_bytes=net), indent=2))

if __name__ == '__main__':
    main()
