"""Measure allocated round6 storage and maintain its dated manifest entry."""
from pathlib import Path
from datetime import datetime, timezone
import json
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BEGIN = '<!-- round6-storage-20261023 begin -->'
END = '<!-- round6-storage-20261023 end -->'


def main():
    stamp = datetime.now(timezone.utc).isoformat()
    groups = [
        ('polarity1d/manuscript/data/runs/figure4_round6_modes_20261023',
         'Targeted modes-only native noise rerun; a₀–a₃₂, |P|, time/circumference and forcing audits. No durable activity profiles. Authoritative evidence for the detector-free persistence test: partial, slope-stencil-sensitive agreement and the failed Goryachev slope screen do not establish a general strong-input explanation. Completion and cap-limited precision are in the figure receipt.'),
        ('polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023',
         'Figures 1–4, preserved earlier supplements, compact plotted CSVs, deterministic susceptibility tests and Fourier/exact-mask analyses. Figures 1–2 are layout-only; new model tests are in Figures 3–4.'),
        ('manuscript_outline/revision_round6_20261023',
         'Verbatim incoming brief, prompt and narrative; starting-outline hashes, standalone legends, integrated reader, response and validation records.'),
        ('manuscript_outline/outline_v12.pdf',
         'Integrated manuscript outline starting from outline_v11. Published migration evidence, reserved WGD experiment and prospective pilot appendix are preserved.'),
        ('.git', 'Versioned scientific source, selected portable evidence and history; bulk mode-only histories remain local outside Git.'),
    ]
    records = []
    for name, purpose in groups:
        path = ROOT / name
        if not path.exists():
            continue
        output = subprocess.check_output(['du', '-sx', '-B1', str(path)], text=True)
        allocated = int(output.split()[0])
        records.append({'path': name, 'allocated_bytes': allocated,
                        'allocated_GiB': allocated / 2**30, 'purpose_and_status': purpose})
    run = ROOT / groups[0][0]
    statuses = list((run / 'workers').glob('c*/seed*/status.json'))
    complete = failed = 0
    for path in statuses:
        value = json.loads(path.read_text())
        complete += bool(value.get('complete'))
        failed += bool(value.get('error'))
    precision_path = run / 'precision/status.json'
    precision = json.loads(precision_path.read_text()) if precision_path.exists() else None
    ledger_path = run / 'allocation_targets.json'
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else {'targets': {}}
    allocated_histories = sum(ledger['targets'].values()) or 4800
    campaign_status = 'ACTIVE' if complete < allocated_histories or not precision or precision.get('early', True) else 'NATIVE_COHORT_COMPLETE'
    result = {'measurement_started_utc': stamp,
        'measurement_finished_utc': datetime.now(timezone.utc).isoformat(),
        'method': 'du -sx -B1; allocated bytes, targeted non-atomic traversal',
        'groups': records, 'mode_seed_status_files': len(statuses),
        'complete_mode_histories': complete, 'failed_mode_histories': failed,
        'allocated_mode_histories': allocated_histories, 'campaign_status': campaign_status,
        'precision_status_available': precision is not None,
        'no_data_removed': True, 'planning_initial_bytes': 4321804800,
        'planning_maximum_bytes': 17287219200}
    receipt = ROOT / 'reports/round6_storage_20261023.json'
    receipt.write_text(json.dumps(result, indent=2) + '\n')
    text = f'''{BEGIN}

## Round6 susceptibility and mode-only revision, 23 October 2026

Measured **{result['measurement_finished_utc']}** using allocated storage
(`du -sx -B1`). These targeted measurements update the named output groups;
the earlier whole-workspace totals remain historical. This is a non-atomic
snapshot. Inventory confers no deletion permission; no original data were removed.

| Path / related output group | Measured allocated storage | Purpose and scientific status |
| --- | ---: | --- |
'''
    for row in records:
        text += f"| `{row['path']}` | {row['allocated_GiB']:.6f} GiB ({row['allocated_bytes']:,} bytes) | {row['purpose_and_status']} |\n"
    text += f'''
The modes-only campaign status is **{campaign_status}**, with
**{complete:,} completed / {allocated_histories:,} allocated histories** and
**{failed:,} failed status records** at this measurement. The prelaunch estimate
was 4,321,804,800 bytes for 24 seeds per condition and 17,287,219,200 bytes at
the 96-seed cap, before compression. Whole-seed precision determines additions.
Final scientific eligibility and unresolved rates remain in the Figure 4 tables.
Deterministic Figure 3 maps, pulse histories and source states are grouped under
the new figure parent above. Round5 and the original compact noise campaign
remain preserved; their dated inventory entries still apply.

[Machine-readable measurement](reports/round6_storage_20261023.json).

{END}
'''
    manifest = ROOT / 'DATA_MANIFEST.md'
    previous = manifest.read_text()
    if BEGIN in previous:
        before, rest = previous.split(BEGIN, 1)
        _, after = rest.split(END, 1)
        updated = before + text + after
    else:
        updated = previous.rstrip() + '\n\n' + text
    manifest.write_text(updated.rstrip() + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
