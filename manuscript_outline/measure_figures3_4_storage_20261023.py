"""Record allocated revision storage without traversing or changing raw data."""
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
GIB = 1073741824
MARKER = 'figures3-4-edge-birth-storage-20261023'


def main():
    runs = ROOT / 'polarity1d/manuscript/data/runs'
    items = [(ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round5_20261023',
              'New Figures 3–4, source tables, deterministic response summaries, certified-event records, exact signed heading accounting, supplements and validation. Derived data only; pre-existing scientific evidence is retained.'),
             (ROOT / 'manuscript_outline/revision_figures3_4_20261023',
              'Standalone figures and legends, pilot observables, archived verbatim prompt/brief, Claude response and preservation hashes, integration and link checks.'),
             (ROOT / 'manuscript_outline/outline_v11.pdf',
              'Integrated outline. Figures 1, 2, 5 and reserved Figure 6 preserve baseline sections and artwork.'),
             (runs / 'signature_figure4_round3_compact_20261022',
              'Pre-existing authoritative compact noise campaign: 260 conditions and 45,696 histories. Read for analysis; no new noisy simulations or changes to original scientific records.'),
             (ROOT / '.git',
              'Versioned scientific source, portable evidence and repository history, including the initial requested checkpoint. Native compact noise histories remain outside Git.')]
    for path in sorted(runs.glob('figure3_round5*')):
        if not path.is_dir():
            continue
        if path.name == 'figure3_round5_edge_pulses_v2_20261023':
            purpose = 'Active frozen deterministic campaign: 120 conditions; 118 polarized responses and two homogeneous Goryachev mean-state exclusions. Half/double amplitude checks; original strict solver results used.'
        else:
            purpose = 'Retained diagnostic preparation/API-failure or canceled numerical-retry record. Not used as scientific response evidence; original frozen files preserved.'
        items.append((path, purpose))
    records = []
    for path, purpose in items:
        assert path.exists(), path
        output = subprocess.check_output(['du', '-sx', '-B1', str(path)], text=True)
        size = int(output.split(maxsplit=1)[0])
        records.append(dict(path=str(path.relative_to(ROOT)), allocated_bytes=size,
                            allocated_GiB=size / GIB, purpose=purpose))
    date = datetime.now(timezone.utc).isoformat()
    report = ROOT / 'reports/figures3_4_round5_storage_20261023.json'
    receipt = dict(measured_utc=date, method='du -sx -B1 separately for each named path',
                   scope='Targeted measurements; not a new whole-workspace inventory. Separate groups can share hard links with earlier editions and are not additive.',
                   no_data_removed=True, records=records)
    report.write_text(json.dumps(receipt, indent=2) + '\n')
    rows = '\n'.join(f'| `{r["path"]}` | {r["allocated_GiB"]:.6f} GiB ({r["allocated_bytes"]:,} bytes) | {r["purpose"]} |'
                     for r in records)
    block = f'''<!-- {MARKER} begin -->

## Figures 3–4 edge response and secondary-domain revision, 23 October 2026

Measured **{date}** using allocated storage (`du -sx -B1`) for each named path. These targeted measurements supersede only their corresponding dated groups; the older whole-workspace totals are historical. Related new outputs are grouped below, including the diagnostic deterministic attempts. No original data were moved or removed, and no new noise campaign was generated. Separate group measurements can overlap through hard links to earlier derived exports; do not sum them as storage growth.

| Path / related output group | Measured allocated storage | Purpose and scientific status |
| --- | ---: | --- |
{rows}

[Machine-readable storage receipt](reports/figures3_4_round5_storage_20261023.json). Scientific limitations are explicit: births and splits are only partially identifiable from compact records; takeover probabilities condition on certification and finite follow-up; signed moment accounting is descriptive. Claude's reviews are preserved verbatim separately from the responses. Inventory confers no deletion permission.

<!-- {MARKER} end -->'''
    manifest = ROOT / 'DATA_MANIFEST.md'
    text = manifest.read_text()
    start = f'<!-- {MARKER} begin -->'
    end = f'<!-- {MARKER} end -->'
    if start in text:
        a = text.index(start)
        b = text.index(end, a) + len(end)
        text = text[:a] + block + text[b:]
    else:
        text = text.rstrip() + '\n\n' + block + '\n'
    manifest.write_text(text)
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
