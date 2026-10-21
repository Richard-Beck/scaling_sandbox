"""Independent point and replication review of primary optical A5 summaries.

This does not reconstruct region tracks or bootstrap confidence intervals.
Those require the separately bound independent all-view native-frame audit.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def review(root, output):
    root = Path(root).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError('An independent review receipt is immutable')
    inventory_path = root / 'inventory.json'
    inventory = read(inventory_path)
    complete_path = root / 'results/analysis/complete.json'
    complete = read(complete_path)
    audit_path = root / 'results/independent_all_views_audit.json'
    audit = read(audit_path)
    inventory_sha = sha(inventory_path)
    cases = inventory['cases']
    if (complete.get('passed') is not True or audit.get('passed') is not True
            or audit.get('complete') is not True
            or complete.get('inventory_sha256') != inventory_sha
            or audit.get('inventory_sha256') != inventory_sha
            or complete.get('cases') != len(cases) or audit.get('cases') != len(cases)):
        raise ValueError('Completed canonical collection and all-case audit required')
    audit_rows = {r['case_id']: r for r in audit['records']}
    if len(audit_rows) != len(cases) or set(audit_rows) != {c['id'] for c in cases}:
        raise ValueError('Independent native audit omits or repeats a canonical case')
    groups = defaultdict(list)
    for case in cases:
        path = root / 'workers' / case['id'] / 'result.json'
        record = read(path)
        checked = audit_rows[case['id']]
        if (record.get('passed') is not True or checked.get('passed') is not True
                or record.get('case') != case
                or record.get('inventory_sha256') != inventory_sha
                or sha(path) != checked['measurement_sha256']
                or record['primary']['sha256'] != checked['primary_sha256']):
            raise ValueError('Canonical primary measurement is not independently bound')
        primary = record['primary']
        if (primary['threshold'] != .5 or primary['smoothed'] is not True
                or record['raw_profile_cadence'] * primary['downsample_factor']
                != record['common_imaging_cadence']):
            raise ValueError('Primary observation threshold, optics or cadence changed')
        summary = primary['summary']
        distribution = summary['region_count_time_distribution']
        values = dict(
            region_count_mean=sum(int(k) * v for k, v in distribution.items()),
            multiple_front_fraction=sum(v for k, v in distribution.items() if int(k) >= 2),
            low_contrast_time_fraction=summary['low_contrast_time_fraction'],
            mean_region_width=summary['mean_region_width'],
            mean_region_width_fraction=summary['mean_region_width_fraction'],
            births_per_time_per_perimeter=summary['observed_births'] / summary['perimeter_time_exposure'])
        key = (case['model'], float(case['C_label']), float(case['f']), case['kind'],
               float(case['sigma']), float(case['native_horizon']))
        groups[key].append(values)
    metrics_path = root / 'results/analysis/metrics.csv'
    if sha(metrics_path) != complete['outputs']['metrics']['sha256']:
        raise ValueError('Completed primary metrics changed')
    errors, checked_points, found = [], 0, set()
    with metrics_path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            if row['primary'] != 'True':
                continue
            key = (row['model'], float(row['C_label']), float(row['f']), row['kind'],
                   float(row['sigma']), float(row['native_horizon']))
            if key not in groups:
                raise ValueError('Unexpected primary condition')
            samples = groups[key]
            statistic = row['statistic']
            if statistic not in samples[0]:
                continue
            pair = (key, statistic)
            if pair in found:
                raise ValueError('Repeated primary condition/statistic')
            found.add(pair)
            finite = [s[statistic] for s in samples if s[statistic] is not None
                      and math.isfinite(s[statistic])]
            expected = math.fsum(finite) / len(finite) if finite else None
            actual = float(row['estimate']) if row['estimate'] else None
            if expected is None:
                if actual is not None:
                    raise AssertionError('Unavailable primary estimate became a value')
            else:
                if actual is None or not math.isclose(actual, expected, rel_tol=2e-12, abs_tol=2e-12):
                    raise AssertionError('Independent primary estimate differs')
                errors.append(abs(actual - expected))
            deterministic = row['kind'] != 'noise'
            if (int(row['independent_seeds']) != (0 if deterministic else len(samples))
                    or int(row['deterministic_controls']) != (1 if deterministic else 0)
                    or int(row['bootstrap_resamples']) != (0 if deterministic else 2000)
                    or int(row['seeds_with_estimand']) != len(finite)
                    or (deterministic and len(samples) != 1)):
                raise AssertionError('Deterministic controls or independent seed counts changed')
            checked_points += 1
    if checked_points != len(groups) * 6 or complete['groups'] != len(groups):
        raise ValueError('All six primary statistics and every canonical group are required')
    receipt = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        scope='Independent point-estimate and replication review of six primary optical A5 statistics; not an independent reconstruction of all tracking or bootstrap intervals.',
        cases=len(cases), groups=len(groups), primary_estimates_checked=checked_points,
        maximum_absolute_error=max(errors, default=0.), inventory_sha256=inventory_sha,
        metrics_sha256=sha(metrics_path), complete_sha256=sha(complete_path),
        all_views_native_audit_sha256=sha(audit_path), reviewer_source_sha256=sha(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(review(args.root, args.output), indent=2))
