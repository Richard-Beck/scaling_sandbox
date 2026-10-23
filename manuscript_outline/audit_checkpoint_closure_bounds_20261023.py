"""Audit Claude's missing-angle concern over the entire preceding cohort."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4/supporting/data'
DEST = ROOT / 'manuscript_outline/revision_figures3_4_20261023'


def rows(name):
    with (SOURCE / name).open() as stream:
        return list(csv.DictReader(stream))


def key(row):
    return row['model'], float(row['C_label']), float(row['f'])


def main():
    coverage = {key(r): r for r in rows('jump_scoring_coverage.csv')}
    audited = []
    for r in rows('drift_decomposition.csv'):
        if float(r['flicker_frames']) != 5 or float(r['MSD_range_fraction']) != .25:
            continue
        kobs = float(r['kobs'])
        if not math.isfinite(float(r['kpred_over_kobs'])) or kobs <= 0:
            continue
        c = coverage[key(r)]
        drift = (2 * math.pi / float(r['C_mean'])) ** 2 * float(r['Ds'])
        low = (drift + float(c['lambda_eff_missing_angle_lower'])) / kobs
        high = (drift + float(c['lambda_eff_missing_angle_upper'])) / kobs
        audited.append(dict(model=r['model'], C_label=r['C_label'], f=r['f'],
                            multiple_front_fraction=r['multiple_front_fraction'],
                            valid_fraction=r['valid_fraction'], scored_fraction=c['jump_scored_fraction'],
                            point_ratio=r['kpred_over_kobs'], missing_angle_lower_ratio=low,
                            missing_angle_upper_ratio=high, unity_within_bounds=low <= 1 <= high,
                            within15_possible=low <= 1.15 and high >= .85))
    with (DEST / 'checkpoint_closure_bounds_audit.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(audited[0]))
        writer.writeheader()
        writer.writerows(audited)
    multiple = [r for r in audited if float(r['multiple_front_fraction']) > .5]
    goryachev = [r for r in audited if r['model'] == 'goryachev']
    result = dict(estimable=len(audited),
                  multiple_front_majority=dict(conditions=len(multiple),
                      unity_within_missing_angle_bounds=sum(r['unity_within_bounds'] for r in multiple),
                      within15_possible=sum(r['within15_possible'] for r in multiple),
                      point_within15=sum(abs(float(r['point_ratio']) - 1) <= .15 for r in multiple)),
                  goryachev=dict(conditions=len(goryachev),
                      point_ratio_min=min(float(r['point_ratio']) for r in goryachev),
                      point_ratio_max=max(float(r['point_ratio']) for r in goryachev),
                      unity_within_bounds=sum(r['unity_within_bounds'] for r in goryachev)),
                  source_sha256={str((SOURCE / n).relative_to(ROOT)):
                                 hashlib.sha256((SOURCE / n).read_bytes()).hexdigest()
                                 for n in ('drift_decomposition.csv', 'jump_scoring_coverage.csv')})
    (DEST / 'checkpoint_closure_bounds_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
