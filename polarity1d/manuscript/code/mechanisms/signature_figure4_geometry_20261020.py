"""Additive geometric qualifiers from independently accepted final ratios.

Time-based excess persistence and short-lag physical mobility are separate
conclusions. Conflicting conclusions remain visible. These labels never alter
the frozen H1--H3 criteria, and an interval including one is compatibility,
not evidence of equivalence.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(4<<20),b''):h.update(b)
    return h.hexdigest()


def boolean(v):
    if v in (True,'True','true'):return True
    if v in (False,'False','false'):return False
    return None


def number(v):
    try:r=float(v)
    except (TypeError,ValueError):return None
    return r if math.isfinite(r) else None


def finite_interval(row,lower,upper):
    lo,hi=number(row.get(lower)),number(row.get(upper))
    return (lo,hi) if lo is not None and hi is not None and 0<=lo<=hi else None


def classify(row,inclusive=False):
    force=boolean(row.get('forced_combined'));geom=number(row.get('actual_geometry_ratio'))
    qualification='included' if inclusive or force is False else 'excluded_forced' if force else 'unavailable_forced_qualification'
    if geom is None or geom<=0:qualification='unavailable_unmatched_or_missing_geometry'
    time=finite_interval(row,'excess_persistence_ci_lower','excess_persistence_ci_upper')
    dx=finite_interval(row,'Dx_ratio_ci_lower','Dx_ratio_ci_upper');point=number(row.get('Dx_large_over_small'))
    time_status='available' if time is not None else 'unavailable_missing_or_unbounded_interval'
    if boolean(row.get('excess_persistence_ci_upper_unbounded')) is True:time_status='unavailable_unbounded_interval'
    if boolean(row.get('evaluable')) is not True:time_status='unavailable_persistence_qualification'
    mobility_status='available' if dx is not None and point is not None and point>0 else 'unavailable_missing_positive_ratio_or_finite_interval'
    if qualification!='included':time_status=mobility_status=qualification
    time_label='unavailable' if time_status!='available' else 'time_exceeds_geometry' if time[0]>1 else 'time_below_geometry' if time[1]<1 else 'time_CI_includes_geometry'
    mobility_label='unavailable' if mobility_status!='available' else 'Dx_decreases_with_size' if dx[1]<1 else 'Dx_increases_with_size' if dx[0]>1 else 'Dx_CI_includes_flat'
    conflict=(time_label,mobility_label) in [('time_exceeds_geometry','Dx_increases_with_size'),('time_below_geometry','Dx_decreases_with_size')]
    if conflict:label='conflicting_time_and_mobility'
    elif time_label=='time_exceeds_geometry' or mobility_label=='Dx_decreases_with_size':label='more_than_geometry'
    elif mobility_label=='Dx_increases_with_size':label='less_than_geometry'
    elif mobility_label=='Dx_CI_includes_flat':label='compatible_with_flat_Dx'
    elif time_label=='time_CI_includes_geometry':label='time_compatible_mobility_unavailable'
    elif time_label=='time_below_geometry':label='time_below_geometry_mobility_unavailable'
    else:label='unavailable'
    return dict(row,forced_inclusive=inclusive,geometric_qualification=qualification,time_geometry_status=time_status,mobility_geometry_status=mobility_status,
        time_geometry_conclusion=time_label,mobility_geometry_conclusion=mobility_label,conflicting_time_mobility=conflict,geometric_qualifier=label,
        geometry_is_not_hypothesis_gate=True,interval_includes_one_interpretation='Compatibility, not equivalence',
        mobility_inference_scope='Positive short-lag physical mobility ratio and finite interval; valid-bootstrap fraction retained separately, not substituted for tau endpoint coverage')


def model_rows(rows):
    grouped=defaultdict(list)
    for r in rows:grouped[(r.get('arm',r.get('analysis_arm')),r['view'],r['model'],r['forced_inclusive'])].append(r)
    out=[]
    for (arm,view,model,inclusive),values in sorted(grouped.items()):
        labels=Counter(r['geometric_qualifier'] for r in values);available=set(labels)-{'unavailable'}
        pattern='unavailable' if not available else next(iter(available)) if len(available)==1 else 'heterogeneous_across_inputs'
        lists={label:sorted(float(r['f']) for r in values if r['geometric_qualifier']==label) for label in sorted(labels)}
        out.append(dict(arm=arm,view=view,model=model,forced_inclusive=inclusive,designed_inputs=len(values),included_inputs=sum(r['geometric_qualification']=='included' for r in values),
            finite_time_intervals=sum(r['time_geometry_status']=='available' for r in values),finite_mobility_intervals=sum(r['mobility_geometry_status']=='available' for r in values),
            conflicting_inputs=sum(r['conflicting_time_mobility'] for r in values),per_model_pattern=pattern,qualifier_counts=dict(labels),qualifier_input_strengths=lists,
            time_conclusion_counts=dict(Counter(r['time_geometry_conclusion'] for r in values)),mobility_conclusion_counts=dict(Counter(r['mobility_geometry_conclusion'] for r in values)),
            geometry_is_not_hypothesis_gate=True,no_equivalence_claim=True))
    return out


def write_csv(path,rows):
    fields=sorted({k for row in rows for k in row})
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for row in rows:writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def export(parent,out):
    complete_path=parent/'complete.json';audit_path=parent/'independent_final_decision_audit.json'
    complete=json.loads(complete_path.read_text());audit=json.loads(audit_path.read_text())
    if not complete.get('passed') or not complete.get('complete') or not audit.get('passed'):raise ValueError('The complete independently accepted final analysis is required')
    if audit.get('producer_complete_sha256')!=digest(complete_path):raise ValueError('The independent final audit belongs to a different producer')
    ratio=parent/'size_ratios.csv'
    if digest(ratio)!=complete['outputs'][ratio.name]:raise ValueError('Accepted final size ratios changed')
    if out.exists():raise FileExistsError('Use a fresh additive geometry export')
    with ratio.open(newline='') as stream:source=list(csv.DictReader(stream))
    rows=[classify(row,inclusive) for inclusive in (False,True) for row in source];models=model_rows(rows);out.mkdir(parents=True);shutil.copy2(Path(__file__),out/Path(__file__).name)
    write_csv(out/'geometric_conditions.csv',rows);write_csv(out/'geometric_models.csv',models)
    status=dict(passed=True,complete=True,source_complete_path=str(complete_path),source_complete_sha256=digest(complete_path),source_independent_audit_path=str(audit_path),source_independent_audit_sha256=digest(audit_path),
        source_size_ratios_path=str(ratio),source_size_ratios_sha256=digest(ratio),source_final_freeze_sha256=complete['source_freeze_sha256'],projection_source_sha256=digest(out/Path(__file__).name),
        outputs={name:digest(out/name) for name in ('geometric_conditions.csv','geometric_models.csv')},condition_rows=len(rows),model_rows=len(models),scope='Additive geometric qualifier only; all308 raw H1-H3 criterion rows and source closure unchanged',
        compatibility_is_not_equivalence=True,unbounded_intervals_do_not_establish_geometric_labels=True,conflicting_time_mobility_conclusions_retained=True,no_native_integrations=True)
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n');return dict(passed=True,complete=True,status_sha256=digest(out/'status.json'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(export(a.parent.resolve(),a.out.resolve()),indent=2))
