"""Project accepted A4/forcing descriptors for concurrent final S11 rendering.

No hypothesis decisions, A5 outcomes or A6 associations are needed. All tables
come from the unchanged frozen final-statistics observation selection and E1
forcing rules, after complete E1--E4 native/estimand gates.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys


def context_rows(common,longest,e1):
    extra=[dict(r,atlas_context='E1_extra') for r in e1 if float(r['sigma']) in (.25,.5,2.)]
    long=[dict(r,atlas_context='E4_long') for r in longest if r['selected_observation_arm']=='preselected_E4']
    if len({(r['model'],float(r['C_label']),float(r['f']),float(r['sigma'])) for r in extra})!=72:
        raise ValueError('S11 E1_extra requires every24 conditions at all three preset extra sigma values')
    if len({(r['model'],float(r['C_label']),float(r['f'])) for r in long})!=52:
        raise ValueError('S11 E4_long requires every52 prespecified long conditions')
    return {'common':common,'E1_extra':extra,'E4_long':long}


def export(bundle,spec_path,out):
    if out.exists():raise FileExistsError('Use a fresh accepted atlas-context export')
    active=json.loads((bundle/'ACTIVE_FINAL_ANALYSIS.json').read_text());freeze=Path(active['active_source_freeze'])
    h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if h(freeze)!=active['active_source_freeze_sha256']:raise ValueError('The active reviewed final source changed')
    definition=json.loads(freeze.read_text());source=Path(definition['source'])
    sys.path.insert(0,str(source.parent));loader=importlib.util.spec_from_file_location('frozen_figure4_final_statistics',source);F=importlib.util.module_from_spec(loader);loader.loader.exec_module(F)
    F.verify_source_freeze(bundle);spec=F.checked_json(spec_path);data=F.assemble_input_data(spec)
    directory=bundle/'results/assembled_observations/common';status_path=directory/'status.json';status=F.checked_json(status_path)
    if not status.get('passed') or not status.get('complete'):raise ValueError('The frozen common observation table must be complete')
    common=F.checked_table(directory,status,'condition_metrics.csv')
    expected=hashlib.sha256(F.csv_bytes([dict(r,analysis_arm='common') for r in data['common']])).hexdigest()
    if status['outputs']['condition_metrics.csv']!=expected:raise ValueError('Atlas common observations differ from accepted observation selection')
    rows=context_rows(common,data['longest'],data['e1']);out.mkdir(parents=True)
    shutil.copy2(Path(__file__),out/Path(__file__).name)
    contexts={'common':dict(path=str(directory/'condition_metrics.csv'),sha256=h(directory/'condition_metrics.csv'),status_path=str(status_path),status_sha256=h(status_path))}
    for name in ('E1_extra','E4_long'):
        folder=out/name;folder.mkdir();table=folder/'condition_metrics.csv';F.write_csv(table,rows[name]);gate=folder/'status.json'
        F.write_json(gate,dict(passed=True,complete=True,context=name,source_freeze_sha256=h(freeze),source_stage_gates={key:value['gates'] for key,value in data['datasets'].items()},input_spec_sha256=h(spec_path),
            outputs={table.name:h(table)},scope='Exact accepted native raw A4 regimes and frozen mean/E1/combined operational flags; no final hypotheses or A5/A6 outcomes'))
        contexts[name]=dict(path=str(table),sha256=h(table),status_path=str(gate),status_sha256=h(gate))
    F.write_json(out/'metrics_plan.json',dict(passed=True,complete=True,contexts=contexts,normalized_stage_roots={name:row['root'] for name,row in spec['stages'].items()},source_freeze_sha256=h(freeze),projection_source_sha256=h(out/Path(__file__).name),input_spec_sha256=h(spec_path),no_native_integrations=True,no_final_hypothesis_selection=True))
    return dict(passed=True,complete=True,path=str(out/'metrics_plan.json'),sha256=h(out/'metrics_plan.json'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,required=True);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    print(json.dumps(export(args.bundle.resolve(),args.spec.resolve(),args.out.resolve()),indent=2))
