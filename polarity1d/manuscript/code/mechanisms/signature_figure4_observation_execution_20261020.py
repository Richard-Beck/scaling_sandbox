"""Run the frozen observation/context phase after complete recovered E4 intake.

This engineering wrapper only binds completed source identities and executes
the already frozen v3 estimands and descriptor projection. It implements no
hypothesis rule and performs no native integration.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())


def run(config_path,expected_config_sha256):
    if sha(config_path)!=expected_config_sha256:raise ValueError('The frozen observation execution configuration changed')
    config=read(config_path)
    if str(Path(__file__).resolve())!=config['execution_source']['path'] or sha(__file__)!=config['execution_source']['sha256']:raise ValueError('Execute the actual frozen observation wrapper')
    bundle=Path(config['bundle']);source_freeze=Path(config['source_freeze']['path'])
    if sha(source_freeze)!=config['source_freeze']['sha256']:raise ValueError('The reviewed final v3 source closure changed')
    closure=read(source_freeze)
    for path,digest in closure['source_sha256'].items():
        if sha(path)!=digest:raise ValueError('The frozen final executable changed')
    native=Path(config['E4_native_root'])
    if sha(native/'manifest.json')!=config['E4_native_manifest_sha256']:raise ValueError('The selected full688 native authority changed')
    gate_path=native/'results/complete_native_numerical_control_prefix_gate_20261020.json';gate=read(gate_path)
    flags=('passed','complete','all_preselected_arms_complete','all_required_prefix_checks_passed','all_required_controls_passed')
    if not all(gate.get(k) is True for k in flags) or gate['native_manifest_sha256']!=config['E4_native_manifest_sha256'] or len(gate['qualified_conditions'])!=52:
        raise ValueError('All52 recovered E4 conditions require complete native/control/prefix qualification')
    stages={}
    for name,source in config['stages'].items():
        root=Path(source['root']);status=root/'results/analysis/status.json';s=read(status)
        if not s.get('passed') or not s.get('complete') or s['stage']!=name:raise ValueError('All declared complete normalized stages are required')
        if source.get('status_sha256') and sha(status)!=source['status_sha256']:raise ValueError('A pinned prior normalized stage changed')
        stages[name]=dict(root=str(root),status_sha256=sha(status))
    if set(stages)!={'E1','E2','E3','E4'}:raise ValueError('The complete observation phase requires E1--E4')
    e4freeze=read(Path(stages['E4']['root'])/'freeze.json');qualification=e4freeze['native_qualification'].get(str(native.resolve()))
    if not qualification or qualification['native_manifest_sha256']!=config['E4_native_manifest_sha256'] or qualification['workers']!=688:
        raise ValueError('Normalized E4 must bind the selected complete688 recovery authority')
    trigger=config['E4_triggers']
    if sha(trigger['path'])!=trigger['sha256']:raise ValueError('The original52-condition E4 selection changed')
    spec=dict(stages=stages,E4_triggers=trigger,E4_native_gate=dict(path=str(gate_path),sha256=sha(gate_path)))
    spec_path=bundle/'inputs/observation_input_spec.json'
    if spec_path.exists():raise FileExistsError('Do not overwrite the immutable first-phase observation specification')
    spec_path.write_text(json.dumps(spec,indent=2)+'\n')
    subprocess.run([sys.executable,'-B',closure['source'],'prepare-observations','--bundle',str(bundle),'--spec',str(spec_path)],check=True)
    checkpoint=config['common_checkpoint']
    if sha(checkpoint['config_path'])!=checkpoint['config_sha256'] or sha(checkpoint['source_path'])!=checkpoint['source_sha256']:
        raise ValueError('The additive early-common checkpoint closure changed')
    subprocess.run([sys.executable,'-B',checkpoint['source_path'],'verify-final','--config',checkpoint['config_path'],'--config-sha256',checkpoint['config_sha256']],check=True)
    projection=read(config['context_projection_freeze']['path'])
    if sha(config['context_projection_freeze']['path'])!=config['context_projection_freeze']['sha256'] or sha(projection['source'])!=projection['source_sha256']:
        raise ValueError('The descriptor projection executable changed')
    target=bundle/'results/atlas_context_projection'
    subprocess.run([sys.executable,'-B',projection['source'],'--bundle',str(bundle),'--spec',str(spec_path),'--out',str(target)],check=True)
    result=dict(passed=True,complete=True,config_sha256=expected_config_sha256,input_spec_path=str(spec_path),input_spec_sha256=sha(spec_path),
        common_status_sha256=sha(bundle/'results/assembled_observations/common/status.json'),longest_status_sha256=sha(bundle/'results/assembled_observations/longest/status.json'),
        common_checkpoint_equality_sha256=sha(bundle/'results/common_observation_checkpoint_20261020/canonical_equality_receipt.json'),
        atlas_metrics_plan_path=str(target/'metrics_plan.json'),atlas_metrics_plan_sha256=sha(target/'metrics_plan.json'),no_native_integrations=True,no_final_hypothesis_decisions=True)
    path=bundle/'results/observation_context_execution_complete.json';path.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--config-sha256',required=True);a=p.parse_args();print(json.dumps(run(a.config.resolve(),a.config_sha256),indent=2))
