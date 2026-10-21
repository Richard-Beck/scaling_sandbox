"""Release recovered E4 derived work only after its full numerical authority."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def launch(controller, config_path):
    controller = Path(controller).resolve()
    config_path = Path(config_path).resolve()
    freeze = read(controller / 'freeze.json')
    for mapping in ('source_sha256', 'plan_sha256', 'script_sha256'):
        for path, digest in freeze[mapping].items():
            if sha(path) != digest:
                raise ValueError('The frozen recovered-intake closure changed')
    if str(Path(__file__).resolve()) not in freeze['source_sha256'] or str(config_path) not in freeze['plan_sha256']:
        raise ValueError('Execute the actual frozen selected-intake source and configuration')
    config = read(config_path)
    native = Path(config['native_root'])
    if sha(native / 'manifest.json') != config['native_manifest_sha256']:
        raise ValueError('The selected full688 native authority changed')
    gate_path = native / 'results/complete_native_numerical_control_prefix_gate_20261020.json'
    gate = read(gate_path)
    flags = ('passed', 'complete', 'all_preselected_arms_complete',
             'all_required_prefix_checks_passed', 'all_required_controls_passed',
             'numerical_integrity_only')
    if (not all(gate.get(name) is True for name in flags) or
            gate.get('scientific_long_arm_qualification') is not False or
            gate['native_manifest_sha256'] != config['native_manifest_sha256'] or
            len(gate['qualified_conditions']) != 52):
        raise ValueError('All52 selected conditions require the complete numerical/control/prefix gate')
    manifest = read(native / 'manifest.json')
    if manifest['stage'] != 'E4' or len(manifest['jobs']) != 688:
        raise ValueError('The complete scientific E4 cohort remains exactly688 slots')
    plan_path = controller / 'inputs/E4_plan.json'
    plan = read(plan_path)
    if plan['native_roots'] != [dict(path=str(native), manifest_sha256=config['native_manifest_sha256'])]:
        raise ValueError('The normalized plan must bind this exact selected native authority')
    receipt = dict(passed=True, complete=True, native_manifest_sha256=config['native_manifest_sha256'],
        numerical_projection_path=str(gate_path), numerical_projection_sha256=sha(gate_path),
        no_scientific_long_arm_qualification=True, all_required_native_controls_and_prefixes_present=True)
    target = controller / 'results/selected_native_release.json'
    if target.exists():
        raise FileExistsError('Do not repeat an accepted recovered-intake release')
    target.write_text(json.dumps(receipt, indent=2) + '\n')
    subprocess.run([sys.executable, '-B', config['stage_execution_source'],
                    '--controller', str(controller), '--plan', str(plan_path)], check=True)
    subprocess.run([sys.executable, '-B', config['followup_source'],
                    '--controller', str(controller)], check=True)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--controller', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(launch(args.controller, args.config), indent=2))
