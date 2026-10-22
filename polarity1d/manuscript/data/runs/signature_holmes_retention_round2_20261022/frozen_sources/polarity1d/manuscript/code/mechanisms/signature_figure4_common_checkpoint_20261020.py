"""Additive common-horizon checkpoint using the reviewed final-v3 estimands.

This releases completed E1/E2/E3 observations before E4, without creating the
canonical full-stage directory or relaxing its gate. Final scientific CSVs must
match this checkpoint byte for byte before its early A6 export is rebound.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def load(config_path, expected):
    if sha(config_path) != expected:
        raise ValueError('Frozen checkpoint configuration changed')
    config = read(config_path)
    own = config['execution_source']
    if str(Path(__file__).resolve()) != own['path'] or sha(__file__) != own['sha256']:
        raise ValueError('Execute the actual frozen checkpoint source')
    closure_path = Path(config['source_freeze']['path'])
    if sha(closure_path) != config['source_freeze']['sha256']:
        raise ValueError('Reviewed final-v3 source closure changed')
    closure = read(closure_path)
    for path, digest in closure['source_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Frozen final-v3 executable changed')
    path = Path(closure['source'])
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location('accepted_final_v3', path)
    final = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(final)
    final.verify_source_freeze(Path(config['bundle']))
    return config, final


def project(config_path, expected):
    config, final = load(config_path, expected)
    if set(config['stages']) != {'E1', 'E2', 'E3'}:
        raise ValueError('The checkpoint requires all three completed prior stages')
    datasets = {name: final.load_stage(row['root'], name, row['status_sha256'])
                for name, row in config['stages'].items()}
    e2 = final.select_sigma_one(datasets['E2'])
    e3 = final.select_sigma_one(datasets['E3'])
    # No E4 observation is fabricated: only the common return value is used.
    common, _, references = final.assemble_metrics(e2['metrics'], e3['metrics'], [], [])
    flags = datasets['E1']['flags']
    required = {final.condition_key(row) for row in flags}
    if len(required) != 24:
        raise ValueError('All reviewed E1 forced-state conditions are required')
    common = final.attach_forced_flags(common, flags, required)
    arrays, groups = final.select_observation_inputs(common, e2, e3, {}, longest=False)
    ratios, tolerances = final.comparison_tables(common, arrays, arm='common')
    out = Path(config['output'])
    out.mkdir(exist_ok=False)
    for name, rows in dict(condition_metrics=[dict(row, analysis_arm='common') for row in common],
                           size_ratios=ratios, noise_tolerance=tolerances).items():
        final.write_csv(out / (name + '.csv'), rows)
    outputs = {path.name: sha(path) for path in out.iterdir() if path.is_file()}
    status = dict(passed=True, complete=True, analysis_arm='common', sigma=1., forcing_sigma=1.,
        checkpoint_only=True, canonical_full_stage_release=False,
        condition_groups=len(groups), source_freeze_sha256=sha(config['source_freeze']['path']),
        rule_freeze_sha256=sha(Path(config['bundle']) / 'rule_freeze.json'),
        input_spec_sha256=expected, source_stage_gates={name: data['gates'] for name, data in datasets.items()},
        outputs=outputs, E3_reference_checks=references,
        scope='Complete E1/E2/E3 common observation checkpoint. No E4, S11 or hypothesis gate is released. The canonical full-stage scientific CSVs must be byte-identical; final A6 must bind the canonical status.',
        no_native_integrations=True, no_final_hypothesis_decisions=True)
    final.write_json(out / 'status.json', status)
    receipt = dict(passed=True, complete=True, configuration_sha256=expected,
        checkpoint_status_path=str(out / 'status.json'), checkpoint_status_sha256=sha(out / 'status.json'),
        outputs=outputs, full_stage_equality_pending=True,
        frozen_helpers=['load_stage', 'select_sigma_one', 'assemble_metrics', 'attach_forced_flags',
                        'select_observation_inputs', 'comparison_tables', 'write_csv'])
    final.write_json(out / 'checkpoint_receipt.json', receipt)
    return receipt


def verify_final(config_path, expected):
    config, final = load(config_path, expected)
    out = Path(config['output'])
    early = read(out / 'status.json')
    canonical = Path(config['bundle']) / 'results/assembled_observations/common'
    accepted = read(canonical / 'status.json')
    if not accepted.get('passed') or not accepted.get('complete'):
        raise ValueError('Complete canonical full-stage common observations are required')
    for name, digest in early['outputs'].items():
        if sha(out / name) != digest or accepted['outputs'].get(name) != digest or sha(canonical / name) != digest:
            raise ValueError('Early and full-stage common scientific tables differ: ' + name)
    result = dict(passed=True, complete=True, checkpoint_status_sha256=sha(out / 'status.json'),
        canonical_status_sha256=sha(canonical / 'status.json'), outputs=early['outputs'],
        all_three_scientific_tables_byte_identical=True,
        final_A6_must_bind_canonical_status=True, no_scientific_rule_change=True)
    target = out / 'canonical_equality_receipt.json'
    if target.exists():
        raise FileExistsError('Do not overwrite the accepted equality receipt')
    final.write_json(target, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['project', 'verify-final'])
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--config-sha256', required=True)
    args = parser.parse_args()
    action = project if args.command == 'project' else verify_final
    print(json.dumps(action(args.config.resolve(), args.config_sha256), indent=2))
