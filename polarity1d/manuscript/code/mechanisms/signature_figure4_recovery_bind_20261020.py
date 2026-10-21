"""Bind prepared E4 derived sources to the exact selected native authority."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

BASE = Path('/home/4473331/projects/scaling_sandbox/polarity1d/manuscript')
RUNS = BASE / 'data/runs'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def bind(expected_manifest, projection_job):
    if not projection_job.isdigit():
        raise ValueError('An exact selected-native projection job ID is required')
    controller = RUNS / 'signature_figure4_recovered_stage_execution_20261020'
    native = RUNS / 'signature_figure4_E4_selected_native_20261020'
    bundle = RUNS / 'signature_figure4_final_statistics_20261020'
    sources = controller / 'frozen_sources'
    if (controller / 'freeze.json').exists():
        raise FileExistsError('Do not rebind or resubmit an immutable recovered controller')
    for path, digest in read(controller / 'source_closure.json')['source_sha256'].items():
        if sha(path) != digest:
            raise ValueError('A prepared reviewed source changed')
    if sha(native / 'manifest.json') != expected_manifest:
        raise ValueError('The declared selected native manifest changed')
    manifest = read(native / 'manifest.json')
    original = read(RUNS / 'signature_figure4_E4_horizon_20261020/manifest.json')
    if manifest['stage'] != 'E4' or len(manifest['jobs']) != 688:
        raise ValueError('The scientific E4 cohort must remain exactly688 slots')
    for selected, prior in zip(manifest['jobs'], original['jobs']):
        if any(selected[key] != prior[key] for key in ('slot', 'model', 'size', 'fraction', 'kind', 'seed', 'sigma', 'protocol')):
            raise ValueError('Source recovery cannot change a scientific native job')
    old_basals = {(row['model'], row['size']): row['input_sha256'] for row in original['basals']}
    if {(row['model'], row['size']): row['input_sha256'] for row in manifest['basals']} != old_basals:
        raise ValueError('Source recovery cannot change a native basal')
    gate = native / 'inputs/e0_interim_gate.json'
    if sha(gate) != '0948e9da7155511056f306fdec406377bd14209df61c79395fae58a9930bb1c5':
        raise ValueError('The delivered prior-stage release gate changed')
    shutil.copy2(Path(__file__), sources / Path(__file__).name)
    inputs = controller / 'inputs'
    if any(inputs.iterdir()):
        raise FileExistsError('Bind only the untouched prepared input directory')
    old_controller = RUNS / 'signature_figure4_late_stage_execution_20261020'
    plan = read(old_controller / 'inputs/E4_plan.json')
    plan.update(native_roots=[dict(path=str(native), manifest_sha256=expected_manifest)],
                gate=dict(path=str(gate), sha256=sha(gate)),
                required_native_master=str(native / 'results/complete_native_independent_audit_20261020.json'))
    write(inputs / 'E4_plan.json', plan)
    shutil.copy2(old_controller / 'inputs/E4_original_seed_plan.json', inputs / 'E4_original_seed_plan.json')
    site_rule = BASE.parent.parent / 'docs/tasks/briefs/20261020_figure4_literal_site_coordinate_addendum_v2.json'
    if sha(site_rule) != '543b78bb11a823bb30791627312e5292c8d7f5c8987ce5cf1af4e9292eb53399':
        raise ValueError('The reviewed literal-v2 route policy changed')
    write(inputs / 'site_E4_plan.json', dict(stage='E4', parent=plan['out'],
        out=str(RUNS / 'signature_figure4_literal_sites_v2_E4_20261020'),
        parent_execution_plan_sha256=sha(inputs / 'E4_plan.json'), rule_path=str(site_rule), rule_sha256=sha(site_rule)))
    early = read(bundle / 'frozen_common_checkpoint_20261020/submission.json')
    observation_source = sources / 'signature_figure4_observation_execution_20261020.py'
    stages = {}
    pins = {'E1': '4e94fa8d1d837ac635b62e1bd1cd75d8295dbd45f4213f74522e34c21c9ffd60',
            'E2': '450836b83b631562c6da7e67a8b261ed5b3d1e3e3705330e128fb043467d1ac6',
            'E3': '401fef2abaa8941a3ff1bffe7d9d19f672759792cd3f19c595ad6e9936bf75ab'}
    for stage in ('E1', 'E2', 'E3', 'E4'):
        row = dict(root=str(RUNS / f'signature_figure4_{stage}_analysis_20261020'))
        if stage in pins:
            if sha(Path(row['root']) / 'results/analysis/status.json') != pins[stage]:
                raise ValueError('An accepted normalized prior stage changed')
            row['status_sha256'] = pins[stage]
        stages[stage] = row
    config = dict(execution_source=dict(path=str(observation_source), sha256=sha(observation_source)),
        bundle=str(bundle), source_freeze=dict(path=str(bundle / 'source_freeze_v3.json'),
            sha256='f37ffb5144144935c1c71e5f0255ec96a0742d143a4f582293be21832cc0791a'),
        E4_native_root=str(native), E4_native_manifest_sha256=expected_manifest, stages=stages,
        E4_triggers=dict(path=str(RUNS / 'signature_figure4_e0_corrected_basals_20261020/results/interim_handoff/e4_trigger_inventory.csv'),
            sha256='fc20b3a8cb68f30e00cabc41c31f76076aa964c453258c3d0f6c90ca2024c8f5'),
        context_projection_freeze=dict(path=str(bundle / 'context_projection_freeze.json'),
            sha256='0ccba8275669e8d621b200c4fc250ad4f1325f305292db6cf00015898ba73025'),
        common_checkpoint={key: early[key] for key in ('config_path', 'config_sha256', 'source_path', 'source_sha256')})
    write(inputs / 'observation_config.json', config)
    original_script = controller / 'launch_original_seeds.slurm'
    observation_script = controller / 'launch_observations.slurm'
    stage_source = sources / 'signature_figure4_stage_execution_20261020.py'
    followup_source = sources / 'signature_figure4_recovery_followups_20261020.py'
    write(inputs / 'followups.json', dict(parent=plan['out'],
        scripts=dict(original_seeds=str(original_script), observations=str(observation_script)),
        site_execution_source=str(sources / 'signature_figure4_site_execution_20261020.py'),
        site_plan=str(inputs / 'site_E4_plan.json')))
    write(inputs / 'intake_config.json', dict(native_root=str(native), native_manifest_sha256=expected_manifest,
        stage_execution_source=str(stage_source), followup_source=str(followup_source)))
    q = shlex.quote
    python = q(sys.executable)
    def script(path, command):
        path.write_text('\n'.join(['#!/bin/bash', '#SBATCH --partition=red', '#SBATCH --qos=small',
            '#SBATCH --cpus-per-task=1', '#SBATCH --mem=8G', '#SBATCH --time=01:00:00',
            f'#SBATCH --output={controller}/logs/{path.stem}_%j.log', 'set -euo pipefail',
            'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1', command]) + '\n')
    script(original_script, f'{python} -B {q(str(stage_source))} --controller {q(str(controller))} --plan {q(str(inputs / "E4_original_seed_plan.json"))}')
    script(observation_script, f'{python} -B {q(str(observation_source))} --config {q(str(inputs / "observation_config.json"))} --config-sha256 {sha(inputs / "observation_config.json")}')
    main = controller / 'launch_selected_normalization.slurm'
    script(main, f'{python} -B {q(str(sources / "signature_figure4_selected_intake_20261020.py"))} --controller {q(str(controller))} --config {q(str(inputs / "intake_config.json"))}')
    freeze = dict(generated_utc=datetime.now(timezone.utc).isoformat(), native_manifest_sha256=expected_manifest,
        source_sha256={str(path): sha(path) for path in sorted(sources.glob('*.py'))},
        plan_sha256={str(path): sha(path) for path in sorted(inputs.glob('*.json'))},
        script_sha256={str(path): sha(path) for path in (main, original_script, observation_script)},
        source_closure_sha256=sha(controller / 'source_closure.json'),
        scientific_stage_and_literal_sources_unchanged=True, actual_final_collector_required_for_followups=True)
    write(controller / 'freeze.json', freeze)
    command = ['sbatch', '--parsable', '--dependency=afterok:' + projection_job,
               '--job-name=F4E4-selected-normalize', str(main)]
    job = subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip().split(';')[0]
    if not job.isdigit():
        raise ValueError('SLURM must return an unambiguous normalized-intake job ID')
    receipt = dict(passed=True, job_id=job, command=command, projection_dependency=projection_job,
        native_manifest_sha256=expected_manifest, execution_freeze_sha256=sha(controller / 'freeze.json'),
        all_followups_queued_by_actual_normalized_collector=True, no_native_integrations=True)
    write(controller / 'results/submission.json', receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--projection-job', required=True)
    args = parser.parse_args()
    print(json.dumps(bind(args.manifest_sha256, args.projection_job), indent=2))
