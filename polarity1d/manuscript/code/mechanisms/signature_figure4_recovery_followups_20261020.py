"""Queue independent derived phases after the actual recovered E4 collector."""
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


def queue(controller):
    controller = Path(controller).resolve()
    freeze = read(controller / 'freeze.json')
    for mapping in ('source_sha256', 'plan_sha256', 'script_sha256'):
        for path, digest in freeze[mapping].items():
            if sha(path) != digest:
                raise ValueError('The frozen derived follow-up closure changed')
    if str(Path(__file__).resolve()) not in freeze['source_sha256']:
        raise ValueError('Execute the actual frozen follow-up source')
    config = read(controller / 'inputs/followups.json')
    pipeline_path = Path(config['parent']) / 'results/pipeline_submission.json'
    pipeline = read(pipeline_path)
    plan_path = controller / 'inputs/E4_plan.json'
    if (pipeline['stage'] != 'E4' or not pipeline['all_arms_submitted'] or
            pipeline['plan_sha256'] != sha(plan_path) or
            pipeline['inventory_sha256'] != sha(Path(config['parent']) / 'inventory.json')):
        raise ValueError('The actual complete recovered normalization pipeline is required')
    collector = pipeline['final_collect_job_id']
    if not str(collector).isdigit():
        raise ValueError('The actual final collector must have an unambiguous job ID')
    target = controller / 'results/followup_submission.json'
    if target.exists():
        raise FileExistsError('Do not duplicate derived follow-up submissions')
    submissions = []
    for role in ('original_seeds', 'observations'):
        command = ['sbatch', '--parsable', '--dependency=afterok:' + str(collector),
                   '--job-name=F4E4-' + role, config['scripts'][role]]
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        job = result.stdout.strip().split(';')[0]
        if not job.isdigit():
            raise ValueError('SLURM must return an unambiguous job ID')
        submissions.append(dict(role=role, job_id=job, command=command))
        target.write_text(json.dumps(dict(parent_final_collect_job_id=collector,
            parent_pipeline_sha256=sha(pipeline_path), submissions=submissions,
            all_followups_submitted=False), indent=2) + '\n')
    subprocess.run([sys.executable, '-B', config['site_execution_source'], 'chain',
                    '--controller', str(controller), '--plan', config['site_plan']], check=True)
    site = read(controller / 'results/E4_chain.json')
    if site['parent_final_collect_job_id'] != collector:
        raise ValueError('Literal-site follow-up must depend on the same actual collector')
    submissions.append(dict(role='literal_sites_v2', **site['submission']))
    result = dict(passed=True, all_followups_submitted=True,
        parent_final_collect_job_id=collector, parent_pipeline_sha256=sha(pipeline_path),
        execution_freeze_sha256=sha(controller / 'freeze.json'), submissions=submissions,
        no_native_integrations=True, no_scientific_gate_relaxed=True)
    target.write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--controller', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(queue(args.controller), indent=2))
