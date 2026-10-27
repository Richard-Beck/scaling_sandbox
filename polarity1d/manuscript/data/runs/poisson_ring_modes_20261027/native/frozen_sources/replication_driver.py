"""Prepare bounded F4 replication with an explicit measurement-scope decision.

The original scientific engine, strict payload and failed Spring receipt remain
unchanged. Preparing/checking runs no integration and grants no scope acceptance.
Submission requires separate root acceptance and independent exact-manifest
preflight for this new payload. There is no automatic submission or resubmission.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import copy
import hashlib
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
sys.dont_write_bytecode = True

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / 'polarity1d/manuscript').is_dir())
MS = ROOT / 'polarity1d/manuscript'
STRICT = MS / 'data/runs/revision_noise_strength_replication_20261019'
DEFAULT = MS / 'data/runs/revision_noise_strength_replication_scoped_20261019'
STRICT_MANIFEST = '0c6820ee06f3fb8f4f182b9a91a60947be7ba065a25908909776ca520bd9f27c'
CORE_SHA = 'f9a5eb848453d25f9992b3358cd94b07cc8103338a95bd0747be1d2f64f90a41'
PILOT_MANIFEST = '795023e22468754dd481e00d09a97bd6d371b15f88b8ea0b75c3a4c55e9c4dce'
SPRING_RECEIPT = '90acb9d1740acaa698a22d12e0d5c385dbe09922c4cc8c60e010e14f6a2a6eb4'
PYTHON = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'
STABLE = ('mean_R', 'mean_A', 'mean_Q', 'mean_contrast_ratio',
          'multiregion_fraction', 'qualified_multiregion_fraction',
          'weak_R_fraction', 'low_contrast_fraction',
          'direction_correlation', 'direction_eligible_fraction')
EXCLUDED = ('multiregion_run_q90', 'qualified_multiregion_run_q90')


def now():
    return datetime.now(timezone.utc).isoformat()


def read(p):
    return json.loads(Path(p).read_text())


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def write(p, value):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    temporary = p.with_suffix(p.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(p)


def core(run):
    path = Path(run) / 'frozen_sources/replication_core.py'
    assert sha(path) == CORE_SHA
    spec = importlib.util.spec_from_file_location('unchanged_replication_core', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PLAN = '''# Bounded common-strength replication with explicit measurement scope

This fresh payload preserves all 360 proposed native jobs from the preceding
strict payload: f=0.003,0.03,0.1, four models, three original sizes, independent
replicate indices2..11. All576 earlier scientific cases and their exact native
hashes remain retained. The target is twelve seeds at these three levels;
original0.001/0.01 already have twelve. Other broad endpoints remain n=2.

Native functions, positive whole-signal lognormal scaling, declared ensemble
means, sigma, correlations, finite forcing knots, basal inputs, solver settings,
horizons and saved cadences are unchanged. The old scientific engine and strict
orchestration core are copied byte-for-byte. There is no new variance-only arm,
range extension, dense companion, mean control or no-input control.

The original selected Spring C0=3,f=1 numerical/observation gate is FAILED:
12 coarse/dense run-q90 comparisons fail;186 other original comparisons pass.
Its original receipt and SHA are immutable. Standard/tight canonical ODE
estimates agree. The failure demonstrates observation-cadence sensitivity of
sampled raw/qualified multiregion run spans on that tested path. It does not
demonstrate a lower-strength failure, nor justify continuous lifetime inference.

This protocol proposes a separate qualification of the ten measurements listed
in STABLE_ESTIMAND_SCOPE.json. For the tested Spring path all180 rows in that
scope pass; all18 run-q90 rows are excluded from that scope, including six ODE
passes. Goryachev's twelve stable-scope unavailable-both comparisons remain
unavailable. The four selected upper-strength tests support this chosen
measurement protocol; they do not establish convergence for every condition.

Primary qualified comparisons cover physical first-moment magnitude R, active
amount ratio A, coherence ratio Q, contrast ratio, raw/contrast-qualified
multiple-region occupancy, weak-R and low-contrast occupancy, and directional
correlation WITH eligible-pair coverage. R=AQ is pointwise; the mean of the
product differs from the product of means. Preserve native physical M0/M1,
actual Spring C(t), basal values, per-seed descriptors and full/first/last-half
windows, including the shared midpoint. Bootstrap whole independent seeds.

All sampled run spans remain explicit cadence-dependent diagnostics with edge
censoring; no continuous lifetimes, event rates or cadence-invariant Spring
duration comparisons are qualified here. q2, peak/prominence sensitivities,
raw-vector correlation, site summaries and selected profiles remain descriptive
evidence; this narrowed acceptance does not claim that the original upper-f
gate directly qualified every such exported quantity. Counts do not establish
competition, site lineage, exclusive phases or a calibrated cross-model input.

Preparation and implementation preflight do not accept this scientific scope.
Submission requires a separate named root decision accepting exactly the scope
hash, plus independent implementation/scope preflight for this exact new
manifest. The failed original receipt remains passed:false. The earlier strict
payload remains blocked and is never overwritten or used for submission here.

Resources remain QoSsmall,1CPU,3GiB,2h per worker, no throttle. Measured first
stage cases predict63.1CPUh and3.83GB additional allocated native output;
30% storage reserve4.98GB. These are estimates, not upper bounds. Preserve
finite/pool/positive-circumference/source guards in every worker. No automatic
resubmission is implemented. No job is submitted by prepare or check.
'''


def prepare(run):
    run = Path(run).resolve()
    assert not run.exists() and run.is_relative_to(MS / 'data/runs')
    assert sha(STRICT / 'manifest.json') == STRICT_MANIFEST
    old = read(STRICT / 'manifest.json')
    for sub in ('inputs', 'frozen_sources', 'protocol', 'workers', 'logs', 'results'):
        (run / sub).mkdir(parents=True)
    sources = {}
    for relative, digest in old['source_sha256'].items():
        source = STRICT / relative
        assert sha(source) == digest
        target_relative = ('frozen_sources/replication_core.py'
                           if relative == 'frozen_sources/replication_driver.py'
                           else relative)
        target = run / target_relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        sources[target_relative] = digest
    shutil.copy2(__file__, run / 'frozen_sources/replication_driver.py')
    sources['frozen_sources/replication_driver.py'] = sha(__file__)
    assert sources['frozen_sources/replication_core.py'] == CORE_SHA
    for b in old['basals']:
        assert sha(STRICT / b['input']) == b['input_sha256']
        shutil.copy2(STRICT / b['input'], run / b['input'])
    shutil.copy2(STRICT / 'protocol/PLAN.md', run / 'protocol/PRIOR_STRICT_PLAN.md')
    shutil.copy2(STRICT / 'protocol/PILOT_PLAN.md', run / 'protocol/PILOT_PLAN.md')
    shutil.copy2(STRICT / old['reuse']['index'], run / old['reuse']['index'])
    (run / 'protocol/PLAN.md').write_text(PLAN)
    receipts = {}
    for model, path in old['numerical_acceptance_required']['receipts'].items():
        data = read(path)
        assert data['manifest_sha256'] == PILOT_MANIFEST and data['model'] == model
        receipts[model] = dict(path=path, sha256=sha(path), original_passed=data['passed'])
    assert receipts['spring_ring']['sha256'] == SPRING_RECEIPT
    assert receipts['spring_ring']['original_passed'] is False
    assert all(v['original_passed'] is True for k, v in receipts.items() if k != 'spring_ring')
    scope = dict(scope_id='f4_broad_strength_stable_estimands_v1',
        pilot_manifest_sha256=PILOT_MANIFEST,
        original_full_gate_spring_passed=False,
        original_spring_receipt_sha256=SPRING_RECEIPT,
        primary_estimands=list(STABLE), excluded_from_qualified_scope=list(EXCLUDED),
        excluded_inference=['continuous lifetimes', 'event rates',
                            'cadence-invariant Spring run-span comparisons'],
        windows=['full', 'first_half', 'last_half'], midpoint_shared=True,
        direction_endpoint_gate='both R>=0.1',
        direction_group_reporting_gate='eligible pairs>=100 and eligible fraction>=0.01; otherwise unavailable',
        chemical_scalar_cadence=5., spring_scalar_cadence=.5,
        chemical_profile_cadence=40., spring_profile_cadence=4.,
        tested_scope='Four selected largest-size f=1 paths, standard/tight canonical and coarse/dense saved observations',
        duration_limit='Observed failure is Spring C0=3,f=1 selected seed only; all-strength spans stay sampled-cadence descriptors without continuous-lifetime claims',
        supplemental_readouts='q2, prominence sensitivity, raw-vector correlation, site summaries and profiles remain descriptive; not directly covered by this narrowed gate',
        scientific_sources_changed=False, forcing_law_changed=False,
        solver_changed=False, observation_cadence_changed=False,
        original_receipts=receipts)
    write(run / 'protocol/STABLE_ESTIMAND_SCOPE.json', scope)
    manifest = copy.deepcopy(old)
    manifest.update(generated_at=now(), experiment='Bounded common-strength replication with explicit measurement scope',
        status='frozen_scope_acceptance_pending_not_submitted', source_sha256=sources,
        protocol_sha256=sha(run / 'protocol/PLAN.md'),
        previous_strict_payload=dict(run=str(STRICT), manifest_sha256=STRICT_MANIFEST,
            immutable=True, all_four_original_pass_guard_retained=True),
        scoped_qualification=dict(scope='protocol/STABLE_ESTIMAND_SCOPE.json',
            scope_sha256=sha(run / 'protocol/STABLE_ESTIMAND_SCOPE.json'),
            original_receipts=receipts,
            acceptance_required='Separate named root acceptance plus exact-manifest independent scoped preflight'))
    manifest['plan']['measurement_scope_changed'] = True
    manifest['plan']['native_jobs_changed_from_strict_payload'] = False
    write(run / 'manifest.json', manifest)
    template = dict(accepted_for_stable_estimands=False,
        review_role='root_scientific_scope_acceptance', accepted_by=None, accepted_at=None,
        manifest_sha256=sha(run / 'manifest.json'),
        stable_estimand_scope_sha256=manifest['scoped_qualification']['scope_sha256'],
        original_spring_receipt_sha256=SPRING_RECEIPT,
        original_full_gate_passed=False, primary_estimands=list(STABLE),
        excluded_from_qualified_scope=list(EXCLUDED),
        rationale='Pending root decision; preparing this template is not acceptance')
    write(run / 'results/STABLE_ESTIMAND_ACCEPTANCE_TEMPLATE.json', template)
    write(run / 'status.json', dict(generated_at=now(), producer='frozen_not_submitted',
        native_jobs_submitted=0, requested_jobs=360, retained_scientific_cases=576,
        scope_acceptance='pending', complete=False))
    return dict(run=str(run), manifest_sha256=sha(run / 'manifest.json'),
                scope_sha256=manifest['scoped_qualification']['scope_sha256'], jobs=360, submitted=False)


def qualification(run, acceptance=None):
    run = Path(run).resolve()
    m = read(run / 'manifest.json')
    q = m['scoped_qualification']
    assert sha(run / q['scope']) == q['scope_sha256']
    scope = read(run / q['scope'])
    assert scope['primary_estimands'] == list(STABLE) and scope['excluded_from_qualified_scope'] == list(EXCLUDED)
    rows = []
    for model, r in q['original_receipts'].items():
        assert sha(r['path']) == r['sha256']
        d = read(r['path'])
        assert d['manifest_sha256'] == PILOT_MANIFEST and d['model'] == model
        assert d['passed'] is r['original_passed']
        gates = d['gate_values']
        selected = [g for g in gates if g['metric'] in STABLE]
        excluded = [g for g in gates if g['metric'] in EXCLUDED]
        assert len(gates) == 198 and len(selected) == 180 and len(excluded) == 18
        assert all(g['metric'] in STABLE + EXCLUDED for g in gates)
        assert all(g['status'] in ('pass', 'unavailable_both') for g in selected)
        if model == 'spring_ring':
            assert r['sha256'] == SPRING_RECEIPT and d['passed'] is False
            assert sum(g['status'] == 'fail' for g in excluded) == 12
        rows.append(dict(model=model, original_whole_gate_passed=d['passed'],
            selected_rows=len(selected), selected_passes=sum(g['status'] == 'pass' for g in selected),
            selected_unavailable_both=sum(g['status'] == 'unavailable_both' for g in selected),
            excluded_rows=len(excluded), excluded_failures=sum(g['status'] == 'fail' for g in excluded),
            receipt=r))
    acceptance_path = Path(acceptance).resolve() if acceptance else run / 'results/STABLE_ESTIMAND_ACCEPTANCE.json'
    a = read(acceptance_path) if acceptance_path.exists() else {}
    accepted = (a.get('accepted_for_stable_estimands') is True
        and a.get('review_role') == 'root_scientific_scope_acceptance'
        and isinstance(a.get('accepted_by'), str) and bool(a['accepted_by'].strip())
        and isinstance(a.get('accepted_at'), str) and bool(a['accepted_at'].strip())
        and a.get('manifest_sha256') == sha(run / 'manifest.json')
        and a.get('stable_estimand_scope_sha256') == q['scope_sha256']
        and a.get('original_spring_receipt_sha256') == SPRING_RECEIPT
        and a.get('original_full_gate_passed') is False
        and a.get('primary_estimands') == list(STABLE)
        and a.get('excluded_from_qualified_scope') == list(EXCLUDED))
    return dict(scoped_numeric_evidence_reconciles=True, original_full_spring_gate_passed=False,
        scope_sha256=q['scope_sha256'], model_evidence=rows,
        scope_accepted=accepted, acceptance_path=str(acceptance_path),
        acceptance_sha256=sha(acceptance_path) if acceptance_path.exists() else None,
        ready_for_scoped_submission=accepted)


def check(run, acceptance=None):
    run = Path(run).resolve()
    m = read(run / 'manifest.json')
    strict = read(STRICT / 'manifest.json')
    assert sha(STRICT / 'manifest.json') == STRICT_MANIFEST
    for key in ('jobs', 'basals', 'native_protocols', 'resources', 'reuse', 'numerical_check'):
        assert m[key] == strict[key], f'Native/reuse change from strict payload: {key}'
    result = core(run).check(run)
    result['original_all_four_numerical_qualification'] = result.pop('numerical_qualification')
    result['scoped_qualification'] = qualification(run, acceptance)
    result.update(native_jobs_identical_to_strict=True, ready_to_submit=False,
        reason='Implementation check only; separate root scope acceptance and exact-manifest independent scoped preflight remain required')
    write(run / 'results/local_preflight.json', result)
    return result


def submit(run, review, acceptance):
    run = Path(run).resolve()
    m = read(run / 'manifest.json')
    r = read(review)
    assert r.get('passed') is True and r.get('manifest_sha256') == sha(run / 'manifest.json'), 'Independent exact-manifest preflight required'
    assert r.get('scoped_qualification_reviewed') is True and r.get('stable_estimand_scope_sha256') == m['scoped_qualification']['scope_sha256'], 'Independent scope review required'
    q = qualification(run, acceptance)
    assert q['ready_for_scoped_submission'], 'Separate named root measurement-scope acceptance required'
    local = read(run / 'results/local_preflight.json')
    assert local['passed'] is True and local['manifest_sha256'] == sha(run / 'manifest.json')
    assert not (run / 'submission.json').exists(), 'No automatic resubmission'
    assert not any((run / 'workers').iterdir()), 'Nonempty worker root requires diagnosis'
    for relative, digest in m['source_sha256'].items():
        assert sha(run / relative) == digest
    assert sha(run / 'protocol/PLAN.md') == m['protocol_sha256']
    args = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
        '--mem=3G', '--time=02:00:00', '--signal=B:USR1@60', '--array=0-359',
        '--job-name=f4-strength-repl', '--output=' + str(run / 'logs/%A_%a.out'),
        '--error=' + str(run / 'logs/%A_%a.err'), '--wrap', 'exec ' + shlex.join([
        PYTHON, str(run / 'frozen_sources/replication_driver.py'), 'worker', '--run', str(run)])
        + ' --slot "$SLURM_ARRAY_TASK_ID"']
    handle = subprocess.check_output(args, text=True).strip().split(';')[0]
    receipt = dict(submitted_at=now(), job_id=handle, manifest_sha256=sha(run / 'manifest.json'),
        preflight_review=str(Path(review).resolve()), preflight_sha256=sha(review),
        scoped_qualification=q, original_full_spring_gate_passed=False, command=args)
    write(run / 'submission.json', receipt)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('prepare', 'check', 'worker', 'submit'))
    p.add_argument('--run', type=Path, default=DEFAULT)
    p.add_argument('--slot', type=int)
    p.add_argument('--review', type=Path)
    p.add_argument('--scope-acceptance', type=Path)
    a = p.parse_args()
    if a.command == 'prepare':
        result = prepare(a.run)
    elif a.command == 'check':
        result = check(a.run, a.scope_acceptance)
    elif a.command == 'worker':
        result = core(a.run).engine(a.run).worker(a.run, a.slot)
    else:
        assert a.review is not None and a.scope_acceptance is not None
        result = submit(a.run, a.review, a.scope_acceptance)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
