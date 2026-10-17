"""A*x response control using exact revision-11 preparation checkpoints."""
from pathlib import Path
import argparse, gzip, json, os, shutil, subprocess, sys
from collections import Counter
import numpy as np
import runner as r
CODE = Path(__file__).resolve().parent
CAMPAIGN = CODE / 'pilots' / 'physical_cue_20261002'
BASE = CODE / 'results' / 'records.json.gz'
OriginalCue = r.Cue

class PhysicalCue(OriginalCue):
    def __call__(self, t, x, L):
        # A*x on the right; A*(L-x) on the left. Preserve crossfade timing.
        return super().__call__(t, x, L) * L

r.Cue = PhysicalCue
original_sources = r.sources

def sources():
    hashes = original_sources()
    for p in [Path(__file__), BASE, CAMPAIGN / 'plot.py']:
        hashes[str(p.relative_to(r.ROOT))] = r.sha(p)
    return hashes
r.sources = sources

def submit(action, array=None, dependency=None):
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '
    import shlex
    command += shlex.join([sys.executable, str(Path(__file__).resolve()), action])
    if array:
        command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    opts = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1', '--mem=3G', '--time=02:00:00', '--job-name=polarity-Ax-' + action, '--chdir=' + str(r.ROOT), '--output=' + str(r.HERE / 'logs' / (action + '_%A_%a.out')), '--error=' + str(r.HERE / 'logs' / (action + '_%A_%a.err'))]
    if array:
        opts.append('--array=' + array)
    if dependency:
        opts.append('--dependency=afterany:' + dependency)
    return subprocess.check_output(opts + ['--wrap', command], text=True).strip().split(';')[0]

def baseline():
    with gzip.open(BASE, 'rt') as f:
        return json.load(f)

def launch():
    CAMPAIGN.mkdir(exist_ok=False)
    plot = (CODE / 'plot.py').read_text().replace('HERE=Path(__file__).resolve().parent', 'HERE=Path(' + repr(str(CODE)) + ')\nimport sys\nsys.path.insert(0,str(HERE))')
    plot = plot.replace('Normalized cues: conditioning, release and response times', 'Response cues A*x; original A*x/L preparations')
    (CAMPAIGN / 'plot.py').write_text(plot)
    r.setup()
    old = baseline()
    assert len(old['preparations']) == len(r.tasks()) == 308
    for prep in old['preparations']:
        checkpoint = r.ROOT / prep['checkpoint']
        assert r.sha(checkpoint) == prep['checkpoint_sha256']
        r.write(r.HERE / 'preparations' / f"task_{prep['task_index']:04d}.json", prep)
    jobs = r.response_plan()
    r.running_status('responses_running')
    response = submit('respond', array=f'0-{len(jobs)-1}')
    final = submit('finalize', dependency=response)
    r.write(r.HERE / 'launch.json', dict(stage='responses_running', response_job_id=response, finalize_job_id=final, response_jobs=len(jobs), launched_at=r.now(), cue='right=A*x, left=A*(L-x), crossfade=20', preparations='Exact revision-11 checkpoints; no preparation rerun'))
    print((r.HERE / 'launch.json').read_text())

def finalize():
    r.verify()
    records = r.collect()
    old = baseline()
    assert records['preparations'] == old['preparations']
    assert len(records['cases']) == 308
    assert all(c['status'] in ['success', 'censored', 'ineligible'] for c in records['cases'])
    records['response_cue'] = dict(right='A*x', left='A*(L-x)', amplitude_by_model=r.PROTOCOL['normalization']['amplitude_by_model'], preparation='Exact Ax/L revision-11 states', source_records_sha256=r.sha(BASE))
    out = r.args.results
    with gzip.open(out / 'records.json.gz', 'wt') as f:
        json.dump(records, f, default=r.convert, allow_nan=False)
    r.write(out / 'run_provenance.json', records['provenance'])
    subprocess.run([sys.executable, str(CAMPAIGN / 'plot.py'), '--results', str(out)], check=True)
    before = {c['task_index']:c for c in old['cases']}
    comparisons = [dict(model=c['model'], assay=c['assay'], length=c['length'], previous_status=before[c['task_index']]['status'], status=c['status'], previous_endpoint=before[c['task_index']].get('endpoint'), endpoint=c.get('endpoint')) for c in records['cases']]
    r.csvwrite(out / 'comparison_to_normalized.csv', comparisons)
    r.write(out / 'validation.json', dict(generated_at=r.now(), passed=True, exact_preparation_records=True, checkpoint_hashes_verified=True, complete_duration_coverage=True, cases=len(records['cases']), finite_trials=sum(len(c.get('tests', [])) for c in records['cases']), source_records_sha256=r.sha(BASE)))
    subprocess.run([sys.executable, '-c', 'import sys; from pathlib import Path; sys.path.insert(0,sys.argv[1]); import validate_regressions as v; v.r.CODE=Path(sys.argv[2]); v.main()', str(CODE), str(CAMPAIGN)], check=True)
    r.write(out / 'status.json', dict(generated_at=r.now(), complete=True, status_counts=dict(Counter(c['status'] for c in records['cases'])), cue='A*x', preparations='Exact original Ax/L states', validation_passed=True))
    stamp = r.now()
    (out / 'index.html').write_text(f'<html><body><h1>A*x response control</h1><p>Generated {stamp}; producer finished. Exact original A*x/L preparation states.</p><a href="length_response_and_settling.pdf">PDF</a><img width="1400" src="length_response_and_settling.png"></body></html>')
    launch = json.loads((r.HERE / 'launch.json').read_text())
    launch.update(stage='finished', finished_at=stamp)
    r.write(r.HERE / 'launch.json', launch)
    print(json.dumps(launch))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=['launch', 'respond', 'finalize', 'status'])
    ap.add_argument('--slot', type=int)
    args = ap.parse_args()
    r.HERE = CAMPAIGN / 'work'
    r.args = argparse.Namespace(results=CAMPAIGN / 'results')
    if args.action == 'launch':
        launch()
    elif args.action == 'respond':
        assert r.sha(r.HERE / 'response_manifest.json') == json.loads((r.HERE / 'response_manifest_sha256.json').read_text())['sha256']
        r.run_response(args.slot)
    elif args.action == 'finalize':
        finalize()
    else:
        print((r.HERE / 'launch.json').read_text())
