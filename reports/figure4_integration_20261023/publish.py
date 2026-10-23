"""Publish the reviewed Figure 4 integration from a session with Git write access.

Default: verify the prepared file hashes and print scope. --publish stages only
listed files and a scoped documentation patch, validates, commits, and pushes.
Run from the repository root with the existing scaling_sandbox Python.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def run(*args, capture=False, **kwargs):
    return subprocess.run(args, cwd=ROOT, check=True, text=True,
                          stdout=subprocess.PIPE if capture else None, **kwargs)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    plan = json.loads((HERE/'publication_plan.json').read_text())
    assert run('git', 'rev-parse', 'HEAD', capture=True).stdout.strip() == plan['base_commit'], 'HEAD changed; review the publication plan again'
    assert not run('git', 'diff', '--cached', '--name-only', capture=True).stdout.strip(), 'Existing staged changes must be handled separately'
    for name, expected in plan['files_sha256'].items():
        assert sha(ROOT/name) == expected, ('File changed since review', name)
    assert sha(HERE/'scoped_documentation.patch') == plan['documentation_patch_sha256']
    run('git', 'apply', '--cached', '--check', str(HERE/'scoped_documentation.patch'))
    print(f"Verified {len(plan['files_sha256'])} explicit files and scoped documentation patch.", flush=True)
    if not args.publish:
        print('Re-run with --publish in a session where .git is writable and GitHub SSH is available.')
        return
    run('git', 'apply', '--cached', str(HERE/'scoped_documentation.patch'))
    names = list(plan['files_sha256'])
    names += [str((HERE/name).relative_to(ROOT)) for name in ('publication_plan.json', 'scoped_documentation.patch')]
    run('git', 'add', '-f', '--pathspec-from-file=-', '--pathspec-file-nul', input='\0'.join(names)+'\0')
    validation = 'manuscript_outline/outline_v9.validation.json'
    run(sys.executable, '-B', 'scripts/verify_manuscript_outline.py',
        'manuscript_outline/outline_v9.json', '--require-tracked', '--report', validation)
    run(sys.executable, '-B', 'manuscript_outline/verify_figure4_round3_20261022.py')
    # The verifiers update only timestamped validation receipts.
    run('git', 'add', '--', validation, 'manuscript_outline/revision_round3_figure4_20261022/verification.json')
    staged = set(run('git', 'diff', '--cached', '--name-only', capture=True).stdout.splitlines())
    assert staged <= set(names) | set(plan['partial_documentation_paths']), 'Unexpected staged paths'
    run('git', 'diff', '--cached', '--check')
    run('git', 'commit', '-m', 'Integrate completed Figure 4 null campaign into manuscript outline')
    commit = run('git', 'rev-parse', 'HEAD', capture=True).stdout.strip()
    upstream = run('git', 'rev-parse', '--abbrev-ref', '--symbolic-full-name', '@{u}', capture=True).stdout.strip()
    remote, branch = upstream.split('/', 1)
    run('git', 'push', remote, f'HEAD:refs/heads/{branch}')
    remote_commit = run('git', 'ls-remote', remote, f'refs/heads/{branch}', capture=True).stdout.split()[0]
    assert commit == remote_commit, (commit, remote_commit)
    result = dict(commit=commit, upstream=upstream, remote_commit=remote_commit, verified=True)
    (HERE/'publication_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
