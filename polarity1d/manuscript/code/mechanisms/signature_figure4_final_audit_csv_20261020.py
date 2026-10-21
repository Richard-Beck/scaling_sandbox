"""Run the unchanged Figure 4 v3 auditor with its accepted CSV capacity.

The published H3 pair field is 145,670 bytes, exceeding Python's default
131,072-byte CSV field limit. The accepted execution recipe raises that limit
to 64 MiB before running the original hash-checked scientific executable.
It changes no audit algorithm, numerical tolerance, input or decision rule.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy
import sys

AUDITOR_NAME = 'signature_figure4_final_audit_20261020.py'
AUDITOR_SHA256 = 'ba2c19889f0d18b6711636200e1a84f8e1c1e8a9eafc73a7e32010641944de64'
CSV_FIELD_LIMIT = 64 * 1024 * 1024


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True,
                        help='Restored final-statistics bundle with its bound inputs and frozen v3 sources.')
    parser.add_argument('--audit-source', type=Path,
                        help='Frozen auditor; defaults to BUNDLE/frozen_sources_v3/' + AUDITOR_NAME)
    parser.add_argument('--check', action='store_true',
                        help='Check the source/path binding and parser capacity without running the audit or writing a receipt.')
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    source = (args.audit_source or bundle / 'frozen_sources_v3' / AUDITOR_NAME).resolve()
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != AUDITOR_SHA256:
        parser.error('The executable differs from the accepted unchanged v3 auditor.')
    complete = json.loads((bundle / 'results/final_analysis/complete.json').read_text())
    if complete['source_sha256'].get(str(source)) != AUDITOR_SHA256:
        parser.error('Use the actual frozen auditor path bound by the producer completion receipt.')
    csv.field_size_limit(CSV_FIELD_LIMIT)
    if args.check:
        print(json.dumps(dict(passed=True,check_only=True,auditor_path=str(source),
            auditor_sha256=digest,csv_field_size_limit=csv.field_size_limit(),
            scientific_auditor_unchanged=True,receipt_written=False),indent=2))
        return
    sys.argv = [str(source), '--bundle', str(bundle)]
    runpy.run_path(str(source), run_name='__main__')


if __name__ == '__main__':
    main()
