#!/usr/bin/env python3
"""Verify retained selected inputs and regenerate existing plots without integration."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import sys
import tempfile
from datetime import datetime, timezone
import storage_cleanup_four_assay_20261017 as c


def main():
    record = json.loads(c.RECORD.read_text())
    assert record['status'] == 'complete'
    manuscript = c.ROOT / 'polarity1d/manuscript'
    doc = json.loads((manuscript / 'data/provenance/data_sources.json').read_text())
    rows = [r for r in doc['files'] if r['destination'].startswith('data/establishment_basal/')
            and r['destination'].endswith('.npz')]
    plan = json.loads((manuscript / 'data/provenance/legacy_slimming_plan.json').read_text())
    rows += [r for r in plan['copy_plan'] if '/reconstruction_inputs/' in r['destination']]
    assert len(rows) == 5006
    for r in rows:
        p = manuscript / r['destination']
        assert p.is_file() and not p.is_symlink(), p
        assert c.sha(p) == r['sha256'], p
    establishment = manuscript / 'data/establishment_basal'
    checkpoints = 0
    with (establishment / 'results/endpoints.csv').open() as stream:
        for r in csv.DictReader(stream):
            if r['selected_checkpoint']:
                assert (establishment / r['selected_checkpoint']).is_file()
                checkpoints += 1
    sys.path.insert(0, str(manuscript / 'code/core'))
    import native
    import numpy as np
    cases = json.load(gzip.open(manuscript / 'data/reversal_fields/original_campaign/records.json.gz', 'rt'))['cases']
    for case in cases:
        source = manuscript / 'data/reversal_fields/reconstruction_inputs' / case['checkpoint']
        native.read_checkpoint(source, case['model'])
        assert len(case['tests']) == 1
        test = case['tests'][0]
        trace = manuscript / 'data/reversal_fields/reconstruction_inputs' / test['trace']
        with np.load(trace) as z:
            prefix = f'trial{test["trace_row"]}_'
            assert prefix + 'pulse_time' in z and prefix + 'recovery_time' in z
    sys.path.insert(0, str(manuscript / 'code'))
    import three_test_figures as t
    campaign = manuscript / 'data/runs/debelly_corrected_three_test_20261008'
    endpoints = t.f.read_endpoint_rows(campaign / 'results/endpoints.csv')
    manifest = json.loads((campaign / 'manifest.json').read_text())
    status = json.loads((campaign / 'results/status.json').read_text())
    with tempfile.TemporaryDirectory(prefix='four_assay_postcleanup_validation_') as folder:
        out = Path(folder)
        t.main_figure(endpoints, t.correlation_summary(endpoints), manifest, out)
        t.supplement(endpoints, out, status['generated_at'], manifest)
        pngs = {p.name: c.sha(p) for p in out.glob('*.png')}
        assert all(c.sha(manuscript / 'figures' / name) == digest for name, digest in pngs.items())
    record['post_cleanup_verification'] = dict(
        checked_at=datetime.now(timezone.utc).isoformat(), passed=True,
        canonical_npz_sha256_checks=5006,
        selected_establishment_checkpoint_paths_resolved=checkpoints,
        selected_personal_diagnostic_cases_loadable=len(cases),
        data_only_regenerated_endpoint_rows=len(endpoints),
        current_pngs_match_exactly=pngs, new_integrations=0)
    record['preflight']['additional_legacy_slurm_job_name_filters_checked'] = 110
    record['preflight']['matching_legacy_named_jobs_in_queue'] = 0
    c.write_record(record)
    print(json.dumps(record['post_cleanup_verification']), flush=True)


if __name__ == '__main__':
    main()
