"""Verify the compact Git codeset, exact restart inputs and plotted readouts.

This checks the shipped snapshot. Recomputing all native trial gates requires
a fresh SLURM campaign and audit_reversal.py, rather than the compact CSVs.
"""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from figure_builder import fit_summary, read_endpoint_rows

BASE = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inventory = read(BASE / 'CODESET_SHA256.json')
    for name, expected in inventory['files'].items():
        path = BASE / name
        if not path.is_file() or sha(path) != expected:
            raise ValueError('Codeset file missing or changed: ' + name)
    establishment = BASE / 'data/establishment_basal'
    with (establishment / 'results/endpoints.csv').open() as stream:
        original = list(csv.DictReader(stream))
    assert len(original) == 3234
    legi = [r for r in original if r['model'] == 'legi']
    assert len(legi) == 462 and all(r['status'] == 'unpolarized' for r in legi)
    models = {'wave_pinning', 'goryachev', 'debelly', 'holmes_model3', 'spring'}
    selected = [r for r in original if r['model'] in models and float(r['amplitude_factor']) == 1]
    assert len(selected) == 330
    checkpoints = 0
    for row in selected:
        if row['status'] != 'polarized':
            continue
        path = establishment / row['selected_checkpoint']
        assert sha(path) == row['selected_checkpoint_sha256']
        with np.load(path, allow_pickle=False) as z:
            assert z['orientation'][-1] < -.05 and z['contrast'][-1] >= .05
        checkpoints += 1
    assert checkpoints == 291
    results = BASE / 'data/runs/corrected_reversal_20261003/results'
    rows = read_endpoint_rows(results / 'endpoints.csv')
    assert len(rows) == 660
    assert Counter(r['status'] for r in rows) == Counter(success=375, censored=207, ineligible=78)
    audit = read(results / 'independent_audit.json')
    validation = read(results / 'validation.json')
    assert audit['passed'] and audit['complete'] and validation['passed']
    assert sha(results / 'independent_audit.json') == validation['independent_audit_sha256']
    audited = {(r['case_index'], r['definition']): r for r in audit['selected_readouts']}
    for row in rows:
        expected = audited[int(row['case_index']), row['definition']]
        assert row['status'] == expected['status'] and row['endpoint'] == expected['endpoint']
        if row['status'] == 'success':
            assert row['selected_field_tail_pass'] == expected['selected_tail_field_pass_0_1_percent']
    fits = fit_summary(rows)
    independently_checked = 0
    for fit in fits:
        if fit['fit_status'] != 'fit':
            continue
        chosen = [r for r in rows if all(r[k] == fit[k] for k in ('model', 'condition', 'definition'))
                  and r['status'] == 'success' and r['selected_field_tail_pass']]
        x = np.array([r['length'] for r in chosen])
        y = np.array([r['endpoint'] for r in chosen])
        slope = np.sum((x-x.mean())*(y-y.mean())) / np.sum((x-x.mean())**2)
        intercept = y.mean() - slope*x.mean()
        midpoint = (x.min()+x.max())/2
        r_squared = 1 - np.sum((y-(slope*x+intercept))**2)/np.sum((y-y.mean())**2)
        np.testing.assert_allclose([fit['slope'], fit['intercept'], fit['normalized_slope'], fit['r_squared']],
                                   [slope, intercept, slope*midpoint/(slope*midpoint+intercept), r_squared],
                                   rtol=1e-9, atol=1e-12)
        independently_checked += 1
    assert independently_checked == 18
    current=BASE/'data/runs/three_test_reversal_20261006'
    current_check={}
    if (current/'results/validation.json').exists():
        import three_test_figures as new
        current_rows=new.f.read_endpoint_rows(current/'results/endpoints.csv')
        assert len(current_rows)==660
        manifest=read(current/'manifest.json');validation=read(current/'results/validation.json')
        audit=read(current/'results/independent_audit.json')
        assert validation['passed'] and audit['passed'] and audit['complete']
        assert sha(current/'results/independent_audit.json')==validation['independent_audit_sha256']
        original=read(BASE/'data/runs/corrected_reversal_20261003/manifest.json')
        preps={(c['model'],c['length']):c for c in original['cases'] if c['condition']=='relative_20pct'}
        for c in manifest['cases']:
            prep=preps[c['model'],c['length']]
            for key in ['checkpoint','checkpoint_sha256','initial_orientation','initial_contrast','eligible']:
                assert c[key]==prep[key]
            factor=64. if c['model']=='goryachev' and c['condition']!='fixed_slope' else 1.
            assert c['response_factor']==factor and c['response_amplitude']==prep['amplitude']*factor
        audited={(c['case_index'],c['definition']):c for c in audit['selected_readouts']}
        for row in current_rows:
            expected=audited[int(row['case_index']),row['definition']]
            assert row['status']==expected['status'] and row['endpoint']==expected['endpoint']
            if row['status']=='success':assert row['selected_field_tail_pass']==expected['selected_tail_field_pass_0_1_percent']
        current_check=new.regression_check(current_rows,new.f.fit_summary(current_rows))
        mechanism=BASE/'data/runs/three_test_mechanisms_20261006/results'
        if (mechanism/'status.json').exists() and read(mechanism/'status.json').get('complete'):
            assert read(mechanism/'spatial_audit.json')['passed']
            cells=read(mechanism/'cell_evidence.json');assert len(cells)==15
            assert {(c['model'],c['geometry']) for c in cells}=={(m,g) for m in models for g in new.TESTS}
            current_check['mechanism_cells']=15
    print(json.dumps(dict(passed=True, files=len(inventory['files']), exact_checkpoints=checkpoints,
                          readouts=len(rows), independently_checked_fits=independently_checked,
                          current_comparison=current_check,
                          scope='Shipped inputs and summary evidence; full trial audit requires rerun.')))


if __name__ == '__main__':
    main()
