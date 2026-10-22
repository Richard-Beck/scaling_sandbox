"""Verify scientific-export provenance and exact outline protection for Figure 4."""
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
from pathlib import Path

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4'
DEST = ROOT/'manuscript_outline/revision_round3_figure4_20261022'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scientific_cohort():
    """Check the portable scientific census without reopening native records."""
    metadata = FIG/'supporting/metadata'
    data = FIG/'supporting/data'
    analysis = read(metadata/'analysis_receipt.json')
    nulls = read(metadata/'null_analysis_receipt.json')
    assert nulls['A']['replicates_per_null_size'] == 1000
    assert nulls['B']['conditions'] == 615
    with (data/'null_B_analytic_agreement.csv').open(newline='') as stream:
        checks = list(csv.DictReader(stream))
    assert len(checks) == 615
    eligible = [r for r in checks if r['status'] == 'precision_met']
    assert all(float(r['SE_log_k']) <= .01 for r in eligible)
    result = dict(null_conditions=615, null_precision_qualified=len(eligible))
    if not analysis.get('completed_campaign'):
        result['model_cohort'] = 'provisional_round2'
        return result
    manifest = read(metadata/'campaign_manifest.json')
    assert sha(metadata/'campaign_manifest.json') == analysis['campaign_manifest_sha256']
    design = {(c['model'], float(c['size']), float(c['fraction'])): c['condition_id']
              for c in manifest['conditions']}
    allocation = read(metadata/'campaign_final_allocation.json')['targets']
    precision = read(metadata/'campaign_final_precision.json')['conditions']
    precision = {r['condition_id']: r for r in precision}
    with (data/'analysis_lag_metrics.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    observed = {(r['model'], float(r['C_label']), float(r['f'])): r for r in rows}
    assert len(design) == len(observed) == len(rows) == 260
    assert set(observed) == set(design)
    assert not analysis['unknown_forcing_conditions']
    for key, row in observed.items():
        cid = design[key]
        n = int(row['n_seeds'])
        assert n == allocation[cid] == precision[cid]['n_seeds']
        assert n in (96, 192, 288, 384)
        if precision[cid]['se_log_k'] is not None:
            assert math.isclose(float(row['SE_log_k']), precision[cid]['se_log_k'],
                                rel_tol=1e-12, abs_tol=1e-14), cid
            assert float(row['analysis_lag']) == precision[cid]['analysis_lag']
            assert float(row['SE_log_k']) <= .07 or n == 384
        elif precision[cid]['analysis_lag'] is not None:
            assert n == 384 and float(row['bootstrap_valid_fraction']) < .975
        else:
            assert row['status'] == 'unavailable_no_analysis_lag'
    result.update(model_cohort='completed_compact_campaign', model_conditions=260,
                  independent_seed_count=sum(allocation.values()),
                  all_precision_cohorts_match=True, unknown_forcing_conditions=0)
    return result


def main():
    prior = read(ROOT/'manuscript_outline/outline_v8.json')
    current = read(ROOT/'manuscript_outline/outline_v9.json')
    integration = read(DEST/'integration.json')
    assert sha(ROOT/integration['baseline']) == integration['baseline_sha256']
    assert sha(ROOT/integration['outline']) == integration['outline_sha256']
    for ident, record in integration['protected_sections'].items():
        assert next(s for s in prior['sections'] if s['id'] == ident) == next(s for s in current['sections'] if s['id'] == ident)
        for path, expected in record['asset_sha256'].items():
            assert sha(ROOT/path) == expected
    assert set(integration['protected_sections']) == {'figure1','figure2','figure3','figure5','figure6'}
    render = read(FIG/'supporting/metadata/render_receipt.json')
    for path, expected in render['inputs'].items():
        assert sha(FIG/path) == expected, path
    for path, expected in render['outputs'].items():
        assert sha(FIG/path) == expected, path
    for path, expected in render.get('documentation', {}).items():
        assert sha(FIG/path) == expected, path
    assert sha(ROOT/render['renderer']) == render['renderer_sha256']
    for path, expected in render.get('code_sha256', {}).items():
        assert sha(ROOT/path) == expected, path
    main_page = PdfReader(FIG/'figure4_main.pdf')
    assert len(main_page.pages) == 1
    page = main_page.pages[0]
    width, height = float(page.mediabox.width)/72, float(page.mediabox.height)/72
    assert abs(width*25.4-180) < 1e-8
    assert abs(height*25.4-render['height_mm']) < 1e-8
    text = page.extract_text()
    for letter in 'ABCDE':
        assert letter in text.split(), letter
    for name in ('Wave-pinning', 'Holmes', 'Spring', 'Goryachev'):
        assert name in text
    assert 'Null models for size-dependent heading persistence' not in text
    assert 'sketch of' not in text and 'what C' not in text
    publication = PdfReader(DEST/'figure4_with_legend.pdf')
    assert len(publication.pages) == 2
    assert publication.pages[0].get_contents().get_data() == page.get_contents().get_data()
    assert sha(DEST/'figure4_with_legend.pdf') == integration['publication_pdf_sha256']
    assert read(ROOT/'manuscript_outline/outline_v9.validation.json')['passed']
    scientific = read(FIG/'supporting/metadata/estimator_verification.json')
    assert scientific['passed']
    assert sha(ROOT/scientific['test_source']) == scientific['test_sha256']
    analysis = read(FIG/'supporting/metadata/analysis_receipt.json')
    assert bool(analysis.get('completed_campaign')) == integration['campaign_complete']
    cohort = scientific_cohort()
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        figure4_status=integration['figure4_status'], campaign_complete=integration['campaign_complete'],
        exact_vector_main_embedded=True, width_mm=180, height_mm=render['height_mm'],
        protected_sections=list(integration['protected_sections']),
        source_hashes_verified=len(render['inputs'])+len(render['outputs'])+len(render.get('documentation', {})),
        scientific_checks=scientific['test_count'], supplement_pages=len(PdfReader(FIG/'figure4_supplement.pdf').pages),
        scientific_cohort=cohort)
    (DEST/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
