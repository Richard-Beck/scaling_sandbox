"""Check preservation, source coordinates, native PDF integration and local links."""
from pathlib import Path
from copy import deepcopy
from datetime import datetime, timezone
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import argparse
import csv
import math
import hashlib
import json
import subprocess

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_round6_20261023'
DEV = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
FIGURES = DEV / 'revision_round6_20261023'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and value:
                self.links.append(value)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--require-complete-mode-campaign', action='store_true')
    args = parser.parse_args()
    starting = read(DEST / 'request/starting_checkpoint.json')
    for name, expected in starting['input_sha256'].items():
        actual = DEST / 'request/revision_brief_round6_20261023.md' if name == 'RevisionBrief.md' else ROOT / name
        assert sha(actual) == expected, name
    original_brief = subprocess.check_output(['git', 'show', starting['baseline_commit'] + ':RevisionBrief.md'], cwd=ROOT)
    assert original_brief == (DEST / 'request/revision_brief_round6_20261023.md').read_bytes()
    assert (ROOT / 'CurrentNarrative.md').read_bytes() == (DEST / 'request/current_narrative_20261023.md').read_bytes()
    baseline = read(ROOT / 'manuscript_outline/outline_v11.json')
    current = read(ROOT / 'manuscript_outline/outline_v12.json')
    for index in (4, 5):
        expected = deepcopy(baseline['sections'][index])
        if index == 4:
            for source in expected['sources']:
                if source.get('path') == 'manuscript_outline/assets/migration_v2/provenance.json':
                    source['label'] = 'Migration display and preserved-input provenance'
        assert current['sections'][index] == expected
    old_pilot = read(ROOT / 'manuscript_outline/revision_figures3_4_20261023/pilot_observables.json')
    new_pilot = read(DEST / 'pilot_observables.json')
    assert len(old_pilot['observables']) == len(new_pilot['observables'])
    for index, (old, new) in enumerate(zip(old_pilot['observables'], new_pilot['observables'])):
        for key in ('observable', 'drift_direction', 'secondary_direction', 'status'):
            assert old[key] == new[key], (index, key)
        if index != 2:
            assert old['operational_definition'] == new['operational_definition'], index
    for key in ('fixed_before_data', 'limitations'):
        assert old_pilot[key] == new_pilot[key], key
    assert new_pilot['interpretation_notes'][:len(old_pilot['interpretation_notes'])] == old_pilot['interpretation_notes']
    assert new_pilot['regime_definition'].startswith(old_pilot['regime_definition'])
    new_section = next(s for s in current['supplementary'] if s.get('id') == 'pilot_observables')
    assert new_section['measurements'] == [r['operational_definition'] for r in new_pilot['observables']]
    assert 'Fig. 4J' in new_pilot['observables'][2]['counterpart']
    source_pairs = [(ROOT / 'polarity1d/manuscript/data/derived/length_correlations.csv',
                     FIGURES / 'figure1/supporting/data/E_length_correlations.csv')]
    for name in ('A_incumbent_organization.csv', 'B_initiation.csv', 'C_replacement_times_all450.csv',
                 'D_gap1_vs_pulse2_all180.csv', 'E_L39_challenger_share_traces.csv'):
        source_pairs.append((DEV / 'revision_round3_20261022/figure2/supporting/data' / name,
                             FIGURES / 'figure2/supporting/data' / name))
    for old, new in source_pairs:
        assert old.read_bytes() == new.read_bytes(), str(new)
    layout_csvs = 0
    for number in (1, 2):
        for path in (FIGURES / f'figure{number}/supporting/data').glob('*.csv'):
            relative_path = str(path.relative_to(ROOT))
            committed = subprocess.check_output(['git', 'show', '209a0c3:' + relative_path], cwd=ROOT)
            assert path.read_bytes() == committed, relative_path
            layout_csvs += 1
    with (FIGURES / 'figure2/supporting/data/C_replacement_times_all450.csv').open() as stream:
        replacement_rows = list(csv.DictReader(stream))
    replacement = {(r['model'], float(r['area']), float(r['width']), float(r['length'])): r
                   for r in replacement_rows}
    reference_folds = []
    reference_models = set()
    for (model, area, width, length), small in replacement.items():
        if area == 1 and length == 19:
            large = replacement[(model, area, width, 39.0)]
            if small['time_estimate'] and large['time_estimate']:
                reference_folds.append(float(large['time_estimate']) / float(small['time_estimate']))
                reference_models.add(model)
    assert len(reference_folds) == 18 and len(reference_models) == 5
    assert min(reference_folds) >= 4.8 and max(reference_folds) <= 12
    assert 'reference cue area 1' in current['sections'][1]['narrative'][0]
    unchanged = subprocess.run(['git', 'diff', '--exit-code', starting['baseline_commit'], '--',
        'author', 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round5_20261023',
        'manuscript_outline/outline_v11.pdf', 'manuscript_outline/outline_v11.json'], cwd=ROOT, capture_output=True, text=True)
    assert unchanged.returncode == 0, unchanged.stdout
    manifest = read(ROOT / 'manuscript_outline/outline_v12.build_manifest.json')
    pdf = PdfReader(ROOT / 'manuscript_outline/outline_v12.pdf')
    cover_text = pdf.pages[0].extract_text()
    assert all(f'Figure {number}' in cover_text for number in range(1, 6))
    note_headings = ('Detailed legend', 'Narrative points', 'Proposed measurements',
                     'Scope, limitations', 'Model definitions', 'Sources and reproducibility')
    orphan_pages = [p['number'] for p in manifest['page_inventory']
                    if p['kind'] == 'figure_notes' and p['text_blocks'] <= 8
                    and not any(h in pdf.pages[p['number']-1].extract_text() for h in note_headings)]
    assert not orphan_pages, orphan_pages
    native_count = 0
    for page in manifest['page_inventory']:
        if page.get('native_vector_pdf'):
            original = PdfReader(ROOT / page['native_vector_pdf']).pages[0]
            integrated = pdf.pages[page['number']-1]
            assert original.get_contents().get_data() == integrated.get_contents().get_data()
            native_count += 1
    sources = 0
    for section in current['sections'][:4]:
        folder = FIGURES / section['id']
        assert (folder / 'scientific_readme.md').read_text().startswith('# Unexpected')
        for source in section.get('main_source_data', []):
            assert sha(ROOT / source['path']) == source['sha256']
            sources += 1
        assert Path(ROOT / section['main_pdf']).is_file()
        assert Path(ROOT / section['supplement_pdf']).is_file()
    links_checked = 0
    for path in DEST.glob('*.html'):
        parser = Links()
        parser.feed(path.read_text())
        for link in parser.links:
            parsed = urlsplit(link)
            if not parsed.scheme and parsed.path:
                target = (path.parent / unquote(parsed.path)).resolve()
                assert target.exists(), (path, link)
                links_checked += 1
    preservation = DEST / 'review_preservation.json'
    f3_audit = read(FIGURES / 'figure3/supporting/metadata/verification.json')
    assert f3_audit['complete'] and f3_audit['AB_data_identical']
    assert f3_audit['worker_source_matches_manifest']
    assert f3_audit['native_stochastic_runs'] == f3_audit['retuned_parameters'] == 0
    for key, name in (('main_pdf_sha256', 'figure3_round6.pdf'),
                      ('main_png_sha256', 'figure3_round6.png'), ('legend_sha256', 'legend.json')):
        assert sha(FIGURES / 'figure3' / name) == f3_audit[key], name
    with (FIGURES / 'figure3/supporting/data/alpha_bridge.csv').open() as stream:
        alpha_rows = list(csv.DictReader(stream))
    alpha_eligible = [r for r in alpha_rows if r['bridge_eligible'].lower() == 'true']
    assert len(alpha_rows) == f3_audit['drift_alpha_requested']
    assert len(alpha_eligible) == f3_audit['drift_alpha_pairs']
    assert sum(r['model'] == 'wave_pinning' for r in alpha_eligible) == f3_audit['wave_pinning_alpha_pairs']
    for row in alpha_eligible:
        assert float(row['fraction']) <= 0.01
        assert abs(float(row['prediction_fraction']) - float(row['fraction'])) < 1e-12
        assert row['matched_input_prediction'].lower() == 'true' and row['extrapolated'].lower() == 'false'
        assert float(row['neighborhood_max_multiple_domain_fraction']) <= 0.05
        assert row['linear_response_eligible'].lower() == 'true'
        assert all(math.isfinite(float(row[key])) for key in ('predicted_alpha', 'observed_alpha'))
    f4_audit = read(FIGURES / 'figure4/supporting/metadata/scientific_validation.json')
    assert f4_audit['passed']
    postreview_data_count = 0
    postreview_archive_count = 0
    if args.require_complete_mode_campaign:
        storage = read(ROOT / 'reports/round6_storage_20261023.json')
        assert storage['campaign_status'] == 'NATIVE_COHORT_COMPLETE'
        assert storage['complete_mode_histories'] == storage['allocated_mode_histories']
        assert storage['failed_mode_histories'] == 0
        precision = read(ROOT / 'polarity1d/manuscript/data/runs/figure4_round6_modes_20261023/precision/status.json')
        assert precision.get('early') is False
        assert len(precision['conditions']) == 200
        assert len({r['condition_id'] for r in precision['conditions']}) == 200
        assert all(r['n_seeds'] in (24, 96) for r in precision['conditions'])
        assert sum(r['n_seeds'] for r in precision['conditions']) == storage['allocated_mode_histories']
        terminal = read(FIGURES / 'figure4/supporting/metadata/stable_artifact_handshake.json')
        assert terminal['stable'] and terminal['no_more_writers']
        assert terminal['controllers_completed'] and terminal['live_producer_stopped']
        assert terminal['native_histories'] == storage['complete_mode_histories']
        assert terminal['native_failures'] == 0
        artifact_hashes = terminal['artifact_sha256']
        presentation_path = FIGURES / 'figure4/supporting/metadata/postreview_presentation_audit.json'
        if presentation_path.exists():
            presentation = read(presentation_path)
            assert presentation['numeric_data_exports_unchanged']
            assert presentation['numeric_data_export_count'] == 62
            for name, expected in artifact_hashes.items():
                assert presentation['prior_artifact_sha256'][name] == expected, name
            artifact_hashes = presentation['artifact_sha256']
            assert terminal['artifact_sha256'].keys() <= artifact_hashes.keys()
            terminal_path = FIGURES / 'figure4/supporting/metadata/stable_artifact_handshake.json'
            recorded_terminal = subprocess.check_output(['git', 'show', '028d9f5:' + str(terminal_path.relative_to(ROOT))], cwd=ROOT)
            assert terminal_path.read_bytes() == recorded_terminal
            export_path = FIGURES / 'figure4/supporting/metadata/data_export_hashes.json'
            original_exports = json.loads(subprocess.check_output(['git', 'show', '028d9f5:' + str(export_path.relative_to(ROOT))], cwd=ROOT))
            assert read(export_path) == original_exports
            archive = ROOT / presentation['reviewed_archive']
            for name, expected in presentation['reviewed_archive_file_sha256'].items():
                original = subprocess.check_output(['git', 'show', '028d9f5:' + str((FIGURES / 'figure4' / name).relative_to(ROOT))], cwd=ROOT)
                assert (archive / name).read_bytes() == original, name
                assert sha(archive / name) == expected, name
                postreview_archive_count += 1
            for name, expected in presentation['presentation_source_sha256'].items():
                assert sha(FIGURES / 'figure4' / name) == expected, name
            freeze = read(FIGURES / 'figure4/supporting/metadata/final_analysis_source_freeze.json')
            for item in freeze['files']:
                assert sha(FIGURES / 'figure4/supporting/metadata/final_analysis_source_snapshot' / item['snapshot']) == item['sha256']
                if item['source'] != str((FIGURES / 'figure4/rebuild_figure4.py').relative_to(ROOT)):
                    assert sha(ROOT / item['source']) == item['sha256']
            for name, expected in original_exports.items():
                assert sha(FIGURES / 'figure4/supporting/data' / name) == expected, name
                postreview_data_count += 1
        for name, expected in artifact_hashes.items():
            assert sha(FIGURES / 'figure4' / name) == expected, name
        evidence = read(FIGURES / 'figure4/supporting/metadata/terminal_evidence_audit.json')
        assert evidence['passed'] and evidence['source_hashes_all_match']
        assert evidence['complete_valid_receipts'] == evidence['unique_seeds'] == terminal['native_histories']
        assert evidence['native_profiles_saved'] is False
    snapshot_count = 0
    for folder, receipt_name, key in ((DEST / 'review_snapshot_initial', 'capture.json', 'sha256'),
            (FIGURES / 'figure4/supporting/review_snapshot_initial', 'SNAPSHOT_HASHES.json', 'files')):
        for relative_name, expected in read(folder / receipt_name)[key].items():
            assert sha(folder / relative_name) == expected, str(folder / relative_name)
            snapshot_count += 1
    if preservation.exists():
        for item in read(preservation)['reviews']:
            assert sha(ROOT / item['path']) == item['sha256'], item['path']
    report = {'generated_utc': datetime.now(timezone.utc).isoformat(), 'passed': True,
        'baseline_commit': starting['baseline_commit'], 'baseline_pdf_preserved': True,
        'round5_files_and_author_documents_preserved': True, 'protected_sections': ['figure5', 'figure6'],
        'figure5_source_label_presentation_update_only': True,
        'pilot_prospective_directions_preserved': True,
        'pilot_references_and_model_conditioning_updated': True,
        'byte_identical_figure1_2_data_tables': len(source_pairs),
        'figure1_2_plotted_csvs_unchanged_by_review_fixes': layout_csvs,
        'figure2_reference_area1_resolved_size_pairs': len(reference_folds),
        'figure2_reference_area1_replacement_fold_range': [min(reference_folds), max(reference_folds)],
        'main_source_hashes_checked': sources, 'native_pdf_pages_preserved': native_count,
        'opening_page_contains_figures1_5_results': True,
        'isolated_source_link_pages': orphan_pages,
        'local_html_links_checked': links_checked,
        'immutable_snapshot_hashes_checked': snapshot_count,
        'figure3_numerical_and_artifact_audit_checked': True,
        'figure3_qualified_exact_input_alpha_pairs': len(alpha_eligible),
        'figure4_numerical_audit_checked': True,
        'figure4_postreview_numerical_exports_verified': postreview_data_count,
        'figure4_postreview_archived_commit_files_verified': postreview_archive_count,
        'terminal_mode_campaign_required_and_checked': args.require_complete_mode_campaign,
        'review_preservation_checked': preservation.exists(),
        'scope': 'Source and packaging audit; model-specific numerical and scientific checks are recorded with Figures 3–4.'}
    (DEST / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
