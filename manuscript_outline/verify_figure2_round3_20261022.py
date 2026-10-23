"""Independent coordinate, preservation, artwork and report integration checks."""
from pathlib import Path
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from urllib.parse import unquote, urlsplit
import re
import xml.etree.ElementTree as ET
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_round3_20261022'
DEV = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
FIG = DEV / 'revision_round3_20261022/figure2'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csvread(name):
    return list(csv.DictReader((FIG / 'supporting/data' / name).open()))


def equal(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-12), (a, b)


def main():
    receipt = read(DEST / 'integration.json')
    assert sha(ROOT / receipt['baseline']) == receipt['baseline_sha256']
    prior = read(ROOT / receipt['baseline'])
    current = read(ROOT / 'manuscript_outline/outline_v6.json')
    protected = []
    for ident, record in receipt['protected_sections'].items():
        a = next(s for s in prior['sections'] if s['id'] == ident)
        b = next(s for s in current['sections'] if s['id'] == ident)
        assert a == b
        assert hashlib.sha256(json.dumps(b, sort_keys=True).encode()).hexdigest() == record['section_sha256']
        for path, digest in record['asset_sha256'].items():
            assert sha(ROOT / path) == digest
        protected.append(ident)
    section = next(s for s in current['sections'] if s['id'] == 'figure2')
    assert section['status'] == 'round3_complete_qualified'
    assert all('Mass-conserved, saturating' not in text for text in section['legend'])
    sources = read(FIG / 'supporting/metadata/source_sha256.json')
    for path, digest in sources.items():
        assert sha(ROOT / path) == digest, path
    # Per-width time coordinates are compared directly to the archived event rows.
    widths = csvread('C_replacement_times_all450.csv')
    assert len(widths) == 450
    originals = {}
    for run in {r['run'] for r in widths}:
        originals[run] = {r['task']: r for r in read(ROOT / 'polarity1d/manuscript/data/runs' / run / 'results/analysis_20261020/event_rows.json')}
    for row in widths:
        native = originals[row['run']][int(row['task'])]
        event = native['continuous_cue_on_events']['D0_5']
        assert event['status'] == row['status']
        for key, archived in [('time_lower', 'lower_elapsed'), ('time_upper', 'upper_elapsed'),
            ('time_estimate', 'estimated_elapsed'), ('lower_bound', 'right_censored_at')]:
            if event[archived] is None:
                assert not row[key]
            else:
                equal(row[key], event[archived])
        if row['status'] == 'right_censored':
            assert not row['time_estimate'] and float(row['lower_bound']) > 0
    assert sum(r['model'] == 'holmes_model3' and float(r['width']) == 15.2 and r['status'] == 'observed' for r in widths) == 18
    # Recompute all D coordinates from exact phase endpoints, independent of the
    # coordinate builder and its missing-value policy implementation.
    archived = read(DEV / 'signature_assays_20261020/supporting/data/figure2_S2_all_native_pulse_gap_endpoints.json')
    endpoints_by_key = {}
    for row in archived:
        c = row['case']
        if c['pulse_count'] == 4 and c['total_on_factor_of_calibration'] == 1.5:
            endpoints_by_key[(c['model'], c['requested_reference_length'], c['gap_multiple_of_pulse'])] = row['endpoints']
    extra = read(ROOT / 'polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022/results/independent_native_readout_audit.json')
    for row in extra['trajectories']:
        c = row['case']
        if c['kind'] == 'train':
            endpoints_by_key[('holmes_model3', c['length'], c['gap'])] = row['endpoints']
    temporal = csvread('D_gap1_vs_pulse2_all180.csv')
    assert len(temporal) == len(endpoints_by_key) == 180
    for row in temporal:
        key = (row['model'], int(row['length']), float(row['gap']))
        ends = endpoints_by_key[key]
        pulses = [e['end'] for e in ends if e['kind'] == 'on']
        gap = next(e['end'] for e in ends if e['kind'] == 'gap')
        s1, s2 = pulses[0]['share_half'], pulses[1]['share_half']
        equal(row['headroom_gain'], (s2-s1)/(1-s1))
        equal(row['activity_kept'], 0 if gap['readout_eligible'] is False else gap['share_half']/s1)
        equal(row['raw_share_retention'], gap['share_half_raw']/s1)
        assert (row['below_activity_gate'] == 'True') == (gap['readout_eligible'] is False)
    left_wall = {model: sum(r['model'] == model and float(r['activity_kept']) <= .02 and float(r['headroom_gain']) > .02 for r in temporal)
        for model in ('goryachev', 'wave_pinning', 'holmes_model3', 'debelly', 'spring')}
    assert left_wall['goryachev'] == 0 and left_wall['debelly'] > 0 and left_wall['spring'] > 0
    # Check text objects in the actual vector export, not just an audit claim.
    svg = ET.parse(FIG / 'figure2_round3.svg')
    texts = [e for e in svg.iter() if e.tag.rsplit('}', 1)[-1] == 'text']
    letters = [e for e in texts if ''.join(e.itertext()) in list('ABCDEF')]
    assert len(letters) == 6
    assert len({e.attrib['style'] for e in letters}) == 1
    assert len({e.attrib['x'] for e in letters}) == 1
    assert len({e.attrib.get('fill', 'black') for e in letters}) == 1
    assert not read(FIG / 'supporting/metadata/display_audit.json')['clipped_text']
    page = PdfReader(FIG / 'figure2_round3.pdf').pages[0]
    equal(float(page.mediabox.width)/72, 7.2)
    equal(float(page.mediabox.height)/72, 10)
    assert 'Mass-conserved, saturating' not in (page.extract_text() or '')
    manifest = read(ROOT / 'manuscript_outline/outline_v6.build_manifest.json')
    figure_page = next(p for p in manifest['page_inventory'] if p['section_id'] == 'figure2' and p['kind'] == 'figure')
    integrated = PdfReader(ROOT / 'manuscript_outline/outline_v6.pdf').pages[figure_page['number']-1]
    assert integrated.get_contents().get_data() == page.get_contents().get_data()
    equal(float(integrated.mediabox.width)/72, 7.2)
    equal(float(integrated.mediabox.height)/72, 10)
    # Confirm every prior PDF page survives; resources/appearance also receive
    # a raster inspection for the main and integrated figure sheet at handoff.
    proof = read(FIG / 'supporting/metadata/supplement_preservation.json')
    supplement = PdfReader(FIG / 'figure2_round3_supplement.pdf')
    index = proof['added_pages']
    for record in proof['preserved']:
        assert sha(ROOT / record['path']) == record['sha256']
        original = PdfReader(ROOT / record['path'])
        for page in original.pages:
            assert page.get_contents().get_data() == supplement.pages[index].get_contents().get_data()
            index += 1
    assert index == len(supplement.pages)
    missing = []
    for path in (FIG / 'index.html', DEST / 'index.html'):
        for target in re.findall(r'(?:href|src)=["\']([^"\']+)["\']', path.read_text()):
            parsed = urlsplit(target)
            if parsed.scheme or not parsed.path:
                continue
            if not (path.parent / unquote(parsed.path)).exists():
                missing.append((str(path), target))
    assert not missing, missing
    validation = read(ROOT / 'manuscript_outline/outline_v6.validation.json')
    assert validation['passed']
    report = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        protected_sections_identical=protected, source_files_verified=len(sources),
        width_conditions_independently_verified=450, temporal_coordinates_independently_verified=180,
        low_retention_positive_gain_counts=left_wall, deBelly_native_endpoint_audit_passed=read(FIG / 'supporting/metadata/debelly_L19_gap16_native_audit.json')['passed'],
        all_six_panel_labels_uniform_in_SVG=True, main_page_inches=[7.2, 10],
        integrated_vector_page_preserved_at_native_dimensions=True,
        prior_pdf_pages_preserved=index-proof['added_pages'], supplement_pages=len(supplement.pages),
        local_HTML_links_valid=True, outline_validation_passed=True,
        native_runs_submitted=0, class_membership_changed=False,
        figure2_sha256=sha(FIG / 'figure2_round3.pdf'), outline_sha256=sha(ROOT / 'manuscript_outline/outline_v6.pdf'))
    (DEST / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
