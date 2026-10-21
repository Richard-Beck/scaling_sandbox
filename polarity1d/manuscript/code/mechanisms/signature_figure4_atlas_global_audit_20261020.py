"""Independently verify complete S11 context, page, archive and lineage coverage."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pypdf import PdfReader


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 ** 2), b''):
            value.update(block)
    return value.hexdigest()


def audit(run):
    run = Path(run).resolve()
    destination = run / 'independent_global_coverage_audit.json'
    if destination.exists():
        raise FileExistsError('Global audit receipts are immutable')
    master = read(run / 'coverage_master.json')
    complete = read(run / 'complete.json')
    binding = read(run / 'physical_identity_binding.json')
    assert master['passed'] and master['all_required_contexts_complete']
    assert complete['coverage_master_sha256'] == sha(run / 'coverage_master.json')
    assert master['physical_identity_binding_sha256'] == sha(run / 'physical_identity_binding.json')
    required = ['common', 'E1_extra', 'E4_long']
    assert master['required_contexts'] == required == binding['required_contexts']
    indices = {entry['atlas_context']: entry for entry in master['page_data_indices']}
    assert len(indices) == len(master['page_data_indices']) == 3 and set(indices) == set(required)
    physical = {}
    pages = {}
    census = {}
    native_values = 0
    shard_bytes = 0
    shard_count = 0
    volumes = 0
    for context in required:
        root = run / 'contexts' / context
        inventory_path = run / 'canonical_inputs' / context / 'inventory.json'
        assert sha(inventory_path) == binding['source_case_inventory_hashes'][context]
        rows = read(inventory_path)['cases']
        assert len(rows) == binding['contexts'][context]['physical_runs']
        for case in rows:
            assert case['physical_run_id'] == case['id'] and case['id'] not in physical
            assert case['independent_stochastic_seed_count'] == int(case['kind'] == 'noise')
            physical[case['id']] = case
        entry = indices[context]
        assert sha(entry['path']) == entry['sha256']
        index = read(entry['path'])
        assert index['summary']['passed'] and index['summary']['float64_roundtrip_equal_every_field_time_actualC']
        assert entry['pages'] == len(index['pages']) == binding['contexts'][context]['expected_pages']
        assert entry['page_and_shard_map'] == index['pages']
        assert read(root / 'results/independent_native_atlas_audit.json')['passed']
        assert read(root / 'results/text_layout_audit.json')['passed']
        assert sha(entry['completion_gate']['path']) == entry['completion_gate']['sha256']
        assert read(entry['completion_gate']['path'])['passed']
        manifest = read(root / 'manifest.json')
        assert [page['slot'] for page in index['pages']] == list(range(len(manifest['jobs'])))
        pdf_pages = {}
        for volume in entry['volumes']:
            path = root / 'results/portable_export/pdf' / volume['file']
            assert sha(path) == volume['sha256'] and path.stat().st_size < 100 * 1024 ** 2
            assert len(PdfReader(path).pages) == volume['pages']
            assert volume['decoded_image_pixels_equal'] and volume['lossless_flate']
            pdf_pages[volume['file']] = set(range(1, volume['pages'] + 1))
            volumes += 1
        selected = 0
        for page in index['pages']:
            slot = page['slot']
            job = manifest['jobs'][slot]
            worker = read(root / 'workers' / f'{slot:04d}.json')
            assert page['condition'] == job['condition'] == worker['condition']
            positions = [worker['selection']['medoid']]
            if worker['selection']['random'] is not None:
                positions.append(worker['selection']['random'])
            expected = [job['cases'][position]['id'] for position in positions]
            assert page['selected_cases'] == expected
            assert len(expected) == (2 if page['condition']['kind'] == 'noise' else 1)
            selected += len(expected)
            with Path(page['selection_CSV']['path']).open() as stream:
                selection = list(csv.DictReader(stream))
            assert [row['case_id'] for row in selection] == expected
            assert sha(page['selection_CSV']['path']) == page['selection_CSV']['sha256']
            assert page['float64_roundtrip_equal'] and not page['sourcefield_resampling'] and not page['sourcefield_quantization']
            for companion in page['companions'].values():
                assert sha(companion['path']) == companion['sha256'] and companion['float64_roundtrip_equal']
                portable = root / 'results/portable_export/page_data' / Path(companion['path']).name
                assert sha(portable) == companion['sha256']
            for shard in page['native_field_CSV_shards']:
                path = Path(shard['path'])
                assert path.is_relative_to(root / 'results/native_field_CSV_archives')
                assert sha(path) == shard['sha256'] and shard['roundtrip_float64_equal']
                assert path.stat().st_size == shard['bytes'] < 90 * 1024 ** 2
                with gzip.open(path, 'rt') as stream:
                    header = next(csv.reader(stream))
                cells = max(record['cells'] for record in page['expected_native_binary_hashes'].values())
                assert header == ['case_index', 'native_frame', 'time', 'actualC'] + [f'active_cell_{i:04d}' for i in range(cells)]
                shard_bytes += shard['bytes']
                shard_count += 1
            native_values += page['native_cell_time_values']
            assert page['pdf_page'] in pdf_pages[page['pdf_volume']]
            pdf_pages[page['pdf_volume']].remove(page['pdf_page'])
            pages[(context, page['pdf_volume'], page['pdf_page'])] = page
        assert all(not remaining for remaining in pdf_pages.values())
        assert selected == index['summary']['selected_histories']
        census[context] = dict(candidates=len(rows), pages=len(index['pages']), selected_histories=selected)
    references = master['E1_common_reference_crosswalk']
    assert len(references) == census['E1_extra']['candidates']
    assert len({row['source_physical_run_id'] for row in references}) == len(references)
    for reference in references:
        branch = physical[reference['source_physical_run_id']]
        common = physical[reference['common_sigma1_physical_run_id']]
        mean = physical[reference['common_sigma0_mean_physical_run_id']]
        zero = physical[reference['common_zero_physical_run_id']]
        assert branch['atlas_context'] == 'E1_extra' and common['atlas_context'] == 'common'
        assert branch['latent_Gaussian_cohort_id'] == common['latent_Gaussian_cohort_id'] == reference['latent_Gaussian_cohort_id']
        assert common['kind'] == 'noise' and common['sigma'] == 1 and mean['kind'] == 'mean' and zero['kind'] == 'zero'
        assert mean['f'] == common['f'] and zero['f'] == 0
        for key in ['model', 'C_label', 'f', 'seed', 'native_horizon']:
            assert branch[key] == common[key]
        for other in [mean, zero]:
            assert other['model'] == common['model'] and other['native_horizon'] == common['native_horizon'] and other['C_label'] == common['C_label']
    prefix_count = 0
    for case in physical.values():
        if case.get('retained_prefix'):
            original = physical[case['prefix_and_topup_relation']['original_physical_run_id']]
            assert original['atlas_context'] == 'common' and case['atlas_context'] == 'E4_long'
            assert case['latent_seed_lineage_id'] == original['latent_seed_lineage_id']
            prefix_count += 1
    for selector in master['longest_page_selector']:
        source = pages[(selector['atlas_context'], selector['pdf_volume'], selector['pdf_page'])]
        assert selector['condition'] == source['condition'] and selector['selected_cases'] == source['selected_cases']
    assert len(physical) == master['independence_and_reuse_census']['physical_runs']
    assert prefix_count == master['independence_and_reuse_census']['E4_retained_prefix_runs']
    result = dict(passed=True, complete=True, generated_utc=datetime.now(timezone.utc).isoformat(),
        coverage_master_sha256=sha(run / 'coverage_master.json'), complete_sha256=sha(run / 'complete.json'),
        physical_identity_binding_sha256=sha(run / 'physical_identity_binding.json'),
        contexts=census, all_required_contexts_complete=True, PDF_volumes=volumes,
        page_count=len(pages), native_cell_time_values=native_values, native_CSV_shards=shard_count,
        all_native_CSV_bytes_hashed=shard_bytes, exact_E1_common_references=len(references),
        exact_E4_prefix_lineages=prefix_count, longest_selector_pages=len(master['longest_page_selector']),
        source_sha256=sha(__file__), native_integrations=0, main_example_selection=False,
        scope='Independent full context/candidate/page/volume/selected-ID coverage, current native CSV bytes and exact physical/Gaussian/prefix reference crosswalks. Numerical native-field roundtrips and decoded-image equality remain bound to their accepted per-context proofs.')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    print(json.dumps(audit(parser.parse_args().run)))
