#!/usr/bin/env python3
"""Compose qualified vector plates, legends and explicit supplementary PDFs.

This launches no integrations. The input manifest explicitly names three figures,
their content JSON, vector plates, captions, supporting PDFs and qualification
receipts. Native profile appendices remain separate, labelled data appendices.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import subprocess
import sys
from urllib.parse import quote, unquote, urlsplit, urlunsplit

from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, TextStringObject

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    spec = json.loads(args.manifest.read_text())
    assert [f['number'] for f in spec['figures']] == [2, 3, 4]
    assert spec['composition_ready'] is True
    expected = spec['inputs_sha256']
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    inputs = {}
    rebased_links = []

    def rebase_links(reader, original):
        """Keep local source annotations valid after merging into the output folder."""
        for page in reader.pages:
            for ref in page.get('/Annots', []):
                action = ref.get_object().get('/A')
                if not action or '/URI' not in action:
                    continue
                old = str(action['/URI'])
                parts = urlsplit(old)
                if parts.scheme or parts.netloc or not parts.path:
                    continue
                target = (original.parent / unquote(parts.path)).resolve()
                assert target.exists(), ('broken source annotation', str(original), old)
                new = urlunsplit(('', '', quote(os.path.relpath(target, out), safe='/'),
                                  parts.query, parts.fragment))
                action[NameObject('/URI')] = TextStringObject(new)
                rebased_links.append(dict(source=str(original), old=old, new=new))

    def source(relative):
        path = (ROOT / relative).resolve()
        assert path.is_file(), path
        name = str(path.relative_to(ROOT))
        digest = sha(path)
        assert name in expected and expected[name] == digest, ('selected input changed', name)
        inputs[name] = digest
        return path

    inputs[str(args.manifest.resolve().relative_to(ROOT))] = sha(args.manifest)
    source(str(Path(__file__).resolve().relative_to(ROOT)))
    sections = []
    captions = []
    for figure in spec['figures']:
        for receipt in figure['qualification_receipts']:
            evidence = json.loads(source(receipt['path']).read_text())
            assert evidence[receipt.get('field', 'passed')] == receipt.get('expected', True), receipt
        content = json.loads(source(figure['content']).read_text())
        section = copy.deepcopy(content['sections'][0] if 'sections' in content else content.get('section', content))
        assert section['legend'] and section['title']
        section['assets'] = []
        section['role'] = ''
        for key in ('narrative', 'measurements', 'limitations', 'supplementary', 'equations'):
            section[key] = []
        # Keep methods/reconstruction links in the HTML index. Appending them
        # to print legends creates otherwise empty continuation pages.
        section['sources'] = []
        section['id'] = f'figure{figure["number"]}'
        section['title'] = f'Figure {figure["number"]}. {figure["label"]}'
        # The renderer accepts mathtext but does not interpret Markdown emphasis.
        section['legend'] = [paragraph.strip().replace('**', '')
                             for block in section['legend']
                             for paragraph in block.split('\n\n')
                             if paragraph.strip() and not paragraph.strip().startswith('# ')]
        for item in figure.get('legend_sources', []):
            source(item['path'])
        sections.append(section)
        captions.append(source(figure['caption']).read_text().strip())
        for name in ('main_pdfs', 'supplementary_pdfs', 'appendix_pdfs'):
            for item in figure.get(name, []):
                source(item['path'])
        source(figure['preview'])
    legend_doc = dict(title='Figures 2–4: standalone legends', date='19 October 2026',
                      stage='Final scientific consolidation', overview=[], sections=sections)
    legend_json = out / 'supporting/metadata/legends.json'
    write_json(legend_json, legend_doc)
    renderer = ROOT / 'scripts/build_manuscript_outline_v2.py'
    for path in (renderer, ROOT / 'scripts/build_manuscript_outline.py'):
        source(str(path.relative_to(ROOT)))
    rendered_legends = out / 'supporting/vectors/legends_render.pdf'
    subprocess.run([sys.executable, str(renderer), '--input', str(legend_json),
                    '--output', str(rendered_legends)], check=True)
    build = json.loads(rendered_legends.with_suffix('.build_manifest.json').read_text())
    assert all(build['validation'][k] for k in (
        'all_local_sources_present', 'all_text_within_page_bounds',
        'all_text_nonoverlapping', 'all_pages_have_extractable_text'))
    legend_reader = PdfReader(rendered_legends)
    rebase_links(legend_reader, rendered_legends)
    legends_only = PdfWriter()
    combined = PdfWriter()
    main_inventory = []
    for figure in spec['figures']:
        first = len(combined.pages) + 1
        for item in figure['main_pdfs']:
            path = source(item['path'])
            reader = PdfReader(path)
            rebase_links(reader, path)
            assert len(reader.pages) == item.get('pages', 1), item
            combined.append(reader, import_outline=False)
        combined.add_outline_item(f'Figure {figure["number"]}', first - 1)
        legend_pages = [p['number'] - 1 for p in build['page_inventory']
                        if p.get('section_id') == f'figure{figure["number"]}']
        assert legend_pages, figure['number']
        for page in legend_pages:
            combined.add_page(legend_reader.pages[page])
            legends_only.add_page(legend_reader.pages[page])
        main_inventory.append(dict(figure=figure['number'], first_page=first,
                                   last_page=len(combined.pages), legend_pages=len(legend_pages)))
    combined.add_metadata({'/Title': 'Figures 2–4 with standalone legends',
                           '/Subject': 'Qualified scientific figure consolidation, 19 October 2026'})
    combined.write(out / 'figures_2_4_main.pdf')
    legends_only.add_metadata({'/Title': 'Figures 2–4: standalone legends'})
    legends_only.write(out / 'legends.pdf')
    supplementary = PdfWriter()
    supplement_inventory = []
    for figure in spec['figures']:
        for item in figure.get('supplementary_pdfs', []):
            first = len(supplementary.pages) + 1
            path = source(item['path'])
            reader = PdfReader(path)
            rebase_links(reader, path)
            assert len(reader.pages) == item['pages'], item
            supplementary.append(reader, import_outline=False)
            supplementary.add_outline_item(item['label'], first - 1)
            supplement_inventory.append(dict(label=item['label'], figure=figure['number'],
                                               first_page=first, last_page=len(supplementary.pages),
                                               source=item['path']))
    assert supplement_inventory
    supplementary.add_metadata({'/Title': 'Figures 2–4: comprehensive supplementary figures'})
    supplementary.write(out / 'figures_2_4_supplement.pdf')
    (out / 'LEGENDS.md').write_text('\n\n'.join(captions) + '\n')

    def href(relative):
        return html.escape(os.path.relpath(source(relative), out), quote=True)

    generated = datetime.now(timezone.utc).isoformat()
    page = ['<!doctype html><html lang="en"><head><meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width,initial-scale=1">',
            '<title>Figures 2–4: consolidated figure set</title>',
            '<style>body{font:17px/1.55 system-ui;max-width:1000px;margin:2rem auto;padding:0 1rem;color:#172d40}img{max-width:100%;height:auto}section{border-top:1px solid #ddd;padding:1.5rem 0}a{color:#126f79}.time{color:#536575}</style></head><body>',
            '<h1>Figures 2–4: consolidated figure set</h1>',
            f'<p class="time">Generated {html.escape(generated)}. Static final assembly; no timed refresh.</p>',
            '<p><a href="figures_2_4_main.pdf">Three main figures with legends</a> · '
            '<a href="figures_2_4_supplement.pdf">Comprehensive supplementary figures</a> · '
            '<a href="LEGENDS.md">Editable legends</a></p>']
    for figure in spec['figures']:
        page += [f'<section><h2>Figure {figure["number"]}: {html.escape(figure["label"])}</h2>',
                 f'<p>{html.escape(figure["summary"])}</p>',
                 f'<a href="{href(figure["preview"])}"><img src="{href(figure["preview"])}" alt="Figure {figure["number"]}"></a>',
                 f'<p><a href="{href(figure["caption"])}">Standalone caption</a></p>']
        for item in figure.get('supplementary_pdfs', []) + figure.get('appendix_pdfs', []):
            page.append(f'<p><a href="{href(item["path"])}">{html.escape(item["label"])}</a> ({item["pages"]} pages)</p>')
        for item in figure.get('legend_sources', []):
            page.append(f'<p><a href="{href(item["path"])}">{html.escape(item["label"])}</a></p>')
        page.append('</section>')
    page.append('</body></html>')
    (out / 'index.html').write_text('\n'.join(page) + '\n')
    record = dict(generated_utc=generated, source_manifest=str(args.manifest.resolve().relative_to(ROOT)),
                  inputs_sha256=inputs, main_page_inventory=main_inventory,
                  supplementary_page_inventory=supplement_inventory,
                  main_pages=len(combined.pages), supplementary_pages=len(supplementary.pages),
                  appendix_policy='Full native-profile appendices remain linked, labelled separate PDFs.',
                  native_integrations_launched=False, vector_main_plates_preserved=True,
                  rebased_local_source_annotations=rebased_links,
                  legends_renderer_validation=build['validation'])
    record['outputs_sha256'] = {str(p.relative_to(out)): sha(p) for p in [
        out / 'figures_2_4_main.pdf', out / 'figures_2_4_supplement.pdf',
        out / 'legends.pdf', out / 'LEGENDS.md', out / 'index.html']}
    write_json(out / 'supporting/metadata/assembly.json', record)
    print(json.dumps({k: v for k, v in record.items() if k not in ('inputs_sha256', 'outputs_sha256')}, indent=2))


if __name__ == '__main__':
    main()
