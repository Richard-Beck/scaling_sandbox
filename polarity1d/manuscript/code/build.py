"""Single entry point for the retained manuscript artifacts; never submits jobs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import subprocess
import sys

BASE = Path(__file__).resolve().parents[1]
CODE = BASE / 'code'
FIGURES = BASE / 'figures'


def run(script, *arguments):
    subprocess.run([sys.executable, str(CODE / script), *arguments], check=True)


def index():
    FIGURES.mkdir(exist_ok=True)
    provenance_path = BASE / 'data/derived/figure_generation_provenance.json'
    provenance = json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    final = provenance.get('status') == 'corrected_basal_matched_reversal_complete'
    main_label = 'Main polarity figure — corrected five-model comparison' if final else 'Corrected protocol preview — response panels pending'
    image_alt = 'Corrected five-model protocol and measured length trends' if final else 'Corrected protocol only; reversal response panels pending'
    items = [
        (main_label, 'main_polarity.png', 'main_polarity.pdf'),
        ('Detailed width / response-definition supplement', 'supplement_width_response.png', 'supplement_width_response.pdf'),
        ('Corrected basal-withdrawal establishment supplement', 'supplement_establishment/index.html', 'supplement_establishment/polarity_establishment_by_model.pdf'),
        ('Interactive reversal fields — personal diagnostic reference', 'diagnostic_reversal_fields/index.html', None),
    ]
    links = []
    for title, target, pdf in items:
        if not (FIGURES / target).is_file():
            continue
        extra = f' · <a href="{pdf}">PDF</a>' if pdf and (FIGURES / pdf).is_file() else ''
        links.append(f'<li><a href="{target}">{html.escape(title)}</a>{extra}</li>')
    stamp = datetime.now(timezone.utc).isoformat()
    sources = {}
    for label, target in [('Corrected five-model reversal', 'runs/corrected_reversal_20261003/results/status.json'),
                          ('Corrected establishment', 'establishment_basal/results/status.json')]:
        path = BASE / 'data' / target
        if path.exists():
            status = json.loads(path.read_text())
            sources[label] = status.get('generated_at', status.get('completed_at', 'See retained status file'))
    source_list = ''.join(f'<li>{html.escape(k)}: {html.escape(str(v))}</li>' for k, v in sources.items())
    (FIGURES / 'index.html').write_text(f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Polarity manuscript figures</title>
<style>body{{max-width:1000px;margin:40px auto;padding:0 24px;font:17px/1.6 system-ui;color:#20252b}}a{{color:#12628c}}img{{max-width:100%}}small{{color:#58626b}}</style>
<h1>Polarity manuscript figures</h1>
<p>Five models: Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring.
The finite right countercue mirrors its corrected finite-patch preparation.
The final main and detailed reversal supplement use the fresh campaign only.</p>
<p><a href="../data/runs/corrected_reversal_20261003/results/index.html">Fresh reversal campaign status and generation time</a></p>
<ul>{''.join(links)}</ul>
<p><a href="../PROTOCOL.md">Protocol distinctions</a> · <a href="../ASSUMPTIONS.md">Assumptions</a> ·
<a href="../RETENTION.md">Retention and deprecation</a> · <a href="CAPTIONS.md">Figure captions</a></p>
<img src="main_polarity.png" alt="{html.escape(image_alt)}">
<p><small>Index generated {stamp}. The index does not submit jobs; producer status and data time appear on the linked campaign page. Refresh manually.</small></p>
<ul>{source_list}</ul></html>''')
    print(f'Wrote {FIGURES / "index.html"}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['render', 'index', 'verify', 'diagnostic'])
    args = parser.parse_args()
    if args.action == 'render':
        run('figure_builder.py')
        index()
    elif args.action == 'index':
        index()
    elif args.action == 'verify':
        run('check_codeset.py' if (BASE / 'CODESET_SHA256.json').exists() else 'verify.py')
    elif args.action == 'diagnostic':
        run('reversal_report/build.py', '--bundle')
        index()


if __name__ == '__main__':
    main()
