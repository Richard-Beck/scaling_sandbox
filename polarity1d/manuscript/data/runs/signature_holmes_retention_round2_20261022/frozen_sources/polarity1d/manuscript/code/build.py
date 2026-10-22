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
    current = provenance.get('status') == 'three_test_shared_checkpoint_complete' and provenance.get('debelly_boundary_method') == 'independent_endpoints'
    final = current
    main_label = 'Main polarity figure — corrected five-model comparison' if final else 'Deprecated comparison — contains flawed centered deBelly results'
    image_alt = 'Corrected five-model protocol and measured length trends' if final else 'Historical comparison containing deprecated centered deBelly simulations'
    items = [
        ('Mechanisms of cell-length-dependent polarity reversal', 'length_mechanisms_20261008/index.html', 'length_mechanisms_20261008/length_mechanisms.pdf'),
        (main_label, 'main_polarity.png', 'main_polarity.pdf'),
        ('Length versus cue duration: five models × three tests', 'supplement_width_response.png', 'supplement_width_response.pdf'),
        ('Initial active material, width and concentration — corrected five-model starts', 'initial_polarity_length/initial_polarity_length.png', 'initial_polarity_length/initial_polarity_length_manuscript.pdf'),
        ('Goryachev length mechanism: weaker countercue activation and stronger incumbent restoring response', 'Goryachev_length_mechanism/index.html', 'Goryachev_length_mechanism/Goryachev_length_mechanism_report.pdf'),
        ('Corrected deBelly MCA zero crossing', 'deBelly_mechanism/debelly_zero_crossing.png', 'deBelly_mechanism/debelly_zero_crossing.pdf'),
        ('Corrected deBelly time to MCA zero crossing versus cell length', 'deBelly_mechanism/debelly_zero_crossing_vs_length.png', 'deBelly_mechanism/debelly_zero_crossing_vs_length.pdf'),
        ('Corrected deBelly advection mechanism', 'deBelly_mechanism/debelly_advection.png', 'deBelly_mechanism/debelly_advection.pdf'),
        ('Historical mechanisms (deBelly centered extension deprecated)', 'three_test_mechanisms/historical_index_20261006.html', 'three_test_mechanisms/historical_three_test_mechanisms_20261006.pdf'),
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
    for label, target in [('Three-test shared-checkpoint reversal', 'runs/debelly_corrected_three_test_20261008/results/status.json'),
                          ('Current length mechanism report', 'runs/length_mechanisms_20261008/report_provenance.json'),
                          ('Mechanism diagnostics', 'runs/three_test_mechanisms_20261006/results/status.json'),
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
<p><strong>deBelly warning:</strong> the old in-loop centered endpoint extension is flawed and deprecated. Do not use its simulations for scientific inference. Corrected figures explicitly identify independent local endpoint commands; historical reports retain deprecated deBelly results.</p>
<p>Five models: Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring.
All three finite challenges restart the same exact relative-20% basal polarized
checkpoint at each length: fixed width 3 (Spring native 0.2), relative width 20%,
and constant physical slope. Goryachev uses one constant stronger patch peak
across all lengths and both patch tests.</p>
<p><a href="../data/runs/debelly_corrected_three_test_20261008/results/index.html">Three-test campaign status and generation time</a></p>
<p><a href="../GRADIENT_REVERSAL.md">Added fixed-slope countergradient protocol</a> ·
<a href="../data/runs/fixed_slope_pilot_20261006/results/index.html">Slope calibration</a> ·
<a href="../data/runs/fixed_slope_reversal_20261006/results/index.html">Fixed-slope comparison status</a> ·
<a href="supplement_gradient_response/index.html">Fixed-slope supplement (after completion)</a></p>
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
        if (BASE / 'data/runs/debelly_corrected_three_test_20261008/results/validation.json').exists():
            run('three_test_figures.py')
        else:
            raise RuntimeError('Corrected deBelly campaign is not audited yet; do not rebuild canonical figures from deprecated centered-model evidence.')
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
