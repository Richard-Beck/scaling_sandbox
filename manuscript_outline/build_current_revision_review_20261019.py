#!/usr/bin/env python3
"""Build the dated current index from the final, explicitly selected figure set."""
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
import hashlib
import html
import json
import os

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019'
OUTPUT = PACKAGE / 'CURRENT_REVIEW.html'
INPUTS = {}


def href(path):
    path = Path(path).resolve()
    assert path.is_file(), path
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return quote(os.path.relpath(path, OUTPUT.parent), safe='/')


def link(label, path):
    return f'<a href="{href(path)}">{html.escape(label)}</a>'


def main():
    generated = datetime.now(timezone.utc).isoformat()
    final = PACKAGE / 'consolidated_figure_set_20261019'
    spec = json.loads((final / 'supporting/metadata/selection.json').read_text())
    page = ['<!doctype html><html lang="en"><head><meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width,initial-scale=1">',
            '<title>Figures 2–4 · current consolidated set</title>',
            '<style>body{font:17px/1.55 system-ui;color:#172d40;max-width:1000px;margin:2rem auto;padding:0 1rem}section{border-top:1px solid #ddd;padding:1rem 0}img{max-width:760px;width:100%;height:auto}a{color:#126f79}.muted{color:#536575}</style></head><body>',
            '<h1>Figures 2–4: consolidated set</h1>',
            f'<p class="muted">Generated {html.escape(generated)}. All scientific producers are terminal. Static snapshot; no timed refresh.</p>',
            '<p>' + ' · '.join([link('Main figures and legends', final / 'figures_2_4_main.pdf'),
                               link('64-page comprehensive supplement', final / 'figures_2_4_supplement.pdf'),
                               link('Reader, methods and profile appendices', final / 'index.html')]) + '</p>',
            '<p>' + link('Selection, interpretation and limits', PACKAGE / 'BEST_EFFORT_FIGURES_20261019.md') + ' · ' + link('Current scientific record', PACKAGE / 'RESULTS.md') + '</p>',
            '<p>All three selected main figures and complete experimental supports have passed independent scientific/source/display checks. Final assembled-package, storage and repository delivery are recorded in ' + link('WORK_STATUS.json', PACKAGE / 'supporting/metadata/WORK_STATUS.json') + '. Earlier prototypes remain dated development history.</p>']
    for figure in spec['figures']:
        page += [f'<section id="figure{figure["number"]}"><h2>Figure {figure["number"]}: {html.escape(figure["label"])}</h2>',
                 f'<p>{html.escape(figure["summary"])}</p>',
                 f'<a href="{href(ROOT / figure["preview"])}"><img src="{href(ROOT / figure["preview"])}" alt="Figure {figure["number"]}"></a>',
                 '<p>' + link('Standalone caption', ROOT / figure['caption']) + '</p></section>']
    page += ['<p>Figures 1, 5 and reserved Figure 6 retain their preceding content. This is a working manuscript figure set, not a completed manuscript or author acceptance.</p>',
             '<p class="muted">Preceding compositions: ' + link('outline v3', ROOT / 'manuscript_outline/outline_v3.pdf') + ' · ' + link('earlier combined review', PACKAGE / 'supporting/vectors/revised_figures_2_4.pdf') + '.</p></body></html>']
    OUTPUT.write_text('\n'.join(page) + '\n')
    result = dict(generated_utc=generated, output=str(OUTPUT.relative_to(ROOT)),
                  output_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
                  inputs_sha256=INPUTS, timed_browser_refresh=False,
                  scope='Final consolidated scientific set; explicit selection and input identities. No integrations launched.')
    destination = PACKAGE / 'supporting/metadata/current_review_entrypoint.json'
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'inputs_sha256'}, indent=2))


if __name__ == '__main__':
    main()
