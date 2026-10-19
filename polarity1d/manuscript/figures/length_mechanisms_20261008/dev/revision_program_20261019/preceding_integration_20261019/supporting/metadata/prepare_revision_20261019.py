#!/usr/bin/env python3
"""Integrate the 19 October Figures 2–4 revisions without changing other figures.

Run with the existing scaling_sandbox Python, after the three figure exporters.
This composes retained assets and reviewed figure content; it runs no simulation.
"""
from __future__ import annotations

import copy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019'
OUTLINE = ROOT / 'manuscript_outline'


def load(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return str(path.resolve().relative_to(ROOT))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def asset_sizes(section):
    for asset in section.get('assets', []):
        path = ROOT / asset['path']
        asset['path'] = rel(path)
        with Image.open(path) as image:
            width, height = image.size
        # Match the existing outline's ten-inch printed figure width and keep
        # the entire image visible on an explicitly sized review sheet.
        # Reserve the heading, short source caption and footer in addition to
        # the image. The renderer independently enforces the printed width.
        asset['page_size'] = [11.8, min(24.0, max(8.0, 3.5 + 10.4 * height / width))]
        asset['minimum_printed_width_inches'] = 10.0
        asset.setdefault('page_title', section['title'])


def main():
    baseline_path = PACKAGE / 'supporting/metadata/integration_baseline.json'
    baseline = load(baseline_path)
    for name, digest in baseline['sha256'].items():
        assert sha(ROOT / name) == digest, ('preserved input changed', name)
    old = load(OUTLINE / 'outline_v2.json')
    document = copy.deepcopy(old)
    document.update(date='19 October 2026', stage='Figures 2–4 revised against the 19 October program')
    document['overview'] = [
        old['overview'][0],
        'Figures 2–4 focus on three questions: which measurements separate explanations for difficult '
        'polarity replacement; how local activity becomes global polarity on a ring; and how physical '
        'site motion, pattern organization and geometry contribute to angular memory. Retained '
        'simulations and targeted distance–duration and Spring-ring experiments supply the evidence. '
        'Detailed assay maps and secondary diagnostics are supporting material.',
        'This is a working manuscript outline. Figures 1 and 5 and the reserved Figure 6 experimental '
        'section retain their preceding content. Activity polarity, modelled replacement, observed '
        'migration and conditional OU velocity memory remain distinct readouts; their biological '
        'connection is an experimental hypothesis.'
    ]
    content_paths = []
    additions = []
    for number in (2, 3, 4):
        path = PACKAGE / f'figure{number}/supporting/metadata/content.json'
        content = load(path)
        content_paths.append(path)
        section = copy.deepcopy(content.get('section', content))
        if number == 3:
            continued_path = PACKAGE / 'figure3/distance_response/supporting/metadata/content_continued.json'
            continued = load(continued_path)
            content_paths.append(continued_path)
            for field in ('assets', 'legend', 'narrative', 'sources', 'supplementary'):
                section.setdefault(field, []).extend(continued.get(field, []))
            section['title'] = 'Figure 3. Geometry, cue placement and duration shape retained polarity'
            section['role'] = (
                'Closed geometry connects local activity to global direction and magnitude; '
                'distance and duration alter retained orientations and regional activity.'
            )
            section['measurements'].append(
                'Measure retained angle, magnitude and local-region activity across characterized '
                'cue distances and durations; report angle alignment and regional depletion as separate endpoints.'
            )
            section['limitations'].append(
                'The duration grids identify shortest tested passes under explicit criteria; '
                'adjacent tested failures and threshold sensitivities remain available, and '
                'no optimized minimum duration is claimed. Spring retains its declared '
                'mean-force global-extension closure and native stretching, dilution and diffusion.'
            )
            section['equations'] = [
                {
                    'math': r'M_k(t)=\int_0^{C(t)} u(s,t)\exp\!\left(\frac{2\pi i k s}{C(t)}\right)\,\mathrm{d}s',
                    'number': '3.1',
                    'description': 'The native active-field angular moments use physical arc quadrature. '
                                   'M₀ is the instantaneous active amount; C(t) is the actual circumference, '
                                   'including Spring expansion and contraction.',
                },
                {
                    'math': r'P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|,\qquad q_k(t)=\frac{|M_k(t)|}{M_0(t)}',
                    'number': '3.2',
                    'description': 'P measures the first moment relative to its initial complex value. '
                                   'The shape moments q₁ and q₂ use instantaneous total active amount; '
                                   'they are distinct from R. Direction requires the declared magnitude gate.',
                },
                {
                    'math': r'\Delta\theta_0=\frac{2\pi d_0}{C_b}',
                    'number': '3.3',
                    'description': 'Prescribed initial-peak-to-cue separation d₀ is mapped to the '
                                   'normal-angle difference using the basal circumference Cb. '
                                   'This cue-placement quantity differs from subsequent activity-site relocation.',
                },
            ]
            section['sources'].extend([
                {'label': 'Spring basal organization, input gain and local recovery',
                 'path': 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/spring_function/README.md'},
                {'label': 'Spring common active-loss response at exact challenged t10 states',
                 'path': 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/spring_challenged_function/README.md'},
                {'label': 'Declared Spring ring model and numerical validation',
                 'path': 'polarity1d/manuscript/SPRING_RING_STATUS_20261019.md'},
            ])
        additions.extend(section.pop('supplementary', []))
        section['id'] = f'figure{number}'
        for field in ('title', 'role', 'assets', 'legend', 'narrative', 'measurements', 'limitations', 'sources'):
            assert section.get(field), (path, field)
        # Exact content identities are in integration.json. The figure notes
        # link the scientific methods and data directly, avoiding a redundant
        # workflow-source link that can create a nearly empty continuation.
        asset_sizes(section)
        document['sections'][number - 1] = section
    for number in (1, 5, 6):
        assert document['sections'][number - 1] == baseline['protected_sections'][f'figure{number}']
    document['supplementary'].extend(additions)
    # The earlier complete deterministic/noisy comparison remains available
    # without carrying its model catalogue into the revised main figures.
    document['supplementary'].append({
        'title': 'Preserved complete mechanism and assay checkpoint',
        'text': [
            'The preceding outline retains all-five mechanism-to-measurement sheets, ring gain/recovery '
            'and finite-cue maps, and noise magnitude, controls and full/first/last-window diagnostics. '
            'Its numerical evidence is preserved; the revised main-figure claims are stated in this edition.'
        ],
        'path': 'manuscript_outline/outline_v2.pdf',
        'sources': [
            {'label': 'Prior all-five mechanism evidence crosswalk', 'path': 'manuscript_outline/mechanism_crosswalk_v2.json'},
            {'label': 'Prior complete ring and noise checkpoint', 'path': 'current manscript figs/vision_aligned_checkpoint.md'},
            {'label': 'Prior physical motion, magnitude and coverage evidence', 'path': 'manuscript_outline/assets/noise_v2/noise_interpretation_v2.content.json'},
        ],
    })
    for item in document['supplementary']:
        if isinstance(item, dict) and item.get('assets'):
            asset_sizes(item)
    document['references'].append({'label': 'Figures 2–4 revision specification', 'path': 'docs/tasks/29_10_19_revision_program.md'})
    write_json(OUTLINE / 'outline_v3.json', document)
    review = copy.deepcopy(document)
    review['title'] = 'Revised Figures 2–4: mechanisms, geometry and polarity dynamics'
    review['subtitle'] = 'Main figures, legends and selected supporting analyses'
    review['sections'] = document['sections'][1:4]
    # Keep this review PDF focused; the full six-figure outline contains the
    # unchanged older supplementary sections.
    review['supplementary'] = additions + [document['supplementary'][-1]]
    write_json(PACKAGE / 'supporting/metadata/revision_review.json', review)
    record = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'source_outline': 'manuscript_outline/outline_v2.json',
        'source_outline_sha256': sha(OUTLINE / 'outline_v2.json'),
        'revised_outline': 'manuscript_outline/outline_v3.json',
        'revised_outline_sha256': sha(OUTLINE / 'outline_v3.json'),
        'figure_content_sha256': {rel(path): sha(path) for path in content_paths},
        'unchanged_sections_verified': ['figure1', 'figure5', 'figure6'],
        'preserved_input_hashes_verified': len(baseline['sha256']),
        'integration_launches_simulations': False,
        'scientific_inputs': 'Retained simulations and separately frozen targeted native campaigns',
    }
    write_json(PACKAGE / 'supporting/metadata/integration.json', record)

    viewer = PACKAGE / 'supporting/viewers/index.html'
    viewer.parent.mkdir(parents=True, exist_ok=True)
    href = lambda path: html.escape(os.path.relpath(path, viewer.parent), quote=True)
    text = ['<!doctype html><html lang="en"><meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            '<title>Revised Figures 2–4</title>',
            '<style>body{max-width:1100px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#172d40}img{width:100%;height:auto}h1,h2{line-height:1.2}a{color:#146d83}.time{color:#536575}section{margin:3rem 0}</style>',
            '<h1>Revised Figures 2–4</h1>',
            '<p><strong>Prototype under revision:</strong> Figure 2 requires substantial reworking after author feedback. The displayed sheet and linked PDF retain the preceding prototype. Figures 3–4 are generated drafts for review; numerical checks do not imply author acceptance.</p>',
            f'<p class="time">Generated {html.escape(record["generated_utc"])}. Static report; native campaigns are separately recorded scientific inputs.</p>',
            f'<p><a href="{href(PACKAGE / "supporting/vectors/revised_figures_2_4.pdf")}">Figures and legends PDF</a> · '
            f'<a href="{href(OUTLINE / "outline_v3.pdf")}">Full revised outline</a></p>']
    for section in document['sections'][1:4]:
        text.append(f'<section><h2>{html.escape(section["title"])}</h2>')
        for paragraph in section['narrative']:
            text.append(f'<p>{html.escape(paragraph)}</p>')
        for asset in section['assets']:
            text.append(f'<a href="{href(ROOT / asset["path"])}"><img src="{href(ROOT / asset["path"])}" alt="{html.escape(asset.get("caption", section["title"]), quote=True)}"></a>')
        text.append('</section>')
    text.append(f'<p><a href="{href(PACKAGE / "supporting/viewers/distance_response.html")}">Interactive distance-response explorer</a> · '
                f'<a href="{href(PACKAGE / "supporting/metadata/INDEPENDENT_REVIEW.md")}">Independent scientific review</a> · '
                f'<a href="{href(PACKAGE / "RESULTS.md")}">Methods, evidence and limitations</a></p></html>')
    viewer.write_text('\n'.join(text) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
