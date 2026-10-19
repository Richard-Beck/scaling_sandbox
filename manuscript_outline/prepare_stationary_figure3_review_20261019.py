#!/usr/bin/env python3
"""Compose the current stationary-cue evidence into a focused review document.

This prepares JSON for the existing scientific PDF renderer. It neither runs
native simulations nor replaces the integrated outline or preceding artwork.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019'


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def source(label: str, path: Path) -> dict:
    assert path.is_file(), path
    return {'label': label, 'path': relative(path)}


def asset(path: Path, title: str, caption: str) -> dict:
    with Image.open(path) as im:
        width, height = im.size
    return {
        'path': relative(path), 'page_title': title, 'caption': caption,
        'page_size': [11.8, min(24.0, max(8.0, 3.5 + 10.4 * height / width))],
        'minimum_printed_width_inches': 10.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--placement-contrast', type=Path, required=True)
    args = parser.parse_args()
    contrast = args.placement_contrast.resolve()
    sizes = PACKAGE / 'figure3/standard_angles_sizes'
    teaching = PACKAGE / 'figure3/standard_angles/vector_concept/teaching_example'
    summary = sizes / 'supporting/data/stationary_sizes_summary.csv'
    validation = json.loads((sizes / 'supporting/metadata/validation.json').read_text())
    assert validation['all_requested_complete'] is True
    assert validation['requested'] == validation['complete'] == 240
    generated = datetime.now(timezone.utc).isoformat()
    section = {
        'id': 'figure3',
        'title': 'Figure 3. Local activity and retained orientation under stationary cues',
        'role': (
            'One stationary-cue experiment links spatial activity to its global vector, '
            'then tests how cue placement, duration and model size affect response and retention. '
            'The following sheets are current evidence for choosing a cohesive main figure.'
        ),
        'assets': [
            asset(teaching / 'holmes_two_state_vector.png',
                  'Figure 3 review: constructing the activity vector',
                  'Measured Holmes profiles illustrate the physical vector sum and the distinction between resultant magnitude and total activity.'),
            asset(sizes / 'response/stationary_sizes_all_angles_overview.png',
                  'Figure 3 review: the complete stationary-cue experiment',
                  'All four models, three qualified sizes, four cue angles and five native durations: withdrawal and later orientation are paired.'),
            asset(contrast, 'Figure 3 review: a controlled placement contrast',
                  'The same prescribed stationary cue is placed at two positions on the same Goryachev basal state; local profiles and complex vectors show the measured responses.'),
        ],
        'legend': [
            'Vector construction. Each boundary segment contributes its real active concentration multiplied by physical segment length and its outward unit normal. The complex first moment stores the summed x and y components. Dividing by its initial complex value rotates and scales the initial vector to one. Its angle describes activity orientation; its magnitude R describes the resultant relative to its initial magnitude. Total active amount M0 is measured separately. In the measured Holmes C40 antipodal example, final total activity is 0.99746 of initial while R is 0.03993: two opposing peaks nearly cancel. Blue arrows group exact native cell contributions into twelve sectors using the same vector scale as the red sum.',
            'Experimental design. Every condition restarts an exact qualified basal state. Chemical circumferences are 20, 30 and 40; Spring resting circumferences are 1, 2 and 3. Cue centres lie at pi/8, pi/4, pi/2 or pi relative to the declared preparation axis. The frozen measured initial first-moment phase defines the small exact offset in plotted P coordinates; centres are never snapped to a different cell. Chemical cue durations are 1, 10, 100, 300 and 1000; Spring durations are 0.1, 1, 10, 100 and 1000. Each is followed by 13000 cue-free native time units.',
            'Complete response map. Rows identify models, columns cue angles and colours the three sizes. Open/dashed curves show withdrawal; filled/solid curves show the later endpoint. Actual angles require R at least 0.1. Weak-resultant endpoints are explicitly marked without assigning them a heading. The curve joins guide the eye between five tested durations. All 240 conditions and all 125 initial cue/core overlaps are retained. The fixed physical full support is 3 for chemical models and 0.2 for Spring; native peaks are 0.1, 6.4, 1 and 4 for Wave-pinning, Goryachev, Holmes and Spring. Clocks and doses are not calibrated across models.',
            'Controlled placement contrast. Goryachev C40 receives a 300-unit stationary cue at 45 or 90 degrees. Both conditions share the full basal state, kinetic parameters and cue shape, strength and duration. The cue arrays differ by an exact 100-cell cyclic shift; both have integral 9.6 and coefficient peak 6.4. These match the prescribed input coefficient; the native activation term also multiplies local inactive substrate, so equal activated biochemical amount is not imposed. Withdrawal headings are 29.59 and 33.48 degrees, but relative resultant magnitudes are 0.924 and 0.721 and both profiles contain two resolved peaks. After release, the 45-degree case progresses to 44.10 degrees while the 90-degree case returns to -0.028 degrees. The original material-region excess approaches 0 and 1, respectively. These are different complex vectors and local states at withdrawal.',
        ],
        'narrative': [
            'Ring measurements distinguish graded retained reorientation from nearly cancelling activity. Across all 180 non-antipodal conditions, the final field is a shifted copy of its own basal pattern with one resolved activity region and recovered resultant magnitude. Intermediate angles therefore reflect a retained spatial shift rather than a weak final activity pattern.',
            'The finite-cue response depends on size within each native protocol. For a 90-degree cue of duration 100, Wave-pinning retains 79.86, 55.36 and 32.52 degrees across circumferences 20, 30 and 40; Goryachev retains 89.09, -0.014 and -0.003 degrees. The complete map preserves other angles and durations rather than reducing them to an arbitrary effective-duration threshold.',
            'Cue-on deflection and later orientation are distinct observations. The same-basal placement pair demonstrates how similar withdrawal angles can accompany different resultant magnitudes, local profiles and subsequent activity redistribution. This supports examining the activity vector and spatial pattern together; it does not identify a unique competition mechanism.',
        ],
        'measurements': [
            'For every condition retain activity direction and resultant magnitude during the cue and release, total active amount, actual local profiles and the fixed original material-region excess. The rich tables include sampled minima/maxima, rebound, and first/retained threshold-entry brackets with finite-horizon censoring. Regional depletion can result from translation, redistribution or local chemistry and is not a site-ancestry measurement.',
            'A boundary activity recording could supply the same weighted sum using measured boundary elements and their actual outward normals. The circle provides those normals exactly in the model. Relating this descriptor to cell movement would require an additional empirical measurement; movement is not inferred here.',
        ],
        'limitations': [
            'The size series keeps physical cue support and native parameters fixed. Physical cue separation, angular footprint, overlap, basal architecture and material totals consequently change with size. These experiments measure that declared protocol response rather than an isolated causal contribution of geometry.',
            'The separate finite-tail field gate passes 228 conditions, is unavailable for 10 compact reused histories and fails for Holmes antipodal duration 300 at circumferences 30 and 40. Values and flags remain displayed. The second weak final resultant is Holmes C30/180 degrees/duration300 (R about 0.0997); the C40/180 degrees/duration1000 teaching example has R about 0.0399. A finite release does not prove asymptotic stability. Repeated samples from one deterministic history are not independent replicates.',
            'The figure composition remains under author discussion. These review sheets preserve the complete systematic evidence and a controlled illustrative pair; they are not an accepted final main-panel layout. The integrated outline PDF still contains the preceding prototype.',
        ],
        'equations': [
            {'math': r'M_1(t)=\sum_j u_j(t)\,\Delta s_j(t)\,[\cos\theta_j+i\sin\theta_j]',
             'number': '3.1',
             'description': 'The activity field u is real. The real and imaginary sums are Cartesian components of the physical activity vector. On each uniformly sampled circle, segment length is actual circumference divided by cell count, including Spring C(t).'},
            {'math': r'P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|,\qquad \phi(t)=\arg P(t)',
             'number': '3.2',
             'description': 'The initial normalization is applied once; it aligns the initial axis and sets the initial magnitude to one. Heading is displayed only above the declared relative-resultant gate.'},
            {'math': r'M_0=\sum_j u_j\Delta s_j,\qquad q_1=\frac{|M_1|}{M_0},\qquad R(t)=\frac{M_0(t)}{M_0(0)}\frac{q_1(t)}{q_1(0)}',
             'number': '3.3',
             'description': 'Relative resultant magnitude combines total activity and directional organization. It can exceed one and differs from current-mass-normalized coherence q1.'},
        ],
        'sources': [
            source('Authoritative complete 240-condition table', summary),
            source('Initial profiles, cue footprints and overlaps', sizes / 'supporting/data/stationary_sizes_initial_context.csv'),
            source('Regional crossing brackets and rebound flags', sizes / 'supporting/data/stationary_sizes_crossings.csv'),
            source('Native archive replay and complete-grid validation', sizes / 'supporting/metadata/validation.json'),
            source('Response-map evidence and enlarged supporting views', sizes / 'response/EVIDENCE.md'),
            source('Exact physical vector contributions in the teaching example', teaching / 'supporting/data/native_cell_contributions.csv'),
            source('Controlled placement-contrast evidence', contrast.parent / 'README.md'),
            source('Independent scientific review', PACKAGE / 'supporting/metadata/INDEPENDENT_REVIEW.md'),
        ],
    }
    document = {
        'title': 'Figure 3: one stationary-cue experiment',
        'subtitle': 'Local activity, complex vectors and retained reorientation',
        'date': '19 October 2026',
        'stage': 'Current evidence for author discussion; main composition open',
        'overview': [
            'The agreed experiment varies stationary cue position and finite duration across three qualified model sizes. This packet connects its design and measurements to the complete results, using an actual vector-construction example and a same-basal placement contrast.',
            'Read the first sheet to interpret the activity vector, the second to compare all 240 responses, and the third to see how different local and vector states evolve after cue withdrawal. All are views of the same stationary-cue experiment.',
            'The numerical experiment and native readouts are complete and independently checked. Main-panel selection and scientific emphasis remain under author discussion. Figure 2 rework and the final integrated manuscript presentation remain separate outstanding work.',
        ],
        'sections': [section], 'supplementary': [], 'references': [],
    }
    output = PACKAGE / 'figure3/supporting/metadata/stationary_review.json'
    output.write_text(json.dumps(document, indent=2, ensure_ascii=False) + '\n')
    record = {
        'generated_utc': generated, 'document': relative(output),
        'document_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'source_summary_sha256': hashlib.sha256(summary.read_bytes()).hexdigest(),
        'native_integrations': 0, 'changes_integrated_outline': False,
        'figure_assets': [a['path'] for a in section['assets']],
    }
    output.with_name('stationary_review_preparation.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
