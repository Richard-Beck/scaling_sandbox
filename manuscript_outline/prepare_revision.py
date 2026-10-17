"""Assemble the critically revised outline without changing the first edition."""
import copy
import json
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads((ROOT / path).read_text())


def source(label, path):
    return dict(label=label, path=path)


def main():
    document = copy.deepcopy(load('manuscript_outline/outline_v1.json'))
    document['stage'] = 'Revised after independent scientific and editorial assessment'
    document['overview'][1] = (
        'Proposed argument: establish model- and stimulus-dependent size effects in finite polarity replacement; '
        'then connect specific contributions to local input response, incumbent restoration, transport, mechanics '
        'and site competition. Existing evolved interventions distinguish some causal contributions; the common '
        'gain/recovery assays provide comparable, mechanism-linked observables rather than a causal partition of '
        'all three duration curves. Supported rings test placement and sustained competition. Replicated noisy '
        'rings measure activity memory, with physical displacement, magnitude and observation coverage alongside '
        'angular descriptors. Published migration observations and conditional OU velocity-memory estimates '
        'supply separate empirical motivation. The new WGD motility assay remains a major pending contribution.'
    )
    document['overview'][2] = (
        'This is a proposed manuscript outline, not a completed manuscript. The order of the mechanism and '
        'existing-migration components remains an editorial choice. Six main figures are numbered for navigation, '
        'with continued sheets where needed for legibility. The second edition adds an all-five mechanism-to-'
        'measurement crosswalk, physical-motion and adaptation diagnostics, and the requested migration curve-'
        'shape display. Mechanism-linked cellular experiments are proposed, and Figure 6 remains wholly unmeasured.'
    )
    sections = {s['id']: s for s in document['sections']}
    crosswalk = load('manuscript_outline/mechanism_crosswalk_v2.json')
    linear = sections['figure2']
    linear['title'] = 'Figure 2. Incumbent organization, local response and model-specific mechanisms of competition'
    linear['assets'][0] = dict(
        path='manuscript_outline/mechanisms_numbered_v2/vision_aligned_figure1_linear.png',
        caption='All-five linear results, with the title renumbered for this outline; scientific pixels and numerical inputs are preserved.'
    )
    linear['assets'].extend(load('manuscript_outline/assets/mechanism_crosswalk_v2.assets.json'))
    for index, asset in enumerate(linear['assets'][1:]):
        asset['page_size'] = [11.8, 11.0 if index == 0 else 14.3]
        asset['page_title'] = 'Figure 2 continued: mechanisms, evidence and measurements'
        asset['minimum_printed_width_inches'] = 10.0
    linear['legend'].append('Continued mechanism-to-measurement sheet. ' + crosswalk['scope'])
    linear['legend'].extend(crosswalk['detailed_legend'])
    linear['narrative'].extend([
        'Mechanistic evidence is model-specific. Wave-pinning has an evolved fixed-dose inactive-mobility rescue; '
        'Holmes does not share that rescue at the tested cap. Spring has a cue-only extension-dependent '
        'deactivation intervention whose effect depends on stimulus geometry. These distinguish transport, '
        'spatial competition and mechanical feedback contributions without defining new counterfactual thresholds.',
        'Corrected deBelly has nonmonotone duration curves, including decreases at L18–21 and L31–33. '
        'Within selected rising ranges, endpoint MCA depletion and common-time accumulation help explain '
        'delayed incumbent Rac decline. L39 requires slower refill in addition to its starting state. A separate '
        'resting-state tension-to-Rho ablation reduces displacement spread by 50.4%; that causal result concerns '
        'resting organization, not an assigned fraction of reversal duration.',
        'Goryachev links lower available inactive concentration to weaker challenger activation and stronger '
        'incumbent repair under specified perturbations. The common active-loss assay, older inactive-bolus '
        'experiment and activation-suppression assay have different clocks and recovery tails and remain separate.'
    ])
    linear['measurements'].extend(
        row['model'] + ': ' + row['measurement'] for row in crosswalk['rows']
    )
    linear['limitations'].append(
        'The crosswalk retains evolved interventions, fixed-state derivatives and temporal/budget associations '
        'at their respective inference levels. It does not assign additive causal fractions to full duration '
        'curves or identify a unique biological model from reporter measurements.'
    )
    seen = {s['path'] for s in linear['sources']}
    linear['sources'].append(source('Mechanism, evidence and corresponding measurement crosswalk',
                                    'manuscript_outline/mechanism_crosswalk_v2.json'))
    for row in crosswalk['rows']:
        for entry in row['sources']:
            if entry['path'] not in seen:
                linear['sources'].append(entry)
                seen.add(entry['path'])

    ring = sections['figure3']
    for asset in ring['assets']:
        asset['path'] = 'manuscript_outline/mechanisms_numbered_v2/' + Path(asset['path']).name
        asset['caption'] = asset['caption'].replace('Delivered checkpoint Figure 2', 'Figure 3')
        asset['caption'] += ' Only the working title numeral is changed from the delivered artwork.'

    noise = sections['figure4']
    revised_noise = load('manuscript_outline/assets/noise_v2/noise_interpretation_v2.content.json')
    noise['title'] = 'Figure 4. Activity memory, physical site motion and forced-state qualification on noisy rings'
    noise['assets'][0] = dict(
        path='manuscript_outline/mechanisms_numbered_v2/vision_aligned_figure3_noise.png',
        caption='Figure 4 memory/time-course sheet: original numerical results with the working title renumbered.'
    )
    noise['assets'].extend(revised_noise['assets'])
    noise['legend'] = ['Memory/time-course sheet: ' + p for p in noise['legend']]
    noise['legend'].extend('Diagnostic/continued sheets: ' + p for p in revised_noise['legend'])
    noise['narrative'] = revised_noise['narrative']
    noise['measurements'] = revised_noise['measurements']
    noise['limitations'] = revised_noise['limitations'] + [
        'Activity-vector memory, conditional physical site displacement and conditional OU velocity memory '
        'remain distinct. No cell-velocity law or calibrated migration-persistence time follows from these rings.'
    ]
    noise['sources'].extend(revised_noise['sources'][2:])
    noise['equations'].append(dict(
        math=r'\Delta\theta=\frac{2\pi\,\Delta s}{C}', number='8',
        description='For the same boundary arclength displacement, angular change decreases with circumference. '
                    'The measured shortest endpoint displacement is conditional on the stated single-site gate; '
                    'it is not an unwrapped trajectory or a cell displacement.'
    ))

    migration = sections['figure5']
    revised_migration = load('manuscript_outline/assets/migration_v2/outline_integration.json')
    migration['assets'] = revised_migration['assets']
    migration['legend'][revised_migration['replace_legend_index']] = (
        revised_migration['replacement_legend'].replace('2.32014 h', 'approximately 2.32 h')
        .replace('1.11406 h', 'approximately 1.11 h')
    )
    migration['narrative'].extend([
        revised_migration['additional_narrative'][0],
        'Lower sampled speed can coexist with greater late-lag displacement in the enlarged-group comparison. '
        'Its causal connection to polarity mechanisms requires matched activity and migration measurements.'
    ])
    migration['limitations'][1] = (
        'The requested observed-point/model-line and curve-shape display is implemented. The original '
        'experimental lag grid remains unknown; neither display samples nor raster-derived slopes supply '
        'independent biological replicates. The observed/model crossover discrepancy remains visible.'
    )
    migration['sources'].extend([
        source('Revised display, exact preserved inputs and direct-user feedback',
               'manuscript_outline/assets/migration_v2/provenance.json'),
        source('Migration display numerical and preservation checks',
               'manuscript_outline/assets/migration_v2/validation.json')
    ])
    document['supplementary'][2]['text'] += (
        ' Revised Figure 4 also displays both strengths, complete half-window diagnostics, physical-site '
        'eligibility and matched C20/C60 native examples; compact source tables support portable rerendering.'
    )
    document['references'].extend([
        source('Independent critical assessment and repair decisions',
               'reports/manuscript_coherence_review_v1_20261017.md'),
        source('Portable original/renumbered mechanism replay and source closure',
               'manuscript_outline/PORTABLE_MECHANISM_RENDER.md'),
        source('Second-edition integration decisions and validation scope',
               'manuscript_outline/REVISION_NOTES.md')
    ])
    # Give dense scientific panels a readable digital foldout rather than
    # scaling them down to the ordinary prose-page width. Continued captions
    # and all detailed legends stay separate from the scientific pixels.
    short_titles = {
        'figure1': 'Figure 1. Five models and three finite challenges',
        'figure2': 'Figure 2. Linear response and mechanistic evidence',
        'figure3': 'Figure 3. Supported-ring competition',
        'figure4': 'Figure 4. Noisy-ring activity memory and physical motion',
        'figure5': 'Figure 5. Migration observations and conditional OU interpretation',
        'figure6': 'Figure 6. RESERVED: new WGD motility results',
    }
    for section in document['sections'] + document['supplementary']:
        for asset in section.get('assets', []):
            with Image.open(ROOT / asset['path']) as im:
                width, height = im.size
            asset.setdefault('page_size', [11.8, max(8.27, round(10.4 * height / width + 3.0, 2))])
            asset.setdefault('page_title', short_titles.get(section['id'], section['title']))
            asset.setdefault('minimum_printed_width_inches', 10.0)
    destination = ROOT / 'manuscript_outline/outline_v2.json'
    destination.write_text(json.dumps(document, ensure_ascii=False, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
