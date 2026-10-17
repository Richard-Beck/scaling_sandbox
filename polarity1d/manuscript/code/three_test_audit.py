"""Independent inputs and complete native-array audit of the three-test suite."""
from pathlib import Path
import copy, csv, json
import numpy as np
import audit_reversal as a


def inputs(run, manifest, package):
    provenance=a.read(run/'provenance.json')
    a.require(a.sha(run/'manifest.json')==provenance['manifest_sha256'],'Manifest changed')
    for name,digest in provenance['sources'].items():
        a.require(a.sha(package/name)==digest and a.sha(run/'frozen_sources'/name)==digest,'Frozen source changed: '+name)
    a.require(set(manifest['models'])==a.MODELS and set(manifest['excluded_models'])=={'otsuji','legi'},'Model scope')
    a.require(manifest['experiment']=='three_geometry_shared_preparation','Experiment identity')
    a.require(manifest['recovery_horizon']==10000. and manifest['sample_interval']==2.5 and manifest['gate_window']==500.,'Observation clocks')
    base=a.read(package/'data/runs/corrected_reversal_20261003/manifest.json')
    originals={(c['model'],c['length']):c for c in base['cases'] if c['condition']=='relative_20pct'}
    slopes=a.read(package/'data/runs/fixed_slope_reversal_20261006/manifest.json')['slope_by_model']
    grids=a.read(package/'data/inputs/reversal_duration_grids.json')['durations_by_model']
    cases=manifest['cases'];initial={}
    expected={(m,L,g) for m,L in originals for g in ['fixed_physical','fixed_relative','fixed_slope']}
    a.require(len(cases)==330 and {(c['model'],c['length'],c['condition']) for c in cases}==expected,'Full 5-model x 3-test coverage')
    a.require([c['case_index'] for c in cases]==list(range(330)),'Case indexing')
    for c in cases:
        old=originals[c['model'],c['length']]
        for key in ['checkpoint','checkpoint_sha256','initial_orientation','initial_contrast','eligible','establishment_case_index','cells','amplitude','amplitude_factor','fixed_scales']:
            a.require(c[key]==old[key],'Shared exact preparation changed: '+key)
        a.require(c['preparation_condition']=='relative_20pct','Different preparation')
        a.require(c['durations']==grids[c['model']],'Base duration coverage changed')
        factor=manifest['goryachev_selection']['selected_factor'] if c['model']=='goryachev' and c['condition']!='fixed_slope' else 1.
        a.close(c['response_factor'],factor,'Length-dependent patch factor')
        a.close(c['response_amplitude'],old['amplitude']*factor,'Length-dependent patch peak')
        if c['condition']=='fixed_slope':a.close(c['gradient_slope'],slopes[c['model']],'Slope changed')
        elif c['condition']=='fixed_physical':
            a.require(c['width_mode']=='physical','Fixed width mode')
            a.close(c['width'],.2 if c['model']=='spring' else 3.,'Fixed physical width changed')
        else:
            a.require(c['width_mode']=='relative' and c['width_fraction_at_reference']==.2,'Relative width changed')
            a.close(c['width'],.2*c['length'],'Nominal relative width')
        if not c['eligible']:continue
        checkpoint=package/c['checkpoint'];a.require(a.sha(checkpoint)==c['checkpoint_sha256'],'Checkpoint digest')
        with np.load(checkpoint) as z:
            a.close(float(z['orientation'][-1]),c['initial_orientation'],'Initial scalar')
            a.require(float(z['orientation'][-1])<-.05 and float(z['contrast'][-1])>=.05,'Eligibility')
            initial[c['case_index']]=a.checkpoint_fields(z,c)
    for path in (run/'work/chunks').glob('r0_*.json'):
        chunk=a.read(path);reuse=chunk.get('reuse_source')
        if not reuse:continue
        c=cases[chunk['case_index']]
        origin=package/reuse['campaign']/reuse['trace']
        a.require(a.sha(origin)==reuse['archive_sha256']==chunk['trace_sha256'],'Reused native archive changed')
        # Positive cached stimuli must be exactly unchanged. Zero is geometry-independent.
        if any(t['duration']>0 for t in chunk['trials']):
            a.require(c['condition']=='fixed_slope' or (c['condition']=='fixed_relative' and c['response_factor']==1.),'Changed stimulus reused')
        original_chunk=a.read(origin.with_suffix('.json'))
        a.require(original_chunk['checkpoint_sha256']==c['checkpoint_sha256'],'Cached checkpoint differs')
        oldrows={t['duration']:t for t in original_chunk['trials']}
        for row in chunk['trials']:
            old=copy.deepcopy(oldrows[row['duration']]);new=copy.deepcopy(row)
            old.pop('trace');new.pop('trace')
            a.require(old==new,'Cached result metadata changed')
    return cases,initial


a.audit_inputs=inputs
if __name__=='__main__':a.main()
