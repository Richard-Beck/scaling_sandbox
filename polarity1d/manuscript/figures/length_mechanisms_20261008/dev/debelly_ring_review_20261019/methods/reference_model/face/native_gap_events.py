"""Fixed-at-birth Rac-gap prototype events; declared component duplicate guard.

One connected Rac>threshold region opens at its highest local maximum if it
contains no existing gap. This is the user's latest chosen birth rule.
Multiple local maxima within the same occupied region are logged separately.
Existing gaps survive component merger until their own face Rac falls below the
threshold. Numerical split faces are retained as ordinary faces after healing.
"""
from __future__ import annotations
from copy import deepcopy
import numpy as np

from peak_detector import periodic_maxima
from ring_core import FIELDS, face_trace, split_cell


def geometry(widths):
    edges = np.r_[0., np.cumsum(widths)]
    return edges, .5 * (edges[:-1] + edges[1:])


def active_components(values, threshold):
    mask = np.asarray(values) > threshold
    size = len(mask)
    if mask.all():
        return [list(range(size))]
    if not mask.any():
        return []
    start = int(np.flatnonzero(~mask)[0])
    groups, current = [], []
    for j in range(1, size + 1):
        index = (start + j) % size
        if mask[index]:
            current.append(index)
        elif current:
            groups.append(current)
            current = []
    return groups


def inventory(widths, state, threshold):
    edges, centers = geometry(widths)
    peaks = periodic_maxima(state['R'], edges[-1], threshold, centers)
    if peaks['status'] in ('roundoff_group_ambiguous', 'ambiguous_transitive_roundoff_plateau'):
        raise RuntimeError('Peak plateau comparison numerically ambiguous; no event convention substituted')
    components = active_components(state['R'], threshold)
    owner = {index: i for i, group in enumerate(components) for index in group}
    records = [dict(cells=group, peaks=[], existing_gap_ids=[]) for group in components]
    for peak in peaks['peaks']:
        if peak['above_threshold'] == 'roundoff_ambiguous':
            raise RuntimeError('Peak at activity threshold is numerically ambiguous')
        if peak['above_threshold'] is True:
            component = owner[peak['sample_indices'][0]]
            records[component]['peaks'].append(peak)
    gaps = []
    for face, gap in state.get('cuts', {}).items():
        value = face_trace(state['R'], widths, face)
        component = None
        adjacent = [(face - 1) % len(widths), face]
        if value >= threshold:
            candidates = [owner[i] for i in adjacent if i in owner]
            if candidates:
                component = candidates[0]
                records[component]['existing_gap_ids'].append(int(gap['id']))
        local_peak = None
        if component is not None and records[component]['peaks']:
            local_peak = max(records[component]['peaks'], key=lambda p: p['maximum'])
        drift = None if local_peak is None else (
            (local_peak['position'] - gap['x'] + .5 * edges[-1]) % edges[-1] - .5 * edges[-1])
        gaps.append(dict(id=int(gap['id']), x=float(gap['x']), face=int(face),
                         Rac_trace=value, component=component,
                         associated_component_peak=None if local_peak is None else local_peak['position'],
                         signed_peak_gap_offset=drift))
    return dict(peaks=peaks['peaks'], components=records, gaps=gaps,
                maximum_Rac=float(np.max(state['R'])), minimum_Rac=float(np.min(state['R'])))


def transition(widths, state, time, threshold, next_id):
    """Close by local trace first; then initiate uncovered active regions."""
    widths = np.asarray(widths, dtype=float)
    state = {**state, 'cuts': deepcopy(state.get('cuts', {}))}
    mass_before = float(np.dot(widths, state['b'] + state['c']))
    records = []
    for face, gap in list(state['cuts'].items()):
        value = face_trace(state['R'], widths, face)
        if value < threshold:
            records.append(dict(event='close', time=float(time), id=int(gap['id']),
                                x=float(gap['x']), Rac_trace=value,
                                old_Uleft=float(gap['Uleft']), old_Uright=float(gap['Uright'])))
            del state['cuts'][face]
    info = inventory(widths, state, threshold)
    births = []
    suppressed = []
    for component in info['components']:
        candidates = sorted(component['peaks'], key=lambda p: (-p['maximum'], p['position']))
        if not candidates:
            continue
        if component['existing_gap_ids']:
            # A drifting maximum still belongs to the same connected region.
            if len(candidates) > 1:
                suppressed.extend(dict(position=p['position'], value=p['maximum'],
                                       reason='distinct maxima share active region containing existing gap')
                                  for p in candidates[1:])
        else:
            tolerance=64*np.finfo(float).eps*max(1.,abs(candidates[0]['maximum']))
            if len(candidates)>1 and abs(candidates[0]['maximum']-candidates[1]['maximum'])<=tolerance:
                raise RuntimeError('Equal distinct highest maxima in gap-free region; no rotationally biased birth selected')
            births.append(candidates[0])
            suppressed.extend(dict(position=p['position'], value=p['maximum'],
                                   reason='only highest maximum in gap-free active region selected')
                              for p in candidates[1:])
    for peak in sorted(births, key=lambda p: p['position']):
        edges, _ = geometry(widths)
        circumference = edges[-1]
        x = float(peak['position'] % circumference)
        distances = abs(edges[:-1] - x)
        distances = np.minimum(distances, circumference-distances)
        face = int(np.argmin(distances))
        tolerance = 64 * np.finfo(float).eps * max(1., circumference)
        if distances[face] > tolerance:
            cell = int(np.searchsorted(edges, x, side='right') - 1)
            fraction = (x - edges[cell]) / widths[cell]
            widths, state, face = split_cell(widths, state, cell, fraction)
        else:
            x=float(edges[face])
        if face in state['cuts']:
            raise RuntimeError('Birth would duplicate an existing physical point')
        displacement = face_trace(state['u'], widths, face)
        current_trace = face_trace(state['R'], widths, face)
        state['cuts'][face] = dict(id=int(next_id), x=x, birth_time=float(time),
                                  Uleft=displacement, Uright=displacement,
                                  raw_left=displacement, raw_right=displacement)
        records.append(dict(event='open', time=float(time), id=int(next_id), x=x,
                            detected_peak_Rac=float(peak['maximum']),
                            inserted_face_Rac=current_trace, initialized_U=displacement,
                            removed_material=0., cells_after=len(widths)))
        next_id += 1
    mass_after = float(np.dot(widths, state['b'] + state['c']))
    if abs(mass_after - mass_before) > 1e-11 * max(1., mass_before):
        raise RuntimeError('Topology event changed MCA material')
    if len(widths) > 2000 or np.min(widths) < 1e-6 * (62. / 402):
        raise RuntimeError('Numerical mesh guard reached; pilot stops rather than merging nearby gaps')
    return widths, state, next_id, dict(events=records,
        distinct_maximum_suppressions=suppressed, material_impulse=mass_after-mass_before,
        prebirth_inventory=info)


def event_preflight():
    widths = np.ones(8)
    state = {q: np.zeros(8) for q in FIELDS}
    state.update(R=np.array([.5,.6,.9,.6,.5,.6,.7,.6]),
                 Q=np.full(8,.7), b=np.full(8,.07), c=np.full(8,.03),cuts={})
    w, a, ident, result = transition(widths,state,0.,.8,1)
    assert len(a['cuts'])==1 and len(w)==9 and result['events'][0]['x']==2.5
    assert abs(np.dot(w,a['b']+a['c'])-.8)<1e-14
    original_x=next(iter(a['cuts'].values()))['x']
    # Put a shifted maximum into the same contiguous active region: no motion/birth.
    a['R']=np.array([.5,.6,.85,.88,.95,.85,.7,.6,.5])
    w2,b,ident2,result2=transition(w,a,1.,.8,ident)
    assert len(b['cuts'])==1 and ident2==ident
    assert next(iter(b['cuts'].values()))['x']==original_x
    # Two maxima within that active region: expose omitted literal peak birth.
    b['R']=np.array([.5,.6,.85,.95,.85,.96,.85,.6,.5])
    _,b2,_,result3=transition(w2,b,2.,.8,ident2)
    assert result3['distinct_maximum_suppressions']
    # Local loss closes a fixed point; there is no imposed minimum lifetime.
    b2['R']=np.full(len(w2),.5)
    _,c,_,result4=transition(w2,b2,3.,.8,ident2)
    assert len(c['cuts'])==0 and result4['events'][0]['event']=='close'
    # Two previously existing gaps may share a merged active component.
    merged={q:np.full(8,.1) for q in FIELDS}
    merged.update(R=np.array([.81,.83,.86,.9,.95,.92,.89,.85]),cuts={
        1:dict(id=1,x=1.,Uleft=0.,Uright=0.),
        6:dict(id=2,x=6.,Uleft=0.,Uright=0.)})
    _, kept, _, last=transition(np.ones(8),merged,0.,.8,3)
    assert len(kept['cuts'])==2 and not last['events']
    # Latest user rule: two maxima within ONE active region create ONE break.
    shared={q:np.full(5,.1) for q in FIELDS}
    shared.update(R=np.array([.6,1.,.9,1.1,.6]),cuts={})
    _,one,_,same=transition(np.ones(5),shared,0.,.8,1)
    assert len(one['cuts'])==1 and same['events'][0]['x']==3.5
    assert same['distinct_maximum_suppressions']
    # Disconnected active regions each create a fixed break.
    shared['R']=np.array([.6,1.,.6,1.1,.6])
    _,two,_,separate=transition(np.ones(5),shared,0.,.8,1)
    assert len(two['cuts'])==2 and len(separate['events'])==2
    _,unchanged,_,again=transition(w,a,0.,.8,ident)
    assert len(unchanged['cuts'])==1 and not again['events']
    shared['R']=np.array([.6,1.,.9,1.,.6])
    try:
        transition(np.ones(5),shared,0.,.8,1)
    except RuntimeError as error:
        assert 'Equal distinct highest maxima' in str(error)
    else:
        raise AssertionError('Equal distinct highest peaks silently selected')
    return dict(success=True, checks=['birth at detected maximum by conservative cell split',
        'fixed gap remains fixed while maximum drifts in its active region',
        'distinct suppressed maximum explicitly flagged', 'local low Rac closes gap',
        'active-region merger does not close existing gaps',
        'two maxima in one connected high region create one break at taller peak',
        'two disconnected high regions each create one break',
        'unchanged state creates no duplicate birth',
        'equal distinct highest peaks report ambiguity'],
        event_material_impulse=result['material_impulse'])


if __name__ == '__main__':
    import json
    print(json.dumps(event_preflight()))
