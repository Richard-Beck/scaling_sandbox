"""Existing-face-only gap placement; native component selection and lifecycle.

Only placement and active-face birth eligibility change. Original numerical
faces never move, split, or disappear. Unresolved active face ties reject.
"""
from copy import deepcopy
import numpy as np
from native_gap_events import geometry, inventory, active_components
from ring_core import face_trace, FIELDS

class FaceLocatorError(RuntimeError):
    def __init__(self, message, locator):
        super().__init__(message)
        self.locator = locator

def locate_face(widths, rac, peak, threshold):
    edges, _ = geometry(widths)
    C = float(edges[-1]); x = float(peak['position'] % C)
    coordinate_tolerance = 64*np.finfo(float).eps*max(1.,C)
    rac_tolerance = 64*np.finfo(float).eps*max(1.,float(np.max(abs(rac))))
    distance = abs((edges[:-1]-x+C/2)%C-C/2)
    minimum = float(distance.min())
    faces = np.flatnonzero(distance <= minimum+coordinate_tolerance).tolist()
    if minimum <= coordinate_tolerance:
        faces = [int(np.argmin(distance))]
    evidence = dict(selected_peak_position=x, distance_tolerance=float(coordinate_tolerance),
        Rac_tolerance=float(rac_tolerance), nearest_distance=minimum, candidates=[],
        rule='nearest existing face; higher native face Rac only for geometric tie')
    for face in faces:
        value = float(face_trace(rac,widths,face))
        classification = 'active' if value > threshold+rac_tolerance else 'inactive' if value < threshold-rac_tolerance or value == threshold else 'threshold_ambiguous'
        evidence['candidates'].append(dict(face=int(face),x=float(edges[face]),distance=float(distance[face]),Rac_trace=value,classification=classification))
    if any(c['classification']=='threshold_ambiguous' for c in evidence['candidates']):
        raise FaceLocatorError('Candidate face at activity-threshold roundoff ambiguity',evidence)
    if not any(c['classification']=='active' for c in evidence['candidates']):
        return dict(evidence,chosen_face=None,chosen_x=None,outcome='deferred_inactive_face')
    candidates=sorted(evidence['candidates'],key=lambda c:-c['Rac_trace'])
    if len(candidates)>1 and abs(candidates[0]['Rac_trace']-candidates[1]['Rac_trace'])<=rac_tolerance:
        raise FaceLocatorError('Equal distinct active candidate face traces; no directional tie convention',evidence)
    chosen=candidates[0]
    return dict(evidence,chosen_face=chosen['face'],chosen_x=chosen['x'],outcome='active_existing_face')

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
    locator_records = []
    deferred_births = []
    for peak in sorted(births, key=lambda p: p['position']):
        location = locate_face(widths, state['R'], peak, threshold)
        row = dict(time=float(time), selected_peak=deepcopy(peak), **location)
        locator_records.append(row)
        if location['chosen_face'] is None:
            deferred_births.append(row)
            continue
        face = int(location['chosen_face'])
        x = float(location['chosen_x'])
        if face in state['cuts']:
            raise RuntimeError('Birth would duplicate an existing physical point')
        displacement = face_trace(state['u'], widths, face)
        current_trace = face_trace(state['R'], widths, face)
        state['cuts'][face] = dict(id=int(next_id), x=x, birth_time=float(time),
                                  Uleft=displacement, Uright=displacement,
                                  raw_left=displacement, raw_right=displacement)
        records.append(dict(event='open', time=float(time), id=int(next_id), x=x,
                            detected_peak_Rac=float(peak['maximum']),
                            selected_peak_position=float(peak['position']),
                            inserted_face_Rac=current_trace, initialized_U=displacement,
                            removed_material=0., cells_after=len(widths), locator=deepcopy(location)))
        next_id += 1
    mass_after = float(np.dot(widths, state['b'] + state['c']))
    if abs(mass_after - mass_before) > 1e-11 * max(1., mass_before):
        raise RuntimeError('Topology event changed MCA material')
    if len(widths) > 2000 or np.min(widths) < 1e-6 * (62. / 402):
        raise RuntimeError('Numerical mesh guard reached; pilot stops rather than merging nearby gaps')
    return widths, state, next_id, dict(events=records,
        distinct_maximum_suppressions=suppressed, material_impulse=mass_after-mass_before,
        prebirth_inventory=info, locator_records=locator_records, deferred_births=deferred_births)

