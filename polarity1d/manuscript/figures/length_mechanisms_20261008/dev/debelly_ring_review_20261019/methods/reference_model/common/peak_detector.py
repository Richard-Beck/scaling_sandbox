#!/usr/bin/env python3
"""Snapshot circular maxima audit, not gap tracking or an evolved model."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
RUN = ROOT / 'polarity1d/manuscript/data/runs/debelly_shared_chemical_arcs_20261017'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def periodic_maxima(values, circumference, threshold=None, positions=None):
    """Maxima of a sampled circular profile, with numerical plateau grouping.

    Piecewise-linear reconstruction through the cell-center samples is used.
    A plateau has one midpoint. No minimum spacing/prominence filter is applied.
    This function deliberately assigns neither site identities nor gap motion.
    """
    values = np.asarray(values, dtype=float)
    size = len(values)
    if size < 3 or not np.all(np.isfinite(values)):
        raise ValueError('Need at least three finite periodic samples')
    if positions is None:
        positions = (np.arange(size) + .5) * circumference / size
    positions = np.asarray(positions, dtype=float)
    if (len(positions) != size or np.any(np.diff(positions) <= 0)
            or positions[0] < 0 or positions[-1] >= circumference):
        raise ValueError('Need increasing cell-center coordinates on [0,C)')
    tolerance = 64 * np.finfo(float).eps * max(1., float(np.max(abs(values))))
    if np.ptp(values) <= tolerance:
        return dict(status='uniform_no_isolated_maximum', tolerance=tolerance, peaks=[])
    boundaries = [i for i in range(size)
                  if abs(values[i] - values[(i - 1) % size]) > tolerance]
    if not boundaries:
        return dict(status='ambiguous_transitive_roundoff_plateau',
                    tolerance=tolerance, peaks=[])
    start = boundaries[0]
    ordered = [(start + j) % size for j in range(size)]
    groups = [[ordered[0]]]
    for index in ordered[1:]:
        if abs(values[index] - values[groups[-1][-1]]) <= tolerance:
            groups[-1].append(index)
        else:
            groups.append([index])
    peaks = []
    ambiguous = False
    for group in groups:
        lo, hi = float(np.min(values[group])), float(np.max(values[group]))
        if hi - lo > tolerance:
            ambiguous = True
            continue
        before = float(values[(group[0] - 1) % size])
        after = float(values[(group[-1] + 1) % size])
        if lo <= max(before, after) + tolerance:
            continue
        first = positions[group[0]]
        span = (positions[group[-1]] - first) % circumference
        position = (first + .5 * span) % circumference
        eligible = None
        if threshold is not None:
            if lo > threshold + tolerance:
                eligible = True
            elif hi < threshold - tolerance:
                eligible = False
            elif lo == hi == threshold:
                eligible = False  # The stipulated strict inequality.
            else:
                eligible = 'roundoff_ambiguous'
        peaks.append(dict(sample_indices=group, position=float(position),
                          minimum=float(lo), maximum=float(hi), above_threshold=eligible))
    return dict(status='roundoff_group_ambiguous' if ambiguous else 'resolved',
                tolerance=tolerance, peaks=sorted(peaks, key=lambda x: x['position']))


def main():
    fixtures = []
    def check(name, values, expected, threshold=.9):
        result = periodic_maxima(values, len(values), threshold)
        actual = len(result['peaks'])
        assert actual == expected, (name, result)
        fixtures.append(dict(name=name, result=result))
        return result
    check('uniform_below', [.5] * 8, 0)
    check('uniform_above_has_no_isolated_maximum', [1.2] * 8, 0)
    check('strict_single_maximum', [.4, .5, .7, 1.2, .8, .6, .5, .3], 1)
    wrapped = check('wrapped_plateau_is_one_maximum', [1.2, 1.2, .7, .5, .4, .5, .7, 1.2], 1)
    assert wrapped['peaks'][0]['sample_indices'] == [7, 0, 1]
    close = check('nearby_peaks_with_resolved_valley_remain_distinct',
                  [.4, 1.2, 1.1, 1.3, .4, .5, .4, .3], 3)
    assert [x['sample_indices'] for x in close['peaks']][:2] == [[1], [3]]
    equal = check('equality_does_not_open', [.3, .4, .9, .4, .3], 1)
    assert equal['peaks'][0]['above_threshold'] is False
    near = check('roundoff_plateau', [.3, 1.2, 1.2 + 2e-15, .3, .2], 1)
    assert near['peaks'][0]['sample_indices'] == [1, 2]
    profile = np.asarray([.4, 1.2, 1.1, 1.3, .4, .5, .4, .3])
    base = periodic_maxima(profile, 8.)
    rotation_errors = []
    for offset in range(len(profile)):
        rotated = periodic_maxima(np.roll(profile, offset), 8.)
        expected = sorted((p['position'] + offset) % 8. for p in base['peaks'])
        observed = [p['position'] for p in rotated['peaks']]
        rotation_errors.append(float(np.max(abs(np.asarray(expected) - observed))))
    assert max(rotation_errors) < 1e-12

    measured = []
    source_hashes = {}
    path = RUN / 'initial_full_loop.npz'
    source_hashes[str(path.relative_to(ROOT))] = sha(path)
    with np.load(path) as d:
        r, position = d['rac'], d['s']
        indices = np.flatnonzero(r >= r.max() - 64 * np.finfo(float).eps)
        measured.append(dict(state='qualified_initial', maximum_Rac=float(r.max()),
                             maximizing_cell_positions=[float(position[i]) for i in indices],
                             detector=periodic_maxima(r, 62.)))
    for slot in range(6):
        path = RUN / f'workers/{slot:02d}.npz'
        source_hashes[str(path.relative_to(ROOT))] = sha(path)
        with np.load(path) as d:
            a, b = d['arc0_rac'][-1], d['arc1_rac'][-1]
            r = np.r_[a, b[::-1]]
            cut = len(a)
            dx = 62. / len(r)
            indices = np.flatnonzero(r >= r.max() - 64 * np.finfo(float).eps)
            measured.append(dict(state=f'worker_{slot:02d}',
                                 condition='cue' if slot % 2 == 0 else 'gap_only',
                                 shorter_arc=cut * dx, time=float(d['time'][-1]),
                                 maximum_Rac=float(r.max()),
                                 maximizing_cell_positions=[float((i + .5) * dx) for i in indices],
                                 old_face_Rac=float(.5 * (r[0] + r[-1])),
                                 challenger_face_Rac=float(.5 * (r[cut - 1] + r[cut])),
                                 detector=periodic_maxima(r, 62.)))
    operator = json.loads((HERE / 'preflight.json').read_text())
    native_midpoint = operator['threshold']
    for row in measured:
        row['native_midpoint_minus_global_maximum'] = native_midpoint - row['maximum_Rac']
    result = dict(success=True, scope='Snapshot maxima only; no peak identity/tracking or dynamic cuts',
                  numerical_tolerance_rule='64 machine eps * max(1,max(abs(Rac)))',
                  native_midpoint_candidate=native_midpoint,
                  threshold_selected=False,
                  fixtures=fixtures, rotation_errors=rotation_errors,
                  archive_measurements=measured, source_sha256=source_hashes,
                  script_sha256=sha(__file__))
    (HERE / 'peak_preflight.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(success=True, fixture_count=len(fixtures),
                          maximum_rotation_error=max(rotation_errors),
                          archive_states_measured=len(measured),
                          native_midpoint_candidate=native_midpoint)))


if __name__ == '__main__':
    main()
