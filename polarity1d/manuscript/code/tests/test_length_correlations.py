"""Statistical edge cases for the scalar-success correlation policy."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import three_test_figures as t


def summaries(monkeypatch, durations, flags=None, extra=()):
    monkeypatch.setattr(t.f, 'MODELS', ['wave_pinning'])
    monkeypatch.setattr(t, 'TESTS', ['fixed_physical'])
    monkeypatch.setattr(t.f, 'GATES', [('orientation', '', '', '')])
    flags = flags if flags is not None else [False] * len(durations)
    rows = [dict(model='wave_pinning', condition='fixed_physical', definition='orientation',
                 length=i + 1, endpoint=value, status='success', selected_field_tail_pass=flag)
            for i, (value, flag) in enumerate(zip(durations, flags))]
    for i, status in enumerate(extra):
        rows.append(dict(rows[0], length=100 + i, endpoint=None, status=status,
                         selected_field_tail_pass=None))
    result = t.correlation_summary(rows)
    t.correlation_check(rows, result)
    return result[0]


def test_scalar_successes_include_failed_and_unknown_flags_but_not_censoring(monkeypatch):
    s = summaries(monkeypatch, [1, 2, 3, 4], [True, False, None, False],
                  ['censored', 'ineligible', 'preparation_invalid'])
    assert s['spearman_rho'] == pytest.approx(1)
    assert s['n_correlation'] == 4 and s['n_total'] == 7
    assert (s['n_field_tail_passed'], s['n_field_tail_failed'], s['n_field_tail_unknown']) == (1, 2, 1)
    assert (s['n_censored'], s['n_ineligible'], s['n_preparation_invalid']) == (1, 1, 1)


@pytest.mark.parametrize('durations,expected', [([4, 3, 2, 1], -1),
                                               ([2, 4, 1, 3], 0),
                                               ([1, 2, 2, 4], np.sqrt(.9))])
def test_negative_zero_and_tied_correlations_remain_available(monkeypatch, durations, expected):
    s = summaries(monkeypatch, durations)
    assert s['correlation_status'] == 'available'
    assert s['spearman_rho'] == pytest.approx(expected)


@pytest.mark.parametrize('durations,status', [([1, 2, 3], 'insufficient_successes'),
                                            ([2, 2, 2, 2], 'constant_duration')])
def test_insufficient_or_constant_data_are_never_colored_as_zero(monkeypatch, durations, status):
    s = summaries(monkeypatch, durations)
    assert s['correlation_status'] == status and s['spearman_rho'] is None
