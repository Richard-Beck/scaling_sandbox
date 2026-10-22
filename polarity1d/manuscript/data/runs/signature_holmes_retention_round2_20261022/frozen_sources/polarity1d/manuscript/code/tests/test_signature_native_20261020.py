"""Regression for qualified basals with unresolved symmetry-breaking seed time."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
from signature_native_20261020 import select_preparation_endpoint


RULE = {"selection": {"minimum_contrast": .05, "minimum_strength": .05,
                       "relative_strength_tolerance": .01}}


def trial(duration, strength, *, verified=True, contrast=.8):
    return dict(duration=duration, verified=verified, final_orientation=-strength,
                cue_aligned_strength=strength, final_contrast=contrast,
                polarity_arrival=1. if verified else None)


def test_tiny_positive_seed_qualifies_endpoint_without_minimum_seed_claim():
    # As observed in Holmes: the zero trial remains homogeneous while every
    # tested positive seed reaches the same fully verified polarized endpoint.
    trials = [trial(0., 0., contrast=0.), trial(.0001953125, .326),
              trial(.000390625, .326), trial(.05, .326), trial(5120., .326)]
    result = select_preparation_endpoint(trials, RULE)
    assert result["selected"]["duration"] == .0001953125
    assert result["selected_relative_duration_bracket"] == 1.
    assert result["preparation_duration_unresolved"] is True
    assert result["upper_optimum_unresolved"] is False


def test_unverified_polarized_endpoints_cannot_qualify():
    with pytest.raises(RuntimeError, match="no independently converged"):
        select_preparation_endpoint([trial(0., 0., contrast=0.),
                                     trial(.1, .326, verified=False)], RULE)


def test_still_increasing_upper_strength_rejects_qualification():
    with pytest.raises(RuntimeError, match="upper observed-strength optimum"):
        select_preparation_endpoint([trial(0., 0., contrast=0.),
                                     trial(.1, .2), trial(1., .3)], RULE)


def test_unverified_upper_trial_is_not_hidden_by_an_earlier_stable_endpoint():
    with pytest.raises(RuntimeError, match="upper observed-strength optimum"):
        select_preparation_endpoint([trial(0., 0., contrast=0.), trial(.1, .326),
                                     trial(1., .326, verified=False)], RULE)
