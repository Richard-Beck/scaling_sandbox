"""A direction-quality gap must not become an invented heading crossing."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure3_heading_20261020 import target_bracket


class Figure3HeadingTests(unittest.TestCase):
    def test_valid_heading_crossing_is_saved_native_bracket(self):
        t = np.arange(5.)
        p = .5*np.exp(1j*np.array([0., .2, .4, .6, .8]))
        result = target_bracket(t, p, .5, 1., .05)
        self.assertEqual(result["status"], "observed_saved_valid_segment_crossing")
        self.assertEqual((result["lower"], result["upper"]), (2., 3.))
        self.assertTrue(result["first_global_crossing_proven"])

    def test_cancelled_vector_is_not_unwrapped_into_a_turn_time(self):
        t = np.arange(4.)
        p = np.array([.5, .1, -.1, -.5], dtype=complex)
        for floor in (1e-10, .05):
            result = target_bracket(t, p, .9*np.pi, 1., floor)
            self.assertEqual(result["status"], "target_attained_when_heading_returns_after_gap")
            self.assertIsNone(result["lower"])
            self.assertIsNone(result["upper"])

    def test_finite_heading_cap_is_a_lower_bound_only_without_quality_gaps(self):
        result = target_bracket(np.arange(3.), .5*np.exp(1j*np.array([0., .2, .4])), .8, 1., .05)
        self.assertEqual(result["status"], "right_censored_saved_heading_at_completed_horizon")
        self.assertEqual(result["lower"], 2.)


if __name__ == "__main__":
    unittest.main()
