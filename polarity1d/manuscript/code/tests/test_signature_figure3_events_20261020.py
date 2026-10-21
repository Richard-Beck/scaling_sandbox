"""Continuous ratio roots versus discontinuous topology measurement boundaries."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure3_events_20261020 import RegionSeedEvent


class RegionEventTests(unittest.TestCase):
    def detector(self, readout):
        placement = dict(initial_10_percent_contour=dict(origin=0.), basal_circumference=1., cue_center_fraction=.4)
        detector = RegionSeedEvent(None, placement, 1., 1.)
        detector.readout = lambda state: readout(float(state[0]))
        return detector

    def two(self, challenger, incumbent=1.):
        return dict(two_region_available=True, state_class="two_resolved_spatially_associated_regions",
                    challenger_excess=challenger, incumbent_excess=incumbent)

    def test_valid_two_region_ratio_has_exact_continuous_root(self):
        readout = lambda t: self.two(.2+1.6*t)
        detector = self.detector(readout)
        result = detector.locate(0., readout(0.), 1., readout(1.), lambda t: np.array([t]))
        self.assertEqual(result["kind"], "continuous_two_region_ratio")
        self.assertTrue(result["continuous_root"])
        self.assertAlmostEqual(result["time"], .5, places=9)
        self.assertAlmostEqual(result["root_value"], 0., places=9)

    def test_challenger_born_above_target_is_not_a_ratio_equality_root(self):
        def readout(t):
            return self.two(2.) if t >= .37 else dict(two_region_available=False, state_class="single_region_incumbent_associated", incumbent_excess=1.)
        detector = self.detector(readout)
        result = detector.locate(0., readout(0.), 1., readout(1.), lambda t: np.array([t]))
        self.assertEqual(result["kind"], "target_attainment_at_resolved_region_birth")
        self.assertFalse(result["continuous_root"])
        self.assertAlmostEqual(result["time"], .37, places=8)
        self.assertLessEqual(result["measurement_boundary_bracket"][1]-result["measurement_boundary_bracket"][0], 1e-9)

    def test_single_region_relocation_is_not_invented_as_handover(self):
        def readout(t):
            return dict(two_region_available=False,
                        state_class="single_region_incumbent_associated" if t < .63 else "single_region_cue_associated", incumbent_excess=1. if t < .63 else 0.)
        detector = self.detector(readout)
        result = detector.locate(0., readout(0.), 1., readout(1.), lambda t: np.array([t]))
        self.assertEqual(result["kind"], "takeover_during_cue_single_region")
        self.assertFalse(result["continuous_root"])
        self.assertAlmostEqual(result["time"], .63, places=8)

    def test_earlier_incumbent_collapse_wins_over_later_ratio(self):
        # Ratio crossing occurs at t=0.65, after the incumbent90% loss at0.6.
        readout = lambda t: self.two(.01+.1*t, .4-.5*t)
        detector = self.detector(readout)
        result = detector.locate(0., readout(0.), .7, readout(.7), lambda t: np.array([t]))
        self.assertEqual(result["kind"], "continuous_incumbent_90_percent_excess_loss")
        self.assertTrue(result["continuous_root"])
        self.assertAlmostEqual(result["time"], .6, places=8)


if __name__ == "__main__":
    unittest.main()
