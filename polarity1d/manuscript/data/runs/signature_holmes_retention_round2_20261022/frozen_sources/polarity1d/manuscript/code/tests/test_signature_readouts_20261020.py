"""Analytic fractional-integration, moving-geometry and fate-readout checks."""
import json
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
from signature_readouts_20261020 import (
    analyze_trace, describe_basal, first_crossing, integrate_profile, measure_profile,
)


class SignatureReadoutTests(unittest.TestCase):
    def test_cell_integral_fractional_support_with_variable_jacobian(self):
        x, edges, active = [.5, 2, 3.5], [0, 1, 3, 4], [1, 4, 2]
        self.assertAlmostEqual(integrate_profile(
            x, active, [(.5, 3.5)], sample_kind="cells", reference_edges=edges), 9.5)
        self.assertAlmostEqual(integrate_profile(
            x, active, [(.5, 3.5)], sample_kind="cells", reference_edges=edges,
            physical_edges=[10, 12, 18, 22]), 29)

    def test_nodal_piecewise_linear_integral_fractional_support(self):
        self.assertAlmostEqual(integrate_profile([0, 1, 3], [0, 2, 0]), 3)
        self.assertAlmostEqual(integrate_profile(
            [0, 1, 3], [0, 2, 0], [(.5, 2)]), 2.25)
        self.assertAlmostEqual(integrate_profile(
            [0, 1, 3], [0, 2, 0], [(.5, 2)], physical_x=[-2, 0, 6]), 6)

    def test_interval_union_does_not_double_count_or_exceed_domain(self):
        self.assertAlmostEqual(integrate_profile(
            [0, 1, 2], [1, 1, 1], [(-4, .8), (.5, 1.5), (1.2, 8)]), 2)
        self.assertEqual(integrate_profile([0, 1], [1, 1], []), 0)

    def test_fractional_midcell_partition_is_additive(self):
        result = measure_profile(
            [.35, 1.05, 1.7], [0, 2, 6], sample_kind="cells",
            reference_edges=[0, .7, 1.4, 2], incumbent_intervals=[(0, .8)],
            challenger_own_intervals=[(.5, 2)], challenger_fixed_intervals=[(.5, 2)])
        self.assertAlmostEqual(result["E_left"], .6)
        self.assertAlmostEqual(result["E_right"], 4.4)
        self.assertAlmostEqual(result["E_total"], 5)
        self.assertAlmostEqual(result["D"], .76)
        self.assertAlmostEqual(result["partition_sum_error"], 0)
        self.assertTrue(result["roi_own_overlap"])
        self.assertFalse(result["roi_own_eligible"])
        self.assertTrue(result["replacement"])

    def test_exact_75_percent_share_and_zone_independence(self):
        common = dict(reference_x=[.25, .75, 1.25, 1.75], active=[2, 0, 6, 0],
                      sample_kind="cells", reference_edges=[0, .5, 1, 1.5, 2],
                      incumbent_intervals=[(0, .5)])
        narrow = measure_profile(**common, challenger_own_intervals=[(1, 1.5)],
                                 challenger_fixed_intervals=[(.25, 2)])
        broad = measure_profile(**common, challenger_own_intervals=[(.25, 2)],
                                challenger_fixed_intervals=[(.25, 2)])
        self.assertEqual(narrow["D"], .5)
        self.assertEqual(narrow["share_half"], .75)
        self.assertEqual(narrow["D"], broad["D"])
        self.assertFalse(narrow["roi_own_overlap"])
        self.assertTrue(broad["roi_own_overlap"])
        self.assertTrue(narrow["roi_fixed_overlap"])

    def test_homogeneous_and_tiny_excess_are_ineligible(self):
        common = dict(reference_x=[0, .5, 1], incumbent_intervals=[(0, .5)],
                      challenger_own_intervals=[(.5, 1)],
                      challenger_fixed_intervals=[(.5, 1)], basal_E_total=1,
                      basal_contrast=1)
        for profile in ([2, 2, 2], [0, 0, 1e-8]):
            with self.subTest(profile=profile):
                result = measure_profile(active=profile, **common)
                self.assertIsNone(result["D"])
                self.assertFalse(result["readout_eligible"])
                self.assertFalse(result["replacement"])
                self.assertIsNotNone(result["ineligibility_reason"])

    def test_low_contrast_ratio_does_not_produce_replacement(self):
        result = measure_profile(
            [0, .5, 1], [1, 1, 1.01], incumbent_intervals=[(0, .5)],
            challenger_own_intervals=[(.5, 1)], challenger_fixed_intervals=[(.5, 1)])
        self.assertEqual(result["D_half_raw"], 1)
        self.assertIsNone(result["D"])
        self.assertFalse(result["replacement"])
        self.assertFalse(result["contrast_eligible"])

    def test_exact_frozen_basal_amount_boundary_is_eligible(self):
        result = measure_profile(
            [0, .5, 1], [0, 0, .04], incumbent_intervals=[(0, .5)],
            challenger_own_intervals=[(.5, 1)], challenger_fixed_intervals=[(.5, 1)],
            basal_E_total=1, basal_contrast=1)
        self.assertAlmostEqual(result["E_total"], .01)
        self.assertTrue(result["readout_eligible"])

    def test_spring_translation_and_uniform_stretch_preserve_share(self):
        args = dict(reference_x=[0, .5, 1], active=[0, 0, 2],
                    incumbent_intervals=[(0, .3)], challenger_own_intervals=[(.6, 1)],
                    challenger_fixed_intervals=[(.5, 1)])
        reference = measure_profile(**args)
        moving = measure_profile(**args, physical_x=[17, 19, 21])
        self.assertEqual(reference["D"], moving["D"])
        self.assertAlmostEqual(moving["E_total"], 4 * reference["E_total"])
        self.assertAlmostEqual(moving["physical_length"], 4)

    def test_basal_halfheight_support_is_fractional_and_frozen(self):
        basal = describe_basal([0, 1, 2, 3], [5, 3, 1, 1])
        self.assertEqual(basal["current_min_floor"], 1)
        self.assertEqual(basal["incumbent_intervals"], [[0., 1.]])
        self.assertAlmostEqual(basal["incumbent_extent_physical"], 1)
        self.assertTrue(basal["peak_in_expected_half"])
        self.assertFalse(basal["multiple_halfheight_regions"])

    def test_crossing_preserves_observed_bracket_and_missing_coverage(self):
        crossing = first_crossing([0, 2, 5], [-1, 0, 1], .5)
        self.assertEqual(crossing["lower_time"], 2)
        self.assertEqual(crossing["upper_time"], 5)
        self.assertEqual(crossing["estimated_time"], 3.5)
        missing = first_crossing([0, 2, 5], [-1, None, 1], .5)
        self.assertIsNone(missing["lower_time"])
        self.assertIsNone(missing["estimated_time"])
        self.assertTrue(missing["left_censored"])
        censored = first_crossing([0, 2, 5], [-1, 0, .2], .5)
        self.assertFalse(censored["observed"])
        self.assertEqual(censored["right_censored_at"], 5)

    def test_trace_freezes_basal_scales_and_has_strict_json_output(self):
        trace = dict(time=[0, 1, 2], x=[0, .5, 1], sample_kind="nodes",
                     active=[[2, 0, 0], [0, 0, 2], [0, 0, 1e-8]],
                     physical_x=[[0, .5, 1], [3, 4, 5], [2, 3, 4]])
        case = dict(challenger_own_intervals=[(.7, 1)], challenger_fixed_intervals=[(.5, 1)])
        result = analyze_trace(trace, case)
        self.assertEqual(result["series"]["D"], [-1., 1., None])
        self.assertIsNone(result["final_d"])
        self.assertFalse(result["final_replacement"])
        self.assertFalse(result["readout_eligible"])
        self.assertIn("amount", result["ineligibility_reason"])
        self.assertEqual(result["series"]["basal_E_total"], [.5, .5, .5])
        json.dumps(result, allow_nan=False)

    def test_physical_width_tracks_stretch_and_exact_inverse_map(self):
        trace = dict(time=[0, 1], xref=[0, .5, 1], sample_kind="nodes",
                     active=[[0, 0, 2], [0, 0, 2]],
                     physical_x=[[0, .5, 1], [10, 11, 14]])
        case = dict(challenger_own_intervals=[(.5, 1)],
                    challenger_fixed_intervals=[(0, 1)],
                    challenger_own_width_physical=1,
                    challenger_fixed_width_physical=4)
        result = analyze_trace(trace, case)
        self.assertEqual(result["challenger_own_intervals_by_frame"][0], [[0., 1.]])
        self.assertAlmostEqual(result["challenger_own_intervals_by_frame"][1][0][0], 5 / 6)
        self.assertAlmostEqual(result["series"]["E_challenger_own"][1], 5 / 3)
        self.assertEqual(result["series"]["D"], [1., 1.])

    def test_left_challenger_and_physical_support(self):
        trace = dict(time=[0], x=[0, .5, 1], active=[[2, 0, 0]],
                     sample_kind="nodes", physical_x=[3, 4, 5])
        case = dict(challenger_side="left", challenger_own_intervals=[(0, .5)],
                    challenger_fixed_intervals=[(0, .5)],
                    challenger_own_width_physical=.5)
        result = analyze_trace(trace, case)
        self.assertEqual(result["challenger_own_intervals_by_frame"], [[[0., .25]]])
        self.assertEqual(result["final_d"], 1)
        self.assertTrue(result["final_replacement"])

    def test_widest_zone_is_frozen_material_while_own_follows_geometry(self):
        trace = dict(time=[0, 1], x=[0, .5, 1], active=[[0, 0, 2], [0, 0, 2]],
                     physical_x=[[0, .5, 1], [10, 11, 14]])
        case = dict(challenger_own_intervals=[(.5, 1)],
                    challenger_fixed_intervals=[(.5, 1)],
                    challenger_own_width_physical=.5,
                    challenger_fixed_width_physical=.5)
        result = analyze_trace(trace, case)
        self.assertEqual(result["challenger_fixed_intervals_by_frame"], [[[.5, 1.]], [[.5, 1.]]])
        self.assertAlmostEqual(result["challenger_own_intervals_by_frame"][1][0][0], 11 / 12)

    def test_explicit_frame_specific_regional_intervals(self):
        trace = dict(time=[0, 1], x=[0, .5, 1], active=[[0, 0, 2], [0, 0, 2]])
        case = dict(challenger_own_intervals=[[[.5, 1]], [[.75, 1]]],
                    challenger_fixed_intervals=[[[0, 1]], [[.5, 1]]])
        result = analyze_trace(trace, case)
        self.assertEqual(result["challenger_own_intervals_by_frame"], [[[.5, 1.]], [[.75, 1.]]])
        self.assertAlmostEqual(result["series"]["E_challenger_own"][0], .5)
        self.assertAlmostEqual(result["series"]["E_challenger_own"][1], .375)
        self.assertEqual(result["series"]["D"], [1., 1.])

    def test_invalid_or_reversed_geometry_fails_explicitly(self):
        for x, y in (([0, 0, 1], [1, 2, 3]), ([0, 1], [1, np.nan])):
            with self.subTest(x=x, y=y), self.assertRaises(ValueError):
                integrate_profile(x, y)
        with self.assertRaises(ValueError):
            integrate_profile([.5], [1], sample_kind="cells")
        with self.assertRaises(ValueError):
            integrate_profile([0, 1], [1, 1], physical_x=[2, 1])


if __name__ == "__main__":
    unittest.main()
