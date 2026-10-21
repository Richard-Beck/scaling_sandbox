"""Geometry remaps preserve domain gaps and censored/topology outcomes."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure3_native_20261020 as native
import signature_figure3_synthesis_20261020 as synthesis


class Figure3SynthesisTests(unittest.TestCase):
    def row(self, slot, gap, value):
        return dict(slot=slot, gap=gap, handover_time=value, classification="external",
            control=None, handover_status="observed" if value is not None else "right_censored_at_completed_native_horizon",
            handover_lower_bound=20. if value is None else None)

    def test_no_interpolation_through_censoring_and_no_antipode_extrapolation(self):
        rows = [self.row(0, 0., 10.), self.row(1, 1., None), self.row(2, 2., 30.)]
        self.assertFalse(synthesis.curve_value(rows, .5)["available"])
        self.assertEqual(synthesis.curve_value(rows, 1.)["lower_bound"], 20.)
        self.assertEqual(synthesis.curve_value(rows, 2.1)["status"], "outside_observed_external_gap_support")
        self.assertEqual(synthesis.curve_value([self.row(0, 0., 10.), self.row(1, 1., 20.)], .5)["value"], 15.)

    def test_measured_width_enters_geometry_and_size_decomposition_adds_exactly(self):
        rows = []
        for model in native.MODELS:
            for size in native.SIZES[model]:
                C = size*(1.02 if model == "spring" else 1.)
                width, support = .4*C, .2 if model == "spring" else 3.
                maximum = C/2-(width+support)/2
                for gap in np.linspace(0., maximum, 9):
                    row = self.row(len(rows), float(gap), 100+3*gap+2*size)
                    row.update(model=model, size_parameter=size, basal_circumference=C,
                        initial_10_percent_width=width, cue_support=support)
                    rows.append(row)
        remapped = synthesis.fixed_angle_remaps(rows)
        for row in remapped:
            if row["theta_degrees"] != 120.:
                continue
            self.assertTrue(row["all_components_available"])
            baseline_size = native.SIZES[row["model"]][0]
            self.assertAlmostEqual(row["matched_gap_size_time_shift"], 2*(row["size_parameter"]-baseline_size), places=10)
            self.assertAlmostEqual(row["geometry_time_shift"]+row["matched_gap_size_time_shift"], row["measured_time"]-row["small_same_angle_time"], places=10)


if __name__ == "__main__":
    unittest.main()
