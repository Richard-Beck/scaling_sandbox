"""Area and discontinuity controls independent of the native solver engines."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / "mechanisms/signature_inputs_20261020.py"
SPEC = importlib.util.spec_from_file_location("signature_inputs", SOURCE)
INPUTS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = INPUTS
SPEC.loader.exec_module(INPUTS)


class CueControls(unittest.TestCase):
    def test_commanded_area_is_identical_on_all_widths_lengths_and_meshes(self):
        for cells in (201, 401):
            for length in (19., 39.):
                x = (np.arange(cells) + .5) * length / cells
                for width in (.95, 1.9, 3.8, 7.6, 15.2):
                    cue = INPUTS.AreaCue(.38, width)
                    diagnostic = cue.profile_diagnostics(x, length)
                    self.assertAlmostEqual(diagnostic["sampled_area"], .38, places=13)
                    self.assertGreater(diagnostic["sampled_peak"], 0)
                    self.assertFalse(diagnostic["support_truncated"])

    def test_spring_area_survives_moving_boundaries_and_reports_truncation(self):
        for length in (.35, .45, .8):
            x = np.linspace(0, length, 201)
            cue = INPUTS.AreaCue(.4, .4, grid_kind="node")
            diagnostic = cue.profile_diagnostics(x, length)
            self.assertAlmostEqual(diagnostic["sampled_area"], .4, places=13)
            self.assertEqual(diagnostic["support_truncated"], length < .4)

    def test_left_and_right_shapes_have_identical_area(self):
        x = (np.arange(201) + .5) * 39 / 201
        right = INPUTS.AreaCue(.38, .95)(0, x, 39)
        left = INPUTS.AreaCue(.38, .95, side="left")(0, x, 39)
        np.testing.assert_allclose(left, right[::-1], rtol=1e-12, atol=1e-12)

    def test_grid_convention_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "declared native grid"):
            INPUTS.AreaCue(.4, .4)(0, np.linspace(0, 1, 201), 1)

    def test_unresolved_cue_footprint_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not resolve"):
            INPUTS.AreaCue(.4, .001)(0, (np.arange(5) + .5) / 5, 1)

    def test_train_exposure_equal_across_counts_on_debelly_lattice(self):
        schedules = [INPUTS.make_train(163.5, 1.5, n, .1, native_step=.5)
                     for n in (2, 4, 8)]
        self.assertEqual(len({s.actual_total_on for s in schedules}), 1)
        for schedule in schedules:
            self.assertAlmostEqual(sum(stop-start for start, stop in schedule.on_intervals),
                                   schedule.actual_total_on)
            self.assertTrue(all(abs(t * 2 - round(t * 2)) < 1e-9
                                for t in schedule.event_times))
            self.assertGreaterEqual(schedule.actual_total_on, schedule.requested_total_on)
            self.assertLess(schedule.actual_total_on-schedule.requested_total_on, 4.)
            self.assertEqual(len(schedule.on_intervals), schedule.pulse_count)

    def test_zero_gap_train_is_exactly_continuous_exposure(self):
        schedule = INPUTS.make_train(7.3, 1.5, 8, 0)
        self.assertAlmostEqual(schedule.end_time, 7.3 * 1.5)
        self.assertTrue(all(kind == "on" for kind, _, _ in schedule.segments))
        self.assertEqual(len(schedule.segments), 8)

    def test_gap_count_excludes_final_recovery(self):
        schedule = INPUTS.make_train(8., 1., 4, 4.)
        gaps = [segment for segment in schedule.segments if segment[0] == "gap"]
        self.assertEqual(len(gaps), 3)
        self.assertAlmostEqual(schedule.end_time, 8. + 3 * 8.)


if __name__ == "__main__":
    unittest.main()
