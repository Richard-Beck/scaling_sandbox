"""Adversarial scientific controls for the signature duration calibration."""
import importlib.util
from pathlib import Path
import sys
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "mechanisms/signature_calibration_20261020.py"
SPEC = importlib.util.spec_from_file_location("signature_calibration", SOURCE)
CAL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CAL
SPEC.loader.exec_module(CAL)


class DurationCalibrationTests(unittest.TestCase):
    def config(self, **changes):
        settings = dict(horizon=32., scan_durations=(1., 2., 4., 8., 16.),
                        absolute_tolerance=.001, relative_tolerance=0.)
        settings.update(changes)
        return CAL.DurationSearch(**settings)

    def test_recovers_analytical_threshold_with_failed_field_tail(self):
        calls = []
        def evaluate(duration):
            calls.append(duration)
            return dict(final_d=.8 if duration >= 7.3 else -.9,
                        field_tail_pass=False)
        result = CAL.calibrate(self.config(), evaluate)
        self.assertTrue(result["train_ready"])
        self.assertLess(result["lower_failing_duration"], 7.3)
        self.assertGreaterEqual(result["upper_committing_duration"], 7.3)
        self.assertLessEqual(result["upper_committing_duration"] -
                             result["lower_failing_duration"], .001)
        self.assertEqual(len(calls), len(set(calls)))
        self.assertFalse(result["field_tail_used_to_define_fate"])

    def test_nonmonotone_dwell_response_blocks_automatic_train(self):
        result = CAL.calibrate(self.config(), lambda t: dict(
            final_d=.6 if 3 <= t <= 10 or t >= 30 else -.6,
            field_tail_pass=True))
        self.assertTrue(result["observed_nonmonotone_commitment"])
        self.assertFalse(result["train_ready"])
        self.assertLess(result["lower_failing_duration"], 3.)
        self.assertGreaterEqual(result["upper_committing_duration"], 3.)

    def test_censoring_has_no_fabricated_threshold(self):
        result = CAL.calibrate(self.config(), lambda t: dict(
            final_d=.49, field_tail_pass=True))
        self.assertEqual(result["status"], "no_sampled_commitment_through_finite_horizon")
        self.assertIsNone(result["upper_committing_duration"])
        self.assertFalse(result["train_ready"])

    def test_no_cue_replacement_invalidates_incumbent(self):
        result = CAL.calibrate(self.config(), lambda t: dict(
            final_d=.7, field_tail_pass=True))
        self.assertEqual(result["status"], "zero_duration_already_replaces_invalid_incumbent_control")
        self.assertFalse(result["train_ready"])

    def test_debelly_half_unit_lattice_is_preserved(self):
        seen = []
        def evaluate(duration):
            seen.append(duration)
            return dict(final_d=.5 if duration >= 7.25 else -.5,
                        field_tail_pass=True, actual_duration=duration)
        result = CAL.calibrate(self.config(native_step=.5), evaluate)
        self.assertEqual((result["lower_failing_duration"],
                          result["upper_committing_duration"]), (7., 7.5))
        self.assertTrue(all(t * 2 == round(t * 2) for t in seen))
        self.assertTrue(result["train_ready"])

    def test_off_lattice_input_is_rejected_before_any_solver_call(self):
        with self.assertRaisesRegex(ValueError, "native lattice"):
            CAL.calibrate(self.config(native_step=.5, scan_durations=(.1,)),
                          lambda t: self.fail("No native call should occur"))

    def test_missing_readout_is_not_replacement_failure(self):
        with self.assertRaises(ValueError):
            CAL.calibrate(self.config(), lambda t: dict(
                final_d=float("nan"), field_tail_pass=True))

    def test_iteration_limit_cannot_enable_train(self):
        result = CAL.calibrate(self.config(max_bisections=1), lambda t: dict(
            final_d=.6 if t >= 7.3 else -.7, field_tail_pass=True))
        self.assertFalse(result["bracket_resolved"])
        self.assertFalse(result["train_ready"])

    def test_explicit_unpolarized_state_is_recorded_separately(self):
        result = CAL.calibrate(self.config(), lambda t: dict(
            final_d=None, field_tail_pass=True, readout_eligible=False,
            ineligibility_reason="positive activity excess below frozen basal fraction"))
        self.assertFalse(result["train_ready"])
        self.assertTrue(all(r["final_d"] is None and not r["replacement"]
                            for r in result["evaluations"]))

    def test_unexplained_missing_dominance_aborts(self):
        with self.assertRaises(TypeError):
            CAL.calibrate(self.config(), lambda t: dict(
                final_d=None, field_tail_pass=True))


if __name__ == "__main__":
    unittest.main()
