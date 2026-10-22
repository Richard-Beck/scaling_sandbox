"""Analytic checks for new S4 inputs and exact sampled-state protocol design."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
from signature_support_20261020 import (
    CenteredAreaCue, DWELL_FACTORS, FOCAL_MODELS, LENGTHS, RelativePeakCue,
    SHARES, actual_dwell, build_plan, position_centers, share_event,
)


class SignatureSupportTests(unittest.TestCase):
    def test_position_cue_support_is_in_right_half_and_area_is_fixed(self):
        for length in LENGTHS:
            for cells in (201, 401):
                x = (np.arange(cells) + .5) * length / cells
                for center in position_centers(length):
                    cue = CenteredAreaCue(.38, 3.8, center)
                    profile, weights, peak = cue.profile(x, length)
                    self.assertAlmostEqual(float(profile @ weights), .38, places=12)
                    self.assertGreater(peak, 0)
                    self.assertGreaterEqual(center - 1.9, length / 2 - 1e-12)
                    self.assertLessEqual(center + 1.9, length + 1e-12)
                    self.assertFalse(cue.profile_diagnostics(x, length)["support_truncated"])

    def test_same_peak_relative_cue_intentionally_changes_area_with_length(self):
        cue = RelativePeakCue(.2)
        integrals = []
        for length in (19, 39):
            x = (np.arange(201) + .5) * length / 201
            diagnostic = cue.profile_diagnostics(x, length)
            self.assertEqual(diagnostic["reconstructed_boundary_peak"], .2)
            integrals.append(diagnostic["sampled_area"])
        self.assertAlmostEqual(integrals[1] / integrals[0], 39 / 19, places=12)

    def test_additional_dwell_preserves_native_debelly_lattice(self):
        self.assertEqual(actual_dwell("debelly", 0), 0)
        self.assertEqual(actual_dwell("debelly", .625), 1)
        self.assertEqual(actual_dwell("debelly", 2), 2)
        self.assertEqual(actual_dwell("wave_pinning", .625), .625)
        with self.assertRaises(ValueError):
            actual_dwell("debelly", -1)

    def test_share_target_selects_exact_upper_sample_and_preserves_nonreach(self):
        measured = dict(series=dict(time=[0., 1., 2.], share_half=[0., None, .8],
                                    D_half_eligible=[True, False, True]))
        event = share_event(measured, .75)
        self.assertEqual(event["upper_sample_index"], 2)
        self.assertEqual(event["actual_upper_sample_share"], .8)
        self.assertFalse(event["native_state_interpolated"])
        self.assertTrue(event["lower_D_unavailable"])
        absent = share_event(measured, .9)
        self.assertEqual(absent["status"], "right_censored")
        self.assertNotIn("upper_sample_index", absent)

    def test_plan_counts_exact_basals_and_reuses_failed_field_tail_without_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calibration, width = root / "calibration", root / "width"
            (calibration / "results").mkdir(parents=True)
            width.mkdir()
            basals, calibrations, width_cases = [], [], []
            for model in FOCAL_MODELS:
                for length in LENGTHS:
                    identifier = f"{model}:L{length}"
                    checkpoint = calibration / (identifier.replace(":", "_") + ".npz")
                    checkpoint.write_bytes(identifier.encode())
                    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
                    basals.append(dict(basal_id=identifier, checkpoint=str(checkpoint),
                        checkpoint_sha256=digest, qualified=True, fixed_scales={"active": 1}))
                    calibrations.append(dict(calibration_id=identifier, basal_sha256=digest,
                        calibration=dict(train_ready=True, upper_committing_duration=2)))
                    for factor in (.5, 1., 2.):
                        task = len(width_cases)
                        width_cases.append(dict(task=task, model=model, requested_reference_length=length,
                                                requested_width_label=3.8, integrated_area_factor=factor))
                        path = width / "workers" / f"{task:04d}" / "result.json"
                        path.parent.mkdir(parents=True)
                        path.write_text(json.dumps(dict(basal=dict(checkpoint_sha256=digest),
                            trial=dict(readout_eligible=True, field_tail_pass=False))))
            (calibration / "results/basals.json").write_text(json.dumps(dict(basals=basals)))
            (calibration / "results/calibrations.json").write_text(json.dumps(dict(calibrations=calibrations)))
            (width / "manifest.json").write_text(json.dumps(dict(cases=width_cases)))
            plan = build_plan(calibration, width)
            self.assertEqual(len(plan["prefix_cases"]), 18)
            self.assertEqual(len(plan["cases"]), 270)
            self.assertEqual(len(plan["reused_stage1_results"]), 54)
            self.assertEqual(sum(c["assay"] == "position" for c in plan["cases"]), 54)
            share = [c for c in plan["cases"] if c["assay"] == "share_withdrawal"]
            self.assertEqual(len(share), 162)
            self.assertEqual({c["target_share"] for c in share}, set(SHARES))
            self.assertEqual({c["dwell_factor_of_calibration"] for c in share}, set(DWELL_FACTORS))
            self.assertTrue(all(not r["field_tail_pass"] for r in plan["reused_stage1_results"]))
            checkpoint = Path(basals[0]["checkpoint"])
            checkpoint.write_bytes(b"changed basal")
            with self.assertRaisesRegex(ValueError, "mismatched exact basal"):
                build_plan(calibration, width)


if __name__ == "__main__":
    unittest.main()
