"""Scientific export gates and censoring remain explicit in Figure 2 assembly."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
from signature_figure2_render_20261020 import (
    condition_key, coverage, load_audited_records, ratio_records, ratio_text, regional_peak_vectorized,
    required_conditions,
)
from signature_figure2_supplement_20261020 import display_indices, exact_binomial_interval


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Figure2RenderTests(unittest.TestCase):
    def test_noise_interval_keeps_zero_eligible_and_finite_seed_uncertainty(self):
        self.assertEqual(exact_binomial_interval(0, 0), (None, None))
        low, high = exact_binomial_interval(12, 12)
        self.assertGreater(low, .7)
        self.assertLess(low, .8)
        self.assertEqual(high, 1.)
        low, high = exact_binomial_interval(0, 12)
        self.assertEqual(low, 0.)
        self.assertGreater(high, .2)
        self.assertLess(high, .3)
    def test_supplement_display_keeps_extrema_and_exact_phase_boundaries(self):
        time = np.linspace(0, 100, 10001)
        values = np.zeros_like(time)
        values[1234], values[8765] = 10, -4
        index = display_indices(time, values, [0, 50, 100], bins=10)
        self.assertTrue({0, 1234, 5000, 8765, 10000}.issubset(set(index)))
        with self.assertRaisesRegex(ValueError, "phase endpoint is missing"):
            display_indices(time, values, [50.00001], bins=10)

    def test_reporter_maximum_retains_cell_overlap_and_nodal_interpolated_edges(self):
        active = np.array([[2., 1., 0.], [1., 3., 2.]])
        cells = regional_peak_vectorized(active, [.5, 1.5, 2.5],
            [[(1.9, 3.)], [(2., 3.)]], sample_kind="cells")
        np.testing.assert_allclose(cells, [1., 1.])
        nodes = regional_peak_vectorized(active, [0., 1., 2.],
            [[(1.5, 2.)], [(1.5, 2.)]], sample_kind="nodes")
        np.testing.assert_allclose(nodes, [.5, 1.5])

    def test_full_factor_grid_is_exactly_450_widths_and_1620_trains(self):
        widths, trains = required_conditions()
        self.assertEqual(len(widths), 450)
        self.assertEqual(len(trains), 1620)
        self.assertFalse(coverage({}, [])["final_export_eligible"])

    def test_absent_calibration_accounts_for_a_planned_train_without_inventing_execution(self):
        key = ("train", "holmes_model3", 39, 4, 4., 1.5)
        report = coverage({}, [], calibration_unavailable={key: dict(status="calibration_unavailable")})
        self.assertEqual(report["audited_train_conditions"], 0)
        self.assertEqual(report["calibration_unavailable_train_conditions"], 1)
        self.assertEqual(report["scientifically_accounted_train_conditions"], 1)
        self.assertEqual(len(report["missing_train_conditions"]), 1619)
        self.assertFalse(report["final_export_eligible"])
        with self.assertRaisesRegex(ValueError, "without executed output"):
            coverage({key: {}}, [], calibration_unavailable={key: {}})

    def test_censored_shape_penalty_cannot_become_finite_or_ranked(self):
        records = {}
        for width, event in ((.95, dict(status="right_censored", right_censored_at=6000)),
                             (7.6, dict(status="observed", lower_elapsed=490, upper_elapsed=500))):
            case = dict(stage="width", model="wave_pinning", requested_reference_length=39,
                        requested_width_label=width, integrated_area_factor=1, checkpoint_sha256="basal")
            records[condition_key(case)] = dict(case=case,
                measurement=dict(continuous_cue_on_events=dict(D0_5=event)))
        row = next(r for r in ratio_records(records) if r["model"] == "wave_pinning" and r["length"] == 39)
        self.assertEqual(ratio_text(row["ratio"]), ">12")
        self.assertIsNone(row["ratio"]["upper"])
        self.assertIsNone(row["ratio"]["estimate_from_interpolation"])
        records[("width", "wave_pinning", 39, .95, 1.)]["case"]["checkpoint_sha256"] = "different"
        with self.assertRaisesRegex(ValueError, "exact frozen basal"):
            ratio_records(records)

    def test_audit_table_requires_digest_and_missing_tasks_remain_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "workers").mkdir()
            directory = root / "results/audit"
            directory.mkdir(parents=True)
            case = dict(task=0, stage="width", model="wave_pinning", requested_reference_length=19,
                        requested_width_label=.95, integrated_area_factor=1)
            (root / "manifest.json").write_text(json.dumps(dict(cases=[case])))
            table = directory / "event_rows.json"
            table.write_text("[]")
            status = dict(producer_state="stopped", manifest_sha256=sha(root / "manifest.json"),
                          outputs={"event_rows.json": sha(table)}, missing_tasks=[0], errors=[],
                          declared_subset_complete=False)
            (directory / "status.json").write_text(json.dumps(status))
            records, receipts = load_audited_records([directory])
            self.assertEqual(records, {})
            self.assertEqual(receipts[0]["missing_tasks"], [0])
            report = coverage(records, receipts)
            self.assertEqual(report["audited_width_conditions"], 0)
            self.assertEqual(len(report["missing_width_conditions"]), 450)
            table.write_text("[{}]")
            with self.assertRaisesRegex(ValueError, "verified digest"):
                load_audited_records([directory])

    def test_manually_passed_review_cannot_override_factor_omissions(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            evidence = root / "evidence.json"
            evidence.write_text("{}")
            review = dict(independent_review_passed=True,
                          panel_completion={p: True for p in "ABCDEF"},
                          stage_completion={str(i): True for i in range(6)},
                          evidence=[dict(stage=str(i), path=str(evidence), sha256=sha(evidence)) for i in range(6)])
            path = root / "review.json"
            path.write_text(json.dumps(review))
            self.assertFalse(coverage({}, [], path)["final_export_eligible"])


if __name__ == "__main__":
    unittest.main()
