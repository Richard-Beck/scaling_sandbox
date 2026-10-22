"""Finite event, censoring and on-disk native/readout provenance checks."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_readouts_20261020 as readouts
from signature_analysis_20261020 import (
    analyze_completed_workers, analyze_worker_output, audit_native_phase_join, pulse_gap_endpoints, replacement_time_ratio,
    threshold_event, width_ratio_table,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class SignatureAnalysisTests(unittest.TestCase):
    def test_cue_dependent_diagnostic_has_one_sided_switch_but_true_state_cannot_jump(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            active = np.array([[1., 0.], [.5, .5], [0., 1.]])
            joined = root / "joined.npz"
            np.savez(joined, time=np.arange(3.), active=active, field_rac=active,
                     field_rac_activation=np.array([[2., 2.], [2., 2.], [0., 0.]]))
            events = []
            for i in range(2):
                path = root / f"phase{i}.npz"
                np.savez(path, time=np.arange(2.), active=active[i:i + 2], field_rac=active[i:i + 2],
                         field_rac_activation=np.full((2, 2), 2. if i == 0 else 0.))
                events.append(dict(index=i, start=i, stop=i + 1, archive=dict(path=str(path))))
            trace = dict(time=np.arange(3.), active=active)
            result = audit_native_phase_join(joined, events, trace, state_fields=("field_rac",))
            self.assertTrue(result["passed"])
            jump = result["phases"][1]["one_sided_diagnostic_boundaries"]["field_rac_activation"]
            self.assertEqual(jump["maximum_absolute_one_sided_difference"], 2)
            self.assertEqual(jump["joined_trace_keeps"], "preceding_phase_endpoint")
            wrong_rac = active[1:].copy()
            wrong_rac[0, 0] += .01
            np.savez(root / "phase1.npz", time=np.arange(2.), active=active[1:], field_rac=wrong_rac,
                     field_rac_activation=np.zeros((2, 2)))
            with self.assertRaisesRegex(ValueError, "differs from native phase: field_rac"):
                audit_native_phase_join(joined, events, trace, state_fields=("field_rac",))

    def test_native_phase_join_requires_all_fields_samples_and_continuity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            active = np.array([[1., 0.], [.5, .5], [0., 1.]])
            full = dict(time=np.arange(3.), active=active, field_active=active,
                        field_inactive=1 - active)
            joined = root / "joined.npz"
            np.savez(joined, **full)
            events = []
            for i in range(2):
                path = root / f"phase{i}.npz"
                np.savez(path, time=np.arange(2.), active=active[i:i + 2],
                         field_active=active[i:i + 2], field_inactive=1 - active[i:i + 2])
                events.append(dict(index=i, start=i, stop=i + 1, archive=dict(path=str(path))))
            audit = audit_native_phase_join(joined, events, full)
            self.assertTrue(audit["passed"])
            self.assertTrue(audit["all_joined_samples_covered"])
            np.savez(root / "phase1.npz", time=np.arange(2.), active=active[1:],
                     field_active=active[1:], field_inactive=np.ones_like(active[1:]))
            with self.assertRaisesRegex(ValueError, "differs from native phase"):
                audit_native_phase_join(joined, events, full)

    def test_task_selection_never_claims_whole_manifest_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            output = root / "results/empty_subset"
            result = analyze_completed_workers(root, output, task_indices=[])
            self.assertEqual(result["planned_workers"], 0)
            self.assertEqual(result["manifest_workers"], 1)
            self.assertFalse(result["full_brief_complete"])
            self.assertFalse(result["all_native_phase_join_audits_pass"])
            with self.assertRaisesRegex(ValueError, "absent from frozen"):
                analyze_completed_workers(root, output, task_indices=[1])

    def test_cue_on_and_recovery_events_are_separate(self):
        series = dict(time=[0, 1, 2, 3], D_half=[-1, -.2, .6, .8],
                      D_half_eligible=[True] * 4)
        cue = threshold_event(series, 0, 1, .5)
        recovery = threshold_event(series, 1, 3, .5)
        self.assertEqual(cue["status"], "right_censored")
        self.assertEqual(cue["right_censored_at"], 1)
        self.assertEqual(recovery["status"], "observed")
        self.assertEqual(recovery["lower_time"], 1)
        self.assertEqual(recovery["upper_time"], 2)
        self.assertEqual(recovery["lower_elapsed"], 0)
        self.assertAlmostEqual(recovery["estimated_elapsed"], .875)

    def test_missing_samples_cannot_be_finite_censoring(self):
        series = dict(time=[0, 1, 2], D_half=[-1, None, -.1],
                      D_half_eligible=[True, True, True])
        result = threshold_event(series, 0, 2, .5)
        self.assertEqual(result["status"], "unscorable_samples")
        self.assertIsNone(result["right_censored_at"])
        series["D_half"][-1] = .8
        result = threshold_event(series, 0, 2, .5)
        self.assertEqual(result["status"], "observed_after_unscorable_samples")
        self.assertIsNone(result["lower_elapsed"])
        incomplete = threshold_event(series, 0, 1.5, .5)
        self.assertEqual(incomplete["status"], "incomplete_window")

    def test_known_low_activity_is_nonreplacement_not_missing_work(self):
        series = dict(time=[0, 1, 2], D_half=[-1, None, -.1],
                      D_half_eligible=[True, False, True])
        result = threshold_event(series, 0, 2, .5)
        self.assertEqual(result["status"], "right_censored")
        self.assertEqual(result["ineligible_samples"], 1)
        self.assertEqual(result["missing_readout_samples"], 0)
        series["D_half"][-1] = .8
        result = threshold_event(series, 0, 2, .5)
        self.assertEqual(result["status"], "observed")
        self.assertEqual(result["lower_time"], 1)
        self.assertTrue(result["lower_D_unavailable"])
        self.assertIsNone(result["estimated_elapsed"])

    def test_initial_activation_is_flagged_without_inventing_latency(self):
        result = threshold_event(dict(time=[0, 1], D_half=[.6, .8]), 0, 1, .5)
        self.assertEqual(result["status"], "already_above_at_interval_start")
        self.assertIsNone(result["estimated_elapsed"])

    def test_observed_ratio_propagates_both_sampled_brackets(self):
        narrow = dict(status="observed", lower_elapsed=90, upper_elapsed=100,
                      estimated_elapsed=95)
        broad = dict(status="observed", lower_elapsed=20, upper_elapsed=25,
                     estimated_elapsed=22.5)
        ratio = replacement_time_ratio(narrow, broad)
        self.assertEqual(ratio["lower"], 3.6)
        self.assertEqual(ratio["upper"], 5)
        self.assertAlmostEqual(ratio["estimate_from_interpolation"], 95 / 22.5)

    def test_censored_ratios_remain_bounds(self):
        narrow = dict(status="right_censored", right_censored_at=6000)
        broad = dict(status="observed", lower_elapsed=1490, upper_elapsed=1500)
        ratio = replacement_time_ratio(narrow, broad)
        self.assertEqual(ratio["lower"], 4)
        self.assertTrue(ratio["lower_strict"])
        self.assertIsNone(ratio["upper"])
        self.assertIsNone(ratio["estimate_from_interpolation"])
        inverse = replacement_time_ratio(broad, narrow)
        self.assertEqual(inverse["lower"], 0)
        self.assertEqual(inverse["upper"], .25)
        self.assertTrue(inverse["upper_strict"])
        absent = replacement_time_ratio(dict(status="missing"), broad)
        self.assertEqual(absent["status"], "unavailable")
        both = replacement_time_ratio(narrow, narrow)
        self.assertFalse(both["informative"])

    def test_gap_retention_has_an_explicit_gain_denominator(self):
        series = dict(time=[0, 1, 2, 3], share_half=[0, .715, .675, None],
                      D_half=[-1, .43, .35, None], E_challenger_half=[0, 4, 3.5, 0])
        events = [dict(index=0, name="00_on", start=0, stop=1),
                  dict(index=1, name="01_gap", start=1, stop=2),
                  dict(index=2, name="recovery", start=2, stop=3)]
        rows = pulse_gap_endpoints(series, events)
        self.assertAlmostEqual(rows[1]["retention_ratio"], .675 / .715)
        self.assertAlmostEqual(rows[1]["D_change"], -.08)
        self.assertAlmostEqual(rows[1]["challenger_amount_change"], -.5)
        self.assertTrue(rows[1]["primary_retention_denominator_eligible"])
        self.assertTrue(rows[1]["primary_retention_success"])
        self.assertAlmostEqual(rows[1]["primary_retention_ratio"], .675 / .715)
        series["share_half"][1] = 0
        self.assertFalse(pulse_gap_endpoints(series, events)[1]["retention_eligible"])

    def test_zero_gap_train_keeps_each_requested_gap_endpoint(self):
        series = dict(time=[0, 1, 2], share_half=[0, .7, .8],
                      D_half=[-1, .4, .6], E_challenger_half=[0, 3, 4])
        events = [dict(index=0, name="00_on", start=0, stop=1),
                  dict(index=1, name="01_on", start=1, stop=2)]
        rows = pulse_gap_endpoints(series, events)
        self.assertEqual([row["kind"] for row in rows], ["on", "gap", "on"])
        self.assertTrue(rows[1]["zero_duration"])
        self.assertEqual(rows[1]["retention_ratio"], 1)
        self.assertEqual(rows[1]["start"], rows[1]["end"])

    def test_known_low_activity_gap_is_retention_loss_with_observed_denominator(self):
        series = dict(time=[0, 1, 2], share_half=[0, .7, None],
                      D_half=[-1, .4, None], E_challenger_half=[0, 3, 0],
                      readout_eligible=[True, True, False])
        events = [dict(index=0, name="00_on", start=0, stop=1),
                  dict(index=1, name="01_gap", start=1, stop=2)]
        gap = pulse_gap_endpoints(series, events)[1]
        self.assertEqual(gap["retention_evidence_state"], "observed_activity_loss")
        self.assertTrue(gap["retention_probability_denominator_eligible"])
        self.assertTrue(gap["observed_retention_loss"])
        self.assertIsNone(gap["retention_ratio"])
        self.assertTrue(gap["primary_retention_denominator_eligible"])
        self.assertFalse(gap["primary_retention_success"])
        series["readout_eligible"][-1] = True
        missing = pulse_gap_endpoints(series, events)[1]
        self.assertEqual(missing["retention_evidence_state"], "missing")
        self.assertFalse(missing["retention_probability_denominator_eligible"])
        self.assertFalse(missing["primary_retention_denominator_eligible"])
        self.assertIsNone(missing["primary_retention_success"])

    def test_primary_retention_requires_nonnegative_gap_D_and_threshold_sensitivity(self):
        series = dict(time=[0, 1, 2], share_half=[.1, .7, .56],
                      D_half=[-.8, .4, .12], readout_eligible=[True] * 3)
        events = [dict(name="00_on", start=0, stop=1), dict(name="01_gap", start=1, stop=2)]
        gap = pulse_gap_endpoints(series, events)[1]
        self.assertTrue(gap["primary_retention_success_by_threshold"]["0.7"])
        self.assertTrue(gap["primary_retention_success_by_threshold"]["0.8"])
        self.assertFalse(gap["primary_retention_success_by_threshold"]["0.9"])
        series["share_half"][-1] = .49
        series["D_half"][-1] = -.02
        self.assertFalse(pulse_gap_endpoints(series, events)[1]["primary_retention_success"])
        series["D_half"][1] = -.1
        self.assertFalse(pulse_gap_endpoints(series, events)[1]["primary_retention_pulse_eligible"])

    def fixture(self, root, schema=1):
        worker = root / "workers/0000"
        folder = worker / "trial"
        folder.mkdir(parents=True)
        source_name = "polarity1d/manuscript/code/mechanisms/signature_readouts_20261020.py"
        source = root / "frozen_sources" / source_name
        source.parent.mkdir(parents=True)
        source.write_bytes(Path(readouts.__file__).read_bytes())
        case = dict(task=0, model="wave_pinning", stage="width", native_length=2.,
                    requested_reference_length=2, cells=3, cue_width=.95, cue_area=.38,
                    requested_width_label=.95, integrated_area_factor=1)
        x = np.array([1 / 3, 1, 5 / 3])
        edges = np.r_[x[0] - .5 * (x[1] - x[0]), x + .5 * (x[1] - x[0])]
        time = np.array([0., 1., 2.])
        active = np.array([[2., 0., 0.], [0., 0., 2.], [0., 0., 2.]])
        trace = dict(time=time, x=x, active=active, sample_kind="cells", reference_edges=edges)
        basal = dict(time=[0], x=x, active=active[:1], sample_kind="cells", reference_edges=edges)
        basal_path = root / "inputs/basal.npz"
        basal_path.parent.mkdir()
        np.savez(basal_path, time=[0], x=x, active=active[:1])
        np.savez(folder / "trace.npz", time=time, x=x, active=active)
        parameters = dict(challenger_own_intervals=[(1.05, 2)],
                          challenger_fixed_intervals=[(0, 2)], challenger_own_width_physical=.95,
                          challenger_side="right", denominator_atol=1e-12, contrast_atol=.05)
        measurement = readouts.analyze_trace(trace, parameters, basal)
        (folder / "readouts.json").write_text(json.dumps(measurement))
        tail = dict(field_tail_pass=True, available=True)
        (folder / "native_tail.json").write_text(json.dumps(tail))
        manifest = dict(cases=[case], source_sha256={source_name: sha(source)}, result_schema=schema)
        (root / "manifest.json").write_text(json.dumps(manifest))
        (root / "freeze.json").write_text(json.dumps(dict(manifest_sha256=sha(root / "manifest.json"))))
        trial = dict(final_d=measurement["final_d"], readout_eligible=True, ineligibility_reason=None,
                     field_tail_pass=True, field_tail=tail, schedule=dict(actual_total_on=1),
                     events=[dict(index=0, name="cue_on", start=0., stop=1.),
                             dict(index=1, name="recovery", start=1., stop=2.)],
                     readouts="readouts.json", full_trace_sha256=sha(folder / "trace.npz"))
        result = dict(case=case, result_schema=schema, trial=trial,
                      basal=dict(checkpoint=str(basal_path), checkpoint_sha256=sha(basal_path)),
                      manifest_sha256=sha(root / "manifest.json"))
        if schema == 2:
            result["artifact_sha256"] = {str(p.relative_to(root)): sha(p) for p in folder.iterdir()}
            trial["artifacts"] = {key: dict(path=str((folder / name).relative_to(root)), sha256=sha(folder / name))
                                  for key, name in (("readouts", "readouts.json"), ("trace", "trace.npz"),
                                                    ("native_tail", "native_tail.json"))}
            trial["trial_directory"] = str(folder.relative_to(root))
            trial["phase_on_intervals"] = [[0., 1.]]
            trial["recovery_interval"] = [1., 2.]
        (worker / "result.json").write_text(json.dumps(result))
        return case, result, folder

    def test_legacy_unhashed_readouts_are_numerically_audited(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, result, folder = self.fixture(root)
            audit = analyze_worker_output(case, result, run_root=root)
            trial = audit["trials"][0]
            self.assertTrue(trial["native_readout_numeric_audit"])
            self.assertFalse(trial["original_readout_hash_available"])
            self.assertEqual(trial["continuous_cue_on_events"]["D0_5"]["status"], "observed")
            self.assertEqual(trial["recovery_events"]["D0_5"]["status"], "already_above_at_interval_start")
            self.assertIn("workers/0000/trial/readouts.json", audit["provenance"]["files_without_original_declared_hash"])
            measurement = json.loads((folder / "readouts.json").read_text())
            measurement["series"]["E_total"][1] += .1
            (folder / "readouts.json").write_text(json.dumps(measurement))
            with self.assertRaisesRegex(ValueError, "Recomputed numeric"):
                analyze_worker_output(case, result, run_root=root)

    def test_strict_schema_verifies_hashes_and_complete_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, result, folder = self.fixture(root, schema=2)
            audit = analyze_worker_output(case, result, run_root=root)
            self.assertTrue(audit["trials"][0]["original_readout_hash_available"])
            (folder / "extra_native.npz").write_bytes(b"unlisted science")
            with self.assertRaisesRegex(ValueError, "lack hash coverage"):
                analyze_worker_output(case, result, run_root=root)
            (folder / "extra_native.npz").unlink()
            (folder / "trace.npz").write_bytes(b"changed native")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                analyze_worker_output(case, result, run_root=root)

    def test_missing_output_is_an_error_not_censoring(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, result, folder = self.fixture(root)
            (folder / "trace.npz").unlink()
            with self.assertRaises(FileNotFoundError):
                analyze_worker_output(case, result, run_root=root)

    def test_width_pairing_retains_missing_conditions_and_basal_identity(self):
        event = dict(status="observed", lower_elapsed=10, upper_elapsed=12)
        def worker(width, task, basal="same"):
            case = dict(model="wave_pinning", requested_reference_length=19,
                        stage="width", integrated_area_factor=1, requested_width_label=width)
            return dict(case=case, task=task, basal_sha256=basal,
                        trials=[dict(label="trial", continuous_cue_on_events=dict(D0_5=event))])
        self.assertEqual(width_ratio_table([worker(.95, 1)])[0]["ratio"]["status"], "unavailable")
        self.assertEqual(width_ratio_table([worker(.95, 1), worker(7.6, 2)])[0]["ratio"]["lower"], 10 / 12)
        with self.assertRaisesRegex(ValueError, "different exact basals"):
            width_ratio_table([worker(.95, 1), worker(7.6, 2, "other")])

    def test_metadata_stream_preserves_subset_and_missing_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, result, folder = self.fixture(root)
            output = root / "results/analysis"
            summary = analyze_completed_workers(root, output)
            self.assertTrue(summary["declared_subset_complete"])
            self.assertFalse(summary["full_width_factor_coverage"])
            self.assertFalse(summary["full_brief_complete"])
            self.assertTrue((output / "workers/0000.json").is_file())
            (root / "workers/0000/result.json").unlink()
            missing = analyze_completed_workers(root, root / "results/missing")
            self.assertEqual(missing["missing_tasks"], [0])
            self.assertFalse(missing["declared_subset_complete"])
            self.assertEqual(json.loads((root / "results/missing/event_rows.json").read_text()), [])


if __name__ == "__main__":
    unittest.main()
