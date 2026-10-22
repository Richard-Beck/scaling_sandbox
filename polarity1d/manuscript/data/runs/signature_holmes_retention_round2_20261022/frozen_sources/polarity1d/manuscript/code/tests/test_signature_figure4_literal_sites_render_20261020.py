"""Bind literal-coordinate displays to every audited native group/window."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_figure4_literal_sites_render_20261020 as D
import signature_figure4_interim_render_20261020 as P


class LiteralDisplayTests(unittest.TestCase):
    def fixture(self, root):
        (root / "results").mkdir(); (root / "workers/0000").mkdir(parents=True); (root / "audits").mkdir()
        source = root / "frozen_source.py"; source.write_text("# reviewed scientific source\n")
        cases = [dict(id="seed1", scalar_sha256="a", profile_sha256="b", native_horizon=32000.)]
        group = dict(task=0, model="spring_ring", C_label=3., f=.1, sigma=1., native_group_sha256="group-sha", noise=cases)
        manifest = root / "manifest.json"
        manifest.write_text(json.dumps(dict(groups=[group], views={view: [] for view in D.VIEWS},
            source_version=2, antipodal_policies=list(D.POLICIES), antipodal_atol_radians=1e-12,
            source_sha256={str(source): P.sha(source)}, input_receipts={})))
        ambiguity_rows = [dict(id="seed1", view=view, method=method, cadence=5., qualified_frames=2,
                               total_frames=2, qualified_adjacent_steps=1, antipodal_steps=1,
                               tracks_before_split=1, tracks_after_split=2)
                          for view in D.VIEWS for method in D.METHODS]
        qualifications = [dict(view=view, method=method, antipodal_policy=policy, cadence=5.,
                               qualified_frames=2, total_frames=2, antipodal_steps=1, seeds_with_antipodal_steps=1)
                          for view in D.VIEWS for method in D.METHODS for policy in D.POLICIES]
        array = root / "workers/0000/arrays.npz"; array.write_bytes(b"immutable independently audited arrays")
        result = root / "workers/0000/result.json"
        result.write_text(json.dumps(dict(passed=True, task=0, manifest_sha256=P.sha(manifest),
            arrays_sha256=P.sha(array), cases=cases, whole_seeds=1, bootstrap_resamples=2000,
            qualifications=qualifications, case_ambiguities=ambiguity_rows)))
        audit = root / "audits/0000.json"
        audit.write_text(json.dumps(dict(passed=True, task=0, source_version=2, manifest_sha256=P.sha(manifest),
            result_sha256=P.sha(result), arrays_sha256=P.sha(array))))
        rows = [dict(model="spring_ring", C_label=3., f=.1, sigma=1., native_horizon=32000.,
                     native_group_sha256="group-sha", view=view, method=method, antipodal_policy=policy,
                     primary_A3=method == D.METHODS[0] and policy == "split", antipodal_steps=1,
                     seeds_with_antipodal_steps=1, qualified_frames=2,
                     n_seeds=1, bootstrap_resamples=2000, lag=5., cadence=5., MSD="nan", coverage=0.)
                for view in D.VIEWS for method in D.METHODS for policy in D.POLICIES]
        curve = root / "results/site_curves.csv.gz"
        with gzip.open(curve, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        ambiguities = root / "results/case_ambiguities.csv.gz"
        expanded = [dict(row, model="spring_ring", C_label=3., f=.1, sigma=1., native_horizon=32000.,
                         native_group_sha256="group-sha") for row in ambiguity_rows]
        with gzip.open(ambiguities, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(expanded[0])); writer.writeheader(); writer.writerows(expanded)
        status = root / "results/status.json"
        status.write_text(json.dumps(dict(passed=True, complete=True, manifest_sha256=P.sha(manifest),
            condition_groups=1, no_native_integrations=True, no_H1_H2_H3_estimand_change=True,
            source_version=2, antipodal_policies=list(D.POLICIES), antipodal_atol_radians=1e-12,
            all_qualified_native_frames_preserved=True, route_choices_are_not_independent_seeds=True,
            methods=list(D.METHODS), rows=len(rows), outputs={curve.name: P.sha(curve), ambiguities.name: P.sha(ambiguities)},
            receipts=[dict(task=0, result_sha256=P.sha(result), audit_sha256=P.sha(audit))])))
        return status

    def test_lost_heading_stays_missing_and_primary_is_literal_heading(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); status = self.fixture(root)
            rows, ambiguities, _ = D.load_literal(root, P.sha(status))
            self.assertEqual(len(rows), 6)
            self.assertEqual(len(ambiguities), 14)
            self.assertTrue(all(row["MSD"] == "nan" for row in rows))
            self.assertEqual(sum(P.boolean(row["primary_A3"]) for row in rows), 1)

    def test_tied_route_is_never_a_primary_or_independent_seed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); status = self.fixture(root)
            claim = json.loads(status.read_text()); claim["route_choices_are_not_independent_seeds"] = False
            status.write_text(json.dumps(claim))
            with self.assertRaisesRegex(ValueError, "complete audited"):
                D.load_literal(root, P.sha(status))

    def test_sidecar_values_cannot_change_even_when_rehashed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); status = self.fixture(root); path = root / "results/case_ambiguities.csv.gz"
            with gzip.open(path, "rt", newline="") as stream: rows = list(csv.DictReader(stream))
            rows[0]["antipodal_steps"] = 0
            with gzip.open(path, "wt", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            claim = json.loads(status.read_text()); claim["outputs"][path.name] = P.sha(path); status.write_text(json.dumps(claim))
            with self.assertRaisesRegex(ValueError, "changed an audited native case value"):
                D.load_literal(root, P.sha(status))

    def test_a_claimed_complete_status_cannot_omit_a_worker_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); status = self.fixture(root)
            payload = json.loads(status.read_text()); payload["receipts"] = []; status.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "own independent receipt"):
                D.load_literal(root, P.sha(status))

    def test_native_array_or_identity_change_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); status = self.fixture(root)
            (root / "workers/0000/arrays.npz").write_bytes(b"changed")
            with self.assertRaises(ValueError):
                D.load_literal(root, P.sha(status))

    def test_canonical_selection_retains_routes_without_adding_overlapping_seeds(self):
        metric = dict(model="wave_pinning", C_label=20, f=.03, view="main", native_horizon=32000,
                      selected_observation_arm="common_E2_augmented_E0", analysis_arm="common", forced_combined="True")
        rows = [dict(metric, antipodal_policy=policy, id="same_native_seed", n_seeds=12) for policy in D.POLICIES]
        curves, cases = D.select_literal_observations([metric], {"E2": rows, "E3": rows}, {"E2": rows, "E3": rows})
        self.assertEqual(len(curves), 3)
        self.assertEqual({row["id"] for row in cases}, {"same_native_seed"})
        self.assertTrue(all(P.boolean(row["display_forced"]) for row in curves))
        with self.assertRaisesRegex(ValueError, "no exact curve"):
            D.select_literal_observations([dict(metric, selected_observation_arm="preselected_E4", analysis_arm="longest")],
                                         {"E2": rows}, {"E2": rows})
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); status = self.fixture(root)
            path = root / "audits/0000.json"; payload = json.loads(path.read_text())
            payload["manifest_sha256"] = "different native group"; path.write_text(json.dumps(payload))
            claim = json.loads(status.read_text()); claim["receipts"][0]["audit_sha256"] = P.sha(path); status.write_text(json.dumps(claim))
            with self.assertRaisesRegex(ValueError, "identities differ"):
                D.load_literal(root, P.sha(status))


if __name__ == "__main__":
    unittest.main()
