"""Receipt-defined publication keeps superseded evidence local and intact."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_figure4_publication_inventory_20261020 as D
import signature_figure4_package_20261020 as Q


class PublicationInventoryTests(unittest.TestCase):
    def record(self, path):
        return dict(path=str(path), sha256=Q.sha(path))

    def fixture(self, root):
        topic = root / "topic"; topic.mkdir()
        data = topic / "final_numeric.csv"; data.write_text("x,y\n1,2\n")
        old = topic / "figure4_S7_E1_v1_20261020.pdf"; old.write_bytes(b"retained layout v1")
        guide = topic / "README.md"; guide.write_text("Current navigation\n")
        content = topic / "content.json"; content.write_text('{"status":"final"}\n')
        candidate = topic / "candidate.json"; candidate.write_text('{"status":"candidate"}\n')
        review = root / "review.json"; review.write_text('{"approved_for_final_export":true}\n')
        final = topic / "final.json"
        final.write_text(json.dumps(dict(passed=True, status="final", final_export_eligible=True,
            artifacts=[self.record(data)], candidate_export=self.record(candidate),
            completion_review=self.record(review), content=self.record(content))))
        archive = root / "immutable_E0"; archive.mkdir()
        (archive / "README.md").write_text("Delivered navigation\n")
        (archive / "accepted.csv").write_text("value\n3\n")
        accepted = topic / "accepted.csv"; accepted.write_bytes((archive / accepted.name).read_bytes())
        receipt = archive / "delivery_receipt.json"
        receipt.write_text(json.dumps(dict(source_topic=str(topic), artifacts={
            "README.md": Q.sha(archive / "README.md"), "accepted.csv": Q.sha(archive / "accepted.csv")})))
        return topic, final, receipt, guide, old

    def test_explicit_allowlist_preserves_and_excludes_superseded_layout(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); topic, final, receipt, guide, old = self.fixture(root)
            before = Q.sha(old)
            report = D.inventory(final, [receipt], topic, root, topic / "publication", [guide])
            published = {row["repo_path"] for row in report["portable_records"]}
            self.assertIn("topic/final_numeric.csv", published)
            self.assertIn("topic/accepted.csv", published)
            self.assertIn("topic/README.md", published)
            self.assertNotIn("topic/figure4_S7_E1_v1_20261020.pdf", published)
            self.assertEqual(Q.sha(old), before)
            self.assertEqual(report["deleted_files"], 0)
            self.assertEqual(len(report["immutable_archived_navigation_retained_local"]), 1)
            self.assertIn(b"topic/publication.json\0", Path(report["git_pathspec"]["path"]).read_bytes())
            self.assertEqual(sum(group['logical_bytes'] for group in report['portable_storage_groups'].values()),report['portable_bytes'])
            self.assertEqual(sum(group['allocated_bytes'] for group in report['portable_storage_groups'].values()),report['portable_allocated_bytes'])
            self.assertGreater(report['portable_storage_groups']['accepted_historical']['files'],0)
            exact_paths=[root/path.decode() for path in Path(report['git_pathspec']['path']).read_bytes().split(b'\0') if path]
            self.assertEqual(report['complete_Git_allowlist_files'],len(exact_paths))
            self.assertEqual(report['complete_Git_allowlist_logical_bytes'],sum(path.stat().st_size for path in exact_paths))
            self.assertEqual(report['complete_Git_allowlist_allocated_bytes'],sum(path.stat().st_blocks*512 for path in exact_paths))

    def test_candidate_cannot_be_labelled_ready_for_git_handoff(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); topic, final, receipt, _, _ = self.fixture(root)
            payload = json.loads(final.read_text()); payload["status"] = "candidate"; final.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "actual final eligible"):
                D.inventory(final, [receipt], topic, root, topic / "publication")

    def test_changed_accepted_scientific_data_cannot_be_silently_omitted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); topic, _, receipt, _, _ = self.fixture(root)
            (topic / "accepted.csv").write_text("value\n4\n")
            with self.assertRaises(ValueError):
                D.accepted_records(receipt)

    def test_mutable_guide_cannot_link_a_superseded_local_only_plot(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); topic, final, receipt, guide, old = self.fixture(root)
            guide.write_text("[Superseded local plot](" + old.name + ")\n")
            with self.assertRaisesRegex(ValueError, "unapproved local-only"):
                D.inventory(final, [receipt], topic, root, topic / "publication", [guide])
            self.assertTrue(old.exists())

    def test_actual_organizer_directory_links_require_public_descendants(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);topic,final,receipt,guide,_=self.fixture(root)
            rows=[]
            for name in ('vectors','data','metadata'):
                folder=topic/'supporting'/name;folder.mkdir(parents=True)
                item=folder/'accepted.txt';item.write_text('Accepted public fixture\n');rows.append(self.record(item))
            obj=json.loads(final.read_text());obj['artifacts'].extend(rows);final.write_text(json.dumps(obj))
            guide.write_text('[vectors](supporting/vectors/) · [data](supporting/data/) · [metadata and audits](supporting/metadata/)\n')
            report=D.inventory(final,[receipt],topic,root,topic/'publication',[guide])
            self.assertEqual(report['public_links_checked'],3)
            (topic/'viewer').mkdir();(topic/'viewer/local_only.html').write_text('Local superseded example')
            guide.write_text('[Viewer](viewer/)\n')
            with self.assertRaisesRegex(ValueError,'local-only directory'):
                D.inventory(final,[receipt],topic,root,topic/'publication',[guide])
            guide.write_text('[Missing](missing/)\n')
            with self.assertRaisesRegex(ValueError,'unapproved local-only or missing'):
                D.inventory(final,[receipt],topic,root,topic/'publication',[guide])


if __name__ == "__main__":
    unittest.main()
