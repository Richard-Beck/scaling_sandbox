"""Exercise scientific coverage and the concrete-candidate review boundary."""
from pathlib import Path
import gzip
import json
import sys
import tempfile
import unittest

from pypdf import PdfWriter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
from signature_figure4_package_20261020 import (
    REQUIRED_EVIDENCE, SECTIONS, assemble, finalize, sha, write_json,validate_atlas_coverage,
    validate_no_clean_mapping_memo,
)


class PackageReviewTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.receipt = self.root / "scientific_receipt.json"
        write_json(self.receipt, dict(passed=True, complete=True))
        self.history = self.root / "historical.json"
        write_json(self.history, dict(original_full_gate_passed=False, stable_scope_accepted=True,
            gate_counts=dict(all_q90_excluded=18, stable_scope_passed=180)))
        self.pdf = self.root / "source.pdf"
        writer = PdfWriter(); writer.add_blank_page(width=200, height=100)
        with self.pdf.open("wb") as stream: writer.write(stream)
        self.png = self.root / "source.png"; self.png.write_bytes(b"fixed reviewed pixel fixture")
        self.csv = self.root / "page.csv"; self.csv.write_text("measurement,status\n1,pass\n")
        self.spec_path = self.root / "inputs.json"
        self.spec = dict(schema="figure4_reviewable_package_v1",
            evidence={key:dict(self.record(self.receipt), assertions=dict(passed=True, complete=True), scope=key)
                      for key in REQUIRED_EVIDENCE},
            sections={key:[dict(pdf=self.record(self.pdf), page=1, csv=self.record(self.csv), title=key,
                               csv_scope="plotted_values",
                               page_metadata=[self.record(self.receipt)] if key=="S11" else [],
                               display_audit=dict(self.record(self.receipt), assertions=dict(passed=True), scope="Displayed-value audit"))]
                      for key in SECTIONS},
            historical_numerics=dict(self.record(self.history), assertions=dict(stable_scope_accepted=True),
                                     scope="Retained failed full gate and accepted stable subset"),
            delivered_interim=self.record(self.pdf),
            main=dict(png=self.record(self.png), pdf=self.record(self.pdf), framing="Reviewed test fixture",no_clean_mapping=False,
                      display_audit=dict(self.record(self.receipt),assertions=dict(passed=True),scope="All plotted main values equal audited sources"),
                      panels=[dict(id="A", kind="native_small_large_kymographs", measurement="Native profiles", source="fixture",
                                   csv=self.record(self.csv),csv_scope="native_plotted_values")]),
            content=dict(title="Figure 4 fixture", role="Exercise export boundaries",
                         legend=["A scientific measurement"], narrative=["A qualified observed result"],
                         measurements=["A declared estimand"], limitations=["A retained limitation"]))
        audit=self.root/"independent_decision_audit.json"
        write_json(audit,dict(passed=True,producer_complete_sha256=sha(self.receipt)))
        self.spec["evidence"]["independent_final_decision_audit"]=dict(self.record(audit),assertions=dict(passed=True),
            scope="Independent decision arithmetic and native/bootstrap/fingerprint replay")
        self.spec["final_statistics"]=dict(producer=self.record(self.receipt))

    def tearDown(self): self.temporary.cleanup()

    def record(self, path): return dict(path=str(path), sha256=sha(path))

    def build(self):
        native=self.spec["sections"]["S11"][0]
        index=self.root/"native_page_index.json"
        raw=dict(pdf_volume=Path(native["pdf"]["path"]).name,pdf_page=native["page"],
            native_field_CSV_shards=native.get("csv_shards") or [native["csv"]],
            selection_CSV=native["page_metadata"][0],companions={})
        write_json(index,dict(summary=dict(passed=True),pages=[raw]))
        master=self.root/"full_atlas_coverage.json"
        write_json(master,dict(passed=True,all_required_contexts_complete=True,
            page_data_indices=[dict(self.record(index),pages=1)]))
        native["atlas_index_sha256"]=sha(index)
        self.spec["atlas_coverage"]=self.record(master)
        write_json(self.spec_path, self.spec)
        return assemble(self.spec_path, self.root / "output", self.root)

    def test_missing_analysis_arm_refuses_candidate(self):
        del self.spec["evidence"]["E4"]
        with self.assertRaisesRegex(ValueError, "Scientific gates"):
            self.build()
        self.assertFalse((self.root / "output").exists())

    def test_candidate_precedes_review_and_finalization_preserves_pixels(self):
        result = self.build()
        metadata = self.root / "output/supporting/metadata"
        self.assertEqual(result["status"], "candidate")
        self.assertEqual(len(result["pages"]), 13)
        self.assertEqual(result["coverage"]["native_atlas_pages"],1)
        self.assertEqual(result["coverage"]["atlas_index_pages"],1)
        self.assertEqual(result["coverage"]["appendix_pages"],1)
        candidate_content=json.loads(Path(result["candidate_content"]["path"]).read_text())
        self.assertEqual(len(candidate_content["supplement_volumes"]),1)
        self.assertEqual(len(candidate_content["main_source_data"]),1)
        self.assertEqual(candidate_content["main_source_data"][0]["rows"],1)
        self.assertNotEqual(candidate_content["supplement_pdf"],candidate_content["supplement_volumes"][0]["path"])
        self.assertFalse((metadata / "content.json").exists())
        artifact_hashes = {record["path"]:record["sha256"] for record in result["artifacts"]}
        candidate_path = metadata / "candidate_export_audit.json"
        review = self.root / "review.json"
        approval = dict(approved_for_final_export=True, candidate_export_audit_sha256=sha(candidate_path),
            candidate_content_sha256=result["candidate_content"]["sha256"], scientific_qualifications_retained=True,
            main_and_supplement_visual_review_passed=True)
        write_json(review, dict(approval, candidate_content_sha256="wrong"))
        with self.assertRaisesRegex(ValueError, "exact candidate"):
            finalize(candidate_path, dict(self.record(review), repo_path="review.json"))
        write_json(review, approval)
        final = finalize(candidate_path, dict(self.record(review), repo_path="review.json"))
        self.assertTrue(final["final_export_eligible"])
        self.assertEqual(json.loads((metadata / "content.json").read_text())["status"], "final")
        self.assertTrue(all(sha(path)==expected for path, expected in artifact_hashes.items()))

    def test_late_changed_csv_cannot_be_finalized_by_stale_review(self):
        result = self.build()
        Path(result["pages"][0]["csv"]["path"]).write_text("measurement,status\n2,changed\n")
        with self.assertRaisesRegex(ValueError, "changed source"):
            finalize(self.root / "output/supporting/metadata/candidate_export_audit.json", {})

    def test_native_value_shards_are_verified_without_copying_or_dead_public_links(self):
        shards=[]
        for index in range(2):
            source=self.root/f"native_part{index}.csv.gz"
            with gzip.open(source,"wt") as stream:stream.write("time,C,u0,u1\n0,20,0.12345678901234566,1\n5,20,0.5,2\n")
            shards.append(dict(self.record(source),rows=2,distribution="HPC_native_archive"))
        page=self.spec["sections"]["S11"][0]
        page.pop("csv");page.update(csv_shards=shards,csv_scope="native_plotted_values")
        result=self.build()
        selected=next(p for p in result["pages"] if p["section"]=="S11")
        self.assertEqual(selected["rows"],4)
        self.assertEqual(len(selected["native_archive_csvs"]),2)
        self.assertEqual(selected["csv_files"],[])
        self.assertFalse(list((self.root/"output").rglob("*.csv.gz")))
        reader=(self.root/"output/figure4_complete_reader.html").read_text()
        self.assertIn("Full native values: HPC archive inventory",reader)
        self.assertNotIn('href="'+str(self.root/"native_part0.csv.gz"),reader)
        with gzip.open(self.root/"native_part1.csv.gz","wt") as stream:stream.write("changed\n")
        with self.assertRaisesRegex(ValueError,"changed source"):
            finalize(self.root/"output/supporting/metadata/candidate_export_audit.json",{})

    def test_atlas_volume_cannot_silently_omit_a_native_page(self):
        atlas=self.root/"two_native_pages.pdf"
        writer=PdfWriter();writer.add_blank_page(width=200,height=100);writer.add_blank_page(width=200,height=100)
        with atlas.open("wb") as stream:writer.write(stream)
        self.spec["sections"]["S11"][0]["pdf"]=self.record(atlas)
        with self.assertRaisesRegex(ValueError,"Every S11 appendix source page"):
            self.build()

    def test_whole_other_atlas_context_cannot_be_omitted(self):
        self.build()
        additional=self.root/"other_context_index.json"
        source=json.loads((self.root/"native_page_index.json").read_text())
        source["pages"][0]["pdf_volume"]="other_context.pdf"
        write_json(additional,source)
        master=self.root/"full_atlas_coverage.json"
        payload=json.loads(master.read_text())
        payload["page_data_indices"].append(dict(self.record(additional),pages=1))
        write_json(master,payload);self.spec["atlas_coverage"]=self.record(master)
        with self.assertRaisesRegex(ValueError,"entire native atlas volume or context"):
            validate_atlas_coverage(self.spec,self.root)

    def test_main_numeric_sources_cannot_be_replaced_by_a_source_label(self):
        self.spec["main"]["panels"][0].pop("csv")
        with self.assertRaisesRegex(ValueError,"exact plotted numeric CSV"):
            self.build()

    def test_no_clean_mapping_cannot_be_final_candidate_without_reviewed_S0_memo(self):
        self.spec['main']['no_clean_mapping']=True
        with self.assertRaisesRegex(ValueError,'mandatory reviewed S0 memo'):
            self.build()

    def test_reviewed_memo_requires_exact_statuses_and_all_displayed_pages(self):
        models=['wave_pinning','goryachev','holmes_model3','spring_ring']
        rows=[dict(arm='common',view='main',forced_inclusive=False,model=model,criterion='H1_model',status='inconclusive') for model in models]
        import csv
        criteria=self.root/'criteria.csv'
        with criteria.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        producer=self.root/'producer.json';write_json(producer,dict(passed=True,complete=True,outputs={'criteria.csv':sha(criteria)}))
        producer_record=self.record(producer)
        independent=self.root/'decision_audit.json';write_json(independent,dict(passed=True,producer_complete_sha256=sha(producer)))
        audit_record=self.record(independent)
        memo=self.root/'reviewed_memo.json'
        obj=dict(approved_for_candidate_framing=True,no_clean_mapping=True,framing=self.spec['main']['framing'],
            producer_complete_sha256=sha(producer),independent_final_decision_audit_sha256=sha(independent),
            overall_interpretation='Qualified fixture statement',retained_qualifications='Finite bounds',
            model_summaries=[dict(model=model,pattern='Pattern',numeric_shortfalls='Deciding bounds',qualifications='Windows',
                                 criterion_refs=[row]) for model,row in zip(models,rows)],
            alternative_framings=[dict(framing=f'Alternative {i}',supporting_panels=['S2'],weaknesses='Finite windows',
                further_runs='Longer observations') for i in range(2)])
        write_json(memo,obj)
        display=self.root/'memo_display.json';write_json(display,dict(passed=True,complete=True,
            source_files={'reviewed_interpretation_memo':self.record(memo)},pdf=dict(self.record(self.pdf),pages=1),
            pages=[dict(pdf_page=1,csv_sha256=sha(self.csv))]))
        self.spec['main']['no_clean_mapping']=True
        self.spec['reviewed_interpretation_memo']=dict(source=self.record(memo),display_audit=self.record(display))
        memo_page=dict(pdf=self.record(self.pdf),page=1,csv=self.record(self.csv),
            display_audit=dict(self.record(display),assertions=dict(passed=True,complete=True),scope='Reviewed S0 memo'))
        self.spec['sections']['S0'].append(memo_page)
        result=validate_no_clean_mapping_memo(self.spec,self.root,producer_record,audit_record)
        self.assertEqual(set(result),{'reviewed_interpretation_memo','reviewed_interpretation_memo_display'})
        self.spec['sections']['S0'].pop()
        with self.assertRaisesRegex(ValueError,'Every mandatory S0 memo page'):
            validate_no_clean_mapping_memo(self.spec,self.root,producer_record,audit_record)
        self.spec['sections']['S0'].append(memo_page)
        obj['model_summaries'][0]['criterion_refs'][0]=dict(rows[0],status='pass');write_json(memo,obj)
        self.spec['reviewed_interpretation_memo']['source']=self.record(memo)
        with self.assertRaisesRegex(ValueError,'different criterion status'):
            validate_no_clean_mapping_memo(self.spec,self.root,producer_record,audit_record)

    def test_changed_original_main_asset_refuses_stale_finalization(self):
        self.build()
        self.png.write_bytes(b"different native example")
        with self.assertRaisesRegex(ValueError,"changed source"):
            finalize(self.root/"output/supporting/metadata/candidate_export_audit.json",{})


if __name__ == "__main__": unittest.main()
