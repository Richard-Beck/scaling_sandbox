"""Final displays cannot mix scientific versions or omit unfavorable views."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_figure4_decisions_render_20261020 as D
import signature_figure4_interim_render_20261020 as P


def criteria():
    return [dict(arm=arm, view=view, forced_inclusive=inclusive, model=model,
                 criterion=criterion, status="inconclusive", reason="Finite observation bounds")
            for arm in D.ARMS for view in D.VIEWS for inclusive in (False, True)
            for model, criterion in D.criterion_order()]


class DecisionDisplayTests(unittest.TestCase):
    def h3_tables(self):
        effect=dict(estimate=1.,ci_lower=.5,ci_upper=None,ci_upper_unbounded=True,
                    independent_small_seeds=12,independent_large_seeds=12,valid_bootstrap_fraction=.9,available=True)
        strong=dict(criterion="H3_strong_input",status="inconclusive",reason="Finite weak bounds")
        e1=dict(criterion="H3_E1_corroboration",status="inconclusive",reason="Outside preset input grid",sigma_checks=[])
        metrics=[dict(arm=arm,model=model,C_label=size,f=f,view=view,tau_lower=0,tau_upper=10,tau_status="bounded",tau_coverage=.8,
                      multiple_fronts_fraction=.7,multiple_fronts_ci_lower=.6,multiple_fronts_ci_upper=.8,
                      low_contrast_fraction=.1,low_contrast_ci_lower=.05,low_contrast_ci_upper=.2,
                      native_horizon=32000,forced_mean=False,forced_E1=False,forced_combined=False)
                 for arm in D.ARMS for model in P.MODELS for size in (1.,3.) for f in (.0001,.1) for view in D.VIEWS]
        candidates=[dict(arm=arm,model=model,view=view,forced_inclusive=inclusive,f=.1,mode="multiple_fronts",
                         mode_descriptor="Native raw regions",strong=json.dumps(strong),E1=json.dumps(e1),fingerprints="{}")
                    for arm in D.ARMS for model in P.MODELS for view in D.VIEWS for inclusive in (False,True)]
        effects=[dict(arm=arm,model=model,view=view,f=.1,statistic=stat,**effect)
                 for arm in D.ARMS for model in P.MODELS for view in D.VIEWS for stat in ("region_count","normalized_contrast")]
        rows=criteria()
        fp=dict(criterion="H3_fingerprint_separation",status="inconclusive",reason="CI overlap",
                components=[dict(criterion=stat,status="inconclusive",first_ci=[.5,None],second_ci=[.5,None])
                            for stat in ("region_count","normalized_contrast")])
        pair=dict(first_model="wave_pinning",first_f=.1,first_mode="multiple_fronts",second_model="goryachev",second_f=.1,
                  second_mode="low_contrast",criterion="H3_distinct_model_pair",status="inconclusive",reason="Existing audited truth",
                  components=[strong,strong,e1,e1,fp])
        for row in rows:
            if row['criterion']=="H3_overall":row['pairs']=json.dumps([pair,pair])
        return {"criteria.csv":rows,"condition_metrics.csv":metrics,"strong_input_candidates.csv":candidates,
                "fingerprint_effects.csv":effects,"E1_condition_metrics.csv":[]}

    def test_all_h3_candidates_and_unbounded_fingerprints_are_visible(self):
        tables=self.h3_tables();rows,_=D.h3_projection(tables)
        self.assertEqual(len(rows),112)
        self.assertIn("∞",D.effect_text(json.loads(rows[0]['display_raw_fingerprint_effects'])['region_count']))
        self.assertIn("[0.6, 0.8]",D.fraction_text(json.loads(rows[0]['display_large_strong_condition']),'multiple_fronts_fraction'))
        self.assertEqual(D.sigma_text(None,None),"not declared")
        tables['strong_input_candidates.csv'].pop()
        with self.assertRaisesRegex(ValueError,"Every predeclared H3"):
            D.h3_projection(tables)

    def test_every_nested_h3_pair_is_projected_without_selecting_a_best_pair(self):
        tables=self.h3_tables();captured=[]
        with patch.object(D,"scientific_table_pages",side_effect=lambda book,rows,cells,headers,widths,label,title,note,**kwargs:captured.append((label,rows,cells))):
            D.h3_pages(None,tables)
        pair_pages=[row for label,rows,_ in captured if label.startswith('S0_H3_pairs_') for row in rows]
        self.assertEqual(len(pair_pages),56)
        self.assertTrue(all(row['status']=='inconclusive' for row in pair_pages))
        self.assertTrue(all('display_first_raw_effects' in row for row in pair_pages))
        candidate_rows=[row for label,rows,_ in captured if label.startswith('S0_H3_candidates_') for row in rows]
        self.assertEqual(len(candidate_rows),112)
        self.assertIn('2 preset distinct-model/mode pairs',"\n".join(D.detail_lines(next(row for row in tables['criteria.csv'] if row['criterion']=='H3_overall'))))

    def test_review_memo_retains_actual_truth_and_no_clean_mapping_alternatives(self):
        memo=dict(approved_for_candidate_framing=True,producer_complete_sha256='producer',
            independent_final_decision_audit_sha256='audit',framing='Reviewed fixture',no_clean_mapping=True,
            overall_interpretation='Qualified test statement',retained_qualifications='Finite windows',
            model_summaries=[dict(model=model,pattern='Observed fixture pattern',qualifications='Observed limits',
                numeric_shortfalls='Finite deciding bounds',criterion_refs=[dict(arm='common',view='main',
                forced_inclusive=False,model=model,criterion='H1_model',status='inconclusive')]) for model in P.MODELS],
            alternative_framings=[dict(framing=f'Alternative {i}',supporting_panels=['S2'],weaknesses='Finite bounds',
                further_runs='Longer observations') for i in range(2)])
        rows=D.reviewed_memo_rows(memo,{'criteria.csv':criteria()},'producer','audit')
        self.assertEqual(len(rows),7)
        self.assertTrue(all(json.loads(row['actual_criterion_records'])[0]['status']=='inconclusive'
                            for row in rows if row['section']=='model'))
        memo['model_summaries'][0]['criterion_refs'][0]['status']='pass'
        with self.assertRaisesRegex(ValueError,'actual audited truth'):
            D.reviewed_memo_rows(memo,{'criteria.csv':criteria()},'producer','audit')
        memo['model_summaries'][0]['criterion_refs'][0]['status']='inconclusive'
        memo['alternative_framings']=[]
        with self.assertRaisesRegex(ValueError,'two or three'):
            D.reviewed_memo_rows(memo,{'criteria.csv':criteria()},'producer','audit')

    def test_unfavorable_late_window_and_forcing_views_cannot_be_dropped(self):
        rows = criteria()
        D.validate_criteria(rows)
        removed = next(i for i, row in enumerate(rows) if row["view"] == "late_quarter" and row["forced_inclusive"])
        with self.assertRaisesRegex(ValueError, "308 raw criteria"):
            D.validate_criteria(rows[:removed] + rows[removed + 1:])
        with self.assertRaisesRegex(ValueError, "308 raw criteria"):
            D.validate_criteria(rows[:-1] + [rows[0]])

    def test_unknown_or_absent_values_do_not_become_negative_or_infinite(self):
        self.assertEqual(D.words(None), "unavailable")
        self.assertEqual(D.words(False), "no")
        self.assertEqual(D.words(float("inf")), "∞")
        self.assertEqual(D.words([None, 3.]), "unavailable, 3")
        self.assertEqual(D.words([]), "none observed")
        rows = criteria(); rows[0]["status"] = "missing"
        with self.assertRaisesRegex(ValueError, "unknown criterion"):
            D.validate_criteria(rows)

    def fixture(self, root):
        out = root / "results/final_analysis"; out.mkdir(parents=True)
        source = root / "scientific_source.py"; source.write_text("# frozen source\n")
        spec = root / "inputs.json"; spec.write_text("{}\n")
        source_inventory = {str(source): P.sha(source)}
        freeze = root / "source_freeze.json"
        freeze.write_text(json.dumps(dict(source_sha256=source_inventory)))
        active = root / "ACTIVE_FINAL_ANALYSIS.json"
        active.write_text(json.dumps(dict(active_source_freeze=str(freeze), active_source_freeze_sha256=P.sha(freeze))))
        outputs = {}
        for name in D.REQUIRED:
            path = out / name
            with path.open("w", newline="") as stream:
                rows = criteria() if name == "criteria.csv" else [dict(scope="complete fixture")]
                writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            outputs[name] = P.sha(path)
        producer = out / "complete.json"
        producer.write_text(json.dumps(dict(passed=True, complete=True, outputs=outputs,
            source_freeze_sha256=P.sha(freeze), source_sha256=source_inventory,
            input_spec_path=str(spec), input_spec_sha256=P.sha(spec))))
        audit = out / "independent_final_decision_audit.json"
        audit.write_text(json.dumps(dict(passed=True, complete=True,
            producer_complete_sha256=P.sha(producer), final_outputs_sha256=outputs)))
        return producer, audit, active

    def test_exact_independent_producer_binding_and_csv_bytes_are_required(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); producer, audit, _ = self.fixture(root)
            tables, _ = D.load_final(root, P.sha(producer), P.sha(audit))
            self.assertEqual(len(tables["criteria.csv"]), 308)
            payload = json.loads(audit.read_text()); payload["producer_complete_sha256"] = "0" * 64
            audit.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "different final producer"):
                D.load_final(root, P.sha(producer), P.sha(audit))
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); producer, audit, _ = self.fixture(root)
            path = root / "results/final_analysis/criteria.csv"
            path.write_text(path.read_text().replace("inconclusive", "pass", 1))
            with self.assertRaises(ValueError):
                D.load_final(root, P.sha(producer), P.sha(audit))

    def test_old_source_revision_is_refused_even_with_valid_old_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); producer, audit, active = self.fixture(root)
            other = root / "source_freeze_new.json"; other.write_text('{"new":"reviewed"}\n')
            active.write_text(json.dumps(dict(active_source_freeze=str(other), active_source_freeze_sha256=P.sha(other))))
            with self.assertRaisesRegex(ValueError, "active reviewed source"):
                D.load_final(root, P.sha(producer), P.sha(audit))


class PrefixRecoveryDisplayTests(unittest.TestCase):
    def test_repair_census_cannot_become_outcome_selection_or_a_relaxed_tolerance(self):
        summary=dict(passed=True,numerical_integrity_only=True,scientific_long_arm_qualification=False,
            response_outcomes_used_for_selection=False,all_failed_originals_retained=True,no_passing_or_live_originals_rerun=True,
            selected_native_failures=0,missing_native_trajectories=0,selected_native_trajectories=3,actual_prefix_repairs=1,
            original_attempted_native_trajectories=3,original_native_integrations_reaching_full_horizon=3,
            original_prefix_validation_passed=2,original_prefix_validation_failed=1,qualified_old_cap_restarts=1,
            passing_originals_preferred=2,total_native_trajectories_including_failed_attempts=4,
            repair_model_counts=dict(wave_pinning=1),prefix_checks=3,bitwise_equal_prefix_checks=2,
            maximum_absolute_prefix_difference=1e-14,unchanged_native_prefix_rule=dict(relative_tolerance=2e-7,absolute_tolerance=2e-9))
        rows=[dict(model="wave_pinning",original_slot="1",prefix_checks="3",bitwise_equal_checks="2",
                   maximum_absolute_difference="1e-14",original_failure="Retained native prefix failure")]
        actual=D.prefix_recovery_rows(summary,rows)
        self.assertEqual(actual[0]["original_failure"],rows[0]["original_failure"])
        self.assertIn("not stationarity",actual[0]["display_scope"])
        with self.assertRaisesRegex(ValueError,"failure-only"):
            D.prefix_recovery_rows(dict(summary,response_outcomes_used_for_selection=True),rows)
        changed=dict(summary,unchanged_native_prefix_rule=dict(relative_tolerance=2e-6,absolute_tolerance=2e-9))
        with self.assertRaisesRegex(ValueError,"changed native prefix tolerance"):
            D.prefix_recovery_rows(changed,rows)


if __name__ == "__main__":
    unittest.main()
