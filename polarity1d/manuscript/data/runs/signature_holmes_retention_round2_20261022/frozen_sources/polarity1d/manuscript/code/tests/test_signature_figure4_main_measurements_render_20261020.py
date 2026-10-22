"""Compact main panels retain bounds, exact rows and source qualifications."""
import math
from pathlib import Path
import sys
import tempfile
import unittest

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure4_main_measurements_render_20261020 as M


class CompactMeasurementTests(unittest.TestCase):
    def fixture(self):
        ratio=dict(arm="common",view="main",model="wave_pinning",f=".03",small_C_label="20",large_C_label="60",
            tau_ratio="",tau_ratio_lower="2.5",tau_ratio_upper="inf",tau_ratio_upper_unbounded="True",
            tau_ratio_ci_lower="2",tau_ratio_ci_upper="inf",tau_ratio_ci_upper_unbounded="True",
            forced_combined="True",native_horizon_small="32000",native_horizon_large="32000",actual_geometry_ratio="9")
        panel=dict(id="B",measurement="persistence_ratio",arm="common",conditions=[dict(model="wave_pinning",f=.03)])
        return dict(**{"size_ratios.csv":[ratio]}),panel

    def test_bounds_are_not_fabricated_points_or_clipped_CI_endpoints(self):
        tables,panel=self.fixture();row=M.selected_rows(tables,panel)[0]
        fig,ax=plt.subplots();M.measured_mark(ax,0,row,"persistence_ratio","blue","o")
        self.assertEqual(row["display_point_status"],"empirical lower bound")
        self.assertEqual(float(ax.lines[0].get_ydata()[0]),2.5)
        self.assertEqual(ax.lines[0].get_marker(),"^")
        self.assertEqual(row["tau_ratio"],"")
        self.assertEqual(row['display_CI_annotation'],'CI∞')
        self.assertEqual(row['display_CI_lower_endpoint_retained'],2.)
        self.assertEqual(float(ax.lines[1].get_ydata()[0]),2.)
        plt.close(fig)
        row.update(tau_ratio="3",tau_ratio_ci_lower="1.7")
        fig,ax=plt.subplots();M.measured_mark(ax,0,row,"persistence_ratio","blue","o")
        self.assertEqual(row["display_point_status"],"native point")
        self.assertEqual(row["display_CI_annotation"],"CI∞")
        self.assertEqual(float(ax.lines[0].get_ydata()[0]),3.)
        self.assertEqual(float(ax.collections[0].get_segments()[0][0,1]),1.7)
        plt.close(fig)

    def test_exact_quotient_identity_and_original_fields_survive_projection(self):
        tables,panel=self.fixture();panel["measurement"]="geometry_normalized_persistence"
        row=M.selected_rows(tables,panel)[0];row.update(excess_persistence="",excess_persistence_ci_lower=".2",excess_persistence_ci_upper="inf")
        point,bounds,ci=M.endpoints(row,panel["measurement"])
        self.assertTrue(math.isnan(point));self.assertEqual(bounds,(2.5/9,math.inf))
        self.assertEqual(row["tau_ratio_lower"],"2.5")
        self.assertEqual(row["display_operationally_forced"],"True")
        self.assertIn("actual_geometry_ratio",row["display_empirical_bound_identity"])

    def test_native_zero_is_explicit_and_never_replaced_by_a_positive_log_point(self):
        row=dict(estimate="0",ci_lower="0",ci_upper=".8",ci_upper_unbounded="False")
        fig,ax=plt.subplots();M.measured_mark(ax,0,row,"fingerprint_region_count","blue","o")
        self.assertEqual(row["estimate"],"0")
        self.assertIn("native zero",row["display_point_status"])
        self.assertEqual(len(ax.lines),0)
        self.assertEqual(ax.texts[0].get_text(),"0\nCI[0,0.8]")
        plt.close(fig)

    def test_missing_or_duplicate_declared_conditions_refuse_and_source_page_fits(self):
        tables,panel=self.fixture();rows=M.selected_rows(tables,panel)
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/"compact.pdf";bounds=M.compact_page(rows,panel,output)
            self.assertTrue(bounds["passed"]);self.assertTrue(output.is_file())
        panel["conditions"][0]["f"]=.1
        with self.assertRaisesRegex(ValueError,"exactly one"):
            M.selected_rows(tables,panel)
        panel["conditions"][0]["f"]=.03;tables["size_ratios.csv"]*=2
        with self.assertRaisesRegex(ValueError,"exactly one"):
            M.selected_rows(tables,panel)

    def test_full_overview_retains_every_model_input_and_bound_without_subset_selection(self):
        tables,_=self.fixture();ratios=[];conditions=[]
        for model in M.P.MODELS:
            for strength in (.001,.03):
                row=dict(tables["size_ratios.csv"][0],model=model,f=strength)
                ratios.append(row)
                conditions.append(dict(arm="common",view="main",model=model,f=strength,C_label=60,regime="multiple_fronts"))
        tables={"size_ratios.csv":ratios,"condition_metrics.csv":conditions}
        panel=dict(measurement="persistence_ratio_overview",arm="common")
        rows=M.selected_rows(tables,panel)
        with tempfile.TemporaryDirectory() as temp:
            bounds=M.overview_page(rows,panel,Path(temp)/"overview.pdf")
            self.assertTrue(bounds["passed"])
        self.assertEqual(len(rows),8)
        self.assertTrue(all(row["display_bound_or_CI_symbol"]=="↑" for row in rows))
        self.assertTrue(all(row["tau_ratio"]=="" for row in rows))
        panel["conditions"]=[dict(model="wave_pinning",f=.03)]
        with self.assertRaisesRegex(ValueError,"favorable subset"):
            M.selected_rows(tables,panel)

    def test_all_views_retains_each_native_bound_and_confidence_endpoint(self):
        tables,panel=self.fixture();original=tables['size_ratios.csv'][0]
        tables['size_ratios.csv']=[dict(original,view=view) for view in M.D.VIEWS]
        panel.update(measurement='persistence_ratio_all_views')
        rows=M.selected_rows(tables,panel)
        self.assertEqual([row['display_view'] for row in rows],list(M.D.VIEWS))
        with tempfile.TemporaryDirectory() as temp:
            self.assertTrue(M.all_views_page(rows,panel,Path(temp)/'views.pdf')['passed'])
        self.assertTrue(all(row['tau_ratio']=='' and row['display_CI_lower_endpoint_retained']==2. for row in rows))
        self.assertTrue(all(row['display_point_status']=='empirical lower bound' for row in rows))
        panel['views']=['main']
        with self.assertRaisesRegex(ValueError,'seven predeclared'):
            M.selected_rows(tables,panel)


if __name__ == "__main__":unittest.main()
