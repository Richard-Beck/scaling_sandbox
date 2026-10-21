"""Prevent silent horizon pooling and gap interpolation in later F4 displays."""
from pathlib import Path
import sys
import unittest

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure4_stage_render_20261020 import select_observation_curves, annotation, presentation_bounds,sigma_metric_confidence
from signature_figure4_interim_render_20261020 import forced, curve, interval_mark
from signature_figure4_interim_render_20261020 import display_rows,coincident_bound_markers


def row(**changes):
    return dict(model="wave_pinning",C_label="20",f=".03",view="main",native_horizon="32000",
                selected_observation_arm="common_E2_augmented_E0",analysis_arm="common",**changes)


class StageDisplayTests(unittest.TestCase):
    def test_native_fraction_confidence_fields_remain_visible_in_sigma_projection(self):
        native=dict(single_front_ci_lower='.11',single_front_ci_upper='.91',
                    multiple_fronts_ci_lower='.12',multiple_fronts_ci_upper='.92',
                    low_contrast_ci_lower='.13',low_contrast_ci_upper='.93')
        self.assertEqual(sigma_metric_confidence(native,'single_fraction'),(.11,.91))
        self.assertEqual(sigma_metric_confidence(native,'multiple_fraction'),(.12,.92))
        self.assertEqual(sigma_metric_confidence(native,'low_contrast_fraction'),(.13,.93))
        self.assertTrue(np.isnan(sigma_metric_confidence({},'multiple_fraction')[0]))

    def test_long_choice_never_falls_back_to_an_available_common_curve(self):
        metric=row();metric.update(selected_observation_arm="preselected_E4",analysis_arm="longest",native_horizon="128000")
        with self.assertRaisesRegex(ValueError,"no exact curve"):
            select_observation_curves([metric],{"E2":[row(lag="5",C=".9")]})
        long=row();long.update(native_horizon="128000",lag="5",C=".9")
        selected=select_observation_curves([metric],{"E4":[long],"E2":[row(lag="5",C=".9")]})
        self.assertEqual(len(selected),1)
        self.assertEqual(selected[0]["analysis_arm"],"longest")

    def test_wrong_horizon_cannot_share_the_same_condition_identity(self):
        with self.assertRaisesRegex(ValueError,"different native horizon"):
            select_observation_curves([row()],{"E2":[dict(row(),native_horizon="128000")]})

    def test_unknown_combined_flag_does_not_become_a_negative_mean_flag(self):
        annotated=annotation(dict(forced="False",forced_combined=None))
        self.assertIsNone(forced(annotated))
        self.assertEqual(annotated["forced"],"False")

    def test_unavailable_lags_break_the_displayed_curve(self):
        fig,ax=plt.subplots()
        curve(ax,[dict(lag="5",C=".9"),dict(lag="10",C=""),dict(lag="15",C=".3")],
              "lag",("C",),color="black",label=None)
        self.assertTrue(np.isnan(ax.lines[0].get_ydata()[1]))
        plt.close(fig)

    def test_inconclusive_estimand_retains_a_known_finite_interval(self):
        fig,ax=plt.subplots()
        state=interval_mark(ax,.1,dict(fhalf_status="inconclusive_interval_overlap",fhalf_lower=".03",fhalf_upper="inf"),
                            "fhalf",color="black")
        self.assertEqual(state,"bounded")
        self.assertEqual(list(ax.lines[0].get_ydata()),[.03])
        plt.close(fig)

    def test_sub_two_cadence_resolution_bound_has_no_fabricated_point(self):
        fig,ax=plt.subplots()
        state=interval_mark(ax,.03,dict(tau_p_status="unresolved_below_two_cadences",
            tau_p="",tau_p_lower="0",tau_p_upper="10",persistence_bound_eligible="True"),
            "tau_p",color="black")
        self.assertEqual(state,"unresolved")
        self.assertEqual(ax.lines[0].get_marker(),"v")
        self.assertEqual(list(ax.lines[0].get_ydata()),[10.])
        plt.close(fig)

    def test_resolved_empirical_point_retains_unbounded_bootstrap_uncertainty(self):
        # Accepted common Gory C60/f.3 and Holmes C20/f.01 exemplars:
        # an empirical crossing exists, but its bootstrap upper is unbounded.
        rows=display_rows([dict(tau='1256.849898429333',tau_status='resolved',tau_ci_lower='0.0',
            tau_ci_upper='',tau_ci_upper_unbounded='True'),dict(tau='6814.045540395658',tau_status='resolved',
            tau_ci_lower='3751.7559472943744',tau_ci_upper='',tau_ci_upper_unbounded='True')])
        fig,ax=plt.subplots()
        for x,row in zip((.3,.01),rows):
            self.assertEqual(interval_mark(ax,x,row,'tau_p',color='black'),'resolved')
        self.assertEqual([t.get_text() for t in ax.texts],['CI[0,∞)','CI∞'])
        self.assertEqual([float(line.get_ydata()[0]) for line in ax.lines],[1256.849898429333,6814.045540395658])
        self.assertEqual(ax.collections[0].get_segments()[0][0,1],3751.7559472943744)
        self.assertEqual(rows[0]['display_tau_p_CI_lower_endpoint'],'0.0')
        self.assertTrue(all(row['display_tau_p_finite_point_upper_unbounded_CI'] for row in rows))
        plt.close(fig)

    def test_coincident_size_bounds_keep_every_original_coordinate_and_identity(self):
        rows=display_rows([dict(model='wave_pinning',C_label=str(size),f='.001',tau_status='right_censored',
            tau_lower='8000',tau_upper='inf',forced='False') for size in (20,60)])
        original=[dict(row) for row in rows];fig,ax=plt.subplots();entries=[]
        for row in rows:
            first=len(ax.lines);state=interval_mark(ax,.001,row,'tau_p',color='black')
            entries.append((row,state,ax.lines[first]))
        coincident_bound_markers(ax,entries,'tau_p')
        self.assertEqual([float(line.get_xdata()[0]) for line in ax.lines],[.001,.001])
        self.assertEqual([float(line.get_ydata()[0]) for line in ax.lines],[8000.,8000.])
        self.assertEqual([line.get_markersize() for line in ax.lines],[4.,6.])
        self.assertEqual(ax.texts[0].get_text(),'×2')
        for row,source in zip(rows,original):
            self.assertTrue(all(row[key]==value for key,value in source.items()))
            self.assertEqual(row['display_tau_p_coincident_bound_count'],2)
            self.assertEqual(row['display_tau_p_coincident_bound_sizes'],'[20.0, 60.0]')
            self.assertTrue(row['display_tau_p_scientific_coordinate_unchanged'])
        plt.close(fig)

    def test_presentation_bounds_distinguish_actual_clipping_from_complete_ticks(self):
        fig,ax=plt.subplots(figsize=(5,3));fig.subplots_adjust(left=.2,right=.8,bottom=.2,top=.8)
        ax.set_xlim(0,8000);ax.set_xticks([0,8000]);self.assertTrue(presentation_bounds(fig)["passed"])
        ax.set_xticks([0,8000,100000]);ax.set_xlim(0,8000)
        self.assertTrue(presentation_bounds(fig)["passed"])
        fig.text(1.01,.5,"clipped text",ha="left")
        check=presentation_bounds(fig);self.assertFalse(check["passed"])
        self.assertTrue(any(row["text"]=="clipped text" and not row["within_page"] for row in check["items"]))
        plt.close(fig)

    def test_footer_overlap_is_rejected_even_when_every_artist_is_on_page(self):
        fig,ax=plt.subplots(figsize=(5,3));fig.subplots_adjust(bottom=.2,top=.8)
        ax.set_xlabel("Lag (native time)")
        footer=fig.text(.2,.08,"A long footer\nSecond line\nThird line",va="bottom",fontsize=10)
        footer.set_gid("document_footer")
        check=presentation_bounds(fig)
        self.assertTrue(check["within_page"]);self.assertFalse(check["header_footer_clear_of_axes"])
        self.assertTrue(check["collisions"])
        footer.set_y(.01);fig.subplots_adjust(bottom=.5)
        self.assertTrue(presentation_bounds(fig)["passed"])
        plt.close(fig)

    def test_document_body_and_axis_annotations_are_included_in_collision_scope(self):
        fig,ax=plt.subplots(figsize=(5,3));ax.set_visible(False)
        header=fig.text(.1,.9,"Document header",va="top");header.set_gid("document_header")
        body=fig.text(.1,.9,"Scientific decision text",va="top")
        check=presentation_bounds(fig)
        self.assertTrue(check["within_page"])
        self.assertFalse(check["header_footer_clear_of_axes"])
        body.set_y(.7);self.assertTrue(presentation_bounds(fig)["passed"])
        ax.set_visible(True);ax.text(1.5,.5,"clipped annotation",transform=ax.transAxes)
        self.assertFalse(presentation_bounds(fig)["within_page"])
        plt.close(fig)


if __name__=="__main__":unittest.main()
