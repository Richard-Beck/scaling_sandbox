"""Keep physical imaging primary and native-factor sensitivity distinct."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_figure4_final_protrusions_render_20261020 as D
import signature_figure4_protrusions_render_20261020 as A


def rows(cadence=5., primary_factor=16, kind="noise"):
    variants = {(th, factor, smooth) for th in (.3, .5, .7) for factor in (1, 2, 4, 8) for smooth in (False, True)}
    variants.add((.5, primary_factor, True))
    return [dict(model="wave_pinning", C_label=60., f=.03, kind=kind, sigma=1. if kind == "noise" else 0.,
                 native_horizon=128000., native_profile_cadence=cadence, threshold=th,
                 downsample_factor=factor, smoothed=smooth, observation_cadence=cadence*factor,
                 primary=(th == .5 and smooth and factor == primary_factor), statistic=statistic,
                 independent_seeds=12 if kind == "noise" else 0,
                 bootstrap_resamples=2000 if kind == "noise" else 0)
            for th, factor, smooth in variants for statistic in A.NAMES]


class FinalA5DisplayTests(unittest.TestCase):
    def test_fine_native_primary_is_a_distinct_25th_variant(self):
        data = rows()
        D.validate_metrics(data, {D.key(data[0])})
        self.assertEqual(len({D.variant(row) for row in data}), 25)
        with self.assertRaisesRegex(ValueError, "one complete physical-cadence primary"):
            D.validate_metrics([row for row in data if row["downsample_factor"] != 16], {D.key(data[0])})

    def test_fixed_native_factor_cannot_replace_common_physical_cadence(self):
        data = rows(primary_factor=2)
        with self.assertRaisesRegex(ValueError, "joined by physical cadence"):
            D.validate_metrics(data, {D.key(data[0])})
        data = rows(); data[0]["observation_cadence"] = 80.
        if data[0]["downsample_factor"] == 16:
            data[0]["observation_cadence"] = 10.
        with self.assertRaisesRegex(ValueError, "different physical observation cadence"):
            D.validate_metrics(data, {D.key(data[0])})

    def test_numerical_companions_and_deterministic_controls_do_not_become_seeds(self):
        data = rows(kind="mean")
        D.validate_metrics(data, {D.key(data[0])})
        data[0]["independent_seeds"] = 1
        with self.assertRaisesRegex(ValueError, "zero independent seeds"):
            D.validate_metrics(data, {D.key(data[0])})

    def test_declared_horizon_and_all_five_size_styles_are_preserved(self):
        data = rows(); identity = D.key(data[0])
        data[0]["native_horizon"] = 32000.
        with self.assertRaisesRegex(ValueError, "outside the exact canonical inventory"):
            D.validate_metrics(data, {identity})
        self.assertEqual([A.size_line_style(i) for i in range(3)], ["-", "--", ":"])
        self.assertEqual(len({repr(A.size_line_style(i)) for i in range(5)}), 5)
        self.assertTrue(A.annotate_forced(dict(display_forced="True"), {}))
        self.assertFalse(A.annotate_forced(dict(display_forced=None), {}))

    def test_event_count_series_never_connect_different_observation_windows(self):
        data=[dict(kind="noise",native_horizon=32000),dict(kind="noise",native_horizon=128000)]
        with patch.object(A,"metric_series") as plot:
            D.metric_series_by_horizon(None,data,count=True)
        self.assertEqual(plot.call_count,2)
        self.assertEqual([call.kwargs["point_marker"] for call in plot.call_args_list],["o","s"])
        self.assertEqual([call.kwargs["label_suffix"] for call in plot.call_args_list],[", H=32000",", H=128000"])
        self.assertTrue(all(len(call.args[1])==1 for call in plot.call_args_list))


if __name__ == "__main__":
    unittest.main()
