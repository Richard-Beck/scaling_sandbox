"""Finite horizons and unavailable E0 conditions survive display assembly."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
from signature_figure4_interim_render_20261020 import (
    condition_key, display_rows, interval_mark, load_inputs, main_rows,
)


class InterimDisplayTests(unittest.TestCase):
    def test_two_censored_times_cannot_become_a_finite_ratio(self):
        source = dict(tau_ratio="", tau_ratio_lower="0", tau_ratio_upper="inf",
                      tau_ratio_ci_upper="inf", tau_ratio_upper_unbounded="True")
        result = display_rows([source])[0]
        self.assertEqual(result["tau_ratio_status"], "unavailable")
        self.assertEqual(result["tau_ratio"], "")
        self.assertEqual(result["tau_ratio_ci_high"], "inf")
        self.assertEqual(result["tau_ratio_upper_unbounded"], "True")

    def test_censored_persistence_plots_the_native_lower_bound_and_ignores_fit(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        row = display_rows([dict(tau_status="right_censored", tau="", tau_lower="8000",
                                  tau_upper="inf", exponential_tau="500000")])[0]
        self.assertEqual(interval_mark(ax, 60, row, "tau_p", color="black"), "bounded")
        self.assertEqual(float(ax.lines[0].get_ydata()[0]), 8000.)
        self.assertEqual(ax.lines[0].get_marker(), "^")
        plt.close(fig)

    def test_subcadence_persistence_remains_a_resolution_limit(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        row = display_rows([dict(tau_status="unresolved_below_two_cadences", tau="",
                                  tau_lower="0", tau_upper="10")])[0]
        self.assertEqual(interval_mark(ax, 20, row, "tau_p", color="black"), "unresolved")
        self.assertEqual(float(ax.lines[0].get_ydata()[0]), 10.)
        self.assertEqual(ax.lines[0].get_marker(), "v")
        plt.close(fig)

    def test_fresh_partial_collector_cannot_export_an_interim(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "results/analysis").mkdir(parents=True)
            (root / "results/analysis/status.json").write_text(json.dumps(dict(passed=True, complete=False)))
            with self.assertRaisesRegex(ValueError, "not declared complete"):
                load_inputs(root)

    def test_spring_alias_keeps_nominal_size_separate_from_physical_size(self):
        row = dict(model="spring", C_label="3", actual_C_mean="3.48", f="0.1", view="main")
        self.assertEqual(condition_key(row), ("spring_ring", 3., .1, "main"))
        self.assertEqual(main_rows([row]), [row])


if __name__ == "__main__":
    unittest.main()
