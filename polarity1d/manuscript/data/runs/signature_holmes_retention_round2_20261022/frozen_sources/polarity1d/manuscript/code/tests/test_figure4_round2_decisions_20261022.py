"""Scientific boundary cases for unchanged S0 truth under new estimands."""
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import figure4_round2_decisions_20261022 as S


def ratio(f, lo=2., hi=4., forced=False, point=3., multi=.8):
    return S.mapped_ratio(dict(f=f, rho_k=point, rho_k_lower=point, rho_k_upper=point,
                              rho_k_ci_lower=lo, rho_k_ci_upper=hi, evaluable=True,
                              coverage_small=.9, coverage_large=.9, forced_combined=forced,
                              large_multiple_front_fraction=multi))


class FrozenCriteria(unittest.TestCase):
    small = dict(fhalf_ci_lower=.001, fhalf_ci_upper=.002)
    large = dict(fhalf_ci_lower=.003, fhalf_ci_upper=.004)

    def h1(self, rows, inclusive=False, large=None):
        return S.F.h1_model(rows, self.small, large or self.large, inclusive=inclusive)["status"]

    def test_gain_at_weak_input_does_not_override_uncertain_later_input(self):
        self.assertEqual(self.h1([ratio(.003), ratio(.05, .7, 1.2)]), "inconclusive")

    def test_forced_reversal_is_excluded_only_from_primary(self):
        rows = [ratio(.003), ratio(.05, .2, .8, forced=True)]
        self.assertEqual(self.h1(rows), "pass")
        self.assertEqual(self.h1(rows, inclusive=True), "fail")

    def test_unknown_forcing_cannot_silently_certify_primary(self):
        self.assertEqual(self.h1([ratio(.003), ratio(.01, forced=None)]), "inconclusive")

    def test_bound_can_certify_gain_without_inventing_point(self):
        row = ratio(.003, 1.2, math.inf, point=None)
        row["tau_ratio_lower"] = 1.5
        row["tau_ratio_upper"] = math.inf
        self.assertEqual(self.h1([row]), "pass")
        self.assertIsNone(row["tau_ratio"])

    def test_overlapping_tolerance_prevents_h1_pass(self):
        self.assertEqual(self.h1([ratio(.003)], large=dict(fhalf_ci_lower=.0015, fhalf_ci_upper=.003)), "inconclusive")

    def test_h2_peak_drop_cannot_bridge_forced_gap(self):
        rows = [ratio(.003, 4., 12., point=8.), ratio(.01, forced=True, point=4.), ratio(.05, .8, 1.3, point=1.)]
        self.assertEqual(S.F.h2_model(rows)["status"], "inconclusive")
        self.assertEqual(S.F.h2_model(rows, inclusive=True)["status"], "pass")


if __name__ == "__main__":
    unittest.main()
