"""Analytic geometry and exact-RHS identities; no native trajectory integration."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_native_20261020 as native
from signature_figure3_campaign_20261020 import exchange_budget


class Figure3GeometryTests(unittest.TestCase):
    def triangular(self, circumference=30., cells=600, radius=2., origin_cell=0):
        x = (np.arange(cells)+.5)*circumference/cells
        p = .3+np.maximum(1.-abs(geometry.arc(x, x[origin_cell], circumference))/radius, 0.)
        return x, p

    def test_periodic_interpolated_edges_and_finite_volume_width(self):
        for origin_cell in (0, 11, 300, 599):
            x, p = self.triangular(origin_cell=origin_cell)
            for fraction in (.1, .5):
                c = geometry.contour(x, 30., p, fraction)
                expected = 2*(1-fraction)
                self.assertAlmostEqual(c["left"], -expected, places=12)
                self.assertAlmostEqual(c["right"], expected, places=12)
                self.assertAlmostEqual(c["weights"].sum(), 2*expected, places=12)

    def test_absolute_external_gap_grid_is_unchanged_across_sizes(self):
        for C in (30., 37.8, 47.6):
            x, p = self.triangular(C, int(round(20*C)))
            rows = geometry.compile_placements("wave_pinning", x, C, p, 3.)
            feasible = [r for r in rows if r["placement_status"] == "feasible"]
            for row in feasible:
                self.assertAlmostEqual(row["short_route"]+row["long_route"], C, places=12)
                self.assertLessEqual(row["short_route"], C/2+1e-12)
                if row["requested_role"] == "absolute_gap_grid":
                    self.assertLess(min(abs(row["edge_gap_10"]-g) for g in geometry.GAPS["wave_pinning"]), 1e-12)
                    self.assertEqual(row["classification_10"], "external")
            antipodes = [r for r in rows if r["requested_role"] == "antipode"]
            self.assertEqual(len(antipodes), 1)
            self.assertAlmostEqual(antipodes[0]["short_route"], C/2, places=12)
            self.assertAlmostEqual(antipodes[0]["long_route"], C/2, places=12)

    def test_infeasible_goryachev_negative_gap_preserved_without_model_abort(self):
        x, p = self.triangular(radius=.25)
        rows = geometry.compile_placements("goryachev", x, 30., p, 3.)
        invalid = [r for r in rows if r["placement_status"] != "feasible"]
        self.assertEqual(len(invalid), 1)
        self.assertEqual(invalid[0]["requested_gap"], -2.)
        self.assertFalse(invalid[0]["claim_eligible"])
        self.assertEqual(invalid[0]["classification_10"], "unavailable")
        self.assertTrue(any(r["requested_role"] == "antipode" for r in rows))
        self.assertFalse(any(r["requested_role"] == "internal_reference" for r in rows))

    def test_spring_uses_physical_basal_circumference_and_material_center(self):
        resting, C = 3.15, 3.19
        n = int(round(100*resting))
        x, p = self.triangular(C, n, radius=.3)
        rows = geometry.compile_placements("spring", x, C, p, .2)
        for row in rows:
            if row["placement_status"] != "feasible":
                continue
            self.assertEqual(row["basal_circumference"], C)
            self.assertAlmostEqual(row["cue_center_fraction"]*C, row["cue_center_physical"], places=12)
            self.assertAlmostEqual(row["target_angle_radians"], 2*np.pi*row["center_distance"]/C, places=12)
        self.assertEqual(native.Adapter("spring", resting).cells, 315)

    def test_low_vector_magnitude_has_coverage_gates_and_no_heading(self):
        q = (np.arange(100)+.5)/100
        z = geometry.activity_vector(np.ones(100), q)
        self.assertFalse(z["direction_available"])
        self.assertIsNone(z["heading_radians"])
        self.assertFalse(z["magnitude_eligibility_005"])
        self.assertFalse(z["magnitude_eligibility_010"])
        self.assertFalse(z["magnitude_eligibility_020"])

    def test_two_region_withdrawal_gap_crosses_periodic_seam(self):
        C, n = 30., 600
        x = (np.arange(n)+.5)*C/n
        first, second = x[0], x[240]
        p = .2+np.maximum(1.-abs(geometry.arc(x, first, C)), 0.)
        p += np.maximum(1.-abs(geometry.arc(x, second, C)), 0.)
        measured = geometry.measured_regions(x, C, p, first, second, .5)
        self.assertTrue(measured["available"])
        self.assertAlmostEqual(measured["center_distance"], 12., places=12)
        self.assertAlmostEqual(measured["measured_gap"], 11., places=12)

    def test_exact_native_rhs_regional_budget_matches_directional_difference(self):
        for key, size, cells in (("wave_pinning", 30., 120), ("spring", 2.5, 100)):
            adapter = native.Adapter(key, size, cells)
            C, q = size, adapter.q
            old, cue = q[0], (q[0]+.4) % 1
            width = .02
            p = .2+1.2*np.exp(-(geometry.arc(q, old, 1.)/width)**2)
            p += 1.0*np.exp(-(geometry.arc(q, cue, 1.)/width)**2)
            state = adapter.initial()
            total = state[:cells]+state[cells:2*cells]
            state[:cells] = p
            state[cells:2*cells] = total-p
            basal_p = .2+1.2*np.exp(-(geometry.arc(q, old, 1.)/width)**2)
            edge = geometry.contour(q*C, C, basal_p, .1)
            placement = dict(basal_circumference=C, cue_center_fraction=cue,
                initial_10_percent_contour=geometry.serial_contour(edge), support=adapter.support)
            budget = exchange_budget(adapter, state, placement)
            self.assertTrue(budget["available"])
            for region in budget["regions"]:
                self.assertAlmostEqual(region["exact_rate"], region["directional_finite_difference"], delta=1e-6)
                terms = sum(region[name] for name in ("chemical_rate", "diffusion_rate", "concentration_dilution_rate", "physical_weight_dilation_rate"))
                self.assertAlmostEqual(region["exact_rate"], terms, places=9)


if __name__ == "__main__":
    unittest.main()
