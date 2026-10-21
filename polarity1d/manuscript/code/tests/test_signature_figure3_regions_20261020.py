"""Analytic whole-region amount/identity tests, no trajectory integration."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_regions_20261020 as regions
import signature_figure3_native_20261020 as native


class CompleteRegionTests(unittest.TestCase):
    def setUp(self):
        self.C, self.n = 30., 600
        self.x = (np.arange(self.n)+.5)*self.C/self.n
        self.old, self.cue = self.x[0], self.x[240]

    def triangle(self, center, radius, height=1.):
        return height*np.maximum(1.-abs(geometry.arc(self.x, center, self.C))/radius, 0.)

    def test_exact_fractional_piecewise_linear_integral_across_seam(self):
        p = .3+self.triangle(self.old, 2.)
        total = regions.excess_integral(p, self.C, 29.173, self.C, .3)
        self.assertAlmostEqual(total, 2., places=12)
        # A cut through the peak must retain its exact linear subcell area.
        partial = regions.excess_integral(p, self.C, self.old, .773, .3)
        self.assertAlmostEqual(partial, .773-.773**2/4, places=12)

    def test_complete_two_region_excess_partitions_full_field(self):
        p = .3+self.triangle(self.old, 2.)+self.triangle(self.cue, 3.)
        result = regions.readout(p, self.C, self.old/self.C, self.cue/self.C)
        self.assertTrue(result["two_region_available"])
        self.assertAlmostEqual(result["incumbent_excess"], 2., places=12)
        self.assertAlmostEqual(result["challenger_excess"], 3., places=12)
        self.assertAlmostEqual(result["challenger_to_incumbent_ratio"], 1.5, places=12)
        self.assertAlmostEqual(result["regional_dominance"], .2, places=12)
        self.assertAlmostEqual(result["full_excess_amount"], 5., places=12)
        for fraction in (.25, .5, .75):
            cut = result["contour_sensitivity"][str(fraction)]
            self.assertAlmostEqual(cut["incumbent_excess"], 2*(1-fraction**2), places=12)
            self.assertAlmostEqual(cut["challenger_excess"], 3*(1-fraction**2), places=12)

    def test_one_relocated_region_is_never_two_region_handover(self):
        p = .3+self.triangle(self.cue, 3.)
        result = regions.readout(p, self.C, self.old/self.C, self.cue/self.C)
        self.assertEqual(result["state_class"], "single_region_cue_associated")
        self.assertFalse(result["two_region_available"])
        self.assertIsNone(result["regional_dominance"])
        self.assertAlmostEqual(result["whole_field_diagnostic"]["dominance"], 1., places=12)

    def test_amount_is_independent_of_beam_width_and_challenger_broadening_retained(self):
        p = .3+self.triangle(self.old, 1.)+self.triangle(self.cue, 4.)
        result = regions.readout(p, self.C, self.old/self.C, self.cue/self.C)
        self.assertAlmostEqual(result["challenger_excess"], 4., places=12)
        self.assertAlmostEqual(result["challenger_to_incumbent_ratio"], 4., places=12)

    def test_rotating_whole_profile_and_axes_preserves_amounts(self):
        p = .3+self.triangle(self.old, 2.)+self.triangle(self.cue, 3.)
        shift = 137
        base = regions.readout(p, self.C, self.old/self.C, self.cue/self.C)
        rotated = regions.readout(np.roll(p, shift), self.C,
            (self.old/self.C+shift/self.n) % 1, (self.cue/self.C+shift/self.n) % 1)
        self.assertAlmostEqual(base["incumbent_excess"], rotated["incumbent_excess"], places=12)
        self.assertAlmostEqual(base["challenger_excess"], rotated["challenger_excess"], places=12)
        self.assertAlmostEqual(base["whole_field_diagnostic"]["dominance"], rotated["whole_field_diagnostic"]["dominance"], places=12)

    def test_frozen_basin_linear_weights_match_fractional_integral(self):
        p = .3+self.triangle(self.old, 2.)+self.triangle(self.cue, 3.)
        for left, width in ((29.173, 7.353), (3.233, 12.456), (18.37, self.C)):
            weights = regions.linear_interval_weights(self.C, self.n, left, width)
            self.assertAlmostEqual(float((p-.3) @ weights),
                regions.excess_integral(p, self.C, left, width, .3), places=12)

    def test_complete_basin_rhs_includes_spring_physical_dilation(self):
        for key, size, cells in (("wave_pinning", 30., 120), ("spring", 2.5, 100)):
            adapter = native.Adapter(key, size, cells)
            old, cue = adapter.q[0], (adapter.q[0]+.4) % 1
            p = .2+1.2*np.exp(-(geometry.arc(adapter.q, old, 1.)/.02)**2)
            p += 1.0*np.exp(-(geometry.arc(adapter.q, cue, 1.)/.02)**2)
            state = adapter.initial()
            total = state[:cells]+state[cells:2*cells]
            state[:cells], state[cells:2*cells] = p, total-p
            proposed = regions.readout(p, size, old, cue)
            budget = regions.frozen_basin_exchange_budget(adapter, state, proposed)
            self.assertTrue(budget["available"])
            for region in budget["regions"]:
                self.assertAlmostEqual(region["exact_rate"], region["directional_finite_difference"], delta=1e-6)
                terms = sum(region[name] for name in ("chemical_rate", "diffusion_rate", "concentration_dilution_rate", "physical_weight_dilation_rate"))
                self.assertAlmostEqual(region["exact_rate"], terms, places=9)
            if key == "spring":
                self.assertGreater(max(abs(r["physical_weight_dilation_rate"]) for r in budget["regions"]), 1e-6)

    def test_high_saddle_boundary_continuous_when_minimum_switches_cells(self):
        profile = np.r_[[2., 1.5, 1.1, 1., 1., 1.1, 1.5, 2.], np.zeros(12)]
        amounts = []
        for perturbation in (-1e-8, 0., 1e-8):
            value = profile.copy()
            value[3] += perturbation
            result = regions.readout(value, 20., .5/20., 7.5/20.)
            self.assertTrue(result["two_region_available"])
            amounts.append(result["challenger_excess"])
        self.assertLess(np.ptp(amounts), 1e-5)

    def test_nonunique_active_flat_saddle_is_explicitly_unavailable(self):
        profile = np.r_[[2., 1.5, 1., 1., 1., 1., 1.5, 2.], np.zeros(12)]
        result = regions.readout(profile, 20., .5/20., 7.5/20.)
        self.assertFalse(result["two_region_available"])
        self.assertEqual(result["state_class"], "two_peaks_with_ambiguous_active_saddle")

    def test_translated_analytic_overlapping_peaks_have_continuous_partition(self):
        amounts, minima = [], []
        for shift in np.linspace(-.08, .08, 161):
            old, cue = .025+shift, 6.325+shift
            p = .2+np.exp(-.5*(geometry.arc(self.x, old, self.C)/2.)**2)
            p += .9*np.exp(-.5*(geometry.arc(self.x, cue, self.C)/2.)**2)
            result = regions.readout(p, self.C, old/self.C, cue/self.C)
            self.assertTrue(result["two_region_available"])
            self.assertAlmostEqual(result["incumbent_excess"]+result["challenger_excess"], result["full_excess_amount"], places=12)
            amounts.append(result["challenger_excess"]-result["incumbent_excess"])
            saddle = max(result["valley_descriptors"], key=lambda v: v["floor"])
            minima.append(int(np.floor(saddle["cell_coordinate"])))
        self.assertGreater(len(set(minima)), 1)
        self.assertLess(float(np.max(abs(np.diff(amounts)))), 1e-4)


if __name__ == "__main__":
    unittest.main()
