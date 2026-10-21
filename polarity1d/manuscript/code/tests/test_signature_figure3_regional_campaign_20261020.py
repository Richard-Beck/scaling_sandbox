"""Primary10% separation differs from the50% sensitivity on the same field."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_native_20261020 as native
import signature_figure3_regional_campaign_20261020 as campaign


class RegionalCampaignTests(unittest.TestCase):
    def test_two_peaks_can_have_merged10_percent_and_separate50_percent_contours(self):
        adapter = native.Adapter("wave_pinning", 30., 600)
        old, cue = adapter.q[0], (adapter.q[0]+.21) % 1
        x = adapter.q*30.
        p = .2+np.exp(-.5*(geometry.arc(x, old*30., 30.)/1.5)**2)
        p += .9*np.exp(-.5*(geometry.arc(x, cue*30., 30.)/1.5)**2)
        state = adapter.initial()
        state[:600] = p
        placement = dict(initial_10_percent_contour=dict(origin=old*30.),
                         basal_circumference=30., cue_center_fraction=cue)
        value = campaign.spatial(adapter, state, placement)
        self.assertTrue(value["two_region_available"])
        primary = campaign.measured_gap(adapter, state, placement, .1)
        sensitivity = campaign.measured_gap(adapter, state, placement, .5)
        self.assertFalse(primary["available"])
        self.assertTrue(sensitivity["available"])
        self.assertAlmostEqual(primary["threshold"], p.min()+.1*np.ptp(p), places=12)
        self.assertNotIn("measured_gap", primary)


if __name__ == "__main__":
    unittest.main()
