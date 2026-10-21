"""Reject corrupted exported A5 identities before scientific plotting."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure4_protrusions_render_20261020 import validate_metric_identities, metric_series


def rows():
    return [dict(model="wave_pinning",C_label="20",f="0" if kind=="zero" else ".03",kind=kind,
                 sigma="1" if kind=="noise" else "0",experiment="E0",threshold=".5",downsample_factor="2",
                 smoothed="True",observation_cadence="80",primary="True",statistic="region_count_mean")
            for kind in ("noise","mean","zero")]


class ProtrusionIdentityTests(unittest.TestCase):
    def test_sink_kind_shadow_cannot_silently_remove_all_primary_noisy_cases(self):
        values=rows()
        self.assertEqual(len(validate_metric_identities(values)),3)
        values[0]["kind"]="survival"
        with self.assertRaisesRegex(ValueError,"condition kinds"):
            validate_metric_identities(values)

    def test_common_primary_cadence_is_independent_of_native_profile_sampling(self):
        values=rows();values[0]["observation_cadence"]="40"
        with self.assertRaisesRegex(ValueError,"common physical optics/cadence"):
            validate_metric_identities(values)

    def test_unidentified_lifetime_point_keeps_its_known_interval(self):
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        row=dict(rows()[0],estimate="",lower_bound="100",upper_bound="200",ci_lower="",ci_upper="")
        metric_series(ax,[row],{})
        segment=ax.collections[0].get_segments()[0]
        self.assertEqual(list(segment[:,1]),[100.,200.])
        self.assertFalse(any(150. in line.get_ydata() for line in ax.lines))
        plt.close(fig)


if __name__=="__main__":unittest.main()
