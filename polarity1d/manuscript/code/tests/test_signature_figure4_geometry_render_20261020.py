"""Geometric displays retain bounds and conflicting raw evidence."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mechanisms"))
import signature_figure4_geometry_render_20261020 as D
import signature_figure4_interim_render_20261020 as P


class GeometricDisplayTests(unittest.TestCase):
    def test_finite_resolution_interval_does_not_create_a_midpoint(self):
        fig, ax = plt.subplots()
        D.ci_mark(ax, .03, dict(excess_persistence="", excess_persistence_ci_lower="0",
            excess_persistence_ci_upper="2"), "excess_persistence", "excess_persistence", color="black", qualified=False)
        self.assertEqual(ax.lines[0].get_marker(), "x")
        self.assertEqual(list(ax.collections[0].get_segments()[0][:, 1]), [0., 2.])
        plt.close(fig)

    def test_unbounded_CI_keeps_its_endpoint_and_unqualified_point_open(self):
        fig, ax = plt.subplots()
        D.ci_mark(ax, .1, dict(excess_persistence="3", excess_persistence_ci_lower="1",
            excess_persistence_ci_upper="", excess_persistence_ci_upper_unbounded="True"),
            "excess_persistence", "excess_persistence", color="black", qualified=False)
        self.assertEqual(list(ax.lines[0].get_ydata()), [3.])
        self.assertEqual(ax.lines[0].get_markerfacecolor(), "none")
        self.assertEqual(ax.lines[1].get_marker(), "^")
        self.assertEqual(list(ax.lines[1].get_ydata()), [1.])
        self.assertEqual(len(ax.collections), 0)
        plt.close(fig)

    def test_raw_negative_mobility_is_retained_without_a_geometric_inference(self):
        fig, ax = plt.subplots()
        D.ci_mark(ax, .1, dict(Dx_large_over_small="-2", Dx_ratio_ci_lower="-3", Dx_ratio_ci_upper=".1"),
            "Dx_large_over_small", "Dx_ratio", color="black", qualified=False)
        self.assertEqual(list(ax.lines[0].get_ydata()), [-2.])
        self.assertEqual(list(ax.collections[0].get_segments()[0][:, 1]), [-3., .1])
        plt.close(fig)

    def test_missing_compatibility_qualification_cannot_silently_become_equivalence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); folder = root / "results/geometric_qualifiers"; folder.mkdir(parents=True)
            status = folder / "status.json"
            status.write_text(json.dumps(dict(passed=True, complete=True, compatibility_is_not_equivalence=False)))
            with self.assertRaisesRegex(ValueError, "scientific qualifications"):
                D.load_geometry(root, P.sha(status))


if __name__ == "__main__":
    unittest.main()
