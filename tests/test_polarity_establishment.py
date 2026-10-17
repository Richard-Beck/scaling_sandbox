"""Scientific selection, persistent-arrival and profile-clustering checks."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'polarity1d/four_assay'))
import polarity_establishment as e


class EstablishmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = e.make_manifest()

    def trial(self, duration, strength, verified=True, contrast=.2):
        return dict(duration=duration, cue_aligned_strength=strength, verified=verified,
                    final_orientation=-strength, final_contrast=contrast, polarity_arrival=2.)

    def test_legi_basal_equilibrium_and_patch_withdrawal(self):
        case = next(c for c in self.manifest['cases'] if c['model'] == 'legi')
        state = e.resting_state(case)
        model = e.r.model_for('legi')
        x = np.linspace(0, case['length'], case['cells'])
        basal = e.withdrawal(case)(0, x, case['length'])
        np.testing.assert_allclose(basal, .1)
        np.testing.assert_allclose(state[2], .5)
        np.testing.assert_allclose(model.reaction(0, x, state, basal), 0, atol=1e-15)
        pulse = e.Patch(case['amplitude'], case['width'], stop=2, background=.1)
        self.assertGreater(pulse(1, x, case['length'])[0], .1)
        np.testing.assert_allclose(pulse(2, x, case['length']), basal)
        # An asymmetric response continues to relax after withdrawal.
        state[2, 0] = .8
        self.assertLess(model.reaction(0, x, state, basal)[2, 0], 0)
        tr = e.r.run('legi', case['length'], state, e.withdrawal(case), 200., n=case['cells'])
        np.testing.assert_allclose(tr.state[2], .5, atol=1e-7)

    def test_arrival_requires_last_exit_and_confirmation(self):
        time = np.array([0., 1., 2., 3., 1003.])
        values = np.array([np.nan, .5, .7, .5, .5])
        self.assertEqual(e.persistent_arrival(time, values, .5, .005), 3.)
        self.assertIsNone(e.persistent_arrival(time[:-1], values[:-1], .5, .005))

    def test_strongest_then_shortest_not_first_polarized(self):
        trials = [self.trial(0, 0), self.trial(1, .1), self.trial(2, .499),
                  self.trial(3, .5), self.trial(4, .3), self.trial(5, .498),
                  self.trial(6, .9, verified=False)]
        passing, best = e.qualifying(trials, self.manifest)
        self.assertEqual(best, .5)
        self.assertEqual([t['duration'] for t in passing], [2, 3, 5])
        # Both disconnected intervals into the strongest band are refined.
        np.testing.assert_allclose(e.threshold_candidates(trials, self.manifest),
                                   [np.sqrt(2), np.sqrt(20)])

    def test_unpolarized_and_insufficient_contrast_not_selected(self):
        self.assertEqual(e.qualifying([self.trial(1, .02)], self.manifest)[0], [])
        self.assertEqual(e.qualifying([self.trial(1, .5, contrast=.001)], self.manifest)[0], [])
        passing, best = e.qualifying([self.trial(1, .9, contrast=.001), self.trial(2, .5)], self.manifest)
        self.assertEqual(best, .5)
        self.assertEqual(passing[0]['duration'], 2)

    def test_profiles_not_scalar_score_define_endpoint_count(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            trials = []
            for i, vector in enumerate([[0., 0.], [1., 0.], [1.005, 0.], [0., 1.]]):
                np.savez(out / f'{i}.npz', endpoint_active=np.array(vector))
                trials.append(dict(duration=float(i), verified=True, trace=f'{i}.npz'))
            count, assignments = e.endpoint_count(out, dict(fixed_scales={'active': 1.}), trials, .01)
            self.assertEqual(count, 3)
            self.assertEqual(assignments[1]['endpoint_cluster'], assignments[2]['endpoint_cluster'])

    def test_no_cue_needed_is_zero_duration(self):
        passed, best = e.qualifying([self.trial(0, .5), self.trial(1, .5)], self.manifest)
        self.assertEqual(passed[0]['duration'], 0)
        self.assertEqual(e.threshold_candidates([self.trial(0, .5), self.trial(1, .5)], self.manifest), [])

    def test_scheduler_arrays_compact_full_grid_and_sparse_resume(self):
        self.assertEqual(list(e.array_batches(list(range(2156)), 100001)),
                         ['0-899', '900-1799', '1800-2155'])
        for spec in e.array_batches(list(range(0, 100000, 3)), 100001):
            self.assertLessEqual(len(spec), 3000)

    def test_complete_link_avoids_chaining_continuous_family(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            trials = []
            for i, value in enumerate([0., .009, .018]):
                np.savez(out / f'{i}.npz', endpoint_active=np.array([value]))
                trials.append(dict(duration=float(i), verified=True, trace=f'{i}.npz'))
            self.assertEqual(e.endpoint_count(out, dict(fixed_scales={'active': 1.}), trials, .01)[0], 2)

    def test_standardized_reference_and_distinct_three_conditions(self):
        cases = self.manifest['cases']
        self.assertEqual(len(cases), 3234)
        self.assertEqual({c['reference_length'] for c in cases if c['model'] != 'spring'}, {15.})
        self.assertEqual({c['reference_length'] for c in cases if c['model'] == 'spring'}, {1.})
        for c in cases:
            if c['width_mode'] == 'physical':
                self.assertAlmostEqual(c['width'], c['reference_length'] * c['width_fraction_at_reference'])
            else:
                self.assertAlmostEqual(c['width'], .2 * c['length'])
        for model in e.r.MODELS:
            self.assertEqual({c['condition'] for c in cases if c['model'] == model},
                             {'fixed_10pct', 'fixed_20pct', 'relative_20pct'})

    def test_relative_width_scales_and_fixed_width_does_not(self):
        for length in [1., 15., 30.]:
            x = np.linspace(0, length, 301)
            left = e.Patch(2., 3., width_mode='relative')(0, x, length)
            right = e.Patch(2., 3., side='right', width_mode='relative')(0, length-x, length)
            np.testing.assert_allclose(left, right, atol=1e-13)
            np.testing.assert_allclose(left, e.Patch(2., .2*length)(0, x, length), atol=1e-13)
        # Spring's moving domain must update the relative support at each call.
        cue = e.Patch(2., .2, width_mode='relative')
        self.assertEqual(cue(0, np.array([.3]), 1.)[0], 0.)
        self.assertGreater(cue(0, np.array([.3]), 2.)[0], 0.)
        self.assertEqual(e.Patch(2., .2)(0, np.array([.3]), 2.)[0], 0.)


if __name__ == '__main__':
    unittest.main()
