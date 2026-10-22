"""Condition census preserves dependencies rather than pooling paths as seeds."""
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_figure4_design_render_20261020 as D

class DesignCensusTests(unittest.TestCase):
    def test_sigma_companions_share_lineage_and_deterministic_controls_add_no_seeds(self):
        case=dict(atlas_context='common',model='wave_pinning',C_label=20.,f=.05,kind='noise',sigma=1.,native_horizon=32000.,latent_Gaussian_cohort_id='same_cohort',independent_stochastic_seed_count=1,scalar_cadence=5.,profile_cadence=5.,origin='native')
        cases=[case,dict(case,atlas_context='E1_extra',sigma=.5),dict(case,kind='mean',sigma=0.,latent_Gaussian_cohort_id=None,independent_stochastic_seed_count=0)]
        rows=D.condition_rows(cases)
        self.assertEqual(sum(row['native_paths'] for row in rows),3)
        self.assertEqual([row['within_condition_stochastic_seed_blocks'] for row in rows if row['kind']=='mean'],[0])
        self.assertTrue(all('dependent across conditions' in row['independence_scope'] for row in rows))
        cases[-1]['independent_stochastic_seed_count']=1
        with self.assertRaisesRegex(ValueError,'Deterministic controls'):
            D.condition_rows(cases)

if __name__=='__main__':unittest.main()
