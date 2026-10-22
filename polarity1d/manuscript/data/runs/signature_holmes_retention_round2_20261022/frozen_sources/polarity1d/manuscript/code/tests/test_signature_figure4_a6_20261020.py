from pathlib import Path
import sys
import unittest
import tempfile
import json
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_figure4_a6_20261020 as A
import signature_figure4_a6_audit_20261020 as Q


def points(models=A.MODELS):
    return pd.DataFrame([dict(model=m,f=f,feature_ratio=x,log_persistence_ratio=-np.log(x))
        for m in models for f,x in zip([.01,.03,.1],[1.,2.,3.])])


class AssociationTests(unittest.TestCase):
    def test_independent_average_rank_audit_handles_ties_and_unavailable_samples(self):
        np.testing.assert_array_equal(Q.average_ranks([5.,1.,1.,3.]),[4.,1.5,1.5,3.])
        result=Q.independent_rank([0,0,1,2],[2,2,1,0],[.01,.02,.03,.04])
        self.assertAlmostEqual(result['spearman'],-1.)
        self.assertEqual(Q.independent_rank([1,1,1],[1,2,3],[.01,.02,.03])['reason'],'constant_predictor')

    def test_independent_source_and_fold_audit_rejects_fabricated_primary_support(self):
        features,response=self.expanded_fixture();p=A.condition_pairs(features,response,experiment='common')
        tables=A.association_tables(p)
        proof=Q.audit_source_points(p,features,response,'common',1.)
        self.assertTrue(proof['all_source_observation_views_retained'])
        Q.audit_ranks(p,tables['cross_model'],tables['within_model'],tables['LOMO_folds'],tables['H2_primary_features'])
        corrupted=tables['H2_primary_features'].copy();corrupted.loc[corrupted.primary_view,'A6_primary_predictor_support']=True
        with self.assertRaisesRegex(ValueError,'Primary H2 criterion'):
            Q.audit_ranks(p,tables['cross_model'],tables['within_model'],tables['LOMO_folds'],corrupted)

    @staticmethod
    def expanded_fixture():
        features=[]
        for size in [30.,60.]:
            for statistic in A.PRIMARY:
                features.append(dict(model='wave_pinning',f=.01,C_label=size,native_horizon=32000.,kind='noise',experiment='common',sigma=1.,
                    threshold=.5,smoothed=True,primary=True,observation_cadence=80.,native_profile_cadence=5. if size==60. else 40.,
                    statistic=statistic,estimate=size/30.,ci_lower=.8,ci_upper=2.2,independent_seeds=8))
        response=dict(model='wave_pinning',f=.01,view='main',small_C_label=30.,large_C_label=60.,analysis_arm='common',sigma=1.,
            native_horizon_small=32000.,native_horizon_large=32000.,evaluable=True,forced_combined=False,
            tau_status_small='resolved',tau_status_large='resolved',tau_ratio=2.,Dtheta_mobility_proxy=2.,mobility_proxy_valid_bootstrap_fraction=1.)
        return pd.DataFrame(features),pd.DataFrame([response])

    def test_four_negative_transfers_pass_the_predeclared_sign(self):
        folds=A.leave_one_model_out(points(),primary=True)
        self.assertEqual(len(folds),4)
        self.assertTrue(all(f['sign_correct'] for f in folds))
        self.assertTrue(A.primary_criterion(-1.,folds)['A6_primary_predictor_support'])

    def test_missing_models_remain_in_four_fold_denominator(self):
        folds=A.leave_one_model_out(points(A.MODELS[:2]),primary=True)
        self.assertEqual(sum(f['status']=='unavailable' for f in folds),2)
        result=A.primary_criterion(-1.,folds)
        self.assertEqual(result['planned_folds'],4)
        self.assertEqual(result['successful_folds'],2)
        self.assertFalse(result['A6_primary_predictor_support'])

    def test_positive_training_and_heldout_signs_do_not_rescue_primary_negative_rule(self):
        p=points();p.log_persistence_ratio*=-1
        self.assertFalse(any(f['sign_correct'] for f in A.leave_one_model_out(p,primary=True)))
        self.assertTrue(all(f['sign_correct'] for f in A.leave_one_model_out(p,primary=False)))

    def test_constant_features_and_insufficient_distinct_input_levels_are_unavailable(self):
        self.assertEqual(A.rank_association([1,1,1],[1,2,3],[.01,.03,.1])['reason'],'constant_predictor')
        self.assertEqual(A.rank_association([1,2,3],[1,2,3],[.01,.01,.03])['reason'],'fewer_than_three_matched_f')

    def test_zero_denominator_is_unavailable_but_zero_numerator_is_measured(self):
        self.assertTrue(np.isnan(A.quotient(1.,0.)[0]))
        self.assertEqual(A.quotient(0.,1.),(0.,'observed'))

    def test_censored_direct_tau_is_not_replaced_by_mobility_proxy(self):
        out=A.outcomes(dict(tau_status_small='right_censored',tau_status_large='resolved',evaluable=True,
            tau_ratio=99.,Dtheta_mobility_proxy=2.,mobility_proxy_valid_bootstrap_fraction=1.))
        self.assertTrue(np.isnan(out['direct_tau'][0]))
        self.assertAlmostEqual(out['Dtheta_proxy'][0],np.log(2.))

    def test_bound_eligible_below_two_cadences_never_creates_a_direct_tau_point(self):
        out=A.outcomes(dict(tau_status_small='unresolved_below_two_cadences',tau_status_large='resolved',evaluable=True,
            tau_ratio=None,tau_ratio_lower=2.,tau_ratio_upper=np.inf,Dtheta_mobility_proxy=3.,mobility_proxy_valid_bootstrap_fraction=1.))
        self.assertTrue(np.isnan(out['direct_tau'][0]))
        self.assertAlmostEqual(out['Dtheta_proxy'][0],np.log(3.))

    def test_native_unit_change_cancels_only_in_within_model_size_quotient(self):
        rows=[]
        for C,x in [(30.,2.),(90.,6.)]:
            for statistic in ['mean_region_width',*A.PRIMARY]:
                rows.append(dict(model='wave_pinning',C_label=C,f=.01,kind='noise',sigma=1.,experiment='E0',threshold=.5,
                    smoothed=True,primary=True,observation_cadence=80.,statistic=statistic,estimate=x,
                    ci_lower=x*.9,ci_upper=x*1.1,independent_seeds=12,seeds_with_estimand=12))
        response=pd.DataFrame([dict(model='wave_pinning',f=.01,view='main',small_C_label=30.,large_C_label=90.,
            tau_status_small='resolved',tau_status_large='resolved',tau_ratio=2.,evaluable=True,forced=False,
            Dtheta_mobility_proxy=3.,mobility_proxy_valid_bootstrap_fraction=1.)])
        f=pd.DataFrame(rows);first=A.condition_pairs(f,response)
        f[['estimate','ci_lower','ci_upper']]*=1e4;second=A.condition_pairs(f,response)
        np.testing.assert_allclose(first.feature_ratio,second.feature_ratio)
        np.testing.assert_allclose(second.small_estimate,first.small_estimate*1e4)

    def test_sink_name_cannot_replace_native_condition_kind(self):
        with self.assertRaisesRegex(ValueError,'condition-kind identity'):
            A.condition_pairs(pd.DataFrame([dict(kind='survival')]),pd.DataFrame())

    def test_missing_primary_noise_view_fails_before_any_association(self):
        features=pd.DataFrame([dict(model='wave_pinning',C_label=30.,f=.01,kind='noise',sigma=1.,experiment='E0',
            threshold=.3,smoothed=False,primary=False,observation_cadence=40.,statistic='region_count_mean',estimate=1.)])
        responses=pd.DataFrame([dict(model='wave_pinning',f=.01,view='main',small_C_label=30.,large_C_label=90.)])
        with self.assertRaisesRegex(ValueError,'primary A5 noise view'):
            A.condition_pairs(features,responses)

    def test_complete_canonical_csv_adapter_preserves_all_four_model_folds(self):
        features,responses=[],[]
        for model in A.MODELS:
            small,large=(2.5,7.5) if model=='spring' else (30.,90.)
            for i,f in enumerate([.01,.03,.1]):
                responses.append(dict(model=model,f=f,view='main',small_C_label=small,large_C_label=large,
                    tau_status_small='resolved',tau_status_large='resolved',tau_ratio=3./(1+i),evaluable=True,forced=False,
                    Dtheta_mobility_proxy=4./(1+i),mobility_proxy_valid_bootstrap_fraction=1.))
                for C,amount in [(small,1.),(large,2.+i)]:
                    for statistic in A.PRIMARY:
                        features.append(dict(model=model,C_label=C,f=f,kind='noise',sigma=1.,experiment='E0',
                            threshold=.5,smoothed=True,primary=True,observation_cadence=8. if model=='spring' else 80.,
                            statistic=statistic,estimate=amount,ci_lower=.9*amount,ci_upper=1.1*amount,
                            independent_seeds=12,seeds_with_estimand=12))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);a5=root/'a5';ratios=root/'ratios';a5.mkdir();ratios.mkdir()
            pd.DataFrame(features).to_csv(a5/'metrics.csv',index=False)
            pd.DataFrame(responses).to_csv(ratios/'size_ratios.csv',index=False)
            (a5/'complete.json').write_text(json.dumps(dict(passed=True,outputs=dict(metrics=dict(sha256=A.sha(a5/'metrics.csv'))))))
            (ratios/'status.json').write_text(json.dumps(dict(passed=True,complete=True,outputs={'size_ratios.csv':A.sha(ratios/'size_ratios.csv')})))
            result=A.analyze(a5,ratios,root/'output')
            self.assertTrue(result['primary_A6_support'])
            complete=A.read(root/'output/complete.json')
            self.assertEqual(len(complete['primary_features']),2)
            self.assertTrue(all(r['successful_folds']==4 for r in complete['primary_features']))
            self.assertEqual(pd.read_csv(root/'output/LOMO_folds.csv').held_out_model.nunique(),4)

    def test_expanded_context_requires_exact_source_horizons_and_arm(self):
        features=pd.DataFrame([dict(model='wave_pinning',f=.01,C_label=C,native_horizon=H)
            for C,H in [(30.,32000.),(60.,64000.)]])
        response=dict(model='wave_pinning',f=.01,small_C_label=30.,large_C_label=60.,analysis_arm='longest',
            native_horizon_small=32000.,native_horizon_large=64000.)
        part,context=A.matched_pair_context(features,response,'longest')
        self.assertEqual(len(part),2)
        self.assertFalse(context['matched_horizons_across_sizes'])
        self.assertEqual(context['native_horizon_large'],64000.)
        with self.assertRaisesRegex(ValueError,'horizon mismatch'):
            A.matched_pair_context(features,{**response,'native_horizon_large':32000.},'longest')
        with self.assertRaisesRegex(ValueError,'analysis_arm'):
            A.matched_pair_context(features,{k:v for k,v in response.items() if k!='analysis_arm'},'longest')

    def test_common_arm_cannot_hide_unequal_size_horizons_or_mixed_rows(self):
        features=pd.DataFrame([dict(model='wave_pinning',f=.01,C_label=C,native_horizon=H)
            for C,H in [(30.,32000.),(60.,64000.)]])
        response=dict(model='wave_pinning',f=.01,small_C_label=30.,large_C_label=60.,analysis_arm='common',
            native_horizon_small=32000.,native_horizon_large=64000.)
        with self.assertRaisesRegex(ValueError,'equal observation horizons'):
            A.matched_pair_context(features,response,'common')
        mixed=pd.concat([features,pd.DataFrame([dict(model='wave_pinning',f=.01,C_label=60.,native_horizon=32000.)])])
        with self.assertRaisesRegex(ValueError,'mixed horizon'):
            A.matched_pair_context(mixed,{**response,'analysis_arm':'longest'},'longest')

    def test_association_tables_never_pool_arms_or_forcing_families(self):
        with self.assertRaisesRegex(ValueError,'Separate A6 runs'):
            A.association_tables(pd.DataFrame(dict(analysis_arm=['common','longest'],sigma=[1.,1.])))
        with self.assertRaisesRegex(ValueError,'Separate A6 runs'):
            A.association_tables(pd.DataFrame(dict(analysis_arm=['common','common'],sigma=[.5,1.])))

    def test_unmatched_physical_cadences_are_retained_without_entering_sensitivity_fits(self):
        features,responses=self.expanded_fixture()
        small=features[(features.C_label==30.)&(features.statistic==A.PRIMARY[0])].copy()
        large=features[(features.C_label==60.)&(features.statistic==A.PRIMARY[0])].copy()
        small['observation_cadence']=40.;large['observation_cadence']=5.
        small['primary']=False;large['primary']=False
        result=A.condition_pairs(pd.concat([features,small,large]),responses,experiment='common')
        unpaired=result[result.observational_pair_status!='matched']
        self.assertEqual(set(unpaired.cadence_multiple),{.5,.0625})
        self.assertTrue(unpaired.feature_ratio.isna().all())
        self.assertTrue(unpaired.small_estimate.notna().any() and unpaired.large_estimate.notna().any())
        fits=A.association_tables(result)['within_model']
        self.assertTrue((fits[fits.imaging_view!='primary'].n==0).all())

    def test_intentional_unmatched_input_is_coverage_not_a_missing_expected_native_job(self):
        features,responses=self.expanded_fixture();features=features[features.C_label==30.]
        responses['comparison_status']='unavailable_unmatched_designed_input';responses['missing_size']=60.
        responses['evaluable']=False;responses['native_horizon_large']=np.nan;responses['forced_combined']=None
        result=A.condition_pairs(features,responses,experiment='common')
        self.assertTrue(result.feature_ratio.isna().all())
        self.assertTrue((result.forcing_status=='unavailable_unknown').all())
        fits=A.association_tables(result)['H2_primary_features']
        self.assertTrue((fits[fits.primary_view].excluded_unknown_forcing==1).all())
        self.assertTrue((fits[fits.primary_view].unavailable_folds==4).all())

    def test_unknown_combined_forcing_is_not_treated_as_nonforced(self):
        features,responses=self.expanded_fixture();responses['forced_combined']=None
        result=A.condition_pairs(features,responses,experiment='common')
        self.assertTrue(result.forced.isna().all())
        fits=A.association_tables(result)
        self.assertTrue((fits['cross_model'][~fits['cross_model'].forced_inclusive].n==0).all())
        self.assertTrue((fits['cross_model'][fits['cross_model'].forced_inclusive].n==1).all())

    def test_forcing_sigma_must_match_ratio_family(self):
        features,responses=self.expanded_fixture();features['sigma']=2.
        with self.assertRaisesRegex(ValueError,'forcing sigma mismatch'):
            A.condition_pairs(features,responses,experiment='common',sigma=2.)

    def test_unequal_horizon_raw_totals_keep_estimates_but_have_no_size_quotient(self):
        features,responses=self.expanded_fixture();features['experiment']='longest'
        features.loc[features.C_label==60.,'native_horizon']=64000.
        responses['analysis_arm']='longest';responses['native_horizon_large']=64000.
        extra=[]
        for statistic in A.RAW_TOTAL_COUNTS:
            copy=features[features.statistic==A.PRIMARY[0]].copy();copy['statistic']=statistic;extra.append(copy)
        result=A.condition_pairs(pd.concat([features,*extra]),responses,experiment='longest')
        raw=result[result.statistic.isin(A.RAW_TOTAL_COUNTS)]
        self.assertTrue(raw.raw_count_exposure_guard.all())
        self.assertTrue(raw.feature_ratio.isna().all())
        self.assertTrue(raw.small_estimate.notna().all() and raw.large_estimate.notna().all())
        self.assertTrue(result[result.statistic.isin(A.PRIMARY)].feature_ratio.notna().all())


if __name__=='__main__':unittest.main()
