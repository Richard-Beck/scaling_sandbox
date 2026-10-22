"""Scientific invariants and independent analytic cases for round-3 estimators."""
import sys
from pathlib import Path
import unittest
import tempfile
import json
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import figure4_round3_estimators_20261022 as E
import figure4_round3_analysis_20261022 as A
import figure4_round3_spring_coordinate_sensitivity_20261022 as S


class Estimators(unittest.TestCase):
    def test_incomplete_campaign_cannot_replace_provisional_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);campaign=root/'campaign';campaign.mkdir();output=root/'figure'
            (campaign/'manifest.json').write_text(json.dumps({'conditions':[]}));(campaign/'status.json').write_text(json.dumps({'complete':False,'status':'RUNNING'}))
            with self.assertRaisesRegex(ValueError,'precision qualification must be COMPLETE'):A.compact(output,campaign)
            self.assertFalse(output.exists())

    def test_old_two_seed_conditions_compare_against_full_new_cohort(self):
        old=dict(arm='common',view='main',model='wave_pinning',C_label='20',f='.01',lag='2000',n_seeds='2',k='.2',k_ci_lower='.1',k_ci_upper='.3')
        key=('wave_pinning',20.,.01)
        def summary(k,n):return dict(lags=np.array([2000.]),k=np.array([k]),kboot=np.full((2000,1),k),n_seeds=n)
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp);(directory/'round2_matched_P_condition_lag_metrics.csv').touch()
            def read(path):
                return [old] if Path(path).name in ['condition_lag_metrics.csv','round2_matched_P_condition_lag_metrics.csv'] else []
            with patch.object(A,'read_csv',side_effect=read):A.export_round2_comparison(directory,{key:summary(.12,12)},{key:summary(.1,96)})
            row=A.read_csv(directory/'round2_comparison.csv')[0]
            self.assertEqual(int(row['n_old_seeds']),2);self.assertEqual(int(row['n_new_seeds']),96)
            self.assertAlmostEqual(float(row['round3_first12_value']),.12);self.assertAlmostEqual(float(row['full_round3_value']),.1)

    def test_fhalf_comparison_preserves_variable_seed_counts_by_input(self):
        def moment(f,n):return dict(arm='common',view='main',model='wave_pinning',C_label='20',f=str(f),lag='2000',n_seeds=str(n),k='.2',k_ci_lower='.1',k_ci_upper='.3')
        def tol(value):return dict(arm='common',view='main',forced_inclusive='True',model='wave_pinning',C_label='20',reference_lag='2000',fhalf=str(value),fhalf_ci_lower='.01',fhalf_ci_upper='.1',status='resolved')
        def summary(n):return dict(lags=np.array([2000.]),k=np.array([.1]),kboot=np.full((2000,1),.1),n_seeds=n)
        keys=[('wave_pinning',20.,f) for f in [.01,.1]]
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            for name in ['round2_matched_P_condition_lag_metrics.csv','round2_matched_P_noise_tolerance.csv']:(directory/name).touch()
            def read(path):
                path=Path(path)
                if path.name in ['condition_lag_metrics.csv','round2_matched_P_condition_lag_metrics.csv']:return [moment(.01,2),moment(.1,12)]
                if path.name=='round3_first12_noise_tolerance.csv':return [tol(.04)]
                if path.name=='round2_matched_P_noise_tolerance.csv':return [tol(.025)]
                if path==directory/'noise_tolerance.csv':return [tol(.03)]
                return [tol(.02)]
            with patch.object(A,'read_csv',side_effect=read):A.export_round2_comparison(directory,{k:summary(12) for k in keys},{keys[0]:summary(96),keys[1]:summary(192)})
            row=next(r for r in A.read_csv(directory/'round2_comparison.csv') if r['measure']=='fhalf_first_downcrossing')
            self.assertEqual(json.loads(row['old_seed_counts_by_f']),{'0.01':2,'0.1':12})
            self.assertEqual(json.loads(row['new_seed_counts_by_f']),{'0.01':96,'0.1':192})
            self.assertEqual(row['n_old_seeds'],'');self.assertAlmostEqual(float(row['full_round3_value']),.03)

    def test_endpoint_cosine_with_missing_middle(self):
        # Endpoint-only eligibility keeps a lag pair even if its interior is invalid.
        angle=np.arange(9)*.37;P=np.exp(1j*angle);valid=np.ones(9,bool);valid[3]=False
        moment=E.seed_moments(P,valid,1.,[2.])
        self.assertEqual(moment[1,0],5)
        self.assertAlmostEqual(moment[0,0]/moment[1,0],np.cos(.74),places=13)

    def test_explicit_moment_over_amount_changes_eligibility(self):
        archive={'M1_real':np.full(5,.5),'M1_imag':np.zeros(5),'M0':np.full(5,10.),'P_real':np.ones(5),'P_imag':np.zeros(5)}
        P=A.normalized_P(archive)
        self.assertEqual(E.seed_moments(P,abs(P)>=.1,1.,[1.])[1,0],0)
        self.assertEqual(E.seed_moments(np.ones(5),np.ones(5,bool),1.,[1.])[1,0],4)

    def test_whole_seed_bootstrap_has_seed_level_support(self):
        moments=np.array([[[90.],[100.],[100.]],[[-.9],[1.],[1.]]])
        summary=E.summarize(moments,[20.,20.],'known_two_seed_case',np.array([1.]),np.array([1.]))
        support=np.unique(np.round(summary['Cboot'][:,0],12))
        np.testing.assert_allclose(support,[-.9,89.1/101,.9],atol=1e-12)

    def test_precision_cli_and_final_compact_use_identical_bootstrap(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp);records=[];moments=[];time=np.arange(801)*5.;target,lags=E.lag_grid('wave_pinning',time[-1]-time[0],5.)
            for seed in range(2):
                theta=np.cumsum(np.random.default_rng(seed+70).normal(0,np.sqrt(.02),len(time)));P=np.exp(1j*theta);path=directory/f'seed{seed}.npz'
                np.savez(path,time=time,P_real=P.real,P_imag=P.imag,circumference=np.array(20.));moments.append(E.seed_moments(P,abs(P)>=.1,5.,lags))
                records.append(dict(condition_id='wave_pinning:C20:f0.01',model='wave_pinning',size=20.,fraction=.01,compact=str(path),status='COMPLETE'))
            index=directory/'index.csv';output=directory/'precision.json';E.write_csv(index,records)
            with patch('sys.argv',['estimator','precision','--input',str(index),'--output',str(output)]):E.main()
            report=json.loads(output.read_text())['conditions'][0];expected=E.summarize(moments,[20.,20.],('compact',('wave_pinning',20.,.01)),target,lags);j=E.analysis_lag(expected)
            self.assertEqual(report['analysis_lag'],float(lags[j]));self.assertAlmostEqual(report['se_log_k'],float(np.std(np.log(expected['kboot'][:,j]),ddof=1)),places=13)

    def test_physical_brownian_msd_recovers_diffusivity(self):
        rng=np.random.default_rng(19);dt=20.;Ds=.002;C=20.
        x=np.r_[0.,np.cumsum(rng.normal(0,np.sqrt(2*Ds*dt),10000))]
        heading=2*np.pi*x/C
        measure,_=E.drift_jump(np.arange(len(x))*dt,C,np.zeros(len(x),int),[{0:a} for a in heading],max_lag=100.)
        self.assertLess(abs(measure['Ds']/Ds-1),.06)

    def test_constant_physical_velocity_has_quadratic_msd(self):
        time=np.arange(501)*2.;velocity=.003;C=20.;heading=2*np.pi*(3+velocity*time)/C
        measure,curve=E.drift_jump(time,C,np.zeros(len(time),int),[{0:a} for a in heading])
        np.testing.assert_allclose([r['MSD'] for r in curve],[(velocity*r['lag'])**2 for r in curve],rtol=1e-10,atol=1e-13)
        self.assertAlmostEqual(measure['MSD_exponent'],2.,places=9)

    def fixture(self):
        groups={}
        for c in E.CHEMICAL_SIZES:
            lags=np.array([250.,500.,1000.,2000.,4000.,8000.,16000.]);k=.4/c**2;C=np.exp(-k*lags)
            groups[c]=dict(lags=lags,C=C,k=np.full(7,k),coverage=np.ones(7),circumference=c,
                           circumference_boot=np.full(2000,c),kboot=np.full((2000,7),k),n_seeds=12)
        return groups

    def test_geometry_alpha_and_full_fixed_neighbors(self):
        groups=self.fixture();r=E.local_alpha(groups,40)
        self.assertEqual(r['status'],'resolved');self.assertAlmostEqual(r['alpha'],2.,places=12)
        del groups[30.]
        self.assertEqual(E.local_alpha(groups,40)['status'],'unavailable_missing_fixed_neighbor')
        self.assertEqual(E.fixed_neighborhood(20),(20.,30.))
        self.assertEqual(E.fixed_neighborhood(60),(50.,60.))

    def test_single_analysis_nearest_half(self):
        # .30 is closer in C than .71, even though the log-log score favors .71.
        self.assertEqual(E.analysis_lag(dict(C=np.array([.30,.71]),coverage=np.array([.1,1.]))),0)

    def test_fhalf_first_downcross_log_strength(self):
        r=E.crossing([.01,.1,1],[.9,.3,.8])
        self.assertAlmostEqual(r['estimate'],.01*10**(2/3))
        self.assertEqual(E.crossing([.01,.1,1],[.9,np.nan,.3])['status'],'bounded_missing_gap')

    def test_no_switch_is_measured_zero_jump_contribution(self):
        rng=np.random.default_rng(31);heading=np.cumsum(rng.normal(0,.02,401));heads=[{0:x} for x in heading]
        r,msd=E.drift_jump(np.arange(401)*5.,20.,np.zeros(401,int),heads)
        self.assertEqual(r['lambda_eff'],0.);self.assertEqual(r['valid_time'],2000.)
        self.assertGreater(r['Ds'],0.);self.assertTrue(np.isnan(r['mean_cos_jump']))

    def test_switch_time_keeps_boundary_intervals(self):
        ids=np.r_[np.zeros(10,int),np.ones(10,int)];heads=[{int(i):float(i)} for i in ids]
        r,_=E.drift_jump(np.arange(20)*2.,20.,ids,heads)
        self.assertEqual(r['valid_time'],38.);self.assertEqual(r['n_switches'],1)
        self.assertAlmostEqual(r['lambda_eff'],(1-np.cos(1))/38)

    def test_initial_short_tenure_is_observed_before_confirmed_switch(self):
        ids=np.r_[np.ones(2,int),np.full(5,2,int)]
        np.testing.assert_array_equal(E.smooth_flicker(ids,5),ids)
        result,_=E.drift_jump(np.arange(len(ids)),20.,ids,[{1:0.,2:1.} for _ in ids])
        self.assertEqual(result['n_switches'],1);self.assertEqual(result['n_scored_switches'],0)
        self.assertAlmostEqual(result['lambda_s'],1/6)

    def test_uniform_dilation_does_not_move_intrinsic_heading(self):
        C=np.linspace(2.,3.,101);ids=np.zeros(101,int);maps=[{0:np.pi/3} for _ in ids]
        values=S.seed_sums(C,ids,maps,1.,20)
        self.assertEqual(float(values['intrinsic_arc'][0].sum()),0.)
        self.assertGreater(float(values['literal_fixed_origin'][0].sum()),0.)

    def test_fixed_circumference_coordinate_conventions_agree(self):
        C=np.full(101,2.);heading=np.linspace(.2,1.4,101);ids=np.zeros(101,int)
        values=S.seed_sums(C,ids,[{0:a} for a in heading],1.,20)
        np.testing.assert_allclose(values['intrinsic_arc'],values['literal_fixed_origin'],atol=1e-12)

    def test_detector_free_does_not_bridge_seed_separator(self):
        P=np.r_[np.exp(1j*np.arange(15)*.01),complex(np.nan),np.exp(1j*np.arange(15)*.02)]
        labels=np.r_[np.zeros(15,int),-1,np.ones(15,int)]
        result=E.wrapped_mixture_mle(P,20.,1.,2,labels)
        self.assertEqual(result['n_increments'],26)

    def test_detector_free_twice_cadence_uses_only_endpoint_magnitude(self):
        P=np.array([1.,1.,.01,1.,1.],complex)
        result=E.wrapped_mixture_mle(P,20.,1.,2,np.zeros(5,int))
        self.assertEqual(result['n_increments'],1)

    def test_wrapped_mixture_recovers_known_diffusion_and_uniform_rate(self):
        rng=np.random.default_rng(51);D=.005;lam=.05;dt=1.;n=12000
        increments=rng.normal(0,np.sqrt(2*D*dt),n);reset=rng.random(n)<-np.expm1(-lam*dt)
        increments[reset]=rng.uniform(-np.pi,np.pi,reset.sum())
        result=E.wrapped_mixture_mle(np.exp(1j*np.r_[0,np.cumsum(increments)]),20.,dt)
        self.assertEqual(result['status'],'resolved')
        self.assertLess(abs(result['Dtheta']/D-1),.08)
        self.assertLess(abs(result['lambda_eff']/lam-1),.15)

    def test_uniform_same_rate_size_curves_identical(self):
        a=A.simulate_null(100,20.,0.,1/7000,24000.,20.,19,True)
        b=A.simulate_null(100,60.,0.,1/7000,24000.,20.,19,True)
        np.testing.assert_array_equal(a,b)

    def test_exact_continuous_reset_acf(self):
        # Independent non-overlapping ensemble endpoint identity E cos = exp(-k t).
        theta=A.simulate_null(40000,20.,.002,.001,1000.,20.,17)
        empirical=np.cos(theta[:,-1]-theta[:,0]).mean()
        analytic=np.exp(-((2*np.pi/20)**2*.002+.001)*1000)
        self.assertLess(abs(empirical-analytic),.012)


if __name__=='__main__':unittest.main()
