from pathlib import Path
import sys
import math
import numpy as np

CODE=Path(__file__).resolve().parents[1]/'mechanisms'
sys.path.insert(0,str(CODE))
import signature_figure4_final_audit_20261020 as Q


def metric(tau=1000,*,regime='single_front',coverage=1.,forced=False):
    return dict(tau_status='resolved',tau_lower=tau,tau_upper=tau,tau_ci_lower=tau,tau_ci_upper=tau,
        tau_coverage=coverage,forced_combined=forced,regime=regime)


def ratio(f=0.01,point=2.,ci=(1.5,2.5),*,forced=False):
    return dict(f=f,tau_ratio=point,tau_ratio_lower=point,tau_ratio_upper=point,tau_ratio_ci_lower=ci[0],tau_ratio_ci_upper=ci[1],
        large_multiple_front_fraction=.8,evaluable=True,coverage_small=1.,coverage_large=1.,forced_combined=forced)


def test_independent_weak_censored_bound_and_low_coverage_rules():
    weak=dict(metric(),tau_status='right_censored',tau_lower=8000,tau_upper=None,tau_upper_unbounded=True)
    assert Q.strong(metric(1000),weak,False)=='pass'
    assert Q.strong(metric(6000),weak,False)=='inconclusive'
    assert Q.strong(metric(1000,coverage=.1),weak,False)=='inconclusive'


def test_independent_h1_requires_disjoint_tolerance_and_all_gains():
    small=dict(fhalf_ci_lower=.001,fhalf_ci_upper=.005);large=dict(fhalf_ci_lower=.01,fhalf_ci_upper=.03)
    assert Q.h1([ratio()],small,large,False)=='pass'
    assert Q.h1([ratio(),ratio(.1,.5,(.4,.8))],small,large,False)=='fail'
    assert Q.h1([],small,large,False)=='inconclusive'


def test_independent_h2_does_not_bridge_a_forced_range():
    rows=[ratio(.01,3.,(2.,4.)),ratio(.03,2.,(1.5,2.5),forced=True),ratio(.1,1.,(.8,1.2))]
    assert Q.h2(rows,False)=='inconclusive'
    assert Q.h2(rows,True)=='pass'


def test_independent_existence_rule_requires_all_preset_sigma_observations():
    rows=[dict(metric(regime='multiple_fronts'),sigma=s) for s in (1.,2.)]
    assert Q.corroborated(rows,'multiple_fronts',False)=='inconclusive'
    rows.extend(dict(metric(),sigma=s) for s in (.25,.5))
    assert Q.corroborated(rows,'multiple_fronts',False)=='pass'


def test_independent_joint_seed_effects_retain_zero_denominator():
    result=Q.independent_effect(np.array([[1,1],[1,1]]),np.array([[3,.25],[3,.25]]))
    assert result['region_count']==(3.,3.) and result['normalized_contrast']==(.25,.25)
    result=Q.independent_effect(np.array([[0,1],[0,1]]),np.array([[3,.25],[3,.25]]))
    assert result['region_count']==(0.,math.inf)


def test_independent_window_changes_are_not_hidden_by_heading_loss():
    assert Q.window(metric(),metric(regime='low_contrast',coverage=0.))=='fail'
    assert Q.window(metric(),metric(coverage=0.))=='inconclusive'
    assert Q.window(metric(),metric(10))=='fail'


def test_independent_fhalf_does_not_replace_weak_censoring_with_a_point():
    assert Q.fhalf([.0001,.03,.1],[8000,2000,1000],[math.inf,2000,1000])==(.0001,.03)
    lower,upper=Q.fhalf([.0001,.03],[10000,1000],[10000,1000])
    assert lower==upper and .0001<lower<.03
    assert Q.fhalf([.0001,.03],[10000,0],[10000,math.inf])==(.0001,math.inf)


def test_independent_bound_eligibility_requires_a_real_covered_crossing():
    row=dict(metric(tau=0),tau=None,tau_status='unresolved_below_two_cadences',tau_upper=10,cadence=5,tau_crossing_bracket=[0,5])
    assert Q.eligible(row)
    assert not Q.eligible(dict(row,tau_coverage=.49))
    assert not Q.eligible(dict(row,tau_crossing_bracket=None))
    assert not Q.eligible(dict(row,tau_status='unavailable_no_contiguous_lag_curve'))
