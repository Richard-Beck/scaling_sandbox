import importlib.util
import math
from pathlib import Path
import sys

import numpy as np

CODE=Path(__file__).resolve().parents[1]/'mechanisms'
sys.path.insert(0,str(CODE))
import signature_figure4_e0_20261020 as E
import signature_figure4_analysis_20261020 as A
import signature_figure4_audit_20261020 as Q
import signature_figure4_inventory_20261020 as I


def synthetic_case():
    time=np.arange(0.,101.)
    phase=.1*time
    P=np.exp(1j*phase)
    P[60:65]=.01
    return dict(row=dict(id='analytic',model='example',C_label=20.,f=.01,kind='noise',seed=1,
                         native_horizon=100.,basal_contrast=10.),
        time=time,P=P,C=np.full(101,20.),contrast=np.full(101,2.),regions=np.ones(101,int))


def test_endpoint_cosine_and_continuous_msd_have_distinct_coverage():
    value=E.case_view(synthetic_case(),'main')
    assert value['cc'][9]>value['mc'][9]
    assert value['metrics']['R_valid_fraction']<1
    assert value['metrics']['single_front_fraction']==1
    assert value['metrics']['single_site_fraction']<1


def test_control_forced_union_does_not_double_count_overlapping_failures():
    case=synthetic_case()
    case['P'][50:81]=.01
    case['contrast'][50:81]=.01
    metric=E.case_view(case,'main')['metrics']
    assert math.isclose(metric['forced_union_fraction'],31/51)
    assert metric['regime']=='low_contrast'


def test_bootstrap_pools_eligible_pairs_not_per_seed_curve_means():
    numbers=np.array([[1.,2.],[3.,0.]])
    counts=np.array([[1,2],[9,0]])
    point,draws=E.pooled_bootstrap(numbers,counts,np.array([[2,0],[0,2]]))
    np.testing.assert_allclose(point,[.4,1.])
    np.testing.assert_allclose(draws[0],[1.,1.])
    assert np.isnan(draws[1,1])


def test_tau_crossing_after_missing_lag_is_unavailable():
    curve=np.array([.99,np.nan,.8,.3])
    assert A.persistence(curve,np.ones(4),5)['status']=='unavailable_lag_gap'
    lower,upper,codes=A.persistence_draws(curve[None,:],5)
    assert codes[0]==3 and lower[0]==0 and math.isinf(upper[0])


def test_missing_weak_condition_after_forced_exclusion_cannot_define_fhalf():
    metrics=[dict(model='m',C_label=1,view='main',f=.0001,n_seeds=2,forced=True,tau_status='resolved'),
             dict(model='m',C_label=1,view='main',f=.01,n_seeds=12,forced=False,tau_status='resolved')]
    row=E.fhalf_row(metrics,{},'m',1,'main',False)
    assert row['fhalf_status']=='unavailable_after_forced_exclusion'
    assert row['weakest_f']==.0001
    assert math.isinf(row['fhalf_ci_upper'])


def test_censored_json_has_explicit_null_and_unbounded_flag(tmp_path):
    E.write_json(tmp_path/'receipt.json',dict(bound=math.inf,unbounded=True))
    value=E.read_json(tmp_path/'receipt.json')
    assert value==dict(bound=None,unbounded=True)


def test_no_evaluable_conditions_does_not_vacuously_support_h1():
    ratios=[dict(model='m',view=v,forced=False,evaluable=False,f=.001,tau_ratio_ci_lower=0.,tau_ratio_ci_upper=math.inf)
            for v in E.VIEWS]
    halves=[dict(model='m',view=v,C_label=C,forced_inclusive=inc,fhalf_ci_lower=.001,fhalf_ci_upper=math.inf)
            for v in E.VIEWS for C in (1,3) for inc in (False,True)]
    rows=E.decisions(ratios,halves)
    assert all(r['H1']=='inconclusive' and r['evaluable_conditions']==0 for r in rows)


def test_independent_coordinate_audit_handles_seam_invalid_gaps_and_dilation():
    P=np.exp(1j*np.array([3.,3.2,3.4,-1.,-2.,-1.5]))
    C=np.array([10.,12.,14.,15.,16.,20.]);valid=np.array([1,1,1,0,1,1],bool)
    theta,arc,direct=Q.coordinates(P,C,valid)
    np.testing.assert_allclose(theta[[0,1,2,4,5]],[0,.2,.4,0,.5])
    np.testing.assert_allclose(arc[[0,1,2,4,5]],np.array([0,2.2,4.8,0,9.])/(2*np.pi))
    np.testing.assert_allclose(direct[[0,1,2,4,5]],np.array([0,2.4,5.6,0,10.])/(2*np.pi))
    num,count=Q.direct_msd(theta,valid,[1,2,3])
    np.testing.assert_array_equal(count,[3,1,0])
    np.testing.assert_allclose(num,[.33,.16,0])


def test_spring_compact_basal_uses_blocked_native_layout(tmp_path):
    active=np.array([1.,3.,2.]);inactive=np.array([20.,21.,22.]);q=np.arange(3)/3
    path=tmp_path/'spring.npz'
    np.savez(path,state=np.r_[active,inactive,4.],material_positions=q)
    result,positions,C=I.basal_profile(path,'spring_ring',3)
    np.testing.assert_array_equal(result,active)
    np.testing.assert_array_equal(positions,q)
    assert C==4
