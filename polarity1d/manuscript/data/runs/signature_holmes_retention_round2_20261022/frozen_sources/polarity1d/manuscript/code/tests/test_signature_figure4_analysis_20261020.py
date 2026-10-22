"""Analytic checks of E0 masks, winding, finite bounds and seed bootstrap."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_figure4_analysis_20261020 as analysis


def test_fft_heading_cosines_match_endpoint_masked_direct_sums():
    theta=np.array([2.9,-3.,-2.5,.8,1.4,-1.]);P=2*np.exp(1j*theta)
    valid=np.array([True,True,False,True,True,True])
    sums,counts,possible=analysis.heading_cosine(P,valid,3)
    for k in range(1,4):
        mask=valid[:-k]&valid[k:]
        assert sums[k-1]==pytest.approx(np.cos(theta[k:]-theta[:-k])[mask].sum(),abs=2e-14)
        assert counts[k-1]==mask.sum() and possible[k-1]==len(theta)-k


def test_zero_usable_heading_pairs_are_missing_rather_than_zero_correlation():
    numerator,count,possible=analysis.heading_cosine(np.ones(5,dtype=complex),np.zeros(5,dtype=bool),2)
    assert np.isnan(analysis.ratio(numerator,count)).all()
    np.testing.assert_array_equal(possible,[4,3])


def test_unwrapped_msd_cannot_bridge_an_invalid_heading_gap():
    P=np.exp(1j*np.array([0.,.2,2.,1.2,1.4]));valid=np.array([True,True,False,True,True])
    x=analysis.angular_coordinate(P,valid);sums,count,possible=analysis.section_msd(x,valid,3)
    assert count[0]==2 and sums[0]==pytest.approx(.08,abs=1e-13)
    assert count[1]==0 and count[2]==0


def test_msd_fft_matches_direct_contiguous_track_squared_displacements():
    x=np.cumsum(np.random.default_rng(2).normal(size=30));valid=np.ones(30,dtype=bool)
    valid[8:11]=False
    sums,count,_=analysis.section_msd(x,valid,10)
    for k in range(1,11):
        pairs=[(x[j+k]-x[j])**2 for j in range(len(x)-k) if valid[j:j+k+1].all()]
        assert sums[k-1]==pytest.approx(sum(pairs),abs=1e-10)
        assert count[k-1]==len(pairs)


def test_arc_increment_gauge_does_not_turn_uniform_dilation_into_site_motion():
    P=np.exp(1j*np.full(5,2.));C=np.arange(10.,15.);valid=np.ones(5,dtype=bool)
    np.testing.assert_allclose(analysis.physical_coordinate(P,C,valid),0.,atol=1e-14)
    P=np.exp(1j*np.arange(5)*.1)
    expected=np.r_[0.,np.cumsum(.1*(C[1:]+C[:-1])/(4*np.pi))]
    np.testing.assert_allclose(analysis.physical_coordinate(P,C,valid),expected,atol=1e-14)


def test_persistence_keeps_short_resolution_and_finite_horizon_censoring():
    short=analysis.persistence(np.array([0.,0.]),np.ones(2),5.)
    assert short['status']=='unresolved_below_two_cadences' and short['upper']==10.
    bounded=analysis.persistence(np.array([.99,.98,.97]),np.ones(3),5.)
    assert bounded['status']=='right_censored' and bounded['lower']==15. and np.isinf(bounded['upper'])
    gap=analysis.persistence(np.array([np.nan,.2]),np.ones(2),5.)
    assert gap['status']=='unavailable_lag_gap'


def test_bootstrap_persistence_intervals_do_not_replace_bounds_with_fitted_points():
    curves=np.array([[.9,.8,.7],[.1,.05,.01],[np.nan,np.nan,np.nan]])
    low,high,status=analysis.persistence_draws(curves,5.)
    np.testing.assert_array_equal(status,[2,1,3])
    assert low[0]==15. and np.isinf(high[0]) and high[1]==10. and np.isinf(high[2])


def test_two_censored_times_do_not_produce_a_ratio_of_their_lower_bounds():
    lo,hi=analysis.interval_ratio(20.,np.inf,10.,np.inf)
    assert lo==0. and np.isinf(hi)


def test_short_lag_diffusivity_is_half_the_signed_free_intercept_slope():
    msd=3.+4.*np.arange(1,11)*.5
    assert analysis.signed_diffusivity(msd,.5)==pytest.approx(2.)
    assert analysis.signed_diffusivity(-msd,.5)==pytest.approx(-2.)


def test_relative_contrast_qualification_supersedes_absolute_native_convenience_flag():
    fractions,label=analysis.regime(np.array([.04,.2,.2,.2]),np.array([1,1,2,2]),1.)
    np.testing.assert_allclose(fractions,[.25,.5,.25,0.])
    assert label=='mixed'


def test_noise_tolerance_uses_log_strength_and_keeps_a_censored_weak_baseline():
    resolved=analysis.conservative_fhalf([.001,.01],[100.,40.],[100.,40.])
    assert resolved['estimate']==pytest.approx(np.exp(np.log(.001)+(5/6)*np.log(10)))
    bounded=analysis.conservative_fhalf([.001,.01],[100.,40.],[np.inf,40.])
    assert bounded['estimate'] is None and bounded['upper']==.01


def test_whole_seed_bootstrap_counts_sum_to_the_number_of_realizations():
    weights=analysis.bootstrap_weights(2,2000,23)
    assert weights.shape==(2000,2)
    assert np.all(weights.sum(axis=1)==2)
    assert set(map(tuple,weights))=={(0,2),(1,1),(2,0)}
