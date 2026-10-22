"""Scientific checks of round-2 estimands, independent of native trajectories."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import revision_round2_figure4_rates_20261022 as R
import signature_figure4_analysis_20261020 as A

def test_endpoint_mask_and_no_gap_interpolation():
    p=np.exp(1j*np.arange(31)*.17);valid=np.ones(31,bool);valid[8:14]=False
    direct=R.seed_lags(p,valid,1.,[2.,5.,17.,2.5])
    n,d,possible=A.heading_cosine(p,valid,20)
    for j,k in enumerate([2,5,17]):
        np.testing.assert_allclose(direct[:,j],[n[k-1],d[k-1],possible[k-1]],atol=1e-12)
    assert np.isnan(direct[:,3]).all()

def test_known_exponential_rate_ratio_and_geometry():
    lag=np.array([250,500,1000,2000,4000,8000.]);a={};b={}
    for x,k,C in [(a,2e-4,20.),(b,2e-5,60.)]:
        x.update(lags=lag,c=np.exp(-k*lag),cov=np.ones(6),Cmean=C,Cboot=np.full(2000,C))
        x['boot']=np.tile(x['c'],(2000,1));x['k']=R.rate(x['c'],lag);x['kboot']=R.rate(x['boot'],lag)
    r=R.ratio_record(a,b,{}, {},0)
    assert r['usable'];np.testing.assert_allclose(r['rho_k'],10);np.testing.assert_allclose(r['rho_geo'],10/9)
    np.testing.assert_allclose(R.secondary(a,b)['secondary_rho_k'],10)

def test_noise_threshold_log_interpolation_and_bounds():
    f=np.array([.001,.01,.1]);r=R.crossing(f,[.9,.6,.4])
    np.testing.assert_allclose(r['estimate'],np.sqrt(.01*.1))
    assert R.crossing(f,[.9,.8,.7])['upper']==np.inf
    assert R.crossing(f,[.4,.3,.2])['lower']==0
    g=R.crossing(f,[.9,np.nan,.4]);assert g['estimate'] is None and g['lower']==.001 and g['upper']==.1

def test_unknown_forced_is_not_unforced():
    assert R.forced_pair('False','') is None
    assert R.forced_pair('False','False') is False
    assert R.forced_pair('True','') is True

def test_undefined_rate_draws_do_not_create_precise_ratio_interval():
    a=dict(lags=np.array([1.]),c=np.array([.8]),boot=np.full((2000,1),.8),cov=np.ones(1),k=np.array([.2]),kboot=np.full((2000,1),.2),Cmean=20.,Cboot=np.full(2000,20.))
    b=dict(lags=np.array([1.]),c=np.array([.9]),boot=np.full((2000,1),.9),cov=np.ones(1),k=np.array([.1]),kboot=np.full((2000,1),.1),Cmean=60.,Cboot=np.full(2000,60.))
    a['kboot'][:100]=np.nan
    result=R.ratio_record(a,b,{}, {},0)
    assert result['rho_k_ci_lower']==0 and result['rho_k_ci_upper']==np.inf
    np.testing.assert_allclose(result['bootstrap_valid_fraction'],.95)
