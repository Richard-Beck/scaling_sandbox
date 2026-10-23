"""Independent numerical checks of the edge/covariance calculation and exports."""
from pathlib import Path
import os,sys,math,json
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
from pypdf import PdfReader
import figure3_round5_edge_20261023 as edge
OUT=edge.OUT;DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'

def run():
    checks={};C=20.;n=400;x=(np.arange(n)+.5)*C/n;u=2+np.cos(2*np.pi*x/C)
    e=edge.edges(u,C);np.testing.assert_allclose(e,[15.,5.],atol=1e-10)
    shifted=2+np.cos(2*np.pi*(x-3.25)/C);es=edge.edges(shifted,C)
    np.testing.assert_allclose(edge.differences(es,e,C),[3.25,3.25],atol=1e-10)
    checks['half_contrast_crossings_and_periodic_translation']=True
    assert edge.edges(np.ones(n),C) is None
    assert edge.edges(2+np.cos(4*np.pi*x/C),C) is None
    checks['homogeneous_and_multiple_regions_rejected']=True
    refinements={}
    for platform,C,ell,tau,knot in [('chemical',20.,1.5,50.,5.),('spring',1.54,.15,5.,.125)]:
        coarse,space=edge.covariance_integral(64,C,ell,tau,knot)
        refined,_=edge.covariance_integral(64,C,ell,tau,knot,phases=48,steps_per_knot=24)
        err=float(np.max(abs(coarse-refined))/np.max(abs(refined)))
        assert err<1e-4,(platform,err)
        # Independent continuous-OU local lognormal identity, a convergent
        # series, is approached by the declared interpolated-knot process.
        exact=tau*sum(1./(math.factorial(m)*m) for m in range(1,35))
        continuous_difference=float(abs(refined[0]-exact)/exact)
        # The native interpolation is intentionally different from continuous OU;
        # use a separate small-knot limit for an analytic check below.
        # Direct quadratic form checks circular FFT normalization and signs.
        rng=np.random.default_rng(20261023);g=rng.normal(size=64)
        matrix=coarse[(np.arange(64)[:,None]-np.arange(64)[None,:])%64]
        direct=float(g@matrix@g);fft=float(g@np.fft.irfft(np.fft.rfft(coarse)*np.fft.rfft(g),n=64))
        np.testing.assert_allclose(direct,fft,rtol=1e-12,atol=1e-11)
        assert np.linalg.eigvalsh(matrix).min()>-1e-9
        refinements[platform]=dict(relative_integration_refinement=err,continuous_OU_log_normal_limit_difference=continuous_difference,direct_quadratic_form=direct)
    limit,_=edge.covariance_integral(16,20.,1.5,50.,.05,phases=8,steps_per_knot=2)
    exact_limit=50.*sum(1./(math.factorial(m)*m) for m in range(1,35))
    limit_error=float(abs(limit[0]-exact_limit)/exact_limit)
    assert limit_error<.001
    checks['small_knot_continuous_OU_log_normal_identity_error']=limit_error
    checks['covariance_quadrature_and_fft']=refinements
    zero,_=edge.covariance_integral(16,20.,1.5,50.,5.,sigma=0,phases=2,steps_per_knot=2)
    np.testing.assert_array_equal(zero,np.zeros(16));checks['zero_sigma_zero_diffusion']=True
    frame=pd.read_csv(DATA/'edge_prediction_comparison.csv');assert len(frame)==120;assert frame.defined.sum()==118
    finite=frame[frame.defined]
    assert finite.maximum_linearity_error.max()<.005
    assert finite.maximum_pulse_gain_error.max()<.025
    assert finite.mean_stationary_rhs.max()<=1e-6
    assert finite.mean_tail_scaled_change.max()<=1e-4
    assert finite.adjoint_residual.max()<1e-6
    np.testing.assert_allclose(finite.Ds_predicted_independent,(finite.D_left+finite.D_right)/4,rtol=1e-12)
    # Prediction must explicitly contain measured pulse gain; this check would
    # fail if the native adjoint amplitudes were accidentally used directly.
    native=pd.read_csv(DATA/'kernel_pulse_validation.csv')
    for r in finite.itertuples():
        q=native[(native.model==r.model)&(native['size']==r.size)&(native.fraction==r.fraction)].sort_values('edge')
        scale=(q.measured_edge_gain/q.predicted_phase_gain).to_numpy()
        np.testing.assert_allclose([r.left_measured_gain_scale,r.right_measured_gain_scale],scale,rtol=1e-12)
        np.testing.assert_allclose(r.Ds_predicted_independent,r.raw_native_Ds_independent*np.mean(scale**2),rtol=1e-12)
    checks['120_terminal_118_defined_two_homogeneous']=True
    checks['measured_pulse_gain_used_for_D']=True
    settling=[]
    for c in edge.read(edge.RUN/'manifest.json')['jobs']:
        p=edge.RUN/'workers'/f"{c['slot']:03d}"/'pulse_responses.npz'
        if not p.is_file():continue
        with np.load(p) as z:
            delta=z['displacement'];settling.append(float(np.max(abs(delta[:,-1]-delta[:,-11]))/np.max(abs(delta[:,-1]))))
    assert max(settling)<.001
    checks['maximum_final_half_tau_relative_response_change']=max(settling)
    checks['maximum_response_scaling_error']=float(finite.maximum_linearity_error.max())
    checks['maximum_native_gain_validation_error']=float(finite.maximum_pulse_gain_error.max())
    manifest=edge.read(edge.RUN/'manifest.json')
    assert edge.sha(edge.RUN/'manifest.json')==edge.read(edge.RUN/'freeze.json')['manifest_sha256']
    for files in ('source_sha256','input_sha256'):
        for name,digest in manifest[files].items():assert edge.sha(edge.RUN/name)==digest,name
    checks['frozen_campaign_sources_and_inputs_unchanged']=True
    # Basal and mean physical pool equality, independently of stationarity.
    import figure4_round3_campaign_20261022 as source
    errors={}
    for platform in ('chemical','spring'):
        native,_=source.native_module(edge.RUN,platform)
        maximum=0.
        for c in manifest['jobs']:
            if c['platform']!=platform:continue
            folder=edge.RUN/'workers'/f"{c['slot']:03d}"
            with np.load(edge.RUN/c['input']) as z:y0=z['state'].copy()
            with np.load(folder/'mean_state.npz') as z:y=z['state'].copy()
            if platform=='chemical':
                grid=native.ring.RingGrid(c['size'],c['cells']);model=native.ring.native_model(c['model'])
                p0=native.ring.pool_amounts(model,grid,y0);p=native.ring.pool_amounts(model,grid,y.reshape(y0.shape))
            else:
                n=c['cells'];p0=np.array([y0[-1]*np.mean(y0[:n]+y0[n:-1])]);p=np.array([y[-1]*np.mean(y[:n]+y[n:-1])])
            assert np.min(y)>=-5e-10
            maximum=max(maximum,float(np.max(abs(p/p0-1))))
        assert maximum<1e-7;errors[platform]=maximum
    checks['mean_preparation_physical_pool_max_relative_error']=errors
    preservation=edge.read(META/'preserved_AB.json');assert preservation['all_identical'];checks['AB_quantitative_data_preserved']=True
    indep=pd.read_csv(DATA/'edge_MSD_curves.csv')
    np.testing.assert_allclose(indep.domain_MSD,(indep.left_MSD+indep.right_MSD+2*indep.edge_cross_product)/4,rtol=1e-12,atol=1e-14)
    assert ((indep.increment_correlation>=-1.000001)&(indep.increment_correlation<=1.000001)).all()
    checks['edge_MSD_covariance_identity']=True
    assert len(PdfReader(OUT/'figure3_round5.pdf').pages)==1
    s=edge.read(META/'supplement_preservation.json');assert len(PdfReader(OUT/'figure3_round5_supplement.pdf').pages)==s['total_pages']
    section=edge.read(OUT/'legend.json')
    def walk_paths(value):
        if isinstance(value,dict):
            for k,v in value.items():
                if k in ('path','native_pdf_page','main_pdf','supplement_pdf','reader'):assert (edge.ROOT/v).is_file(),v
                else:walk_paths(v)
        elif isinstance(value,list):
            for v in value:walk_paths(v)
    walk_paths(section);checks['deliverable_PDFs_and_legend_links']=True
    result=dict(generated_utc=edge.now(),passed=True,checks=checks)
    edge.write(META/'numerical_checks.json',result);print(json.dumps(result,indent=2));return result
if __name__=='__main__':run()
