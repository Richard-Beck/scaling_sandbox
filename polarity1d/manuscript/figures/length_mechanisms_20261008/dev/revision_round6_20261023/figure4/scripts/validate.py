"""Independent analytic and complete-export checks; no new native noise runs."""
from common import *
import numpy as np
from mode_analysis import fourier,sufficient,aggregate_modes
from exemplars import exact_decomposition
from pypdf import PdfReader

def main():
    checks={}
    # Analytic continuous first harmonic sampled at native cell centres.
    N=400;theta=2*np.pi*(np.arange(N)+.5)/N;t=np.arange(401.)*5.;heading=.001*t
    p=1+.4*np.cos(theta[None,:]-heading[:,None]);a=fourier(p)
    error=max(np.max(abs(a[:,0]-1)),np.max(abs(a[:,1]-.2*np.exp(1j*heading))),np.max(abs(a[:,2:])))
    assert error<2e-14;checks['analytic_Fourier_max_error']=float(error)
    # Irregularly thinned exact endpoints must still reproduce the analytic
    # phase correlation of a constant-amplitude rotating harmonic.
    take=np.unique(np.linspace(0,len(t)-1,271,dtype=int));s=sufficient(t[take],a[take],'wave_pinning')
    eligible=s['counts']>0
    complex_corr=s['num'][eligible,1]/s['den'][eligible,1]
    phase_corr=s['phase_num'][eligible]/s['phase_den'][eligible]
    analytic=np.cos(.001*s['lags'][eligible])
    err=float(max(np.max(abs(complex_corr-analytic)),np.max(abs(phase_corr-analytic))))
    assert err<2e-12;checks['irregular_endpoint_phase_complex_correlation_error']=err
    # Exponential complex correlations with known mode-rate q² scaling.
    lags=np.array([5.,10.,20.,40.,80.,160.,320.,640.,1280.]);g=np.array([0,.0004,.0016,.0036])
    counts=np.full(len(lags),10000.);num=np.exp(-lags[:,None]*g)*counts[:,None];den=np.ones_like(num)*counts[:,None]
    seed=dict(lags=lags,target=lags,num=num,den=den,counts=counts,phase_num=num[:,1],phase_den=counts,spectrum=np.array([1,.04,.01,.005]),frames=10001)
    rates,curves=aggregate_modes([seed]*4,dict(condition_id='analytic',model='wave_pinning',size=20.,fraction=.01),identity='analytic_validation')
    err=max(abs(r['gamma']-g[r['n']]) for r in rates[1:]);assert err<1e-12
    checks['analytic_exponential_rate_max_error']=float(err)
    # Uniform background has zero first moment. One complete peak basin covers
    # all excess, so exact native-mass moment reconstruction must close by algebra.
    exact=exact_decomposition(t,p,np.full(len(t),20.),p[0],100.)
    full=next(r for r in exact if r['convention']=='exact_full_mask')
    assert abs(full['remainder_fraction'])<1e-11
    assert max(r['max_closure_error'] for r in exact)<1e-11
    checks['rigid_translation_exact_native_mass_R']=full['remainder_fraction']
    checks['synthetic_all_convention_loss_closure']=max(r['max_closure_error'] for r in exact)
    # Full export checks inspect independent outcomes rather than copying a
    # numerical implementation into a second test.
    ex=csvread(DATA/'exemplar_exact_masks.csv');assert ex
    maxerr=max(numrow(r,'max_closure_error') for r in ex);assert maxerr<1e-10
    checks['actual_exact_mask_max_heading_closure_error']=maxerr
    complete=list((OUT/'exemplar_analysis').glob('c*/complete.json'));assert len(complete)==260
    assert all(read(p).get('algorithm_sha256')==digest(Path(__file__).with_name('exemplars.py')) for p in complete)
    assert all(read(p).get('mode_algorithm_sha256')==digest(Path(__file__).with_name('mode_analysis.py')) for p in complete)
    checks['all260_exemplar_algorithms_match']=True
    closure=csvread(DATA/'H_closure_intervals.csv');assert len(closure)==260
    assert all(numrow(r,'kpred_lower')<=numrow(r,'kpred_upper')+1e-15 for r in closure)
    assert len({numrow(r,'f') for r in closure})==16
    checks['closure_all260_conditions_all16_inputs']=True
    for r in csvread(DATA/'J_takeover_selection.csv'):assert numrow(r,'minimum_births')>=200
    assert all(numrow(r,'distance_bin')>0 for r in csvread(DATA/'J_takeover_points.csv'))
    checks['takeover_total_birth_gate_and_positive_distance_plot']=True
    checks['main_pdf_pages']=len(PdfReader(OUT/'figure4_main.pdf').pages)
    checks['supplement_pdf_pages']=len(PdfReader(OUT/'figure4_supplement.pdf').pages)
    write(META/'scientific_validation.json',dict(passed=True,checks=checks))
    print(json.dumps(checks,indent=2))

def numrow(r,k):return float(r[k])
if __name__=='__main__':main()
