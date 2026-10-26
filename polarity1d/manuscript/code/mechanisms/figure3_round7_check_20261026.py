"""Integrity, scientific numeric and display checks for the round7 Figure3."""
from pathlib import Path
import os,sys,json,hashlib
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
from scipy.optimize._numdiff import approx_derivative
from pypdf import PdfReader
import figure3_conduction_20261026 as campaign
ROOT=campaign.ROOT;OUT=campaign.OUT;DATA=campaign.DATA;META=campaign.META

def main():
    d=pd.read_csv(DATA/'conduction_sweep.csv');checks={}
    assert len(d)==100 and d.defined.all()
    assert d.groupby('model').size().to_dict()==dict(holmes_model3=30,spring_ring=35,wave_pinning=35)
    assert d.rhs_residual.max()<1e-6 and d.adjoint_residual.max()<1e-7
    assert d.translation_eigenvalue.abs().max()<1e-6 and d.right_translation_overlap.min()>.95
    assert d.negative_min.min()>0 and np.isfinite(d.gain).all()
    checks['all100_native_conditions_completed']=True;checks['max_stationary_RHS']=float(d.rhs_residual.max());checks['max_adjoint_residual']=float(d.adjoint_residual.max());checks['min_translation_overlap']=float(d.right_translation_overlap.min())
    native,_=campaign.source.native_module(campaign.edge.RUN,'spring');n=200;closure=native.Closure(2.,n)
    # Independent centered difference in C verifies the additional inactive
    # diffusion contribution to the complete Spring Jacobian.
    with np.load(DATA/'conduction_states_spring_ring.npz') as z:y=z['D1_C2'].copy()
    n=(len(y)-1)//2;C=y[-1];v=y[n:-1];D=1.;step=1e-6
    lap=np.roll(v,1)-2*v+np.roll(v,-1)
    num=((D-10)*lap*(n/(C+step))**2-(D-10)*lap*(n/(C-step))**2)/(2*step)
    exact=-2*(D-10)*lap*(n/C)**2/C
    np.testing.assert_allclose(num,exact,rtol=1e-8,atol=1e-10)
    checks['modified_Spring_diffusion_length_derivative_relative_error']=float(np.max(abs(num-exact))/np.max(abs(exact)))
    front=json.loads((META/'isolated_front.json').read_text());assert abs(front['maxwell_residual'])<1e-12 and front['quadrature_relative_change']<1e-6
    checks['independent_front_equal_area_and_quadrature_refinement']=True
    # Pilot values are rounded descriptive comparisons, never solver targets.
    records=[]
    def add(model,D,C,pilot,status='rounded_numeric'):
        z=d[(d.model==model)&(d.D==D)&(d['size']==C)].iloc[0]
        records.append(dict(model=model,D=D,size=C,our_gain=z.gain,pilot_gain=pilot,pilot_status=status,our_over_pilot_minus_one=(z.gain/pilot-1 if pilot is not None else None)))
    for D,vals in [(2.5,[1.69,1.07,None]),(5,[2.03,1.48,1.16]),(10,[2.28,1.88,1.59]),(20,[2.45,2.19,1.98]),(40,[2.54,2.4,2.27]),(100,[2.6,2.54,2.48]),(1000,[2.640,2.641,2.634])]:
        for C,val in zip([20,40,60],vals):add('wave_pinning',D,C,val,'domain_lost' if val is None else 'rounded_numeric')
    for D,vals in [(5,[.351,.207]),(12.5,[.416,.31]),(25,[.444,.376]),(50,[.46,.42]),(200,[.473,.462]),(5000,[.477,.478])]:
        for C,val in zip([20,60],vals):add('holmes_model3',D,C,val)
    for D,vals in [(1,[.276,.19]),(3,[.294,.259]),(10,[.302,.299]),(100,[.305,.318]),(1000,[.305,.32])]:
        for C,val in zip([1,3],vals):add('spring_ring',D,C,val)
    pd.DataFrame(records).to_csv(DATA/'pilot_gain_comparison.csv',index=False)
    checks['pilot_comparisons_exported_not_fitted']=len(records)
    provenance=json.loads((META/'provenance.json').read_text())
    for p,h in provenance['source_hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    checks['retained_sources_unchanged']=len(provenance['source_hashes'])
    section=json.loads((OUT/'legend.json').read_text());assert not section.get('supplement_pdf');assert all((ROOT/a['path']).is_file() and (ROOT/a['native_pdf_page']).is_file() for a in section['assets'])
    pdf=PdfReader(OUT/'figure3_round7.pdf');assert len(pdf.pages)==1;assert max(float(pdf.pages[0].mediabox.width)/72,float(pdf.pages[0].mediabox.height)/72)<=24
    svg=(OUT/'figure3_round7.svg').read_text();assert all('>'+x+'<' in svg for x in 'ABCDEF')
    assert (OUT/'scientific_readme.md').read_text().startswith('# Unexpected')
    checks['main_formats_one_page_no_supplement_and_letters']=True
    out=dict(generated_utc=campaign.source.now(),passed=True,checks=checks)
    campaign.write(META/'validation.json',out)
    provenance['output_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.rglob('*') if p.is_file() and p.suffix!='.log' and p.name!='provenance.json'}
    provenance['code_sha256'][str(Path(__file__).relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    campaign.write(META/'provenance.json',provenance);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
