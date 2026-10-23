"""Lightweight derivative checks; no native simulation or repeated covariance solve."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,re
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
from pypdf import PdfReader
from matplotlib.mathtext import MathTextParser
import figure3_round5_review_20261023 as revision
base=revision.base;OUT=revision.OUT;DATA=revision.DATA;META=revision.META;ROOT=revision.ROOT

def main():
    checks={};history=base.read(META/'reviewed_49754ab_preservation.json')
    for row in history['reviewed_main_and_supplement']+history['unchanged_scientific_data']:
        path=ROOT/row['path'];assert base.edge.sha(path)==row['sha256'],path
        previous=subprocess.run(['git','show','49754ab:'+row['path']],cwd=ROOT,capture_output=True,check=True).stdout
        assert hashlib.sha256(previous).hexdigest()==row['sha256'],row['path']
    checks['reviewed_artifacts_and_scientific_data_match_commit49754ab']=True
    source=pd.read_csv(DATA/'edge_prediction_comparison.csv');normal=pd.read_csv(DATA/'size_normalized_prediction_comparison.csv');areas=pd.read_csv(DATA/'normalized_integrated_edge_response.csv')
    assert len(source)==len(normal)==120;assert normal.defined.sum()==115
    undefined=normal[~normal.defined];assert (undefined.model=='goryachev').all() and (undefined.fraction==.01).all() and len(undefined)==5
    for r in normal.itertuples():
        q=source[(source.model==r.model)&(source.fraction==r.fraction)];reference=q.sort_values('size').iloc[0];s=q[q['size']==r.size].iloc[0]
        assert r.reference_size==reference['size']
        np.testing.assert_allclose(r.observed_ratio,s.observed_Ds/reference.observed_Ds,rtol=1e-12)
        if r.defined:np.testing.assert_allclose(r.shared_translation_predicted_ratio,s.Ds_predicted_full_covariance/reference.Ds_predicted_full_covariance,rtol=1e-12)
        else:assert np.isnan(r.shared_translation_predicted_ratio)
    checks['same_input_fixed_smallest_size_normalization']=True
    checks['undefined_Goryachev_reference_not_substituted']=True
    for model in revision.MODELS:
        q=areas[areas.model==model].sort_values('size');np.testing.assert_allclose(q.integrated_response_over_smallest,q.integrated_displacement_per_pulse_dose/q.iloc[0].integrated_displacement_per_pulse_dose,rtol=1e-12)
    checks['C_inset_normalization']=True
    recipe=base.read(META/'review_render_recipe.json')
    assert recipe['C_inset']['y_limits']==[0,1.1]
    assert recipe['D_absolute']['filled']=='Ds_predicted_full_covariance' and recipe['D_absolute']['open']=='Ds_predicted_independent'
    assert recipe['D_absolute']['parameter_fits']==0
    checks['shared_inset_axes_both_absolute_candidates_no_fit']=True
    q=normal[(normal.model=='wave_pinning')&(normal['size']==60)&(normal.fraction==.001)].iloc[0]
    np.testing.assert_allclose(q.shared_translation_predicted_ratio,.4939004889085995,rtol=1e-12);np.testing.assert_allclose(q.observed_ratio,.47324320218752636,rtol=1e-12)
    checks['WP_C60_C20_f001']=dict(predicted=float(q.shared_translation_predicted_ratio),observed=float(q.observed_ratio))
    section=base.read(OUT/'legend.json');assert len(section['equations'])==2
    assert [x['number'] for x in section['equations']]==['3.1','3.2']
    parser=MathTextParser('path')
    for equation in section['equations']:parser.parse('$'+equation['math']+'$',dpi=100)
    legend=' '.join(section['legend']);assert 'boundary)..' not in legend
    for phrase in ('rejected independent-edge','0.49 predicted and 0.47 measured','different size grids','vary photoactivation offset'):
        assert phrase in legend,phrase
    assert len(legend.split())<650
    checks['concise_scientific_legend_two_assays_and_scope']=dict(words=len(legend.split()),equations=2)
    def paths(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k in ('path','native_pdf_page','main_pdf','supplement_pdf','reader'):assert (ROOT/v).is_file(),v
                else:paths(v)
        elif isinstance(x,list):
            for value in x:paths(value)
    paths(section)
    pdf=PdfReader(ROOT/section['main_pdf']);assert len(pdf.pages)==1
    assert abs(float(pdf.pages[0].mediabox.width)/72-16.2)<1e-8 and abs(float(pdf.pages[0].mediabox.height)/72-15.5)<1e-8
    preservation=base.read(META/'review_supplement_preservation.json');assert len(PdfReader(ROOT/section['supplement_pdf']).pages)==preservation['total_pages']
    for row in preservation['preserved']:assert base.edge.sha(ROOT/row['path'])==row['sha256']
    checks['PDF_dimensions_pages_and_local_sources']=True
    previous_checks=base.read(META/'numerical_checks.json');assert previous_checks['passed']
    checks['prior_numerical_validation_retained_no_repeated_native_work']=True
    artifacts=[OUT/f'{revision.STEM}.{ext}' for ext in ('pdf','png','svg')]+[OUT/f'{revision.STEM}_supplement.pdf',OUT/'legend.json',OUT/'legend.md',OUT/'METHODS_RESULTS_review.md',OUT/'index.html',DATA/'size_normalized_prediction_comparison.csv',DATA/'normalized_integrated_edge_response.csv']
    receipt=dict(generated_utc=base.edge.now(),passed=True,checks=checks,artifacts=[dict(path=base.rel(p),sha256=base.edge.sha(p)) for p in artifacts],native_reruns=0,retuned_parameters=0)
    base.write(META/'review_validation.json',receipt);print(json.dumps(receipt,indent=2));return receipt
if __name__=='__main__':main()
