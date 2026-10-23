"""Audit corrected estimands, exact retained evidence and publication integration."""
import csv
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'polarity1d/manuscript/code/mechanisms'))
import figure4_round3_estimators_20261022 as E
import figure4_round4_reanalysis_20261023 as A

FIG=A.OUT;DATA=FIG/'supporting/data';META=FIG/'supporting/metadata'
DEST=ROOT/'manuscript_outline/revision_figure4_scoring_20261023'


def read(p):return json.loads(p.read_text())
def rows(p):
    with p.open(newline='') as stream:return list(csv.DictReader(stream))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(r):return r['model'],float(r['C_label']),float(r['f'])
def primary(r):return float(r['flicker_frames'])==5 and float(r['MSD_range_fraction'])==.25
def close(a,b):assert np.isclose(float(a),float(b),rtol=2e-12,atol=1e-14,equal_nan=True),(a,b)


def main():
    baseline=read(META/'baseline_round3_sources.json')
    for path,expected in baseline['source_sha256'].items():assert sha(ROOT/path)==expected,('preceding evidence changed',path)
    receipt=read(META/'scoring_revision_receipt.json')
    for path,expected in receipt['source_sha256'].items():assert sha(ROOT/path)==expected,path
    assert sha(A.CAMPAIGN/'results/compact_index.csv')==receipt['compact_index_sha256']
    source=rows(DATA/'drift_decomposition.csv');mainrows=[r for r in source if primary(r)]
    oldrows=[r for r in rows(ROOT/baseline['baseline']/'supporting/data/drift_decomposition.csv') if primary(r)]
    old={key(r):r for r in oldrows};tracked={key(r):r for r in mainrows}
    assert len(source)==1300 and len(mainrows)==len(old)==260
    for r in source:
        counted=int(r['n_switches']);scored=int(r['n_scored_switches']);exposure=float(r['valid_time'])
        assert 0<=scored<=counted
        expected=float(r['lambda_s'])*(1-float(r['mean_cos_jump'])) if scored else 0. if counted==0 and exposure>0 else math.nan
        close(r['lambda_eff'],expected)
        expected=(2*np.pi/float(r['C_mean']))**2*float(r['Ds'])+expected
        close(r['kpred'],expected)
        close(r['kpred_over_kobs'],expected/float(r['kobs']) if float(r['kobs'])>0 else math.nan)
        if not primary(r):continue
        for field in ['Ds','MSD_exponent','lambda_s','n_switches','n_scored_switches','mean_cos_jump','valid_time','MSD_max_lag','kobs','analysis_lag','n_seeds']:
            close(r[field],old[key(r)][field])
        model,size,f=key(r);ref=tracked[(model,1. if model=='spring_ring' else 20.,f)]
        den=float(ref['Ds']);ratio=float(r['Ds'])/den if den>0 else math.nan
        close(r['Ds_over_smallest'],ratio)
        with np.load(DATA/'bootstrap'/A.filename(key(r))) as z:
            for field in ['Ds','lambda_eff','kpred','kpred_over_kobs']:
                lo,hi=E.percentile(z[field],True);close(r[field+'_ci_lower'],lo);close(r[field+'_ci_upper'],hi)
            if key(r)==key(ref):lo=hi=1.
            else:
                with np.load(DATA/'bootstrap'/A.filename(key(ref))) as rz:
                    with np.errstate(invalid='ignore',divide='ignore'):lo,hi=E.percentile(z['Ds']/rz['Ds'],True)
            close(r['Ds_over_smallest_ci_lower'],lo);close(r['Ds_over_smallest_ci_upper'],hi)
    assert sum(int(r['n_switches']) for r in mainrows)==receipt['counted_replacements']==915757
    assert sum(int(r['n_scored_switches']) for r in mainrows)==receipt['scored_replacements']==321423
    assert sum(math.isfinite(float(r['lambda_eff'])) for r in mainrows)==receipt['jump_finite']==245
    fits=rows(DATA/'detector_free_decomposition.csv')
    original_fits=rows(DATA/'detector_free_mle.csv')
    assert len(fits)==len(original_fits)==520
    for r,original in zip(fits,original_fits):
        assert all(r[k]==v for k,v in original.items())
        predicted=float(r['Dtheta'])+float(r['lambda_eff']) if r['status']=='resolved' else math.nan
        close(r['kpred'],predicted)
        close(r['kpred_over_kobs'],predicted/float(r['kobs']) if float(r['kobs'])>0 else math.nan)
    alpha={key(r):float(r.get('alpha') or 'nan') for r in rows(DATA/'local_alpha.csv')}
    for audit in receipt['closure_by_regime']:
        cohort=mainrows if audit['method']=='scorable_events' else [r for r in fits if int(r['detector_free_step'])==1]
        cohort=[r for r in cohort if ('multiple_front_majority' not in audit['regime'] or float(r['multiple_front_fraction'])>.5)
                and ('alpha_below_2' not in audit['regime'] or alpha[key(r)]<2)]
        ratios=[r for r in cohort if math.isfinite(float(r['kpred_over_kobs']))]
        assert len(cohort)==audit['conditions'] and len(ratios)==audit['estimable']
        assert sum(abs(float(r['kpred_over_kobs'])-1)<=.15 for r in ratios)==audit['within_15_percent']
    for r in rows(DATA/'jump_scoring_coverage.csv'):
        lo=float(r['lambda_eff_missing_angle_lower']);hi=float(r['lambda_eff_missing_angle_upper']);value=float(r['lambda_eff'])
        assert 0<=lo<=hi
        if math.isfinite(value):assert lo-1e-14<=value<=hi+1e-14
    for name in ['A_acf','B_alpha','B_folds','C_alpha','C_folds','D_acf']:
        assert sha(DATA/f'plotted_{name}.csv')==sha(ROOT/baseline['baseline']/f'supporting/data/plotted_{name}.csv'),('A-D changed',name)
    points=rows(DATA/'plotted_E.csv')
    assert sum(r['quantity']=='kpred_over_kobs' and r['method']=='scorable_events' for r in points)==154
    assert sum(r['quantity']=='kpred_over_kobs' and r['method']=='detector_free' for r in points)==172
    assert sum(r['quantity']=='lambda_eff' and r['method']=='scorable_events' and float(r['value'])==0 for r in points)==108
    all_ratio_values=[float(r['value']) for r in points if r['quantity']=='kpred_over_kobs']
    assert min(all_ratio_values)>.003 and max(all_ratio_values)<1e4
    render=read(META/'render_receipt.json')
    for category in ['inputs','outputs','documentation']:
        for path,expected in render[category].items():assert sha(FIG/path)==expected,('render',category,path)
    for path,expected in render['code_sha256'].items():assert sha(ROOT/path)==expected,('source',path)
    integration=read(DEST/'integration.json');document=read(ROOT/integration['outline']);prior=read(ROOT/integration['baseline'])
    assert sha(ROOT/integration['outline'])==integration['outline_sha256']
    assert sha(ROOT/integration['baseline'])==integration['baseline_sha256']
    for ident,record in integration['protected_sections'].items():
        assert next(s for s in prior['sections'] if s['id']==ident)==next(s for s in document['sections'] if s['id']==ident)
        for path,expected in record['asset_sha256'].items():assert sha(ROOT/path)==expected
    assert len(PdfReader(DEST/'figure4_with_legend.pdf').pages)==2
    page=PdfReader(FIG/'figure4_main.pdf').pages[0]
    assert np.allclose([float(page.mediabox.width)/72*25.4,float(page.mediabox.height)/72*25.4],[180,500])
    assert len(PdfReader(FIG/'figure4_supplement.pdf').pages)>30
    assert '42 estimable' not in ' '.join(next(s for s in document['sections'] if s['id']=='figure4')['legend'])
    result=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,
                no_new_simulations=True,preceding_exports_unchanged=True,all_A_to_D_coordinates_identical=True,
                all_other_figure_sections_and_artwork_identical=True,checked_decomposition_rows=1300,
                checked_whole_seed_bootstrap_conditions=260,tracked_closure_ratios=154,detector_free_closure_ratios=172,
                zero_jump_terms_visible=108,all_finite_ratio_points_inside_main_axis=True,
                event_based_multiple_front_within15='0/33',detector_free_multiple_front_within15='3/33',
                render_receipt_sha256=sha(META/'render_receipt.json'),integration_sha256=sha(DEST/'integration.json'))
    write=DEST/'verification.json';write.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
