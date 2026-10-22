"""Independent arithmetic checks for post hoc H1 and the partial integration."""
import csv,json,math
csv.field_size_limit(50_000_000)
from pathlib import Path
from datetime import datetime,timezone
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022'
DEST=ROOT/'manuscript_outline/revision_round2_20261022'
def read(p):return json.loads(p.read_text())
def rows(p):return list(csv.DictReader(p.open()))
def num(r,k):return float(r[k]) if r.get(k) else math.nan
def main():
    p=FIG/'figure4/supporting/data';calls=rows(p/'S0_posthoc_H1.csv');ratios=rows(p/'size_ratios.csv');tols=rows(p/'noise_tolerance.csv')
    assert len(calls)==140
    key=lambda r:(r['arm'],r['view'],r['forced_inclusive'],r['model'])
    assert len({key(r) for r in calls})==140
    independently_checked=0
    for call in calls:
        if call['model']=='all':continue
        arm,view,inc,model=key(call)
        ts=sorted([r for r in tols if key(r)==key(call)],key=lambda r:num(r,'C_label'))
        cutoff=num(ts[-1],'fhalf')
        rr=[r for r in ratios if (r['arm'],r['view'],r['model'])==(arm,view,model) and num(r,'f')<cutoff]
        assert json.loads(call['retained_inputs'])==[num(r,'f') for r in rr]
        if not math.isfinite(cutoff):assert call['status']=='inconclusive';continue
        eligible=[r for r in rr if r['evaluable']=='True' and (inc=='True' or r['forced_combined']=='False')]
        unknown=inc=='False' and any(r['evaluable']=='True' and r['forced_combined'] not in ('True','False') for r in rr)
        reverse=any(num(r,'rho_k_ci_upper')<1 for r in eligible)
        gains=bool(eligible) and not unknown and all(num(r,'rho_k_ci_lower')>1 for r in eligible)
        tolpass=num(ts[-1],'fhalf_ci_lower')>num(ts[0],'fhalf_ci_upper')
        tolfail=num(ts[-1],'fhalf_ci_upper')<=num(ts[0],'fhalf_ci_lower')
        status='fail' if reverse or tolfail else 'pass' if gains and tolpass else 'inconclusive'
        assert call['status']==status,(key(call),call['status'],status)
        independently_checked+=1
    for call in calls:
        if call['model']!='all':continue
        component=[r for r in calls if key(r)[:3]==key(call)[:3] and r['model']!='all']
        passed=sum(r['status']=='pass' for r in component);unknown=sum(r['status']=='inconclusive' for r in component)
        expected='pass' if passed>=3 else 'fail' if passed+unknown<3 else 'inconclusive'
        assert call['status']==expected
    maincalls=[r for r in calls if r['model']=='all' and r['view']=='main' and r['forced_inclusive']=='False']
    assert len(maincalls)==2 and all(r['status']=='pass' and r['frozen_rate_status']=='inconclusive' for r in maincalls)
    retired=rows(p/'S0_retired_gates.csv');a6=rows(p/'S0_retired_A6.csv')
    assert len(retired)==168 and len(a6)==28 and all(r['status']=='retired' for r in retired+a6)
    summary=read(FIG/'figure4/supporting/metadata/main_summary.json')
    assert summary['case']=='A' and summary['post_hoc'] and 'post hoc' in summary['legend'].lower()
    for model in ('wave_pinning','holmes_model3','spring_ring'):
        tt=sorted([r for r in tols if (r['arm'],r['view'],r['model'],r['forced_inclusive'])==('common','main',model,'False')],key=lambda r:num(r,'C_label'))
        assert num(tt[-1],'fhalf_ci_lower')>num(tt[0],'fhalf_ci_upper')
    # Geometry normalization is an outcome, not part of the H1 decision.
    for r in rows(p/'main_size_ratios.csv'):
        assert r['model']!='goryachev'
        if math.isfinite(num(r,'rho_k')):assert math.isclose(num(r,'rho_geo'),num(r,'rho_k')/num(r,'geometry_expectation'),rel_tol=1e-12)
    baseline=read(ROOT/'manuscript_outline/outline_v4.json');draft=read(DEST/'outline_v5_in_progress.json')
    for ident in ('figure1','figure5','figure6'):assert next(s for s in baseline['sections'] if s['id']==ident)==next(s for s in draft['sections'] if s['id']==ident)
    assert 'IN PROGRESS' in draft['stage']
    assert read(FIG/'figure2/legend.json')['status']=='in_progress_pending_Holmes_SLURM'
    ext=rows(FIG/'figure3/optional_extension_candidates.csv');assert len(ext)==9
    out=dict(passed=True,generated_utc=datetime.now(timezone.utc).isoformat(),posthoc_H1_rows=140,independent_model_arithmetic_checks=independently_checked,independent_overall_checks=28,retired_gate_rows=len(retired)+len(a6),H1_main=[(r['arm'],r['frozen_rate_status'],r['status']) for r in maincalls],protected_sections=['figure1','figure5','figure6'],optional_extension_count=9,formal_outline_v5_integration_deferred=True,figure4_supplement_pages=len(PdfReader(FIG/'figure4/figure4_supplement.pdf').pages))
    (DEST/'followup_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
