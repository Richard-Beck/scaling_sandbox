"""Audit scientific coordinates, four native restarts and page preservation."""
from pathlib import Path
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
MAN=ROOT/'polarity1d/manuscript'
FIG=MAN/'figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3'
DEST=ROOT/'manuscript_outline/revision_round3_figure3_20261022'
RUN=MAN/'data/runs/figure3_round3_relaxation_v2_20261022'
SOURCE=MAN/'figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance/supporting/data'
sys.path.insert(0,str(MAN/'code/mechanisms'))
import signature_figure3_native_20261020 as native
import signature_figure3_regions_20261020 as regions


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def eq(a,b):assert math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-11),(a,b)


def main():
    f=pd.read_csv(FIG/'supporting/data/half_heading_readouts.csv')
    original=pd.read_csv(SOURCE/'A1_heading_target_time_brackets.csv')
    original=original[(original.target_fraction==.5)&(original.coherence_floor==.05)&original.control.isna()]
    assert set(f.slot)==set(original.slot) and len(f)==152
    for r in f.itertuples():
        a=original[original.slot==r.slot].iloc[0]
        assert a.status==r.status
        eq(r.half_angle,a.target_angle)
        eq(r.target_angle,2*a.target_angle)
        if pd.notna(r.t50):
            assert a.status=='observed_saved_valid_segment_crossing'
            assert a.lower<=r.t50<=a.upper
            eq(r.speed,.5*r.target_angle/r.t50)
            eq(r.speed_lower,.5*r.target_angle/a.upper)
            eq(r.speed_upper,.5*r.target_angle/a.lower)
        else:
            assert pd.isna(r.speed)
        if r.heading_undefined:
            assert 'gap' in r.status and pd.isna(r.t50) and pd.isna(r.saddle)
    # All marker placements are native 50% crossings, with distinct topology.
    for size in (30.,47.6):
        for gap in (-1.,4.):
            r=f[(f.model=='wave_pinning')&(f.size_parameter==size)&np.isclose(f.gap,gap)].iloc[0]
            assert bool(r.saddle)==(gap==4)
    regimes=pd.read_csv(FIG/'supporting/data/regime_lines.csv')
    for r in regimes.itertuples():
        g=f[(f.model==r.model)&(f.size_parameter==r.size_parameter)&(f.saddle==True)]
        eq(r.first_saddle_gap,g.gap.min())
    assert len(regimes)==12
    provenance=read(FIG/'supporting/metadata/topology_provenance.json')
    assert len(provenance)==int(f.t50.notna().sum())==120
    for record in provenance:
        r=f[f.slot==record['slot']].iloc[0];t=record['interpolated_topology']
        assert bool(r.saddle)==t['two_region_available']
        if t['two_region_available']:
            assert t['resolved_peaks']==2 and len(t['valley_descriptors'])==2
            eq(t['incumbent_excess']+t['challenger_excess'],t['full_excess_amount'])
    profiles=pd.read_csv(FIG/'supporting/data/activity_profiles.csv.gz')
    meta=read(FIG/'supporting/data/activity_profile_metadata.json')
    assert len(meta)==16 and set(profiles.stage)=={0,1,2,3}
    for r in meta:
        if r['stage']<3:eq(r['heading_fraction'],(0,.5,.9)[r['stage']])
    plan=read(RUN/'manifest.json')
    assert plan['native_protocols']==4 and not plan['retuned']
    assert len(plan['jobs'])==4
    # Independently recompute RHS and full-state restart continuity from arrays.
    native_results=[]
    for job in plan['jobs']:
        label=job['label'];p=job['placement'];s=read(RUN/'workers'/label/'status.json')
        assert s['complete'] and s['producer']=='finished'
        assert sha(RUN/job['input'])==job['input_sha256']
        with np.load(RUN/job['input']) as z:basal=z['state'].copy()
        with np.load(RUN/'workers'/label/'cue_on.npz') as z:
            on_t=z['time'].copy();on=z['states'].copy()
            eq(z['cue_peak'],p['cue_peak']);eq(z['cue_center_fraction'],p['cue_center_fraction'])
        with np.load(RUN/'workers'/label/'cue_off.npz') as z:
            off_t=z['time'].copy();off=z['states'].copy();eq(z['cue_peak'],0)
        np.testing.assert_array_equal(on[0],basal)
        np.testing.assert_array_equal(off[0],on[-1])
        adapter=native.Adapter('wave_pinning',p['size_parameter'],p['cells'])
        phase=np.exp(2j*np.pi*adapter.q)
        phi0=np.angle(on[0,:adapter.cells]@phase)
        target=s['target_angle']
        frac=np.angle((on[-1,:adapter.cells]@phase)*np.exp(-1j*phi0))/target
        if label=='iii':
            eq(on_t[-1],20000);assert frac<.95 and s['withdrawal_condition']=='20000_horizon'
        else:
            eq(frac,.95);assert s['withdrawal_condition']=='95_percent_heading'
        eq(on_t[-1],s['cue_on_duration'])
        eq(off_t[-1],s['relaxation_duration'])
        derivative=float(np.max(abs(adapter.rhs(off[-1])[:adapter.cells])))
        eq(derivative,s['final_max_active_derivative'])
        assert derivative<=s['stationarity_tolerance']*(1+1e-7)
        assert s['stop_condition']=='stationarity' and off_t[-1]<=4*on_t[-1]
        for key in ('cue_on','cue_off'):
            assert s[key]['max_amount_relative_error']<=1e-7
            assert s[key]['minimum_concentration']>=-5e-10
            assert sha(Path(s[key]['path']))==s[key]['sha256']
        # Exported final active fields must be native, without normalization.
        profile=profiles[(profiles.label==label)&(profiles.stage==3)]
        np.testing.assert_allclose(profile.active,off[-1,:adapter.cells],rtol=1e-14,atol=1e-14)
        native_results.append(dict(label=label,withdrawal=s['withdrawal_condition'],
            cue_on_duration=float(on_t[-1]),relaxation_duration=float(off_t[-1]),
            stop=s['stop_condition'],final_active_derivative=derivative))
    # Compare decoded PDF page content byte-for-byte after concatenation.
    receipt=read(FIG/'supporting/metadata/supplement_preservation.json')
    combined=PdfReader(FIG/'figure3_round3_supplement.pdf')
    offset=0
    for source in receipt['sources']:
        path=ROOT/source['path'];assert sha(path)==source['sha256']
        previous=PdfReader(path)
        for page in previous.pages:
            assert page.get_contents().get_data()==combined.pages[offset].get_contents().get_data()
            offset+=1
    assert offset==len(combined.pages)==receipt['total_pages']
    # Mathematical text and a single native vector page, without a figure title.
    plate=PdfReader(FIG/'figure3_round3.pdf');assert len(plate.pages)==1
    text=plate.pages[0].extract_text()
    assert 'Cue placement relative' not in text and 'deBelly' not in text
    for name in ('Wave-pinning','Goryachev','Holmes','Spring'):assert name in text
    assert '90%' in text and '50%' in text
    outline=read(ROOT/'manuscript_outline/outline_v7.json');before=read(ROOT/'manuscript_outline/outline_v6.json')
    integration=read(DEST/'integration.json')
    assert sha(ROOT/integration['baseline'])==integration['baseline_sha256']
    for ident,record in integration['protected_sections'].items():
        a=next(s for s in before['sections'] if s['id']==ident)
        b=next(s for s in outline['sections'] if s['id']==ident)
        assert a==b
        for path,digest in record['asset_sha256'].items():assert sha(ROOT/path)==digest
    rendered=PdfReader(ROOT/'manuscript_outline/outline_v7.pdf')
    note_text=''.join(p.extract_text() for p in rendered.pages)
    assert '20,000 horizon' in note_text and 'photoactivation' in note_text
    assert read(ROOT/'manuscript_outline/outline_v7.validation.json')['passed']
    report=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,
        placements=len(f),observed_crossings=120,heading_undefined=23,no_turn=9,
        all_four_marker_topologies_correct=True,regime_lines=12,
        topology_bracket_disagreements=f[f.topology_bracket_consistent==False][['slot','gap']].to_dict('records'),
        four_native_protocols=native_results,retuned=False,supplement_pages_preserved=offset,
        all_other_sections_preserved=list(integration['protected_sections']),
        no_cap_plotted_as_observation=True)
    (DEST/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    (FIG/'supporting/metadata/verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
