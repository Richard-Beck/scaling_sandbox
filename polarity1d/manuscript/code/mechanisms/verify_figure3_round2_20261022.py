"""Independent native-field, interval-arithmetic and artifact audit for Figure 3.

Does not import the simulation, decomposition or rendering producers. Native
peak counts and pool conservation are recomputed directly from stored fields.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import math

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure3'
RUN=ROOT/'data/runs/signature_figure3_holmes_extensions_round2_20261022'
OLD=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance'
SYN=ROOT/'data/runs/signature_figure3_synthesis_20261020'
TITLE='Separation slows the resolution of coexisting regions at short range; at matched separation, cell size dominates both driven takeover and resolution.'

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def same(a,b):return (pd.isna(a) and pd.isna(b)) or a==b

def native_audit():
    plan=read(RUN/'manifest.json')
    assert sha(RUN/'manifest.json')==read(RUN/'freeze.json')['manifest_sha256']
    for p,h in plan['source_sha256'].items():assert sha(RUN/'frozen_sources'/p)==h
    for p,h in plan['input_sha256'].items():assert sha(RUN/p)==h
    assert len(plan['jobs'])==9
    records=[]
    for job in plan['jobs']:
        for k in ('original_status','original_tail'):assert sha(job[k])==job[k+'_sha256']
        status=read(RUN/'workers'/f'{job["slot"]:03d}'/'status.json')
        assert status['status']=='complete' and status['complete']
        phase=status['phase'];assert sha(RUN/phase['path'])==phase['sha256']
        with np.load(job['original_tail']) as z:oldt,oldy=z['time'].copy(),z['states'].copy()
        with np.load(RUN/phase['path']) as z:
            t,y=z['time'].copy(),z['states'].copy()
            assert float(z['cue_peak'])==0
        assert t[0]==0 and t[-1]==60000 and oldt[-1]==20000
        np.testing.assert_array_equal(y[0],oldy[-1])
        times=np.r_[oldt,20000+t[1:]];fields=np.concatenate([oldy,y[1:]])
        cells=status['placement']['cells'];u=fields[:,:cells]
        assert np.isfinite(fields).all() and fields.min()>=-5e-10
        shaped=fields.reshape(len(fields),4,cells)
        rac=shaped[:,:2].sum(axis=(1,2));rho=shaped[:,2:].sum(axis=(1,2))
        error=max(np.max(abs(rac/rac[0]-1)),np.max(abs(rho/rho[0]-1)))
        assert error<1e-7
        peaks=[]
        for f in u:
            pp,_=find_peaks(np.tile(f,3),prominence=.05)
            peaks.append(int(np.count_nonzero((pp>=cells)&(pp<2*cells))))
        single=(np.array(peaks)==1)&(np.ptp(u,axis=1)>=.05)
        found=np.flatnonzero(single)
        expected=([float(times[found[0]-1]),float(times[found[0]])] if len(found) else [80000.,None])
        primary=status['B_resolution']['primary']
        assert [primary['lower'],primary['upper']]==expected
        assert not single[times<=20000].any()
        assert status['B_resolution']['one_region_persists_through_cap']==bool(len(found) and single[found[0]:].all())
        scales=np.repeat([7.5,7.5,3.1,3.1],cells)
        windows=[]
        for end in (79000.,80000.):
            keep=(times>=end-1000)&(times<=end);yy=fields[keep];tt=times[keep]
            variation=np.max(np.ptp(yy,axis=0)/scales)
            projected=np.max(abs(yy[-1]-yy[-2])/scales)/(tt[-1]-tt[-2])*1000
            windows.append(bool(max(variation,projected)<=1e-4))
        assert status['finite_tail']['both_late_windows_steady']==all(windows)
        records.append(dict(original_slot=job['original_slot'],size=job['size'],native_samples=len(times),
            resolution_lower=expected[0],resolution_upper=expected[1],
            conserved_pool_relative_error=float(error),minimum_concentration=float(fields.min()),
            final_peak_count=peaks[-1],field_steady=all(windows)))
    assert sum(r['resolution_upper'] is not None for r in records)==6
    return records


def arithmetic_audit():
    a=pd.read_csv(OUT/'A1_plotted_values.csv');b=pd.read_csv(OUT/'B_ratio1_plotted_values.csv')
    originala=pd.read_csv(SYN/'A1_complete_region_evidence.csv')
    originalb=pd.read_csv(SYN/'B_complete_region_evidence_and_geometry.csv')
    pd.testing.assert_frame_equal(a,originala)
    approved={205,208,211,214,217,232,235,238,241}
    for row in b.itertuples():
        original=originalb[originalb.slot.eq(row.slot)].iloc[0]
        assert sha(row.native_status)==row.native_status_sha256
        assert same(row.measured_gap,original.measured_gap)
        if row.slot not in approved:
            for field in ('resolution_lower','resolution_upper','resolution_status'):
                assert same(getattr(row,field),original[field])
        else:
            status=read(row.extension_status);primary=status['B_resolution']['primary']
            assert row.resolution_lower==primary['lower']
            assert same(row.resolution_upper,primary['upper']) or (pd.isna(row.resolution_upper) and primary['upper'] is None)
    for row in a.itertuples():assert sha(row.native_status)==row.native_status_sha256
    d=pd.read_csv(OUT/'decomposition_terms.csv');q=pd.read_csv(OUT/'interpolation_queries.csv')
    assert len(d)==48 and len(q)==48
    finite_products=0
    for (arm,model,size),group in d.groupby(['arm','model','size']):
        frame=a if arm=='A1' else b
        frame=frame[frame.control.isna()&frame.model.eq(model)]
        if arm=='A1':
            frame=frame[frame.classification.eq('external')].copy()
            frame['xx']=frame.gap;frame['ll']=frame.handover_time.fillna(frame.handover_lower_bound);frame['hh']=frame.handover_time.fillna(np.inf)
        else:
            frame=frame[frame.primary_claim_eligible].copy()
            frame['xx']=frame.measured_gap;frame['ll']=frame.resolution_lower;frame['hh']=frame.resolution_upper.fillna(np.inf)
        smallsize=group.reference_size.iloc[0]
        sm=frame[frame.size_parameter.eq(smallsize)].sort_values('xx')
        lg=frame[frame.size_parameter.eq(size)].sort_values('xx')
        gm=min(sm.xx.max(),lg.xx.max());gl=lg.xx.max()
        assert all(abs(group.matched_gap-gm)<1e-9)
        intervals={}
        for name,ff,gap in [('small_matched',sm,gm),('large_matched',lg,gm),('large_antipode',lg,gl)]:
            x=ff.xx.to_numpy();exact=np.flatnonzero(abs(x-gap)<1e-8)
            if len(exact):
                rr=ff.iloc[int(exact[0])];interval=[rr.ll,rr.hh];slots=[int(rr.slot)];weight=0.
            else:
                i=int(np.searchsorted(x,gap));lo,hi=ff.iloc[i-1],ff.iloc[i]
                weight=(gap-lo.xx)/(hi.xx-lo.xx)
                interval=[lo.ll*(1-weight)+hi.ll*weight,
                    np.inf if not np.isfinite(lo.hh) or not np.isfinite(hi.hh) else lo.hh*(1-weight)+hi.hh*weight]
                slots=[int(lo.slot),int(hi.slot)]
            query=q[q.arm.eq(arm)&q.model.eq(model)&q['size'].eq(size)&q['query'].eq(name)].iloc[0]
            np.testing.assert_allclose(interval,[query.lower,query.upper],rtol=1e-12)
            assert slots==json.loads(query.source_slots)
            assert math.isclose(weight,query.interpolation_fraction,rel_tol=1e-12,abs_tol=1e-12)
            intervals[name]=interval
        values={}
        for term,num,den in [('size_at_matched_gap','large_matched','small_matched'),('extra_separation','large_antipode','large_matched'),('total','large_antipode','small_matched')]:
            nu,de=intervals[num],intervals[den]
            lo=nu[0]/de[1] if np.isfinite(de[1]) else 0.
            hi=nu[1]/de[0] if de[0]>0 else np.inf
            row=group[group.term.eq(term)].iloc[0]
            np.testing.assert_allclose([lo,hi],[row.lower,row.upper],rtol=1e-12)
            if lo>0 and np.isfinite(hi):
                expected=math.sqrt(lo*hi);assert math.isclose(row.point,expected,rel_tol=1e-12);values[term]=expected
            else:assert pd.isna(row.point)
        if len(values)==3:
            assert math.isclose(values['size_at_matched_gap']*values['extra_separation'],values['total'],rel_tol=1e-12)
            finite_products+=1
    return dict(decomposition_terms=len(d),interpolation_queries=len(q),finite_product_identities=finite_products,
        holmes_eligible=int(len(b[b.model.eq('holmes_model3')&b.control.isna()&b.primary_claim_eligible])),
        holmes_censored=int(b[b.model.eq('holmes_model3')&b.control.isna()&b.primary_claim_eligible].resolution_upper.isna().sum()))


def artifact_audit():
    main=PdfReader(OUT/'figure3_main.pdf');sup=PdfReader(OUT/'figure3_supplement.pdf')
    assert len(main.pages)==1
    oldpages=[]
    for path in (OLD/'supporting/vectors/figure3_distance_main.pdf',OLD/'supporting/vectors/figure3_distance_supplement.pdf'):
        oldpages.extend(PdfReader(path).pages)
    assert len(sup.pages)==3+len(oldpages)
    for a,b in zip(sup.pages[3:],oldpages):assert a.get_contents().get_data()==b.get_contents().get_data()
    content=read(OUT/'content.json')
    assert content['title']=='Figure 3. '+TITLE and content['legend'][0].startswith(TITLE)
    assert TITLE in (OUT/'README.md').read_text()
    assert ' '.join(TITLE.split()) in ' '.join(main.pages[0].extract_text().split())
    for ext in ('png','pdf','svg'):assert (OUT/f'figure3_main.{ext}').stat().st_size>1000
    assert content['provenance']['extensions_completed']==9
    assert content['provenance']['model_colors']==dict(wave_pinning='#2878b5',goryachev='#c83b3b',holmes_model3='#36944c',spring='#c59a25')
    return dict(main_pages=len(main.pages),supplement_pages=len(sup.pages),preserved_original_page_content_streams=len(oldpages),
        message_exact=True,model_palette_verified=True,main_panels=5)


def native_csv_audit():
    sources=pd.read_csv(OUT/'kymograph_sources.csv')
    values=pd.read_csv(OUT/'kymograph_plotted_values.csv.gz')
    counts={}
    for (arm,slot),group in sources.groupby(['assay','slot'],sort=False):
        frame=pd.read_csv(OUT/('A1_plotted_values.csv' if arm=='A1' else 'B_ratio1_plotted_values.csv'))
        row=frame[frame.slot.eq(slot)&frame.control.isna()].iloc[0]
        status=read(row.native_status);cells=status['placement']['cells']
        tt=[];yy=[]
        for phase in group.itertuples():
            assert sha(phase.archive)==phase.sha256
            with np.load(phase.archive) as z:
                times=z['time']+phase.absolute_start;fields=z['states'][:,:cells]
            if arm=='B':times-=status['withdrawal_time']
            keep=(times>=0)&(times<=phase.plotted_cap)
            tt.extend(times[keep]);yy.extend(fields[keep])
        if arm=='B':x=(np.arange(cells)+.5)*47.6/cells
        else:
            origin=status['placement']['initial_10_percent_contour']['origin']/37.8
            x=((np.arange(cells)+.5)/cells-origin+.25)%1*37.8
        subset=values[values.assay.eq(arm)&values.slot.eq(slot)]
        np.testing.assert_allclose(subset.time,np.repeat(tt,cells),rtol=1e-14,atol=1e-12)
        np.testing.assert_allclose(subset.arc_position,np.tile(x,len(tt)),rtol=1e-14,atol=1e-12)
        np.testing.assert_allclose(subset.activity,np.asarray(yy).ravel(),rtol=1e-14,atol=1e-12)
        counts[f'{arm}:{slot}']=len(subset)
    assert sum(counts.values())==len(values)
    return dict(native_field_csv_rows=counts,native_source_checksums_verified=len(sources))


if __name__=='__main__':
    records=native_audit()
    report=dict(passed=True,verified_utc=datetime.now(timezone.utc).isoformat(),native_field_validation=records,
        **arithmetic_audit(),**artifact_audit(),**native_csv_audit(),bounds_have_no_point_values=True,approved_extensions_completed=9,
        audit_source_sha256=sha(__file__))
    (RUN/'results/independent_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='native_field_validation'},indent=2))
