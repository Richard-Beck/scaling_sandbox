"""Freeze deterministic Figure 4 D examples as portable compact records.

This preparation step reads immutable round-2 native fields once. The normal
figure rebuild subsequently reads only these small portable NPZ exports.
"""
from __future__ import annotations

import gzip
import json
import math
from pathlib import Path
import numpy as np

import figure4_round3_estimators_20261022 as E
from figure4_round3_render_20261022 import (OUT,ROOT,MODELS,NAMES,n,finite,
    read_csv,write_csv,write_json,sha,equivalent_size,primary)


def prepare(out=OUT,campaign=None):
    data,meta=out/'supporting/data',out/'supporting/metadata'
    lag=primary(read_csv(data/'condition_lag_metrics.csv'))
    tols=primary(read_csv(data/'noise_tolerance.csv'))
    analysis=primary(read_csv(data/'analysis_lag_metrics.csv'))
    manifest_path=Path(campaign)/'manifest.json' if campaign else ROOT/'data/runs/signature_figure4_E2_analysis_20261020/arms/sigma_1p0/analysis_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    compact_index=read_csv(Path(campaign)/'results/compact_index.csv') if campaign else []
    compact_events=read_csv(data/'seed0_switch_events.csv') if campaign else []
    source='completed_compact_campaign' if campaign else 'provisional_round2'
    index=[];selection=[];sources=[]
    for model in MODELS:
        small,large=((1.,3.) if model=='spring_ring' else (20.,60.))
        ref=200. if model=='spring_ring' else 2000.
        endpoint=[r for r in lag if r['model']==model and n(r,'C_label') in (small,large) and math.isclose(n(r,'lag'),ref)]
        fs=sorted({n(r,'f') for r in endpoint})
        halves={c:n(next((r for r in tols if r['model']==model and math.isclose(n(r,'C_label'),c) and math.isclose(n(r,'reference_lag'),ref)),{}),'fhalf') for c in (small,large)}
        midpoint=math.sqrt(halves[small]*halves[large]) if all(finite(v) and v>0 for v in halves.values()) else math.nan
        qualified=[]
        for f in fs:
            a=next((r for r in endpoint if n(r,'f')==f and n(r,'C_label')==small),{})
            b=next((r for r in endpoint if n(r,'f')==f and n(r,'C_label')==large),{})
            if n(a,'C')<.4 and n(b,'C')>.6:qualified.append(f)
        candidates=qualified or fs
        if not candidates:
            selection.append(dict(model=model,status='unavailable_no_ref_correlations'))
            continue
        if not qualified and not finite(midpoint):
            selection.append(dict(model=model,status='unavailable_no_qualifying_input_or_finite_tolerance_midpoint',selection_reason='No designed input satisfies the reference-correlation inequalities and an endpoint noise tolerance is unresolved; no input is substituted'))
            continue
        chosen=min(candidates,key=lambda f:abs(math.log(f/midpoint))) if finite(midpoint) else min(qualified)
        reason=('C_small(reference)<0.4 and C_large(reference)>0.6; nearest log tolerance midpoint among qualifying designed inputs' if finite(midpoint) else 'C_small(reference)<0.4 and C_large(reference)>0.6; lowest qualifying designed input, endpoint tolerance midpoint unavailable') if qualified else 'No designed input satisfies the two reference-correlation inequalities; nearest designed input to the log midpoint of endpoint noise tolerances'
        a=next((r for r in analysis if r['model']==model and n(r,'C_label')==small and math.isclose(n(r,'f'),chosen)),{})
        k=n(a,'k');requested=6/k if k>0 else math.nan
        horizon=3200. if model=='spring_ring' else 32000.
        if campaign:
            pair_rows=[r for r in compact_index if r['model']==model and n(r,'size') in (small,large) and math.isclose(n(r,'fraction'),chosen) and int(r['replicate'])==0 and r['status'].lower()=='complete']
            if pair_rows:horizon=min(n(r,'horizon',default=horizon) for r in pair_rows)
        start=horizon*.5;duration=min(requested,horizon*.5) if finite(requested) else horizon*.5
        stop=start+duration
        selection.append(dict(model=model,f=chosen,reference_lag=ref,selection_reason=reason,window_requested=requested,window_actual=duration,window_start=start,window_stop=stop,window_clipped=not finite(requested) or requested>duration,source=source,seed_policy='declared replicate0; no seed search' if campaign else 'lowest numeric existing seed; first declared replicate, no seed search',small_fhalf=halves[small],large_fhalf=halves[large],kobs_small=k,analysis_lag=n(a,'analysis_lag')))
        for c in (small,large):
            if campaign:
                row=next(r for r in compact_index if r['model']==model and math.isclose(float(r['size']),c) and math.isclose(float(r['fraction']),chosen) and int(r['replicate'])==0 and r['status'].lower()=='complete')
                case=dict(profile_path=row['exemplar'],seed=row['seed'],profile_sha256=sha(Path(row['exemplar'])),id=row['condition_id'])
            else:
                group=next(g for g in manifest['groups'] if g['model']==model and math.isclose(float(g['C_label']),c) and math.isclose(float(g['f']),chosen))
                case=min(group['noise'],key=lambda r:int(r['seed']))
            path=Path(case['profile_path'])
            if sha(path)!=case['profile_sha256']:raise ValueError('Native example hash changed: '+str(path))
            with np.load(path,allow_pickle=False) as archive:
                fulltime=archive['time']
                record_start=float(fulltime[0]) if campaign else start
                record_stop=record_start+duration
                keep=(fulltime>=record_start)&(fulltime<=record_stop)
                t=fulltime[keep];field=archive['primary'][keep].astype(np.float32)
                circ=np.broadcast_to(archive['circumference'],fulltime.shape)[keep].astype(np.float32)
                positions=archive['material_positions'].astype(float) if 'material_positions' in archive else (np.arange(field.shape[1])+.5)/field.shape[1]
                space_edges=np.r_[0,(positions[:-1]+positions[1:])/2,1]*float(circ.mean())
            if len(t)>2000:
                take=np.unique(np.rint(np.linspace(0,len(t)-1,2000)).astype(int));t=t[take];field=field[take];circ=circ[take]
            if field.shape[1]>400:
                # Physical-bin block averages, never upsampling or spatial interpolation.
                bins=np.array_split(np.arange(field.shape[1]),400)
                space_edges=space_edges[np.r_[0,[ix[-1]+1 for ix in bins]]]
                field=np.column_stack([field[:,ix].mean(1) for ix in bins]).astype(np.float32)
            switch_times=[];track_source=None
            roots=[ROOT/'data/runs/signature_figure4_regions_common_20261020',ROOT/'data/runs/signature_figure4_regions_E2_20261020']
            record_path=next((r/'workers'/case['id']/'result.json' for r in roots if (r/'workers'/case['id']/'result.json').exists()),None) if not campaign else None
            if campaign:
                switch_times=[n(r,'time') for r in compact_events if r['model']==model and n(r,'C_label')==c and math.isclose(n(r,'f'),chosen) and record_start<=n(r,'time')<=record_stop]
                track_source=Path(row['compact'])
            if record_path:
                record=json.loads(record_path.read_text())
                raw=next((s for s in record['sensitivities'] if s['threshold']==.5 and s['downsample_factor']==1 and not s['smoothed']),None)
                if raw:
                    track_source=Path(raw['path'])
                    with gzip.open(track_source,'rt') as stream:frames=json.load(stream)['result']['frames']
                    ids=E.smooth_flicker([f['dominant_track_id'] if f['dominant_track_id'] is not None else -1 for f in frames],5)
                    # Missing accepted-region centroids split tenures, exactly as the estimator.
                    for i,f in enumerate(frames):
                        headings={int(r['track_id']):r['heading'] for r in f['regions'] if r['heading'] is not None}
                        if int(ids[i]) not in headings:ids[i]=-1
                    sections=E.tenure_sections(ids)
                    for (aa,bb),(cc,dd) in zip(sections[:-1],sections[1:]):
                        when=frames[cc]['time']
                        if bb==cc and ids[aa]!=ids[cc] and start<=when<=stop:switch_times.append(when)
            name=f'model_D_{model}_C{c:g}.npz'
            np.savez_compressed(data/name,time=t.astype(np.float32),field=field,circumference=circ,space_edges=space_edges.astype(np.float32))
            index.append(dict(file=name,model=model,C_label=c,C_mean=float(circ.mean()),f=chosen,seed=case['seed'],replicate_index=0,source=source,native_window_start=float(t[0]),native_window_stop=float(t[-1]),switch_times=json.dumps(switch_times),switch_status='counted_from_compact_tracks' if campaign else 'counted_from_retained_raw_tracks' if track_source else 'unavailable_no_retained_raw_tracks',native_profile_path=str(path),native_profile_sha256=case['profile_sha256'],native_track_path=str(track_source) if track_source else '',native_track_sha256=sha(track_source) if track_source else '',shape=json.dumps(list(field.shape))))
            sources.append(dict(model=model,C_label=c,f=chosen,seed=case['seed'],profile_path=str(path),profile_sha256=case['profile_sha256'],output=name,output_sha256=sha(data/name)))
    write_csv(data/'model_D_kymographs.csv',index)
    write_csv(data/'model_D_selection.csv',selection)
    write_json(meta/'model_D_provenance.json',dict(source=source,stage_manifest=str(manifest_path),stage_manifest_sha256=sha(manifest_path),selection=selection,sources=sources,seed_search=False))
    return selection


if __name__=='__main__':prepare()
