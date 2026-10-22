"""Spring-only fixed-origin centroid sensitivity; intrinsic arc stays primary.

Both conventions use identical observed tenures and the primary fit cutoff.
The literal coordinate includes dilation of the fixed native material origin;
its result therefore depends on the initial angular branch of each tenure.
"""
from __future__ import annotations
import argparse
import gzip
import json
import math
from pathlib import Path
import numpy as np
from scipy.fft import rfft,irfft,next_fast_len
import figure4_round3_estimators_20261022 as E
import figure4_round3_analysis_20261022 as A


def sums(x,maxstep):
    x=np.asarray(x,float);x=x-x.mean();n=len(x);m=min(maxstep,n-1)
    numerator=np.zeros(maxstep);count=np.zeros(maxstep)
    if m<1:return numerator,count
    size=next_fast_len(2*n-1);fourier=rfft(x,size);cross=irfft(fourier*np.conj(fourier),size)
    prefix=np.r_[0.,np.cumsum(x*x)];lag=np.arange(1,m+1)
    numerator[:m]=np.maximum(prefix[n]-prefix[lag]+prefix[n-lag]-2*cross[lag],0.);count[:m]=n-lag
    # Independent direct arithmetic checks both ends of the FFT range.
    for k in sorted({1,m}):
        truth=float(np.sum((x[k:]-x[:-k])**2))
        if not np.isclose(numerator[k-1],truth,rtol=2e-8,atol=1e-18+1e-11*np.sum(x*x)):
            raise ValueError('Literal-coordinate FFT MSD failed direct arithmetic check')
    return numerator,count


def seed_sums(C,ids,maps,dt,maxstep):
    ids=E.smooth_flicker(ids,5);heading=np.array([maps[i].get(int(identity),math.nan) for i,identity in enumerate(ids)])
    ids[~np.isfinite(heading)]=-1
    result={name:np.zeros((2,maxstep)) for name in ['intrinsic_arc','literal_fixed_origin']}
    for start,stop in E.tenure_sections(ids):
        theta=np.unwrap(heading[start:stop]);circ=np.asarray(C[start:stop]);literal=theta*circ/(2*np.pi)
        arc=np.r_[0.,np.cumsum(np.diff(theta)*.5*(circ[1:]+circ[:-1])/(2*np.pi))]
        local={}
        for name,coordinate in [('intrinsic_arc',arc),('literal_fixed_origin',literal)]:
            num,count=sums(coordinate,maxstep);result[name][0]+=num;result[name][1]+=count;local[name]=np.array([num,count])
        if np.ptp(circ)<1e-10:
            np.testing.assert_allclose(local['intrinsic_arc'],local['literal_fixed_origin'],rtol=1e-7,atol=1e-12)
    return result


def fit(numerator,count,lag):
    with np.errstate(invalid='ignore',divide='ignore'):
        msd=numerator/count;valid=np.isfinite(msd);value=np.nansum(msd*lag,axis=-1)/(2*np.sum(valid*lag**2,axis=-1))
    return np.where(valid.sum(axis=-1)>=2,value,np.nan)


def histories(source,campaign,key):
    if source=='compact':
        rows=[r for r in A.read_csv(campaign/'results/compact_index.csv') if (r['model'],float(r['size']),float(r['fraction']))==key]
        for row in sorted(rows,key=lambda r:int(r['replicate'])):
            with np.load(row['compact'],allow_pickle=False) as z:
                offsets=z['region_offset'];identities=z['region_id'];head=z['region_heading']
                maps=[{int(identities[j]):float(head[j]) for j in range(int(a),int(b)) if np.isfinite(head[j])} for a,b in zip(offsets[:-1],offsets[1:])]
                yield z['circumference'].astype(float),z['dominant_track_id'].astype(int),maps,float(np.median(np.diff(z['time'])))
        return
    stage='E3' if key[1] in [1.5,2.5] else 'E2'
    manifest=json.loads((A.ROOT/f'data/runs/signature_figure4_{stage}_analysis_20261020/arms/sigma_1p0/analysis_manifest.json').read_text())
    g=next(g for g in manifest['groups'] if (g['model'],float(g['C_label']),float(g['f']))==key)
    for case in g['noise']:
        candidates=[A.ROOT/f'data/runs/{root}/workers'/case['id']/'result.json' for root in ['signature_figure4_regions_common_20261020','signature_figure4_regions_E2_20261020']]
        path=next((p for p in candidates if p.exists()),None)
        if path is None:continue
        record=json.loads(path.read_text());raw=next(s for s in record['sensitivities'] if s['threshold']==.5 and s['downsample_factor']==1 and not s['smoothed'])
        with gzip.open(raw['path'],'rt') as stream:frames=json.load(stream)['result']['frames']
        if A.digest(raw['path'])!=raw['sha256']:raise ValueError('Raw Spring track archive changed')
        maps=[{int(r['track_id']):float(r['heading']) for r in f['regions'] if r['heading'] is not None} for f in frames]
        yield np.array([f['circumference'] for f in frames]),np.array([f['dominant_track_id'] if f['dominant_track_id'] is not None else -1 for f in frames]),maps,float(frames[1]['time']-frames[0]['time'])


def run(output,source='round2',campaign=None):
    data=output/'supporting/data';primary=[r for r in A.read_csv(data/'drift_decomposition.csv') if r['model']=='spring_ring' and r['flicker_frames']=='5' and r['MSD_range_fraction']=='0.25'];rows=[]
    for row in primary:
        key=('spring_ring',float(row['C_label']),float(row['f']));maxlag=float(row.get('MSD_max_lag') or 0);parts=[];dt=None
        for C,ids,maps,cadence in histories(source,campaign,key):
            if dt is None:dt=cadence
            if not np.isclose(cadence,dt):raise ValueError('Mixed Spring cadence within condition')
            parts.append(seed_sums(C,ids,maps,dt,int(round(maxlag/dt))))
        base=dict(model=key[0],C_label=key[1],f=key[2],source=source,n_seeds=len(parts),MSD_max_lag=maxlag,source_cadence=dt,
                  primary_coordinate_convention='intrinsic_arc_integral_midpoint_C_dtheta_over_2pi',sensitivity_coordinate_convention='literal_unwrapped_theta_times_actual_C_over_2pi_fixed_native_origin')
        if not parts or maxlag<=0:rows.append(dict(base,status='unavailable_primary_MSD_fit_range'));continue
        arrays={name:np.array([p[name] for p in parts]) for name in ['intrinsic_arc','literal_fixed_origin']};n=len(parts)
        identity=('compact',key) if source=='compact' else key;weights=np.random.default_rng(E.stable_seed(identity)).multinomial(n,np.full(n,1/n),size=2000);lag=np.arange(1,arrays['intrinsic_arc'].shape[-1]+1)*dt
        estimates={};boots={}
        for name,a in arrays.items():estimates[name]=float(fit(a[:,0].sum(0),a[:,1].sum(0),lag));boots[name]=fit(weights@a[:,0],weights@a[:,1],lag)
        with np.errstate(invalid='ignore',divide='ignore'):ratio=estimates['literal_fixed_origin']/estimates['intrinsic_arc'];ratio_boot=boots['literal_fixed_origin']/boots['intrinsic_arc']
        if math.isfinite(float(row['Ds'])) and not np.isclose(estimates['intrinsic_arc'],float(row['Ds']),rtol=1e-7,atol=1e-15):raise ValueError('Sensitivity tenures do not reproduce primary intrinsic-arc diffusivity')
        lo,hi=E.percentile(ratio_boot,True);dlo,dhi=E.percentile(boots['literal_fixed_origin'],True)
        rows.append(dict(base,intrinsic_arc_Ds=estimates['intrinsic_arc'],literal_fixed_origin_Ds=estimates['literal_fixed_origin'],literal_fixed_origin_Ds_ci_lower=dlo,literal_fixed_origin_Ds_ci_upper=dhi,
                         literal_over_arc_Ds=ratio,literal_over_arc_Ds_ci_lower=lo,literal_over_arc_Ds_ci_upper=hi,bootstrap_valid_fraction=float(np.isfinite(ratio_boot).mean()),status='resolved' if np.isfinite(ratio) else 'unavailable_primary_MSD_fit'))
    E.write_csv(data/'spring_coordinate_sensitivity.csv',rows);return rows


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=A.DEFAULT_OUT);p.add_argument('--source',choices=['round2','compact'],default='round2');p.add_argument('--campaign',type=Path);args=p.parse_args()
    if args.source=='compact' and args.campaign is None:p.error('--source compact requires --campaign')
    run(args.output,args.source,args.campaign)
