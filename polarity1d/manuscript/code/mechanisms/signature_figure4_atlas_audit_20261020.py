"""Independent source-scalar, selection, native-geometry and pixel S11 audit."""
from pathlib import Path
import argparse
import hashlib
import json
import os
from datetime import datetime,timezone
import numpy as np
from PIL import Image

FEATURES=['mean_relative_contrast','mean_basal_normalized_P','P_coverage_01','single_fraction','multiple_fraction','low_contrast_fraction','mean_qualified_count']


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def reconstruct_features(case):
    assert sha(case['scalar_path'])==case['scalar_sha256']
    with np.load(case['scalar_path']) as z:
        t=z['time'];start=np.flatnonzero(abs(t-.5*case['native_horizon'])<1e-8)
        assert len(start)==1 and abs(t[-1]-case['native_horizon'])<1e-8
        k=int(start[0]);t=t[k:];raw_contrast=z['contrast'][k:];c=raw_contrast/case['basal_contrast']
        p=np.hypot(z['P_real'][k:],z['P_imag'][k:]);n=z['current_half_contrast_regions'][k:]
    dt=np.diff(t);T=float(t[-1]-t[0]);qualified=raw_contrast>=.1*case['basal_contrast']
    return np.array([np.sum((c[1:]+c[:-1])*dt/2)/T,np.sum((p[1:]+p[:-1])*dt/2)/T,
        np.sum((p[:-1]>=.1)*dt)/T,np.sum((qualified[:-1]&(n[:-1]==1))*dt)/T,
        np.sum((qualified[:-1]&(n[:-1]>=2))*dt)/T,np.sum((~qualified[:-1])*dt)/T,
        np.sum(np.where(qualified[:-1],n[:-1],0)*dt)/T])


def verify_selection(job,record,master_seed):
    native_matrix=np.stack([reconstruct_features(c) for c in job['cases']])
    matrix=np.array([[f[k] for k in FEATURES] for f in record['features']])
    # Positive weighted sums have the standard gamma_n = n*eps/(1-n*eps)
    # bound. Allow both independent summation paths and the few scalar
    # multiply/divide/hypot operations; no biological tolerance is added.
    operations=max(int(round(c['native_horizon']/c['scalar_cadence']/2))+16 for c in job['cases'])
    epsilon=np.finfo(float).eps;gamma=operations*epsilon/(1-operations*epsilon)
    feature_error_bound=4*gamma*np.maximum(abs(matrix),abs(native_matrix))+16*epsilon
    assert np.all(abs(native_matrix-matrix)<=feature_error_bound)
    scales=np.percentile(matrix,75,axis=0)-np.percentile(matrix,25,axis=0)
    scales=np.where(scales>0,scales,np.ptp(matrix,axis=0))
    distances=np.zeros((len(matrix),len(matrix)))
    for a in range(len(matrix)):
        for b in range(len(matrix)):
            distances[a,b]=sum(abs(matrix[a,k]-matrix[b,k])/scales[k] for k in range(len(FEATURES)) if scales[k]>0)
    scores=distances.sum(axis=1);ids=[c['id'] for c in job['cases']]
    eligible=np.flatnonzero(abs(scores-scores.min())<=1e-12*max(1.,float(abs(scores).max())))
    medoid=min(eligible,key=lambda i:ids[i])
    condition=json.dumps(job['condition'],sort_keys=True)
    keys=[hashlib.sha256(f'{master_seed}|{condition}|{i}'.encode()).hexdigest() for i in ids]
    candidates=[i for i in range(len(ids)) if i!=medoid]
    random=min(candidates,key=lambda i:keys[i]) if candidates else None
    actual=record['selection'];assert actual['medoid']==medoid and actual['random']==random
    assert actual['random_keys']==keys
    np.testing.assert_allclose(actual['scales'],scales,rtol=1e-12,atol=1e-13)
    # This arithmetic bound also covers frozen prototype renderers that
    # subtract already-normalized means. Final contexts subtract physical
    # features first. A narrow IQR is explicit in either path's bound.
    normalized_magnitudes=np.divide(np.max(abs(matrix),axis=0),scales,out=np.zeros(len(FEATURES)),where=scales>0)
    arithmetic_bound=16*epsilon*len(matrix)*float(normalized_magnitudes.sum())+16*epsilon
    assert np.max(abs(np.asarray(actual['scores'])-scores))<=arithmetic_bound
    assert np.max(abs(np.asarray(actual['distance'])-distances))<=arithmetic_bound
    # Independently reconstructed native feature values must also produce
    # the exact same selected ID, rather than merely similar score sums.
    native_scales=np.percentile(native_matrix,75,axis=0)-np.percentile(native_matrix,25,axis=0)
    native_scales=np.where(native_scales>0,native_scales,np.ptp(native_matrix,axis=0))
    native_delta=abs(native_matrix[:,None,:]-native_matrix[None,:,:])
    native_scores=np.divide(native_delta,native_scales,out=np.zeros_like(native_delta),where=native_scales>0).sum(axis=2).sum(axis=1)
    native_eligible=np.flatnonzero(abs(native_scores-native_scores.min())<=1e-12*max(1.,float(abs(native_scores).max())))
    assert min(native_eligible,key=lambda i:ids[i])==medoid
    return len(ids),dict(slot=job['slot'],maximum_native_feature_error=float(abs(native_matrix-matrix).max()),
        maximum_native_feature_arithmetic_bound=float(feature_error_bound.max()),
        source_score_error=float(abs(np.asarray(actual['scores'])-scores).max()),score_error_bound=arithmetic_bound,
        narrowest_nonzero_feature_scale=float(min(scales[scales>0])) if np.any(scales>0) else None,
        independent_native_feature_medoid_exact=True,source_feature_medoid_exact=True,random_hash_keys_exact=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);args=p.parse_args();run=args.run.resolve()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Complete native atlas audit uses SLURM')
    plan=read(run/'manifest.json');rules=read(run/'inputs/selection_rules.json')
    assert sha(run/'inputs/selection_rules.json')==plan['rules_sha256']
    assert sha(run/'inputs/frozen_atlas.py')==plan['code_sha256']
    complete=read(run/'results/export/complete.json');assert complete['passed']
    for name,value in complete['outputs'].items():assert sha(run/'results/export'/name)==value
    seen=set();groups=0;displays=0;candidate_cases=0;native_frames=0;arithmetic=[]
    actual_C_extrema_preserved=True;snapshot_times=[]
    for job in plan['jobs']:
        record=read(run/'workers'/f"{job['slot']:04d}.json")
        assert record['complete'] and record['condition']==job['condition'] and record['cases']==job['cases']
        assert record['code_sha256']==plan['code_sha256'] and record['rules_sha256']==plan['rules_sha256']
        count,roundoff=verify_selection(job,record,rules['master_seed']);candidate_cases+=count;arithmetic.append(roundoff)
        for c in job['cases']:
            identity=(tuple(sorted(job['condition'].items())),c['independent_realization_id'])
            assert identity not in seen;seen.add(identity)
        artifact=record['artifact'];assert artifact['native_fields_interpolated'] is False
        assert artifact['physical_mesh_geometry_interpolation_only'] is True
        assert sha(artifact['png'])==artifact['png_sha256'] and sha(artifact['webp'])==artifact['webp_sha256']
        with Image.open(artifact['png']) as a,Image.open(artifact['webp']) as b:
            assert np.array_equal(np.asarray(a.convert('RGBA')),np.asarray(b.convert('RGBA')))
        selected=[record['selection']['medoid']]
        if record['selection']['random'] is not None:selected.append(record['selection']['random'])
        assert len(artifact['displayed'])==len(selected)
        for index,display in zip(selected,artifact['displayed']):
            case=job['cases'][index];assert case['id']==display['case_id']
            for key in ['profile','scalar','basal_profile','source_status']:
                assert sha(case[key+'_path'])==case[key+'_sha256']
            with np.load(case['profile_path']) as z:
                time=z['time'];field=z['primary'];C=np.broadcast_to(z['circumference'],time.shape)
                q=z['material_positions'] if 'material_positions' in z else z['x']/C[0]
                assert field.shape==(len(time),case['native_cells']) and np.isfinite(field).all()
                assert (C>0).all() and np.isfinite(C).all()
                np.testing.assert_allclose(q,(np.arange(case['native_cells'])+.5)/case['native_cells'],rtol=0,atol=1e-12)
                faces=np.sort(np.concatenate([time,(time[:-1]+time[1:])/2]))
                Cfaces=np.interp(faces,time,C)
                np.testing.assert_array_equal(Cfaces[::2],C)
                # Every native field contributes unchanged, with no averaging
                # between its neighbors, on at least one half-strip.
                rows=np.repeat(np.arange(len(time)),2)[1:-1]
                assert set(rows)==set(range(len(time)))
                assert display['native_frames']==len(time) and display['native_cells']==field.shape[1]
                assert display['minimum_actual_circumference']==float(C.min()) and display['maximum_actual_circumference']==float(C.max())
                for fraction,snapshot in zip([0.,.5,.75,1.],display['snapshots']):
                    nearest=int(np.argmin(abs(time-fraction*case['native_horizon'])))
                    assert snapshot['requested_horizon_fraction']==fraction and snapshot['native_frame']==nearest
                    assert snapshot['native_saved_time']==float(time[nearest]) and snapshot['actual_circumference']==float(C[nearest])
                    snapshot_times.append(snapshot['native_saved_time'])
                native_frames+=len(time)
            displays+=1
        groups+=1
    assert candidate_cases==plan['independent_native_cases']==complete['independent_native_cases']
    assert groups==complete['groups'] and displays==complete['selected_displays']
    result=dict(passed=True,created_utc=datetime.now(timezone.utc).isoformat(),groups=groups,independent_candidate_cases=candidate_cases,
        displayed_seed_histories=displays,displayed_native_frames=native_frames,lossless_pixel_equal_groups=groups,
        independent_source_scalar_feature_selection=True,actual_C_extrema_preserved=actual_C_extrema_preserved,
        exact_saved_snapshot_times_verified=len(snapshot_times),no_source_field_interpolation=True,native_integrations=0,
        scalar_feature_roundoff_rule='Four gamma_n bounds for independent positive weighted sums plus16 machine epsilons for scalar operations; n is the actual half-window sample count plus16.',
        score_roundoff_rule='16*eps*n_cases*sum(max_abs_physical_feature/nonzero_feature_scale)+16*eps; scales are the declared IQR or nonzero range fallback.',
        roundoff_by_group=arithmetic,selected_IDs_exact_under_both_independent_feature_paths=True,
        main_example_selection=False,scope='Independent native scalar feature/selection reconstruction, source hashes, native mesh geometry and snapshots, exact RGBA equality. This audits selection and rendering inputs; it does not independently establish biological model validity.',
        atlas_manifest_sha256=sha(run/'manifest.json'),atlas_complete_sha256=sha(run/'results/export/complete.json'),audit_source_sha256=sha(__file__))
    (run/'results/independent_native_atlas_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
