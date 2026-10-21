import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_figure4_site_20261020 as S
import signature_figure4_site_audit_20261020 as Q
import signature_figure4_analysis_20261020 as A
import signature_figure4_e0_20261020 as E


def test_literal_position_retains_uniform_dilation_at_fixed_material_heading():
    C=np.array([2.,3.,4.]);M=np.full(3,1j);valid=np.ones(3,bool)
    np.testing.assert_allclose(S.literal_coordinate(M,C,valid),C/4)
    np.testing.assert_allclose(Q.coordinate(M,C,valid),C/4)
    np.testing.assert_allclose(A.physical_coordinate(M,C,valid),0.)


def test_literal_unwrap_keeps_material_origin_and_resets_branch_after_invalid_gap():
    phase=np.array([1.9,2.1,0.,.3,.5])*np.pi;valid=np.array([True,True,False,True,True]);C=np.full(5,2*np.pi)
    x=S.literal_coordinate(np.exp(1j*phase),C,valid)
    np.testing.assert_allclose(x[valid],phase[valid]);assert np.isnan(x[2])
    np.testing.assert_allclose(x,Q.coordinate(np.exp(1j*phase),C,valid),equal_nan=True)
    _,counts,_=A.section_msd(x,valid,3);np.testing.assert_array_equal(counts,[2,0,0])


def test_literal_map_uses_physical_moment_instead_of_rotated_normalized_P(tmp_path):
    row=case(tmp_path,0);tracks=S.native_tracks(row);t,M,C,valid=tracks[S.METHODS[0]]
    assert np.angle(M[0])==pytest.approx(.6)
    assert S.literal_coordinate(M,C,valid)[0]==pytest.approx(.6*C[0]/(2*np.pi))


def case(folder,seed):
    p=folder/f'case{seed}';p.mkdir();t=np.arange(1601.);phase=.6+.001*t+seed*.03;M=np.exp(1j*phase);M[210:215]=0;P=M/M[0];C=2.5+.1*np.sin(t/100)
    np.savez(p/'scalar.npz',time=t,P_real=P.real,P_imag=P.imag,M1_real=M.real,M1_imag=M.imag,circumference=C,contrast=np.ones(len(t)),current_half_contrast_regions=np.ones(len(t)))
    pt=t[::4];q=np.arange(16)/16;field=1+np.cos(2*np.pi*q[None,:]-phase[::4,None]);np.savez(p/'profile.npz',time=pt,circumference=C[::4],material_positions=q,primary=field)
    np.savez(p/'basal.npz',primary=field[0]);(p/'status.json').write_text('{"complete":true}')
    r=dict(id=f'test{seed}',model='spring_ring',C_label=2.5,f=.03,kind='noise',seed=seed,scalar_cadence=1.,profile_cadence=4.,native_horizon=1600.,basal_contrast=2.)
    for prefix,name in [('scalar','scalar.npz'),('profile','profile.npz'),('basal_profile','basal.npz'),('source_status','status.json')]:r[prefix+'_path']=str(p/name);r[prefix+'_sha256']=Q.digest(p/name)
    return r


def fixture_root(tmp_path):
    root=tmp_path/'literal';root.mkdir();(root/'workers').mkdir();(root/'audits').mkdir();(root/'results').mkdir()
    cases=[case(tmp_path,i) for i in range(2)]
    group=dict(task=0,model='spring_ring',C_label=2.5,f=.03,sigma=1.,source_root='synthetic',source_stage='test',native_group_sha256='analytic',noise=cases)
    E.write_json(root/'manifest.json',dict(groups=[group],methods=S.METHODS,views=E.VIEWS,antipodal_policies=S.POLICIES,antipodal_atol_radians=S.ANTIPODAL_ATOL,input_receipts={},source_sha256={str(Path(module.__file__).resolve()):Q.digest(module.__file__) for module in (S,Q,A,E)}))
    return root


def test_literal_worker_independent_native_audit_and_compressed_export(tmp_path):
    root=fixture_root(tmp_path);assert S.worker(root,0)['passed'];assert Q.audit(root,0)['passed'];assert S.collect(root)['complete']
    status=E.read_json(root/'results/status.json');assert status['outputs']['site_curves.csv.gz']==Q.digest(root/'results/site_curves.csv.gz')
    import gzip,csv
    with gzip.open(root/'results/site_curves.csv.gz','rt') as stream:rows=list(csv.DictReader(stream))
    assert {r['method'] for r in rows}==set(S.METHODS)
    assert {r['view'] for r in rows}==set(E.VIEWS)
    assert {r['antipodal_policy'] for r in rows}==set(S.POLICIES)
    assert all(r['primary_A3']=='True' for r in rows if r['method']=='heading_literal_actualC' and r['antipodal_policy']=='split')
    assert all(r['primary_A3']=='False' for r in rows if r['antipodal_policy']!='split')


def test_independent_literal_audit_rejects_changed_interior_MSD(tmp_path):
    root=fixture_root(tmp_path);S.worker(root,0);p=root/'workers/0000/arrays.npz'
    with np.load(p) as z:values={k:z[k] for k in z.files}
    values['main__heading_literal_actualC__split__MSD'][4]+=.1;np.savez_compressed(p,**values)
    result=E.read_json(p.parent/'result.json');result['arrays_sha256']=Q.digest(p);E.write_json(p.parent/'result.json',result)
    with pytest.raises(AssertionError,match='Direct native coordinate MSD'):Q.audit(root,0)


def test_invalid_circumference_does_not_create_a_literal_coordinate():
    with pytest.raises(ValueError):S.literal_coordinate(np.ones(2,complex),np.array([1.,np.nan]),np.ones(2,bool))


def test_antipodal_split_preserves_both_frames_but_removes_route_ambiguous_pairs():
    phase=np.pi*np.array([.2,1.2,1.3]);signal=np.exp(1j*phase);C=np.array([3.,4.,5.]);valid=np.ones(3,bool)
    x,cuts,ties=S.route_data(signal,C,valid,'split')
    assert np.isfinite(x).all() and valid.sum()==3
    np.testing.assert_array_equal(ties,[False,True,False])
    n,c,p=S.split_msd(x,valid,cuts,2);np.testing.assert_array_equal(c,[1,0]);np.testing.assert_array_equal(p,[2,1])
    q,qcuts,qties=Q.route(signal,C,valid,'split');np.testing.assert_allclose(x,q);np.testing.assert_array_equal(cuts,qcuts);np.testing.assert_array_equal(ties,qties)
    dn,dc=Q.direct(q,valid,[1,2],qcuts);np.testing.assert_allclose(n,dn);np.testing.assert_array_equal(c,dc)
    base=S.pair_counts(valid,np.zeros(3,bool),2);np.testing.assert_array_equal(base-c,[1,1])


def test_antipodal_routes_are_paired_sensitivities_with_independent_audit_agreement():
    signal=np.exp(2j*np.pi*np.array([.4,.9,.0,.5,.6]));C=np.linspace(3,4,5);valid=np.ones(5,bool)
    curves=[]
    for policy in S.POLICIES:
        x,cuts,ties=S.route_data(signal,C,valid,policy);q,qcuts,qties=Q.route(signal,C,valid,policy)
        np.testing.assert_allclose(x,q);np.testing.assert_array_equal(cuts,qcuts);np.testing.assert_array_equal(ties,qties)
        n,c,_=S.split_msd(x,valid,cuts,4);dn,dc=Q.direct(q,valid,np.arange(1,5),qcuts);np.testing.assert_allclose(n,dn,rtol=1e-10,atol=1e-10);np.testing.assert_array_equal(c,dc)
        curves.append(n)
    assert not np.allclose(curves[1],curves[2])
