from pathlib import Path
import sys
import math
import json
import numpy as np
import pytest

CODE=Path(__file__).resolve().parents[1]/'mechanisms'
sys.path.insert(0,str(CODE))
import signature_figure4_stage_analysis_20261020 as S
import signature_figure4_e0_20261020 as E
import signature_figure4_audit_20261020 as Q


def row(sigma=1.,seed=1,kind='noise',origin='E1',horizon=200.):
    return dict(id=f'{kind}_{sigma}_{seed}',model='m',C_label=20.,f=.02 if kind!='zero' else 0.,kind=kind,seed=seed,
        sigma=sigma,origin=origin,native_horizon=horizon,scalar_cadence=1.,profile_cadence=4.,latent_gaussian_sha256=f'latent_{seed}')


def test_sigma_bootstrap_keys_are_paired_by_native_seed_block():
    rows=[row(sigma,seed) for sigma in (.25,.5,1.,2.) for seed in range(12)]
    rows+=[row(0.,0,'mean'),row(0.,0,'zero')]
    selected=S.select_stage(rows,'E1')
    groups=[S.condition_groups(selected,'E1',s)[0][0] for s in (.25,.5,1.,2.)]
    assert len({g['bootstrap_pairing_key'] for g in groups})==1
    assert all([r['seed'] for r in g['noise']]==list(range(12)) for g in groups)


def test_mismatched_sigma_seed_is_not_silently_paired():
    rows=[row(sigma,seed) for sigma in (.25,.5,1.,2.) for seed in range(12)]
    rows[0]['seed']=100
    with pytest.raises(ValueError,match='common latent'):S.select_stage(rows,'E1')


def test_matching_integer_seeds_with_changed_gaussian_knots_are_rejected():
    rows=[row(sigma,seed) for sigma in (.25,.5,1.,2.) for seed in range(12)]
    rows[0]['latent_gaussian_sha256']='different'
    with pytest.raises(ValueError,match='Gaussian-knot hashes'):S.select_stage(rows,'E1')


def test_horizon_matched_controls_are_required():
    rows=[row(seed=s) for s in (1,2)]+[row(0.,0,'mean',horizon=800.),row(0.,0,'zero')]
    with pytest.raises(ValueError,match='horizon-matched'):S.condition_groups(rows,'E2',1.)


def test_e4_does_not_pool_original_short_observations():
    rows=[row(origin='E0'),row(origin='E4',horizon=800.)]
    selected=S.select_stage(rows,'E4')
    assert len(selected)==1 and selected[0]['native_horizon']==800.


def test_split_e1_native_runs_cannot_hide_an_entire_missing_condition():
    base=[row(0,0,'zero')]
    selection=dict(e1_strengths={'m':[.02,.03]})
    only_ready=[row(origin='E1')]
    with pytest.raises(ValueError,match='full declared condition coverage'):S.declared_coverage(only_ready,'E1',base,selection)
    complete=only_ready+[dict(only_ready[0],f=.03)]
    S.declared_coverage(complete,'E1',base,selection)


def test_e4_literal_trigger_set_is_not_capped_by_rough_run_estimate():
    trigger=[dict(model='m',C_label=20,f=.02,original_horizon=200,trigger=True),dict(model='m',C_label=20,f=.03,original_horizon=200,trigger=True)]
    with pytest.raises(ValueError,match='full declared condition coverage'):S.declared_coverage([row(origin='E4',horizon=800.)],'E4',[],triggers=trigger)


def test_duplicate_independent_identity_is_an_error():
    a=row();a.update(scalar_sha256='a',profile_sha256='b')
    b=dict(a,scalar_sha256='c',profile_sha256='d')
    with pytest.raises(ValueError,match='observation-companion'):S.deduplicate([a,b])


def test_finer_observation_companion_counts_once_after_exact_scalar_check(tmp_path):
    p=tmp_path/'scalar.npz';np.savez(p,time=[0.,1.],P_real=[1.,1.],P_imag=[0.,0.],contrast=[1.,1.],circumference=20.,current_half_contrast_regions=[1,1])
    a=row();a.update(scalar_path=str(p),scalar_sha256='a',profile_sha256='b',profile_cadence=40.)
    b=dict(a,id='fine',profile_sha256='fine',profile_cadence=5.,observation_companion_of='retained')
    selected,replaced=S.deduplicate([a,b])
    assert len(selected)==1 and selected[0]['id']=='fine' and len(replaced)==1


def test_additional_forced_state_requires_observed_mean_agreement():
    metrics={s:dict(regime='multiple_fronts') for s in (0.,.25,.5,1.)}
    assert S.abnormal_persistence(metrics)['forced_E1'] is True
    metrics[0.]['regime']='single_front'
    result=S.abnormal_persistence(metrics)
    assert result['forced_E1'] is False and result['low_sigma_abnormal_persistence'] is True
    metrics={s:dict(regime='single_front') for s in (0.,.25,.5,1.)}
    assert S.abnormal_persistence(metrics)['forced_E1'] is False
    del metrics[0.]
    assert S.abnormal_persistence(metrics)['forced_E1'] is None


def test_interim_delivery_gate_is_explicit(tmp_path):
    gate=tmp_path/'release.json';gate.write_text(json.dumps(dict(e0_interim_delivered=False)))
    with pytest.raises(ValueError,match='Deliver'):S.release_check(gate)


def test_no_stronger_input_substitutes_for_missing_declared_weak_baseline():
    m=[dict(model='m',C_label=20.,view='main',f=f,n_seeds=12,forced=False,tau_status='resolved') for f in (.001,.01)]
    result=E.fhalf_row(m,{},'m',20.,'main',False)
    assert result['fhalf_status']=='unavailable_missing_declared_weak_input'
    assert result['fhalf_lower']==.0001 and math.isinf(result['fhalf_upper'])


def test_future_native_census_keeps_blocked_spring_basal_exact(tmp_path):
    run=tmp_path/'native';(run/'results').mkdir(parents=True);(run/'inputs').mkdir();folder=run/'workers/0000';folder.mkdir(parents=True)
    checkpoint=run/'inputs/spring.npz';active=np.array([1.,3.,2.]);q=np.arange(3)/3
    np.savez(checkpoint,state=np.r_[active,[10.,11.,12.],4.],material_positions=q)
    basal=dict(model='spring_ring',size=3.,input='inputs/spring.npz',input_sha256=E.digest(checkpoint))
    job=dict(slot=0,model='spring_ring',size=3.,fraction=.02,kind='noise',seed=1,sigma=1.,protocol=dict(horizon=3200.,scalar_interval=.5,profile_interval=.5))
    E.write_json(run/'manifest.json',dict(stage='E2',basals=[basal],jobs=[job]))
    E.write_json(run/'results/status.json',dict(complete=True,missing_slots=[],failed_slots=[],completed=1,total=1))
    status=folder/'status.json';E.write_json(status,dict(complete=True))
    scalar=folder/'trace.npz';profile=folder/'primary.npz';np.savez(scalar,time=[0.,1.]);np.savez(profile,primary=active[None,:])
    case=dict(model='spring_ring',size=3.,fraction=.02,kind='noise',seed=1,sigma=1.,case_id='E2_0000',
        trace=str(scalar),trace_sha256=E.digest(scalar),profile=str(profile),profile_sha256=E.digest(profile),scalar_interval=.5,profile_interval=.5,horizon=3200.,
        new_gaussian_path=False,retained_prefix=dict(old_seed=1,full_prefix_sha256='native-prefix'),paired_sigma1_source=dict(worker='exact-source'))
    E.write_json(folder/'result.json',dict(job=job))
    case.update(native_result_sha256=E.digest(folder/'result.json'),basal_input_sha256=E.digest(checkpoint))
    E.write_json(run/'results/case_census.json',dict(native_run_manifest_sha256=E.digest(run/'manifest.json'),cases=[case]))
    rows=S.native_rows(run,tmp_path/'analysis')
    with np.load(rows[0]['basal_profile_path']) as z:np.testing.assert_array_equal(z['primary'],active)
    assert rows[0]['actualC0']==4 and rows[0]['basal_contrast']==2
    assert rows[0]['new_gaussian_path'] is False and rows[0]['retained_prefix']==case['retained_prefix']
    assert rows[0]['paired_sigma1_source']==case['paired_sigma1_source'] and rows[0]['native_metadata']['retained_prefix']==case['retained_prefix']


def test_variable_group_worker_export_and_independent_audit_roundtrip(tmp_path):
    root=tmp_path/'analysis';(root/'inputs').mkdir(parents=True);(root/'results').mkdir();(root/'analysis_workers').mkdir()
    cases=[];time=np.arange(801.);q=np.arange(8)/8;basal=2+np.cos(2*np.pi*q);bp=root/'inputs/basal.npz'
    np.savez(bp,primary=basal,material_positions=q,circumference=20.)
    for i,kind in enumerate(('noise','noise','mean','zero')):
        folder=root/f'input_{i}';folder.mkdir();phase=.02*time if kind=='noise' else np.zeros(len(time));P=np.exp(1j*phase)
        scalar=folder/'trace.npz';profile=folder/'primary.npz';status=folder/'status.json';status.write_text('{"complete":true}')
        np.savez(scalar,time=time,P_real=P.real,P_imag=P.imag,M1_real=P.real,M1_imag=P.imag,circumference=20.,contrast=np.full(len(time),2.),current_half_contrast_regions=np.ones(len(time),int))
        pt=time[::4];np.savez(profile,time=pt,primary=2+np.cos(2*np.pi*q[None,:]-phase[::4,None]),x=q*20.,circumference=20.)
        r=row(seed=i+1 if kind=='noise' else 0,kind=kind,origin='E2',horizon=800.);r['id']=str(i)
        for role,path in [('scalar',scalar),('profile',profile),('source_status',status),('basal_profile',bp)]:r[role+'_path']=str(path);r[role+'_sha256']=E.digest(path)
        r.update(basal_contrast=2.,actualC0=20.,native_cells=8)
        cases.append(r)
    groups,_=S.condition_groups(cases,'E2',1.)
    E.write_json(root/'inventory.json',dict(cases=cases))
    E.write_json(root/'inputs/estimands.json',{})
    E.write_json(root/'results/native_identity_verification.json',dict(passed=True))
    sources={str(Path(x.__file__).resolve()):E.digest(x.__file__) for x in (E,E.A,Q)}
    manifest=dict(groups=groups,views=E.VIEWS,stage='E2',sigma=1.,scope='Synthetic test only',census_counts=dict(noise=2,mean=1,zero=1),source_sha256=sources,
        inventory_sha256=E.digest(root/'inventory.json'),estimands_sha256=E.digest(root/'inputs/estimands.json'),native_identity_sha256=E.digest(root/'results/native_identity_verification.json'))
    E.write_json(root/'analysis_manifest.json',manifest);E.write_json(root/'collector_freeze_v3.json',dict(source_sha256=sources))
    E.analyze_group(root,0)
    audited=Q.audit(root,0)
    assert audited['passed'] and audited['noisy_realizations']==2
    exported=E.collect(root)
    assert exported['passed'] and exported['groups']==1
    gate=E.read_json(root/'results/analysis/status.json')
    assert gate['condition_groups']==1 and gate['noisy_realizations']==2
    assert E.read_json(root/'analysis_workers/0000/result.json')['case_metrics'][14]['per_case_tau_status']=='right_censored'


def test_paired_sigma_aggregation_has_one_deterministic_mean_and_no_forced_healthy_front(tmp_path):
    root=tmp_path/'stage';root.mkdir();(root/'results').mkdir();shared=root/'shared';shared.mkdir()
    q=np.arange(8)/8;time=np.arange(801.);basal=2+np.cos(2*np.pi*q)
    bp=shared/'basal.npz';np.savez(bp,primary=basal,material_positions=q,circumference=20.)
    sources={str(Path(x.__file__).resolve()):E.digest(x.__file__) for x in (S,E,E.A,Q,S.I)}
    def make_case(folder,kind,seed,sigma):
        folder.mkdir();phase=.01*sigma*time*(1+seed/12) if kind=='noise' else np.zeros(len(time));P=np.exp(1j*phase)
        scalar=folder/'trace.npz';profile=folder/'primary.npz';status=folder/'status.json';status.write_text('{"complete":true}')
        np.savez(scalar,time=time,P_real=P.real,P_imag=P.imag,M1_real=P.real,M1_imag=P.imag,circumference=20.,contrast=np.full(len(time),2.),current_half_contrast_regions=np.ones(len(time),int))
        np.savez(profile,time=time[::4],primary=2+np.cos(2*np.pi*q[None,:]-phase[::4,None]),x=q*20,circumference=20.)
        r=row(sigma,seed,kind,horizon=800.);r['id']=f'{kind}_{seed}' if kind!='noise' else f'{sigma}_{seed}'
        for role,path in [('scalar',scalar),('profile',profile),('source_status',status),('basal_profile',bp)]:r[role+'_path']=str(path);r[role+'_sha256']=E.digest(path)
        r.update(basal_contrast=2.,actualC0=20.,native_cells=8)
        return r
    controls=[make_case(shared/k,k,0,0.) for k in ('mean','zero')];arms=[]
    for sigma in (.25,.5,1.,2.):
        child=root/f'sigma_{sigma}';child.mkdir();(child/'inputs').mkdir();(child/'results').mkdir();(child/'analysis_workers').mkdir()
        cases=[make_case(child/f'input_{seed}','noise',seed,sigma) for seed in range(12)]+controls
        groups,_=S.condition_groups(cases,'E1',sigma)
        E.write_json(child/'inventory.json',dict(cases=cases));E.write_json(child/'inputs/estimands.json',{});E.write_json(child/'results/native_identity_verification.json',dict(passed=True))
        m=dict(groups=groups,views=E.VIEWS,stage='E1',sigma=sigma,scope='Synthetic only',census_counts=dict(noise=12,mean=1,zero=1),source_sha256=sources,
            inventory_sha256=E.digest(child/'inventory.json'),estimands_sha256=E.digest(child/'inputs/estimands.json'),native_identity_sha256=E.digest(child/'results/native_identity_verification.json'))
        E.write_json(child/'analysis_manifest.json',m);E.write_json(child/'collector_freeze_v3.json',dict(source_sha256=sources))
        E.analyze_group(child,0);Q.audit(child,0);E.collect(child);Q.collect(child)
        arms.append(dict(path=str(child),sigma=sigma,condition_groups=1))
    E.write_json(root/'results/arms.json',dict(stage='E1',arms=arms))
    E.write_json(root/'freeze.json',dict(source_sha256=sources,E1_forced_rule=dict(root_approved=True,rule='same_abnormal_majority_at_sigma1_both_lower_and_mean')))
    assert S.collect(root)['passed']
    import csv
    flags=list(csv.DictReader((root/'results/analysis/E1_forced_flags.csv').open()))
    assert all(r['forced_E1']=='False' for r in flags)
    control=list(csv.DictReader((root/'results/analysis/control_metrics.csv').open()))
    assert len(control)==14 and all(r['bootstrap_resamples']=='0' and r['n_seeds']=='0' for r in control)
    compare=list(csv.DictReader((root/'results/analysis/sigma_comparisons.csv').open()))
    assert len(compare)==21 and all(r['paired_seed_blocks']=='12' for r in compare)
    for r in compare:
        expected=float(r['sigma'])**2
        assert math.isclose(float(r['Dtheta_ratio']),expected,rel_tol=1e-8)
        assert math.isclose(float(r['Dtheta_ratio_ci_lower']),expected,rel_tol=1e-8)
        assert math.isclose(float(r['Dtheta_ratio_ci_upper']),expected,rel_tol=1e-8)


def test_original_seed_sensitivity_is_a_separate_long_arm_not_new_realizations(tmp_path,monkeypatch):
    monkeypatch.setattr(S.I,'PACKAGE',tmp_path)
    parent=tmp_path/'data/runs/parent';(parent/'results/analysis').mkdir(parents=True)
    rows=[dict(row(seed=i,origin='E4',horizon=800),new_gaussian_path=i>=2,retained_prefix={'old_seed':i} if i<2 else None) for i in range(12)]
    rows+=[row(0,0,'mean',origin='E4',horizon=800),row(0,0,'zero',origin='E4',horizon=800)]
    E.write_json(parent/'inventory.json',dict(cases=rows))
    E.write_json(parent/'results/analysis/status.json',dict(passed=True,complete=True,stage='E4',outputs={}))
    E.write_json(parent/'freeze.json',dict(stage='E4',inventory_sha256=E.digest(parent/'inventory.json'),release_path='release',release_sha256='sha',estimands_path='estimands',estimands_sha256='sha'))
    out=tmp_path/'data/runs/originals'
    result=S.prepare_original_seed_sensitivity(out,parent)
    assert result['stage']=='E4_original_seed_sensitivity' and result['counts']==dict(noise=2,mean=1,zero=1)
    inventory=E.read_json(out/'inventory.json')
    assert {r['seed'] for r in inventory['cases'] if r['kind']=='noise'}=={0,1}
    assert all(r['native_horizon']==800 for r in inventory['cases'])
    assert E.read_json(out/'freeze.json')['stage']=='E4_original_seed_sensitivity'
    rows[0].pop('new_gaussian_path');E.write_json(parent/'inventory.json',dict(cases=rows))
    f=E.read_json(parent/'freeze.json');f['inventory_sha256']=E.digest(parent/'inventory.json');E.write_json(parent/'freeze.json',f)
    with pytest.raises(ValueError,match='explicitly distinguish'):S.prepare_original_seed_sensitivity(tmp_path/'data/runs/missing',parent)


def test_native_engineering_gate_covers_exact_frozen_workers_and_unchanged_results(tmp_path):
    run=tmp_path/'native';(run/'results').mkdir(parents=True);folder=run/'workers/0000';folder.mkdir(parents=True)
    E.write_json(run/'manifest.json',dict(jobs=[dict(slot=0)]))
    E.write_json(folder/'result.json',dict(native='declared'));E.write_json(folder/'status.json',dict(complete=True))
    audit=run/'results/audit_0000.json';E.write_json(audit,dict(passed=True,slot=0,native_manifest_sha256=E.digest(run/'manifest.json'),native_result_sha256=E.digest(folder/'result.json'),native_status_sha256=E.digest(folder/'status.json'),hashed_native_file_snapshots={}))
    gate=dict(passed=True,full_native_run_covered=True,missing_tasks=[],failed_tasks=[],native_manifest_sha256=E.digest(run/'manifest.json'),receipts=[dict(slot=0,path=str(audit),sha256=E.digest(audit),passed=True)])
    E.write_json(run/'results/complete_native_independent_audit_20261020.json',gate)
    assert S.require_complete_native_gate(run)['workers']==1
    gate['receipts']=[];E.write_json(run/'results/complete_native_independent_audit_20261020.json',gate)
    with pytest.raises(ValueError,match='every frozen worker'):S.require_complete_native_gate(run)
    gate['receipts']=[dict(slot=0,path=str(audit),sha256=E.digest(audit),passed=True)];E.write_json(run/'results/complete_native_independent_audit_20261020.json',gate)
    E.write_json(folder/'result.json',dict(native='changed-after-audit'))
    with pytest.raises(ValueError,match='Immutable input hash changed'):S.require_complete_native_gate(run)
