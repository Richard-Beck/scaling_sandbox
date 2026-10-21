from pathlib import Path
import sys
import math
import numpy as np
import pytest

CODE=Path(__file__).resolve().parents[1]/'mechanisms'
sys.path.insert(0,str(CODE))
import signature_figure4_final_statistics_20261020 as F


def metric(*,f=.03,tau=20,lower=None,upper=None,status='resolved',coverage=1.,regime='single_front',forced=False,C=60,view='main',sigma=1.):
    lo=tau if lower is None else lower;hi=tau if upper is None else upper
    return dict(model='wave_pinning',C_label=C,f=f,view=view,sigma=sigma,tau_status=status,tau=tau if status=='resolved' else None,
        tau_lower=lo,tau_upper=hi,tau_ci_lower=lo,tau_ci_upper=hi,tau_coverage=coverage,regime=regime,forced=forced,Dtheta_ci_lower=.01,Dtheta_ci_upper=.02,Dtheta=.015,actual_C_mean=C)


def ratio(f,point=2.,ci=(1.5,2.5),multiple=.8,forced=False,evaluable=True):
    return dict(f=f,tau_ratio=point,tau_ratio_lower=point,tau_ratio_upper=point,tau_ratio_ci_lower=ci[0],tau_ratio_ci_upper=ci[1],
        large_multiple_front_fraction=multiple,forced=forced,evaluable=evaluable,coverage_small=1.,coverage_large=1.)


def feature(name,rho=-.8,success=3,unknown=0,outcome='direct_tau'):
    return dict(statistic=name,pooled_spearman=rho,successful_folds=success,unavailable_folds=unknown,planned_folds=4,
        primary_view=True,primary_H2_feature=True,imaging_view='primary',outcome=outcome,forced_inclusive=False)


def test_unbounded_json_interval_is_not_missing_or_zero():
    assert F.interval(dict(tau_lower=8000,tau_upper=None,tau_upper_unbounded=True),'tau')==(8000,math.inf)
    assert F.interval(dict(tau_lower=8000,tau_upper=None),'tau') is None


def test_no_evaluable_conditions_cannot_vacuously_support_h1():
    r=F.h1_model([],None,None)
    assert r['status']=='inconclusive'


def test_h1_needs_noise_tolerance_and_all_ratio_cis():
    rows=[ratio(.01),ratio(.03)]
    small=dict(fhalf_ci_lower=.01,fhalf_ci_upper=.02);large=dict(fhalf_ci_lower=.04,fhalf_ci_upper=.05)
    assert F.h1_model(rows,small,large)['status']=='pass'
    assert F.h1_model(rows,small,dict(fhalf_ci_lower=.015,fhalf_ci_upper=.05))['status']=='inconclusive'
    rows.append(ratio(.1,point=.5,ci=(.3,.8)))
    assert F.h1_model(rows,small,large)['status']=='fail'


def test_forced_inclusive_is_a_separate_h1_result():
    rows=[ratio(.01),ratio(.03,point=.5,ci=(.3,.8),forced=True)]
    s=dict(fhalf_ci_lower=.01,fhalf_ci_upper=.02);l=dict(fhalf_ci_lower=.04,fhalf_ci_upper=.05)
    assert F.h1_model(rows,s,l)['status']=='pass'
    assert F.h1_model(rows,s,l,inclusive=True)['status']=='fail'


def test_unknown_required_forced_flag_is_not_nonforced():
    row=dict(forced=False,E1_required=True,forced_E1=None)
    assert F.forced(row) is None and not F.include(row,False)


def test_h2_peak_drop_requires_all_designed_range_points():
    rows=[ratio(.01,3.,(2,4)),ratio(.03,2.,(1.5,2.5)),ratio(.1,1.,(.8,1.2))]
    assert F.h2_model(rows)['status']=='pass'
    rows[1]['forced']=True
    assert F.h2_model(rows)['status']=='inconclusive'
    assert F.h2_model(rows,inclusive=True)['status']=='pass'


def test_h2_can_use_certified_censored_peak_without_inventing_a_point():
    peak=ratio(.01,None,(3.,math.inf));peak.update(tau_ratio_lower=3.,tau_ratio_upper=math.inf)
    end=ratio(.1,1.,(.8,1.2))
    assert F.h2_model([peak,end])['status']=='pass'


def test_a6_gate_needs_declared_direct_feature_and_three_of_four_folds():
    rows=[feature(n) for n in F.PRIMARY_A6]
    assert F.a6_gate(rows)['status']=='pass'
    assert F.a6_gate([feature(n,success=2,unknown=2) for n in F.PRIMARY_A6])['status']=='inconclusive'
    assert F.a6_gate([feature(n,success=2,unknown=0) for n in F.PRIMARY_A6])['status']=='fail'
    assert F.a6_gate([feature(n,outcome='Dtheta_proxy') for n in F.PRIMARY_A6])['status']=='inconclusive'


def test_missing_a6_fold_cannot_be_a_success():
    rows=[feature(n,success=2,unknown=1) for n in F.PRIMARY_A6]
    assert F.a6_gate(rows)['status']=='inconclusive'


def test_strong_input_uses_censored_weak_lower_bound():
    weak=metric(f=.0001,tau=None,lower=8000,upper=math.inf,status='right_censored')
    assert F.strong_input(metric(tau=2000),weak)['status']=='pass'
    assert F.strong_input(metric(tau=5000),weak)['status']=='inconclusive'


def test_contrast_loss_does_not_create_a_tau_upper_bound():
    weak=metric(f=.0001,tau=None,lower=8000,upper=math.inf,status='right_censored')
    low=metric(tau=None,lower=0,upper=math.inf,status='unavailable_no_contiguous_lag_curve',coverage=0.,regime='low_contrast')
    assert F.strong_input(low,weak)['status']=='inconclusive'


def test_strong_input_must_keep_actual_weak_baseline_and_small_size_coverage():
    weak=metric(f=.0001,tau=8000)
    assert F.strong_input(metric(tau=1000),dict(weak,f=.001))['status']=='inconclusive'
    assert F.strong_input(metric(tau=1000),weak,small_strong=metric(coverage=.1))['status']=='inconclusive'


def test_e1_existence_rule_reports_all_preset_sigmas():
    rows=[metric(regime='multiple_fronts',sigma=s) for s in (.25,.5,1.,2.)]
    rows[0]['regime']='single_front';rows[1]['regime']='single_front'
    result=F.e1_mode_corroboration(rows,'multiple_fronts')
    assert result['status']=='pass' and len(result['sigma_checks'])==3
    assert [r['sigma'] for r in result['sigma_checks'] if r['same_abnormal_majority']]==[2.]


def test_e1_missing_sigma_is_missing_work_not_selected_success():
    rows=[metric(regime='multiple_fronts',sigma=s) for s in (1.,2.)]
    assert F.e1_mode_corroboration(rows,'multiple_fronts')['status']=='inconclusive'
    rows.append(dict(rows[-1],sigma=.5,f=.1))
    with pytest.raises(ValueError,match='same model'):F.e1_mode_corroboration(rows,'multiple_fronts')


def test_fingerprint_frames_are_time_weighted_within_one_seed():
    answer=F.time_weighted_seed_summary([0,1,4],[1,3,1000],[2,4,1000],2)
    np.testing.assert_allclose(answer,[2.5,1.75])


def test_whole_seed_size_effects_have_joint_weights_and_no_frame_replication():
    result=F.fingerprint_effects([[1,1],[1,1]],[[3,.25],[3,.25]],draws=2000)
    assert result['region_count']['estimate']==3 and result['region_count']['ci_lower']==3
    assert result['normalized_contrast']['ci_upper']==.25
    assert result['region_count']['independent_small_seeds']==2


def test_zero_count_denominator_is_unavailable_not_pseudocounted():
    result=F.fingerprint_effects([[0,1],[0,1]],[[3,.25],[3,.25]])
    assert result['region_count']['available'] is False
    assert result['region_count']['estimate'] is None and math.isinf(result['region_count']['ci_upper'])


def test_distinct_modes_must_be_in_different_models_and_both_fingerprints_separate():
    pass_=F.decision('test','pass','analytic')
    a=dict(model='a',f=.03,mode='multiple_fronts',strong=pass_,E1=pass_,fingerprints=F.fingerprint_effects([[1,1],[1,1]],[[3,.8],[3,.8]]))
    b=dict(model='b',f=.03,mode='low_contrast',strong=pass_,E1=pass_,fingerprints=F.fingerprint_effects([[1,1],[1,1]],[[.2,.1],[.2,.1]]))
    assert F.h3_overall([a,b])['status']=='pass'
    assert F.h3_overall([a,dict(b,model='a')])['status']=='inconclusive'
    b['fingerprints']['normalized_contrast']=a['fingerprints']['normalized_contrast']
    assert F.h3_overall([a,b])['status']=='inconclusive'


def test_disjoint_censored_window_intervals_are_sensitive():
    a=metric(tau=None,lower=8000,upper=math.inf,status='right_censored');b=metric(tau=1000)
    assert F.interval_and_label_compatibility(a,b,label='window')['status']=='fail'


def test_majority_change_remains_visible_with_missing_heading_coverage():
    a=metric();b=metric(tau=None,lower=0,upper=math.inf,status='unavailable',coverage=0,regime='low_contrast')
    result=F.interval_and_label_compatibility(a,b,label='window')
    assert result['status']=='fail' and result['regime_changed'] and result['persistence_comparison_unavailable']


def test_missing_window_is_not_a_robustness_pass():
    views={v:metric(view=v) for v in F.VIEWS};del views['late_quarter']
    assert F.window_qualification(views)['status']=='inconclusive'


def test_e4_can_uncensor_a_criterion_without_common_support():
    yes=F.decision('x','pass','test');unknown=F.decision('x','inconclusive','censored')
    assert F.support_label(unknown,yes,yes,common_only_censored=True,horizon=yes)=='qualified_long_horizon_only'
    assert F.support_label(yes,yes,yes,horizon=yes)=='robust_common_and_long'
    no=F.decision('x','fail','changed majority')
    assert F.support_label(unknown,yes,no,common_only_censored=True)=='sensitive_support'


def test_preselected_e4_has_no_silent_common_fallback():
    base=[metric(C=20),metric(C=60)];triggers=[dict(model='wave_pinning',C_label=60,f=.03,trigger=True)]
    with pytest.raises(ValueError,match='fallback'):F.assemble_metrics(base,[],[],triggers)
    long=[dict(metric(C=60),long_arm_qualified=False)]
    with pytest.raises(ValueError,match='qualification'):F.assemble_metrics(base,[],long,triggers)


def test_common_assembly_keeps_dense_inputs_and_only_new_e3_sizes():
    base=[metric(C=20),metric(C=60),metric(C=60,f=.005)]
    e3=[metric(C=20),metric(C=30),metric(C=50),metric(C=60)]
    long=[dict(metric(C=60),tau=100,long_arm_qualified=True)]
    triggers=[dict(model='wave_pinning',C_label=60,f=.03,trigger=True)]
    common,longest,checks=F.assemble_metrics(base,e3,long,triggers)
    assert len(common)==5 and len(checks)==2
    assert next(r for r in longest if r['C_label']==60 and r['f']==.03)['tau']==100
    assert next(r for r in longest if r['f']==.005)['selected_observation_arm']=='sole_common_arm_not_selected_for_E4'


def test_overall_counts_keep_all_four_models():
    p=F.decision('x','pass','test');u=F.decision('x','inconclusive','missing')
    result=F.at_least('overall',[p,p,u,u],3)
    assert result['status']=='inconclusive' and result['total']==4


def test_inclusive_a6_keeps_its_separate_nonprimary_view_schema():
    rows=[dict(feature(name),forced_inclusive=True,primary_view=False) for name in F.PRIMARY_A6]
    assert F.a6_gate(rows,inclusive=True)['status']=='pass'
    assert F.a6_gate(rows)['status']=='inconclusive'


def test_a_missing_majority_label_is_not_window_compatibility():
    assert F.interval_and_label_compatibility(metric(),dict(metric(),regime=None),label='window')['status']=='inconclusive'
    assert F.interval_and_label_compatibility(metric(),dict(metric(),regime='low_contrast',tau_ci_lower=None),label='window')['status']=='fail'


def test_physical_fingerprint_window_uses_overlap_not_frame_replication():
    frames=[dict(time=t,region_count=n,contrast=c) for t,n,c in [(0,1,2),(2,3,4),(8,1000,1000)]]
    np.testing.assert_allclose(F.physical_seed_summary(frames,2,1,5),[2.5,1.75])
    with pytest.raises(ValueError,match='complete observed coverage'):F.physical_seed_summary(frames,2,-1,5)


def test_final_forced_join_keeps_missing_required_qualification_explicit():
    m=metric();flags=[dict(model=m['model'],C_label=60,f=.03,view='main',forced_E1=None)]
    row=F.attach_forced_flags([m],flags,{F.condition_key(m)})[0]
    assert row['E1_required'] and row['forced_combined'] is None
    with pytest.raises(ValueError,match='missing work'):F.attach_forced_flags([m],[],{F.condition_key(m)})
    row=F.attach_forced_flags([m],[],set())[0]
    assert not row['E1_required'] and row['forced_combined'] is False


def test_long_qualification_binds_every_selected_condition():
    row=metric();triggers=[dict(row,trigger=True)]
    gate=dict(passed=True,complete=True,all_preselected_arms_complete=True,all_required_prefix_checks_passed=True,all_required_controls_passed=True,qualified_conditions=[row])
    assert F.qualify_long_rows([row],gate,triggers)[0]['long_arm_qualified']
    gate['qualified_conditions']=[]
    with pytest.raises(ValueError,match='preselected'):F.qualify_long_rows([row],gate,triggers)


def test_tolerance_uses_whole_seed_intervals_and_does_not_invent_weak_baseline():
    rows=[metric(f=.0001,tau=100),metric(f=.1,tau=10)]
    arrays={F.condition_key(r):{'main__tau_boot_lower':np.full(2000,r['tau']),'main__tau_boot_upper':np.full(2000,r['tau'])} for r in rows}
    result=F.noise_tolerance(rows,arrays,'wave_pinning',60.,'main')
    assert result['fhalf_status']=='resolved' and result['fhalf_ci_lower']==result['fhalf_ci_upper']
    assert F.noise_tolerance(rows[1:],arrays,'wave_pinning',60.,'main')['fhalf_status']=='unavailable_missing_declared_weak_input'


def test_unqualified_heading_cannot_supply_a_finite_tolerance_bound():
    rows=[metric(f=.0001,tau=10000),metric(f=.1,tau=10,coverage=.1)]
    arrays={F.condition_key(r):{'main__tau_boot_lower':np.full(2000,r['tau']),'main__tau_boot_upper':np.full(2000,r['tau'])} for r in rows}
    result=F.noise_tolerance(rows,arrays,'wave_pinning',60.,'main')
    assert result['fhalf_status']!='resolved' and math.isinf(result['fhalf_ci_upper'])


def analytic_full_arm():
    metrics=[];arrays={};e1=[];summaries={}
    for model in F.MODELS:
        for C in (20.,60.):
            for f in (.0001,.03,.1):
                values={}
                for view in F.VIEWS:
                    tau=(10000 if C==20 else 20000) if f==.0001 else (4000 if C==20 else 1000)
                    mode='single_front' if f==.0001 else 'multiple_fronts' if model in ('wave_pinning','holmes_model3') else 'low_contrast'
                    r=dict(metric(f=f,tau=tau,regime=mode,C=C,view=view),model=model,n_seeds=2,multiple_fronts_fraction=.8 if mode=='multiple_fronts' else 0.)
                    metrics.append(r)
                    for suffix in ('tau_boot_lower','tau_boot_upper'):values[view+'__'+suffix]=np.full(2000,tau)
                    for suffix in ('Dtheta_boot','Cmean_boot'):values[view+'__'+suffix]=np.ones(2000)
                    summaries[F.metric_key(r)]=np.tile([1.,1.] if C==20 else [3.,.8] if mode=='multiple_fronts' else [.2,.1],(2,1))
                    if f!=.0001:
                        for sigma in (.25,.5,1.,2.):e1.append(dict(r,sigma=sigma))
                arrays[(model,C,f)]=values
    features=[dict(feature(name),forced_inclusive=inclusive,primary_view=not inclusive) for inclusive in (False,True) for name in F.PRIMARY_A6]
    return metrics,arrays,e1,summaries,features


def test_full_arm_keeps_all_models_windows_forced_views_and_candidate_inputs():
    metrics,arrays,e1,summaries,features=analytic_full_arm()
    result=F.evaluate_arm(metrics,arrays,e1,summaries,features,arm='analytic')
    assert len(result['criteria'])==154
    assert len(result['strong_input_candidates'])==112
    assert len(result['window_qualification'])==24
    h3=[r for r in result['criteria'] if r['criterion']=='H3_overall']
    assert len(h3)==14 and all(r['status']=='pass' for r in h3)
    assert len(result['fingerprint_effects'])==168
    assert all(r['independent_small_seeds']==2 for r in result['fingerprint_effects'])


def test_missing_declared_window_is_execution_error_not_censoring():
    metrics,arrays,e1,summaries,features=analytic_full_arm()
    with pytest.raises(ValueError,match='missing work'):F.evaluate_arm(metrics[1:],arrays,e1,summaries,features,arm='analytic')


def test_longest_unselected_intermediate_size_retains_e3_source():
    base=[metric(C=20),metric(C=60)];e3rows=[metric(C=30),metric(C=50)]
    common,longest,_=F.assemble_metrics(base,e3rows,[],[])
    def dataset(rows):return dict(arrays={F.condition_key(r):r for r in rows},groups={F.condition_key(r):r for r in rows})
    arrays,groups=F.select_observation_inputs(longest,dataset(base),dataset(e3rows),dataset([]),longest=True)
    assert {k[1] for k in arrays}=={20,30,50,60}


def test_comparison_export_keeps_an_unmatched_designed_input_explicit():
    metrics,arrays,_,_,_=analytic_full_arm()
    extra=dict(metrics[0],f=.015)
    metrics.append(extra);arrays[F.condition_key(extra)]=arrays[F.condition_key(metrics[0])]
    ratios,_=F.comparison_tables(metrics,arrays,arm='analytic')
    row=next(r for r in ratios if r['f']==.015)
    assert not row['evaluable'] and row['comparison_status']=='unavailable_unmatched_designed_input'


def test_ratio_proxy_and_geometry_use_actual_circumference_and_seed_draws():
    sm=dict(metric(C=20),n_seeds=2,multiple_fronts_fraction=.8,actual_C_mean=18,Dtheta=2)
    lg=dict(metric(C=60),n_seeds=2,multiple_fronts_fraction=.8,actual_C_mean=54,Dtheta=.5)
    def arrays(row):return {'main__tau_boot_lower':np.full(2000,20),'main__tau_boot_upper':np.full(2000,20),'main__Dtheta_boot':np.full(2000,row['Dtheta']),'main__Cmean_boot':np.full(2000,row['actual_C_mean'])}
    result=F.bootstrap_ratio_row(sm,lg,arrays(sm),arrays(lg))
    assert result['actual_geometry_ratio']==9 and result['Dtheta_mobility_proxy']==4
    assert result['Dtheta_mobility_proxy_ci_lower']==4 and result['Dx_large_over_small']==2.25


def test_csv_contract_keeps_nullable_flags_and_nested_unbounded_intervals():
    import csv,io,json
    rows=[dict(forced_combined=None,sigma=1.,bounds=dict(ci_upper=math.inf,ci_upper_unbounded=True),label='native raw; optical fingerprint')]
    decoded=list(csv.DictReader(io.StringIO(F.csv_bytes(rows).decode())))[0]
    assert decoded['forced_combined']=='' and decoded['sigma']=='1.0'
    assert json.loads(decoded['bounds'])==dict(ci_upper=None,ci_upper_unbounded=True)


def audit_spec(tmp_path,value,name='audit.json'):
    import json
    path=tmp_path/name;path.write_text(json.dumps(value))
    return dict(path=str(path),sha256=F.sha(path))


def test_a5_independent_gate_binds_inventory_all_controls_and_measured_hashes(tmp_path):
    accepted={key:dict(sha256=key+'-measurement') for key in ('noise','mean','zero')}
    records=[dict(case_id=key,passed=True,measurement_sha256=row['sha256']) for key,row in accepted.items()]
    master=dict(passed=True,inventory_sha256='selected-inventory',cases=3,records=records)
    assert len(F.validate_A5_audits([audit_spec(tmp_path,master)],'selected-inventory',list(accepted),accepted)[0])==3
    with pytest.raises(ValueError,match='different canonical inventory'):
        F.validate_A5_audits([audit_spec(tmp_path,dict(master,inventory_sha256='unrelated-E0'))],'selected-inventory',list(accepted),accepted)
    with pytest.raises(ValueError,match='case count disagrees|every canonical noisy case and control'):
        F.validate_A5_audits([audit_spec(tmp_path,dict(master,cases=3,records=records[:-1]))],'selected-inventory',list(accepted),accepted)
    changed=[dict(r) for r in records];changed[0]['measurement_sha256']='changed-result'
    with pytest.raises(ValueError,match='measured result changed'):
        F.validate_A5_audits([audit_spec(tmp_path,dict(master,records=changed))],'selected-inventory',list(accepted),accepted)


def test_legacy_a5_master_requires_inventory_bound_hashed_shards_with_exact_coverage(tmp_path):
    accepted={key:dict(sha256=key) for key in ('noise','control')}
    shard=dict(passed=True,inventory_sha256='selected',cases=2,records=[dict(case_id=key,passed=True,measurement_sha256=key) for key in accepted])
    master=dict(passed=True,native_cases=2,shards=[audit_spec(tmp_path,shard,'shard.json')])
    assert len(F.validate_A5_audits([audit_spec(tmp_path,master)],'selected',list(accepted),accepted)[0])==2
    disagree=dict(master,records=[dict(case_id=key,passed=True,measurement_sha256='wrong') for key in accepted])
    with pytest.raises(ValueError,match='master records disagree'):
        F.validate_A5_audits([audit_spec(tmp_path,disagree)],'selected',list(accepted),accepted)
    shard['inventory_sha256']='old-E0';master['shards']=[audit_spec(tmp_path,shard,'shard.json')]
    with pytest.raises(ValueError,match='not bound to the selected inventory'):
        F.validate_A5_audits([audit_spec(tmp_path,master)],'selected',list(accepted),accepted)
    with pytest.raises(ValueError,match='lacks exact canonical case records'):
        F.validate_A5_audits([audit_spec(tmp_path,dict(passed=True,inventory_sha256='selected',cases=2))],'selected',list(accepted),accepted)


def finite_resolution_bound():
    return dict(metric(tau=None,lower=0,upper=10,status='unresolved_below_two_cadences'),cadence=5.,tau_crossing_bracket='[0,5]')


def test_genuine_sub_two_cadence_crossing_is_a_bound_without_a_point_tau():
    strong=finite_resolution_bound();weak=metric(f=.0001,tau=None,lower=8000,upper=math.inf,status='right_censored')
    assert F.persistence_eligible(strong) and strong['tau'] is None
    result=F.strong_input(strong,weak)
    assert result['status']=='pass' and result['strong_bounds']==(0.,10.) and result['strong_point_tau'] is None


@pytest.mark.parametrize('change',[dict(tau_coverage=.49),dict(tau_status='unavailable_no_contiguous_lag_curve'),dict(tau_crossing_bracket=None),dict(tau_crossing_bracket='[0,100]'),dict(tau_upper=math.inf),dict(tau=3.)])
def test_finite_resolution_bound_never_rescues_missing_crossing_coverage_or_invents_a_point(change):
    assert not F.persistence_eligible(dict(finite_resolution_bound(),**change))


def test_a_resolution_bound_as_ratio_denominator_preserves_unbounded_upper():
    sm=dict(finite_resolution_bound(),C_label=20.,actual_C_mean=20.,n_seeds=12,multiple_fronts_fraction=.8)
    lg=dict(metric(tau=100),n_seeds=12,multiple_fronts_fraction=.8)
    def arrays(lower,upper,C):return dict(main__tau_boot_lower=np.full(2000,lower),main__tau_boot_upper=np.full(2000,upper),main__Dtheta_boot=np.ones(2000),main__Cmean_boot=np.full(2000,C))
    row=F.bootstrap_ratio_row(sm,lg,arrays(0.,10.,20.),arrays(100.,100.,60.))
    assert row['evaluable'] and row['tau_ratio'] is None and row['tau_ratio_lower']==10 and math.isinf(row['tau_ratio_upper'])
