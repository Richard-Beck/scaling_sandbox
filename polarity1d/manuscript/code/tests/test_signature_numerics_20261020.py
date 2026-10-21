"""Conservation, controller-clock and event-preservation numerical regressions."""
from dataclasses import dataclass
import math
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_numerics_20261020 as numerical


def test_fv_overlap_preserves_sharp_profile_amount_and_positivity():
    old=np.r_[np.zeros(121),np.ones(80)]
    new=numerical.fv_remap(old,401)
    np.testing.assert_allclose(new.mean(),old.mean(),rtol=0,atol=2e-15)
    assert new.min()>=0 and new.max()<=1+1e-14
    assert new[0]==0 and new[-1]==pytest.approx(1.)


def test_each_holmes_species_and_both_pools_survive_fv_transfer():
    old=np.random.default_rng(42).uniform(.01,3.,(4,201))
    original=old.copy();new,receipt=numerical.remap_state('holmes_model3',old)
    np.testing.assert_array_equal(old,original)
    np.testing.assert_allclose(new.mean(axis=1),old.mean(axis=1),rtol=2e-14,atol=1e-14)
    assert new.min()>=0 and new.shape==(4,401)
    assert set(receipt['pools'])=={'Rac','Rho'}
    assert receipt['exact_main_checkpoint_replaced'] is False


def test_spring_transfer_preserves_physical_amounts_and_raw_endpoints():
    q=np.linspace(0.,1.,201);u=.4+.2*np.cos(np.pi*q);v=2.-u
    old=np.empty(404);old[:-2:2]=u;old[1:-2:2]=v;old[-2:]=(-7.,-4.5)
    new,receipt=numerical.remap_state('spring',old)
    np.testing.assert_array_equal(new[-2:],old[-2:])
    for component in (0,1):
        np.testing.assert_allclose(np.trapezoid(new[component:-2:2],np.linspace(0.,1.,401))*2.5,
            np.trapezoid(old[component:-2:2],q)*2.5,rtol=2e-14,atol=1e-14)
    assert np.min(new[:-2])>0 and len(receipt['pools'])==2


def test_debelly_preserves_all_eleven_components_and_exact_controller_values():
    rng=np.random.default_rng(7)
    old={name:rng.uniform(.01,1.,201) for name in numerical.DB_NAMES}
    old['displacement']=np.linspace(-2.,3.,201)
    old['velocity']=np.linspace(-.4,.1,201)
    for i,name in enumerate(numerical.DB_CONTROLLERS):old[name]=np.full(201,-.3+i*.01)
    new,receipt=numerical.remap_state('debelly',old)
    assert set(new)==set(old)==set(numerical.DB_NAMES)
    assert all(value.shape==(401,) for value in new.values())
    for name in numerical.DB_CONTROLLERS:assert np.all(new[name]==old[name][0])
    np.testing.assert_allclose((new['mca_bound']+new['mca_unbound']).mean(),
        (old['mca_bound']+old['mca_unbound']).mean(),rtol=1e-14,atol=1e-14)
    assert receipt['companion_only'] is True


def test_negative_densities_and_missing_debelly_state_are_not_repaired_silently():
    with pytest.raises(ValueError,match='Negative'):
        numerical.remap_state('wave_pinning',np.array([[-1.,1.],[2.,2.]]))
    with pytest.raises(ValueError,match='eleven'):
        numerical.remap_state('debelly',{'rac':np.ones(201)})


def test_quarter_step_preserves_physical_filter_clock_while_literal_factor_does_not():
    factor=numerical.rate_preserving_relaxation(.98,.25,.5)
    assert numerical.rate_preserving_relaxation(.98,.5,.5)==.98
    assert factor*factor==pytest.approx(.98,rel=1e-15)
    assert math.log(factor)/.25==pytest.approx(math.log(.98)/.5,rel=1e-14)
    # Retaining literal 0.98 at every quarter step doubles the decay rate.
    assert math.log(.98)/.25==pytest.approx(2*math.log(.98)/.5)


def test_controller_override_restores_the_native_factory_after_failure():
    @dataclass(frozen=True)
    class Parameters:
        boundary_relaxation:float=.98
        rac_deactivation:float=.04
    baseline=SimpleNamespace(parameters=Parameters())
    original=lambda key:baseline
    fake=SimpleNamespace(model_for=original,n=SimpleNamespace(ds=SimpleNamespace(
        DeBelly=lambda p:SimpleNamespace(parameters=p))))
    with pytest.raises(RuntimeError,match='test exit'):
        with numerical.quarter_step_controller(fake):
            configured=fake.model_for('debelly').parameters
            assert configured.boundary_relaxation==numerical.rate_preserving_relaxation()
            assert configured.rac_deactivation==.04
            raise RuntimeError('test exit')
    assert fake.model_for is original and baseline.parameters.boundary_relaxation==.98


def test_observation_coarsening_retains_non_grid_phase_boundaries_and_tail_end():
    times=np.unique(np.r_[np.arange(0.,10.3,.5),10.3,10.3+np.arange(.5,10.,.5),20.3])
    events=[dict(start=0.,stop=10.3),dict(start=10.3,stop=20.3)]
    indices=numerical.observation_indices(times,events,2.5)
    assert {0.,10.3,20.3}.issubset(set(times[indices]))
    assert indices[-1]==len(times)-1
    with pytest.raises(ValueError,match='absent'):
        numerical.observation_indices(times,events,.3)


def test_overlapping_crossing_brackets_are_compared_before_midpoint_precision():
    first=dict(observed=True,lower_time=10.,upper_time=12.5,left_censored=False)
    second=dict(observed=True,lower_time=12.,upper_time=15.,left_censored=False)
    comparison=numerical.event_comparison(first,second,tolerance=0.)
    assert comparison['agrees'] and comparison['brackets_overlap']
    assert not numerical.event_comparison(first,dict(observed=False))['agrees']


def test_fate_or_field_failures_cannot_be_hidden_by_small_final_D_difference():
    baseline=dict(readout_eligible=True,final_d=.49,field_tail_pass=True)
    assert numerical.comparison_classification(baseline,
        dict(baseline,final_d=.51),{'agrees':True})=='replacement_fate_changed'
    assert numerical.comparison_classification(baseline,
        dict(baseline,field_tail_pass=False),{'agrees':True})=='finite_field_tail_status_changed'
    assert numerical.comparison_classification(baseline,
        dict(baseline,maximum_endpoint_scaled_field_difference=.06),{'agrees':True})=='endpoint_native_field_difference_exceeds_declared_tolerance'


def test_retained_collection_controls_are_required_and_read_from_a_bound_sidecar(tmp_path):
    rows=[];controls=[]
    for model in numerical.MODELS:
        for length in numerical.LENGTHS:
            checkpoint=tmp_path/f'{model}_{length}.npz';checkpoint.write_bytes(f'{model}:{length}'.encode())
            identifier=f'{model}:L{length}'
            rows.append(dict(basal_id=identifier,model=model,requested_reference_length=length,
                qualified=True,checkpoint=str(checkpoint),checkpoint_sha256=numerical.sha(checkpoint)))
            controls.append(dict(basal_id=identifier,passed=True,no_cue_control_pass=True,
                                 spontaneous_D_replacement=False))
    path=tmp_path/'basals.json';path.write_text(json.dumps(dict(basals=rows)))
    with pytest.raises(ValueError,match='Missing or ambiguous'):
        numerical._source_basals([path])
    sidecar=tmp_path/'control_flags.json';sidecar.write_text(json.dumps(dict(controls=controls)))
    read=numerical._source_basals([path])
    assert len(read)==10
    assert all(row['source_control_sha256']==numerical.sha(sidecar) for row in read.values())
    controls[-1]['passed']=False;sidecar.write_text(json.dumps(dict(controls=controls)))
    with pytest.raises(ValueError,match='lacks a passed'):
        numerical._source_basals([path])


def test_temporal_companion_pool_gate_keeps_redistribution_distinct_from_loss():
    fields=dict(rac=np.ones((2,4)),rho=np.ones((2,4)),
        mca_bound=np.array([[.2,.2,.2,.2],[.1,.2,.2,.3]]),
        mca_unbound=np.full((2,4),.8))
    trace=SimpleNamespace(fields=fields)
    gate=numerical._debelly_pool_gate(trace,1.)
    assert gate['passed'] and gate['maximum_relative_MCA_pool_error']<1e-14
    fields['mca_unbound'][1]+=1e-5
    assert not numerical._debelly_pool_gate(trace,1.)['passed']
    fields['mca_unbound'][1]-=1e-5;fields['rac'][1,0]=-1e-5
    assert not numerical._debelly_pool_gate(trace,1.)['passed']


def test_missing_native_work_cannot_be_relabelled_as_an_unavailable_calibration(tmp_path):
    run=tmp_path/'calibration';source=run/'workers/0003/result.json';source.parent.mkdir(parents=True)
    result=dict(basal=dict(checkpoint_sha256='exact_basal'),manifest_sha256='manifest',
        calibration=dict(train_ready=False,status='no_sampled_commitment_through_finite_horizon'))
    source.write_text(json.dumps(result));digest=numerical.sha(source)
    (run/'results').mkdir();receipt=run/'results/artifact_verification.json'
    receipt.write_text(json.dumps(dict(workers=[dict(worker_result_sha256=digest,passed=True,
        all_scientific_artifacts_hashed=True)])))
    row=dict(model='holmes_model3',requested_reference_length=39,pulse_count=4,
        gap_multiple_of_pulse=4.,total_on_factor_of_calibration=1.5,
        calibration_result_path=str(source),calibration_result_sha256=digest,
        exact_basal_sha256='exact_basal',calibration_manifest_sha256='manifest',
        calibration_artifact_hashes_verified=True,train_executed=False,calibration_duration=None)
    path=tmp_path/'unavailable.json';collection=dict(native_missing_or_failed=False,conditions=[row])
    path.write_text(json.dumps(collection));assert numerical._unavailable_train_sources(path)==[row]
    collection['native_missing_or_failed']=True;path.write_text(json.dumps(collection))
    with pytest.raises(ValueError,match='Missing or failed execution'):
        numerical._unavailable_train_sources(path)


def test_final_fate_agreement_cannot_hide_a_gap_sign_switch():
    def measured(values):
        return dict(series=dict(time=[0.,1.,2.],D_half=values,D_half_raw=values,
            D_half_eligible=[True]*3,share_half=[(v+1)/2 for v in values],
            contrast_fraction_of_basal=[1.]*3,amount_fraction_of_basal=[1.]*3,
            ineligibility_reason=['']*3))
    reference=measured([-.9,-.6,.9]);companion=measured([-.9,.2,.901])
    result=numerical.phase_endpoint_comparison(reference,companion,
        [dict(name='gap',stop=1.),dict(name='recovery',stop=2.)])
    assert result['maximum_absolute_eligible_D_difference']==pytest.approx(.8)
    assert not result['all_available_challenger_bias_signs_agree']
    assert not result['all_eligible_endpoint_D_within_tolerance']
