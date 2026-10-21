from pathlib import Path
import json
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_figure4_geometry_20261020 as G


def row(**changes):
    return dict(dict(arm='common',view='main',model='model',f=.01,forced_combined=False,evaluable=True,actual_geometry_ratio=9.,excess_persistence_ci_lower=1.1,excess_persistence_ci_upper=1.5,Dx_large_over_small=.6,Dx_ratio_ci_lower=.5,Dx_ratio_ci_upper=.8),**changes)


def test_geometry_preserves_conflicting_time_and_mobility_conclusions():
    r=G.classify(row(Dx_large_over_small=2.,Dx_ratio_ci_lower=1.5,Dx_ratio_ci_upper=2.5))
    assert r['time_geometry_conclusion']=='time_exceeds_geometry'
    assert r['mobility_geometry_conclusion']=='Dx_increases_with_size'
    assert r['geometric_qualifier']=='conflicting_time_and_mobility'


def test_unbounded_intervals_are_unavailable_and_CI_including_one_is_not_equivalence():
    r=G.classify(row(excess_persistence_ci_upper=float('inf'),Dx_ratio_ci_lower=.5,Dx_ratio_ci_upper=1.5))
    assert r['time_geometry_conclusion']=='unavailable'
    assert r['mobility_geometry_conclusion']=='Dx_CI_includes_flat'
    assert r['geometric_qualifier']=='compatible_with_flat_Dx'
    assert 'not equivalence' in r['interval_includes_one_interpretation']


def test_forced_and_unknown_qualifications_are_distinct_from_inclusive_sensitivity():
    assert G.classify(row(forced_combined=True))['geometric_qualification']=='excluded_forced'
    assert G.classify(row(forced_combined=None))['geometric_qualification']=='unavailable_forced_qualification'
    assert G.classify(row(forced_combined=None),True)['geometric_qualifier']=='more_than_geometry'


def test_finite_Dx_interval_does_not_override_unresolved_point_ratio():
    r=G.classify(row(Dx_large_over_small=None))
    assert r['mobility_geometry_conclusion']=='unavailable'


def test_model_summary_exposes_input_heterogeneity_without_selecting_one_strength():
    values=[G.classify(row()),G.classify(row(f=.1,excess_persistence_ci_lower=.2,excess_persistence_ci_upper=.8,Dx_large_over_small=2.,Dx_ratio_ci_lower=1.5,Dx_ratio_ci_upper=2.5))]
    m=G.model_rows(values)[0]
    assert m['per_model_pattern']=='heterogeneous_across_inputs'
    assert m['qualifier_input_strengths']=={'more_than_geometry':[.01],'less_than_geometry':[.1]}


def test_geometry_export_rejects_unrelated_passed_audit_and_binds_accepted_sources(tmp_path):
    parent=tmp_path/'accepted';parent.mkdir();G.write_csv(parent/'size_ratios.csv',[row()])
    complete=dict(passed=True,complete=True,outputs={'size_ratios.csv':G.digest(parent/'size_ratios.csv')},source_freeze_sha256='a'*64)
    (parent/'complete.json').write_text(json.dumps(complete));audit=parent/'independent_final_decision_audit.json'
    audit.write_text(json.dumps(dict(passed=True,producer_complete_sha256='unrelated')))
    with pytest.raises(ValueError,match='different producer'):G.export(parent,tmp_path/'geometry')
    audit.write_text(json.dumps(dict(passed=True,producer_complete_sha256=G.digest(parent/'complete.json'))))
    assert G.export(parent,tmp_path/'geometry')['complete']
    status=json.loads((tmp_path/'geometry/status.json').read_text())
    assert status['source_size_ratios_sha256']==G.digest(parent/'size_ratios.csv')
    assert status['source_independent_audit_sha256']==G.digest(audit)
