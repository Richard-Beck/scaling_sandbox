import json
import numpy as np
import corrected_reversal as c


def test_refinement_covers_every_sampled_fail_to_pass_transition():
    trials=[dict(duration=d,orientation_success=gate,full_polarity_success=gate) for d,gate in [(0.,False),(1.,True),(2.,False),(3.,True),(4.,False)]]
    extra=c.refinement(trials)
    assert .5 in extra and 2.5 in extra
    assert all(d>0 for d in extra)
    assert 3.5 not in extra


def test_resolved_bracket_needs_no_more_refinement():
    trials=[dict(duration=d,orientation_success=gate,full_polarity_success=gate) for d,gate in [(0.,False),(99.5,False),(100.,True)]]
    assert c.refinement(trials)==[]


def test_duration_zero_cannot_supply_a_response_endpoint():
    case=dict(case_index=0,establishment_case_index=0,model='wave_pinning',length=15.,condition='fixed_20pct',width_mode='physical',width=3.,amplitude=.1,initial_orientation=-.4,eligible=True)
    # A zero control that reverses is an invalid preparation, not instant response.
    trial=dict(duration=0.,orientation_success=True,full_polarity_success=True,no_cue_preparation_valid=False)
    rows,records=c.rows_and_records(dict(cases=[case]),{0:[trial]},True)
    assert all(row['status']=='preparation_invalid' and row['endpoint'] is None for row in rows)


def test_frozen_duration_grid_includes_the_corrected_full_bounds():
    data=c.read(c.GRID)['durations_by_model']
    for model in c.MODELS:
        ds=data[model]
        assert ds[0]==0. and np.all(np.diff(ds)>0)
        hi=50. if model=='spring' else 2000. if model=='debelly' else 5120.
        assert ds[-1]==hi
