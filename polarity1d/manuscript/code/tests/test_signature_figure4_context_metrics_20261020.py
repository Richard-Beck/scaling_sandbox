from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
import signature_figure4_context_metrics_20261020 as M


def rows():
    e1=[dict(model='model',C_label=i,f=.1,sigma=sigma,forced_combined=None) for i in range(24) for sigma in (.25,.5,1.,2.)]
    longest=[dict(model='model',C_label=i,f=.1,selected_observation_arm='preselected_E4',forced_combined=None) for i in range(52)]
    longest.append(dict(model='model',C_label=80,f=.1,selected_observation_arm='common_new_E3_size',forced_combined=False))
    return e1,longest


def test_contexts_keep_all_extra_sigmas_all52long_and_unknown_qualification():
    e1,longest=rows();result=M.context_rows([{'retained_common':True}],longest,e1)
    assert len(result['E1_extra'])==72 and len(result['E4_long'])==52
    assert {r['sigma'] for r in result['E1_extra']}=={.25,.5,2.}
    assert all(r['forced_combined'] is None for r in result['E1_extra']+result['E4_long'])
    assert result['common']==[{'retained_common':True}]


def test_contexts_reject_missing_preset_arms_and_prespecified_long_conditions():
    e1,longest=rows()
    with pytest.raises(ValueError,match='three preset'):M.context_rows([],longest,e1[1:])
    with pytest.raises(ValueError,match='every52'):M.context_rows([],longest[1:],e1)
