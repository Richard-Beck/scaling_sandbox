from pathlib import Path
import numpy as np
import pytest
import native as r
import polarity_establishment as e


def test_basal_legi_state_is_an_equilibrium():
    case=dict(model='legi',length=15.,cells=31,background_stimulus=.1)
    state=e.resting_state(case)
    model=r.model_for('legi')
    grid=r.Grid1D(15.,31)
    derivative=model.reaction(0.,grid.x,state,e.withdrawal(dict(case,width=3.))(0.,grid.x,15.))
    np.testing.assert_allclose(state,np.tile([[.1],[.1],[.5]],(1,31)),atol=1e-15)
    np.testing.assert_allclose(derivative,0.,atol=1e-15)


@pytest.mark.parametrize('length',[3.,15.,45.])
@pytest.mark.parametrize('mode',['physical','relative'])
def test_establishment_mirror_width_and_basal_withdrawal(length,mode):
    x=np.linspace(0,length,10001)
    left=e.Patch(.4,3.,width_mode=mode,fraction=.2,background=.1,stop=10.)
    right=e.Patch(.4,3.,side='right',width_mode=mode,fraction=.2,background=.1,stop=10.)
    np.testing.assert_allclose(left(0,x,length),right(0,length-x,length),atol=1e-15)
    np.testing.assert_allclose(left(10,x,length),.1,atol=1e-15)
    np.testing.assert_allclose(left(0,np.array([0.]),length),.5,atol=1e-15)


@pytest.mark.parametrize('length',[1.,3.,15.,45.])
def test_establishment_patch_integral_with_truncation(length):
    x=np.linspace(0,length,50001)
    width=3.;amplitude=.1;background=.1
    pulse=e.Patch(amplitude,width,background=background)
    area=.5*(min(length,width)+width/np.pi*np.sin(np.pi*min(length,width)/width))
    np.testing.assert_allclose(np.trapezoid(pulse(0.,x,length)-background,x),amplitude*area,rtol=1e-8)


def test_guard_protects_all_retained_data():
    for n in ['linear_v11','width_reversal','establishment_basal','reversal_fields']:
        with pytest.raises(ValueError):r.protect_retained_output(r.DATA/n/'work')
    r.protect_retained_output(r.DATA/'runs/fresh')


def test_runtime_imports_are_inside_core():
    import polarity1d
    assert Path(polarity1d.__file__).resolve().is_relative_to(r.CODE)
    for model in [r.wave_pinning,r.otsuji,r.goryachev,r.legi,r.db,r.sp,r.ho]:
        assert Path(model.__file__).resolve().is_relative_to(r.CODE)
