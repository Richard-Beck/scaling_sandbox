import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from gradient_reversal import LinearGradient


def test_physical_slope_survives_length_changes_and_growth():
    cue=LinearGradient(.3,background=.1)
    for length in [1.,2.7,15.,39.]:
        x=np.linspace(0,length,201);y=cue(1,x,length)
        np.testing.assert_allclose(np.diff(y)/np.diff(x),.3)
        assert y[0]==.1
        assert np.isclose(y[-1],.1+.3*length)
        assert np.isclose(np.trapezoid(y-.1,x),.3*length**2/2)


def test_withdrawal_retains_background_and_exact_cutoff():
    cue=LinearGradient(.2,background=.1,stop=5.)
    x=np.linspace(0,15,201)
    assert cue.event_times==(5.,)
    assert cue(4.999,x,15)[-1]>.1
    np.testing.assert_array_equal(cue(5.,x,15),np.full_like(x,.1))
    np.testing.assert_array_equal(cue(6.,x,15),np.full_like(x,.1))


@pytest.mark.parametrize('slope',[0.,-1.,np.nan,np.inf])
def test_invalid_slope_rejected(slope):
    with pytest.raises(ValueError):LinearGradient(slope)
