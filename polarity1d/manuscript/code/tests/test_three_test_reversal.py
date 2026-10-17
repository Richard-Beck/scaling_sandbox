"""Scientific input invariants: changing length cannot retune patch amplitude."""
import copy
import importlib.util
from pathlib import Path
import sys
import numpy as np

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE))
import three_test_reversal as t


def test_patch_peak_is_length_independent_and_support_obeys_geometry():
    base=next(c for c in t.base_cases() if c['model']=='goryachev' and c['eligible'])
    for geometry in t.TESTS[:2]:
        for length in [5.,15.,45.]:
            c=copy.deepcopy(base);c['length']=length
            c=t.make_case(c,geometry,64.)
            cue=t.response_cue(c)
            width=3. if geometry=='fixed_physical' else .2*length
            np.testing.assert_allclose(cue(0,np.array([length,length-width/2,length-width]),length),[6.4,3.2,0.],atol=1e-12)
            assert c['checkpoint']==base['checkpoint'] and c['checkpoint_sha256']==base['checkpoint_sha256']


def test_gradient_keeps_physical_slope_and_shared_checkpoint():
    base=next(c for c in t.base_cases() if c['model']=='goryachev' and c['eligible'])
    slopes=[]
    for length in [5.,15.,45.]:
        c=copy.deepcopy(base);c['length']=length;c=t.make_case(c,'fixed_slope')
        signal=t.response_cue(c)(0,np.array([0.,length]),length)
        slopes.append(signal[-1]/length)
        assert c['checkpoint_sha256']==base['checkpoint_sha256']
    np.testing.assert_allclose(slopes,[slopes[0]]*3)
