import numpy as np
from polarity1d.debelly import DeBelly, DeBellyParameters, simulate_debelly
from polarity1d.protocols import LocalizedPulse


def test_independent_endpoints_recover_one_sided_case_and_reflect():
    kwargs=dict(length=21,cells=105,end_time=80,step=.5,output_interval=2)
    results=[]
    for mode,side in [('published_one_sided','left'),
                      ('reflection_symmetric_extension','left'),
                      ('reflection_symmetric_extension','right')]:
        results.append(simulate_debelly(
            DeBelly(DeBellyParameters(assay_mode=mode,center_boundary_motion=False)),
            protocol=LocalizedPulse(amplitude=3.6,start=10,stop=50,taper=20,
                                    width_fraction=.2,sides=(side,)),**kwargs))
    original,left,right=results
    for field in original.fields:
        np.testing.assert_allclose(left.fields[field],original.fields[field],atol=1e-12)
    assert left.fields['left_boundary_displacement'].min()<-.1
    np.testing.assert_allclose(left.fields['right_boundary_displacement'],0,atol=1e-12)
    assert right.fields['right_boundary_displacement'].max()>.1
    np.testing.assert_allclose(right.fields['left_boundary_displacement'],0,atol=1e-12)
    np.testing.assert_allclose(left.fields['rac'],right.fields['rac'][:,::-1],atol=1e-10)
    np.testing.assert_allclose(left.fields['displacement'],-right.fields['displacement'][:,::-1],atol=1e-10)
    np.testing.assert_allclose(left.fields['velocity'],-right.fields['velocity'][:,::-1],atol=1e-10)
