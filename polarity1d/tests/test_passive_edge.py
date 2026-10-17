import numpy as np
from polarity1d.passive_edge import solve_passive_edge
from polarity1d.debelly import DeBelly,DeBellyParameters,simulate_debelly
from polarity1d.substrate_clutch import FixedSubstrateClutch
from polarity1d.protocols import LocalizedPulse


def test_loaded_edge_moves_and_both_boundary_conditions_hold():
    u,v,speed,ub,vb,res=solve_passive_edge(
        old_u=np.zeros(30),old_boundary=0.,active_displacement=-.02,
        active_velocity=-.05,passive_side='right',friction=np.full(30,8.),
        stiffness=100.,viscosity=10000.,active_stress=100.,rho=np.full(30,.69),
        resting_rho=.69,dx=.2,dt=.1,memory=np.zeros(30),drag=np.full(30,.5))
    assert ub[0]==-.02 and vb[0]==-.05
    assert abs(ub[1])>1e-6 and abs(vb[1])>1e-6
    np.testing.assert_allclose(vb[1],ub[1]/.1,atol=1e-10)
    np.testing.assert_allclose(100*2*(ub[1]-u[-1])/.2+10000*2*(vb[1]-v[-1])/.2,0,atol=1e-7)
    assert abs(np.sum(.5*(v+speed))*.2)<1e-10
    assert abs(res['global_force_residual'])<1e-7
    # Independently evaluate both interior force balances, including face BCs.
    def laplacian(field,boundary):
        extended=np.r_[2*boundary[0]-field[0],field,2*boundary[1]-field[-1]]
        return np.diff(extended,n=2)/.2**2
    np.testing.assert_allclose(8*(u/.1-v),100*laplacian(u,ub),atol=1e-8)
    np.testing.assert_allclose(10000*laplacian(v,vb)-8*(v-u/.1)-.5*(v+speed),0,atol=1e-8)


def test_passive_rest_and_reflection():
    common=dict(length=6,cells=30,end_time=20,step=.1,output_interval=1)
    def model(side):
        return DeBelly(DeBellyParameters(assay_mode='reflection_symmetric_extension',
                        center_boundary_motion=False,passive_boundary=side))
    rest=simulate_debelly(model('right'),substrate_clutch=FixedSubstrateClutch(),**common)
    for name in ('displacement','velocity','cell_translation'):
        np.testing.assert_allclose(rest.fields[name],0,atol=1e-10)
    runs=[]
    for active,passive in [('left','right'),('right','left')]:
        runs.append(simulate_debelly(model(passive),substrate_clutch=FixedSubstrateClutch(),
            protocol=LocalizedPulse(amplitude=3.6,start=1,stop=15,taper=3,width_fraction=.2,sides=(active,)),**common))
    left,right=runs
    for name in ('displacement','velocity','cell_translation'):
        np.testing.assert_allclose(left.fields[name],-right.fields[name][:,::-1],atol=1e-8)
    np.testing.assert_allclose(left.fields['rac'],right.fields['rac'][:,::-1],atol=1e-9)
    assert abs(left.fields['right_boundary_displacement']).max()>1e-4
    assert left.metadata['maximum_boundary_residuals']['passive_stress_residual']<1e-7
    assert left.metadata['maximum_boundary_residuals']['global_force_residual']<1e-7
