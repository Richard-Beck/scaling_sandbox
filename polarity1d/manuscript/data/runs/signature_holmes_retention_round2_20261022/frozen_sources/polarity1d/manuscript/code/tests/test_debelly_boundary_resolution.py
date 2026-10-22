"""Boundary-physics regression checks for the separate corrected adapter."""
from pathlib import Path
import sys
from dataclasses import replace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mechanisms'))
from debelly_boundary_resolution import model,n,base,run

def asymmetric(side,L=39.):
 state=base(L);state['rac'][0 if side=='left' else -1]*=2
 return state

def test_inactive_endpoint_is_not_commanded_by_opposite_rac():
 tr=run('independent_endpoints',39.,asymmetric('right'),n.ConstantInput(0),.5,interval=.5)
 np.testing.assert_allclose(tr.fields['left_boundary_displacement'][-1],0,atol=1e-12)
 assert tr.fields['right_boundary_displacement'][-1,0]>.07

def test_independent_extension_recovers_published_single_active_end():
 L=39.;state=asymmetric('left');corrected=run('independent_endpoints',L,state,n.ConstantInput(0),20.,interval=1)
 published=n.ds.simulate_debelly(n.ds.DeBelly(replace(model('independent_endpoints').parameters,assay_mode='published_one_sided')),protocol=n.ConstantInput(0),length=L,cells=201,end_time=20,output_interval=1,step=.5,initial_checkpoint=state)
 for k in published.fields:np.testing.assert_allclose(corrected.fields[k],published.fields[k],atol=1e-11,rtol=1e-10)

def test_independent_extension_is_reflection_equivariant():
 a=run('independent_endpoints',39.,asymmetric('left'),n.ConstantInput(0),20.,interval=1)
 b=run('independent_endpoints',39.,asymmetric('right'),n.ConstantInput(0),20.,interval=1)
 for k in ['rac','rho','mca_bound','mca_unbound','tension']:np.testing.assert_allclose(a.fields[k],b.fields[k][:,::-1],atol=1e-11,rtol=1e-10)
 for k in ['displacement','velocity']:np.testing.assert_allclose(a.fields[k],-b.fields[k][:,::-1],atol=1e-11,rtol=1e-10)

def test_legacy_centering_changes_internal_mechanics_not_only_origin():
 state=asymmetric('right');a=run('legacy_centered',39.,state,n.ConstantInput(0),.5,interval=.5);b=run('independent_endpoints',39.,state,n.ConstantInput(0),.5,interval=.5)
 np.testing.assert_allclose(a.fields['boundary_separation'],b.fields['boundary_separation'],atol=1e-12)
 assert np.max(abs(a.fields['tension'][-1]-b.fields['tension'][-1]))>1
 shift=.5*(b.fields['left_boundary_displacement'][-1,0]+b.fields['right_boundary_displacement'][-1,0]);assert np.max(abs(a.fields['displacement'][-1]-(b.fields['displacement'][-1]-shift)))>.03

def test_display_centering_preserves_strain_without_changing_dynamics():
 tr=run('independent_endpoints',39.,asymmetric('right'),n.ConstantInput(0),20.,interval=1)
 u=tr.fields['displacement'];c=.5*(tr.fields['left_boundary_displacement']+tr.fields['right_boundary_displacement']);display=u-c
 np.testing.assert_allclose(np.diff(u,axis=1),np.diff(display,axis=1),atol=1e-14)
