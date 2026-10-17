"""Force closure, fixed-site memory, and symmetry checks for substrate clutches."""
import numpy as np
import pytest
from polarity1d.debelly import DeBelly, DeBellyParameters, simulate_debelly, _dirichlet_elliptic
from polarity1d.protocols import LocalizedPulse
from polarity1d.substrate_clutch import FixedSubstrateClutch, FixedSubstrateParameters


def test_site_memory_stays_on_substrate_and_new_sites_are_unbound():
    clutch = FixedSubstrateClutch(FixedSubstrateParameters(spacing=0.05))
    clutch.initialize(np.linspace(0, 2, 11))
    clutch.n[:] = np.linspace(0.1, 0.9, len(clutch.n))
    clutch.q[:] = 0.2 * clutch.n
    ids, _, _, n, q = clutch._contact(np.linspace(0.5, 2.5, 11))
    overlap = ids <= clutch.ids[-1]
    np.testing.assert_allclose(n[overlap], clutch.n[10:])
    np.testing.assert_allclose(q[overlap], clutch.q[10:])
    assert np.all(n[~overlap] == 0) and np.all(q[~overlap] == 0)


def test_turnover_releases_extension_and_new_bonds_are_unstretched():
    clutch = FixedSubstrateClutch(FixedSubstrateParameters(spacing=0.05))
    edges = np.linspace(0, 2, 11)
    clutch.initialize(edges)
    clutch.n[:] = 0.2
    clutch.q[:] = 0.04 * np.where(np.arange(len(clutch.n)) < len(clutch.n)//2, -1, 1)
    old_q = clutch.q.copy()
    dt = 0.1
    off = clutch.parameters.off_rate * np.exp(0.2 / clutch.parameters.rupture_force)
    # A clamped, stationary cortex isolates the population update.
    def stationary(diagonal, viscosity, rhs, dx, left, right):
        return np.zeros_like(rhs)
    clutch.advance(edges=edges, dx=0.2, dt=dt, diagonal=np.ones(10), viscosity=1,
                   rhs=np.zeros(10), boundary=(0, 0), solve=stationary)
    np.testing.assert_allclose(clutch.q, old_q / (1 + dt*off), atol=1e-14)
    np.testing.assert_allclose(clutch.n, (0.2+dt*0.1)/(1+dt*(0.1+off)))
    assert np.all(np.abs(clutch.q/clutch.n) < 0.2)


def test_coupled_force_balance_and_velocity_frame():
    clutch = FixedSubstrateClutch()
    edges = np.linspace(0, 4, 21)
    clutch.initialize(edges)
    diagonal = np.full(20, 8.)
    rhs = np.linspace(-2, 3, 20)
    v = clutch.advance(edges=edges, dx=0.2, dt=0.25, diagonal=diagonal,
                       viscosity=10, rhs=rhs, boundary=(-0.1, 0.05), solve=_dirichlet_elliptic)
    traction = clutch.fields['clutch_traction']
    np.testing.assert_allclose(v, _dirichlet_elliptic(diagonal, 10, rhs-traction, 0.2, -0.1, 0.05), atol=1e-12)
    assert abs(np.sum(traction)*0.2) < 1e-12
    assert 0 <= clutch.n.min() <= clutch.n.max() <= 1


@pytest.mark.parametrize('centered',[True,False])
def test_unstimulated_cell_does_not_move_and_mirrored_cue_reverses_motion(centered):
    model = DeBelly(DeBellyParameters(assay_mode='reflection_symmetric_extension',center_boundary_motion=centered))
    kwargs = dict(length=6, cells=30, end_time=12, step=0.1, output_interval=2)
    rest = simulate_debelly(model, substrate_clutch=FixedSubstrateClutch(), **kwargs)
    np.testing.assert_allclose(rest.fields['cell_translation'], 0, atol=1e-12)
    trajectories = []
    for side in ('left', 'right'):
        trajectories.append(simulate_debelly(model, substrate_clutch=FixedSubstrateClutch(),
            protocol=LocalizedPulse(amplitude=3.6, start=1, stop=9, taper=2, width_fraction=.2, sides=(side,)), **kwargs))
    left, right = trajectories
    assert np.max(np.abs(left.fields['cell_translation'])) > 1e-5
    np.testing.assert_allclose(left.fields['cell_translation'], -right.fields['cell_translation'], atol=1e-9)
    np.testing.assert_allclose(left.fields['rac'], right.fields['rac'][:, ::-1], atol=1e-9)
    np.testing.assert_allclose(left.fields['cortical_velocity_lab'],
        left.fields['velocity'] + left.fields['cell_velocity'])
    assert left.metadata['mca_conservation_error'] < 1e-12
    if not centered:
        np.testing.assert_allclose(left.fields['right_boundary_displacement'],0,atol=1e-12)
        assert left.fields['left_boundary_displacement'].min()<0
