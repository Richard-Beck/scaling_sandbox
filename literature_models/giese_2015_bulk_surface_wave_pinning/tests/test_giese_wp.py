from __future__ import annotations

import numpy as np

from giese_wp import (GieseSolver, Stimulus, half_mass_pf,
                      make_circular_system, make_gmsh_circular_system,
                      membrane_metrics)


def test_circle_and_membrane_measures_converge_to_geometry():
    system = make_circular_system(3.0, target_h=0.25)
    assert abs(system.area - np.pi * 1.5**2) / (np.pi * 1.5**2) < 0.01
    assert abs(system.perimeter - 2 * np.pi * 1.5) / (2 * np.pi * 1.5) < 0.01
    assert system.max_edge_length <= 0.25


def test_shared_reaction_flux_conserves_mass_and_transfers_oppositely():
    system = make_circular_system(3.0, target_h=0.5)
    solver = GieseSolver(system)
    u0, v0 = solver.initial_state()
    dt = 0.01
    reaction = solver.reaction_load(u0, v0, 0.0, (Stimulus(0.44, 0.0),))
    u1, v1, _ = solver.step(u0, v0, 0.0, dt, (Stimulus(0.44, 0.0),))
    membrane_change = np.ones(u1.size) @ system.membrane_mass @ (u1 - u0)
    bulk_change = np.ones(v1.size) @ system.bulk_mass @ (v1 - v0)
    assert reaction.sum() > 0.0
    assert np.isclose(membrane_change, dt * reaction.sum(), atol=1e-13)
    assert np.isclose(bulk_change, -dt * reaction.sum(), atol=1e-13)
    assert np.isclose(solver.total_mass(u1, v1), solver.total_mass(u0, v0), atol=1e-13)


def test_stimulus_has_exact_physical_membrane_coverage():
    for diameter in (1.5, 3.0, 4.5):
        system = make_circular_system(diameter, target_h=0.25)
        solver = GieseSolver(system)
        _, v = solver.initial_state()
        v[:] = 1.0
        stimulus = Stimulus(1.0, 0.0, membrane_fraction=0.05)
        load = solver.stimulus_reaction_load(v, stimulus, 0.0)
        assert np.isclose(load.sum() / system.perimeter, 0.05, atol=1e-14)


def test_no_stimulus_stays_spatially_homogeneous_and_conserves_mass():
    solver = GieseSolver(make_circular_system(3.0, target_h=0.5))
    result = solver.run(20.0, dt=0.1, output_interval=1.0)
    assert max(row["relative_amplitude"] for row in result.metrics) < 1e-4
    assert np.max(np.abs(result.total_mass / result.total_mass[0] - 1.0)) < 1e-11


def test_published_pf_uses_component_around_global_maximum():
    u = np.array([1.0, 3.0, 3.0, 1.0, 1.0, 2.5, 1.0, 1.0])
    metrics = membrane_metrics(u, np.ones_like(u))
    # Global peak component occupies two of eight equal boundary edges/vertices.
    assert np.isclose(metrics["PF"], 0.75)
    assert metrics["n_above_mean_clusters"] == 2


def test_half_mass_pf_is_zero_for_uniform_and_positive_for_localized_mass():
    edges = np.ones(64)
    assert abs(half_mass_pf(np.ones(64), edges)) < 1e-3
    localized = np.ones(64)
    localized[:8] = 20.0
    assert half_mass_pf(localized, edges) > 0.7


def test_gmsh_mesh_is_irregular_and_refined_at_membrane():
    system = make_gmsh_circular_system(3.0, boundary_h=0.08, bulk_h=0.24)
    facets = system.mesh.facets
    lengths = np.linalg.norm(
        system.mesh.p[:, facets[1]] - system.mesh.p[:, facets[0]], axis=0
    )
    boundary_lengths = lengths[system.mesh.boundary_facets()]
    assert lengths.max() > 2.0 * boundary_lengths.mean()
    assert np.std(lengths) > 0.02
    assert abs(system.area - np.pi * 1.5**2) / (np.pi * 1.5**2) < 0.005
