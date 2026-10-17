from __future__ import annotations

import numpy as np

from bulk_surface import FixedGeometry, SimulationConfig
from bulk_surface.holmes import HolmesParameters, holmes_direct_case


def test_ellipse_geometry_preserves_measures_and_constants():
    geometry = FixedGeometry.ellipse(3.0, 1.5, dx=0.25, n_surface=96)
    assert np.isclose(geometry.area, np.pi * 3.0 * 1.5, rtol=1e-12)
    assert np.allclose(
        geometry.trace_bulk(np.full(geometry.n_bulk, 2.75)), 2.75, atol=1e-14
    )
    assert np.max(np.abs(geometry.bulk_stiffness @ np.ones(geometry.n_bulk))) < 1e-12
    assert np.max(
        np.abs(geometry.surface_stiffness @ np.ones(geometry.n_surface))
    ) < 1e-12


def test_holmes_direct_model_has_nine_fields_and_conserves_each_family():
    geometry = FixedGeometry.ellipse(2.5, 1.25, dx=0.4, n_surface=48)
    parameters = HolmesParameters(dissociation_rate=0.5)
    case = holmes_direct_case(geometry, parameters=parameters)
    assert case.initial.bulk.shape == (3, geometry.n_bulk)
    assert case.initial.surface.shape == (6, geometry.n_surface)
    result = case.simulator.run(
        case.initial,
        SimulationConfig(end_time=0.1, dt=0.002, output_interval=0.02),
    )
    for family in ("cdc42", "rac", "rho"):
        total = result.pool_totals[f"{family}_total"]
        assert np.max(np.abs(total / total[0] - 1.0)) < 2e-11
    assert np.min(result.bulk) >= 0.0
    assert np.min(result.surface) >= 0.0


def test_holmes_cartesian_stimulus_is_stronger_at_positive_x():
    geometry = FixedGeometry.ellipse(2.5, 1.25, dx=0.4, n_surface=48)
    case = holmes_direct_case(geometry)
    stimulus = case.simulator.reactions.stimulus.values_at_membrane
    assert stimulus[np.argmax(geometry.membrane_points[:, 0])] == np.max(stimulus)
    assert stimulus[np.argmin(geometry.membrane_points[:, 0])] == np.min(stimulus)
    assert np.isclose(np.min(stimulus), 0.2)
    assert np.isclose(np.max(stimulus), 0.4)
