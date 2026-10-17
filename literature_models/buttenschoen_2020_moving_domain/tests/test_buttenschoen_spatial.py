from pathlib import Path
import sys

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))

from buttenschoen import (  # noqa: E402
    InitialCondition,
    MechanicalParameters,
    SimulationParameters,
    SimulationResult,
    figure9_regimes,
    initial_state,
    scale_wgd_condition,
    simulate,
)
from buttenschoen.model import endpoint_forces  # noqa: E402
from buttenschoen.render import (  # noqa: E402
    plot_kymograph,
    plot_length_centroid_comparison,
    plot_soft_protrusive_mechanism,
    plot_wgd_companion,
)


def test_force_sign_flips_between_protrusive_and_contractile():
    active = np.asarray([0.1, 0.2, 2.0])
    protrusive = endpoint_forces(active, MechanicalParameters(coupling="protrusive"))
    contractile = endpoint_forces(active, MechanicalParameters(coupling="contractile"))
    assert protrusive[0] < 0 < protrusive[1]
    assert contractile == pytest.approx(tuple(-value for value in protrusive))


def test_initial_state_has_requested_physical_amount():
    mechanics = MechanicalParameters(rest_length=50.0)
    initial = InitialCondition(total_amount=10.0)
    _, state = initial_state(mechanics, initial, 101)
    radius = state[-1]
    concentration_sum = state[:101] + state[101:202]
    amount = radius*(2.0/101)*np.sum(concentration_sum)
    assert amount == pytest.approx(initial.total_amount)


@pytest.mark.parametrize("coupling", ["protrusive", "contractile"])
def test_mass_is_conserved_while_cell_moves(coupling):
    result = simulate(
        mechanics=MechanicalParameters(coupling=coupling, elastic_modulus=500.0),
        initial=InitialCondition(total_amount=25.0),
        simulation=SimulationParameters(n_cells=41, end_time=20.0,
                                        output_interval=1.0,
                                        rtol=1e-8, atol=1e-10),
    )
    assert np.ptp(result.length) > 1e-4
    assert np.max(result.relative_mass_error) < 2e-7


def test_figure9_presets_are_explicit_and_distinct():
    regimes = figure9_regimes()
    assert set(regimes) == {
        "contractile_stiff", "contractile_soft",
        "protrusive_stiff", "protrusive_soft",
    }
    assert regimes["contractile_stiff"][0].elastic_modulus == 1500
    assert regimes["contractile_soft"][1].total_amount == 10
    assert regimes["protrusive_soft"][0].elastic_modulus == 250


def test_kymograph_renderer(tmp_path):
    result = simulate(
        simulation=SimulationParameters(n_cells=21, end_time=2.0,
                                        output_interval=1.0),
    )
    output = plot_kymograph(result, tmp_path / "kymograph.png")
    assert output.is_file()
    assert output.stat().st_size > 1000


def test_saved_metadata_reports_conservation_and_endpoint_metrics(tmp_path):
    result = simulate(
        simulation=SimulationParameters(n_cells=21, end_time=2.0,
                                        output_interval=1.0),
    )
    metadata = result.save(tmp_path)
    import json
    saved = json.loads(metadata.read_text())
    assert saved["maximum_relative_mass_error"] < 1e-5
    assert saved["minimum_length"] <= saved["maximum_length"]
    assert 0 <= saved["final_active_contrast"] <= 1


def test_saved_result_round_trips(tmp_path):
    result = simulate(
        simulation=SimulationParameters(n_cells=21, end_time=2.0,
                                        output_interval=1.0),
    )
    result.save(tmp_path)
    restored = SimulationResult.load(tmp_path)
    assert restored.mechanics == result.mechanics
    assert restored.initial_condition == result.initial_condition
    assert restored.simulation == result.simulation
    assert np.array_equal(restored.times, result.times)
    assert np.allclose(restored.active, result.active)
    assert np.allclose(restored.inactive, result.inactive)
    assert np.allclose(restored.centroid, result.centroid)
    assert np.allclose(restored.radius, result.radius)


def test_wgd_scaling_laws_and_initial_concentration_are_matched():
    mechanics = MechanicalParameters(
        rest_length=50, elastic_modulus=250, cell_viscosity=60,
        endpoint_friction=100000, maximum_force=50000,
    )
    initial = InitialCondition(total_amount=10, centroid=0)
    scaled_mechanics, scaled_initial = scale_wgd_condition(mechanics, initial, 2)
    scale = 2**(1/3)
    assert scaled_mechanics.rest_length == pytest.approx(scale*50)
    assert scaled_mechanics.elastic_modulus == pytest.approx(scale*250)
    assert scaled_mechanics.cell_viscosity == pytest.approx(scale*60)
    assert scaled_mechanics.endpoint_friction == pytest.approx(scale**2*100000)
    assert scaled_mechanics.maximum_force == pytest.approx(scale**2*50000)
    assert scaled_initial.total_amount == pytest.approx(scale*10)
    assert scaled_initial.centroid == pytest.approx(0)
    assert (scaled_initial.total_amount/scaled_mechanics.rest_length
            == pytest.approx(initial.total_amount/mechanics.rest_length))


def test_all_figure9_conditions_begin_at_centroid_zero():
    for mechanics, initial in figure9_regimes().values():
        assert initial.centroid == 0
        scaled_mechanics, scaled_initial = scale_wgd_condition(
            mechanics, initial, 2
        )
        assert scaled_initial.centroid == 0
        _, baseline_state = initial_state(mechanics, initial, 21)
        _, scaled_state = initial_state(scaled_mechanics, scaled_initial, 21)
        assert baseline_state[-2] == 0
        assert scaled_state[-2] == 0


def test_wgd_pair_conserves_mass_under_motion():
    mechanics, initial = figure9_regimes()["protrusive_soft"]
    mechanics, initial = scale_wgd_condition(mechanics, initial, 2)
    result = simulate(
        mechanics=mechanics,
        initial=initial,
        simulation=SimulationParameters(n_cells=31, end_time=10,
                                        output_interval=1,
                                        rtol=1e-8, atol=1e-10),
    )
    assert np.ptp(result.length) > 1e-4
    assert np.max(result.relative_mass_error) < 2e-7


def test_wgd_companion_renderer(tmp_path):
    simulation = SimulationParameters(n_cells=15, end_time=1,
                                      output_interval=1)
    results = {}
    for name, (mechanics, initial) in figure9_regimes().items():
        scaled = scale_wgd_condition(mechanics, initial, 2)
        results[name] = {
            "1x": simulate(mechanics=mechanics, initial=initial,
                           simulation=simulation),
            "2x": simulate(mechanics=scaled[0], initial=scaled[1],
                           simulation=simulation),
        }
    output = plot_wgd_companion(results, tmp_path/"companion.png")
    assert output.is_file()
    assert output.stat().st_size > 1000

    mechanism = plot_soft_protrusive_mechanism(
        results["protrusive_soft"], tmp_path/"mechanism.png"
    )
    comparison = plot_length_centroid_comparison(
        results, tmp_path/"length_centroid.png"
    )
    assert mechanism.is_file()
    assert mechanism.stat().st_size > 1000
    assert comparison.is_file()
    assert comparison.stat().st_size > 1000
