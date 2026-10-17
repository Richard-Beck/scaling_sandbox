import json

import numpy as np

from migration1d.focused_motile_length import (
    CenteredHillForce,
    all_tasks,
    make_model,
    mechanics_parameters,
    node_count,
    parser,
    run_task,
)
from migration1d.forces import EndpointState
from migration1d.model import endpoint_forces
from migration1d.tension import GlobalTensionFeedback, SpatialTensionFeedback


def arguments(*extra):
    return parser().parse_args([
        "--chemistry-class", "single_species",
        *extra,
    ])


def test_production_design_has_240_unique_tasks_per_chemistry():
    tasks = all_tasks(arguments())
    assert len(tasks) == 240
    assert len(set(tasks)) == len(tasks)
    assert {task[2] for task in tasks} == {"soft", "stiff"}
    assert {task[3] for task in tasks} == {
        "single", "dual_equal", "dual_biased", "reversal"
    }
    assert {task[4] for task in tasks} == {0.0, 5.0, 20.0}


def test_design_rejects_more_than_three_tension_strengths():
    with np.testing.assert_raises_regex(ValueError, "between one and three"):
        all_tasks(arguments("--tension-strengths", "0,1,2,3"))


def test_mechanical_cases_are_one_hundred_fold_apart_and_length_normalized():
    for length in (5.0, 20.0, 50.0):
        soft = mechanics_parameters(length, "soft")
        stiff = mechanics_parameters(length, "stiff")
        assert np.isclose(stiff.spring_stiffness / soft.spring_stiffness, 100.0)
        assert np.isclose(soft.spring_stiffness * length, 1.0)
        assert np.isclose(stiff.spring_stiffness * length, 100.0)
        assert np.isclose(soft.dashpot_viscosity * length, 20.0)


def test_node_rule_preserves_minimum_and_maximum_spacing():
    assert node_count(5.0, 81, 0.25) == 81
    assert node_count(20.0, 81, 0.25) == 81
    assert node_count(50.0, 81, 0.25) == 201


def test_single_species_force_is_balanced_at_rest():
    args = arguments("--lengths", "20")
    model = make_model(args, 20.0, "soft", "single")
    position = np.linspace(0.0, 20.0, 11)
    membrane = np.repeat(model.chemistry.initial_membrane[:, None], 11, axis=1)
    cytoplasm = np.repeat(model.chemistry.initial_cytoplasm[:, None], 11, axis=1)
    np.testing.assert_allclose(
        endpoint_forces(model, 0.0, position, membrane, cytoplasm), 0.0,
        atol=1e-14,
    )


def test_feedback_mode_constructs_requested_tension_model():
    global_model = make_model(arguments(), 20.0, "soft", "single", alpha_t=5.0)
    assert isinstance(global_model.tension_feedback, GlobalTensionFeedback)

    spatial_args = arguments(
        "--tension-mode", "spatial",
        "--spatial-tension-diffusivity", "10",
        "--spatial-tension-relaxation-time", "10",
        "--spatial-tension-generation-rate", "0.01",
    )
    spatial_model = make_model(
        spatial_args, 20.0, "soft", "single", alpha_t=20.0
    )
    feedback = spatial_model.tension_feedback
    assert isinstance(feedback, SpatialTensionFeedback)
    assert feedback.alpha_t == 20.0
    assert feedback.characteristic_length == 10.0


def test_centered_force_translates_toward_the_active_pole():
    force = CenteredHillForce(baseline=0.2, gain=1.0, half_level=1.0)
    position = np.linspace(0.0, 10.0, 11)
    membrane = np.full((1, 11), 0.2)
    membrane[0, -1] = 2.0
    inactive = np.ones((1, 11))
    states = [
        EndpointState(0.0, side, position, membrane, inactive, ("active_gtpase",))
        for side in ("left", "right")
    ]
    assert force(states[0]) == 0.0
    assert force(states[1]) > 0.0


def test_short_task_writes_valid_motile_artifacts(tmp_path, monkeypatch):
    monkeypatch.setenv("MIGRATION1D_OUTPUT_DIR", str(tmp_path))
    args = arguments(
        "--lengths", "5",
        "--mechanical-cases", "soft",
        "--protocols", "single",
        "--duration", "2",
        "--output-interval", "0.5",
        "--max-step", "0.1",
        "--minimum-nodes", "21",
        "--stimulus-start", "0.5",
        "--stimulus-stop", "1.5",
        "--reversal-start", "1.6",
        "--reversal-stop", "1.9",
    )
    json_path, trajectory_path = run_task(args)
    payload = json.loads(json_path.read_text())
    trajectory = np.load(trajectory_path)
    assert payload["task_key"] == ["single_species", 5.0, "soft", "single", 0.0]
    assert payload["row"]["minimum_length_um"] > 0.0
    assert trajectory["position"].shape == (5, 21)
    assert np.ptp(trajectory["center"]) > 0.0
    assert trajectory["elastic_tension"].shape == (5,)
    np.testing.assert_array_equal(trajectory["tension_added_deactivation_rate"], 0.0)
    assert payload["row"]["maximum_active_gtpase_conservation_error"] < 1e-7
