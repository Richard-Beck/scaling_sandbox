import numpy as np

from migration1d.single_species.wgd import (
    parser,
    resting_state,
    scaled_parameters,
    transition_class,
)
from migration1d.single_species.wgd_extended import (
    all_tasks as extended_tasks,
    domain_transition,
    parser as extended_parser,
    reversal_transition,
)


def test_wgd_null_scales_geometry_but_preserves_resting_concentrations():
    args = parser().parse_args(["--mode", "shard"])
    baseline_mu, baseline_lambda, baseline_h = 3.0, 40.0, 2.0
    baseline, u0, v0, scale = scaled_parameters(
        args, baseline_mu, baseline_lambda, baseline_h, 2.0
    )

    assert np.isclose(scale, 2.0 ** (1.0 / 3.0))
    assert np.isclose(baseline.domain_length, scale * args.length)
    assert np.isclose(baseline.dimensionless_length, scale * baseline_lambda)
    assert np.isclose(baseline.cytoplasm_depth, scale * baseline_h)
    assert np.isclose(baseline.total_abundance, u0 + baseline.cytoplasm_depth * v0)
    np.testing.assert_allclose(resting_state(baseline), (u0, v0), rtol=1e-8)


def test_transition_class_highlights_polarity_changes():
    assert transition_class("single_polarized", "dual_polarized") == "single_to_dual"
    assert transition_class("dual_polarized", "single_polarized") == "dual_to_single"
    assert transition_class("homogeneous_low", "single_polarized") == "gained_polarity"
    assert transition_class("single_polarized", "homogeneous_high") == "lost_polarity"
    assert transition_class("homogeneous_low", "homogeneous_low") == "unchanged"


def test_extended_wgd_screen_has_one_array_task_per_simulation():
    args = extended_parser().parse_args(["--mode", "task"])
    tasks = extended_tasks(args)
    assert len(tasks) == 8685
    assert sum(task[0] == "dual" for task in tasks) == 6084
    assert sum(task[0] == "reversal" for task in tasks) == 2601


def test_extended_transitions_highlight_requested_changes():
    assert domain_transition(1, 2) == "one_to_two"
    assert domain_transition(2, 1) == "two_to_one"
    assert domain_transition(0, 2) == "gained_regions"
    assert reversal_transition("retained_original_pole", "reversed") == "became_reversible"
    assert reversal_transition("reversed", "lost_polarity") == "lost_reversibility"
