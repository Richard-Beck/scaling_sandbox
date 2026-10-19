"""Small independent algebraic checks; no trajectory integration or fitting."""
from __future__ import annotations

from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np

from spring_ring_closure_20261017 import Closure as OriginalClosure, native
from spring_ring_model_20261019 import (
    Closure, CuePhase, run_protocol, physical_raised_cosine,
    polarity_vector, _active_vector_history,
)


def check_model():
    # A small grid permits a full independent finite-difference Jacobian check,
    # including entries far from each concentration's nearest neighbours.
    closure = Closure(2.0, 24)
    original = OriginalClosure(2.0, 24)
    n = closure.cells
    rng = np.random.default_rng(20261019)
    state = np.r_[rng.uniform(0.65, 1.4, n), rng.uniform(0.3, 0.8, n), 2.043]
    constant_cue = rng.uniform(0.0, 4.0, n)
    np.testing.assert_array_equal(closure.initial(), original.initial())
    np.testing.assert_array_equal(closure.rhs(state, constant_cue), original.rhs(state, constant_cue))
    assert closure.amount(state) == original.amount(state)

    def cue_at(y):
        return physical_raised_cosine(closure, y, 0.97, 4.0, 0.31)

    def physical_rhs(y):
        return closure.rhs(y, cue_at(y)[0])

    comparisons = []
    for label, function, jacobian in (
        ("state_independent_profile", lambda y: closure.rhs(y, constant_cue),
         closure.jacobian(state, constant_cue)),
        ("fixed_physical_width_profile", physical_rhs,
         closure.jacobian(state, *cue_at(state))),
    ):
        numerical = np.empty_like(jacobian)
        # Several step sizes expose truncation/roundoff behaviour rather than
        # selecting agreement from only one finite-difference step.
        step_results = []
        for relative_step in (1e-4, 1e-5, 1e-6):
            for column in range(len(state)):
                step = relative_step * max(1.0, abs(float(state[column])))
                above, below = state.copy(), state.copy()
                above[column] += step
                below[column] -= step
                numerical[:, column] = (function(above) - function(below)) / (2.0*step)
            residual = np.abs(jacobian - numerical)
            step_results.append(dict(step=relative_step,
                max_absolute_error=float(residual.max()),
                max_scaled_error=float(np.max(residual / (1.0 + np.abs(jacobian))))))
        assert step_results[-1]["max_scaled_error"] < 5e-7
        comparisons.append(dict(profile=label, finite_difference=step_results))

    # Differentiate the physical-amount conservation identity independently.
    # grad(M).J + rhs.T.Hess(M) must vanish: M=C*mean(u+v).
    derivative = physical_rhs(state)
    jacobian = closure.jacobian(state, *cue_at(state))
    amount_gradient = np.r_[np.full(2*n, state[-1]/n), np.mean(state[:-1]) * 2.0]
    amount_hessian_rhs = np.r_[np.full(2*n, derivative[-1]/n), np.sum(derivative[:-1])/n]
    conservation_derivative = amount_gradient @ jacobian + amount_hessian_rhs
    assert np.max(np.abs(conservation_derivative)) < 1e-10
    assert abs(amount_gradient @ derivative) < 1e-10

    # A global coupling is physically present and cannot be omitted from a
    # sparse local stencil: u_0 alters dilution of distant inactive v_(n/2).
    global_entry = float(jacobian[n+n//2, 0])
    assert abs(global_entry) > 1e-8

    # Test full covariance of the derivative and Jacobian, including the cue
    # circumference derivative, rather than only a uniform concentration.
    shift = 7
    permutation = np.r_[np.roll(np.arange(n), shift), n+np.roll(np.arange(n), shift), 2*n]
    shifted = state[permutation]
    shifted_cue = physical_raised_cosine(closure, shifted, (0.97 + shift/n) % 1.0, 4.0, 0.31)
    rhs_rotation = closure.rhs(shifted, shifted_cue[0]) - derivative[permutation]
    jac_rotation = closure.jacobian(shifted, *shifted_cue) - jacobian[np.ix_(permutation, permutation)]
    assert np.max(np.abs(rhs_rotation)) < 1e-10
    assert np.max(np.abs(jac_rotation)) < 1e-10

    # Compare complete homogeneous derivatives with the native interval model.
    uniform = np.r_[np.full(n, 1.15), np.full(n, 0.72), 2.045]
    native_state = np.empty(2*n + 2)
    native_state[:-2:2], native_state[1:-2:2] = uniform[:n], uniform[n:-1]
    native_state[-2:] = 0.0, uniform[-1]
    native_derivative = native.rhs(0.0, native_state, L0=2.0)
    expected_native = np.r_[native_derivative[:-2:2], native_derivative[1:-2:2],
                            native_derivative[-1] - native_derivative[-2]]
    uniform_residual = float(np.max(np.abs(closure.rhs(uniform, 0.0) - expected_native)))
    assert uniform_residual < 1e-12

    # The first circular moment must distinguish amplitude from angle and
    # remain equivariant under cyclic changes in the chosen angular origin.
    vector = polarity_vector(closure, state)
    rotated_vector = polarity_vector(closure, shifted)
    expected_vector = complex(vector["x"], vector["y"]) * np.exp(2j*np.pi*shift/n)
    vector_residual = abs(complex(rotated_vector["x"], rotated_vector["y"]) - expected_vector)
    assert vector_residual < 1e-14
    assert polarity_vector(closure, uniform)["direction_radians"] is None
    with np.errstate(divide="raise", invalid="raise"):
        zero_active = uniform.copy()
        zero_active[:n] = 0.0
        assert polarity_vector(closure, zero_active) == dict(
            x=0.0, y=0.0, magnitude=0.0, direction_radians=None)
        zero_history_vector = _active_vector_history(closure, np.stack([zero_active[:n], state[:n]]))
        assert zero_history_vector[0] == 0.0j
        assert np.isfinite(zero_history_vector).all()
        assert abs(zero_history_vector[1] - complex(vector["x"], vector["y"])) < 1e-15

    # A seam-crossing physical cue is smooth at support edges and agrees with
    # direct periodic distances. Check its C derivative independently.
    signal, derivative_cue = cue_at(state)
    assert signal[0] > 0 and signal[-1] > 0
    h = 1e-6
    above, below = state.copy(), state.copy()
    above[-1] += h
    below[-1] -= h
    cue_derivative_error = float(np.max(np.abs((cue_at(above)[0] - cue_at(below)[0])/(2*h) - derivative_cue)))
    assert cue_derivative_error < 1e-8

    return dict(generated_at=datetime.now(timezone.utc).isoformat(), passed=True,
        trajectory_integrations=0, unchanged_rhs_exact=True,
        native_source_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        source_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), Path(__file__).with_name("spring_ring_model_20261019.py"),
                         Path(__file__).with_name("spring_ring_closure_20261017.py"))},
        jacobian_checks=comparisons,
        physical_amount_derivative_residual=float(amount_gradient @ derivative),
        differentiated_conservation_residual=float(np.max(np.abs(conservation_derivative))),
        distant_inactive_dilution_jacobian_entry=global_entry,
        rotated_rhs_residual=float(np.max(np.abs(rhs_rotation))),
        rotated_jacobian_residual=float(np.max(np.abs(jac_rotation))),
        uniform_native_full_derivative_residual=uniform_residual,
        polarity_vector_rotation_residual=float(vector_residual),
        zero_active_pool_vector_and_history_pass=True,
        seam_crossing_cue_C_derivative_residual=cue_derivative_error,
        scope="Operator and Jacobian identities only; no integration, spatial convergence, basin or biological calibration claim")


def smoke_protocol(output):
    """Small SLURM-only regression for the reusable solver's zero active pool.

    This tests API/runtime behaviour, not polarization or mesh resolution. The
    combined protocol must reproduce two separately restarted calls exactly.
    """
    closure = Closure(2.0, 24)
    initial = closure.initial()
    initial[:closure.cells] = 0.0
    initial[closure.cells:-1] = 2.0
    phases = [CuePhase(0.1, 4.0, closure.material_positions[0], 0.2, 0.05),
              CuePhase(0.2, output_interval=0.1)]
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        combined = run_protocol(closure, phases, initial_state=initial)
        first = run_protocol(closure, phases[:1], initial_state=initial)
        second = run_protocol(closure, phases[1:], initial_state=first["states"][-1])
    np.testing.assert_array_equal(combined["states"][0], initial)
    np.testing.assert_array_equal(combined["states"][:len(first["time"])], first["states"])
    np.testing.assert_array_equal(combined["states"][len(first["time"])-1:], second["states"])
    np.testing.assert_array_equal(combined["time"], [0.0, 0.05, 0.1, 0.1 + 0.1, 0.1 + 0.2])
    for name in ("polarity_x", "polarity_y", "polarity_magnitude"):
        assert combined[name][0] == 0.0
        assert np.isfinite(combined[name]).all()
    assert np.isnan(combined["polarity_direction_radians"][0])
    assert np.isfinite(combined["polarity_direction_radians"][1:]).all()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = output.with_suffix(".npz")
    np.savez_compressed(fields, **{key: value for key, value in combined.items()
                                 if isinstance(value, np.ndarray)})
    result = dict(generated_at=datetime.now(timezone.utc).isoformat(), passed=True,
        slurm_job_id=__import__("os").environ["SLURM_JOB_ID"], cells=closure.cells,
        exact_split_restart_continuity=True, zero_active_pool_vector=[0.0, 0.0],
        zero_active_pool_magnitude=0.0, zero_active_pool_direction=None,
        warnings="RuntimeWarnings promoted to errors; none occurred",
        max_saved_amount_relative_error=combined["max_saved_amount_relative_error"],
        minimum_saved_concentration=combined["minimum_saved_concentration"],
        full_field_archive=fields.name,
        full_field_sha256=hashlib.sha256(fields.read_bytes()).hexdigest(),
        source_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), Path(__file__).with_name("spring_ring_model_20261019.py"),
                         Path(__file__).with_name("spring_ring_closure_20261017.py"), Path(native.__file__))},
        scope="Zero-active initial API regression and exact staged restart; no preparation, polarity-stability or grid-resolution claim")
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--slurm-smoke-output", type=Path)
    args = parser.parse_args()
    result = smoke_protocol(args.slurm_smoke_output) if args.slurm_smoke_output else check_model()
    print(json.dumps(result, indent=2, allow_nan=False))
